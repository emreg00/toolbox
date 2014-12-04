import urllib2, os, cPickle, json
import configuration

try:
    CONFIG = configuration.Configuration() 
    LINCS_USER_KEY = CONFIG.get("LINCS_API_KEY") # change this value with custom API key
except:
    print "Warning: Using default LINCS API key"
    LINCS_USER_KEY = 'lincsdemo' 


def main():
    #drug = "InChIKey=PSGAAPLEWMOORI-PEINSRQWSA-N" 
    drug = "InChIKey=WIGIZIANZCJQQY-UHFFFAOYSA-N" # Glimepiride
    cell_line = None #"MCF7"
    get_drug_signature(drug, cell_line)
    return


def get_data(command, parameter, parameter2=None):
    if command == "list":
        parameter = parameter.replace(" ", "+")
        txt = 'pertinfo?q={"inchi_key":"%s"}' % (parameter)
        # example: InChIKey=PSGAAPLEWMOORI-PEINSRQWSA-N
    elif command == "get":
        if parameter2 is None:
            txt = 'siginfo?q={"pert_id":"%s"}' % (parameter)
        else:
            txt = 'siginfo?q={"pert_id":"%s","cell_id":"%s"}' % (parameter, parameter2)
        # example: BRD-K82216340
    else:
        raise ValueError("Unknown command: " + command)
    url = 'http://api.lincscloud.org/a2/%s&user_key=%s' % (txt, LINCS_USER_KEY)
    #print url
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    while True:
	try:
	    response = json.load(response)
	    break
	except:
	    print "Problem with response:", parameter, parameter2
	    response = urllib2.urlopen(req)
    return response


def get_drug_signature(drug, cell_line, probe_to_genes = None):
    #try:
    response = get_data("list", drug)
    #except urllib2.HTTPError:
    if len(response) > 1:
        print "Warning: multiple perturbations"
        for row in response:
            pert_id = row["pert_id"]
            print pert_id
    elif len(response) == 0:
        print "No info for", drug
    	return []
    pert_id = response[0]["pert_id"]
    response = get_data("get", pert_id, cell_line)
    values = []
    perturbation_times = []
    for row in response:
        #print row["pert_time"], row["cell_id"]
	if "pert_time" not in row:
	    continue
        perturbation_times.append(int(row["pert_time"]))
    if len(perturbation_times) > 0:
	perturbation_time = min(perturbation_times) 
    for row in response:
	if "pert_time" in row:
	    p_time = int(row["pert_time"])
	else:
	    p_time = None
        #if p_time is not None and perturbation_time != p_time: # to use only the shortet/longest perturbation(s) 
        #    continue 
	cell_type = row["cell_id"]
        values_up = row["up100_full"] 
        values_down = row["dn100_full"] 
        #values_up = row["up100_bing"] 
        #values_down = row["dn100_bing"] 
        #values_up = row["up50_lm"]
        #values_down = row["dn50_lm"]
        #print len(values_up), len(values_down)
	if probe_to_genes is not None:
	    genes_up = []
	    genes_down = []
	    for drug_probes, drug_genes in [(values_up, genes_up), (values_down, genes_down)]:
		for probe in drug_probes:
		    if probe in probe_to_genes:
			for gene in probe_to_genes[probe]:
			    if gene == "":
				continue
			    drug_genes.append(gene)
	else:
	    genes_up = values_up
	    genes_down = values_down
        values.append((genes_up, genes_down, cell_type, p_time))
    return values


def get_probe_mapping(probe_mapping_file, id_type="geneid"):
    """
	Parses trimmed HGU133A annotation file from GEO (GPL96)
	to get genes
	id_type: symbol / geneid
    """
    probe_to_genes = {}
    f = open(probe_mapping_file )
    line = f.readline()
    for line in f:
	probe, gene, geneid = line.strip("\n").split("\t")
	if id_type == "symbol":
	    for word in gene.split("///"):
		word = word.strip()
		probe_to_genes.setdefault(probe, []).append(word)
	elif id_type == "geneid":
	    words = geneid.split("///")
	    words.sort(lambda x,y: cmp(int(x), int(y)))
	    word = words[0]
	    probe_to_genes.setdefault(probe, []).append(word)
	else:
	    raise ValueError("Unknown id_type " + id_type)
    f.close()
    return probe_to_genes


def get_cmap_connectivity_score(disease_up, disease_down, drug_genes, measure_type="connectivity"):
    """
	average_rank: ranks of signature genes in profile
	connectivity: cmap connectivity score
	correlation: spearman correlation
    """
    s = None
    if measure_type == "connectivity":
        ks_up = get_cmap_ks_score(disease_up, drug_genes)
        ks_down = get_cmap_ks_score(disease_down, drug_genes)
        if ks_up >= 0 and ks_down >= 0:
            s = 0
        elif ks_up < 0 and ks_down < 0:
            s = 0
        else:
            s = ks_up - ks_down
    elif measure_type == "average_rank":
        ks_up = get_cmap_rank_score(disease_up, drug_genes)
        ks_down = get_cmap_rank_score(disease_down, drug_genes)
	if ks_up is None or ks_down is None:
	    s = 0.0
	else:
	    s = ks_down - ks_up 
    elif measure_type == "correlation":
        drug_genes_set = set(drug_genes) 
        rank_disease = []
        disease_genes = []
        disease_down.reverse()
        for i, gene in enumerate(disease_up + disease_down):
            if gene in drug_genes_set:
                rank_disease.append(i+1)
                disease_genes.append(gene)
        gene_to_rank = dict([ (gene, i+1) for i, gene in enumerate(drug_genes) ])
        rank_drug = [ gene_to_rank[gene] for gene in disease_genes ]
        if len(rank_drug) == 0 or len(rank_disease) == 0:
            return 0
        s, pval = stat_utilities.correlation(rank_disease, rank_drug, cor_type = "spearman")
	ks_up = pval
    else:
        raise ValueError("Unknown measure type! " + measure_type)
    return s, ks_up


def get_cmap_ks_score(genes_signature, genes_profile):
    """
        genes_signature: ordered list of genes in the signature, aka tags (based on correlation or differential expression)
        genes_profile: ordered list of (up and down) genes in the instance/perturbation profile
    """
    genes_signature = set(genes_signature)
    t = float(len(set(genes_profile) & genes_signature))
    n = float(len(genes_profile))
    j = 0
    values_a = []
    values_b = []
    for i, gene in enumerate(genes_profile):
        if gene in genes_signature:
            j += 1
            v_j = i + 1
            a = (j / t) - (v_j / n)
            b = (v_j / n) - ((j-1) / t) 
            values_a.append(a)
            values_b.append(b)
    if len(values_a) == 0 or len(values_b) == 0: 
        return 0.0
    a = max(values_a)
    b = max(values_b)
    if a >= b:
        ks = a
    else:
        ks = -b
    return ks


def get_cmap_rank_score(genes_signature, genes_profile):
    """
        genes_signature: ordered list of genes in the signature, aka tags (based on correlation or differential expression)
        genes_profile: ordered list of (up and down) genes in the instance/perturbation profile
    """
    genes_signature = set(genes_signature)
    values_e = []
    for i, gene in enumerate(genes_profile):
        if gene in genes_signature:
            v_j = i + 1
	    values_e.append(v_j) 
    if len(values_e) == 0: 
	ks = None
    else:
	ks = sum(values_e) / float(len(values_e)) 
    return ks 


if __name__ == "__main__":
    main()


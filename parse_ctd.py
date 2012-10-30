dir_name = "/home/emre/arastirma/netzcore/data/ctd/"

def main():
    #disease_to_values = parse_CTD()
    #get_inference_scores_for_disease(disease_to_values, "alzheimer")
    get_nc_scores_for_disease("alzheimer")
    get_nc_scores_for_disease("diabetes")
    get_nc_scores_for_disease("aids")
    return

def parse_CTD():
    disease_to_genes = {}
    disease_to_values = {}
    for line in open(dir_name + "CTD_genes_diseases.tsv"):
	if line.startswith("#"):
	    continue
	words = line.strip().split("\t")
	gene = words[0]
	disease = words[2].lower()
	disease = disease.split()[0]
	if len(disease.split("/"))>1:
	    disease = disease.split("/")[0]

	disease = disease.rstrip(",")
	disease_to_genes.setdefault(disease, set()).add(gene)

	is_direct = 0
	if words[4] != "":
	    is_direct = 1
	
	score = "NA"
	if words[6] != "":
	    score = float(words[6])

	disease_to_values.setdefault(disease, []).append((gene, score, is_direct))

    for disease, genes in disease_to_genes.iteritems():
	print disease, len(genes)
	#print len(set(zip(*disease_to_values[disease])[0]))
												      
    for disease, values in disease_to_values.iteritems():                                               
	genes = set()
	genes_direct = set()
	for gene, score, is_direct in values:
	    if is_direct == 1:
		genes_direct.add(gene)
	    genes.add(gene)
	f = open(dir_name + "associations/" + disease + "_direct.txt", 'w')
	for gene in genes_direct:
	    f.write("%s\n" % gene)
	f.close()
	f = open(dir_name + "associations/" + disease + ".txt", 'w')
	for gene in genes:
	    f.write("%s\n" % gene)
	f.close()
    return disease_to_values

def get_netcombo_scores(disease):
    node_to_score = {}
    for line in open("../data/output/biana_no_tap_relevance/new_omim_%s/nc3/node_scores.sif.genesymbol" % disease):
	gene, score = line.strip().split()
	node_to_score[gene] = score
    return node_to_score

def get_inference_scores_for_disease(disease_to_values, disease_name):
    node_to_score = get_netcombo_scores()
    f = open(dir_name + "%s.dat" % disease_name, 'w')
    f.write("gene i_score is_direct nc_score\n")
    for disease, values in disease_to_values.iteritems():                                               
	if disease.startswith(disease_name):                                                           
	    gene_to_score = {}
	    gene_to_direct = {}
	    for gene, score, is_direct in values:
		if gene in gene_to_score:
		    if is_direct == 1:
			gene_to_direct[gene] = is_direct
		    if gene_to_score[gene] != score:
			print "Warning: iscore inconsitency", gene, score, gene_to_score[gene]
			if gene_to_score[gene] == "NA":
			    gene_to_score[gene] = score
		else:
		    gene_to_direct[gene] = is_direct
		    gene_to_score[gene] = score

	    for gene, is_direct in gene_to_direct.iteritems():
		if gene in node_to_score:
		    f.write("%s %s %d %s\n" % (gene, str(gene_to_score[gene]), is_direct, node_to_score[gene]))
    f.close()
    return

def get_nc_scores_for_disease(disease):
    node_to_score = get_netcombo_scores(disease)
    if disease == "aids":
	disease = "hiv"
    genes = set([ line.strip() for line in open(dir_name + "associations/%s.txt" % disease) ])
    genes_direct = set([ line.strip() for line in open(dir_name + "associations/%s_direct.txt" % disease) ])
    f = open(dir_name + "%s_nc.dat" % disease, 'w')
    f.write("gene in_ctd is_direct nc_score\n")
    for gene, score in node_to_score.iteritems():
	in_ctd = 0
	is_direct = 0
	if gene in genes:                                               
	    in_ctd = 1
	if gene in genes_direct:
	    is_direct = 1
	f.write("%s %d %d %s\n" % (gene, in_ctd, is_direct, score))
    f.close()
    return
 
if __name__ == "__main__":
    main()


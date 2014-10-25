import urllib
import urllib2
from bs4 import BeautifulSoup
import cPickle

# Taken from Durek and Walter http://www.biomedcentral.com/1752-0509/2/100
CURRENCY_METABOLITES = set(["C00001", "C00002", "C00003", "C00004", "C00005", "C00006", "C00007", "C00008", "C00009", "C00010", "C00011", "C00013", "C00014", "C00015", "C00016", "C00018", "C00019", "C00020", "C00021", "C00023", "C00027", "C00028", "C00030", "C00034", "C00035", "C00038", "C00044", "C00050", "C00055", "C00061", "C00063", "C00070", "C00075", "C00076", "C00080", "C00105", "C00112", "C00113", "C00115", "C00120", "C00125", "C00126", "C00138", "C00139", "C00144", "C00175", "C00194", "C00205", "C00238", "C00291", "C01352"])


def main():
    out_file = "/home/emre/data/interactome/METABOLIC/network.txt"
    #output_interactions(out_file)
    #print get_disease_info("H00021")
    #print get_disease_info("H00410")
    #print get_drug_disease_mapping(["D06402"])
    return


def get_drug_disease_mapping(drugbank_to_kegg_id, dump_file):
    if os.path.exists(dump_file):
	drug_to_diseases = cPickle.load(open(dump_file))
	return drug_to_diseases 
    #! Get kegg_ids in drugbank parser
    #kegg_ids if multiple keggs reduce(lambda x,y: x | y, drugbank_to_kegg_id.values())
    kegg_ids = set(drugbank_to_kegg_id.values()) 
    kegg_drug_to_phenotype_and_mesh_ids = get_kegg_drug_mesh_mapping(kegg_ids)
    drug_to_diseases = {} # (mesh_id, mesh_term, may_treat) 
    #flag = False
    for drugbank_id, kegg_id in drugbank_to_kegg_id.iteritems():
	#if drugbank_id == "DB04575":
	#    flag = True
	#if flag == False: 
	#    continue
	print drugbank_id, kegg_id
        if kegg_id not in kegg_drug_to_phenotype_and_mesh_ids:
            continue
        values = kegg_drug_to_phenotype_and_mesh_ids[kegg_id]
        print values
	for disease, dui in values:
	    disease = disease.lower()
            val = 1
	    drug_to_diseases.setdefault(drugbank_id, []).append((disease, dui, val))
    cPickle.dump(drug_to_diseases, open(dump_file, 'w'))
    raise ValueError("Not implemented!")
    return



def get_data(command, parameter, parameter2=None):
    """
        get drug / disease: get/dr:D00001 or get/ds:H00001
        get all drug-disease links: link/drug/disease
        get all diseases: list/disease
    """
    url = 'http://rest.kegg.jp/%s/%s' % (command, parameter) #list/reaction
    if parameter2 is not None: # command == "link"
        url += '/%s' % parameter2
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    for line in response:
        yield line


def get_kegg_drug_mesh_mapping(kegg_ids = None):
    # If no drug is given get all drugs 
    if kegg_ids is None:
        kegg_ids = []
        for line in get_data("list", "drug"):
            drug, description = line.strip("\n").split("\t")
            kegg_ids.append(drug)
        print len(kegg_ids), kegg_ids[:5]
    # Get drug disease mapping
    kegg_drug_to_diseases = {}
    for line in get_data("link", "disease", "drug"):
        drug, disease = line.strip("\n").split()
        drug = drug[len("dr:"):] # strip dr:
        disease = disease[len("ds:"):] # strip ds:
        kegg_drug_to_diseases.setdefault(drug, set()).add(disease)
    print len(kegg_drug_to_diseases), kegg_drug_to_diseases.items()[:5]
    # Get corresponding mesh ids of diseases 
    diseases = reduce(lambda x,y: x | y, kegg_drug_to_diseases.values())
    print len(diseases), list(diseases)[:5]
    kegg_disease_to_phenotype_and_mesh_ids = {}
    for kegg_id in diseases:
        phenotype, mesh_ids = get_disease_info(kegg_id)
        kegg_disease_to_phenotype_and_mesh_ids[kegg_id] = (phenotype, mesh_ids)
    kegg_drug_to_phenotype_and_mesh_ids = {}
    # Map disease id to mesh for each drug
    for kegg_id in kegg_ids:
        if kegg_id not in kegg_drug_to_diseases:
            continue
        for kegg_disease in kegg_drug_to_diseases[kegg_id]:
            if kegg_disease not in kegg_disease_to_phenotype_and_mesh_ids:
                continue
            phenotype, mesh_ids = kegg_disease_to_phenotype_and_mesh_ids[kegg_disease]
            for mesh_id in mesh_ids:
                kegg_drug_to_phenotype_and_mesh_ids.setdefault(kegg_id, set()).add((phenotype, mesh_id))
    print len(kegg_drug_to_phenotype_and_mesh_ids), kegg_drug_to_phenotype_and_mesh_ids.items()[:5]
    return kegg_drug_to_phenotype_and_mesh_ids


def get_disease_info(kegg_disease):
    mesh_ids = None
    flag = False
    for line in get_data("get", kegg_disease):
        if line.startswith("NAME"):
            phenotype = " ".join(line.strip("\n").split()[1:])
            continue
        elif line.startswith("DBLINKS"):
            flag = True
        elif not line.startswith(" "):
            flag = False
            continue
        else:
            if flag == False:
                continue
        #print line,
        words = line.strip("\n").split()
        if words[0] == "DBLINKS":
            db = words[1]
            terms = words[2:]
        else:
            db = words[0]
            terms = words[1:]
        #print dummy, db, term
        if db == "MeSH:":
            mesh_ids = set(terms)
    return phenotype, mesh_ids


def output_interactions(out_file):
    reactions = get_reactions()
    reaction_to_values = {}
    for reaction_id in reactions:
        enzymes, lefts, rights = get_reaction_info(reaction_id)
        #print reaction_id, enzymes
        #print lefts, rights
        lefts = set(lefts) - CURRENCY_METABOLITES
        rights = set(rights) - CURRENCY_METABOLITES
        if len(lefts) + len(rights) == 0:
            continue
        geneids = set() 
        for enzyme in enzymes:
            geneids |= set(get_enzyme_info(enzyme))
        #print geneids
        if len(geneids) == 0:
            continue
        reaction_to_values[reaction_id] = (geneids, lefts, rights)
    reaction_ids = reaction_to_values.keys()
    print len(reaction_ids)
    interactions = set()
    for i, reaction_id1 in enumerate(reaction_ids):
        for j, reaction_id2 in enumerate(reaction_ids):
            if i <= j:
                (geneids1, lefts1, rights1) = reaction_to_values[reaction_id1]
                (geneids2, lefts2, rights2) = reaction_to_values[reaction_id1]
                if len(lefts1 & rights2) > 0 or len(lefts2 & rights1) > 0:
                    for geneid1 in geneids1:
                        for geneid2 in geneids2:
                            interactions.add(tuple(sorted((geneid1, geneid2))))
    print len(interactions)
    f = open(out_file, 'w')
    for geneid1, geneid2 in interactions:
        f.write("%s\t%s\n" % (geneid1, geneid2))
    f.close()
    return


def get_reactions():
    reactions = []
    for line in get_data("list", "reaction"):
        reaction_id = line.strip("\n").split()[0]
        reactions.append(reaction_id)
        #if len(reactions) > 200:
        #    break #!
    return reactions


def get_reaction_info(reaction_id):
    lefts = []
    rights = []
    enzymes = []
    for line in get_data("get", reaction_id):
        if line.startswith("EQUATION"):
            words = line.strip("\n").split()
            passed = False
            for word in words:
                if word.startswith("C"):
                    if passed:
                        rights.append(word)
                    else:
                        lefts.append(word)
                elif word == "<=>":
                    passed = True
        if line.startswith("ENZYME"):
            words = line.strip("\n").split()
            for word in words[1:]:
                enzymes.append(word)
    return enzymes, lefts, rights


def get_enzyme_info(ec_id, prefix="HSA:"):
    geneids = []
    try:
        for line in get_data("get", ec_id):
            if line.startswith("GENES") or line.lstrip().startswith(prefix):
                idx = line.find(prefix)
                if idx != -1:
                    words = line[idx+5:].strip().split()
                    for word in words:
                        e_idx = word.rfind("(")
                        geneid = word[:e_idx]
                        geneids.append(geneid)
    except urllib2.HTTPError:
        print "No info for", ec_id
        return geneids
    return geneids


if __name__ == "__main__":
    main()


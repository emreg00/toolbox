
import TsvReader
import re, os


def main():
    dir_name = "../data/"
    uniprot_dir_name = "uniprot_2011_Nov_2"
    mim_to_traits = get_mim_to_traits(dir_name + "disease/omim/morbidmap")
    print len(mim_to_traits)
    print mim_to_traits["600807"]
    #merge_uniprot_chromosome_files(uniprot_dir_name)
    #disease_to_genes, disease_to_loci = get_disease_gene_mapping(dir_name)
    #get_candidates_by_loci_matching(disease_to_genes, disease_to_loci, dir_name, uniprot_dir_name)
    #get_all_genes_in_morbidmap("diabetes_type_2_morbidmap.txt", "diabetes_type_2_genes.txt")
    #check_all_loci_in_morbidmap(dir_name + os.sep + "morbidmap")
    #get_disease_similarity_matrix(dir_name, dir_name + "similarity.dat")
    return

def get_mim_to_traits(file_name):
    mim_to_traits = {}
    trait_exp = re.compile("(.*), (\d{6})")
    f = open(file_name)
    for line in f:
	words = line.strip().split("|")
	#if words[0][0] == "?":
	#    continue
	m = trait_exp.match(words[0])
	if m:
	    mim = m.group(2)
	    trait = m.group(1).strip("{}")
	    print mim, trait
	else:
	    continue
	mim_to_traits.setdefault(mim, []).append(trait)
    f.close()
    return mim_to_traits


def get_disease_genes(base_dir, genes_to_be_considered=None, top_percentage=None):
    disease_to_genes = {}
    for file_name in os.listdir(base_dir):
	if not file_name.startswith("omim_"):
	    continue
	disease = file_name[5:-4]
	print disease,
	scores_and_genes = []
	f = open(base_dir + file_name)
	if top_percentage is None:
	    genes = set([ line.strip().split()[0] for line in f ])
	else:
	    for i, line in enumerate(f):
		words = line.strip().split()
		if len(words)>1:
		    gene = words[0]
		    score = float(words[1])
		else:
		    gene = words[0]
		    score = 9999999 - i
		scores_and_genes.append((score, gene))
	    scores_and_genes.sort()
	    scores_and_genes.reverse()
	    k = int(round(top_percentage * len(scores_and_genes) / 100.0))
	    genes = set(zip(*scores_and_genes[:k])[1])
	f.close()
	if genes_to_be_considered is not None:
	    genes = genes & genes_to_be_considered
	print len(genes)
	disease_to_genes[disease] = genes
    return disease_to_genes

def get_disease_similarity_matrix(base_dir, output_file, genes_to_be_considered=None, top_percentage=None):
    get_disease_genes(base_dir, genes_to_be_considered, top_percentage)
    diseases = disease_to_genes.keys()
    diseases.sort()
    f = open(output_file, "w")
    f.write("%s\n" % " ".join(diseases))
    for disease1 in diseases:
	f.write("%s" % disease1)
	for disease2 in diseases:
	    genes1 = disease_to_genes[disease1]
	    genes2 = disease_to_genes[disease2]
	    f.write(" %d" % len(genes1 & genes2))
	f.write("\n")
    f.close()
    return

def check_all_loci_in_morbidmap(file_name):
    f = open(file_name)
    for line in f:
	words = line.strip().split("|")
	tokenize_loci(words[3])
    f.close()
    return


def get_all_genes_in_morbidmap(file_name, out_file_name):
    f = open(file_name)
    genes = set()
    for line in f:
	words = line.strip().split("|")
	if words[0][0] == "?":
	    continue
	for word in words[1].split(", "):
	    genes.add(word.strip()) 
    f.close()
    f = open(out_file_name, 'w')
    [ f.write("%s\n" % gene) for gene in genes ]
    f.close()
    return

def merge_uniprot_chromosome_files(uniprot_dir_name):
    f_merged = open(uniprot_dir_name + os.sep + "uniprot_loci.txt", 'w')
    for file_name in sorted(os.listdir(uniprot_dir_name)):
	f = open(uniprot_dir_name + os.sep + file_name)
	prev_line = None
	flag = False
	for line in f:
	    if flag == False:
		if prev_line is not None and prev_line.strip().startswith("name") and line.strip().startswith("____"):
		    flag = True
	    elif line.strip() == "":
		break
	    else:
		words = line.split()
		f_merged.write("%s\t%s\n" % (words[0], words[1]))
		for word in words[2:]:
		    if word[0] == "[" and word[-1] == "]":
			f_merged.write("%s\t%s\n" % (word[1:-1], words[1]))
	    prev_line = line
	f.close()
    f_merged.close()
    return


def get_disease_gene_mapping(dir_name):

    network_genes = set([ line.strip() for line in open("bppi_all.txt") ])
    network_genes &= set([ line.strip() for line in open("entrez_all.txt") ])
    network_genes &= set([ line.strip() for line in open("goh_all.txt") ])

    f = open(dir_name + os.sep + "morbidmap")
    disorder_to_genes = {}
    disorder_to_loci = {}
    diseases = set()
    for line in f:
	words = line.strip().split("|")
	if words[0][0] == "?":
	    continue
	disorder = words[0].lstrip("[{").rstrip("}]")
	disorder = disorder.lower()
	disorder_words = disorder.split()
	if disorder_words[0] in ("von", "van", "h."):
	    disease = "".join(disorder_words[:2])
	else:
	    disease = disorder_words[0]
	diseases.add(disease)
	for word in words[1].split(", "):
	    disorder_to_genes.setdefault(disorder, set()).add(word.strip()) 
	disorder_to_loci.setdefault(disorder, set()).add(words[3].strip()) 
    f.close()

    disease_to_genes = {}
    disease_to_loci = {}
    for disorder, genes in disorder_to_genes.iteritems():
	#print disorder, genes
	disorder_words = disorder.split()
	if disorder_words[0] in ("von", "van", "h."):
	    disease = "".join(disorder_words[:2])
	else:
	    disease = disorder_words[0]
	disease = disease.replace("-","_").replace("/","_").rstrip(",")
	a = disease_to_genes.setdefault(disease, set())
	disease_to_genes[disease] = a.union(genes) 
	a = disease_to_loci.setdefault(disease, set())
	disease_to_loci[disease] = a.union(disorder_to_loci[disorder]) 
	for disorder_word in disorder_words:
	    for separator in ("-", "/"):
		words = disorder_word.split(separator)
		if len(words) > 1 and words[1] in diseases:
		    disease = words[1]
		    disease = disease.replace("-","_").replace("/","_").rstrip(",")
		    a = disease_to_genes.setdefault(disease, set())
		    disease_to_genes[disease] = a.union(genes) 
		    a = disease_to_loci.setdefault(disease, set())
		    disease_to_loci[disease] = a.union(disorder_to_loci[disorder]) 
    print disease_to_genes.keys()[:10]
    print len(disorder_to_genes), len(disease_to_genes)

    selected_disease_to_genes = {}
    selected_disease_to_loci = {}
    for disease, genes in disease_to_genes.iteritems():
	if len(network_genes&genes) > 5:
	    print disease, len(genes), len(network_genes&genes)
	    selected_disease_to_genes[disease] = genes #network_genes&genes
	    selected_disease_to_loci[disease] = disease_to_loci[disease]


    print selected_disease_to_loci.keys()[:10]
    print selected_disease_to_loci.values()[:10]
    print len(selected_disease_to_genes)

    diseases = selected_disease_to_genes.keys()
    diseases.sort()
    print diseases
    for disease in diseases:
	name = "new_omim_" + disease + ".txt" 
	f = open(dir_name + os.sep + name, 'w')
	genes = disease_to_genes[disease]
	[ f.write("%s\n" % g) for g in genes ]
	f.close()
    return selected_disease_to_genes, selected_disease_to_loci

def tokenize_loci(loci):
    ch, start, end = None, None, None
    chr_exp = re.compile("([\dXY]\d{0,1})([pq].*)")
    band_exp = re.compile("[pq]\d{1,2}(\.\d{1,2}){0,1}")
    words = loci.split("-")
    m = chr_exp.match(words[0])
    inconsistency = False
    if words[0][-3:] == "cen":
	ch = words[0][:-3]
	inconsistency = True
    elif m:
	ch = m.group(1)
	band = m.group(2)
	m2 = band_exp.match(band)
	if m2:
	    start = band
	else:
	    inconsistency = True
    else:
	inconsistency = True
    if len(words) > 1:
	band = words[1]
	m2 = band_exp.match(band)
	if m2:
	    end = band
	else:
	    inconsistency = True
    if start is None and end is not None:
	start = end 
    elif end is None and start is not None:
	end = start
    if inconsistency:
	if ch is None:
	    print "Ignoring:", loci, ch, start, end
    return ch, start, end

def get_genetic_interval(loci):
    import math
    loci_exp = re.compile("([\dXY]\d{0,1})(.*)")
    loci_exp2 = re.compile("(Chr\.){0,1}([\dXY]\d{0,1})$")
    ch, start, end = None, None, None
    m = loci_exp2.match(loci)
    if m:
	ch = m.group(2)
	return ch, start, end
    loci = loci.replace("ter","99")
    loci = loci.replace("cen","q0")
    m = loci_exp.match(loci)
    if m:
	ch = m.group(1) 
	band = m.group(2)
	if band == "":
	    return ch, start, end
	words = band.split("-")
	if len(words) > 1:
	    start = words[0]
	    end = words[1]
	    if end[0] == "p":
		end = float(end.replace("p", "-"))
	    elif end[0] == "q":
		end = float(end.replace("q", ""))
	else:
	    start = words[0]
	    if start == "q":
		start = "q0"
		end = 99.0
	    elif start == "p":
		start = "q0"
		end = -99.0
    else:
	print "Misformatted loci", loci
    start_org = start
    if start[0] == "p":
	start = float(start.replace("p", "-"))
    elif start[0] == "q":
	start = float(start.replace("q", ""))
    if end is None:
	if abs(start) > 10:
	    if start_org.find(".") > -1: 
		end = math.ceil(abs(start))-0.01
	    else: #if start == int(start):
		end = abs(start) + 0.99
	    if start < 0:
		end = -end
	else:
	    end = abs(start) + 0.99
	    if start < 0:
		end = -end
    if abs(start) > 10:
	start /= 10
    if abs(end) > 10:
	end /= 10
    if end < start:	
	temp = start
	start = end
	end = temp
    if start < 0:
	val = str(start)
	val = val[1:]
	if start == int(start):
	    val = val[:-2]
	if len(val) == 1:
	    start = float("-"+val+".999")
	elif len(val) == 3:
	    start = float("-"+val+"99")
	elif len(val) == 4:
	    start = float("-"+val+"9")
    if end > 0:
	val = str(end)
	if end == int(end):
	    val = val[:-2]
	if len(val) == 1:
	    end = float(val+".999")
	elif len(val) == 3:
	    end = float(val+"99")
	elif len(val) == 4:
	    end = float(val+"9")

    return ch, start, end

def check_loci_consistency(loci, uniprot_loci):
    inside = False
    if uniprot_loci == loci:
	inside = True
    else:
	try:
	    ch, start, end = get_genetic_interval(loci)
	    u_ch, u_start, u_end = get_genetic_interval(uniprot_loci)
	except:
	    print "Problem in getting interval:", loci, uniprot_loci
	    return None
	print ch, start, end
	print u_ch, u_start, u_end
	if u_ch is None:
	    print "u_ch is None", loci, uniprot_loci
	if ch == u_ch:
	    if start is None or u_start is None:
		inside = True
	    else:
		if (u_start >= start and u_start <= end):
		    inside = True
		elif (u_end >= start and u_end <= end):
		    inside = True
		elif (start >= u_start and start <= u_end):
		    inside = True
		elif (end >= u_start and end <= u_end):
		    inside = True
    return inside

def get_candidates_by_loci_matching(disease_to_genes, disease_to_loci, dir_name, uniprot_dir_name, via="biomart"):

    temp_file = ".biomart_results.txt"

    disease_to_candidates = {}
    if via == "biomart":
	# Checking all genes under the linkage interval of disease
	for disease, vals in disease_to_loci.iteritems():
	    for loci in vals:
		chromosome, start, end = tokenize_loci(loci)
		os.system("~/bin/R -f get_candidates.r --slave --args %s %s %s %s" % (chromosome, start, end, temp_file))
		for line in open(temp_file):
		    disease_to_candidates.setdefault(disease, set()).add(line.strip())
    elif via == "uniprot":
	gene_to_loci = {}
	f = open(uniprot_dir_name + os.sep + "uniprot_loci.txt")
	for line in f.readlines():
	    gene, loci = line.strip().split()
	    if gene.strip() == "-":
		continue
	    if gene_to_loci.has_key(gene):
		if loci != gene_to_loci[gene]:
		    print "Warning: inconsistent loci duplication for", gene, loci
		continue
	    loci = loci.replace(",",".")
	    loci = loci.rstrip(".")
	    gene_to_loci[gene] = loci
	f.close()

	# Checking all genes under the linkage interval of disease
	for disease, vals in disease_to_loci.iteritems():
	    for loci in vals:
		for gene, uniprot_loci in gene_to_loci.iteritems():
		    if check_loci_consistency(loci, uniprot_loci):
			disease_to_candidates.setdefault(disease, set()).add(gene)

    for disease, genes in disease_to_genes.iteritems():
	disease_to_candidates[disease] |= genes

    for disease, vals in disease_to_candidates.iteritems():
	f = open(dir_name + os.sep + "candidates/omim_%s.txt" % disease, "w")
	for gene in vals:
	    f.write("%s\n" % (gene)) 
	f.close()
    return


if __name__ == "__main__":
    main()
    #print check_loci_consistency("19q13.32", "19q13.13")
    #print check_loci_consistency("19p13.2", "19q13.13")
    #print check_loci_consistency("19q13.2", "19q13.13")



dir_name = "../data/gad/"

disease_to_genes = {}
f = open(dir_name + "all.txt")
f.readline()
f.readline()
f.readline()
for line in f:
    words = line.strip().split("\t")

    try:
	words[1]
    except:
	#print "Skipping:", line,
	continue

    if words[1] != "Y":
	continue
    gene = words[8]
    if gene == "":
	continue
    disease = words[2].lower()
    if disease == "":
	continue

    disease = disease.replace("|", "; ")
    diseases = disease.split("; ")

    try:
	disease = diseases[0].split()[0]
    except:
	#print "Skipping:", line,
	continue
    if disease.startswith("alzheimer"):
	disease = "alzheimer"
    if disease.startswith("aging"):
	disease = "aging"
    disease = disease.rstrip(",")
    if len(disease.split("/")) > 1:
	disease = disease.split("/")[0]
    disease_to_genes.setdefault(disease, set()).add(gene)

    if len(diseases) > 1:
	if diseases[1] == "":
	    continue
	disease = diseases[1].split()[0]
	if disease.startswith("alzheimer"):
	    disease = "alzheimer"
	if disease.startswith("aging"):
	    disease = "aging"
	disease = disease.rstrip(",")
	if len(disease.split("/")) > 1:
	    disease = disease.split("/")[0]
	disease_to_genes.setdefault(disease, set()).add(gene)

f.close()

for disease, genes in disease_to_genes.iteritems():
    print disease, len(genes)
    
for disease, genes in disease_to_genes.iteritems():
    f = open(dir_name + "associations/" + disease + ".txt", 'w')
    #if disease.startswith("alzheimer"):
    for gene in genes:
	f.write("%s\n" % gene)
    f.close()


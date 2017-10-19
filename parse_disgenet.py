import csv

def get_disgenet_genes(file_name):
    disease_to_genes = {}
    disease_to_sources = {} # not keeping sources of individual associations
    f = open(file_name)
    reader = csv.DictReader(filter(lambda row: row[0]!='#', f), delimiter='\t')
    for row in reader:
	disease = row["diseaseName"].lower()
	sources = disease_to_sources.setdefault(disease, set())
	sources |= set(row["sourceId"].lower().split(","))
	disease_to_genes.setdefault(disease, set()).add(row["geneId"])
    f.close()
    return disease_to_genes, disease_to_sources



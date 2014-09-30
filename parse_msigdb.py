
def main():
    file_name = "/home/emre/arastirma/data/ontology/msigdb/c2.all.v4.0.entrez.gmt"
    pathway_to_geneids, geneid_to_pathways = get_msigdb_info(file_name, prefix=None)
    print pathway_to_geneids["biocarta_alk_pathway"]
    return


def get_msigdb_info(pathway_file, prefix=None):
    pathway_to_geneids = {}
    geneid_to_pathways = {}
    for line in open(pathway_file):
	words = line.strip().split("\t")
	pathway = words[0].lower()
	if prefix is not None:
	    if not pathway.startswith(prefix):
		continue
	geneids = words[2:]
	pathway_to_geneids[pathway] = set(geneids)
	for geneid in geneids:
	    geneid_to_pathways.setdefault(geneid, set()).add(pathway)
    return pathway_to_geneids, geneid_to_pathways


if __name__ == "__main__":
    main()


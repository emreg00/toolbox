
def main():
    base_dir = "/home/emre/arastirma/data/drug/cmap"
    desc_file = base_dir + "/cmap_instances_02.csv"
    matrix_file = base_dir + "/rankMatrix.txt"
    probe_mapping_file = base_dir + "/probe_mapping.txt"
    drug_to_top_geneids = get_cmap_info(desc_file, probe_mapping_file, matrix_file)
    return

def get_cmap_info(desc_file, probe_mapping_file, matrix_file):
    # Get instance info
    f = open(desc_file)
    line = f.readline()
    words = line.strip().split("\t")
    header_to_idx = dict((word.lower(), i) for i, word in enumerate(words))
    drugs = set()
    instance_to_values = {} # drug, concentration, duration, cell_type
    drug_to_instances = {} 
    for line in f:
	words = line.strip().split("\t")
	drug = words[header_to_idx["cmap_name"]].lower()
	instance = words[header_to_idx["instance_id"]]
	concentration = float(words[header_to_idx["concentration (m)"]])
	duration = float(words[header_to_idx["duration (h)"]])
	cell_type = words[header_to_idx["cell2"]]
	instance_to_values[instance] = (drug, concentration, duration, cell_type)
	drug_to_instances.setdefault(drug, []).append(instance)
    f.close()
    # Get instance with highest concentration & shortest duration
    instance_to_drug = {}
    for drug, instances in drug_to_instances.iteritems():
	max_concentration = None
	min_duration = None
	selected_instance = None
	for instance in instances:
	    drug_, concentration, duration, cell_type = instance_to_values[instance]
	    if max_concentration is None:
		max_concentration = concentration
	    if min_duration is None:
		min_duration = duration
	    if duration < min_duration:
		min_duration = duration
		selected_instance = instance
	    elif duration == min_duration:
		if concentration >= max_concentration:
		    max_concentration = concentration
		    selected_instance = instance
	if selected_instance is None:
	    raise ValueError("None instance %s %s" % (drug, instances))
	instance_to_drug[selected_instance] = drug
    # Get probe geneid mapping
    f = open(probe_mapping_file)
    line = f.readline()
    probe_to_geneid = {}
    for line in f:
	probe, geneid = line.strip().split("\t")
	probe_to_geneid[probe] = geneid
    f.close()
    # Get expression ranks
    f = open(matrix_file)
    line = f.readline()
    words = line.strip().split("\t")
    #header_to_idx = dict((word.lower(), i) for i, word in enumerate(words))
    header_values = words
    instance_to_ranks = {}
    geneids = []
    for line in f:
	words = line.strip().split("\t")
	probe = words[0]
	if probe not in probe_to_geneid:
	    geneid = "*%s" % probe
	else:
	    geneid = probe_to_geneid[probe]
	geneids.append(geneid)
	for rank, instance in zip(words[1:], header_values):
	    instance_to_ranks.setdefault(instance, []).append(int(rank))
    f.close()
    geneids = numpy.array(geneids)
    # Get top geneids for each instance
    drug_to_top_geneids = {}
    for instance, ranks in instance_to_ranks.iteritems():
	if instance not in instance_to_drug:
	    continue
	indices = numpy.argsort(ranks)
	up_geneids = []
	down_geneids = []
	i = 0
	for geneid in geneids[indices][:250]:
	    if i >= 50:
		break
	    if geneid[0] == "*":
		continue
	    up_geneids.append(geneid)
	    i += 1
	i = 0
	for geneid in reversed(geneids[indices][-250:]):
	    if i >= 50:
		break
	    if geneid[0] == "*":
		continue
	    down_geneids.append(geneid)
	    i += 1
	drug_to_top_geneids[instance_to_drug[instance]] = (up_geneids, down_geneids)
    return drug_to_top_geneids


if __name__ == "__main__":
    main()


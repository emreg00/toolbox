
def main():
    import network_utilities
    module_file = "modules.txt"
    g =  network_utilities.create_network_from_sif_file("../data/interactions.sif")
    get_modules_of_graph(g, "mcl", output_file=module_file, inflation=1.7)
    return

def get_modules_from_file(output_file):
    f = open(output_file)
    modules = []
    for line in f:
	words = line.strip().split("\t")
	modules.append(words)
    f.close()
    return modules

def get_modules_of_graph(sub_graph, module_detection_type, output_file=None, inflation=1.7):
    if module_detection_type == "connected":
	import network_utilities
	modules = network_utilities.get_connected_components(sub_graph, return_as_graph_list=True)
    elif module_detection_type == "mcl":
	from os import system
	if output_file is None:
	    temp_file_name = ".temp_module_file.txt89734234"
	else:
	    temp_file_name = output_file
	f = open(temp_file_name, 'w')
	nodes = set()
	for node1, node2, data in sub_graph.edges(data=True):
	    nodes.add(node1)
	    nodes.add(node2)
	    f.write("%s\t%s\t%f\n" % (node1, node2, data))
	for node in sub_graph.nodes():
	    if node not in nodes:
		f.write("%s\n" % node)
	f.close()
	# Optimum inflation parameter was 1.7-1.8 in a recent comparison paper
	system("mcl %s --abc -I %f -o %s 2>> %s" % (temp_file_name, inflation, temp_file_name + ".mcl", temp_file_name + ".err"))
	f = open(temp_file_name + ".mcl")
	modules = []
	for line in f:
	    words = line.strip().split("\t")
	    modules.append(words)
	f.close()
    else:
	raise ValueError("Unrecognized module detection type")
    #print len(modules), map(len, modules)
    return modules

if __name__ == "__main__":
    main()



import argparse
import network_utilities, wrappers

def main():                         
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--network_file') #, required=True) 
    parser.add_argument('-s', '--nodes_from') #, required=True) 
    parser.add_argument('-t', '--nodes_to') #, required=True) 
    parser.add_argument('-o', '--out_file') #, required=True)
    parser.add_argument('-n', '--n_random', type=int, default=1000) 
    parser.add_argument('-m', '--min_bin_size', type=int, default=100) 
    parser.add_argument('-x', '--n_seed', type=int, default=452456) 
    parser.add_argument('-f', '--parameter_file', type=str, default=None) 
    args = parser.parse_args()
    if args.parameter_file is not None:
	nodes_from = None
	nodes_to = None
	network_file = None 
	n_random = None
	min_bin_size = None
	n_seed = None
	out_file = None
	arguments = open(args.parameter_file).read().split()
	#print arguments
	i = 0
	if arguments[0] == "":
	    i = 1
	while i < len(arguments):
	    if arguments[i] == "-e":
		network_file = arguments[i+1]
	    if arguments[i] == "-s":
		nodes_from = arguments[i+1].split(",")
	    if arguments[i] == "-t":
		nodes_to = arguments[i+1].split(",")
	    if arguments[i] == "-o":
		out_file = arguments[i+1]
	    if arguments[i] == "-m":
		min_bin_size = int(arguments[i+1])
	    if arguments[i] == "-n":
		n_random = int(arguments[i+1])
	    if arguments[i] == "-x":
		n_seed = int(arguments[i+1])
	    #print arguments[i], arguments[i+1]
	    i += 2
    else:
	nodes_from = args.nodes_from.split(",")
	nodes_to = args.nodes_to.split(",")
	network_file = args.network_file
	n_random = args.n_random
	min_bin_size = args.min_bin_size
	n_seed = args.n_seed
	out_file = args.out_file
    network = network_utilities.create_network_from_sif_file(network_file, use_edge_data = False, delim = None, include_unconnected=True)
    #print args
    print network_file, nodes_from, nodes_to, n_random, min_bin_size, n_seed, out_file
    values = wrappers.calculate_proximity(network, nodes_from = nodes_from, nodes_to = nodes_to, n_random = n_random, min_bin_size = min_bin_size, seed = n_seed)
    if values is not None: # not in network
	d, z, (m, s) = values
	#print z, d, (m, s)
	open(out_file, 'w').write("%f %f %f %f\n" % (z, d, m, s))
    return


if __name__ == "__main__":
    main()


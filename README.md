# toolbox

Toolbox is a repository encapsulating various scripts used in my research on the analysis of disease and drug related biological data sets. 
It contains generic utilities for data processing (e.g., parsing, network-based analysis, etc, ...). 

The code here has been developed during the analysis of data in various projects such as
- [BIANA](http://github.com/emreg00/biana) ([@javigx2](https://twitter.com/javigx2) was the lead developer)
- [GUILD](http://github.com/GUILD)
<- PEPPER: PErsonalized Perturbation ProfilER>
- Proximity: A method to calculate distances between two groups of nodes in the network while correcting for degree biases (e.g., incompleteness or study bias).

The package mainly consists of two types of files:
- parser_{resource_name_to_be_parsed}.py
- {type_of_data/software}_utilities.py

For instance, [parse_drugbank.py](parse_drugbank.py) contains methods to parse DrugBank data base (v.3) XML dump 
and [network_utilities.py](network_utilities.py) contains methods related to network generation and analysis. 

### Proximity

For calculating proximity, the inputs are: node_set_1, node_set_2 and network. 
The nodes in the network are binned such that the nodes in the same bin have similar degrees. 
Next, random nodes matching the number and the degree of the nodes in the node sets are chosen.
The average distance from the nodes in one set to the other is then calculated and compared to the 
random expectation (the distances observed in random groups).
Below are the relevant methods in [network_utilities.py](network_utilities.py) for calculating the proximity.

- Creating node bins by degree:

  bins = get_degree_binning(network, min_bin_size)
  
- Selecting nodes randomly matching the degrees of the nodes in the given set:

  nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware)

- Calculating average distance from the nodes in nodes_from to the closest node in the nodes_to (lengths is the dictionary containing all pairwise shortest path distances, distance="closest", parameters={}):

  d = get_separation(network, lengths, nodes_from, nodes_to, distance, parameters)





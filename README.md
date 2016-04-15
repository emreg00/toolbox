# toolbox

Toolbox is a repository encapsulating various scripts used in my research on the analysis of disease and drug related biological data sets. 
It contains generic utilities for data processing (e.g., parsing, network-based analysis, etc, ...). 

The code here has been developed during the analysis of data in various projects such as
- [BIANA](http://github.com/emreg00/biana) ([@javigx2](https://twitter.com/javigx2) was the lead developer)
- [GUILD](http://github.com/emreg00/guild)
- Proximity: A method to calculate distances between two groups of nodes in the network while correcting for degree biases (e.g., incompleteness or study bias).

The package mainly consists of two types of files:
- parser_{resource_name_to_be_parsed}.py
- {type_of_data/software}_utilities.py

For instance, [parse_drugbank.py](parse_drugbank.py) contains methods to parse DrugBank data base (v.3) XML dump 
and [network_utilities.py](network_utilities.py) contains methods related to network generation and analysis. 

## Proximity

### Proximity analysis
To replicate the analysis in the paper please refer to [proximity](http://github.com/emreg00/proximity) repository.

### Proximity calculation

See `calculate_proximity` method in [wrappers.py](wrappers.py)  for calculating proximity:

`calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, n_random=1000, min_bin_size=100, seed=452456)`

For instance, to calculate the proximity from (A, C) to (B, D, E) in a toy network (given below):

```python
>>> from toolbox import network_utilities, wrappers
>>> file_name = "toy.sif"
>>> network = wrappers.get_network(file_name, only_lcc = True)
>>> nodes_from = ["A", "C"]
>>> nodes_to = ["B", "D", "E"]
>>> d, z, (mean, sd) = wrappers.calculate_proximity(network, nodes_from, nodes_to, min_bin_size = 2)
>>> print (d, z, (mean, sd))
(1.0, 0.97823676194805476, (0.75549999999999995, 0.24993949267772786))
>>>
```

Toy network (toy.sif):
```
A 1 B
A 1 C
A 1 D
A 1 E
A 1 F
A 1 G
A 1 H
B 1 C
B 1 D
B 1 I
B 1 J
C 1 K
D 1 E
D 1 I
E 1 F
```

The inputs are the two groups of nodes and the network. 
The nodes in the network are binned such that the nodes in the same bin have similar degrees. 
For real networks, use a larger `min_bin_size` (e.g., 100). 
The random nodes matching the number and the degree of the nodes in the node sets are chosen
using these bins.
The average distance from the nodes in one set to the other is then calculated and compared to the 
random expectation (the distances observed in random groups).


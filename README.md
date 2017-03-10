# toolbox

Toolbox is a repository encapsulating various scripts used in my research on the analysis of disease and drug related biological data sets. 
It contains generic utilities for data processing (e.g., parsing, network-based analysis, proximity, etc, ...).

Contents 

* [Background](#background)
* [Parsers](#parsers)
* [Wrappers](#wrappers)
* [Proximity](#proximity)
* [Citation](#citation)


## Background

The code here has been developed during the analysis of data in various projects such as
- [BIANA](http://github.com/emreg00/biana) ([@javigx2](https://twitter.com/javigx2) was the lead developer): Biological data analysis and network integration
- [GUILD](http://github.com/emreg00/guild): Network-based disease-gene prioritization
- [Proximity](http://github.com/emreg00/proximity): A method to calculate distances between two groups of nodes in the network while correcting for degree biases (e.g., incompleteness or study bias)

The package mainly consists of two types of files:
- parser_{resource_name_to_be_parsed}.py
- {type_of_data/software}_utilities.py

For instance, [parse_drugbank.py](parse_drugbank.py) contains methods to parse DrugBank data base (v.3) XML dump 
and [network_utilities.py](network_utilities.py) contains methods related to network generation and analysis. 


## Parsers

Parsers available for the following resources:

- [Clinical trials](http://clinicaltrials.gov)
- [CMAP](https://www.broadinstitute.org/cmap)
- [CTD](http://ctdbase.org)
- [DailyMed](http://dailymed.nlm.nih.gov)
- [DrugBank](http://www.drugbank.ca)
- [Drugrepurposing.info](http://drugrepurposing.info)
- [GAD](https://geneticassociationdb.nih.gov)
- [GeneCards](http://www.genecards.org)
- [GDSC](http://www.cancerrxgene.org)
- [GO](http://geneontology.org)
- [KEGG API](www.genome.jp/kegg)
- [KEGG Brite](http://www.genome.jp/kegg/brite.html)
- [LabeledIn](http://www.ncbi.nlm.nih.gov/pubmed/25220766)
- [LINCS API](http://www.lincsproject.org)
- [MEDI](https://medschool.vanderbilt.edu/cpm/center-precision-medicine-blog/medi-ensemble-medication-indication-resource)
- [Metab2Mesh](http://metab2mesh.ncibi.org)
- [MeSH](http://www.ncbi.nlm.nih.gov/mesh)
- [MSigDB](http://www.broadinstitute.org/msigdb)
- [NCBI Gene](http://www.ncbi.nlm.nih.gov/gene)
- [NDFRT API](https://rxnav.nlm.nih.gov/NdfrtAPIs.html)
- [OpenFDA API](https://open.fda.gov)
- [Orphanet](http://www.orpha.net)
- [OMIM](http://www.omim.org)
- [SIDER](http://sideeffects.embl.de)
- [SNOMEDCT](https://www.nlm.nih.gov/healthit/snomedct)
- [STITCH](http://stitch.embl.de)
- [STRING](http://string-db.org)
- [UMLS](https://www.nlm.nih.gov/research/umls)
- [UniProt](http://www.uniprot.org)

The parsers are provided "as is" and might not work due to updates on the data format of these resources. Please contact me for suggestions, bug reports and enquiries.


## External packages (Optional)

Some functions in toolbox rely on the following packages. The package will load properly but certain functionality might not be available.

* [GenRev](https://bioinfo.uth.edu/GenRev.html) for Steiner-tree algorithm implementation
* [negex](https://github.com/chapmanbe/negex) for identifying negatiation in textual medical records
* [FuncAssoc client](https://gist.github.com/yy/869845) for connecting FuncAssociate server (incorporated in the repository)


## Wrappers

[wrappers.py](wrappers.py) provides an easy to use interface to various methods I commonly use. It is continuously under development. Currently it contains methods to 
- Map UniProt, ENTREZ ids and gene symbols
- Creating networkx network from file 
- Calculating proximity
- Calculating functional enrichment using [FuncAssociate](http://llama.mshri.on.ca/funcassociate/) API


## GUILD

See below for python interface to run GUILD (assumes it is properly compiled and accessible at executable_path) using A and C as seeds and a toy network:

```python
>>> from toolbox import wrappers
>>> file_name = "toy.sif"
>>> network = wrappers.get_network(file_name, only_lcc = True)
>>> nodes = set(network.nodes())
>>> seeds = ["A", "C"]
>>> node_to_score = dict((node, 1) for node in seeds)
>>> name = "sample_run"
>>> output_dir = "./"
>>> wrappers.run_guild(name, node_to_score, nodes, file_name, output_dir, executable_path)
```

After this command input node score file "sample_run.node" and output node score file "sample_run.ns" will be created in the current directory.


## Proximity

### Proximity analysis
To replicate the analysis in the paper please refer to [proximity](http://github.com/emreg00/proximity) repository.

### Proximity calculation

See `calculate_proximity` method in [wrappers.py](wrappers.py)  for calculating proximity:

`calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, n_random=1000, min_bin_size=100, seed=452456)`

For instance, to calculate the proximity from (A, C) to (B, D, E) in a toy network (given below):

```python
>>> from toolbox import wrappers
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
The proximity is not symmetric (if nodes_from and nodes_to are swapped, 
the results would be different, see below for details).
The nodes in the network are binned such that the nodes in the same bin have similar degrees. 
For real networks, use a larger `min_bin_size` (e.g., 10, 25, 50, 100, see below for choosing the bin size). 
The random nodes matching the number and the degree of the nodes in the node sets are chosen
using these bins.
The average distance from the nodes in one set to the other is then calculated and compared to the 
random expectation (the distances observed in random groups).

### Proximity calculation considerations

* From/to nodes: Note that proximity, by definition, is not symmetric and the order of `nodes_from` and 
`nodes_to` makes a difference. If you do not have an intrinsic relationship between 
the two sets of nodes (the proximity from one to the other, such as from drugs to diseases),
you can use the node set smaller in size.

* Degree binning: For random selection of nodes with similar degrees to those in the original node sets, 
proximity uses binning of similar degree-nodes. This is because, high-degree nodes in the network are
less common and the randomization algorithm would choose always the same set of nodes if binning is not used.
That being said, if the bins are too large, the nodes within the bin are not a good representative of 
the original nodes (spanning many nodes with different degrees). Accordingly, the bins should contain enough
number of nodes that would allow a representative random sampling. In each bin the nodes with
degree higher degree (e.g., k+1, if k is the degree of the nodes in the current bin) are added iteratively
till `min_bin_size` is reached. For instance, `min_bin_size` can be chosen to be at least twice 
as large as the max number of nodes in the node sets.


## Citation

* If you use biomedical data base parsers or proximity related methods please cite: Guney E, Menche J, Vidal M, Barab&aacute;si AL. Network-based in silico drug efficacy screening. Nat. Commun. 7:10331 doi: 10.1038/ncomms10331 (2016). [link](http://www.nature.com/ncomms/2016/160201/ncomms10331/full/ncomms10331.html)

* If you use GUILD related methods please cite: 
Guney E, Oliva B. Exploiting Protein-Protein Interaction Networks for Genome-Wide Disease-Gene Prioritization. PLoS ONE 7(9): e43557 (2012). [link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043557)


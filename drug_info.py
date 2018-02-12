from toolbox import configuration, wrappers, stat_utilities 
from toolbox import parse_drugbank, parse_sider_v4, parse_stitch
from toolbox import TsvReader, network_utilities
import os, cPickle, numpy, random #, time,
try:
    from indigo import indigo
except:
    print "Indigo not found, chemical similarity will not be available!" 


def get_drug_info(parameters, parser=None):
    # Get drugbank data
    if parser is None:
	parser = get_drugbank(parameters)
    # Get drug target geneids for selected drugs (e.g., approved)
    drug_to_geneids = get_drug_targets(parameters, parser)
    selected_drugs = get_drugs_by_type(parameters, parser)
    drug_to_geneids = dict((drug, geneids) for drug, geneids in drug_to_geneids.iteritems() if drug in selected_drugs)
    return drug_to_geneids


def get_drugbank(parameters):
    file_name = parameters.get("drugbank_file")
    dump_file = file_name + ".pcl"
    if os.path.exists(dump_file):
	parser = cPickle.load(open(dump_file))
    else:
	parser = parse_drugbank.DrugBankXMLParser(file_name)
	parser.parse()
	cPickle.dump(parser, open(dump_file, 'w'))
    return parser


###### Drugbank related ######

def get_drug_targets(parameters, parser=None, id_type="geneid"):
    """
    id_type: uniprot | geneid | symbol
    """
    if parser is None:
	parser = get_drugbank(parameters)
    target_types = set(parameters.get("target_type").split("|"))
    drug_to_uniprots = parser.get_targets(target_types, only_paction=parameters.get_boolean("only_paction"))
    #print len(drug_to_uniprots), drug_to_uniprots.items()[:5]
    uniprot_ids = reduce(lambda x,y: x|y, drug_to_uniprots.values())
    if id_type == "uniprot":
	drug_to_targets = drug_to_uniprots
    elif id_type == "geneid":
	uniprot_to_id = wrappers.get_uniprot_to_geneid(parameters.get("uniprot_file"), uniprot_ids)
    elif id_type == "symbol":
	uniprot_to_id = wrappers.get_uniprot_to_symbol(parameters.get("uniprot_symbol_file"), uniprot_ids)
    else:
	raise ValueError("Unknown id type: %s" % id_type)
    drug_to_targets = {}
    for drug, uniprots in drug_to_uniprots.iteritems():
	for uniprot in uniprots:
	    if uniprot in uniprot_to_id:
		drug_to_targets.setdefault(drug, set()).add(uniprot_to_id[uniprot])
    #print len(drug_to_targets), drug_to_targets.items()[:5]
    #print drug_to_uniprots["DB00536"]
    #print drug_to_targets["DB00536"]
    return drug_to_targets


def get_drugbank_ids_for_drugs(drugs, parameters=None, parser=None, check_synonyms=True, use_text_matching=False):
    if parser is None:
	if parameters is None:
	    raise ValueError("One of the parser or parameters arguments are required!")
	parser = get_drugbank(parameters)
    name_to_drug, synonym_to_drug = parser.get_synonyms(selected_drugs=None, only_synonyms=False)
    drug_to_db_id = {}
    not_in_db = set()
    for drug in drugs:
	name = drug.lower()
	if name in name_to_drug:
	    drug = name_to_drug[name]
	    drug_to_db_id[name] = drug
	elif check_synonyms:
	    if name in synonym_to_drug:
		drug = synonym_to_drug[name]
		if name in drug_to_db_id:
		    print "Ignoring synonym for already matched id", name, drug_to_db_id[name], drug
		else:
		    drug_to_db_id[name] = drug
	if name not in drug_to_db_id:
	    if use_text_matching == True:
		found = False
		for drug, db_name in parser.drug_to_name.iteritems():
		    db_name = db_name.lower()
		    if db_name.find(name) != -1:
			found = True
			drug_to_db_id[db_name] = drug
			print "Id found by text matching", name, drug
		if not found:
		    not_in_db.add(name)
	    else:
		not_in_db.add(name)
    print "Not in DrugBank:", not_in_db
    return drug_to_db_id


def get_drugs_by_type(parameters, parser=None):
    if parser is None:
	parser = get_drugbank()
    drug_type = parameters.get("drug_type")
    groups_to_include = set([drug_type])
    groups_to_exclude = set([]) #["withdrawn"])
    if drug_type == "all":
	selected_drugs = set(parser.drug_to_name.keys())
	return selected_drugs 
    elif drug_type == "all_but_withdrawn":
	groups_to_include = set(["approved", "experimental", "investigational", "nutraceutical", "illicit"]) # , "withdrawn"
	groups_to_exclude = set(["withdrawn"])
    elif drug_type == "withdrawn":
	groups_to_exclude = set([])
    selected_drugs = parser.get_drugs_by_group(groups_to_include = groups_to_include, groups_to_exclude = groups_to_exclude) 
    return selected_drugs 


def get_drug_smiles_by_target(parameters, parser=None):
    if parser is None:
	parser = get_drugbank(parameters)
    drug_to_geneids = get_drug_targets(parameters, parser)
    smiles_to_geneids = {}
    geneid_to_smile_strings = {}
    for drug, smiles in parser.drug_to_smiles.iteritems():
	if drug in drug_to_geneids:
	    geneids = drug_to_geneids[drug]
	    for geneid in geneids:
		for words in (smiles.split("\n"), smiles.split("<br"), smiles.split(";"), smiles.split()):
		    if len(words) > 1:
			print "Potential multiple smiles", smiles 
		    elif len(smiles) == 0:
			print "Empty smiles", smiles 
		smiles_to_geneids.setdefault(smiles, set()).add(geneid)
		geneid_to_smile_strings.setdefault(geneid, set()).add(smiles)
    return smiles_to_geneids, geneid_to_smile_strings


def get_drug_drug_interactions(parameters, drug_names, out_file=None):
    parser = drug_info.get_drugbank(parameters)
    drug_to_db_id = drug_info.get_drugbank_ids_for_drugs(drug_names, parameters=None, parser=parser, check_synonyms=True, use_text_matching=False)
    db_id_to_interactions = parser.drug_to_interactions
    #print len(db_id_to_interactions), db_id_to_interactions.items()[:3]
    print len(drug_names), len(drug_to_db_id)
    drugs_new = []
    for drug in drug_names:
	if drug not in drug_to_db_id:
	    continue
	db_id = drug_to_db_id[drug]
	if db_id not in db_id_to_interactions:
	    continue
	drugs_new.append(drug)
    drugs_new.sort()
    if out_file is not None:
	f = open(out_file, 'w')
	f.write("Drug\t%s\n" % "\t".join(drugs_new))
	for drug in drugs_new:
	    values = ["0"] * len(drugs_new)
	    for i, drug2 in enumerate(drugs_new):
		if drug != drug2:
		    if drug_to_db_id[drug2] in db_id_to_interactions[drug_to_db_id[drug]]:
			values[i] = "1"
	    f.write("%s\t%s\n" % (drug, "\t".join(values)))
	f.close()
    return drug_to_db_id, db_id_to_interactions


###### Drug similarity related ######

def get_smiles_similarity(smiles1, smiles2, fp_type = "sim", metric = "tanimoto"):
    """
    fp_type: sim | sub
    metric: tanimoto | tversky
    """
    if len(smiles1) == 0 or len(smiles2) == 0:
	return None
    ind = indigo.Indigo()
    m = ind.loadMolecule(smiles1) 
    m.aromatize()
    fp = m.fingerprint(fp_type)
    m2 = ind.loadMolecule(smiles2) 
    m2.aromatize() # Aromatize molecules in case they are not in aromatic form
    fp2 = m2.fingerprint(fp_type) # Calculate similarity between "similarity" fingerprints
    d = ind.similarity(fp, fp2, metric) 
    return d


def get_target_similarity(targets1, targets2, target_to_occurrences=None):
    """
    Weighted jaccard, if target_to_occurrences (diseases / side effects) is not None
    """
    if len(targets1) == 0 or len(targets2) == 0:
	return None
    targets_common = targets1 & targets2
    if target_to_occurrences is None:
	#d = len(targets_common) / float(max(len(targets1), len(targets2))) # ~worse
	d = len(targets_common) / float(len(targets1|targets2)) 
    else:
	d = 0.0
	for t in targets_common:
	    d += 1.0/len(target_to_occurrences[t])
	d /= len(targets1|targets2) # max(len(targets1), len(targets2))
    return d


def get_target_ppi_similarity(targets1, targets2, network):
    if len(targets1) == 0 or len(targets2) == 0:
	return None
    vals = []
    for target1 in targets1:
	for target2 in targets2:
	    d = network_utilities.get_shortest_path_length_between(network, target1, target2)
	    vals.append(d)
    d = numpy.exp(-numpy.mean(vals))
    return d


def get_drug_similarity(drug_to_values, method="target", network=None, dump_file=None):
    """
    values = targets or smiles
    """
    if dump_file is not None and os.path.exists(dump_file):
	drug_to_drug_similarity = cPickle.load(open(dump_file))
	return drug_to_drug_similarity 
    if network is None and method == "target-ppi":
	raise ValueError("Network is required for target-ppi")
    drug_to_drug_similarity = {}
    drugs = drug_to_values.keys()
    for i, drug1 in enumerate(drugs):
	for j, drug2 in enumerate(drugs):
	    if i >= j:
		continue
	    #comb = tuple(sorted([drug1, drug2]))
	    val1 = drug_to_values[drug1]
	    val2 = drug_to_values[drug2]
	    d = None
	    if method == "target":
		d = get_target_similarity(val1, val2)
	    elif method == "target-ppi":
		d = get_target_ppi_similarity(val1, val2, network)
	    elif method == "chemical":
		d = get_smiles_similarity(val1, val2)
	    else:
		raise ValueError("Uknown method: %s" % method)
	    drug_to_drug_similarity.setdefault(drug1, {})[drug2] = d
	    drug_to_drug_similarity.setdefault(drug2, {})[drug1] = d
    if dump_file is not None:
	cPickle.dump(drug_to_drug_similarity, open(dump_file,'w'))
    return drug_to_drug_similarity


###### Chemical similarity based target prediction (~SEA) related ######

def get_chemical_similarity_based_target_predictions(parameters, smiles_list, cutoff=0.9, method="fishers"):
    """
    method: any_smiles / at_least_one_above (inherit targets of any matching smiles w.r.t. cutoff) | majority_above (inherit the target majority of whose smiles is above cutoff) | all_above (inherit the target all of whose smiles are above cutoff)
    """
    parser = get_drugbank(parameters)
    smiles_to_geneids, geneid_to_smiles_strings = get_drug_smiles_by_target(parameters, parser)
    #print len(geneid_to_smiles_strings), geneid_to_smiles_strings.items()[:5]
    all_smiles = reduce(lambda x,y: x|y, geneid_to_smiles_strings.values())
    smiles_to_smiles_similarity = {}
    for smiles1 in smiles_list:
	smiles_to_smiles_similarity[smiles1] = {}
	for smiles2 in all_smiles:
	    try:
		d = get_smiles_similarity(smiles1, smiles2, fp_type="sim", metric="tanimoto") 
	    except:
		#print smiles1, smiles2 # chirality not possible
		continue
	    smiles_to_smiles_similarity[smiles1][smiles2] = d 
    smiles_to_geneids_predicted = {}
    if method == "fishers":
	smiles_to_query_smiles = {}
	for smiles in all_smiles:
	    for smiles_query, smiles_to_value in smiles_to_smiles_similarity.iteritems():
		try:
		    d = smiles_to_value[smiles]
		except:
		    continue
		if d >= cutoff:
		    smiles_to_query_smiles.setdefault(smiles, set()).add(smiles_query)
	print "Drugs with matching smiles:", len(smiles_to_query_smiles)
	smiles_to_geneids_predicted = get_side_effect_targets_fishers(smiles_to_geneids, smiles_to_query_smiles, cutoff=float(parameters.get("fdr_cutoff")), correct_pvalues=True)
    elif method == "any_smiles" or method == "at_least_one_above":
	for smiles in smiles_list:
	    for smiles2, d in smiles_to_smiles_similarity[smiles].iteritems():
		if d >= cutoff:
		    print smiles2, d
		    geneids = smiles_to_geneids_predicted.setdefault(smiles, set())
		    geneids |= smiles_to_geneids[smiles2]
    elif method in ("majority_above", "all_above"):
	for smiles in smiles_list:
	    for geneid, smiles_strings in geneid_to_smiles_strings.iteritems():
		values = []
		for smiles2 in smiles_strings:
		    try:
			#print smiles2 # chirality not possible
			d = smiles_to_smiles_similarity[smiles][smiles2]
			values.append(d >= cutoff)
		    except:
			continue
		if len(values) == 0:
		    continue
		if method == "majority_above":
		    if sum(values) / float(len(values)) > 0.5: #any(values):
		    
			smiles_to_geneids_predicted.setdefault(smiles, set()).add(geneid)
		elif method == "all_above":
		    if all(values):
			smiles_to_geneids_predicted.setdefault(smiles, set()).add(geneid)
	    # If smiles is among known smiles, add known targets
	    if smiles in smiles_to_geneids:
		smiles_to_geneids_predicted[smiles] |= smiles_to_geneids[smiles]
    else:
	raise ValueError("Unknown method: %s!" % method)
    return smiles_to_geneids_predicted


###### Side effect targets related ######

def get_side_effect_targets(parameters, source = "sider"):
    if source.startswith("sider"):
	dump_file = parameters.get("sider_dir") + "/side_effect_to_targets.pcl" 
	if source != "sider":
	    dump_file += "." + source
    elif source == "offsides":
	dump_file = parameters.get("offsides_dir") + "/side_effect_to_targets.pcl" 
    else:
	raise ValueError("Uknown source: %s" % source)
    if os.path.exists(dump_file):
	side_effect_to_targets = cPickle.load(open(dump_file))
	return side_effect_to_targets
    # Get drugs and their targets
    drug_to_geneids = get_drug_info(parameters)
    # Get side effect info
    if source.startswith("sider"):
	drug_to_side_effects = get_drug_side_effects(parameters, source) 
    elif source == "offsides":
	drug_to_side_effects = get_offsides(parameters)
    #print len(drug_to_side_effects), drug_to_side_effects.items()[:5]
    # Side effect protein target sets w.r.t. FDR <=0.2
    side_effect_to_targets = get_side_effect_targets_fishers(drug_to_geneids, drug_to_side_effects, cutoff=float(parameters.get("fdr_cutoff")), correct_pvalues=True)
    cPickle.dump(side_effect_to_targets, open(dump_file,'w'))
    return side_effect_to_targets


def get_drug_side_effect_subset(parameters, drug_to_side_effects):
    n_fold = int(parameters.get("n_fold"))
    random.seed(int(parameters.get("random_seed")))
    pairs = [ (drug, side_effect) for drug, side_effects in drug_to_side_effects.iteritems() for side_effect in side_effects ]
    random.shuffle(pairs)
    values = []
    n = len(pairs) / n_fold 
    #n = len(pairs) / 10 # to check 1-fold of a possible 10-fold
    for i in xrange(n_fold):
    #for i, j in [(0, 9*n), (9*n, len(pairs))]:
	drug_to_side_effects_sub = {}
	for drug, side_effect in pairs[i*n:(i+1)*n]:
	#for drug, side_effect in pairs[i:j]:
	    drug_to_side_effects_sub.setdefault(drug, set()).add(side_effect)
	values.append(drug_to_side_effects_sub)
    return values


def get_side_effect_target_symbols(parameters, source = "sider", output_file=None):
    # Get gene id - name mapping
    geneid_to_names, name_to_geneid = wrappers.get_geneid_symbol_mapping(parameters.get("id_mapping_file"))
    # Get side effect targets
    side_effect_to_targets = get_side_effect_targets(parameters, source)
    # Create side effect target mapping file
    if output_file is not None:
	f = open(file_name, 'w') 
    side_effect_to_genes = {}
    for side_effect, targets in side_effect_to_targets.iteritems(): 
	values = [] 
	for target in targets:
	    if target in geneid_to_names:
		values.extend(list(geneid_to_names[target]))
	values.sort()
	side_effect_to_genes[side_effect] = set(values)
	#print side_effect, len(targets), len(values)
	#values = targets
	if output_file is not None:
	    f.write("\t%s\t%s\n" % (side_effect, "\t".join(values)))
    if output_file is not None:
	f.close()
    return side_effect_to_genes


def get_side_effect_targets_fishers(drug_to_geneids, drug_to_side_effects, cutoff=0.2, correct_pvalues=True): # min_n_drug=5, 
    """
    cutoff: p-value or fdr cutoff (0.2)
    correct_pvalues: apply multiple hypothesis testing (True)
    (obselete) min_n_drug: Consider only side effects that are associated with at least n drugs (5) 
    """
    # Get side effect to drugs
    drugs_all = set()
    side_effect_to_drugs = {}
    for drug, side_effects in drug_to_side_effects.iteritems():
	if drug not in drug_to_geneids:
	    continue
	for side_effect in side_effects:
	    side_effect_to_drugs.setdefault(side_effect, set()).add(drug)
	drugs_all.add(drug)
    #print len(side_effect_to_drugs), side_effect_to_drugs.items()[:5]
    # Get target to drugs
    target_to_drugs = {}
    for drug, geneids in drug_to_geneids.iteritems():
	if drug not in drug_to_side_effects:
	    continue
	for geneid in geneids:
	    target_to_drugs.setdefault(geneid, set()).add(drug)
	if drug not in drugs_all:
	    print "Side effect info but no target info:", drug
    #print len(target_to_drugs), target_to_drugs.items()[:5]
    # Get side effect to targets
    side_effect_to_targets = {}
    n_less_than_five = 0
    for side_effect, drugs_se in side_effect_to_drugs.iteritems():
	#if len(drugs_se) < min_n_drug: 
	#    n_less_than_five += 1
	#    continue
	values = []
	for target, drugs_target in target_to_drugs.iteritems():
	    tp = len(drugs_se & drugs_target)
	    fp = len(drugs_target) - tp
	    fn = len(drugs_se) - tp
	    tn = len(drugs_all) - (tp + fp + fn)
	    oddsratio, pvalue = stat_utilities.fisher_exact(tp, fp, fn, tn, alternative="greater")
	    #if target == "19" and side_effect == "nausea": 
	    #	print side_effect, oddsratio, pvalue
	    if correct_pvalues:
		values.append((pvalue, target))
	    else:
		if pvalue <= cutoff:
		    side_effect_to_targets.setdefault(side_effect, set()).add(target)
	if correct_pvalues:
	    pvalues_new = stat_utilities.correct_pvalues_for_multiple_testing(zip(*values)[0])
	    for i, pvalue in enumerate(pvalues_new):
		if pvalue <= cutoff:
		    side_effect_to_targets.setdefault(side_effect, set()).add(values[i][1])
    #print len(side_effect_to_drugs), n_less_than_five, len(side_effect_to_targets)
    return side_effect_to_targets


###### Sider related ######

def get_sider(parameters):
    indication_file = parameters.get("sider_dir") + "/meddra_all_label_indications.tsv.gz" 
    side_effect_file = parameters.get("sider_dir") + "/meddra_all_label_se.tsv.gz" 
    dump_file = parameters.get("sider_dir") + "/sider.pcl"
    pubchem_to_indications, pubchem_to_side_effects = parse_sider_v4.get_sider_info(indication_file, side_effect_file, dump_file)
    return pubchem_to_indications, pubchem_to_side_effects 


def get_drug_side_effects(parameters, source=None):
    dump_file = parameters.get("sider_dir") + "/drug_side_effects.pcl" 
    if os.path.exists(dump_file):
	drugbank_id_to_side_effects = cPickle.load(open(dump_file))
	if source is not None:
	    if source == "sider1":
		drugbank_id_to_side_effects = get_drug_side_effect_subset(parameters, drugbank_id_to_side_effects)[0]
	    elif source == "sider2":
		drugbank_id_to_side_effects = get_drug_side_effect_subset(parameters, drugbank_id_to_side_effects)[1]
	return drugbank_id_to_side_effects 
    # Get sider info
    pubchem_to_indications, pubchem_to_side_effects = get_sider(parameters)
    #print "SIDER"
    #print len(pubchem_to_indications), pubchem_to_indications.items()[:5]
    #print len(pubchem_to_side_effects), pubchem_to_side_effects.items()[:5]
    # Get stitch drugbank mapping
    pubchem_to_drugbank_ids, pubchem_to_target_to_score = get_stitch(parameters)
    #print "STITCH"
    #print len(pubchem_to_drugbank_ids), pubchem_to_drugbank_ids.items()[:5]
    #print len(pubchem_to_target_to_score), pubchem_to_target_to_score.items()[:5]
    # Get drugbank pubchem mapping
    parser = get_drugbank(parameters)
    pubchem_to_drugbank_ids_db = {}
    for drug_to_values in (parser.drug_to_pubchem, parser.drug_to_pubchem_substance):
	for drugbank_id, pubchem in drug_to_values.iteritems():
	    pubchem_to_drugbank_ids_db.setdefault(pubchem, set()).add(drugbank_id)
    #print "DB"
    #print len(pubchem_to_drugbank_ids_db), pubchem_to_drugbank_ids_db.items()[:5]
    drugbank_id_to_side_effects = {}
    for pubchem, side_effects in pubchem_to_side_effects.iteritems(): 
	if pubchem not in pubchem_to_side_effects:
	    continue
	db_ids1 = set()
	db_ids2 = set()
	if pubchem in pubchem_to_drugbank_ids:
	    db_ids1 = pubchem_to_drugbank_ids[pubchem]
	    #drugbank_ids |= db_ids1
	if pubchem in pubchem_to_drugbank_ids_db:
	    db_ids2 = pubchem_to_drugbank_ids_db[pubchem]
	    #drugbank_ids |= db_ids2
	drugbank_ids = db_ids2
	if len(db_ids2) == 0: # No mapping in drugbank
	    if len(db_ids1) == 1: # accept pubchem mapping
		drugbank_ids = db_ids1
	    else: # ignore (the right pubchem will come from the one from drugbank)
		continue 
	for drugbank_id in drugbank_ids:
	    for side_effect in pubchem_to_side_effects[pubchem]:
		drugbank_id_to_side_effects.setdefault(drugbank_id, set()).add(side_effect)
	#if len(db_ids1 & db_ids2) == 0:
	#    print pubchem, len(db_ids1 & db_ids2), len(db_ids1 | db_ids2), db_ids1, db_ids2
    cPickle.dump(drugbank_id_to_side_effects, open(dump_file,'w'))
    return drugbank_id_to_side_effects


##### STITCH related #####

def get_stitch(parameters):
    base_dir = parameters.get("stitch_dir") + "/" 
    chemicals_file = base_dir + "chemicals.inchikeys.v4.0.1.tsv.gz"
    alias_file = base_dir + "chemical.aliases.v4.0.tsv.gz"
    inchikey_file = base_dir + "chemicals.inchikeys.v4.0.1.tsv.gz"
    links_file = base_dir + "9606.protein_chemical.links.detailed.v4.0.tsv.gz"
    gene_mapping_file = parameters.get("string_dir") + "/entrez_gene_id.vs.string.v10.28042015.tsv"
    dump_file = base_dir + "stitch.pcl"
    cid_to_drugbank_ids, cid_to_target_to_score = parse_stitch.get_stitch_info(chemicals_file, alias_file, inchikey_file, links_file, gene_mapping_file, dump_file, species_prefix = "9606")
    return cid_to_drugbank_ids, cid_to_target_to_score


##### Other side effect resoureces (OFFSIDES / NUGENT) related #####

def get_offsides(parameters):
    # Get offsides data
    parser = TsvReader.TsvReader(parameters.get("offsides_file"), delim="\t", inner_delim = None, quotation='"')
    header_to_idx, cid_to_values = parser.read(fields_to_include = ["stitch_id", "event", "pvalue", "bg_correction", "sider", "future_aers", "medeffect"], keys_to_include = None, merge_inner_values = False)
    #print len(cid_to_values) #, cid_to_values.items()[:3]
    # Get stitch drugbank mapping
    pubchem_to_drugbank_ids, pubchem_to_target_to_score = get_stitch(parameters)
    # Get drugbank side effect mapping
    drug_to_side_effects = {}
    not_in_db = set()
    side_effects = set()
    for cid, values in cid_to_values.iteritems():
	side_effects |= set(map(lambda x: x.lower(), zip(*values)[0]))
    n_side_effects = len(side_effects)
    #print "Number of side effects:", n_side_effects
    for cid, values in cid_to_values.iteritems():
	cid = cid[3:]
	if cid.startswith("1"):
	    cid = "%s" % (abs(int(cid)) - 100000000)
	else:
	    cid = "%s" % abs(int(cid)) 
	if cid not in pubchem_to_drugbank_ids:
	    not_in_db.add(cid)
	    continue
	drugs = pubchem_to_drugbank_ids[cid]
	for side_effect, pval, bg_correction, sider, future_aers, medeffect in values:
	    #print side_effect, pval, bg_correction
	    # Convert all lower / upper case starting side effect words to Sider format
	    #side_effect = " ".join(map(lambda x: x[0].lower() + x[1:], side_effect.split(" ")))
	    side_effect = side_effect.lower()
	    side_effect = side_effect[0].upper() + side_effect[1:] 
	    #if n_side_effects*float(pval) <= 0.2: # and float(bg_correction) <= 0.001:
	    #if future_aers == "1":
	    if medeffect == "1": # and sider != "1": 
		for drug in drugs:
		    drug_to_side_effects.setdefault(drug, set()).add(side_effect)
    #print len(drug_to_side_effects), drug_to_side_effects.items()[:3]
    #print "Drugs not in drugbank mapping:", len(not_in_db) #, not_in_db
    return drug_to_side_effects 


def get_nugent(parameters):
    # Get nugent data
    f = open(parameters.get("nugent_file"))
    f.readline()
    name_to_side_effect_counts = {}
    for line in f:
	tweet_id, drugs, side_effects = line.strip("\n").split("\t")
	if side_effects == "":
	    continue
	for drug in drugs.split("|"):
	    d = name_to_side_effect_counts.setdefault(drug, {})
	    #d |= set(map(lambda x: x[0].upper() + x[1:], side_effects.split("|")))
	    for side_effect in side_effects.split("|"):
		side_effect = side_effect[0].upper() + side_effect[1:]
		e = d.setdefault(side_effect, 0)
		d[side_effect] = e + 1
    # Get drugbank side effect mapping
    parser = get_drugbank(parameters)
    name_to_drug, synonym_to_drug = parser.get_synonyms(selected_drugs=None, only_synonyms=False)
    drug_to_side_effects = {}
    not_in_db = set()
    for name, side_effect_to_count in name_to_side_effect_counts.iteritems():
	# Get drugbank name mapping
	if name in name_to_drug:
	    drug = name_to_drug[name]
	elif name in synonym_to_drug:
	    drug = synonym_to_drug[name]
	else:
	    not_in_db.add(name)
	    continue
	count_total = 0.0
	for side_effect, count in side_effect_to_count.iteritems():
	    count_total += count
	for side_effect, count in side_effect_to_count.iteritems():
	    #if count / count_total > 0.5: 
	    if count >= 5:
		drug_to_side_effects.setdefault(drug, set()).add(side_effect)
    print len(drug_to_side_effects), drug_to_side_effects.items()[:3]
    print "Drugs not in drugbank mapping:", len(not_in_db) #, not_in_db
    return drug_to_side_effects 



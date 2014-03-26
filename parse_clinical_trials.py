
def main():
    return

def get_ctid_to_conditions(base_dir):
    condition_file = base_dir + "/conditions.txt"
    condition_file2 = base_dir + "/condition_browse.txt"
    # Get conditions
    ctid_to_conditions = {} 
    f = open(condition_file)
    f.readline()
    for line in f:
	words = line.strip().split("|")
	ctid = words[1]
	condition = words[2].lower()
	ctid_to_conditions.setdefault(ctid, set()).add(condition)
    f.close() 
    return ctid_to_conditions
    f = open(condition_file2)
    f.readline()
    for line in f:
	words = line.strip().split("|")
	ctid = words[1]
	condition = words[2].lower()
	ctid_to_conditions.setdefault(ctid, set()).add(condition)
    f.close() 
    return ctid_to_conditions


def get_ctid_to_phase(base_dir):
    study_file = base_dir + "/clinical_study.txt"
    # Get phase information
    ctid_to_phase = {}
    f = open(study_file)
    line = f.readline()
    words = line.strip().split("|")
    header_to_idx = dict((word.lower(), i) for i, word in enumerate(words))
    ctid_to_phase = {}
    text = None
    while line:
	line = f.readline()
	if line.startswith("NCT"):
	    if text is not None:
		words = text.strip().split("|")
		ctid = words[0]
		try:
		    phase = words[header_to_idx["phase"]]
		except:
		    print words
		    return
		if phase.strip() != "":
		    ctid_to_phase[ctid] = phase
	    text = line
	else:
	    text += line
    f.close() 
    return ctid_to_phase


def get_interventions(base_dir): 
    #ctid_to_drugs = {}
    drug_to_ctids = {}
    intervention_file = base_dir + "/interventions.txt"
    f = open(intervention_file)
    f.readline()
    prev_row = 0
    for line in f:
	words = line.strip().split("|")
	try:
	    row = int(words[0])
	    if row != prev_row + 1:
		continue
	except:
	    continue
	prev_row += 1
	try:
	    ctid = words[1]
	except:
	    print words
	intervention = words[2]
	drug = words[3].lower()
	drug = drug.decode("ascii", errors="ignore").encode("ascii")
	drug = drug.strip("\"'")
	if intervention != "Drug":
	    continue
	drug_to_ctids.setdefault(drug, set()).add(ctid)
	#ctid_to_drugs.setdefault(ctid, set()).add(drug)
	#conditions = drug_to_interventions.setdefault(drug, set())
	#conditions |= ctid_to_conditions[ctid]
    f.close()
    intervention_file = base_dir + "/intervention_browse.txt"
    f = open(intervention_file)
    f.readline()
    for line in f:
	words = line.strip().split("|")
	row = int(words[0])
	ctid = words[1]
	drug = words[2].lower()
	drug = drug.decode("ascii", errors="ignore").encode("ascii")
	drug = drug.strip("\"'")
	drug_to_ctids.setdefault(drug, set()).add(ctid)
	#ctid_to_drugs.setdefault(ctid, set()).add(drug)
    f.close()
    intervention_file = base_dir + "/intervention_other_names.txt"
    f = open(intervention_file)
    f.readline()
    for line in f:
	words = line.strip().split("|")
	row = int(words[0])
	ctid = words[1]
	drug = words[3].lower()
	drug = drug.decode("ascii", errors="ignore").encode("ascii")
	drug = drug.strip("\"'")
	drug_to_ctids.setdefault(drug, set()).add(ctid)
	#ctid_to_drugs.setdefault(ctid, set()).add(drug)
    f.close()
    #print set(ctid_to_phase.values())
    return drug_to_ctids #ctid_to_drugs 


def get_drug_to_interventions(base_dir, name_to_drug, synonym_to_drug):
    drug_to_ctids = get_interventions(base_dir)
    ctid_to_conditions = get_ctid_to_conditions(base_dir)
    #ctid_to_phase = get_ctid_to_phase(base_dir) 
    # Get interventions
    drug_to_interventions = {}
    for drug, ctids in drug_to_ctids.iteritems():
	drugbank_id = None
	if name_to_drug is None:
	    drugbank_id = drug
	else:
	    if drug in name_to_drug:
		drugbank_id = name_to_drug[drug]
	    elif drug in synonym_to_drug:
		drugbank_id = synonym_to_drug[drug]
	    else:
		continue
	values = set()
	for ctid in ctids:
	    #if ctid_to_phase[ctid] != "Phase 3":  
	    #	continue
	    values |= ctid_to_conditions[ctid]
	if len(values) == 0:
	    continue
	drug_to_interventions.setdefault(drugbank_id, values)
    #print drug_to_interventions["voriconazole"]
    #for ctid in drug_to_ctids["voriconazole"]:
    #	print ctid, ctid_to_conditions[ctid], ctid_to_phase[ctid]
    return drug_to_interventions 


def get_frequent_interventions(drug_to_interventions):
    condition_to_count = {}
    for drug, interventions in drug_to_interventions.iteritems():
	for condition in interventions:
	    if condition in condition_to_count:
		condition_to_count[condition] += 1
	    else:
		condition_to_count[condition] = 1
    values = []
    for condition, count in condition_to_count.iteritems():
	values.append((count, condition))
    values.sort()
    values.reverse()
    #print values[:50]
    return values

if __name__ == "__main__":
    main()


import TsvReader
import sql_utilities

def main():
    file_name = "/home/emre/data/target/chembl/chembl_drugtargets-18_17_33_12.txt"
    drug_to_target_id = get_molecule_to_target_id(file_name)
    print drug_to_target_id["CHEMBL502"]
    file_name = "/home/emre/data/target/chembl/target-18_15_28_11.txt"
    target_to_uniprots = get_target_id_to_uniprots(file_name, organism="Homo sapiens")
    print target_to_uniprots["CHEMBL220"]
    print target_to_uniprots["CHEMBL2095233"]
    #drugs = ["CHEMBL502", "CHEMBL1678", "CHEMBL2171854"]
    drugs = drug_to_target_id.keys()
    drug_to_targets = retrieve_targets_from_database(drugs)
    print len(drug_to_targets)
    file_name = "/home/emre/data/target/chembl/targets_retrieved.txt"
    f = open(file_name, 'w')
    f.write("ChemblId\tUniprotId\tType\n")
    for drug, targets in drug_to_targets.iteritems():
	if drug in drug_to_target_id:
	    target_id = drug_to_target_id[drug]
	    if target_id in target_to_uniprots:
		for target in target_to_uniprots[target_id]:
		    f.write("%s\t%s\t%s\n" % (drug, target, "Target"))
	for target in targets:
	    f.write("%s\t%s\t%s\n" % (drug, target, "Activity"))
    f.close()
    return


def get_molecule_to_target_id(file_name):
    parser = TsvReader.TsvReader(file_name, delim="\t", inner_delim=None)
    column_to_index, id_to_values = parser.read(fields_to_include=["MOLECULE_CHEMBL_ID", "TARGET_CHEMBL_ID"], keys_to_include=None, merge_inner_values=False)
    # MOLECULE_NAME, CANONICAL_SMILES, ACTION_TYPE
    # name  = values[0][column_to_index["molecule_name"]])
    drug_to_target = {}
    for drug, values in id_to_values.iteritems():
	if len(values[0]) > 1:
	    print values
	drug_to_target[drug] = values[0][0]
    return drug_to_target 


def get_target_id_to_uniprots(file_name, organism="Homo sapiens"):
    parser = TsvReader.TsvReader(file_name, delim="\t", inner_delim=None)
    column_to_index, id_to_values = parser.read(fields_to_include=["CHEMBL_ID", "PROTEIN_ACCESSION", "ORGANISM"], keys_to_include=None, merge_inner_values=False)
    #TID     PREF_NAME       
    target_to_uniprots = {}
    for target, values in id_to_values.iteritems():
	if len(values[0]) > 2:
	    print values
	if values[0][column_to_index["organism"]] != organism:
	    continue
	target_to_uniprots[target] = set(values[0][column_to_index["protein_accession"]].split(", "))
    return target_to_uniprots 


def retrieve_targets_from_database(drugs):
    drug_to_targets = {}
    dbhost = None #"127.0.0.1"
    dbuser = "emre" 
    dbname = "chembl_23"
    cursor, db = sql_utilities.get_db_cursor(dbhost, dbname, dbuser, dbpass=None)
    tables = ["docs D", "molecule_dictionary M", "activities A", "assays S", "target_dictionary T", "target_components C", "component_sequences CC"]
    columns = ["M.chembl_id", "CC.accession"] 
    fixed_conditions = [["T.tax_id", "=", "9606"]]
    join_conditions = [["M.molregno", "=", "A.molregno"], ["D.doc_id", "=", "A.doc_id"], ["A.assay_id", "=", "S.assay_id"], ["S.tid", "=", "T.tid"], ["C.tid", "=", "T.tid"], ["C.component_id", "=", "CC.component_id"], ["M.chembl_id", "IN", "(\'%s\')" % "\',\'".join(drugs)]] 
    sql = sql_utilities.get_select_sql_query(tables, columns = columns, fixed_conditions = fixed_conditions, join_conditions = join_conditions, group_conditions=None, distinct_columns=False)
    #print sql
    results = sql_utilities.get_values_from_database(cursor, sql)
    if results:
	for i, row in enumerate(results):
	    if i < 5:
		print row
	    drug_to_targets.setdefault(row[0], set()).add(row[1])
    sql_utilities.close_db_cursor(cursor, db)
    return drug_to_targets

   

if __name__ == "__main__":
    main()


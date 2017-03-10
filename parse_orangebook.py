#import os, cPickle
from datetime import datetime

def main():
    name = "mesna"
    base_dir = "/home/eguney/data/orangebook/"
    name_to_ids, id_to_exclusivity, id_to_patent = get_orangebook_data(base_dir)
    dates = sorted(id_to_exclusivity["001"])
    print name_to_ids[name], dates[0], dates[-1] 
    dates = sorted(id_to_patent["001"])
    print name_to_ids[name], dates[0], dates[-1] 
    return 


def get_orangebook_data(base_dir):
    mapping_file = base_dir + "products.txt"
    name_to_ids = get_product_mapping(mapping_file)
    mapping_file = base_dir + "exclusivity.txt"
    id_to_values_exc = get_exclusivity_mapping(mapping_file)
    mapping_file = base_dir + "patent.txt"
    id_to_values_pat = get_patent_mapping(mapping_file)
    return name_to_ids, id_to_values_exc, id_to_values_pat 


def get_product_mapping(mapping_file, delim = "~"):
    """ 
    Ingredient~DF;Route~Trade_Name~Applicant~Strength~Appl_Type~Appl_No~Product_No~TE_Code~Approval_Date~RLD~Type~Applicant_Full_Name
    BUDESONIDE~AEROSOL, FOAM;RECTAL~UCERIS~VALEANT PHARMS INTL~2MG/ACTUATION~N~205613~001~~Oct 7, 2014~Yes~RX~VALEANT PHARMACEUTICALS INTERNATIONAL
    """
    name_to_ids = {}
    f = open(mapping_file)
    header = f.readline()
    for line in f:
        words = line.strip().split(delim)
        name = words[2].lower()
        product_id = words[7]
	date = words[9]
	#date = datetime.strptime(date, "%b %d, %Y") # "Approved prior to"
	#applicant = words[10]
	#if name in name_to_id:
	#    print "Overwriting", name, name_to_id[name], product_id
	#name_to_id[name] = product_id
	name_to_ids.setdefault(name, set()).add(product_id)
    f.close()
    return name_to_ids
    
	
def get_exclusivity_mapping(mapping_file, delim = "~"):
    """ 
    Appl_Type~Appl_No~Product_No~Exclusivity_Code~Exclusivity_Date
    N~008372~008~ODE~Oct 15, 2017
    """
    id_to_values = {}
    f = open(mapping_file)
    header = f.readline()
    for line in f:
        words = line.strip().split(delim)
        product_id = words[2]
        code = words[3]
	date = words[4]
	date = datetime.strptime(date, "%b %d, %Y")
	#if product_id in id_to_values:
	#    print "Overwriting", product_id, id_to_values[product_id]
        id_to_values.setdefault(product_id, set()).add((date, code))
    f.close()
    return id_to_values
    

def get_patent_mapping(mapping_file, delim = "~"):
    """ 
    Appl_Type~Appl_No~Product_No~Patent_No~Patent_Expire_Date_Text~Drug_Substance_Flag~Drug_Product_Flag~Patent_Use_Code~Delist_Flag
    N~013718~002~5872147~Dec 5, 2017~~~U-585~
    """
    id_to_values = {}
    f = open(mapping_file)
    header = f.readline()
    for line in f:
        words = line.strip().split(delim)
        product_id = words[2]
        code = words[3]
	date = words[4]
	date = datetime.strptime(date, "%b %d, %Y")
	#if product_id in id_to_values:
	#    print "Overwriting", product_id, id_to_values[product_id]
        id_to_values.setdefault(product_id, set()).add((date, code))
    f.close()
    return id_to_values
    

if __name__ == "__main__":
    main()


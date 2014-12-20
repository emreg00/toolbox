import parse_ncbi, parse_drugbank


def main():
    base_dir = "/home/emre/arastirma/data/drug/sensitivity/gdsc/"
    #file_name = base_dir + "gdsc_compounds_conc_w5.csv"
    #compound_to_concentrations = get_compounds(file_name)
    #print len(compound_to_concentrations), compound_to_concentrations["GNF-2"]
    geneid_to_names, name_to_geneid = parse_ncbi.get_geneid_symbol_mapping(base_dir + "../../../proteome/ncbi/geneid_to_symbol.txt")
    file_target = base_dir + "gdsc_en_output_w5.csv"
    file_response = base_dir + "gdsc_manova_output_w5.csv"
    get_gsdc_info(file_target, file_response, name_to_geneid) 
    return


def get_gsdc_info(file_target, file_response, name_to_geneid, name_to_drug=None, drug_to_geneids=None):
    log_info = False
    # Parse files
    compound_to_targets = get_targets(file_target)
    #print len(compound_to_targets), compound_to_targets["GNF-2"]
    compound_to_gene_to_values = get_drug_response(file_response)
    #print len(compound_to_gene_to_values), len(compound_to_gene_to_values["GNF-2"]), compound_to_gene_to_values["GNF-2"]["FGFR3"]
    # Map genes to geneids 
    ##name_to_geneid["PI3KB"] = "5287"
    ##name_to_geneid["SRC"] = "6714"
    name_to_geneid["PDGFR"] = "5159" # PDGFRB
    name_to_geneid["ABL"] = "25" # ABL1
    #name_to_geneid["ATM"] = "472"
    name_to_geneid["FAK"] = "5747" # PTK2
    name_to_geneid["MTORC1"] = "84335" # AKT1S1
    name_to_geneid["MTORC2"] = "253260" # RICTOR
    name_to_geneid["BCLW"] = "599" # BCL2L2
    name_to_geneid["BCLXL"] = "598" # BCL2L1
    name_to_geneid["MEK2"] = "5605" # MAP2K2
    name_to_geneid["MEK1"] = "5604" # MAP2K1
    name_to_geneid["RAF"] = "5894" # RAF1
    name_to_geneid["HER2"] = "2064" # ERBB2
    name_to_geneid["CHK1"] = "1111" # CHEK1
    name_to_geneid["CHK2"] = "11200" # CHEK2
    name_to_geneid["FAM123B"] = name_to_geneid["AMER1"]
    name_to_geneid["MYCL1"] = name_to_geneid["MYCL"]
    compound_to_geneids = {}
    compounds_not_found = {}
    for compound, targets in compound_to_targets.iteritems():
        compound = compound.lower()
        targets = targets.replace(" AND ", " : ")
        if targets in ("SRC FAMILY ", "SRC-FAMILY"):
            targets = "SRC"
        elif targets == "P53-MDM2 INTERACTION":
            targets = "TP53 : MDM2"
        elif targets == "BCL2 FAMILY":
            targets = "BCL2"
        elif targets == "SERUM/GLUCOCORTICOID REGULATED KINASE 1 ":
            targets = "SGK1"
        elif targets == "CASPASE 3 ACTIVATOR":
            targets = "CASP3"
        elif targets == "TBK1 :  PDK1 :  IKK ":
            targets = "TBK1 : PDK1 : IKBKB"
        elif targets == "AURORA B":
            targets = "AURKB"
        elif targets == "AMPK AGONIST":
            targets = "PRKAA1 : PRKAA2"
        elif targets == "SPLEEN TYROSINE KINASE":
            targets = "SYK"
        elif targets == "ATM (IC50 13 NM) (ATR >>10 MM)":
            targets = "ATM : ATR"
        elif targets == "PRKC":
            targets = "PRKCA : PRKCB :  PRKCD : PRKCE : PRKCG : PRKCI"
        elif targets == "PROLYL4HYDROXYLASE.":
            targets = "P4HA1"
        elif targets == "DIHYDROFOLATE REDUCTASE (DHFR)":
            targets = "DHFR"
        elif targets in ("BCR-ABL", "BCRABL ONLY (ALLOSTERIC NONATP COMPETETIVE)"):
            targets = "BCR : ABL1"
        elif targets == "G-SECRETASE":
            targets = "APP : PSENEN : APH1A : APH1B"
        else:
            idx = targets.find("/")
            if idx != -1:
                prefix = targets[:idx-1]
                targets = prefix + targets[idx-1] + " : " + prefix + targets[idx+1]
                targets = targets.replace(" ", "")
        targets = targets.split(":")
        geneids = set()
        for target in targets:
            target = target.strip().replace("-", "")
            if target in name_to_geneid:
                geneids.add(name_to_geneid[target])
            else:
                compounds_not_found.setdefault(compound, set()).add(target)
        if compound == "bibw2992":
            geneids.add("2066")
        if compound in name_to_drug:
            drug = name_to_drug[compound]
            if drug in drug_to_geneids:
                geneids |= drug_to_geneids[drug]
        if len(geneids) == 0:
            continue
        compound_to_geneids[compound] = geneids
    compound_to_gene_to_values_mod = {}
    gene_to_geneid = {}
    genes_not_found = set()
    for compound, gene_to_values in compound_to_gene_to_values.iteritems():
        compound = compound.lower()
        for gene, values in gene_to_values.iteritems():
            if gene not in name_to_geneid:
                genes_not_found.add(gene)
            else:
                compound_to_gene_to_values_mod.setdefault(compound, {})[gene] = values
                gene_to_geneid[gene] = name_to_geneid[gene]
    compound_to_gene_to_values = compound_to_gene_to_values_mod 
    if log_info:
        print "Not found:"
        for k, v in compounds_not_found.iteritems():
            print k, v
        print len(compound_to_geneids), compound_to_geneids["erlotinib"]
        print genes_not_found
        print "Identical targets:"
        for i, compound in enumerate(compound_to_geneids):
            geneids = compound_to_geneids[compound]
            for j, compound2 in enumerate(compound_to_geneids):
                if i < j:
                    geneids2 = compound_to_geneids[compound2]
                    if geneids == geneids2: #len(geneids & geneids2) > 0:
                        print compound, compound2, geneids & geneids2
    return compound_to_geneids, compound_to_gene_to_values, gene_to_geneid


def get_compounds(file_name):
    compound_to_concentrations = {}
    f = open(file_name)
    line = f.readline()
    for line in f:
	words = line.strip().split(",")
        compound, c_min, c_max = words
        compound_to_concentrations[compound] = (float(c_min), float(c_max))
    f.close()
    return compound_to_concentrations


def get_targets(file_name):
    compound_to_targets = {}
    f = open(file_name)
    line = f.readline()
    for line in f:
	words = line.strip().split(",")
        compound =  words[2]
        if compound in compound_to_targets:
            continue
        #print compound
        targets = words[3].upper()
        compound_to_targets[compound] = targets 
    f.close()
    return compound_to_targets


def get_drug_response(file_name):
    compound_to_gene_to_values = {}
    f = open(file_name)
    line = f.readline()
    for line in f:
	words = line.strip().split(",")
        compound = words[0]
        gene = words[1]
        p_value = float(words[6])
        effect = float(words[7])
        compound_to_gene_to_values.setdefault(compound, {})[gene] = (p_value, effect)
    f.close()
    return compound_to_gene_to_values


if __name__ == "__main__":
    main()



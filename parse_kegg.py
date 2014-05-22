import urllib
import urllib2
from bs4 import BeautifulSoup
from diseasome.src import diseasome

# Taken from Durek and Walter http://www.biomedcentral.com/1752-0509/2/100
CURRENCY_METABOLITES = set(["C00001", "C00002", "C00003", "C00004", "C00005", "C00006", "C00007", "C00008", "C00009", "C00010", "C00011", "C00013", "C00014", "C00015", "C00016", "C00018", "C00019", "C00020", "C00021", "C00023", "C00027", "C00028", "C00030", "C00034", "C00035", "C00038", "C00044", "C00050", "C00055", "C00061", "C00063", "C00070", "C00075", "C00076", "C00080", "C00105", "C00112", "C00113", "C00115", "C00120", "C00125", "C00126", "C00138", "C00139", "C00144", "C00175", "C00194", "C00205", "C00238", "C00291", "C01352"])


def main():
    out_file = "/home/emre/data/interactome/METABOLIC/network.txt"
    output_interactions(out_file)
    return


def output_interactions(out_file):
    reactions = get_reactions()
    reaction_to_values = {}
    for reaction_id in reactions:
        enzymes, lefts, rights = get_reaction_info(reaction_id)
        #print reaction_id, enzymes
        #print lefts, rights
        lefts = set(lefts) - CURRENCY_METABOLITES
        rights = set(rights) - CURRENCY_METABOLITES
        if len(lefts) + len(rights) == 0:
            continue
        geneids = set() 
        for enzyme in enzymes:
            geneids |= set(get_enzyme_info(enzyme))
        #print geneids
        if len(geneids) == 0:
            continue
        reaction_to_values[reaction_id] = (geneids, lefts, rights)
    reaction_ids = reaction_to_values.keys()
    print len(reaction_ids)
    interactions = set()
    for i, reaction_id1 in enumerate(reaction_ids):
        for j, reaction_id2 in enumerate(reaction_ids):
            if i <= j:
                (geneids1, lefts1, rights1) = reaction_to_values[reaction_id1]
                (geneids2, lefts2, rights2) = reaction_to_values[reaction_id1]
                if len(lefts1 & rights2) > 0 or len(lefts2 & rights1) > 0:
                    for geneid1 in geneids1:
                        for geneid2 in geneids2:
                            interactions.add(tuple(sorted((geneid1, geneid2))))
    print len(interactions)
    f = open(out_file, 'w')
    for geneid1, geneid2 in interactions:
        f.write("%s\t%s\n" % (geneid1, geneid2))
    f.close()
    return


def get_data(command, parameter):
    url = 'http://rest.kegg.jp/%s/%s' % (command, parameter) #list/reaction
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    for line in response:
        yield line


def get_reactions():
    reactions = []
    for line in get_data("list", "reaction"):
        reaction_id = line.strip("\n").split()[0]
        reactions.append(reaction_id)
        #if len(reactions) > 200:
        #    break #!
    return reactions


def get_reaction_info(reaction_id):
    lefts = []
    rights = []
    enzymes = []
    for line in get_data("get", reaction_id):
        if line.startswith("EQUATION"):
            words = line.strip("\n").split()
            passed = False
            for word in words:
                if word.startswith("C"):
                    if passed:
                        rights.append(word)
                    else:
                        lefts.append(word)
                elif word == "<=>":
                    passed = True
        if line.startswith("ENZYME"):
            words = line.strip("\n").split()
            for word in words[1:]:
                enzymes.append(word)
    return enzymes, lefts, rights


def get_enzyme_info(ec_id, prefix="HSA:"):
    geneids = []
    try:
        for line in get_data("get", ec_id):
            if line.startswith("GENES") or line.lstrip().startswith(prefix):
                idx = line.find(prefix)
                if idx != -1:
                    words = line[idx+5:].strip().split()
                    for word in words:
                        e_idx = word.rfind("(")
                        geneid = word[:e_idx]
                        geneids.append(geneid)
    except urllib2.HTTPError:
        print "No info for", ec_id
        return geneids
    return geneids


if __name__ == "__main__":
    main()


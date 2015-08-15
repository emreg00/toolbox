import urllib2, os, cPickle
from bs4 import BeautifulSoup
import time
from xml.etree.ElementTree import iterparse, tostring
import parse_uniprot

NS = ""

def main():
    base_dir = "/home/emre/data/disease/orphanet"
    file_name = base_dir + "/gene_associations.xml"
    out_file = base_dir + "/disease_to_genes.tsv"
    disease_to_genes, gene_to_uniprot = get_disease_to_genes(file_name)
    uniprot_file = base_dir + "/../../proteome/uniprot/idmapping.tab"
    uniprot_ids = gene_to_uniprot.values() #reduce(lambda x,y: x+y, gene_to_uniprot.values())
    uniprot_to_geneid = parse_uniprot.get_uniprot_to_geneid_from_idmapping_file(uniprot_file, uniprot_ids)
    output_disease_genes(file_name, out_file, uniprot_to_geneid)
    return


def output_disease_genes(file_name, out_file, uniprot_to_geneid=None):
    #geneid_to_names, name_to_geneid = diseasome.get_geneid_symbol_mapping() 
    disease_to_genes, gene_to_uniprot = get_disease_to_genes(file_name)
    f = open(out_file, 'w') 
    count = 0
    all_geneids = set()
    for disease, genes in disease_to_genes.iteritems():
	#f.write("\t%s\t%s\n" % (disease, "\t".join(genes)))
        if uniprot_to_geneid is None:
            geneids = genes #set([ name_to_geneid[gene] for gene in genes if gene in name_to_geneid ])
        else:
            geneids = set()
            for gene in genes:
                if gene in gene_to_uniprot:
                    uniprot = gene_to_uniprot[gene]
                    if uniprot in uniprot_to_geneid:
                        geneid = uniprot_to_geneid[uniprot]
                        try:
                            geneid = int(geneid)
                        except:
                            print "Problem with:", geneid
                        geneids.add(geneid) 
	if len(geneids) > 5:
	    count += 1
        all_geneids |= geneids
        #print disease, geneids
	f.write("\t%s\t%s\n" % (disease, "\t".join(map(str, geneids))))
    f.close()
    print "N disease > 5", count, "N all geneids", len(all_geneids)
    return


def get_disease_to_genes(file_name):
    disease_to_genes = {}
    gene_to_synonyms = {}
    gene_to_uniprot = {}
    context = iterparse(file_name, ["start", "end"])
    context = iter(context)
    event, root = context.next()
    state_stack = [ root.tag ]
    phenotype, gene, source = None, None, None
    for (event, elem) in context:
        if event == "start":
            state_stack.append(elem.tag)
            if elem.tag == NS+"Name":
                gene = None
            if elem.tag == NS+"Source":
                source = None
        elif event == "end":
            if elem.tag == NS+"Name":
                if state_stack[-2] == NS+"Disorder":
                    phenotype = elem.text.encode("ascii", "replace")
            elif elem.tag == NS+"Symbol":
                if state_stack[-4] == NS+"DisorderGeneAssociationList" and state_stack[-2] == "Gene": # This is always the case
                    gene = elem.text
                    disease_to_genes.setdefault(phenotype, set()).add(gene)
            elif elem.tag == NS+"Synonym":
                if state_stack[-5] == NS+"DisorderGeneAssociationList" and state_stack[-3] == "Gene": 
                    synonym = elem.text
                    gene_to_synonyms.setdefault(gene, set()).add(synonym)
            elif elem.tag == NS+"Source":
                if state_stack[-6] == NS+"DisorderGeneAssociationList" and state_stack[-4] == "Gene": 
                    source = elem.text
            elif elem.tag == NS+"Reference": 
                if state_stack[-6] == NS+"DisorderGeneAssociationList" and state_stack[-4] == "Gene": 
                    if source == "SwissProt":
                        if gene in gene_to_uniprot:
                            if gene_to_uniprot[gene] != elem.text:
                                # Only GNAS has inconsistency but the geneids are the same 
                                print "Uniprot inconsistency:", gene, gene_to_uniprot[gene], elem.text
                        gene_to_uniprot[gene] = elem.text
            # Does not work
            elif False and elem.tag == NS+"ExternalReference": #NS+"Source":
                if state_stack[-5] == NS+"DisorderGeneAssociationList" and state_stack[-3] == "Gene": 
                    elem_inner = elem.find(NS+"Source")
                    print tostring(elem_inner)
                    if elem_inner.text == "SwissProt":
                        elem_inner = elem.find(NS+"Reference")
                        gene_to_uniprot[gene] = elem.text
            elem.clear()
            state_stack.pop()
    root.clear()
    print len(disease_to_genes), len(gene_to_synonyms), len(gene_to_uniprot)
    print disease_to_genes["Lymphangioleiomyomatosis"], gene_to_synonyms["TSC1"], gene_to_uniprot["TSC1"]
    return disease_to_genes, gene_to_uniprot


# Does not work due to the unknown tags in xml file
def get_disease_to_genes_bs4(file_name):
    disease_to_genes = {}
    xml_doc = open(file_name) 
    soup = BeautifulSoup(xml_doc, "xml")
    for tag in soup.find_all('Disorder'): #, recursive=False):
        phenotype = None
        for tag_inner in tag.children:
            print tag_inner.name
            if tag_inner.name == "DisorderGeneAssociationList":
                for child in tag_inner.DisorderGeneAssociation.Gene.children:
                    #print child.name
                    if child.name == "Symbol":
                        gene = child.string.encode("ascii", "ignore")
                    if child.name == "ExternalReferenceList":
                        for child_inner in child.children:
                            if child_inner.Source.string == "SwissProt":
                                uniprot = child_inner.Reference.string
            elif tag_inner.name == "Name":
                phenotype = tag_inner.string.encode("ascii", "ignore")
        if phenotype is not None:
            disease_to_genes.setdefault(phenotype, set()).add(gene)
    print len(disease_to_genes)
    return disease_to_genes


if __name__ == "__main__":
    main()


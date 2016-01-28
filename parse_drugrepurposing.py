#import urllib2, os, cPickle, re, time
from bs4 import BeautifulSoup
#from toolbox import text_utilities

def main():
    #drug = "montelukast" 
    #disease = "type 2 diabetes mellitus"
    #disease = "asthma"
    base_dir = "~/Dropbox/sbnb/"
    name, indication = read_repurposing_data(base_dir + "Drug Repurposing Info.html")
    print name
    print indication
    return 


def read_repurposing_data(file_name):
    drug_to_values = {}
    #name = None
    #indication = [] 
    html_doc = open(file_name)
    soup = BeautifulSoup(html_doc, "xml")
    #name = tag.strong.string.encode("ascii", "ignore")
    flag = False
    for tag in soup.find_all('tr'):
        #print tag.name 
	if flag:
	    break
	for tag_p in tag.find_all_next():
	    try:
		if tag_p.name == "td": 
		    if tag_p['class'] is not None and tag_p['class'].startswith("recaptcha"):
			flag = True
			break
	    except:
		print tag.name, tag_p.name
		continue
	    txt = tag_p.get_text(" ").strip() # " " For separating headers from text
	    if txt == "":
		continue
	    txt = " ".join(txt.split()) 
	    txt = txt.encode("ascii", "ignore")
	    if txt == "":
		continue
	    #for i, tag_p in enumerate(tag.next_siblings):
    return drug_to_values #name, indication 
    
	
if __name__ == "__main__":
    main()
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


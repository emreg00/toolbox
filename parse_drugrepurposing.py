##############################################################################
# DrugRepurposing.info parser
#
# eg 22/01/2016
##############################################################################

from bs4 import BeautifulSoup

def main():
    base_dir = "~/Dropbox/sbnb/"
    drug_to_values = read_repurposing_data(base_dir + "Drug Repurposing Info.html")
    print len(drug_to_values), drug_to_values 
    return 


def read_repurposing_data(file_name):
    drug_to_values = []
    html_doc = open(file_name)
    soup = BeautifulSoup(html_doc, "lxml")
    flag = False
    for i, tag in enumerate(soup.find_all('tr')):
	#print i, tag.name 
	values = []
	for tag_p in tag.descendants:
	    #print tag_p.name 
	    if tag_p.name == "td": 
		if tag_p.get('class') is not None and tag_p['class'][0].startswith("recaptcha"):
		    flag = True
		    break
	    if tag_p.name == "input":
		val = tag_p.get('value')
		values.append(val)
	    if tag_p.name == "option" and tag_p.get('selected') is not None:
		val = str(tag_p.get_text())
		#val = val.encode("ascii", "ignore")
		values.append(val)
	if flag: # or i>10: 
	    break
	print values
	if len(values) > 0:
	    drug_to_values.append(values[:-1])
    return drug_to_values  
    
	
if __name__ == "__main__":
    main()
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 



def convert_to_R_string(txt, mapping=[(" ", "."), (",", ""), ("'", "")]):
    #txt = txt.replace(",","").replace(" ", ".").replace("'", "")
    for a, b in mapping:
	txt = txt.replace(a,b)
    return txt


def tokenize_disease_name(disease, exact=True):
    disease = disease.lower()
    disease_mod = disease.replace(" and ", ", ")
    disease_mod = disease.replace("-", ", ")
    phrases = disease_mod.split(",")
    values = []
    for phrase in phrases:
	inner_values = []
	words = phrase.strip().split()
	for i, token in enumerate(words):
	    if token.endswith("'s") or token.endswith("^s") :
		token = token[:-2]
	    if i == len(words) - 1:
		if token[-1] == "s":
		    token = token[:-1]
	    if token in ("disease", "disorder", "syndrome"):
		continue
	    inner_values.append(token)
	if exact:
	    values.append(" ".join(inner_values))
	else:
	    values += inner_values
    return values




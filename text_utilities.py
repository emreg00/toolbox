
try:
    from toolbox.external.negex import negex #sortRules, negTagger
except:
    print "Import error: Negex. Using keyword matching instead"


KEYWORDS_NEGATIVE = [ " not ", " no ", " except ", " exception ", " inappropriate ", " without ", " absence " ]
KEYWORDS_SYMPTOMATIC = [ " protect", " maint", " manage", " symptom", " relie", " palliati", " alleviat" ]


def get_negex_rules(file_name):
    f = open(file_name)
    rules = negex.sortRules(f.readlines())
    return rules


def is_negated(txt, phrase, rules = None):
    negative = False
    if rules is not None:
        tagger = negex.negTagger(sentence = txt, phrases = [ phrase ], rules = rules, negP = False)
	negative = tagger.getNegationFlag() == "negated"
	#txt_tagged = tagger.getNegTaggedSentence()
    else:
	negative, i = in_keywords(txt, KEYWORDS_NEGATIVE)
    return negative


def is_symptomatic(txt):
    return in_keywords(txt, KEYWORDS_SYMPTOMATIC)


def in_keywords(txt, keywords):
    flag = False
    for i, keyword in enumerate(keywords):
	idx = txt.find(keyword)
	if idx != -1:
	    flag = True
	    break
    return flag, i


def convert_to_R_string(txt):
    txt = txt.replace_chars(txt, mapping=[(" ", "."), (",", ""), ("'", ""), ("-", "."), ("/", ".")])
    return txt


def replace_chars(txt, mapping=[(" ", "_"), (",", ""), ("'", ""), ("-", "_"), ("/", "_")]):
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




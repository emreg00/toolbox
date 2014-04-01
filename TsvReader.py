
from FormattedFileProcessor import FormattedFileProcessor

class TsvReader(FormattedFileProcessor):
    """
	Read/process TSV (tab seperated) formatted files
    """
    def __init__(self, input_file_name, delim="\t", inner_delim=None, quotation=None):
	FormattedFileProcessor.__init__(self, input_file_name=input_file_name, input_type="tsv", delim = delim, inner_delim = inner_delim, quotation = quotation)
	return

    def read(self, fields_to_include, keys_to_include=None, merge_inner_values=False):
	"""
	    Read input file in "delim" seperated format line by line into a dictionary.
	    Returns two dictionaries: (column_header, index) dictionary and (first_field_included, other_field_values) dictionary. 
	    Values of second dictionary are of type list corresponding to included_field_values (as a list, that is it is a list of lists)
	    If merge_inner_values is True these included_field_values are expanded using "inner_delim" as a list.

	    fields_to_include: columns that would be included in the dictionary. If None all columns are included.
	    keys_to_include: use only lines whose value of the first column is inside this set. Set None for including all the lines in the file.
	    merge_inner_values: If True expand values in a field as a list using "inner_delim". Upon its usage, issuing a reduce(lambda x,y: x+y, vals) statement is suggested on the vals of dictionary.
	    
	    Ex:
	    id	uniprotaccession
	    81234   P12345, Q23W42
	    65747   A12342
	    81234   Q23W42

	    without merge_inner_values:
	    >>> { '81234': [["P12345, Q23W42"], ["Q23W42"]], '65747': [["A12342"]] }

	    with merge_inner_values: 
	    >>> { '81234': [["P12345", "Q23W42", "Q23W42"]], '65747': [["A12342"]] }

	    >>> reduce(lambda x,y: x+y, [["P12345", "Q23W42"], ["Q23W42"]])
	    >>> ['P12345', 'Q23W42', 'Q23W42']
	"""
	file = open(self.input_file_name)
	# Read header
	line = file.readline()
	cols = [ c.lower().strip(self.quotation) for c in line.strip("\n").split(self.delim) ]
	if fields_to_include is None:
	    first_column = cols[0]
	    fields_to_include = []
	    fields_to_include.extend(cols)
	else:
	    fields_to_include = [ f.lower() for f in fields_to_include ]
	    first_column = fields_to_include[0]
	fields_to_include.remove(first_column) # this will be stored as key
	columns = dict(zip(cols, range(len(cols))))
	if merge_inner_values:
	    if self.inner_delim is None:
		raise Exception("merge_inner_values requires that inner_delim is defined!")
	#print "Start:", columns 
	id_to_values = {}
	line_prev = line
	line = file.readline()
	#vals = line.strip().split(self.delim)
	#print vals
	if isinstance(keys_to_include, list):
	    keys_to_include = set(keys_to_include)
	while line:
	    try:
		vals = line.rstrip("\n").split(self.delim)
		id_ = vals[columns[first_column]].strip(self.quotation)
		if keys_to_include is None or id_ in keys_to_include:
		    new_vals = []
		    if merge_inner_values:
			for f in fields_to_include:
			    new_vals.extend(map(lambda x: x.strip().strip(self.quotation), vals[columns[f]].split(self.inner_delim)))
		    else:
			new_vals = [vals[columns[f]].strip(self.quotation) for f in fields_to_include]
		    id_to_values.setdefault(id_, []).append(new_vals)
	    except:  
		print line_prev, line 
		import traceback
		traceback.print_exc()
		break
	    line_prev = line
	    line = file.readline()
	file.close()
	#print columns, "\n", id_to_values.items()[0]
	column_to_index = dict(zip(fields_to_include, range(len(fields_to_include))))
	#print "End:", columns
	return column_to_index, id_to_values


    def process(self, out_method, fields_to_include, overwrite_keys, keys_to_include):
	"""
	    ! Now sort of OBSOLETE, use read (above) instead !
	    Read and process an input file line by line. If out_method is None a dictionary storing read lines are returned.
	    out_method: method to output columns in current line on the fly in tsv format
	    fields_to_include: columns that would be included in the dictionary or processed with the function
	    overwrite_keys: allows overwriting keys (in case of entries with duplicate primary column values in the file) returning values as list in the dictionary
			    if False, returns list of values as list in the dictionary (each list element corresponding to value list of distinct entries)
	    keys_to_include: use only lines whose value of the first column is inside this set. Set None for including all the lines in the file
	"""
	file = open(self.input_file_name)
	line = file.readline()
	cols = [ c.lower() for c in line.strip("\n").split('\t') ]
	if fields_to_include is None:
	    first_column = cols[0]
	else:
	    fields_to_include = [ f.lower() for f in fields_to_include ]
	    first_column = fields_to_include[0]
	columns = dict(zip(cols, range(len(cols))))
	#print "Start:", columns
	id_to_value = {}
	i=0
	line_prev = line
	line = file.readline()
	vals = line.strip("\n").split('\t')
	#print vals
	while line:
	    try:
		vals = line.strip("\n").split('\t')
		id_ = vals[columns[first_column]]
		if keys_to_include is None or id_ in keys_to_include:
		    if out_method is None:
			if fields_to_include is None:
			    if overwrite_keys:
				id_to_value[id_] = vals
			    else:
				id_to_value.setdefault(id_, []).append(vals)
			else:
			    if overwrite_keys:
				id_to_value[id_] = [ vals[columns[f]] for f in fields_to_include]
			    else:
				id_to_value.setdefault(id_, []).append( [vals[columns[f]] for f in fields_to_include] )
		    else:
			out_method("%s\n" % "\t".join([ vals[columns[f]] for f in fields_to_include ]))
	    except: #Exception, e: 
	    	#print "In: ", __file__, e, vals
		print line_prev, line 
		import traceback
		traceback.print_exc()
	    i+=1
	    #if i > 20:
	    #	break
	    line_prev = line
	    line = file.readline()
	file.close()
	#print columns, "\n", id_to_value.items()[0]
	if out_method is None:
	    if fields_to_include is not None:
		cols2 = []
		for c in cols:
		    if c in fields_to_include:
			cols2.append(c)
		columns = dict(zip(cols2, range(len(cols2))))
	    #print "End:", columns
	    return columns, id_to_value
	else:
	    return



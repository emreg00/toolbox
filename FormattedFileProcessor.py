
class FormattedFileProcessor(object):
    """
	Wrapper around reading/processing various formatted input streams
    """
    allowed_formats = set(["sif", "fasta", "tsv"])
    def __init__(self, input_file_name, input_type, delim, inner_delim=None):
	"""
	    Initialize an object of this class storing 
	    input_file_name: file name to be read/processed
	    input_type: could be sif, tsv, fasta, etc..
	"""
	self.input_file_name = input_file_name
	self.input_type = input_type
	self.delim = delim
	self.inner_delim = inner_delim
	if not self.input_type in self.allowed_formats:
	    raise Exception("Unrecognized input type")
	return

    #def read(self, fields_to_include=None, overwrite_keys=True, keys_to_include=None):
    def read(self, fields_to_include=None, keys_to_include=None):
	"""
	    Read the file into a dictionary and return cloumns names and value dictionary
	"""
	raise Exception("Call method of FormattedInputProcesor abstract class")
	#return self.process(out_method = None, fields_to_include = fields_to_include, overwrite_keys = overwrite_keys, keys_to_include = keys_to_include)

    def process(self, out_method, fields_to_include, overwrite_keys):
	"""
	    Read and process an input file line by line. If out_method is None a dictionary storing read lines are returned.
	    out_method: method to output columns in current line on the fly in the input_type format
	    fields_to_include: columns that would be included in the dictionary or processed with the function
	"""
	raise Exception("Call method of FormattedInputProcesor abstract class")
	return


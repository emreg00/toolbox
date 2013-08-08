import ConfigParser #, os

#CONFIG = configuration.Configuration(CONFIG_FILE, CONFIG_SECTION) 

class Configuration():

    CONFIG_FILE = "default.ini"
    CONFIG_SECTION = "DEFAULT"

    def __init__(self, config_file=None, config_section=None):
	#print os.getcwd()
	if config_file is None:
	    config_file = Configuration.CONFIG_FILE
	if config_section is None:
	    config_section = Configuration.CONFIG_SECTION
	self.config_section = config_section
	self.config = ConfigParser.SafeConfigParser()
	self.config.read(config_file)
	return

    def get(self, var_name):
	return self.config.get(self.config_section, var_name) 

# CONFIG = configuration.get_configuration(CONFIG_FILE, CONFIG_SECTION) 
#def get_configuration(config_file=CONFIG_FILE, config_section=CONFIG_SECTION, var_name=None):
#    config = ConfigParser.SafeConfigParser()
#    config.read(CONFIG_FILE)
#    if var_name is None:
#	return config
#    else:
#	return config.get(CONFIG_SECTION, var_name) 


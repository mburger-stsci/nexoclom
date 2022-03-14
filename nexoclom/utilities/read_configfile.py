import os
from nexoclom.utilities.exceptions import ConfigfileError


DEFAULT_DATABASE = 'thesolarsystemmb'
DEFAULT_PORT = 5432

class NexoclomConfig:
    """Configure external resources used in the model.
    The following parameters can be saved in the file `$HOME/.nexoclom`.
    * savepath = <path where output files are saved>
    * database = <name of the postgresql database to use> (*optional*)
    * port = <port for postgreSQL server to use> (*optional*)
    
    If savepath is not present, an exception is raised
    """
    def __init__(self, configfile=None):
        if configfile is None:
            # print('config', os.environ.get('NEXOCLOMCONFIG', 'Not Set'))
            configfile = os.environ.get('NEXOCLOMCONFIG', os.path.join(
                os.environ['HOME'], '.nexoclom'))
        else:
            pass
        print(f'Using configuration file {configfile}')
        self.configfile = configfile
        
        config = {}
        if os.path.isfile(configfile):
            # Read the config file into a dict
            for line in open(configfile, 'r'):
                if '=' in line:
                    key, value = line.split('=')
                    config[key.strip()] = value.strip()
                else:
                    pass
        else:
            pass

        self.savepath = config.get('savepath', None)
        if self.savepath is None:
            raise ConfigfileError(configfile, self.savepath)
        elif not os.path.exists(self.savepath):
            os.makedirs(self.savepath)
        else:
            pass
        
        self.database = config.get('database', DEFAULT_DATABASE)
        
        if 'port' not in config:
            self.port = DEFAULT_PORT
        else:
            self.port = int(config['port'])
            
        for key, value in config.items():
            if key not in self.__dict__:
                self.__dict__[key] = value
            else:
                pass
            
    def __repr__(self):
        return self.__dict__.__repr__()
    
    def __str__(self):
        return self.__dict__.__str__()

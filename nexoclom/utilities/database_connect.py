import psycopg2
from nexoclom.utilities.read_configfile import read_configfile


DEFAULT_DATABASE = 'thesolarsystemmb'
DEFAULT_PORT = 5432

def database_connect(return_con=True):
    """Wrapper for psycopg.connect() that determines which database and port to use.

    :return:
    :param database: Default = None to use value from config file
    :param port: Default = None to use value from config file
    :param return_con: False to return database name and port instead of connection
    :return: Database connection with autocommit = True unless return_con = False
    """
    config = read_configfile()
    database = config.get('database', DEFAULT_DATABASE)
    port = config.get('port', DEFAULT_PORT)

    if return_con:
        con = psycopg2.connect(dbname=database, port=port)
        con.autocommit = True

        return con
    else:
        return database, port
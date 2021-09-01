import os
import psycopg2


nexoclom_tables = ['geometry_with_time', 'geometry_without_time',
                   'surface_int_constant', 'surface_int_map',
                   'surface_int_tempdependent', 'forces',
                   'spatdist_uniform', 'spatdist_surfmap',
                   'spatdist_spot', 'spatdist_fittedoutput',
                   'speeddist_gaussian', 'speeddist_maxwellian',
                   'speeddist_sputtering', 'speeddist_flat',
                   'speeddist_fittedoutput', 'speeddist_user',
                   'angdist_isotropic', 'options', 'outputfile',
                   'modelimages', 'uvvsmodels']


def database_connect(database=None, port=None, return_con=True):
    """Wrapper for psycopg2.connect() that determines which database and port to use.

    :return:
    :param database: Default = None to use value from config file
    :param port: Default = None to use value from config file
    :param return_con: False to return database name and port instead of connection
    :return: Database connection with autocommit = True unless return_con = False
    """
    configfile = os.path.join(os.environ['HOME'], '.nexoclom')
    config = {}
    if os.path.isfile(configfile):
        for line in open(configfile, 'r').readlines():
            if '=' in line:
                key, value = line.split('=')
                config[key.strip()] = value.strip()
            else:
                pass

        if (database is None) and ('database' in config):
            database = config['database']
        else:
            pass

        if (port is None) and ('port' in config):
            port = int(config['port'])
        else:
            pass
    else:
        pass

    if database is None:
        database = 'thesolarsystemmb'
    else:
        pass

    if port is None:
        port = 5432
    else:
        pass

    if return_con:
        con = psycopg2.connect(database=database, port=port)
        con.autocommit = True

        return con
    else:
        return database, port

def export_database():
    # Get database name and port
    database, port = database_connect(return_con=False)

    if not os.path.exists('backup'):
        os.makedirs('backup')
    else:
        pass
    
    for table in nexoclom_tables:
        print(f'Backup up {table}')
        savef = os.path.join('backup', f'nexoclom_{table}.sql')
        os.system(f'pg_dump -p {port} -t {table} {database} > {savef}')
        
def import_database():
    # Get database name and port
    database, port = database_connect(return_con=False)
    
    for table in nexoclom_tables:
        print(f'Importing {table}')
        savef = os.path.join('backup', f'nexoclom_{table}.sql')
        os.system(f'psql -c "drop table {table};" thesolarsystemmb')
        os.system(f'psql {database} < {savef}')

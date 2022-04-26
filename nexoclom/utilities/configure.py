"""Create and read configuration file, create necessary database tables."""
import os
import pandas as pd
from nexoclom.utilities import NexoclomConfig
from nexoclom import __file__ as basefile


basepath = os.path.dirname(basefile)


def configure_nexoclom(verbose=False):
    # Create the database if necessary
    config = NexoclomConfig(verbose=verbose)
    config.verify_database_running()
    
    # Validate nexoclom output tables
    with config.database_connect() as con:
        cur = con.cursor()
        cur.execute('select table_name from information_schema.tables')
        tables = [r[0] for r in cur.fetchall()]
        cur.execute('''SELECT distinct pg_type.typname AS enum_type
                       FROM pg_type JOIN pg_enum
                       ON pg_enum.enumtypid = pg_type.oid;''')
        tables.extend([r[0] for r in cur.fetchall()])
        
        with open(os.path.join(basepath, 'data', 'schema.sql'), 'r') as sqlfile:
            done = False
            while not done:
                line = sqlfile.readline()
                nextline = ''
                if 'CREATE' in line:
                    # table_to_test = line[len('CREATE TABLE '):-3]
                    table_to_test = line.split()[2].lower()

                    if table_to_test in tables:
                        # Need to verify schema
                        pass
                    else:
                        # Create the table if it isn't there
                        query = line
                        nextline = sqlfile.readline()
                        while (nextline.strip()) and ('DONE' not in nextline):
                            query += nextline
                            nextline = sqlfile.readline()
                        print(query)
                        cur.execute(query)
                done = ('DONE' in nextline) or ('DONE' in line)
    return config

def configure_solarsystem():
    # Make a pickle file with the planetary constants
    pklfile = os.path.join(basepath, 'data', 'PlanetaryConstants.pkl')
    if not os.path.exists(pklfile):
        constfile = os.path.join(basepath, 'data', 'PlanetaryConstants.dat')
        constants = pd.read_csv(constfile, sep=':', comment='#')
        constants.columns = [col.strip() for col in constants.columns]
        constants['Object'] = constants['Object'].apply(lambda x: x.strip())
        constants['orbits'] = constants['orbits'].apply(lambda x: x.strip())
        constants.to_pickle(pklfile)
    else:
        pass

    # Do something with the naif_ids.
    # pklfile = os.path.join(basepath, 'data', 'naif_ids.pkl')
    # if not os.path.exists(pklfile):
        # naiffile = os.path.join(basepath, 'naif_ids.dat')

def configure_atomicdata():
    # Make gvalue table
    gvalue_file = os.path.join(basepath, 'data', 'g-values', 'g-values.pkl')
    if not os.path.exists(gvalue_file):
        from nexoclom.atomicdata.initialize_atomicdata import make_gvalue_table
        make_gvalue_table()
    else:
        pass

    # Make the photorates table
    photorates_file = os.path.join(basepath, 'data', 'Loss', 'photorates.pkl')
    if not os.path.exists(photorates_file):
        from nexoclom.atomicdata.initialize_atomicdata import make_photorates_table
        make_photorates_table()
        
def configure(verbose=False):
    config = configure_nexoclom(verbose=True)
    configure_solarsystem()
    configure_atomicdata()
    
    return config

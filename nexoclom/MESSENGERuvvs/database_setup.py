import os
import os.path
import sys
import glob
import psycopg2
import types


def database_connect(database=None, port=None, return_con=True):
    """Return a database connection to saved atomic data
    Wrapper for ``psycopg2.connect()`` that determines database and port to use.

    **Parameters**

    database
        Database to connect to. If not given, it must be supplied in
        the $HOME/.nexoclom configuration file.

    port
        Port the database server uses. If not given, it must be supplied in
        the $HOME/.nexoclom configuration file.

    return_con
        False to return database name and port instead of connection.
        Default = True

    **Returns**

    Database connection with autocommit = True unless return_con = False

    **Examples**
    ::

        >>> from atomicdataMB import database_connect
        >>> database, port = database_connect(return_con=False)
        >>> print(f'database = {database}; port = {port}')
        database = thesolarsystemmb; port = 5432
        >>> with database_connect() as con:
        ...     cur = con.cursor()
        ...     cur.execute('SELECT DISTINCT species from gvalues')
        ...     species = cur.fetchall()
        >>> species = [s[0] for s in species]
        >>> print(species)
        ['Ca', 'OH', 'O', 'Ti', 'C', 'Mg+', 'Na', 'Mg', 'H', 'Mn', 'He',
         'Ca+', 'K', 'S']

    """
    configfile = os.path.join(os.environ['HOME'], '.nexoclom')
    config = {}
    if os.path.isfile(configfile):
        for line in open(configfile, 'r'):
            key, value = line.split('=')
            config[key.strip()] = value.strip()
        
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
    
def messenger_database_setup(force=False):
    """Setup the database from SQL database dump files.
    Repopulates the database using a SQL backup rather than the original
    IDL save files. See :doc:`database_fields` for a description of the
    tables and fields used by MESSENGERuvvs.
    
    **Parameters**
    
    force
        If True, deletes old database tables and remakes them.
        Default is False, which only creates the tables if necessary.
        
    **Returns**
    
    No output.
    
    """
    # Get database name and port
    database, port = database_connect(return_con=False)

    if ((isinstance(sys.modules['psycopg2'], types.ModuleType)) and
        ('test' not in database)):
        # Read in current config file if it exists
        configfile = os.path.join(os.environ['HOME'], '.nexoclom')
        datapath = None
        if os.path.isfile(configfile):
            for line in open(configfile, 'r').readlines():
                key, value = line.split('=')
                if key.strip() == 'datapath':
                    datapath = value.strip()
                else:
                    pass
        else:
            pass

        if datapath is None:
            datapath = input('What is the path to the MESSENGER data? ')
            with open(configfile, 'a') as f:
                f.write(f'datapath = {datapath}\n')
        else:
            pass

        # Verify database is running
        status = os.popen('pg_ctl status').read()
        if 'no server running' in status:
            os.system(f'pg_ctl start -D $HOME/.postgres/main '
                      f'-l $HOME/.postgres/logfile -o "-p {port}"')
        else:
            pass

        # Create MESSENGER database if necessary
        with database_connect(database='postgres') as con:
            cur = con.cursor()
            cur.execute('select datname from pg_database')
            dbs = [r[0] for r in cur.fetchall()]

            if database not in dbs:
                print(f'Creating database {database}')
                cur.execute(f'create database {database}')
            else:
                pass

        # Create the MESSENGER tables if necessary
        with database_connect() as con:
            cur = con.cursor()
            cur.execute('select table_name from information_schema.tables')
            tables = [r[0] for r in cur.fetchall()]

            mestables = ['capointing', 'cauvvsdata', 'caspectra',
                         'mgpointing', 'mguvvsdata', 'mgspectra',
                         'napointing', 'nauvvsdata', 'naspectra',
                         'mesmercyear']
            there = [m in tables for m in mestables]

            if (False in there) or force:
                # Delete any tables that may exist
                for mestab in mestables:
                    if mestab in tables:
                        cur.execute(f'drop table {mestab}')
                    else:
                        pass

                # Import the dumped tables
                datafiles = glob.glob(datapath+'/UVVS*sql')
                for dfile in datafiles:
                    print(f'Loading {os.path.basename(dfile)}')
                    os.system(f'psql -d {database} -p {port} -f {dfile}')
            else:
                pass

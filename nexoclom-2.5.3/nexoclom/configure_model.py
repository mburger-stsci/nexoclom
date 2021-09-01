"""Create and read configuration file, create necessary database tables."""
import os
import os.path
from .database_connect import database_connect
import sys
import types


def configfile():
    """Configure external resources used in the model.

    The following parameters can be saved in the file `$HOME/.nexoclom`.
    
    * savepath = <path where output files are saved>
    
    * datapath = <path where MESSENGER data is kept>
    
    * database = <name of the postgresql database to use> (*optional*)
    
    * port = <port for postgreSQL server to use> (*optional*)
    
    If savepath and datapath are not present, user is prompted to enter them.
    """
    # Determne the savepath
    cfile = os.path.join(os.environ['HOME'], '.nexoclom')
    savepath = None
    if os.path.isfile(cfile):
        for line in open(cfile, 'r').readlines():
            key, value = line.split('=')
            if key.strip() == 'savepath':
                savepath = value.strip()
            else:
                pass
    else:
        pass

    if savepath is None:
        savepath = input('Where should outputs be saved: ')
        with open(cfile, 'a') as f:
            f.write(f'savepath = {savepath}\n')

        # Create save directory if necessary
        if not os.path.isdir(savepath):
            try:
                os.makedirs(savepath)
            except:
                assert 0, f'Could not create directory {savepath}'

    # Create the database if necessary
    database, port = database_connect(return_con=False)
    with database_connect(database='postgres') as con:
        cur = con.cursor()
        cur.execute('select datname from pg_database')
        dbs = [r[0] for r in cur.fetchall()]

        if database not in dbs:
            print(f'Creating database {database}')
            cur.execute(f'create database {database}')
        else:
            pass

    return savepath


def verify_output_tables():
    """Create the database tables used by nexoclom to save output."""
    with database_connect() as con:
        cur = con.cursor()
        cur.execute('select table_name from information_schema.tables')
    tables = [r[0] for r in cur.fetchall()]
    
    with open(os.path.join(os.path.dirname(__file__),
                           'data',
                           'schema.sql'), 'r') as sqlfile:
        done = False
        while not done:
            line = sqlfile.readline()
            nextline = ''
            if 'TABLE' in line:
                table_to_test = line[len('CREATE TABLE '):-3]
                if table_to_test in tables:
                    # Need to verify schema
                    pass
                else:
                    query = line
                    nextline = sqlfile.readline()
                    while (nextline.strip()) and ('DONE' not in nextline):
                        query += nextline
                        nextline = sqlfile.readline()
                    print(query)
                    cur.execute(query)
            done = ('DONE' in nextline) or ('DONE' in line)


def configure_model():
    """Ensure the database and configuration file are set up for nexoclom.
    
    **Parameters**
    
    No parameters.
    
    **Returns**
    
    No output.
    """
    if isinstance(sys.modules['psycopg2'], types.ModuleType):
        # Get database name and port
        database, port = database_connect(return_con=False)

        # Verify database is running
        status = os.popen('pg_ctl status').read()
        if 'no server running' in status:
            os.system(f'pg_ctl start -D $HOME/.postgres/main '
                      f'-l $HOME/.postgres/logfile -o "-p {port}"')
        else:
            pass
        
        with database_connect() as con:
            # Check whether SSObject has been created
            cur = con.cursor()
            query = """
                    SELECT exists
                    (select 1 from pg_type where typname = 'ssobject');
                    """
            cur.execute(query)
            result = cur.fetchall()

        if not result[0][0]:
            from solarsystemMB import create_SSObject
            create_SSObject()

        else:
            pass
        
        verify_output_tables()

from .database_connect import database_connect


def configure_database():
    """Set up the database tables used by nexoclom.
    Tables used are defined in :doc:`database_tables`.
    """
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
        
    with database_connect() as con:
    

"""Backup the MESSENGERuvvs database tables."""
import os
from .database_setup import database_connect


def databasebackups():
    """Backup the MESSENGERuvvs database tables.
    
    Dump the MESSENGERuvvs data into SQL files that can be restored if
    necessary. Tables that are backed-up are: cauvvsdata, capointing,
    mguvvsdata, mgpointing, nauvvsdata, napointing, mesmercyear.
    
    This function takes no arguments. The path to save the database dumps
    must be specified in $HOME/.nexoclom in the format
    ``datapath = <path to data>``. The default database and port are
    ``thesolarsystemmb`` and ``5432``. These can also be specified in the
    .nexoclom file.
    """
    
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
    assert datapath is not None, 'Undefined datapath.'

    # Get database name and port
    database, port = database_connect(return_con=False)

    mestables = ['capointing', 'cauvvsdata', 'caspectra',
                 'mgpointing', 'mguvvsdata', 'mgspectra',
                 'napointing', 'nauvvsdata', 'naspectra',
                 'mesmercyear']

    for table in mestables:
        print(f'Backing up {table}')
        savef = os.path.join(datapath, f'UVVS_{table}.sql')
        os.system(f"pg_dump -p {port} -t {table} {database} > {savef}")

if __name__ == '__main__':
    databasebackups()

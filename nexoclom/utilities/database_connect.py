import os
import psycopg
import subprocess
from nexoclom.utilities.read_configfile import NexoclomConfig


def verify_database_running(config=None, configfile=None):
    if config is None:
        config = NexoclomConfig(configfile=configfile)
    else:
        pass
    
    # verify database is running
    proc = subprocess.run('pg_ctl status', capture_output=True, shell=True)
    if 'no server running' in str(proc.stdout):
        subprocess.run(f'pg_ctl -o "-p {config.port}" start -l {os.environ["PGDATA"]}/logfile',
                       shell=True)
        return 'Started Database'
    else:
        return 'Database Already Running'


def database_connect(configfile=None):
    """Wrapper for psycopg.connect() that determines which database and port to use.

    :return:
    :param database: Default = None to use value from config file
    :param port: Default = None to use value from config file
    :param return_con: False to return database name and port instead of connection
    :return: Database connection with autocommit = True unless return_con = False
    """
    config = NexoclomConfig(configfile=configfile)
    verify_database_running(config)
    
    con = psycopg.connect(dbname=config.database, port=config.port)
    con.autocommit = True

    return con
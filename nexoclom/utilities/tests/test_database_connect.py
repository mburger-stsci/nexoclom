import subprocess
from nexoclom.utilities import database_connect
from nexoclom.utilities.database_connect import verify_database_running
import pytest


@pytest.mark.utilities
def test_database_connect():
    # Verify that it can start the postgres server
    proc = subprocess.run('pg_ctl stop', capture_output=True, shell=True)
    assert (('server stopped' in str(proc.stdout)) or
            ('Is server running') in  str(proc.stdout))
    
    with database_connect() as con:
        assert con.autocommit
        assert con.info.dbname == 'thesolarsystemmb_test'

    # Test database connection when server already running
    with database_connect() as con:
        assert con.autocommit
        assert con.info.dbname == 'thesolarsystemmb_test'

@pytest.mark.utilities
def test_verify_database_running():
    proc = subprocess.run('pg_ctl stop', capture_output=True, shell=True)
    assert (('server stopped' in str(proc.stdout)) or
            ('Is server running') in str(proc.stdout))
    
    # Database stopped; needs to start
    output = verify_database_running()
    assert output == 'Started Database'
    
    # Database already running
    output = verify_database_running()
    assert output == 'Database Already Running'


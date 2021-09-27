import os
from ..database_connect import database_connect

def test_database_connect():
    # Test database connection
    with database_connect() as con:
        assert con.autocommit
        
    # Test picks right defaults
    assert database_connect(return_con=False) == ('thesolarsystemmb', 5432)
    
    # Test non-defaults
    assert database_connect(database='doesntexist', port=1111,
                            return_con=False) == ('doesntexist', 1111)
    
    # Test no config file
    os.rename(os.path.join(os.environ['HOME'], '.nexoclom'), '_temp')
    with database_connect() as con:
        assert con.autocommit
    os.rename('_temp', os.path.join(os.environ['HOME'], '.nexoclom'))
    
    

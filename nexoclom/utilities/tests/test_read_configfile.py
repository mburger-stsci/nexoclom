from nexoclom.utilities.read_configfile import NexoclomConfig
import os
import shutil
import pytest

@pytest.mark.utilities
def test_read_configfile():
    # Test with environment variable
    if os.path.exists('/Users/mburger/Work/Research/modeloutputs_test'):
        shutil.rmtree('/Users/mburger/Work/Research/modeloutputs_test')
    result = {'savepath':'/Users/mburger/Work/Research/modeloutputs_test',
              'database':'thesolarsystemmb_test',
              'port':5432}
    config = NexoclomConfig()
    assert config.savepath == result['savepath']
    assert config.database == result['database']
    assert config.port == result['port']

    # Test with defaults
    del os.environ['NEXOCLOMCONFIG']
    result = {'savepath':'/Users/mburger/Work/Research/modeloutputs',
              'database':'thesolarsystemmb',
              'port':5432}
    config = NexoclomConfig()
    assert config.savepath == result['savepath']
    assert config.database == result['database']
    assert config.port == result['port']
    os.environ['NEXOCLOMCONFIG'] = os.path.join(os.environ['HOME'], '.nexoclom_test')

    # Test with different port
    result = {'savepath':'/Users/mburger/Work/Research/modeloutputs2',
              'database':'thesolarsystemmb',
              'port':1234}
    config = NexoclomConfig(os.path.join(os.environ['HOME'],
                                         '.nexoclom_test2'))
    assert config.savepath == result['savepath']
    assert config.database == result['database']
    assert config.port == result['port']
    
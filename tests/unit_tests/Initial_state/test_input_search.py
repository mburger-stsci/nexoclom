import os
import pytest
import random
import glob
from nexoclom import Input, Output, __path__


# inputfiles = ('Ca.reference.input', 'Na.reference.input',
#               'Ca.surfacemap.maxwellian.input', 'Na.surfacemap.maxwellian.input')
basepath = os.path.dirname(__path__[0])
inputfiles = glob.glob(os.path.join(basepath, 'tests', 'test_data',
                                    'inputfiles', '*.input'))
n_its = [random.randint(1, 10) for _ in inputfiles]

@pytest.mark.initial_state
@pytest.mark.parametrize('inputfile, n_it', zip(inputfiles, n_its))
def test_input_search(inputfile, n_it):
    # inputfile_ = os.path.join(basepath, 'tests', 'test_data', 'inputfiles',
    #                           inputfile)
    # inputs = Input(inputfile_)
    print(inputfile)
    
    inputs = Input(inputfile)
    if inputs.spatialdist.type == 'surface map':
        inputs.spatialdist.coordinate_system = 'solar-fixed'
    else:
        pass
    
    # Remove anything that might be there
    inputs.delete_files()
    
    outputs = [Output(inputs, 1000, run_model=False) for _ in range(n_it)]
    outputfiles = set(output.filename for output in outputs)

    search_results = inputs.search()
    assert len(search_results[0]) == n_it, 'Did not find all the files'
    assert set(search_results[1]) == outputfiles, 'Did not find correct files'
    assert search_results[2] == 1000 * n_it
    assert search_results[3] == 1000 * n_it
    for outputfile in outputfiles:
        assert os.path.exists(outputfile), f'{outputfile} not found'

    # Cleanup
    inputs.delete_files()
    for outputfile in outputfiles:
        assert not os.path.exists(outputfile), f'{outputfile} not deleted'

if __name__ == '__main__':
    for inputfile, n_it in zip(inputfiles, n_its):
        test_input_search(inputfile, n_it)

import sys
import pickle
from nexoclom.data_simulation.compute_iteration import compute_iteration


datafile = sys.argv[1]
with open(datafile, 'rb') as file:
    self, outputfile, scdata = pickle.load(file)
 
iteration_result = compute_iteration(self, outputfile, scdata)
savefile = datafile.replace('inputs', 'outputs').replace('.pkl', '.txt')
with open(savefile, 'w') as file:
    file.write(iteration_result.modelfile)

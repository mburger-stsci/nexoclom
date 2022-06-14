import sys
import pickle
from nexoclom import Output


datafile = sys.argv[1]
npackets = int(sys.argv[2])
with open(datafile, 'rb') as file:
    inputs = pickle.load(file)

Output(inputs, npackets)

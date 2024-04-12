import pandas as pd

gvalues = pd.read_pickle('g-values.pkl')
gvalues = gvalues[gvalues.species != 'Mg']
gvalues = gvalues[gvalues.species != 'Na']
gvalues = gvalues[gvalues.species != 'Ca']

species = (('Mg', 2852), ('Na', 5897), ('Na', 5891), ('Ca', 4227))
for sp, lam in species:
    data = pd.read_csv(f'{sp}.Killen2022.csv')
    data['species'] = sp

    data['wavelength'] = lam
    data['refpoint'] = 0.352
    data['filename'] = f'{sp}.Killen2022.csv'
    data['reference'] = 'Killen et al. 2022'
    data.rename(columns={str(lam): 'gvalue'}, inplace=True)
    gvalues = pd.concat([gvalues, data])


from inspect import currentframe, getframeinfo
frameinfo = getframeinfo(currentframe())
print(frameinfo.filename, frameinfo.lineno)
from IPython import embed; embed()
import sys; sys.exit()

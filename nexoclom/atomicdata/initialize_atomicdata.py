"""Read atomicdata from the text files and save as pandas dataframes
"""
import os
import glob
from nexoclom import __file__ as basefile
import pandas as pd

basepath = os.path.dirname(basefile)

def make_gvalue_table():
    ref = 'Killen et al. (2009)'
    datafiles = glob.glob(os.path.join(basepath, 'data', 'g-values', '*.dat'))
    
    gvalues = pd.DataFrame(columns=['species', 'wavelength', 
                                    'velocity', 'gvalue', 'refpoint', 
                                    'filename', 'reference'])
    
    for datafile in datafiles:
        # Determine the species
        species = os.path.basename(datafile).split('.')[0]

        with open(datafile) as f:
            # Determine the reference point
            refpt_str = f.readline().strip()
        refpt = float(refpt_str.split('=')[1])

        gvalue_species = pd.read_csv(datafile, sep=':', skiprows=1)
        wavelengths = [float(wave) for wave in gvalue_species.columns[1:]]
        gvalue_species.columns = ['vel'] + wavelengths

        for wave in wavelengths:
            print(species, wave)
            for _, row in gvalue_species.iterrows():
                newrow = {'species': species,
                          'wavelength': wave,
                          'velocity': row['vel'],
                          'gvalue': row[wave],
                          'refpoint': refpt,
                          'filename': datafile,
                          'reference': ref}
                gvalues.loc[len(gvalues)] = newrow
        
    gvalue_file = os.path.join(basepath, 'data', 'g-values', 'g-values.pkl')
    print(gvalue_file)
    gvalues.to_pickle(gvalue_file)
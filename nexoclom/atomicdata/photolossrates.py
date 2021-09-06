"""``photolossrates`` - Determine photoionization and photodissociation rates"""
import os, os.path
import glob
import pandas as pd
import astropy.units as u
from .database_connect import database_connect


def make_photo_table():
    """ Creates and populates photorates database table.
    Creates a database table called photorates. Fields in the table:
    
    filename
        Filename in project tree containing the data; used only for
        populating the database
        
    reference
        Source of the photoionization or photodissociation rate
        
    species
        Atomic or molecular species
        
    reaction
        Photoionization or photodissociation reaction
        
    If multiple reaction rates are found for a reaction, user is prompted
    to choose the best one. Most of the reactions are in:
    Huebner & Mukherjee (2015), Astrophys. Space Sci., 195, 1-294.
    
    **Parameters**
    
    None
    
    **Returns**
    
    No output
    
    """
    con = database_connect()

    # drop the old table if necessary
    cur = con.cursor()
    cur.execute('select table_name from information_schema.tables')
    tables = [r[0] for r in cur.fetchall()]
    if 'photorates' in tables:
        cur.execute('''DROP TABLE photorates''')
    else:
        pass

    # Make the photorates table
    cur.execute('''CREATE TABLE photorates (
                     filename text,
                     reference text,
                     species text,
                     reaction text,
                     kappa float,
                     bestvalue boolean)''')

    photodatafiles = glob.glob(os.path.join(os.path.dirname(__file__), 'data',
                                            'Loss', 'Photo', '*.dat'))

    for f in photodatafiles:
        print(f'  {f}')
        ref = ''
        for line in open(f):
            if 'reference' in line.lower():
                ref = line.split('//')[0].strip()
            # elif 'datatype' in line.lower():
            #     dtype = line.split('//')[0].strip()
            # elif 'reactype' in line.lower():
            #     rtype = line.split('//')[0].strip()
            # elif 'ratecoefunits' in line.lower():
            #     un = line.split('//')[0].strip()
            elif len(line.split(':')) == 4:
                parts = line.split(':')
                sp = parts[0].strip()
                reac = parts[1].strip()
                kappa = parts[2].strip()

                cur.execute('''INSERT INTO photorates values(
                                   %s, %s, %s, %s, %s, %s)''',
                            (f, ref, sp, reac, kappa, False))

    # Look for duplicates
    cur.execute('SELECT DISTINCT reaction from photorates')
    temp = cur.fetchall()
    rlist = [t[0] for t in temp]
    for r in rlist:
        print(r)
        cur.execute('SELECT reference from photorates where reaction=%s',
                    (r, ))
        if cur.rowcount > 1:
            temp = cur.fetchall()
            refs = [a[0] for a in temp]
            print('Reaction = {}'.format(r))
            for i, a in enumerate(refs):
                print('({}) {}'.format(i, a))
            q = 0
            cur.execute('''UPDATE photorates
                           SET bestvalue=True
                           WHERE reaction=%s and reference=%s''',
                        (r, refs[q]))
        else:
            cur.execute('''UPDATE photorates
                           SET bestvalue=True
                           WHERE reaction=%s''',
                        (r, ))
    con.close()


class PhotoRate:
    r"""Determine photoreactions and photorates for a species.

    **Parameters**

    species
        Species to compute rates for.

    aplanet
        Distance from the Sun. Default is 1 AU. Given as either a numeric
        type or an astropy quantity with length units.

    **Class Attributes**

    species
        Species

    aplanet
        Distance from the Sun; astropy quantity with units AU

    rate
        Reaction rate; astropy quantity with units s^{-1}. Rate is the sum
        of all possible reactions for the species.

    reactions
        Pandas dataframe with columns for reaction and rate (in s^{-1}) for
        each reaction for the species. This can be used to determine the
        products produced by photolysis and photoionization.
        
    **Example**
    ::
        >>> from atomicdataMB import PhotoRate
        >>> kappa = PhotoRate('Na', 0.33)
        >>> print(kappa)
        Species = Na
        Distance = 0.33 AU
        Rate = 6.666666666666666e-05 1 / s
        >>> print(kappa.rate)
        6.666666666666666e-05 1 / s
        >>> print(kappa.reactions)
                       reaction                  kappa
        0  Na, photon -> Na+, e  6.666666666666666e-05
        >>> kappa = PhotoRate('H_2O')
        >>> print(kappa)
        Species = H_2O
        Distance = 1.0 AU
        Rate = 1.2056349999999999e-05 1 / s
        >>> print(kappa.reactions)
                             reaction     kappa
        0      H_2O, photon -> H_2, O  5.97e-07
        1       H_2O, photon -> OH, H  1.03e-05
        2     H_2O, photon -> H, H, O  7.54e-07
        3   H_2O, photon -> H, OH+, e  5.54e-08
        4   H_2O, photon -> OH, H+, e  1.31e-08
        5    H_2O, photon -> H_2O+, e  3.31e-07
        6  H_2O, photon -> H_2, O+, e  5.85e-09

    """
    def __init__(self, species, aplanet_=1.*u.AU):
        with database_connect() as con:
            prates = pd.read_sql(
                f'''SELECT reaction, kappa
                    FROM photorates
                    WHERE species='{species}' and bestvalue=True''', con)

        try:
            aplanet_.value
            aplanet = aplanet_
        except:
            aplanet = aplanet_*u.AU

        self.species = species
        self.aplanet = aplanet
        a0 = 1*u.AU

        # Photo rate adjusted to proper heliocentric distance
        if len(prates) == 0:
            print('No photoreactions found')
            self.reactions = None
            self.rate = 1e-30/u.s
        else:
            prates['kappa'] = prates['kappa'].apply(
                lambda k: k * (a0/aplanet)**2)
            self.reactions = prates
            self.rate = prates['kappa'].sum()/u.s

    def __str__(self):
        output = (f'Species = {self.species}\n'
                  f'Distance = {self.aplanet}\n'
                  f'Rate = {self.rate}')
        return output

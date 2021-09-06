"""Create and read configuration file, create necessary database tables."""
import os
import numpy as np
import pandas as pd
import psycopg2
import glob
from astropy.io import ascii
from .database_connect import database_connect


curpath = os.path.dirname(__file__)
basepath = os.path.sep + os.path.join(*curpath.split(os.path.sep)[:-1])


def configure_nexoclom():
    # Create the database if necessary
    database, port = database_connect(return_con=False)
    with psycopg2.connect(database=database, port=port) as con:
        con.autocommit = True
        cur = con.cursor()
        cur.execute('select datname from pg_database')
        dbs = [r[0] for r in cur.fetchall()]

        if database not in dbs:
            print(f'Creating database {database}')
            cur.execute(f'create database {database}')
        else:
            pass

    # Validate nexoclom output tables
    with database_connect() as con:
        cur = con.cursor()
        cur.execute('select table_name from information_schema.tables')
    tables = [r[0] for r in cur.fetchall()]
    
    with open(os.path.join(basepath, 'data', 'schema.sql'), 'r') as sqlfile:
        done = False
        while not done:
            line = sqlfile.readline()
            nextline = ''
            if 'TABLE' in line:
                # table_to_test = line[len('CREATE TABLE '):-3]
                table_to_test = line.split()[2]
                if table_to_test in tables:
                    # Need to verify schema
                    pass
                else:
                    # Create the table if it isn't there
                    query = line
                    nextline = sqlfile.readline()
                    while (nextline.strip()) and ('DONE' not in nextline):
                        query += nextline
                        nextline = sqlfile.readline()
                    print(query)
                    cur.execute(query)
            done = ('DONE' in nextline) or ('DONE' in line)

def configure_solarsystem():
    # Make a pickle file with the planetary constants
    pklfile = os.path.join(basepath, 'data', 'PlanetaryConstants.pkl')
    if not os.path.exists(pklfile):
        constfile = os.path.join(basepath, 'data', 'PlanetaryConstants.dat')
        constants = pd.read_csv(constfile, sep=':', comment='#')
        constants.columns = [col.strip() for col in constants.columns]
        constants['Object'] = constants['Object'].apply(lambda x: x.strip())
        constants['orbits'] = constants['orbits'].apply(lambda x: x.strip())
        constants.to_pickle(pklfile)
    else:
        pass

    # Do something with the naif_ids.
    # pklfile = os.path.join(basepath, 'data', 'naif_ids.pkl')
    # if not os.path.exists(pklfile):
        # naiffile = os.path.join(basepath, 'naif_ids.dat')

def configure_atomicdata():
    # Make gvalue table
    with database_connect() as con:
        cur = con.cursor()
        cur.execute('select table_name from information_schema.tables')
    tables = [r[0] for r in cur.fetchall()]
    if 'gvalues' not in tables:
        # Create the table
        with database_connect() as con:
            cur = con.cursor()
            cur.execute('''CREATE TABLE gvalues (
                            filename text,
                            reference text,
                            species text,
                            refpt float, -- AU
                            wavelength float, -- A
                            velocity float[], -- km/s
                            g float[])''')  # 1/s

        # Look up the gvalue datafiles
        datafiles = glob.glob(os.path.join(basepath, 'data', 'g-values', '*.dat'))
        ref = 'Killen et al. (2009)'
    
        for d in datafiles:
            # Determine the species
            f = os.path.basename(d)
            sp = f.split('.')[0]

            # Determine reference point for the file
            with open(d, 'r') as f:
                # Determine the reference point
                astr = f.readline().strip()
                a = float(astr.split('=')[1])

                # Determine the wavelengths
                ww = f.readline().strip().split(':')[1:]
                wavestr = [w.strip() for w in ww]

            # Read in the data table
            data = ascii.read(d, delimiter=':', header_start=1)

            # make the vel array
            vel = np.array(data['vel'])
            q = np.argsort(vel)
            vel = vel[q]

            # Make an array of g-values for each wavelength and add the row
            for w in wavestr:
                gg = np.array(data[w])
                gg = gg[q]
                wave = float(w.strip())
                print(d, sp, wave)

                with database_connect() as con:
                    cur = con.cursor()
                    cur.execute('''INSERT into gvalues values (
                                %s, %s, %s, %s, %s, %s, %s)''',
                                (d, ref, sp, a, wave, list(vel), list(gg)))
    else:
        pass

    # Make the photorates table
    if 'photorates' not in tables:
        # Make the photorates table
        with database_connect() as con:
            cur = con.cursor()
            cur.execute('''CREATE TABLE photorates (
                            filename text,
                            reference text,
                            species text,
                            reaction text,
                            kappa float,
                            bestvalue boolean)''')

        photodatafiles = glob.glob(os.path.join(basepath, 'data', 'Loss', 
                                                'Photo', '*.dat'))

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

                    with database_connect() as con:
                        cur = con.cursor()
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
            

def configure():
    configure_nexoclom()
    configure_solarsystem()
    configure_atomicdata()
    # configure_MESSENGERuvvs()
    

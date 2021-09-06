import os
import numpy as np
import pandas as pd
import pickle
import glob
from scipy import io
from astropy.time import Time
from astropy import units as u
from solarsystemMB import SSObject, planet_geometry
from .database_setup import database_connect


def merc_year(datatime=None, initialize=False):
    """Insert/read start date for each Mercury year from database.

    This creates and reads from database table *MESmercyear*
    """
    
    tstart = Time('2011-03-18T00:00:00', format='isot', scale='utc')
    tend = Time('2015-04-30T23:59:59', format='isot', scale='utc')
    
    if initialize:
        times_ = np.arange(tstart.jd, tend.jd)
        times = [Time(t, format='jd', scale='utc') for t in times_]
        
        taa = np.ndarray((len(times),))*u.rad
        for i, t in enumerate(times):
            time = Time(t, format='jd', scale='utc')
            geo = planet_geometry(time, 'Mercury')
            taa[i] = geo['taa']
        
        styear = [times[0]]
        for a, b, c in zip(taa[0:-1], taa[1:], times[1:]):
            if a > b:
                styear.append(c)
                print(c.iso)
        endyr = [*styear[1:], tend]
        
        with database_connect() as con:
            cur = con.cursor()
            try:
                cur.execute('DROP table MESmercyear')
            except:
                pass
            
            print('creating MESmercyear')
            cur.execute('''CREATE table MESmercyear
                             (yrnum int PRIMARY KEY,
                              yrstart timestamp,
                              yrend timestamp)''')
            for i, d in enumerate(zip(styear, endyr)):
                cur.execute(f'''INSERT into MESmercyear
                                values ({i}, '{d[0].iso}', '{d[1].iso}')''')
    else:
        pass
    
    if datatime is not None:
        with database_connect() as con:
            yrnum = pd.read_sql('''SELECT * from MESmercyear''', con)
        
        myear = np.ndarray((len(datatime),), dtype=int)
        for i, yr in yrnum.iterrows():
            q = (datatime > yr.yrstart) & (datatime < yr.yrend)
            myear[q] = yr.yrnum
        
        return myear
    else:
        return None
    
    
def initialize_MESSENGERdata(datapath):
    """Store data from IDL summary files in a database.
    The IDL summary files were provided by Aimee Merkel.
    
    Two tables are created for each species (Ca, Na, and Mg): *xxuvvsdata* and
    *xxuvvspointing* where xx is the species. See :doc:`database_fields` for
    a description of these tables and fields.
    
    **Parameters**
    
    datapath
        Path to the IDL summary files
        
    **Returns**
    
    No output.
    """
    mercury = SSObject('Mercury')
    
    # Add to the database
    with database_connect() as con:
        cur = con.cursor()
        cur.execute('select table_name from information_schema.tables')
        tables = [r[0] for r in cur.fetchall()]

        mestables = ['capointing', 'cauvvsdata', 'mgpointing',
                     'mguvvsdata', 'napointing', 'nauvvsdata',
                     'caspectra', 'naspectra', 'mgspectra',
                     'uvvsmodels_oribt', 'uvvsmodels_query']
        
        # Delete any tables that may exist
        for mestab in mestables:
            if mestab in tables:
                cur.execute(f'drop table {mestab}')
            else:
                pass
            
        # print('creating MESmercyear table')
        # merc_year(initialize=True)

        print('creating UVVS tables')
        spec = ['Ca', 'Na', 'Mg']
        for sp in spec:
            # Table with spectrum information
            print(f'Creating {sp}uvvsdata')
            cur.execute(f'''CREATE table {sp}uvvsdata (
                               unum SERIAL PRIMARY KEY,
                               species text,
                               frame text,
                               UTC timestamp,
                               orbit int,
                               merc_year int,
                               taa float,
                               rmerc float,
                               drdt float,
                               subslong float,
                               g float,
                               radiance float,
                               sigma float)''')
            
            # Table with MESSENGER geometry and UVVS pointing
            print(f'Creating {sp}pointing')
            cur.execute(f'''CREATE table {sp}pointing (
                               pnum SERIAL PRIMARY KEY,
                               x float,
                               y float,
                               z float,
                               xbore float,
                               ybore float,
                               zbore float,
                               obstype text,
                               obstype_num int,
                               xtan float,
                               ytan float,
                               ztan float,
                               rtan float,
                               alttan float,
                               longtan float,
                               lattan float,
                               loctimetan float,
                               slit text)''')  # Not including slit corners
            
            # Table with spectra
            print(f'Creating {sp}spectra')
            cur.execute(f'''CREATE table {sp}spectra (
                                snum SERIAL PRIMARY KEY,
                                wavelength float[],
                                calibrated float[],
                                raw float[],
                                dark float[],
                                solarfit float[])''')

    savfiles = glob.glob(datapath+'/*_temp.pkl')
    savfiles = sorted(savfiles)
    for oldfile in savfiles:
        # realfile = oldfile.replace('.sav', '_temp.pkl')
        # newfile = oldfile.replace('.sav', '.pkl')
        newfile = oldfile.replace('_temp', '')
        print(f'{oldfile}\n{newfile}\n***')
        # data = io.readsav(oldfile, python_dict=True)
        # data = pickle.load(open(realfile, 'rb'))
        data = pickle.load(open(oldfile, 'rb'))

        kR = u.def_unit('kR', 1e3*u.R)
        Rmerc = u.def_unit('R_Mercury', mercury.radius)
        nm = u.def_unit('nm', 1e-9*u.m)
        
        npts = len(data['orb_num'])
        species = os.path.basename(oldfile)[0:2].lower()
        
        # Determine UT for each spectrum
        t_iso = ['{}:{}:{}'.format('20'+time[0:2].decode('utf-8'),
                                   time[2:5].decode('utf-8'),
                                   time[6:].decode('utf-8'))
                 for time in data['step_utc_time']]
        UTC = Time(t_iso, format='yday')
        
        # Orbit number for each data spectrum
        orbit = np.array([int(o) for o in data['orb_num']])
        
        # determine Mercury year
        myear = merc_year(UTC)
        rmerc = (np.sqrt(np.sum(data['planet_sun_vector_tg']**2,
                                axis=1))*u.km).to(u.AU)
        
        radiance = data[f'{species.lower()}_tot_rad_kr']
        sigma = radiance/data[f'{species.lower()}_tot_rad_snr']
        
        # Spacecraft position and boresight in MSO
        xyz = np.ndarray((npts, 3))
        bore = np.ndarray((npts, 3))
        corn0 = np.ndarray((npts, 3))
        corn1 = np.ndarray((npts, 3))
        corn2 = np.ndarray((npts, 3))
        corn3 = np.ndarray((npts, 3))
        for i in np.arange(npts):
            xyz[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                               data['planet_sc_vector_tg'][i, :]
                               )/mercury.radius.value
            bore[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                data['boresight_unit_vector_center_tg'][i, :])
            corn0[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c1_tg'][i, :])
            corn1[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c2_tg'][i, :])
            corn2[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c3_tg'][i, :])
            corn3[i, :] = np.dot(data['mso_rotation_matrix'][i, :, :],
                                 data['boresight_unit_vector_c4_tg'][i, :])
        
        xcorner = np.array([corn0[:, 0], corn1[:, 0],
                            corn2[:, 0], corn3[:, 0]]).transpose()
        ycorner = np.array([corn0[:, 1], corn1[:, 1],
                            corn2[:, 1], corn3[:, 1]]).transpose()
        zcorner = np.array([corn0[:, 2], corn1[:, 2],
                            corn2[:, 2], corn3[:, 2]]).transpose()
        
        # Determine tangent point
        t = -np.sum(xyz*bore, axis=1)
        tanpt = xyz+bore*t[:, np.newaxis]
        rtan = np.linalg.norm(tanpt, axis=1)
        
        slit = np.array(['Surface' if s == 0
                         else 'Atmospheric'
                         for s in data['slit']])
        obstype = np.array(
            [str(ob).replace('b', '').replace("'", '').strip()
             for ob in data['obs_typ']])
        
        # Add in the spectra
        spectra = data[species.lower()+'_rad_kr']
        wavelength = data['wavelength']
        raw = data['orig']
        try:
            corrected = data['fully_corr_cr']
        except:
            corrected = data['corr']
        dark = data['dark']
        solarfit = data['sol_fit']
        
        ndata = pd.DataFrame(
            {'species': species,
             'frame': 'MSO',
             'UTC': UTC,
             'orbit': orbit,
             'merc_year': myear,
             'TAA': data['true_anomaly']*np.pi/180.,
             'rmerc': rmerc.value,
             'drdt': data['rad_vel'],
             'subslong': data['subsolar_longitude']*np.pi/180.,
             'g': data['gvals']/u.s,
             'radiance': radiance,
             'sigma': sigma,
             'x': xyz[:, 0]*Rmerc,
             'y': xyz[:, 1]*Rmerc,
             'z': xyz[:, 2]*Rmerc,
             'xbore': bore[:, 0], 'ybore': bore[:, 1], 'zbore': bore[:, 2],
             'xcorn1': xcorner[:, 0], 'xcorn2': xcorner[:, 1],
             'xcorn3': xcorner[:, 2], 'xcorn4': xcorner[:, 3],
             'ycorn1': ycorner[:, 0], 'ycorn2': ycorner[:, 1],
             'ycorn3': ycorner[:, 2], 'ycorn4': ycorner[:, 3],
             'zcorn1': zcorner[:, 0], 'zcorn2': zcorner[:, 1],
             'zcorn3': zcorner[:, 2], 'zcorn4': zcorner[:, 3],
             'obstype': obstype,
             'obstype_num': data['obs_typ_num'],
             'xtan': tanpt[:, 0], 'ytan': tanpt[:, 1],
             'ztan': tanpt[:, 2], 'rtan': rtan,
             'alttan': data['target_altitude_set'][:, 0],
             'minalt': data['minalt'],
             'longtan': data['target_longitude_set'][:, 0]*np.pi/180,
             'lattan': data['target_latitude_set'][:, 0]*np.pi/180,
             'loctimetan': data['obs_solar_localtime'],
             'slit': slit})
        ndata.fillna(-999, inplace=True)
        
        spectra = [spectra[i,:] for i in range(spectra.shape[0])]
        wavelength = [wavelength[i,:] for i in range(wavelength.shape[0])]
        raw = [raw[i,:] for i in range(raw.shape[0])]
        corrected = [corrected[i,:] for i in range(corrected.shape[0])]
        dark = [dark[i,:] for i in range(dark.shape[0])]
        solarfit = [solarfit[i,:] for i in range(solarfit.shape[0])]
        spectra = pd.DataFrame(
            {'spectra': spectra,
             'wavelength': wavelength,
             'raw': raw,
             'corrected': corrected,
             'dark': dark,
             'solarfit': solarfit})
        
        # save this for later
        with open(newfile, 'wb') as f:
            pickle.dump(ndata, f, pickle.HIGHEST_PROTOCOL)
        with open(newfile.replace('.pkl', '_spectra.pkl'), 'wb') as f:
            pickle.dump(spectra, f, pickle.HIGHEST_PROTOCOL)
        
        print('Inserting UVVS data')
        with database_connect() as con:
            print(f'Saving {species} Data')
            for i, dpoint in ndata.iterrows():
                cur.execute(f'''INSERT into {species}uvvsdata (
                                    species, frame, UTC, orbit, merc_year,
                                    taa, rmerc, drdt, subslong, g, radiance,
                                    sigma) values (
                                    '{dpoint.species}',
                                    '{dpoint.frame}',
                                    '{dpoint.UTC.iso}',
                                    {dpoint.orbit},
                                    {dpoint.merc_year},
                                    {dpoint.TAA},
                                    {dpoint.rmerc},
                                    {dpoint.drdt},
                                    {dpoint.subslong},
                                    {dpoint.g},
                                    {dpoint.radiance},
                                    {dpoint.sigma})''')
                cur.execute(f'''INSERT into {species}pointing (
                                    x, y, z, xbore, ybore, zbore,
                                    obstype, obstype_num, xtan, ytan, ztan,
                                    rtan, alttan, longtan, lattan,
                                    loctimetan, slit) values (
                                    {dpoint.x},
                                    {dpoint.y},
                                    {dpoint.z},
                                    {dpoint.xbore},
                                    {dpoint.ybore},
                                    {dpoint.zbore},
                                    '{dpoint.obstype}',
                                    {dpoint.obstype_num},
                                    {dpoint.xtan},
                                    {dpoint.ytan},
                                    {dpoint.ztan},
                                    {dpoint.rtan},
                                    {dpoint.alttan},
                                    {dpoint.longtan},
                                    {dpoint.lattan},
                                    {dpoint.loctimetan},
                                    '{dpoint.slit}')''')
                
            print(f'Saving {species} Spectra')
            for i, spec in spectra.iterrows():
                cur.execute(f'''INSERT into {species}spectra (wavelength,
                                    calibrated, raw, dark, solarfit) values (
                                    %s, %s, %s, %s, %s)''',
                            (spec.wavelength.tolist(), spec.spectra.tolist(),
                             spec.raw.tolist(), spec.dark.tolist(),
                             spec.solarfit.tolist()))

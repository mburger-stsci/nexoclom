import os, os.path
import psycopg2
from . import database


def configfile(setconfig=False):
    """Configure external resources used in the model

    Paths are saved in $HOME/.nexoclom
    * savepath = path where output files are saved
    * database = name of the postgresql database to use
    """

    # Read in current config file if it exists
    configfile = os.path.join(os.environ['HOME'], '.nexoclom')
    config = {}
    if os.path.isfile(configfile):
        for line in open(configfile, 'r').readlines():
            if '=' in line:
                key, value = line.split('=')
                config[key.strip()] = value.strip()
    else:
        setconfig = True

    if setconfig:
        # Prompt user for paths
        if 'savepath' in config:
            oldfile = config['savepath']
            savepath_ = input(f'Path to save files [{oldfile}]: ')
            savepath = oldfile if savepath_ == '' else savepath_
        else:
            savepath = input(f'Path to save model output: ')

        if 'datapath' in config:
            oldfile = config['datapath']
            datapath_ = input(f'Path to data files [{oldfile}]: ')
            datapath = oldfile if datapath_ == '' else datapath_
        else:
            datapath = input(f'Path to data files: ')

        # make sure the database exists
        with psycopg2.connect(host='localhost', database='postgres') as con:
            con.autocommit = True
            cur = con.cursor()
            try:
                cur.execute(f'CREATE database {database}')
            except:
                pass

        # Create save directory if necessary
        if not os.path.isdir(savepath):
            try:
                os.makedir(savepath)
            except:
                assert 0, f'Could not create directory {savepath}'

        # Create the data directory if necessary
        if not os.path.isdir(datapath):
            try:
                 os.makedir(datapath)
            except:
                assert 0, f'Could not create directory {datapath}'

        # Save the config file
        try:
            cfile = open(configfile, 'w')
        except:
            assert 0, f'Could not open {configfile}'

        cfile.write(f'savepath = {savepath}\n')
        cfile.write(f'datapath = {datapath}\n')
        cfile.close()
    else:
        savepath = config['savepath']
        datapath = config['datapath']

    return savepath, datapath


def set_up_output_tables(con):
    cur = con.cursor()

    # drop tables if necessary
    tables = ['outputfile', 'geometry', 'sticking_info', 'forces',
              'spatialdist', 'speeddist', 'angulardist', 'options',
              'modelimages', 'uvvsmodels']
    for tab in tables:
        try:
            cur.execute(f'''DROP table {tab}''')
        except:
            con.rollback()

    # create outputfile table
    cur.execute('''CREATE TABLE outputfile (
                       idnum SERIAL PRIMARY KEY,
                       filename text UNIQUE,
                       npackets bigint,
                       totalsource float,
                       creationtime timestamp NOT NULL)''')
    print('Created outputfile table')

    # create geometry table
    cur.execute('''CREATE TABLE geometry (
                       geo_idnum bigint PRIMARY KEY,
                       planet SSObject,
                       StartPoint SSObject,
                       objects SSObject ARRAY,
                       starttime timestamp,
                       phi real ARRAY,
                       subsolarpt point,
                       TAA float)''')
    print('Created geometry table')

    # Create sticking_info table
    cur.execute('''CREATE TABLE sticking_info (
                       st_idnum bigint PRIMARY KEY,
                       stickcoef float,
                       tsurf float,
                       stickfn text,
                       stick_mapfile text,
                       epsilon float,
                       n float,
                       tmin float,
                       emitfn text,
                       accom_mapfile text,
                       accom_factor float)''')
    print('Created sticking_info table')

    # create forces table
    cur.execute('''CREATE TABLE forces (
                       f_idnum bigint PRIMARY KEY,
                       gravity boolean,
                       radpres boolean)''')
    print('Created forces table')

    # create spatialdist table
    cur.execute('''CREATE TABLE spatialdist (
                       spat_idnum bigint PRIMARY KEY,
                       type text,
                       exobase float,
                       use_map boolean,
                       mapfile text,
                       longitude float[2],
                       latitude float[2])''')
    print('Created spatialdist table')

    # create table speeddist
    cur.execute('''CREATE TABLE speeddist (
                       spd_idnum bigint PRIMARY KEY,
                       type text,
                       vprob float,
                       sigma float,
                       U float,
                       alpha float,
                       beta float,
                       temperature float,
                       delv float)''')
    print('Created speeddist table')

    # create table angulardist
    cur.execute('''CREATE TABLE angulardist (
                       ang_idnum bigint PRIMARY KEY,
                       type text,
                       azimuth float[2],
                       altitude float[2],
                       n float)''')
    print('Created angulardist table')

    ## Skipping perturbvel and plasma_info for now

    # create table options
    cur.execute('''CREATE TABLE options (
                       opt_idnum bigint PRIMARY KEY,
                       endtime float,
                       resolution float,
                       at_once boolean,
                       atom text,
                       lifetime float,
                       fullsystem boolean,
                       outeredge float,
                       motion boolean,
                       streamlines boolean,
                       nsteps int)''')
    print('Created options table')

    # Create table for model images
    cur.execute('''CREATE TABLE modelimages (
                       idnum SERIAL PRIMARY KEY,
                       out_idnum bigint,
                       quantity text,
                       origin text,
                       dims float[2],
                       center float[2],
                       width float[2],
                       subobslongitude float,
                       subobslatitude float,
                       mechanism text,
                       wavelength text,
                       filename text)''')
    print('Created modelimages table')

    # # Create table for MESSENGER comparison
    cur.execute('''CREATE TABLE uvvsmodels (
                       idnum SERIAL PRIMARY KEY,
                       out_idnum bigint,
                       quantity text,
                       orbit int,
                       dphi float,
                       mechanism text,
                       wavelength text,
                       filename text)''')
    print('Created uvvsmodels table')

def configure_model(setcfg=None, setdb=None):
    if setcfg is None:
        cfgfile = input('Reset the model configuration file? (y/n) ')
        setcfg = True if cfgfile.lower() in ('y', 'yes') else False
    else:
        pass
    database, savepath, datapath = configfile(setconfig=setcfg)

    if setdb is None:
        cfgdb = input('Reset the modeloutputs database? (y/n) ')
        setdb = True if cfgfile.lower() in ('y', 'yes') else False
    else:
        pass

    if setdb:
        with psycopg2.connect(host='localhost', database=database) as con:
            con.autocommit = True
            set_up_output_tables(con)

# Installation of Neutral Cloud Model 

### Packages created by Matthew Burger that need to be installed:
1. nexoclom
2. MESSENGERuvvs

------------------

**Installation process**

1. Install Anaconda Python (version >= 3.8):
   1. Download the installer from:
           https://www.anaconda.com/distribution/
   2. double-click the installer to install. This installs anaconda python in
           $HOME/anaconda3 and does not need sysadmin privileges.
   3. Verify it works: Open a new terminal window and start `ipython`. You should
see something like this:
```
(base) [sunra m🍔 /~/]$ ipython
Python 3.8.8 (default, Apr 13 2021, 12:59:45) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.32.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: 
```
NOTE: I think Anaconda python likes the bash shell, but there
are probably ways around that.

2. Create a new python environment with the model.
   1. Obtain the file nexoclom_environment.yml
   2. Create the envirnoment:
           > conda env create -f nexoclom_environment.yml
       (c) To use this environment run:
           > conda activate nexoclom
           None of this will work if the correct environment is not active.
       (d) To turn it off run:
           (nexoclom) > conda deactivate

(3) Create the .nexoclom file
    (a) In your home directory create a file called .nexoclom with the
    following lines:
savepath = <fullpath>/modeloutputs
datapath = <fullpath>/ModelData
database = thesolarsystemmb
mesdatapath = <fullpath>/UVVSData
mesdatabase = messengeruvvsdb

<fullpath> does not need to be the same in all lines, but the directories all
need to be valid.

(3) Initialize the postgres server if necessary:
    (a) In your .bashrc or .bash_profile file (the file that runs when you
    start a terminal window) add the line:
        export PGDATA=/Users/mburger/.postgres/main
    (This step technically isn't needed because the environment variable gets
    set when you activate the environment).
    (b) > initdb -D $PGDATA
    (c) > pg_ctl -l $PGDATA/logfile start
    (d) > createdb <username>
        * Find <username> with > echo $USER
    (e) > createdb thesolarsystemmb
        * This needs to match database in the .nexoclom file
    (g) > createdb messengeruvvsdb
        * This needs to match mesdatabase in the .nexoclom file

(4) Configure the MESSENGER UVVS database
    (a) Download the MESSENGERdata package. Get the link from Matt.
    (b) Put the file in the mesdatapath directory and untar it.
        > tar -xvzf Level1.tar.gz
    (c) Then run:
    (nexoclom)$ ipython
    Python 3.8.13 | packaged by conda-forge | (default, Mar 25 2022, 06:06:49)
    Type 'copyright', 'credits' or 'license' for more information
    IPython 8.2.0 -- An enhanced Interactive Python. Type '?' for help.

    In [1]: from MESSENGERuvvs import initialize_MESSENGERdata

    In [2]: initialize_MESSENGERdata()

    This will take a while to run (hours probably).

(5)

(6) To install updates, run:
    (nexoclom) pip install --upgrade nexoclom
    (nexoclom) pip install --upgrade MESSENGERuvvs

"""Read the model inputs from a file and create the Input object.

The Input object is build from smaller objects defining different model
options.

Geometry
    Defines the Solar System geometry for the Input.

SurfaceInteraction
    Defines the surface interactions.

Forces
    Set which forces act on model particles.

SpatialDist
    Define the initial spatial distribution of particles.

SpeedDist
    Define the initial speed distribution of particles.

AngularDist
    Define the initial angular distribtuion of particles.

Options
    Configure other model parameters
"""
import os
import os.path
import numpy as np
import pandas as pd
from astropy.time import Time
from .Output import Output
from ..utilities import database_connect
from ..utilities.read_configfile import read_configfile
from .input_classes import (Geometry, SurfaceInteraction, Forces, SpatialDist,
                            SpeedDist, AngularDist, Options)
from .ModelImage import ModelImage

class Input:
    def __init__(self, infile):
        """Read the input options from a file.

        **Parameters**
        
        infile
            Plain text file containing model input parameters. See
            :doc:`inputfiles` for a description of the input file format.

        **Class Attributes**

        * geometry
        
        * surface_interaction
        
        * forces
        
        * spatialdist
        
        * speeddist
        
        * angulardist
        
        * options
        
        """
        # Read the configuration file
        config = read_configfile()
        self._savepath = config['savepath']

        # Read in the input file:
        self._inputfile = infile
        params = []
        if os.path.isfile(infile):
            # Remove everything in the line after a comment character
            for line in open(infile, 'r'):
                if ';' in line:
                    line = line[:line.find(';')]
                elif '#' in line:
                    line = line[:line.find('#')]
                else:
                    pass
                    
                if line.count('=') == 1:
                    param_, val_ = line.split('=')
                    if param_.count('.') == 1:
                        sec_, par_ = param_.split('.')
                        params.append((sec_.casefold().strip(),
                                       par_.casefold().strip(),
                                       val_.casefold().strip()))
                    else:
                        pass
                else:
                    pass
        else:
            raise FileNotFoundError(infile)
            
        def extract_param(tag):
            return {b:c for (a,b,c) in params if a == tag}

        self.geometry = Geometry(extract_param('geometry'))
        self.surfaceinteraction = SurfaceInteraction(extract_param(
            'surfaceinteraction'))
        self.forces = Forces(extract_param('forces'))
        self.spatialdist = SpatialDist(extract_param('spatialdist'))
        self.speeddist = SpeedDist(extract_param('speeddist'))
        self.angulardist = AngularDist(extract_param('angulardist'))
        self.options = Options(extract_param('options'))
        
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        result = (self.geometry.__str__() + '\n' +
                  self.surfaceinteraction.__str__() + '\n' +
                  self.forces.__str__() + '\n' +
                  self.spatialdist.__str__() + '\n' +
                  self.speeddist.__str__() + '\n' +
                  self.angulardist.__str__() + '\n' +
                  self.options.__str__())
        
        return result

    def search(self):
        """ Search the database for previous model runs with the same inputs.
        See :doc:`searchtolerances` for tolerances used in searches.
        
        **Parameters**
        
        No parameters.
        
        **Returns**
        
        * A list of filenames corresponding to the inputs.
        
        * Number of packets contained in those saved outputs.
        
        * Total modeled source rate.
        """
        geo_id = self.geometry.search()
        sint_id = self.surfaceinteraction.search()
        for_id = self.forces.search()
        spat_id = self.spatialdist.search()
        spd_id = self.speeddist.search()
        ang_id = self.angulardist.search()
        opt_id = self.options.search()
        
        if None in [geo_id, sint_id, for_id, spat_id, spd_id, ang_id, opt_id]:
            return [], [], 0., 0.
        else:
            query = f'''SELECT idnum, filename, npackets, totalsource
                        FROM outputfile
                        WHERE geo_type = '{self.geometry.type}' and
                              geo_id = {geo_id} and
                              sint_type = '{self.surfaceinteraction.sticktype}' and
                              sint_id = {sint_id} and
                              force_id = {for_id} and
                              spatdist_type = '{self.spatialdist.type}' and
                              spatdist_id = {spat_id} and
                              spddist_type = '{self.speeddist.type}' and
                              spddist_id = {spd_id} and
                              angdist_type = '{self.angulardist.type}' and
                              angdist_id = {ang_id} and
                              opt_id = {opt_id}'''
            with database_connect() as con:
                result = pd.read_sql(query, con)
            
            return (result.idnum.to_list(), result.filename.to_list(),
                    result.npackets.sum(), result.totalsource.sum())

    def run(self, npackets, packs_per_it=None, overwrite=False, compress=True):
        """Run the nexoclom model with the current inputs.
        
        **Parameters**
        
        npackets
            Number of packets to simulate
        
        packs_per_it
            Maximum number of packets to run at one time. Default = 1e5 in
            constant step-size mode; 1e6 in adaptive step-size mode.
        
        overwrite
            Erase any files matching the current inputs that exist.
            Default = False
            
        compress
            Remove packets with frac=0 from the outputs to reduce file size.
            Default = True
            
        **Outputs**
        
        Nothing is returned, but model runs are saved and cataloged.
        """
        t0_ = Time.now()
        print(f'Starting at {t0_}')

        # Determine how many packets have already been run
        if overwrite:
            self.delete_files()
            totalpackets = 0
        else:
            _, outputfiles, totalpackets, _ = self.search()
            print(f'Found {len(outputfiles)} files with {totalpackets} '
                  'packets.')

        npackets = int(npackets)
        ntodo = npackets - totalpackets
        
        if ntodo > 0:
            if (packs_per_it is None) and (self.options.step_size == 0):
                packs_per_it = 1000000
            elif packs_per_it is None:
                packs_per_it = (1e8 * self.options.step_size /
                                self.options.endtime.value)
            else:
                pass
            packs_per_it = int(np.min([ntodo, packs_per_it]))
            
            # Determine how many iterations are needed
            nits = int(np.ceil(ntodo/packs_per_it))
            
            print('Running Model')
            print(f'Will complete {nits} iterations of {packs_per_it} packets.')

            for _ in range(nits):
                tit0_ = Time.now()
                print(f'Starting iteration #{_+1} of {nits}')

                # Create an output object
                Output(self, packs_per_it, compress=compress)
                tit1_ = Time.now()
                print(f'Completed iteration #{_+1} in '
                      f'{(tit1_ - tit0_).sec} seconds.')
        else:
            pass

        t2_ = Time.now()
        dt_ = (t2_-t0_).sec
        if dt_ < 60:
            dt_ = f'{dt_} sec'
        elif dt_ < 3600:
            dt_ = f'{dt_/60} min'
        else:
            dt_ = f'{dt_/3600} hr'
        print(f'Model run completed in {dt_} at {t2_}.')

    def produce_image(self, format_, filenames=None, overwrite=False):
        return ModelImage(self, format_, filenames=filenames,
                          overwrite=overwrite)
    
    def delete_files(self):
        """Delete output files and remove them from the database.

        **Parameters**

        filelist
            List of files to remove. This can be found with Inputs.search()

        **Returns**

        No outputs.

        """
        idnum, filelist, _, _ = self.search()
        with database_connect() as con:
            cur = con.cursor()
            
            for i, f in zip(idnum, filelist):
                # Remove from database and delete file
                print(f'Deleting outputfile {os.path.basename(f)}')
                cur.execute('''DELETE FROM outputfile
                               WHERE idnum = %s''', (i,))
                if os.path.exists(f):
                    os.remove(f)

                # Delete any model images that depend on this output
                cur.execute('''SELECT idnum, filename FROM modelimages
                               WHERE out_idnum = %s''', (i,))
                for mid, mfile in cur.fetchall():
                    print(f'Deleting model image {os.path.basename(mfile)}')
                    cur.execute('''DELETE from modelimages
                                   WHERE idnum = %s''', (mid,))
                    if os.path.exists(mfile):
                        os.remove(mfile)
                
                # Delete any uvvsmodels that depend on this output
                cur.execute('''SELECT idnum, filename FROM uvvsmodels
                               WHERE out_idnum = %s''', (i,))
                for mid, mfile in cur.fetchall():
                    print(f'Deleting uvvsmodel result {os.path.basename(mfile)}')
                    cur.execute('''DELETE from uvvsmodels
                                   WHERE idnum = %s''', (mid,))
                    if os.path.exists(mfile):
                        os.remove(mfile)

                # Remove the fitted uvvsmodels that depend on this output
                uvvsmods = pd.read_sql(
                    f'''SELECT filename FROM uvvsmodels
                        WHERE unfit_idnum={i}''', con)
                for _, uvvsmod in uvvsmods.iterrows():
                    print(f'Deleting fitted uvvsmodel result ' +
                          os.path.basename(uvvsmod.filename))
                    cur.execute(f'''DELETE FROM uvvsmodels
                                    WHERE filename='{uvvsmod.filename}';''')
                    if os.path.exists(uvvsmod.filename):
                        os.remove(uvvsmod.filename)

                # Delete any fitted outputs that depend on this output
                spatdist_id = pd.read_sql(
                    f'''SELECT idnum FROM spatdist_fittedoutput
                        WHERE unfit_outid = {i}''', con)
                for num in spatdist_id.idnum.to_list():
                    # Remove the outputfile
                    fitoutfile = pd.read_sql(
                        f'''SELECT idnum, filename FROM outputfile
                            WHERE spatdist_type = 'fitted output' and
                                  spatdist_id = {num}''', con)

                    for _, row in fitoutfile.iterrows():
                        print(f'Deleting fitted output ' +
                              os.path.basename(row.filename))
                        cur.execute(f'''DELETE FROM outputfile
                                        WHERE filename='{row.filename}';''')
                        if os.path.exists(row.filename):
                            os.remove(row.filename)
                        
                        # Remove the modelimages
                        modims = pd.read_sql(
                            f'''SELECT filename FROM modelimages
                                WHERE out_idnum={row.idnum}''', con)
                        for _, modim in modims.iterrows():
                            print(f'Deleting model image ' +
                                  os.path.basename(modim.filename))
                            cur.execute(f'''DELETE FROM modelimages
                                            WHERE filename='{modim.filename}';''')
                            if os.path.exists(modim.filename):
                                os.remove(modim.filename)
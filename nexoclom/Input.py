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
import sys
import pandas as pd
import logging
from astropy.time import Time
from .Output import Output
from .configure_model import configfile
from .database_connect import database_connect
from .input_classes import (Geometry, SurfaceInteraction, Forces, SpatialDist,
                            SpeedDist, AngularDist, Options)
from .produce_image import ModelImage
from .delete_files import delete_files


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
        
        **Class Methods**

        * findpackets()
        
        * run(npackets, packs_per_it=None, overwrite=False, compress=True)
        
        * produce_image(format, filenames=None)
        
        """
        # Read the configuration file
        self._savepath = configfile()

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

        if extract_param('geometry'):
            self.geometry = Geometry(extract_param('geometry'))
        else:
            assert 0, 'Need to define default action.'
            
        self.surfaceinteraction = SurfaceInteraction(extract_param(
            'surfaceinteraction'))
        
        self.forces = Forces(extract_param('forces'))
        
        self.spatialdist = SpatialDist(extract_param('spatialdist'))
        if self.spatialdist.type == 'surface map':
            if self.spatialdist.coordsystem == 'default':
                if self.geometry.startpoint == self.geometry.planet.object:
                    self.spatialdist.coordsystem = 'solar-fixed'
                else:
                    self.spatialdist.coordsystem = 'planet-fixed'
            else:
                pass
        
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

    def findpackets(self):
        """ Search the database for previous model runs with the same inputs.
        See :doc:`searchtolerances` for tolerances used in searches.
        
        **Parameters**
        
        No parameters.
        
        **Returns**
        
        * A list of filenames corresponding to the inputs.
        
        * Number of packets contained in those saved outputs.
        
        * Total modeled source rate.
        """
        georesult = self.geometry.search(startlist=None)
        if georesult is not None:
            surfintresult = self.surfaceinteraction.search(startlist=georesult)
        else:
            return [], 0, 0

        if surfintresult is not None:
            forceresult = self.forces.search(startlist=surfintresult)
        else:
            return [], 0, 0

        if forceresult is not None:
            spatresult = self.spatialdist.search(startlist=forceresult)
        else:
            return [], 0, 0

        if spatresult is not None:
            spdresult = self.speeddist.search(startlist=spatresult)
        else:
            return [], 0, 0

        if spdresult is not None:
            angresult = self.angulardist.search(startlist=spdresult)
        else:
            return [], 0, 0

        if angresult is not None:
            finalresult = self.options.search(startlist=angresult)
        else:
            return [], 0, 0

        if finalresult is not None:
            result_ = [str(s) for s in finalresult]
            resultstr = f"({', '.join(result_)})"
            with database_connect() as con:
                result = pd.read_sql(
                    f'''SELECT filename, npackets, totalsource
                        FROM outputfile
                        WHERE idnum in {resultstr}''', con)
            npackets = result.npackets.sum()
            totalsource = result.totalsource.sum()

            return result.filename.to_list(), npackets, totalsource
        else:
            return [], 0, 0

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
        # Configure the logger
        # Note: The logfile name will be changed to match the outputfile name
        #       and stored in the Output object.
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        log_file_handler = logging.FileHandler('log.out', 'w')
        logger.addHandler(log_file_handler)
        out_handler = logging.StreamHandler(sys.stdout)
        logger.addHandler(out_handler)
        fmt = logging.Formatter('%(levelname)s: %(msg)s')
        log_file_handler.setFormatter(fmt)
        out_handler.setFormatter(fmt)

        t0_ = Time.now()
        logger.info(f'Starting at {t0_}')
        
        if len(self.geometry.planet) != 1:
            logger.error('Gravity and impact check not working for '
                          'planets with moons.')
            sys.exit()
            
        # Determine how many packets have already been run
        # outputfiles, totalpackets, _ = self.findpackets()
        # logger.info(f'Found {len(outputfiles)} files with {totalpackets} '
        #              'packets.')
        #
        # if (overwrite) and (totalpackets > 0):
        #     # delete files and remove from database
        #     delete_files(outputfiles)
        #     totalpackets = 0
        # else:
        #     pass
        totalpackets = 0
        
        npackets = int(npackets)
        ntodo = npackets - totalpackets
        
        if ntodo > 0:
            if packs_per_it is None:
                packs_per_it = (100000
                                if self.options.step_size > 0
                                else int(1e6))
            else:
                pass
            packs_per_it = min(ntodo, packs_per_it)
            
            # Determine how many iterations are needed
            nits = ntodo//packs_per_it + 1
            
            logger.info('Running Model')
            logger.info(f'Will complete {nits} iterations of {packs_per_it} '
                         'packets.')
            
            for _ in range(nits):
                tit0_ = Time.now()
                logger.info(f'Starting iteration #{_+1} of {nits}')
                
                # Create an output object
                Output(self, packs_per_it, compress=compress,
                       logger=logger)
                # Just run and save the model when output is created
                # No reason to explicitly call run
                
                tit1_ = Time.now()
                logger.info(f'Completed iteration #{_+1} in '
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
        logger.info(f'Model run completed in {dt_} at {t2_}.')
        out_handler.close()
        
    def produce_image(self, format_, filenames=None):
        return ModelImage(self, format_, filenames=filenames)

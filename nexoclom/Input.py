"""Read the model inputs from a file and create the Input object.

The Input object is build from smaller objects defining different model
options.

Geometry
    Defines the Solar System geometry for the Input.

StickingInfo
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
import pandas as pd
from .configure_model import configfile
from .database_connect import database_connect
from .input_classes import (Geometry, StickingInfo, Forces, SpatialDist,
                            SpeedDist, AngularDist, Options)
from .modeldriver import modeldriver
from .produce_image import ModelImage


class Input:
    def __init__(self, infile):
        """Read the input options from a file.

        **Parameters**
        
        infile
            Plain text file containing model input parameters. See
            :doc:`inputfiles` for a description of the input file format.

        **Class Attributes**

        * geometry
        
        * sticking_info
        
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

        print(extract_param('geometry'))
        self.geometry = Geometry(extract_param('geometry'))
        self.sticking_info = StickingInfo(extract_param('sticking_info'))
        self.forces = Forces(extract_param('forces'))
        self.spatialdist = SpatialDist(extract_param('spatialdist'))
        self.speeddist = SpeedDist(extract_param('speeddist'))
        self.angulardist = AngularDist(extract_param('angulardist'),
                                       self.spatialdist)
        self.options = Options(extract_param('options'), self.geometry.planet)
        
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        result = (self.geometry.__str__() + '\n' +
                  self.sticking_info.__str__() + '\n' +
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
            stickresult = self.sticking_info.search(startlist=georesult)
        else:
            return [], 0, 0

        if stickresult is not None:
            forceresult = self.forces.search(startlist=stickresult)
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
        
        modeldriver(self, npackets, packs_per_it=packs_per_it,
                    overwrite=overwrite, compress=compress)

    def produce_image(self, format, filenames=None):
        return ModelImage(self, format, filenames=filenames)

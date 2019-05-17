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
import os, os.path
import numpy as np
import psycopg2
import pandas as pd
from astropy.io import ascii
import astropy.units as u
from .produce_image import ModelImage
from .modeldriver import modeldriver
from .configure_model import configfile
from .input_classes import (Geometry, StickingInfo, Forces,
                            SpatialDist, SpeedDist, AngularDist, Options)


class Input():
    def __init__(self, infile):
        """Read the input options from a file.

        **Parameters**
        infile
            Plain text file containing model input parameters.

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
        * run(npackets, overwrite=False, compress=True)
        * produce_image(format_, filenames=None)
        """
        # Read the configuration file
        self.__savepath, self.__database = configfile()

        # Read in the input file:
        self.__inputfile = infile
        if os.path.isfile(infile):
            data = ascii.read(infile, delimiter='=', comment=';',
                              data_start=0, names=['Param', 'Value'])
        else:
            assert 0, 'Input file {} does not exist.'.format(infile)

        section = [d.split('.')[0].casefold() for d in data['Param']]
        param = [d.split('.')[1].casefold() for d in data['Param']]
        values = [v.split(';')[0].strip().casefold()
                  if ';' in v else v.casefold() for v in data['Value']]

        # Extract the geometry parameters
        gparam = {b:c for (a,b,c) in zip(section, param, values)
                  if a == 'geometry'}
        self.geometry = Geometry(gparam)

        # Extract the sticking information
        sparam = {b:c for a,b,c in zip(section, param, values)
                  if a == 'sticking_info'}
        self.sticking_info = StickingInfo(sparam)

        # Extract the forces
        fparam = {b:c for a,b,c in zip(section, param, values)
                  if a == 'forces'}
        self.forces = Forces(fparam)

        # Extract the spatial distribution
        sparam = {b:c for a,b,c in zip(section, param, values)
                  if a == 'spatialdist'}
        self.spatialdist = SpatialDist(sparam)

        # Extract the speed distribution
        sparam = {b:c for a,b,c in zip(section, param, values)
                  if a == 'speeddist'}
        self.speeddist = SpeedDist(sparam)

        # Extract the angular distribution
        aparam = {b:c for a,b,c in zip(section, param, values)
                  if a == 'angulardist'}
        self.angulardist = AngularDist(aparam, self.spatialdist)

        # Extract the options
        oparam = {b:c for a,b,c in zip(section, param, values)
                  if a == 'options'}
        self.options = Options(oparam, self.geometry.planet)

    def __str__(self):
        print(self.geometry)
        print(self.sticking_info)
        print(self.forces)
        print(self.spatialdist)
        print(self.speeddist)
        print(self.angulardist)
        print(self.options)
        return ''

    def findpackets(self):
        '''
        Search the database for identical inputs
        '''
        georesult = self.geometry.search(self.__database, startlist=None)
        stickresult = self.sticking_info.search(self.__database,
                                                startlist=georesult)
        forceresult = self.forces.search(self.__database,
                                         startlist=stickresult)
        spatresult = self.spatialdist.search(self.__database,
                                             startlist=forceresult)
        spdresult = self.speeddist.search(self.__database,
                                          startlist=spatresult)
        angresult = self.angulardist.search(self.__database,
                                            startlist=spdresult)
        finalresult = self.options.search(self.__database,
                                          startlist=angresult)

        if finalresult is None:
            return [], 0, 0
        else:
            result_ = [str(s) for s in finalresult]
            resultstr = f"({', '.join(result_)})"
            with psycopg2.connect(database=self.__database) as con:
                result = pd.read_sql(
                    f'''SELECT filename, npackets, totalsource
                        FROM outputfile
                        WHERE idnum in {resultstr}''', con)
            npackets = result.npackets.sum()
            totalsource = result.totalsource.sum()

            return result.filename.to_list(), npackets, totalsource

    def run(self, npackets, overwrite=False, compress=True):
        modeldriver(self, npackets, overwrite, compress)

    def produce_image(self, format_, filenames=None):
        return ModelImage(self, format_, filenames=filenames)

import numpy as np
import astropy.units as u
from scipy.spatial import KDTree

from nexoclom.math import rotation_matrix, Histogram2d
from nexoclom.data_simulation.ModelResult import ModelResult
from nexoclom.particle_tracking.Output import Output
from nexoclom import __path__ as basepath

import bokeh.plotting as bkp
from bokeh.palettes import Inferno256
from bokeh.models import (HoverTool, ColumnDataSource, ColorBar,
                          LogColorMapper, LogTicker, LinearColorMapper)
from bokeh.io import curdoc, export_png
from bokeh.themes import Theme


class ModelDensity:
    def __init__(self, inputs, xpts, ypts, zpts, dr=0.05):
        """ Create Images from model results.
        This Assumes the model has already been run.
        
        Parameters
        ==========
        inputs
            An Input object
        
        params
            A dictionary with format information or a path to a formatfile.
            
        filenames
            A filename or list of filenames to use. Default = None is to
            find all files created for the inputs.
            
        overwrite
            If True, deletes any images that have already been computed.
            Default = False
            
        Notes
        =====
        Required parameters
            xpts, ypts, zpts
        
        Optional parameters
            dr
        
        """
        self.type = 'density'
        
        self.origin = inputs.geometry.planet
        self.unit = u.def_unit('R_' + self.origin.object, self.origin.radius)
        
        self.density = np.zeros(len(xpts))
        self.packets = np.zeros(len(xpts))
        self.totalsource = 0
        
        self.dr = dr*self.unit
        self.Vpix = (4/3/np.pi*self.dr**3).to(u.cm**3)

        self.outid, self.outputfiles, _, _ = inputs.search()
        # self.outputfiles = [self.outputfiles[0]]

        for fname in self.outputfiles:
            print(f'Output filename: {fname}')
            output = Output.restore(fname)
            
            packets = output.X
            data_pts = packets[['x', 'y', 'z']].values
            tree = KDTree(data_pts)
            
            sample_pts = np.stack([xpts, ypts, zpts]).to(self.unit).transpose().value
            pts = tree.query_ball_point(sample_pts, self.dr.value)
            mapped = list(map(lambda x: (packets.iloc[x].frac.sum(), len(x)), pts))
            density = np.array([m[0] for m in mapped])
            packets = np.array([m[1] for m in mapped])
            
            # image_, packets_, totalsource_ = image_step(self, fname, overwrite)
            self.density += density
            self.packets += packets
            self.totalsource += output.totalsource

        mod_rate = self.totalsource / output.inputs.options.endtime.value
        self.atoms_per_packet = 1e23 / mod_rate
        self.sourcerate = 1e23 / u.s
        self.density = self.density * self.atoms_per_packet/self.Vpix

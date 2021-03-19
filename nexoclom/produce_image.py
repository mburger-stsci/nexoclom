import os.path
import numpy as np
import numpy.matlib as npmat
import pandas as pd
import pickle
import astropy.units as u
from solarsystemMB import SSObject
from mathMB import rotation_matrix
from .ModelResults import ModelResult
from .database_connect import database_connect
from .input_classes import InputError
from .Output import Output


class ModelImage(ModelResult):
    def __init__(self, inputs, format, filenames=None, overwrite=False):
        """ Create Images from model results.
        This Assumes the model has already been run.
        
        Parameters
        ==========
        inputs
            An Input object
        
        format
            A dictionary with format information or a path to a formatfile.
            
        filenames
            A filename or list of filenames to use. Default = None is to
            find all files created for the inputs.
            
        overwrite
            If True, deletes any images that have already been computed.
            Default = False
        """
        super().__init__(inputs, format, filenames=filenames)
        self.type = 'image'
        
        if 'origin' in self.format:
            self.origin = SSObject(self.format['origin'])
        else:
            self.origin = inputs.geometry.planet
            
        self.unit = u.def_unit('R_' + self.origin.object,
                               self.origin.radius)

        if 'dims' in self.format:
            dimtemp = self.format['dims'].split(',')
            self.dims = [int(dimtemp[0]), int(dimtemp[1])]
        else:
            self.dims = [800, 800]

        if 'center' in self.format:
            centtemp = self.format['center'].split(',')
            self.center = [float(centtemp[0])*self.unit,
                           float(centtemp[1])*self.unit]
        else:
            self.center = [0*self.unit, 0*self.unit]

        if 'width' in self.format:
            widtemp = self.format['width'].split(',')
            self.width = [float(widtemp[0])*self.unit,
                          float(widtemp[1])*self.unit]
        else:
            self.width = [8*self.unit, 8*self.unit]

        if 'subobslongitude' in self.format:
            self.subobslongitude = float(self.format['subobslongitude'])*u.rad
        else:
            self.subobslongitude = 0.*u.rad

        if 'subobslongitude' in self.format:
            self.subobslatitude = float(self.format['subobslatitude'])*u.rad
        else:
            self.subobslatitude = np.pi*u.rad

        self.image = np.zeros(self.dims)
        self.packet_image = np.zeros(self.dims)
        self.blimits = None
        immin = tuple(c - w/2 for c, w in zip(self.center, self.width))
        immax = tuple(c + w/2 for c, w in zip(self.center, self.width))
        scale = tuple(w/(d-1) for w, d in zip(self.width, self.dims))
        self.Apix = (scale[0]*scale[1]).to(u.cm**2)

        self.xaxis = np.linspace(immin[0], immax[0], self.dims[0])
        self.zaxis = np.linspace(immin[1], immax[1], self.dims[1])

        for i, fname in enumerate(self.filenames):
            # Search to see if its already been done
            print(f'Output filename: {fname}')
            image_, packets_ = self.restore(fname, overwrite=overwrite)

            if image_ is None:
                image_, packets_ = self.create_image(fname)
                print(f'Completed image {i+1} of {len(self.filenames)}')
            else:
                print(f'Image {i+1} of {len(self.filenames)} '
                       'previously completed.')

            self.image += image_
            self.packet_image += packets_

        self.image *= self.atoms_per_packet

    def save(self, fname, image, packets):
        # Determine the id of the outputfile
        idnum = int(os.path.basename(fname).split('.')[0])

        # Insert the image into the database
        if self.quantity == 'radiance':
            mech = ', '.join(sorted([m for m in self.mechanism]))
            wave_ = sorted([w.value for w in self.wavelength])
            wave = ', '.join([str(w) for w in wave_])
        else:
            mech = None
            wave = None

        width = [w.value for w in self.width]
        center = [c.value for c in self.center]
        
        with database_connect() as con:
            cur = con.cursor()
            cur.execute(f'''INSERT into modelimages (out_idnum, quantity,
                                origin, dims, center, width, subobslongitude,
                                subobslatitude, mechanism, wavelength,
                                filename)
                            VALUES (%s, %s, %s, %s::INT[2],
                                    %s::DOUBLE PRECISION[2],
                                    %s::DOUBLE PRECISION[2],
                                    %s, %s, %s, %s, 'temp')''',
                        (idnum, self.quantity, self.origin.object,
                         self.dims, center, width,
                         self.subobslongitude.value, self.subobslatitude.value,
                         mech, wave))

        idnum_ = pd.read_sql(f'''SELECT idnum
                                FROM modelimages
                                WHERE filename = 'temp';''', con)
        assert len(idnum_) == 1
        idnum = int(idnum_.idnum[0])

        savefile = os.path.join(os.path.dirname(fname), f'image.{idnum}.pkl')
        with open(savefile, 'wb') as f:
            pickle.dump((image, packets), f)
        cur.execute(f'''UPDATE modelimages
                        SET filename=%s
                        WHERE idnum = %s''', (savefile, idnum))
        con.close()

    def restore(self, fname, overwrite=False):
        # Determine the id of the outputfile
        if self.quantity == 'radiance':
            mech = ("mechanism = '" +
                    ", ".join(sorted([m for m in self.mechanism])) +
                    "'")
            wave_ = sorted([w.value for w in self.wavelength])
            wave = ("wavelength = '" +
                    ", ".join([str(w) for w in wave_]) +
                    "'")
        else:
            mech = 'mechanism is NULL'
            wave = 'wavelength is NULL'

        with database_connect() as con:
            idnum_ = pd.read_sql(f'''SELECT idnum
                                FROM outputfile
                                WHERE filename='{fname}' ''', con)
            oid = idnum_.idnum[0]

            result = pd.read_sql(
                f'''SELECT filename FROM modelimages
                    WHERE out_idnum = {oid} and
                          quantity = '{self.quantity}' and
                          origin = '{self.origin.object}' and
                          dims[1] = {self.dims[0]} and
                          dims[2] = {self.dims[1]} and
                          center[1] = {self.center[0].value} and
                          center[2] = {self.center[1].value} and
                          width[1] = {self.width[0].value} and
                          width[2] = {self.width[1].value} and
                          subobslongitude = {self.subobslongitude.value} and
                          subobslatitude = {self.subobslatitude.value} and
                          {mech} and
                          {wave}''', con)

        if (len(result) == 1) and overwrite:
            if os.path.exists(result.filename[0]):
                os.remove(result.filename[0])
            with database_connect() as con:
                cur = con.cursor()
                cur.execute('''DELETE FROM modelimages
                               WHERE filename = %s''', (result.filename[0],))
            image, packets = None, None
        elif len(result) == 1:
            image, packets = pickle.load(open(result.filename[0], 'rb'))
        elif len(result) == 0:
            image, packets = None, None
        else:
            raise RuntimeError('ModelImage.restore',
                               'Should not be able to get here.')

        return image, packets

    def create_image(self, fname):
        # Determine the proper frame rotation
        M = self.image_rotation()

        # Load data in solar reference frame
        output = Output.restore(fname)
        packets = output.X
        if self.origin != self.inputs.geometry.planet:
            super().transform_reference_frame(packets)
        
        packets['radvel_sun'] = (packets['vy'] +
                                 output.vrplanet.to(self.unit/u.s).value)

        # packet positions in an array
        pts_sun = packets[['x', 'y', 'z']].values

        # Rotate to observer frame
        pts_obs = np.array(np.matmul(M, pts_sun.transpose()).transpose())

        # Determine which packets are not blocked by planet
        rhosqr_obs = np.linalg.norm(pts_obs[:, [0, 2]], axis=1)
        inview = (rhosqr_obs > 1) | (pts_obs[:,1] < 0)
        packets['frac'] *= inview

        # Which packets are in sunlight
        rhosqr_sun = np.linalg.norm(pts_sun[:, [0, 2]], axis=1)
        out_of_shadow = (rhosqr_sun > 1) | (pts_sun[:,1] < 0)

        # Packet weighting
        self.packet_weighting(packets, out_of_shadow, output.aplanet)
        packets['weight'] /= self.Apix

        # Edges of each pixel (bin) for the histogram
        dx = (self.xaxis[1]-self.xaxis[0])/2.
        bx = np.append(self.xaxis-dx, self.xaxis[-1]+dx)
        dz = (self.zaxis[1]-self.zaxis[0])/2.
        bz = np.append(self.zaxis-dz, self.zaxis[-1]+dz)

        pts_obs = pts_obs.transpose()
        image, _, _ = np.histogram2d(pts_obs[0,:], pts_obs[2,:],
                                     weights=packets['weight'], bins=(bx, bz))
        packim, _, _ = np.histogram2d(pts_obs[0,:], pts_obs[2,:],
                                       bins=(bx, bz))

        self.save(fname, image, packim)
        del output

        return image, packim

    def display(self, savefile='image.png', limits=None, show=True):
        import bokeh.plotting as bkp
        from bokeh.palettes import Inferno256
        from bokeh.models import (HoverTool, ColumnDataSource, ColorBar,
                                  LogColorMapper, LogTicker)
        from bokeh.io import curdoc
        from bokeh.themes import Theme

        if self.unit.__str__() == 'R_Mercury':
            ustr = 'R_M'
        else:
            ustr = 'R_obj'
            
        if self.quantity == 'radiance':
            runit = 'kR'
            rname = 'Radiance'
        elif self.quantity == 'column':
            runit = 'cm-2'
            rname = 'Column'
        else:
            assert 0

        tooltips = [('x', '$x{0.1f} ' + ustr),
                    ('y', '$y{0.1f} ' + ustr),
                    (rname, '@image ' + runit)]
        
        curdoc().theme = Theme(os.path.join(os.path.dirname(__file__),
                                            'data', 'bokeh.yml'))

        low = np.min(self.image[self.image > 0])
        color_mapper = LogColorMapper(palette=Inferno256, low=low,
                                      high=np.max(self.image))

        x0 = np.min(self.xaxis.value)
        y0 = np.min(self.zaxis.value)
        dw = np.max(self.xaxis.value) - np.min(self.xaxis.value)
        dh = np.max(self.zaxis.value) - np.min(self.zaxis.value)

        fig = bkp.figure(plot_width=1000, plot_height=1000,
                         title=f'{self.inputs.options.species} {rname}',
                         x_axis_label=f'Distance ({ustr})',
                         y_axis_label=f'Distance ({ustr})',
                         x_range=[np.min(self.xaxis.value),
                                  np.max(self.xaxis.value)],
                         y_range=[np.min(self.zaxis.value),
                                  np.max(self.zaxis.value)],
                         tooltips=tooltips)

        fig.image(image=[self.image], x=x0, y=y0, dw=dw, dh=dh,
                  color_mapper=color_mapper)
        xc = np.cos(np.linspace(0, 2*np.pi, 1000))
        yc = np.sin(np.linspace(0, 2*np.pi, 1000))
        fig.patch(xc, yc, fill_color='yellow')
        
        bkp.output_file(savefile)
        if show:
            bkp.show(fig)
        else:
            bkp.save(fig)
            
        assert 0

        # from astropy.visualization import (PercentileInterval, LogStretch,
        #                                    ImageNormalize)
        #
        # # Determine limits if none given
        # if limits is None:
        #     interval = PercentileInterval(95)
        #     self.blimits = interval.get_limits(self.image[self.image > 0])
        # elif len(limits) == 2:
        #     self.blimits = limits
        # else:
        #     assert 0, 'Problem with the display limits'

#        norm = ImageNormalize(self.image, stretch=LogStretch(),
#                              vmin=self.blimits[0], vmax=self.blimits[1])

        # # Make the colorbar
        # if self.quantity == 'column':
        #     clabel = f'$N_{{ {self.inputs.options.species} }}\ cm^{{-2}}$'
        # else:
        #     clabel = f'$I_{{ {self.inputs.options.species} }} R$'
        # cbar = fig.colorbar(im, shrink=0.7, label=clabel)

        # Put Planet's disk in the middle
        # xc, yc = (np.cos(np.linspace(0, 2*np.pi, 1000)),
        #           np.sin(np.linspace(0, 2*np.pi, 1000)))
        # ax.fill(xc, yc, 'y')

        return fig
    
    def image_rotation(image):
        slong = image.subobslongitude
        slat = image.subobslatitude
        
        pSun = np.array([0., -1., 0.])
        pObs = np.array([np.sin(slong)*np.cos(slat),
                         -np.cos(slong)*np.cos(slat),
                         np.sin(slat)])
        if np.array_equal(pSun, pObs):
            M = npmat.identity(3)
        else:
            costh = np.dot(pSun, pObs)/np.linalg.norm(pSun)/np.linalg.norm(pObs)
            theta = np.arccos(np.clip(costh, -1, 1))
            axis = np.cross(pSun, pObs)
            M = rotation_matrix(theta, axis)
        
        #M = np.transpose(M)
        return M



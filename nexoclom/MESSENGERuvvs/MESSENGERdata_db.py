"""MESSENGER UVVS data class"""
import numpy as np
import pandas as pd
import pickle
import copy
from astropy import units as u

from nexoclom import Input, LOSResult
import mathMB
from .database_setup import database_connect
from .plot_methods import plot_bokeh, plot_plotly, plot_fitted


class InputError(Exception):
    """Raised when a required parameter is not included."""
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class MESSENGERdata:
    """Retrieve MESSENGER data from database.
    Given a species and set of comparisons, retrieve MESSSENGER UVVS
    data from the database. The list of searchable fields is given at
    :doc:`database_fields`.
    
    Returns a MESSENGERdata Object.
    
    **Parameters**
    
    species
        Species to search. This is required because the data from each
        species is stored in a different database table.
        
    query
        A SQL-style list of comparisons.
        
    load_spectra
        Load the individual spectra as well as the radiance data.
        Default = False
        
    The data in the object created is extracted from the database tables usingu
    the query:
    
    ::
    
        SELECT *
        FROM <species>uvvsdata, <species>pointing
        WHERE <query>
    
    See examples below.
    
    **Class Atributes**
    
    species
        The object can only contain a single species.
        
    frame
        Coordinate frame for the data, either MSO or Model.
        
    query
        SQL query used to search the database and create the object.
    
    data
        Pandas dataframe containing result of SQL query. Columns in the
        dataframe are the same as in the database except *frame* and
        *species* have been dropped as they are redundant. If models have been
        run, there are also columns in the form modelN for the Nth model run.
        
    taa
        Median true anomaly for the data in radians.
        
    model_label
        If *N* models have been run, this is a dictionary in the form
        `{'model0':label0, ..., 'modelN':labelN}` containing descriptions for
        the models.
        
    model_strength
        If *N* models have been run, this is a dictionary in the form
        `{'model0':strength0, ..., 'modelN':strengthN}` containing modeled
        source rates in units of :math:`10^{23}` atoms/s.
        
    **Examples**
    
    1. Loading data
    
    ::
    
        >>> from MESSENGERuvvs import MESSENGERdata
        
        >>> CaData = MESSENGERdata('Ca', 'orbit = 36')
        
        >>> print(CaData)
        Species: Ca
        Query: orbit = 36
        Frame: MSO
        Object contains 581 spectra.
        
        >>> NaData = MESSENGERdata('Na', 'orbit > 100 and orbit < 110')
        
        >>> print(NaData)
        Species: Na
        Query: orbit > 100 and orbit < 110
        Frame: MSO
        Object contains 3051 spectra.
        
        >>> MgData = MESSENGERdata('Mg',
                'loctimetan > 5.5 and loctimetan < 6.5 and alttan < 1000')
        
        >>> print(len(MgData))
        45766
        
    2. Accessing data.
    
    * The observations are stored within the MESSENGERdata object in a
      `pandas <https://pandas.pydata.org>`_ dataframe attribute called *data*.
      Please see the `pandas documentation <https://pandas.pydata.org>`_ for
      more information on how to work with dataframes.
    
    ::
    
        >>> print(CaData.data.head(5))
                                 utc  orbit  merc_year  ...  loctimetan         slit               utcstr
        unum                                            ...
        3329 2011-04-04 21:24:11.820     36          0  ...   14.661961  Atmospheric  2011-04-04T21:24:11
        3330 2011-04-04 21:25:08.820     36          0  ...   12.952645  Atmospheric  2011-04-04T21:25:08
        3331 2011-04-04 21:26:05.820     36          0  ...   12.015670  Atmospheric  2011-04-04T21:26:05
        3332 2011-04-04 21:27:02.820     36          0  ...   12.007919  Atmospheric  2011-04-04T21:27:02
        3333 2011-04-04 21:27:59.820     36          0  ...   12.008750  Atmospheric  2011-04-04T21:27:59
        
        [5 rows x 29 columns]

    * Individual observations can be extracted using standard Python
      slicing techniques:
     
    ::
        
        >>> print(CaData[3:8])
        Species: Ca
        Query: orbit = 36
        Frame: MSO
        Object contains 5 spectra.

        >>> print(CaData[3:8].data['taa'])
        unum
        3332    1.808107
        3333    1.808152
        3334    1.808198
        3335    1.808243
        3336    1.808290
        Name: taa, dtype: float64

    3. Modeling data
    
    ::
    
        >>> inputs = Input('Ca.spot.Maxwellian.input')
        >>> CaData.simulate_data_from_inputs(inputs, 1e5, label='Model 1')
        >>> inputs.speeddist.temperature /= 2.  # Run model with different temperature
        >>> CaData.simulate_data_from_inputs(inputs, 1e5, label='Model 2')
        
    4. Plotting data
    
    ::
    
        >>> CaData.plot('Ca.orbit36.models.html')
    
    5. Exporting data to a file
    
    ::
    
        >>> CaData.export('modelresults.csv')
        >>> CaData.export('modelresults.html', columns=['taa'])

    
    """
    def __init__(self, species, comparisons=None, load_spectra=False):
        allspecies = ['Na', 'Ca', 'Mg']
        if species in allspecies:
            self.species = species
        else:
            raise InputError('MESSENGERdata.__init__', f'{species} not valid.')
        
        self.frame = None
        self.query = None
        self.data = None
        self.taa = None
        self.app = None  # Plotly app
        self.model_info = {}
        
        if comparisons is None:
            # Return list of queryable fields
            with database_connect() as con:
                columns = pd.read_sql(
                    f'''SELECT * from {species}uvvsdata, {species}pointing
                        WHERE 1=2''', con)
            print('Available fields are:')
            for col in columns.columns:
                print(f'\t{col}')
        else:
            # Run the query and try to make the object
            query = f'''SELECT * from {species}uvvsdata, {species}pointing
                        WHERE unum=pnum and ({comparisons})
                        ORDER BY unum'''
            try:
                with database_connect() as con:
                    data = pd.read_sql(query, con)
            except Exception:
                print('Requested data:')
                print(query)
                raise InputError('MESSENGERdata.__init__',
                                 'Problem with comparisons given.')

            if len(data) > 0:
                self.frame = data.frame[0]
                self.query = self._format_query(comparisons)
                self.taa = np.median(data.taa)*u.rad
                
                data.drop(['species', 'frame'], inplace=True, axis=1)
                data.loc[data.alttan < 0, 'alttan'] = 1e10
                data.set_index('unum', inplace=True)
                
                if load_spectra:
                    specquery = f'''SELECT *
                                    FROM {species}spectra
                                    WHERE snum in {data.index.to_list()}'''
                    specquery = specquery.replace('[', '(').replace(']', ')')
                    with database_connect() as con:
                        spectra = pd.read_sql(specquery, con)
                        
                    spectra.set_index('snum', inplace=True)
                    for r, row in spectra.iterrows():
                        spectra.loc[r, 'wavelength'] = np.array(row.wavelength)
                        spectra.loc[r, 'raw'] = np.array(row.raw)
                        spectra.loc[r, 'solarfit'] = np.array(row.solarfit)
                        spectra.loc[r, 'dark'] = np.array(row.dark)
                        spectra.loc[r, 'calibrated'] = np.array(row.calibrated)
                    self.data = data.merge(spectra, left_index=True,
                                           right_index=True)
                else:
                    self.data = data
            else:
                print(query)
                print('No data found')
                
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        result = ('MESSENGER UVVS Data Object\n'
                  f'Species: {self.species}\n'
                  f'Query: {self.query}\n'
                  f'Frame: {self.frame}\n'
                  f'Object contains {len(self)} spectra.')
        return result

    def __len__(self):
        return len(self.data)

    def __getitem__(self, q_):
        if isinstance(q_, int):
            q = slice(q_, q_+1)
        elif isinstance(q_, slice):
            q = q_
        elif isinstance(q_, pd.Series):
            q = np.where(q_)[0]
        else:
            raise TypeError

        new = copy.deepcopy(self)
        new.query = self.query + ' [Subset]'
        new.taa = self.taa
        new.data = self.data.iloc[q].copy()
        new.model_info = self.model_info

        return new

    def __iter__(self):
        for i in range(len(self.data)):
            yield self[i]

    def keys(self):
        """Return all keys in the object, including dataframe columns"""
        keys = list(self.__dict__.keys())
        keys.extend([f'data.{col}' for col in self.data.columns])
        return keys

    @staticmethod
    def _format_query(comparison):
        """Try to impose some regularity on the query."""
        query = comparison.lower()
    
        # Remove some extra spaces
        chars = ['=', '>', '<', '!=', '>=', '<=']
        for char in chars:
            query = query.replace(f'{char} ', char)
            query = query.replace(f' {char}', char)

        return query

    def set_frame(self, frame=None):
        """Convert between MSO and Model frames.

        More frames could be added if necessary.
        If Frame is not specified, flips between MSO and Model."""
        if (frame is None) and (self.frame.casefold() == 'MSO'.casefold()):
            frame = 'Model'
        elif (frame is None) and (self.frame.casefold() == 'MODEL'.casefold()):
            frame = 'MSO'
        else:
            pass

        allframes = ['MODEL'.casefold(), 'MSO'.casefold()]
        if frame.casefold() not in allframes:
            print('{} is not a valid frame.'.format(frame))
        elif frame == self.frame:
            pass
        elif ((self.frame.casefold() == 'MSO'.casefold()) and
              (frame.casefold() == 'MODEL'.casefold())):
            # Convert from MSO to Model
            self.data.x, self.data.y = self.data.y.copy(), -self.data.x.copy()
            self.data.xbore, self.data.ybore = (self.data.ybore.copy(),
                                                -self.data.xbore.copy())
            self.data.xtan, self.data.ytan = (self.data.ytan.copy(),
                                              -self.data.xtan.copy())
            self.frame = 'Model'
        elif ((self.frame.casefold() == 'MODEL'.casefold()) and
              (frame.casefold() == 'MODEL'.casefold())):
            self.data.x, self.data.y = -self.data.y.copy(), self.data.x.copy()
            self.data.xbore, self.data.ybore = (-self.data.ybore.copy(),
                                                self.data.xbore.copy())
            self.data.xtan, self.data.ytan = (-self.data.ytan.copy(),
                                              self.data.xtan.copy())
            self.frame = 'MSO'
        else:
            assert 0, 'You somehow picked a bad combination.'
            

    def compute_radiance_from_map(self, unfitnum, sourcemap_, label=None,
                                  save_packets=False):
        oldkey = f'model{unfitnum:02d}'
        oldinfo = self.model_info[oldkey]
        
        if isinstance(sourcemap_, str) and sourcemap_.endswith('.pkl'):
            with open(sourcemap_, 'rb') as mfile:
                sourcemap = pickle.load(mfile)
        elif isinstance(sourcemap_, str) and sourcemap_.endswith('.sav'):
            from scipy.io import readsav
            sourcemap_ = readsav(sourcemap_)
            sourcemap = {'longitude':sourcemap_['longitude'][0]*u.rad,
                         'latitude':sourcemap_['latitude'][0]*u.rad,
                         'abundance':sourcemap_['map'][0].transpose(),
                         'coordinate_system':str(sourcemap_['coordinate_system'][0])}
        elif isinstance(sourcemap_, dict):
            sourcemap = sourcemap_
        else:
            raise InputError('MESSENGERuvvs.compute_radiance_from_map: '
                             'Improper sourcemap format')
        
        # Create the LOSResult object
        model_result = LOSResult(self, 'radiance')
    
        # Compute fit to data
        model_result.recompute_radiance(unfitnum, sourcemap)
    
        # Attach the model_result to the data
        modnum = len(self.model_info)
        modkey = f'model{modnum:02d}'
        npackkey = f'npackets{modnum:02d}'
        maskkey = f'mask{modnum:02d}'
    
        self.data[modkey] = model_result.radiance / 1e3  # Convert to kR
        self.data[npackkey] = model_result.npackets
        strength, goodness_of_fit, mask = mathMB.fit_model(self.data.radiance.values,
                                                           self.data[modkey].values,
                                                           self.data.sigma.values,
                                                           fit_method=oldinfo['fit_method'],
                                                           masking=oldinfo['masking'],
                                                           altitude=self.data.alttan)
        strength *= u.def_unit('10**23 atoms/s', 1e23 / u.s)
        self.data[modkey] = self.data[modkey] * strength.value
        self.data[maskkey] = mask
    
        if label is None:
            label = modkey.capitalize()
        else:
            pass
    
        model_info = {'inputs':model_result.inputs,
                      'fit_method':oldinfo['fit_method'],
                      'masking':oldinfo['masking'],
                      'goodness-of-fit':goodness_of_fit,
                      'strength':strength,
                      'label':label,
                      'outputfiles':model_result.outputfiles,
                      'fitted': False,
                      'modelfiles':model_result.modelfiles,
                      'sourcemap':model_result.sourcemap}
    
        self.model_info[modkey] = model_info
        print(f'Model strength for {label} = {strength}')

    def determine_source(self, unfitnum, label=None, weight_method='scaling'):
        oldkey = f'model{unfitnum:02d}'
        oldinfo = self.model_info[oldkey]
        
        # Create the LOSResult object
        model_result = LOSResult(self, 'radiance')
        
        # Compute fit to data
        model_result.determine_source_from_data(unfitnum,
                                                weight_method=weight_method)
        
        # Attach the model_result to the data
        modnum = len(self.model_info)
        modkey = f'model{modnum:02d}'
        npackkey = f'npackets{modnum:02d}'
        maskkey = f'mask{modnum:02d}'

        self.data[modkey] = model_result.radiance / 1e3  # Convert to kR
        self.data[npackkey] = model_result.npackets
        strength, goodness_of_fit, mask = mathMB.fit_model(self.data.radiance.values,
                                                           self.data[modkey].values,
                                                           self.data.sigma.values,
                                                           fit_method=oldinfo['fit_method'],
                                                           masking=oldinfo['masking'],
                                                           altitude=self.data.alttan)
        strength *= u.def_unit('10**23 atoms/s', 1e23 / u.s)
        self.data[modkey] = self.data[modkey] * strength.value
        self.data[maskkey] = mask

        if label is None:
            label = modkey.capitalize()
        else:
            pass

        model_info = {'inputs':model_result.inputs,
                      'fit_method': oldinfo['fit_method'],
                      'masking': oldinfo['masking'],
                      'goodness-of-fit':goodness_of_fit,
                      'strength':strength,
                      'label':label,
                      'outputfiles':model_result.outputfiles,
                      'fitted': True,
                      'modelfiles':model_result.modelfiles,
                      'sourcemap':model_result.sourcemap}

        self.model_info[modkey] = model_info
        print(f'Model strength for {label} = {strength}')

    def model(self,
              inputs,
              npackets,
              quantity='radiance',
              fit_method='chisq',
              dphi=1*u.deg,
              overwrite=False,
              masking=None,
              label=None,
              packs_per_it=None,
              save_packets=False):
        """Run the nexoclom model with specified inputs and fit to the data.
        
        ** Parameters**
        inputs
            A nexoclom Input object or a string with the path to an inputs file.
            
        npackets
            Number of packets to run.
            
        fit_method
            Allows user to specify the quantity to be minimized when fitting
            the model to the data. Default = chisq (chi-squared
            minimization). Other option is 'difference' which minimizes the
            sum of the differences of the data and the model ignoring the
            data uncertainty.
            
        dphi
            Angular size of the view cone. Default = 1 deg. Must be
            given as an astropy quantity.
            
        overwrite
            Set to True to erase any previous model runs with these inputs.
            Default = False
        
        masking
            Allows user to specify which data points are included in the fit
            to the model. Default = None
            * middleX: Use the middle X% of values. For example, middle50
              excludes the faintest and brightest 25% spectra
            
            * minsnrX: specifies the minimum signal-to-noise ratio to use.
              minsnr3 exludes any spectra with SNR < 3
            
            * siglimitX: Excludes any data where the model is not within X σ
              of the data in an initial fit to the data.
        
        start_from_fit
            Uses a previous fit to the data as the starting point for the model
            run. Creates a source map based on the fit. Default = False
            * If start_from_fit = True, the inputs must be the label of the
              model run that was fit to the data. Must have:
              data.model_info[inputs]['fitted'] = True
        
        start_from_fit_options
            A dictionary with the fields:
            * smooth: True if the derived surface map should be smoothed
              before running. Default = False
            * nlonbins: Number of longitude bins in grid. Default = 72
            * nlatbins: Number of latitude bins in grid. Default = 36
            * nvelbins: Number of velocity bins: Default = 100
            * mapfile: File to which the derived source should be saved.
              Should end in '.pkl'. Default='mapfile_temp.pkl'
            * The speed distribution needs to be defined. It can be defined
              using the standard SpeedDistribution fields or with
              type = 'from fit' which approximates the derived speed
              distribution. Default is 'from fit'
        """
        # Put data into model units
        self.set_frame('Model')
        
        # Create the LOSResult object
        model_result = LOSResult(self, quantity, dphi=dphi)
        
        # simulate the data
        inputs.options.fitted = False
        model_result.simulate_data_from_inputs(inputs, npackets, overwrite,
                                               packs_per_it, save_packets=save_packets)
        
        # Attach the model_result to the data
        modnum = len(self.model_info)
        modkey = f'model{modnum:02d}'
        npackkey = f'npackets{modnum:02d}'
        maskkey = f'mask{modnum:02d}'

        self.data[modkey] = model_result.radiance / 1e3 # Convert to kR
        self.data[npackkey] = model_result.npackets
        strength, goodness_of_fit, mask = mathMB.fit_model(self.data.radiance.values,
                                                           self.data[modkey].values,
                                                           self.data.sigma.values,
                                                           fit_method=fit_method,
                                                           masking=masking,
                                                           altitude=self.data.alttan)
        strength *= u.def_unit('10**23 atoms/s', 1e23 / u.s)
        self.data[modkey] = self.data[modkey] * strength.value
        self.data[maskkey] = mask
    
        if label is None:
            label = modkey.capitalize()
        else:
            pass
    
        model_info = {'inputs': model_result.inputs,
                      'fit_method': fit_method,
                      'masking': masking,
                      'goodness-of-fit': goodness_of_fit,
                      'strength': strength,
                      'label': label,
                      'outputfiles': model_result.outputfiles,
                      'fitted': False,
                      'modelfiles': model_result.modelfiles,
                      'sourcemap': model_result.sourcemap}

        self.model_info[modkey] = model_info
        print(f'Model strength for {label} = {strength}')

    def plot(self, filename=None, plot_method='plotly', show=False):
        if plot_method == 'plotly':
            app = plot_plotly(self, filename)
            if show:
                app.run_server(debug=False)
            else:
                pass
            self.app = app
        elif plot_method == 'bokeh':
            plot_bokeh(self, filename, show)
        else:
            print('Not a valid plotting method')
            
    def plot_fitted_model(self, filestart='fitted', show=False,
                          make_frames=False, smooth=False, savepng=False):
        plot_fitted(self, filestart, show, make_frames, smooth=smooth,
                    savepng=savepng)

    def export(self, filename, columns=('utc', 'radiance')):
        """Export data and models to a file.
        **Parameters**
        
        filename
            Filename to export model results to. The file extension determines
            the format. Formats available: csv, pkl, html, tex
            
        columns
            Columns from the data dataframe to export. Available columns can
            be found by calling the `keys()` method on the data object.
            Default = ['utc', 'radiance'] and all model result columns. Note:
            The default columns are always included in the output
            regardless of whether they are specified.
        
        **Returns**
        
        No outputs.
        
        """
        columns_ = list(columns)
        if len(self.model_info) != 0:
            columns_.extend(self.model_info.keys())
        else:
            pass
        
        # Make sure radiance is in there
        if 'radiance' not in columns_:
            columns_.append('radiance')
        else:
            pass

        # Make sure UTC is in there
        if 'utc' not in columns_:
            columns_.append('utc')
        else:
            pass

        if len(columns_) != len(set(columns_)):
            columns_ = list(set(columns_))
        else:
            pass

        for col in columns_:
            if col not in self.data.columns:
                columns_.remove(col)
            else:
                pass

        subset = self.data[columns_]
        if filename.endswith('.csv'):
            subset.to_csv(filename)
        elif filename.endswith('.pkl'):
            subset.to_pickle(filename)
        elif filename.endswith('.html'):
            subset.to_html(filename)
        elif filename.endswith('.tex'):
            subset.to_latex(filename)
        else:
            print('Valid output formats = csv, pkl, html, tex')

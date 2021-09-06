import os
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve

import bokeh.plotting as bkp
from bokeh.models import (HoverTool, Whisker, CDSView, BooleanFilter,
                          ColorBar, LinearColorMapper)
from bokeh.models.tickers import SingleIntervalTicker, DatetimeTicker, FixedTicker
from bokeh.layouts import column, gridplot
from bokeh.palettes import Set1, Turbo256
from bokeh.io import export_png, curdoc
from bokeh.themes import Theme

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash.dependencies as ddep
import plotly.io as pio
import plotly.graph_objects as go


WIDTH, HEIGHT = 1000, 500
FONTSIZE, NUMFONTSIZE = '16pt', '12pt'

BOKEH_THEME_FILE = os.path.join(os.path.dirname(__file__), 'data', 'bokeh.yml')


def color_generator(pallette=Set1[9]):
    i = -1
    while True:
        i = (i + 1) % len(pallette)
        yield pallette[i]


def plot_bokeh(self, filename=None, show=True, savepng=False):
    curdoc().theme = Theme(BOKEH_THEME_FILE)
    
    if filename is not None:
        if not filename.endswith('.html'):
            filename += '.html'
        else:
            pass
        bkp.output_file(filename)
    else:
        pass
    
    # Format the date correction
    self.data['utcstr'] = self.data['utc'].apply(
            lambda x:x.isoformat()[0:19])
    
    # Put the dataframe in a useable form
    self.data['lower'] = self.data.radiance-self.data.sigma
    self.data['upper'] = self.data.radiance+self.data.sigma
    self.data['lattandeg'] = self.data.lattan*180/np.pi
    
    mask = self.data.alttan != self.data.alttan.max()
    if np.any(mask):
        m = self.data[self.data.alttan != self.data.alttan.max()].alttan.max()
    else:
        m = 1e10
    col = np.interp(self.data.alttan, np.linspace(0, m, 256),
                    np.arange(256)).astype(int)
    self.data['color'] = [Turbo256[c] for c in col]
    source = bkp.ColumnDataSource(self.data)
    
    # Tools
    tools = ['pan', 'box_zoom', 'wheel_zoom', 'xbox_select',
             'hover', 'reset', 'save']
    
    # tool tips
    tips = [('index', '$index'),
            ('UTC', '@utcstr'),
            ('Radiance', '@radiance{0.2f} kR'),
            ('LTtan', '@loctimetan{2.1f} hr'),
            ('Lattan', '@lattandeg{3.1f} deg'),
            ('Alttan', '@alttan{0.f} km')]
    
    # Make the figure
    width, height = 1200, 600
    fig0 = bkp.figure(plot_width=width, plot_height=height,
                      x_axis_type='datetime',
                      title=f'{self.species}, {self.query}',
                      x_axis_label='UTC',
                      y_axis_label='Radiance (kR)',
                      y_range=[0, self.data.radiance.max()*1.5],
                      tools=tools, active_drag="xbox_select")
    
    # plot the data
    dplot = fig0.circle(x='utc', y='radiance', size=7, color='black',
                        legend_label='Data', hover_color='yellow',
                        source=source, selection_color='orange')
    fig0.line(x='utc', y='radiance', color='black', legend_label='Data',
              source=source)
    fig0.xaxis.ticker = DatetimeTicker(num_minor_ticks=5)
    
    # Add error bars
    fig0.add_layout(Whisker(source=source, base='utc', upper='upper',
                            lower='lower'))
    renderers = [dplot]
    
    # Plot the model
    col = color_generator()
    modplots, maskedplots = [], []
    for modkey, result in self.model_result.items():
        c = next(col)
        # fig0.line(x='utc', y=modkey, source=source,
        #           legend_label=result.label, color=c)
        
        maskkey = modkey.replace('model', 'mask')
        mask = (self.data[maskkey]).to_list()
        view = CDSView(source=source, filters=[BooleanFilter(mask)])
        modplots.append(fig0.circle(x='utc', y=modkey, size=9, color=c,
                                    source=source, legend_label=result.label,
                                    view=view))
        maskkey = modkey.replace('model', 'mask')
        mask = np.logical_not(self.data[maskkey]).to_list()
        view = CDSView(source=source, filters=[BooleanFilter(mask)])
        maskedplots.append(fig0.circle(x='utc', y=modkey, size=9,
                                       source=source, line_color=c,
                                       fill_color='yellow', view=view,
                                       legend_label=result.label +
                                                    ' (Data Point Not Used)'))
        renderers.extend(modplots)
        renderers.extend(maskedplots)
    
    datahover = HoverTool(tooltips=tips, renderers=renderers)
    fig0.add_tools(datahover)
    
    ##############
    # Plot tangent point
    color_mapper = LinearColorMapper(palette="Turbo256", low=0, high=m)
    
    width, height = 1200, 600
    tools = ['pan', 'box_zoom', 'wheel_zoom', 'box_select',
             'hover', 'reset', 'save']
    fig1 = bkp.figure(plot_width=width, plot_height=height,
                      title=f'Tangent Point Location',
                      x_axis_label='Local Time (hr)',
                      y_axis_label='Latitude (deg)',
                      x_range=[0, 24],
                      y_range=[-90, 90], tools=tools,
                      active_drag="box_select")
    tanplot = fig1.circle(x='loctimetan', y='lattandeg', size=5,
                          selection_color='orange', hover_color='purple',
                          source=source, color='color')
    fig1.xaxis.ticker = SingleIntervalTicker(interval=6,
                                             num_minor_ticks=6)
    fig1.yaxis.ticker = SingleIntervalTicker(interval=45,
                                             num_minor_ticks=3)
    color_bar = ColorBar(color_mapper=color_mapper, title='Altitude (km)',
                         label_standoff=12, border_line_color=None,
                         location=(0, 0))
    fig1.add_layout(color_bar, 'right')
    datahover = HoverTool(tooltips=tips, renderers=[tanplot])
    fig1.add_tools(datahover)
    
    grid = column(fig0, fig1)
    
    if filename is not None:
        bkp.output_file(filename)
        if savepng:
            export_png(grid, filename=filename.replace('.html', '.png'))
        else:
            pass
        bkp.save(grid)
    else:
        pass
    
    if show:
        bkp.show(grid)
    
    return fig0, fig1


def plot_plotly(self, filename=None):
    pio.templates.default = 'plotly_white'
    
    # Format the date correction
    self.data['utcstr'] = self.data['utc'].apply(
            lambda x:x.isoformat()[0:19])
    
    # Convert latitude to degrees
    lattandeg = self.data.lattan*180/np.pi
    
    ######################
    # Make the data figure
    data_figure = go.Figure()
    data_figure.update_layout(title=f'{self.species}, {self.query}',
                              xaxis_title='UTC',
                              yaxis_title='Radiance (kR)',
                              font={'family':'Arial', 'size':12},
                              hovermode='x',
                              clickmode='event+select',
                              )
    data_figure.update_xaxes(rangeslider_visible=True)

    datafig = go.Scatter(name='Data',
                         x=self.data.utc,
                         y=self.data.radiance,
                         error_y={'type':'data',
                                  'array':self.data.sigma.values},
                         mode='lines+markers',
                         marker={'color':'black', 'size':3},
                         line={'color':'black'},
                         )
    data_figure.add_trace(datafig)
    dmax = self.data.radiance.max()
    
    # Add in the model lines
    colors = {'model0':'red',
              'model1':'green',
              'model2':'blue',
              'model3':'cyan',
              'model4':'purple'}
    for modkey, result in self.model_result.items():
        maskkey = modkey.replace('model', 'mask')
        marker_line_color = self.data[maskkey].apply(
                lambda x:colors[modkey] if x else 'yellow')
        # marker = self.data[maskkey].apply(lambda x: 0 if x else 100)
        modfig = go.Scatter(name=result.label,
                            x=self.data.utc,
                            y=self.data[modkey],
                            mode='markers+lines',
                            marker={'size':3, 'color':colors[modkey],
                                    'line':{'color':marker_line_color,
                                            'width':1}},
                            line={'color':colors[modkey]},
                            )
        data_figure.add_trace(modfig)
        dmax = np.max([dmax, self.data[modkey].max()])
    
    dataline = go.Scatter(name='Current Spectrum',
                          x=[self.data.iloc[0].utc, self.data.iloc[0].utc],
                          y=[0, 1e5],
                          mode='lines',
                          line={'color':'purple'})
    data_figure.add_trace(dataline)
    data_figure.update_yaxes(range=[0, dmax*1.2])
    if filename is not None:
        data_figure.write_html(filename)

    ##########################
    # Make the tangent point plot
    amax = self.data[self.data.alttan != self.data.alttan.max()].alttan.max()
    tan_figure = go.Figure()
    tan_figure.update_layout(title=f'Tangent Point',
                             xaxis_title='Local Time (hr)',
                             yaxis_title='Latitude (deg)',
                             font={'family':'Arial', 'size':12},
                             clickmode='event+select',
                             hovermode='closest',
                             )
    tan_figure.update_xaxes(range=[0, 24])
    tan_figure.update_yaxes(range=[-90, 90])
    tanfig = go.Scatter(name='TanPoint',
                        x=self.data.loctimetan,
                        y=lattandeg,
                        mode='markers',
                        marker={'size':3,
                                'color':self.data.alttan,
                                'cmin':0, 'cmax':amax,
                                'showscale':True,
                                'colorbar':{'title':'Altitude (km)',
                                            'ticks':'inside'}},
                        )
    tan_figure.add_trace(tanfig)
    
    ##########################
    # Make the Spectrum figure
    if 'wavelength' in self.data:
        # Make the spectrum plot
        spectrum_figure = go.Figure()
        spectrum_figure.update_layout(
                title=f'Spectrum: {self.data.iloc[0].utcstr}',
                xaxis_title='Wavelength (nm)',
                yaxis_title='Counts',
                font={'family':'Arial', 'size':12})
        
        rawfig = go.Scatter(name='Raw',
                            x=self.data.iloc[0].wavelength,
                            y=self.data.iloc[0].raw,
                            mode='markers+lines',
                            marker={'color':'blue'},
                            line={'color':'blue'})
        spectrum_figure.add_trace(rawfig)
        
        solarfig = go.Scatter(name='Solar+Dark',
                              x=self.data.iloc[0].wavelength,
                              y=(self.data.iloc[0].solarfit +
                                 self.data.iloc[0].dark),
                              mode='markers+lines',
                              marker={'color':'red'},
                              line={'color':'red'})
        spectrum_figure.add_trace(solarfig)
        
        correctfig = go.Scatter(name='Corrected',
                                x=self.data.iloc[0].wavelength,
                                y=(self.data.iloc[0].raw -
                                   self.data.iloc[0].solarfit -
                                   self.data.iloc[0].dark),
                                mode='markers+lines',
                                marker={'color':'black'},
                                line={'color':'black'})
        spectrum_figure.add_trace(correctfig)
    else:
        spectrum_figure = None
    
    ##################
    # Make the dash app
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    app = dash.Dash(f'MESSENGER UVVS {self.species} {self.query}',
                    external_stylesheets=external_stylesheets)
    
    app.layout = html.Div([
        html.H1(f'MESSENGER UVVS {self.query}',
                style={'textAlign':'center'}),
        
        html.Div([
            dcc.Slider(id='UTCSlider',
                       min=0, max=len(self),
                       step=1,
                       value=0,
                       updatemode='drag'),
            html.Div(id='UTCValue'),
            ]),
        
        html.Table([
            html.Tr([
                html.Td(dcc.Graph(id='DataGraph',
                                  figure=data_figure),
                        colSpan=2),
                ]),
            html.Tr([
                html.Td(dcc.Graph(id='TangentGraph')),
                html.Td(dcc.Graph(id='SpectrumGraph',
                                  figure=spectrum_figure))
                ]),
            ]),
        
        html.Div(id='SelectionData', hidden=True)
        ])
    
    @app.callback(
            [ddep.Output('UTCValue', 'children'),
             ddep.Output('DataGraph', 'figure'),
             ddep.Output('TangentGraph', 'figure'),
             ddep.Output('SpectrumGraph', 'figure')],
            [ddep.Input('UTCSlider', 'value'),
             ddep.Input('SelectionData', 'children')])
    def update_figures_from_slider(slider_value, selection_data):
        ctx = dash.callback_context
        triggered = ctx.triggered[0]['prop_id']
        
        if triggered == 'UTCSlider.value':
            # Update UTCValue
            utcvalue = f'UTC = {self.data.iloc[slider_value].utcstr}'
            
            # Update the DataGraph
            data_figure.update_traces(x=[self.data.iloc[slider_value].utc,
                                         self.data.iloc[slider_value].utc],
                                      selector={'name':'Current Spectrum'})
            data_figure.update_yaxes(range=[0, dmax*1.2])
            
            # Update the TangentPlot
            marker_size = np.zeros(len(self))+3
            marker_size[slider_value] = 10
            marker_symbol = ['circle' for _ in range(len(self))]
            marker_symbol[slider_value] = 'x'
            tan_figure.update_traces(marker={'size':marker_size,
                                             'symbol':marker_symbol},
                                     selector={'name':'TanPoint'})
            
            # Update spectrum figure
            if 'wavelength' in self.data:
                spectrum_figure.update_traces(
                        x=self.data.iloc[slider_value].wavelength,
                        y=self.data.iloc[slider_value].raw,
                        selector={'name':'Raw'})
                spectrum_figure.update_traces(
                        x=self.data.iloc[slider_value].wavelength,
                        y=(self.data.iloc[slider_value].solarfit +
                           self.data.iloc[slider_value].dark),
                        selector={'name':'Solar+Dark'})
                spectrum_figure.update_traces(
                        x=self.data.iloc[slider_value].wavelength,
                        y=(self.data.iloc[slider_value].raw -
                           self.data.iloc[slider_value].solarfit -
                           self.data.iloc[slider_value].dark),
                        selector={'name':'Corrected'})
            else:
                pass
        elif triggered == 'SelectionData.children':
            if selection_data is None:
                selection_data = [0]
            else:
                pass
            
            # Update UTCValue
            utcvalue = f'UTC = {self.data.iloc[min(selection_data)].utcstr}'
            
            # Update the DataGraph
            data_figure.update_traces(x=[self.data.iloc[min(selection_data)].utc,
                                         self.data.iloc[min(selection_data)].utc],
                                      selector={'name':'Current Spectrum'})
            data_figure.update_yaxes(range=[0, dmax*1.2])
            
            marker_size = np.zeros(len(self))+3
            marker_size[selection_data] = 10
            marker_symbol = np.array(['circle' for _ in range(len(self))])
            marker_symbol[selection_data] = 'x'
            
            data_figure.update_traces(marker={'size':marker_size,
                                              'symbol':marker_symbol},
                                      selector={'name':'Data'})
            tan_figure.update_traces(marker={'size':marker_size,
                                             'symbol':marker_symbol},
                                     selector={'name':'TanPoint'})
        
        else:
            utcvalue = 0
        
        return utcvalue, data_figure, tan_figure, spectrum_figure
    
    @app.callback(
            ddep.Output('UTCSlider', 'value'),
            [ddep.Input('DataGraph', 'clickData'),
             ddep.Input('TangentGraph', 'clickData')])
    def update_slider_value(dataclick, tanclick):
        ctx = dash.callback_context
        triggered = ctx.triggered[0]['prop_id']
        
        if triggered == 'DataGraph.clickData':
            return dataclick['points'][0]['pointIndex']
        elif triggered == 'TangentGraph.clickData':
            return tanclick['points'][0]['pointIndex']
        else:
            return 0
    
    @app.callback(
            ddep.Output('SelectionData', 'children'),
            [ddep.Input('DataGraph', 'selectedData'),
             ddep.Input('TangentGraph', 'selectedData')])
    def update_selection_data(dataselect, tanselect):
        ctx = dash.callback_context
        triggered = ctx.triggered[0]['prop_id']
        
        if triggered == 'DataGraph.selectedData':
            selected_points = [point['pointIndex']
                               for point in dataselect['points']]
            return selected_points
        elif triggered == 'TangentGraph.selectedData':
            selected_points = [point['pointIndex']
                               for point in tanselect['points']]
            return selected_points
        else:
            return None
    
    return app


def make_fitted_plot(self, result, filestart='fitted', show=True, ut=None,
                     smooth=False, savepng=False):
    curdoc().theme = Theme(BOKEH_THEME_FILE)
    
    if smooth:
        kernel = Gaussian2DKernel(x_stddev=1)
        source = convolve(result['abundance'], kernel, boundary='wrap')
        packets = convolve(result['p_available'], kernel, boundary='wrap')
    else:
        source = result['abundance']
        packets = result['p_available']
    
    # Tools
    tools = ['save']

    local_time = (result['longitude'].value * 12/np.pi + 12) % 24
    arg = np.argsort(local_time[:-1])
    source, packets = source[arg,:], packets[arg,:]

    # Distribution of available packets
    fig0 = bkp.figure(plot_width=WIDTH, plot_height=HEIGHT,
                      title=f'{self.species}, {self.query}, Available Packets',
                      x_axis_label='Local Time (hr)',
                      y_axis_label='Latitude (deg)',
                      x_range=[0, 24], y_range=[-90, 90],
                      tools=tools)
    fig0.title.text_font_size = FONTSIZE
    fig0.xaxis.axis_label_text_font_size = FONTSIZE
    fig0.yaxis.axis_label_text_font_size = FONTSIZE
    fig0.xaxis.major_label_text_font_size = NUMFONTSIZE
    fig0.yaxis.major_label_text_font_size = NUMFONTSIZE
    fig0.xaxis.ticker = FixedTicker(ticks=[0, 6, 12, 18, 24])
    fig0.yaxis.ticker = FixedTicker(ticks=[-90, 45, 0, 45, 90])
    
    fig0.image(image=[packets.transpose()],
               x=0, y=-90, dw=24, dh=180, palette='Spectral11')
    
    # Distribution of packets used in the final model
    fig1 = bkp.figure(plot_width=WIDTH, plot_height=HEIGHT,
                      title=f'{self.species}, {self.query}, Packets Used',
                      x_axis_label='Local Time (hr)',
                      y_axis_label='Latitude (deg)',
                      x_range=[0, 24], y_range=[-90, 90],
                      tools=tools)
    fig1.title.text_font_size = FONTSIZE
    fig1.xaxis.axis_label_text_font_size = FONTSIZE
    fig1.yaxis.axis_label_text_font_size = FONTSIZE
    fig1.xaxis.major_label_text_font_size = NUMFONTSIZE
    fig1.yaxis.major_label_text_font_size = NUMFONTSIZE
    fig1.xaxis.ticker = FixedTicker(ticks=[0, 6, 12, 18, 24])
    fig1.yaxis.ticker = FixedTicker(ticks=[-90, 45, 0, 45, 90])
    
    fig1.image(image=[source.transpose()],
               x=0, y=-90, dw=24, dh=180, palette='Spectral11')
    
    fig2 = bkp.figure(plot_width=WIDTH, plot_height=HEIGHT,
                      title=f'{self.species}, {self.query}, Speed Distribution',
                      x_axis_label='Speed (km/s)',
                      y_axis_label='Relative Number',
                      y_range=[0, 1.2],
                      tools=tools)
    fig2.title.text_font_size = FONTSIZE
    fig2.xaxis.axis_label_text_font_size = FONTSIZE
    fig2.yaxis.axis_label_text_font_size = FONTSIZE
    fig2.xaxis.major_label_text_font_size = NUMFONTSIZE
    fig2.yaxis.major_label_text_font_size = NUMFONTSIZE
    
    fig2.line(x=result['velocity'][:-1], y=result['v_available'],
              legend_label='Packets Available', color='red')
    fig2.line(x=result['velocity'][:-1], y=result['vdist'],
              legend_label='Packets Used', color='blue')
    
    # Full orbit time series
    # Format the date correction
    self.data['utcstr'] = self.data['utc'].apply(
            lambda x:x.isoformat()[0:19])
    
    # Put the dataframe in a useable form
    self.data['lower'] = self.data.radiance-self.data.sigma
    self.data['upper'] = self.data.radiance+self.data.sigma
    self.data['lattandeg'] = self.data.lattan*180/np.pi
    
    m = self.data[self.data.alttan != self.data.alttan.max()].alttan.max()
    col = np.interp(self.data.alttan, np.linspace(0, m, 256),
                    np.arange(256)).astype(int)
    self.data['color'] = [Turbo256[c] for c in col]
    source = bkp.ColumnDataSource(self.data)
    
    # Tools
    tools = ['pan', 'box_zoom', 'wheel_zoom', 'xbox_select',
             'hover', 'reset', 'save']
    
    # tool tips
    tips = [('index', '$index'),
            ('UTC', '@utcstr'),
            ('Radiance', '@radiance{0.2f} kR'),
            ('LTtan', '@loctimetan{2.1f} hr'),
            ('Lattan', '@lattandeg{3.1f} deg'),
            ('Alttan', '@alttan{0.f} km')]
    
    # Make the radiance figure
    title_ = f'{self.species}, {self.query}'
    if ut is not None:
        title_ += f', UTC = {ut.isoformat()}'
    else:
        pass
    
    fig3 = bkp.figure(plot_width=WIDTH, plot_height=HEIGHT,
                      x_axis_type='datetime',
                      title=title_,
                      x_axis_label='UTC',
                      y_axis_label='Radiance (kR)',
                      y_range=[0, self.data.radiance.max()*1.5],
                      tools=tools, active_drag="xbox_select")
    fig3.title.text_font_size = FONTSIZE
    fig3.xaxis.axis_label_text_font_size = FONTSIZE
    fig3.yaxis.axis_label_text_font_size = FONTSIZE
    fig3.xaxis.major_label_text_font_size = NUMFONTSIZE
    fig3.yaxis.major_label_text_font_size = NUMFONTSIZE
    
    # plot the data
    dplot = fig3.circle(x='utc', y='radiance', size=7, color='black',
                        legend_label='Data', hover_color='yellow',
                        source=source, selection_color='orange')
    fig3.line(x='utc', y='radiance', color='black', legend_label='Data',
              source=source)
    fig3.xaxis.ticker = DatetimeTicker(num_minor_ticks=5)
    
    # Add error bars
    fig3.add_layout(Whisker(source=source, base='utc', upper='upper',
                            lower='lower'))
    renderers = [dplot]
    
    # Plot the model
    col = color_generator()
    modplots, maskedplots = [], []
    for modkey, result in self.model_result.items():
        if result.fitted:
            c = next(col)
            fig3.line(x='utc', y=modkey, source=source,
                      legend_label=result.label, color=c)

            maskkey = modkey.replace('model', 'mask')
            mask = self.data[maskkey].to_list()
            view = CDSView(source=source, filters=[BooleanFilter(mask)])
            modplots.append(fig3.circle(x='utc', y=modkey, size=7, color=c,
                                        source=source, legend_label=result.label,
                                        view=view))
            
            maskkey = modkey.replace('model', 'mask')
            mask = np.logical_not(self.data[maskkey]).to_list()
            view = CDSView(source=source, filters=[BooleanFilter(mask)])
            maskedplots.append(fig3.circle(x='utc', y=modkey, size=7,
                                           source=source, line_color=c,
                                           fill_color='yellow', view=view,
                           legend_label=result.label + '(Data Point Not Used)'))
            renderers.extend(modplots)
            renderers.extend(maskedplots)
        
        if ut is not None:
            yr = fig3.y_range
            fig3.line(x=[ut, ut], y=[0, 1e5], color='purple')
            fig3.y_range = yr
        else:
            pass
    
    datahover = HoverTool(tooltips=tips, renderers=renderers)
    fig3.add_tools(datahover)
    
    grid = gridplot([[fig3, fig2], [fig0, fig1]])
    
    # Save png version
    if savepng:
        export_png(grid, filename=filestart+'.png')
    else:
        pass
    
    bkp.output_file(filestart+'.html')
    bkp.save(grid)  # html files not needed
    
    if show:
        bkp.show(grid)
    else:
        pass
    
    return grid


def plot_fitted(self, filestart=None, show=True, make_frames=False,
                smooth=False, savepng=False):
    # Check that the path exists
    if ((os.path.dirname(filestart) != '') and
        (not os.path.exists(os.path.dirname(filestart)))):
        os.makedirs(os.path.dirname(filestart))
    
    # Compute and plot final source results
    nfitted = sum(self.model_result[result].fitted for result in self.model_result)
    for key, result in self.model_result.items():
        if result.fitted:
            final_result = self.model_info[key]['sourcemap']
            if nfitted > 1:
                filestart_ = f'{filestart}_{key}'
            else:
                filestart_ = filestart
            make_fitted_plot(self, final_result, filestart=filestart_, show=show,
                             smooth=smooth, savepng=savepng)
    
            if make_frames:
                assert 0
            else:
                pass
        else:
            pass

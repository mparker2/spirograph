"""
A script for drawing nice spirals
"""

from itertools import chain, repeat
from operator import attrgetter
import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection

# initial parameter values:
tmax0 = 12 * np.pi
tstep0 = 0.1
R0 = 20
r0 = 10
p0 = 10
col0 = 1
a0 = 1
b0 = 0.5
c0 = 0


class PiString(str):
    '''
    object which can be string formatted (internally by matplotlib)
    with number and displays the number in terms of pi.
    '''
    def __mod__(self, number):
        return u'{:.1f}\u03C0'.format(number / np.pi)


def spiro_linefunc(t, R, r, p):
    '''spiralograph x transformation'''
    x = (R - r) * np.cos(t) + p * np.cos(((R - r) / r) * t)
    y = (R - r) * np.sin(t) - p * np.sin(((R - r) / r) * t)
    return x, y


def spiro_linewidths(t, a, b, c):
    '''calc variable linewidth'''
    return np.sin(t * a) * b + c


def segments(x, y):
    '''
    convert x and y points into segments to plot as individual lines
    using matplotlib LineCollection
    '''
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


class SpiroGraph(object):

    '''
    Spirograph drawer with matplotlib slider widgets to change parameters.
    Parameters of line are:

        R: The radius of the big circle
        r: The radius of the small circle which rolls along the inside of the
           bigger circle
        p: distance from centre of smaller circle to point in the circle where
           the pen hole is.
        tmax: the angle through which the smaller circle is rotated to draw the
              spirograph
        tstep: how often matplotlib plots a point
        a, b, c: parameters of the linewidth equation.
    '''

    # kwargs for each of the matplotlib sliders
    slider_kwargs = (
        {'label': 't_max', 'valmin': np.pi, 'valmax': 200 * np.pi,
         'valinit': tmax0, 'valfmt': PiString()},
        {'label': 't_step', 'valmin': 0.01,
         'valmax': 10, 'valinit': tstep0},
        {'label': 'R', 'valmin': 1, 'valmax': 200, 'valinit': R0},
        {'label': 'r', 'valmin': 1, 'valmax': 200, 'valinit': r0},
        {'label': 'p', 'valmin': 1, 'valmax': 200, 'valinit': p0},
        {'label': 'colour', 'valmin': 0, 'valmax': 1, 'valinit': 1},
        {'label': 'width_a', 'valmin': 0.5, 'valmax': 10, 'valinit': 1},
        {'label': 'width_b', 'valmin': 0, 'valmax': 10, 'valinit': 0},
        {'label': 'width_c', 'valmin': 0, 'valmax': 10, 'valinit': 0.5})

    rbutton_kwargs = (
        {'labels': ('black', 'white'), 'activecolor': 'white', 'active': 0},
        {'labels': ('solid', 'variable'), 'activecolor': 'white', 'active': 0})

    def __init__(self, colormap, figsize=(7, 10)):
        self.colormap_name = colormap
        self.variable_color = False
        # Use ScalarMappable to map full colormap to range 0 - 1
        self.colormap = ScalarMappable(cmap=colormap)
        self.colormap.set_clim(0, 1)

        # set up main axis onto which to draw spirograph
        self.figsize = figsize
        plt.rcParams['figure.figsize'] = figsize
        self.fig, self.mainax = plt.subplots()
        plt.subplots_adjust(bottom=0.3)
        title = self.mainax.set_title('Spirograph Drawer!',
                                      size=20,
                                      color='white')
        self.text = [title, ]
        # set up slider axes
        self.slider_axes = [plt.axes([0.25, x, 0.65, 0.015])
                            for x in np.arange(0.05, 0.275, 0.025)]
        # same again for radio buttons
        self.rbutton_axes = [plt.axes([0.025, x, 0.1, 0.15])
                             for x in np.arange(0.02, 0.302, 0.15)]
        # use log scale for tstep slider
        self.slider_axes[1].set_xscale('log')
        # turn off frame, ticks and tick labels for all axes
        for ax in chain(self.slider_axes, self.rbutton_axes, [self.mainax, ]):
            ax.axis('off')
        # use axes and kwargs to create list of sliders/rbuttons
        self.sliders = [Slider(ax, **kwargs)
                        for ax, kwargs in zip(self.slider_axes,
                                              self.slider_kwargs)]
        self.rbuttons = [RadioButtons(ax, **kwargs)
                         for ax, kwargs in zip(self.rbutton_axes,
                                               self.rbutton_kwargs)]
        self.update_figcolors()

        # set up initial line
        self.t = np.arange(0, tmax0, tstep0)
        x, y = spiro_linefunc(self.t, R0, r0, p0)
        self.linecollection = LineCollection(
            segments(x, y),
            linewidths=spiro_linewidths(self.t, a0, b0, c0),
            color=self.colormap.to_rgba(col0))
        self.mainax.add_collection(self.linecollection)

        # creates the plot and connects sliders to various update functions
        self.run()

    def update_figcolors(self, bgcolor='black'):
        '''
        function run by background color radiobutton. Sets all labels, text,
        and sliders to foreground color, all axes to background color
        '''
        fgcolor = 'white' if bgcolor == 'black' else 'black'
        self.fig.set_facecolor(bgcolor)
        self.mainax.set_axis_bgcolor(bgcolor)
        for ax in chain(self.slider_axes, self.rbutton_axes):
            ax.set_axis_bgcolor(bgcolor)

        # set fgcolor elements to black or white, mostly elements of sliders
        for item in chain(map(attrgetter('label'), self.sliders),
                          map(attrgetter('valtext'), self.sliders),
                          map(attrgetter('poly'), self.sliders),
                          self.text,
                          *map(attrgetter('labels'), self.rbuttons)):
            item.set_color(fgcolor)

        self.update_radiobutton_colors()
        plt.draw()

    def update_linewidths(self, *args):
        '''
        function run by a, b and c parameter sliders. Sets width of each line
        in linecollection according to sine function
        '''
        a, b, c = (s.val for s in self.sliders[6:])
        self.linecollection.set_linewidths(spiro_linewidths(self.t, a, b, c))
        plt.draw()

    def update_linecolors(self, *args):
        '''
        function run by color slider and indirectly by variable/solid color
        radiobutton. Updates colors of each line in linecollection using the
        set colormap.
        '''
        # get current color value (a value between 1 and 0)
        col_val = self.sliders[5].val
        if not self.variable_color:
            # if solid color, convert color value to rgb and set the color
            self.linecollection.set_color(self.colormap.to_rgba(col_val))
        else:
            # create values between 0 and 1 for each line segment
            colors = (self.t / max(self.t)) + col_val
            # use color value to roll colors
            colors[colors > 1] -= 1
            self.linecollection.set_color(
                [self.colormap.to_rgba(i) for i in colors])
        plt.draw()

    def update_lineverts(self, *args):
        '''
        function run by R, r, p, tmax and tstep sliders to update line vertices
        '''
        tmax, tstep, R, r, p = (s.val for s in self.sliders[:5])
        self.t = np.arange(0, tmax, tstep)
        x, y = spiro_linefunc(self.t, R, r, p)
        self.linecollection.set_verts(segments(x, y))
        # change axis limits to pad new line nicely
        self.mainax.set(xlim=(min(x) - 5, max(x) + 5),
                        ylim=(min(y) - 5, max(y) + 5))
        plt.draw()

    def update_linecolor_setting(self, val):
        '''
        function run by solid/variable colour slider, alters variable_color
        attribute then calls update_linecolors
        '''
        if val == 'variable':
            self.variable_color = True
        elif val == 'solid':
            self.variable_color = False
        # need to update radiobutton colors here.
        self.update_radiobutton_colors()
        self.update_linecolors()

    def update_radiobutton_colors(self):
        '''
        makes radiobutton colors correct even on a changing axis background
        '''
        bgcolor = self.rbuttons[0].value_selected
        fgcolor = 'white' if bgcolor == 'black' else 'black'
        for i, b in enumerate(self.rbuttons):
            # find out index of the active button
            active_idx = self.rbutton_kwargs[i]['labels'].index(
                b.value_selected)
            # set button colors accordingly
            b.circles[not active_idx].set_color(bgcolor)
            b.circles[active_idx].set_color(fgcolor)

    def run(self):
        '''
        set up slider functions
        '''
        verts_func = self.update_lineverts
        colors_func = self.update_linecolors
        widths_func = self.update_linewidths
        # create iterable of slider funcs to zip with sliders
        slider_update_funcs = chain(repeat(verts_func, 5),
                                    [colors_func, ],
                                    repeat(widths_func, 3))
        # set slider on_changed functions
        for s, f in zip(self.sliders, slider_update_funcs):
            s.on_changed(f)

        self.rbuttons[0].on_clicked(self.update_figcolors)
        self.rbuttons[1].on_clicked(self.update_linecolor_setting)
        plt.show()

if __name__ == '__main__':
    # get name of the colormap
    a = argparse.ArgumentParser()
    a.add_argument('-c', '--colormap',
                   required=False,
                   default='gist_ncar',
                   help='matplotlib colormap to use')
    args = a.parse_args()
    SpiroGraph(args.colormap)

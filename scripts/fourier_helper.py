# This helper was co-created by Olivier Mesnard and Natalia Clementi

import numpy
from IPython.display import HTML
import ipywidgets
from matplotlib import animation, pyplot


class Circle:
    """Class to create a circle object"""

    def __init__(self, R=0.5, xc=0.0, yc=0.0):
        """Define a circle

        Arguments:
            R {float}  -- Radius of a circle (default: {0.5})
            xc {float} -- x-coordinate of the circle's center (default: {0.0})
            yc {float} -- y-coordinate of the circle's center (default: {0.0})
        """

        self.R = R
        self.xc, self.yc = xc, yc

    def coordinates(self, theta=None, num=50):
        """Get the coordinates on the circle

        Arguments:
        theta {array, float} -- Array/float of angles in radians
                                (default: {None})
        num {int}            -- Number of points to discretize theta if theta
                                is None (default: {50})
        Returns:
        x {numpy.ndarray of floats} -- x-coordinates of the points
                                       on the circle (1D array).
        y {numpy.ndarray of floats} -- y-coordinates of the points on
                                       the circle (1D array).
        """
        if theta is None:
            theta = numpy.linspace(0.0, 2 * numpy.pi, num=num)
        x = self.xc + self.R * numpy.cos(theta)
        y = self.yc + self.R * numpy.sin(theta)
        return x, y


class Element:
    """An Element contains a circle and function."""

    def __init__(self, f, shift_x=0.0, shift_y=0.0):
        """Initialize the frequency and create its associated circle.

        Arguments
        ----------
        f : float
            The frequency.
        shift_x : float (optional)
            Horizontal shift of the circle (x-coordinate of the center);
            default: 0.0.
        shift_y : float (optional)
            Vertical shift of the circle (y-coordinate of the center);
            default: 0.0.

        """
        self.f = f
        self.circle = Circle(R=1 / (self.f), xc=shift_x, yc=shift_y) 
        self.label = r'$1\sin({0}\theta)/{0}$'.format(self.f)
    
    def signal(self, theta):
        """Compute the signal.

        Arguments
        ----------
        theta : float or numpy.ndarray of floats
                Angle(s) in radians.

        Returns
        -------
        s : float or numpy.ndarray of floats
            The signal.

        """
        s = self.circle.R * numpy.sin(self.f * theta)
        return s
    
    def hline(self, t, xmax=None):
        """Get the position of the horizontal line.

        Parameters
        ----------
        t    : float
               Angle in radians.
        xmax : float (optional)
               Right-hand side x-coordinate of the horizontal line;
               default: None

        Returns
        -------
        pos : tuple
              Position of the horizontal line.
        """

        xmin, y = self.circle.coordinates(theta=self.f * t)
        if xmax is None:
            xmax = t
        return ((xmin, xmax), (y, y))

    def rline(self, t):
        """Get the position of the radial line from the center.
        
        Parameters
        ----------
        t : float
            Angle in radians.

        Returns
        -------
        pos : tuple
            Position of the horizontal line.
        """
        x, y = self.circle.coordinates(theta=self.f * t)
        return ((self.circle.xc, x), (self.circle.yc, y))


def create_elements_center(freq, x_shift, y_shift):
    """ It creates center elements to produce animation
    for different waves and respective circles based on
    their frequency.

    Arguments
    ---------
    freq: float or list. Frequency of desired waves.
    x_shift: float, x-shift of the circle (x-coord of center)
    y_shift: float, y-shift of the circle (y-coord of center)

    Returns
    -------
    elems: list with the elements all center on same coordinates.
    """

    elems = [Element(f, shift_x=x_shift, shift_y=y_shift) for f in freq]
    
    return elems


def define_angles(periods, steps_per_period):
    """Creates array of angles that goes from 0 to 2*pi*periods, and contains
    periods*steps_per_period number of steps.

    Arguments
    ---------
    periods: float, number of periods.
    steps_per_period: int, number of steps per period.

    Returns
    -------
    theta_arr: numpy array, array of angles.
    """

    steps = periods * steps_per_period
    theta_arr = numpy.linspace(0.0, 2*numpy.pi*periods, num=steps)

    return theta_arr


def create_init_fig_center(elements, theta_arr):
    """ Creates initial figure needed for animation. In this case the 
    figure is the one that has all the circles centered. 

    Arguments
    ---------
    elements: list of objects, list of objects created with the function 
    create_elements_center. 
    theta: array, array of angles created with define_angles.

    Returns
    -------
    Plot with initial figure needed to plot animation. 
    """
    pyplot.ioff() #to avoid displaying figure.

    fig, ax = pyplot.subplots(figsize=(14.0, 6.0))
    ax.grid()
    ax.set_xlabel(r'$\theta$')
    ax.axvline(0.0, color='black')
    signals, hlines, rlines = [], [], []
    for i, elem in enumerate(elements):
        color = 'C' + str(i)
        ax.plot(*elem.circle.coordinates(), color=color, linewidth=2.0)
        hlines.append(ax.plot(*elem.hline(theta_arr[0]), color=color, linestyle='--', marker='o')[0])
        rlines.append(ax.plot(*elem.rline(theta_arr[0]), color=color)[0])
        signals.append(ax.plot(theta_arr[:1], elem.signal(theta_arr[:1]),
                            label=elem.label, color=color, linewidth=2.0)[0])
    ax.legend(ncol=len(elements), loc='lower center', prop={'size': 12})
    ax.axis('scaled', adjustable='box')
    ax.set_xlim(-4.0, 14.0)
    ax.set_ylim(-2.0, 2.0)

    return fig, ax, hlines, rlines, signals


def update_figure_center(n, anim_dict, theta, display_fig=False):
    """Update the figure at a given step.

    Parameters
    ----------
    n : integer
        The step index.
    anim_dict : dict, contains elements needed to create animation.
    theta : numpy.ndarray of floats
        All the angles.
    display_fig : boolean (optional)
        Set True to display the fig with ipywidgets.interactive;
        default: False.

    """
    for i, elem in enumerate(anim_dict['elems']):
        anim_dict['hlines'][i].set_data(elem.hline(theta[n]))
        anim_dict['rlines'][i].set_data(elem.rline(theta[n]))
        anim_dict['signals'][i].set_data(theta[:n + 1], elem.signal(theta[:n + 1]))
    if display_fig:
        display(anim_dict['fig'])
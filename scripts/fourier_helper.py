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
        theta {array} -- Array of angles in radians (default: {None})
        num {int}     -- Number of points to discretize theta if theta is None
                         (default: {50})
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
        shift : float (optional)
            Horizontal shift of the circle (x-coordinate of the center);
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

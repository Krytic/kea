"""
Histogram classes to contain event rate data and allow for easy plotting

Author: Max Briel
"""
import numpy as np
import matplotlib.pyplot as plt
import kea.helpers

class histogram:
    """A histogram which can contains data and can be manipulated.
    Either **xlow**, **xup**, and **nr_bins** is given or **edges**

    As per any histogram, the upper edges are non inclusive, except for the
    last bin.

    Parameters
    ----------
    edges : array
        An array with items defining the edges.
    values: array
        An array of floats defining the values of each bin.

    Attributes
    ----------

    values : array
        The values contained in each bin.
    nr_bins : int
        the number of bins in the histogram.
    edges : array
        An array of bin edges.
    bin_width : array
        The width of each bin.
    """
    def __init__(self, edges, values):
        if len(edges) != len(values)+1:
            raise Exception("values need to be one shorter than the lin_edges")
        self.values = np.array(values)
        self.nr_bins = len(self.values)
        self.edges = np.array(edges)
        self.bin_width = np.diff(edges)

    def __str__(self):
        return str(self.values)

    def __len__(self):
        return len(self.values)

    def __getitem__(self, key):
        return self.values[key]

    def __setitem__(self, key, value):
        self.values[key] = value

    def __repr__(self):
        return self.values

    def __mul__(self, other):
        out = self.copy()
        out.values = self.values * other
        return out

    def __div__(self, other):
        out = self.copy()
        out.values = self.values / other
        return out

    def __truediv__(self, other):
        out = self.copy()
        out.values = self.values / other
        return out

    def __idiv__(self, other):
        out.values = self.values / other
        return self

    def __itruediv__(self, div):
        self.values = self.values / div
        return self

    def __imul__(self, mul):
        self.values = self.values * mul
        return self

    def copy(self):
        """Creates a copy of the hisogram

        Returns
        -------
        histogram
            An exact copy of the histogram

        """
        return histogram(self.edges, self.values)

    def get_bin(self, x):
        """Returns the bin number at value **x**

        Parameters
        ----------
        x : float
            value where you want to know the bin number

        Returns
        -------
        int
            The bin number

        """
        return kea.helpers._numba_get_bin(x,self.edges)

    def get_bin_center(self, i):
        """Returns the center of the bin

        Parameters
        ----------
        i : int
            Bin number

        Returns
        -------
        float
            The center of bin *i*

        """
        return (self.edges[i] + self.edges[i+1])/2


    def plot(self, *argv, **kwargs):
        """
        Plot the histogram. matplotlib.pyplot arguments can be passed on too
        """
        plt.hist(self.edges[:-1], self.edges, weights=self.values, histtype=u'step', *argv, **kwargs)
        return None

    def integral(self, x1, x2):
        """Returns the integral of the histogram between **x1** and **x2**.

        Parameters
        ----------
        x1 : float
            lower bound of the integration
        x2 : float
            upper bound of the integration

        Returns
        -------
        float
            The integral between **x1** and **x2**

        """
        return kea.helpers._numba_integral(x1, x2, self.edges, self.values, self.bin_width)




class histogram_BPASS(histogram):
    """ Container for the BPASS data to reside and make it possible to perform basic
    operations on the enclosed data.

    This class inherits all function of histogram, but overrides a few to allow
    for easier manipulation of linear and logarithmic data.

    The internal binning is linear in Gyr, but plotting can be done in linear and
    logarithmic binning.
    """
    def __init__(self, values):
        super().__init__(self.getLinEdges(), values)

    def integral(self, x1, x2):
        """ Returns the integral of the histogram between **x1** and **x2**.
            It return an integral in units :math:`events/yr/M_odot`

        Parameters
        ----------
        x1 : float
            lower bound of the integration
        x2 : float
            upper bound of the integration

        Returns
        -------
        float
            The integral between **x1** and **x2**

        """
        return histogram.integral(self, x1, x2) *1e9


    def plotLog(self, *argv, **kwargs):
        """ Plot the histogram on a logirthmic axis
        """
        _ = plt.hist(self.getLogEdges()[:-1], self.getLogEdges(), weights=self.values, histtype=u'step', *argv, **kwargs)
        _ = plt.xlabel("log(age/yr)")
        return None

    def plotLin(self, *argv, **kwargs):
        """Plot the histogram on a linear axis.
        """
        _ = plt.hist(self.getLinEdges()[:-1], self.getLinEdges(), weights=self.values, histtype=u'step', *argv, **kwargs)
        _ = plt.xlabel("age/Gyr")
        return None

    def getLogBins(self):
        """Returns the middle points of all bins, except the first bin, which
        runs from 0 to :math:`10^{6.05}`, but returns 6.

        Returns
        -------
        array
            The middle points of all bins in logarithmic space

        """
        return np.linspace(6,11,51)

    def getLinBins(self):
        """Returns the middle points of all bins, except the first bin, which
        runs from 0 to :math:`10^{6.05}``, but returns :math:`10^6`.

        Returns
        -------
        array
            The middle points of all bins in linear space

        """
        return 10**np.linspace(6, 11, 51)

    def getLogEdges(self):
        """Returns the edges in logarithmic space, except for the first bin,
        where the lowest edge is set to 5.95 instead of 0 for plotting purposes.

        Returns
        -------
        array
            The bin edges in logarithmic space

        """
        return np.linspace(5.95, 11.05, 52)


    def getLinEdges(self):
        """Returns the edges in linear space


        Returns
        -------
        array
            The bin edges in linear space

        """
        x = [0.0]
        x.extend(10**np.linspace(6.05, 11.05, 51)/1e9)
        return x

    def copy(self):
        """Create a copy of the BPASS histogram

        Returns
        -------
        BPASS_hist
            A copy of the BPASS historam.
        """

        return BPASS_hist(self.values)

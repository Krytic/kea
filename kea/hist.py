import numpy as np
import matplotlib.pyplot as plt


class histogram:
    """
    A class to make a simple histogram and fill it with data.
    It creates a histogram based on a lower and upper limit
    and a number of bins within those limits.

    The lowest bin and the highest bin are used to catch
    overflow or underflow values! Make sure your histogram
    covers your whole space.
    The upper edges are non inclusive.
    """
    def __init__(self,xlow=None, xup=None, nr_bins=None, edges=None):
        if xlow != None and xup != None and nr_bins != None:
            self._xlow = xlow
            self._xup = xup
            self._nr_bins = nr_bins
            self._bin_edges = np.linspace(xlow, xup, nr_bins+1)

        elif isinstance(edges, type([])) or isinstance(edges, type(np.array([]))):
            self._xlow = edges[0]
            self._xup = edges[-1]
            self._nr_bins = len(edges)-1
            self._bin_edges = np.array(edges)
        else:
            raise Exception("Not given the correct input")

        self._values = np.zeros(self._nr_bins)
        self.lower_edges = self._bin_edges[:-1]
        self.upper_edges = self._bin_edges[1:]

    def __len__(self):
        return len(self._values)

    def __str__(self):
        return str(self._values)

    def __repr__(self):
        return f"The bins: {self._bin_edges}\n The values: {self._values}"

    def __mul__(self, other):
        """
        Multiply the values with the given number or multiply each value with the
        value from the numpy array
        """
        out = self.copy()
        out._values = self._values * other
        return out

    def __div__(self, other):
        """
        Divides the values with the given number or divides each value with the
        value from the numpy array
        """
        out = self.copy()
        out._values = self._values / other
        return out

    def __truediv__(self, other):
        """
        Divides the values with the given number or divides each value with the
        value from the numpy array
        """
        out = self.copy()
        out._values = self._values / other
        return out

    def copy(self):
        out = histogram(xlow=self._xlow, xup=self._xup, nr_bins=self._nr_bins)
        out._values = self._values
        return out

    def Fill(self, x, w=1):
        """
        Add data to the histogram.
        x is the x-axis position. If no weight is given 1 is used.
        """
        def f(f, g):
            if f >= self._xup:
                 self._values[self._nr_bins-1] += g
            elif f <= self._xlow:
                self._values[0] += g
            else:
                bin_nr = np.where(self.lower_edges <= f)[0][-1]
                self._values[bin_nr] += g

        if not isinstance(x, type(0.0)):
            if w==1:
                for i in range(0, len(x)):
                    f(x[i], 1)
            elif len(x) != len(w):
                raise Exception("weights needs to be as long as x")
            else:
                for i in range(0, len(x)):
                    f(x[i], w[i])
        else:
            f(x, w)

    def plot(self, *argv, **kwargs):
        """
        Plot the histogram
        """
        _ = plt.hist(self._bin_edges[:-1], self._bin_edges, weights=self._values,histtype=u'step', *argv, **kwargs)
        return None

    def getBinContent(self, bin_nr):
        """
        Return the value of the given bin
        """
        return self._values[bin_nr]

    def getNBins(self):
        """
        Return the number of bins in the histogram
        """
        return self._nr_bins

    def getValues(self):
        """
        Return all the values
        """
        return self._values

    def getBinWidth(self, i):
        """
        return the bin width of the given bin
        """
        return self.upper_edges[i] - self.lower_edges[i]

    def getBinCenter(self, i):
        return (self.upper_edges[i] + self.lower_edges[i])/2

    def getBin(self, x):
        """
        Get the bin number of the data value
        """
        if x < self._bin_edges[0] or x > self._bin_edges[-1]:
            raise Exception("x outside of range")
        out = np.where(x >= self._bin_edges)[0][-1]
        if out == self._nr_bins:
            out = out-1
        return out

    def getBinEdges(self):
        return self._bin_edges


    def sum(self, x1, x2):
        """
        Perform a binwise summation of the values in between the given values
        """
        if x1 >= x2:
            raise Exception("x1 should be larger than x2")
        if x1 < self._bin_edges[0]:
            print("WARNING: lower limit is below lowest bin edge")
        if x2 > self._bin_edges[-1]:
            print("WARNING: higher limit is above the highest bin edge")
        lower_bin = self.getBin(x1)
        upper_bin = self.getBin(x2)
        if lower_bin == upper_bin:
            bin_width = self.getBinWidth(lower_bin)
            return self.getBinContent(lower_bin) * (x2 - x1) / bin_width
        else:
            total = 0
            # get lower bin part
            bin_width = self.getBinWidth(lower_bin)
            total += self.getBinContent(lower_bin) * (self.upper_edges[lower_bin] - x1)/bin_width

            # get upper bin part
            bin_width = self.getBinWidth(upper_bin)
            total += self.getBinContent(upper_bin) * (x2 - self.lower_edges[upper_bin])/bin_width

            # get the parts in between if they are there
            if (lower_bin + 1) != upper_bin:
                for i in range(lower_bin+1, upper_bin):
                    total += self._values[i]

            return total

    def integral(self, x1, x2):
        """
        Perform a binwise integration between x1 and x2.
        """
        if x1 >= x2:
            raise Exception("x1 should be larger than x2")
        if x1 < self._bin_edges[0]:
            print("WARNING: lower limit is below lowest bin edge")
        if x2 > self._bin_edges[-1]:
            print("WARNING: higher limit is above the highest bin edge")
        lower_bin = self.getBin(x1)
        upper_bin = self.getBin(x2)
        if lower_bin == upper_bin:
            bin_width = self.getBinWidth(lower_bin)
            return self.getBinContent(lower_bin) * (x2 - x1)
        else:
            total = 0
            # get lower bin part
            bin_width = self.getBinWidth(lower_bin)
            total += self.getBinContent(lower_bin) * (self.upper_edges[lower_bin] - x1)

            # get upper bin part
            bin_width = self.getBinWidth(upper_bin)
            total += self.getBinContent(upper_bin) * (x2 - self.lower_edges[upper_bin])

            # get the parts in between if they are there
            if (lower_bin + 1) != upper_bin:
                for i in range(lower_bin+1, upper_bin):
                    total += self._values[i] * self.getBinWidth(i)

            return total

class BPASS_hist(histogram):
    """
    A container for the BPASS data to reside and make it possible to perform basic
    operations on the enclosed data.

    The internal bining is linear, because it allows for easier comparison to other data.
    But logarithmic binning can be plotted and requested.
    """
    def __init__(self):
        super().__init__(edges=self.getLinEdges())

    def Fill(self, x, w=1,ty=None):
        """
        Add data to the BPASS plot.

        ty is "log" or "lin". If none is given, log is assumed.
        """

        if (ty == "log") or (ty == None):
            x = 10**x /1e9

        return histogram.Fill(self, x, w)

    def integral(self, x1, x2):
        """
        Performs the integration between the values and corrects for the binning
        being in Gyr instead of yrs.
        """
        return histogram.integral(self, x1, x2) *1e9


    def plotLog(self, *argv, **kwargs):
        """
        Plot the dat in a logarithmic way
        """
        _ = plt.hist(self.getLogEdges()[:-1], self.getLogEdges(), weights=self._values, histtype=u'step', *argv, **kwargs)
        _ = plt.xlabel("log(age/yr)")
        return None

    def plotLin(self, *argv, **kwargs):
        """
        Plot the data on a linear scale.
        """
        _ = plt.hist(self.getLinEdges()[:-1], self.getLinEdges(), weights=self._values, histtype=u'step', *argv, **kwargs)
        _ = plt.xlabel("age/Gyr")
        return None

    def getLogBins(self):
        """
        returns the middle point of the bin. Except for the first bin, which runs from 0 to 10**6.05,
        but return 6.
        """
        return np.linspace(6,11,51)

    def getLinBins(self):
        """
        return the middle point of the bin in linspace. Except for the first bin,
         which runs from 0 to 10**6.05,but return 10**6.
        """
        return 10**np.linspace(6, 11, 51)

    def getLogEdges(self):
        """
        Returns the bin edges in log space
        """
        return np.linspace(5.95, 11.05, 52)


    def getLinEdges(self):
        """
        returns the bin edges in linear space
        """
        x = [0.0]
        x.extend(10**np.linspace(6.05, 11.05, 51)/1e9)
        return x

    def copy(self):
        out = BPASS_hist()
        out._values = self._values
        return out

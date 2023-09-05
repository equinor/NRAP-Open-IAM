# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import h5py
import numpy as np
from numpy import ma
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from matplotlib import animation
from matk.ordereddict import OrderedDict
from copy import copy
import warnings
from math import atan2, degrees
import pandas as pd

np.set_printoptions(threshold=sys.maxsize)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class Mesh2D():
    def __init__(self, xs, ys, As=None):
        """
        Constructor method of MeshGrid object

        :param xs: Sorted, unique x values of mesh
        :type xs: array(float)
        :param ys: Sorted, unique y values of mesh
        :type ys: array(float)
        :param As: Cell areas, either a float if all areas are identical or
            a numpy meshgrid containing cell areas.
        :type As: float or array(float)
        :param sparse: If True a sparse grid is returned in order to conserve
            memory. Default is False.
        :type sparse : bool, optional
        :returns: NRAP-Open-IAM MeshGrid object
        """
        # TODO: We could automatically figure out cell areas.
        # Automatically figuring out interior cell areas would be easy,
        # however, boundary cell areas would require assumptions.
        self.xs = xs
        self.ys = ys
        self.xm, self.ym = np.meshgrid(xs, ys)
        # Transpose xm and ym matrices to have size (len(xs), len(ys))
        self.xm = self.xm.T
        self.ym = self.ym.T
        self.shape = self.xm.shape

        if isinstance(As, float):
            self.Am = np.ones(self.shape)*As
        elif not As is None and As.shape == self.shape:
            self.Am = As
        else: self.Am = None

        # Dictionary of variables keyed by name (e.g., "pressure")
        self.variables = {}

    def nans_like_data(self):
        """
        Return an array of nans with correct dimensions to add data with add_data method
        """
        return np.zeros(self.shape)*np.nan

    def zeros_like_data(self):
        """
        Return an array of zeros with correct dimensions to add data with add_data method
        """
        return np.zeros(self.shape)

    def add_variable(self, name):
        """
        Add a variable. Data associated with the variable at multiple time steps can then be added
        """
        self.variables[name] = Variable(name, self)
        return self.variables[name]

class Variable():
    def __init__(self, name, parent):
        """
        Constructor method of Variable object

        :param name: Variable name
        :type name: string
        :param parent: Parent object
        :type parent: Mesh2D object
        :returns: NRAP-Open-IAM Variable object
        """
        self.name = name
        self._parent = parent
        # Ordered dictionary of datasets keyed by time
        self.datasets = OrderedDict()
        #TODO: could add interpolator here. You would then need to specify
        # the time you want the interpolation.
        # self.interpolators = {}
    @property
    def xs(self):
        return self._parent.xs

    @property
    def ys(self):
        return self._parent.ys

    @property
    def xm(self):
        return self._parent.xm

    @property
    def ym(self):
        return self._parent.ym

    @property
    def Am(self):
        return self._parent.Am

    @property
    def shape(self):
        return self._parent.shape

    @property
    def times(self):
        return np.array(list(self.datasets.keys()))

    def nans_like_data(self):
        """
        Return an array of nans with correct dimensions to add data with add_data method
        """
        return np.zeros(self.shape)*np.nan

    def zeros_like_data(self):
        """
        Return an array of zeros with correct dimensions to add data with add_data method
        """
        return self._parent.zeros_like_data()

    def add_dataset(self, data, time):
        """
        Add a data set

        :param data: Meshgrid of data to add to Variable object
        :type data: Numpy Meshgrid
        :param time: Time associated with data
        :type time: float
        """
        if not data.shape == self.shape:
            print("Error: Shape of data must be {}".format(self.shape))
        self.datasets[time] = DataSet(data, time, self)

    def plume_areas(self, thresh):
        '''
        Calculate plume area

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of plume areas
        '''
        return np.array([d.plume_area(thresh) for d in list(self.datasets.values())])

    def plume_areas_dt(self, thresh):
        '''
        Calculate the derivative of plume area with respect to time

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of plume areas
        '''
        return np.gradient(
            np.array([d.plume_area(thresh) for d in list(self.datasets.values())]),
            self.times)

    def plume_centroids(self, thresh):
        '''
        Calculate plume centroid (https://en.wikipedia.org/wiki/Image_moment).

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of plume centroids (x,y pairs)
        '''
        return np.array([d.plume_centroid(thresh) for d in list(self.datasets.values())])

    def plume_centroids_dt(self, thresh):
        '''
        Calculate the derivative of plume centroid with respect to time
        (https://en.wikipedia.org/wiki/Image_moment).

        :param thresh: Cutoff value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of the derivative of plume centroids with time
            (x,y pairs)
        '''
        dcdt = np.gradient(
            np.array([d.plume_centroid(thresh) for d in list(self.datasets.values())]),
            self.times, axis=0)
        # Convert nans to zeros
        dcdt[np.isnan(dcdt)] = 0
        return dcdt

    def mobility(self, thresh):
        '''
        Calculate the mobility of the plume.

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy arrays of plume speed, direction relative to positive x axis
        '''
        return self.effective_velocity(thresh)

    def effective_velocity(self, thresh):
        '''
        Calculate the effective velocity of the plume.

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy arrays of plume speed, direction relative to positive x axis
        '''
        # Calculate x and y centroid diffs between times
        diffs = np.diff(
            np.array([d.plume_centroid(thresh) for d in list(self.datasets.values())]),
            axis=0)
        # Calculate vector angle from positive x axis
        angles = [atan2(y, x)*180/np.pi for x, y in diffs]
        # Calculate squared values diffs
        sqrd_diffs = np.power(diffs, 2)
        # Calculate hypotenuse of x and y diffs (distance traveled by plume)
        lngths = np.sqrt(np.add(sqrd_diffs[:, 0], sqrd_diffs[:, 1]))
        lngths[np.isnan(lngths)] = 0
        dldt = np.zeros_like(lngths)
        for i, dl, dt in zip(range(len(lngths)), lngths, np.diff(self.times)):
            if dt > 0.:
                dldt[i] = dl / dt
        # Add zero for t=0
        dldt = np.insert(dldt, [0], 0.)
        angles = np.insert(angles, [0], 0.)
        return dldt, angles

    def spreading(self, thresh):
        # Calculate spreading in primary and secondary directions
        dldt = self.plume_spreads_dt(thresh)
        # Calclulate primary and secondary angles
        angles = self.plume_angles(thresh)
        # Find indices for max spreading coefficient
        maxind = np.argmax(np.abs(dldt), axis=1)
        # Return max spreading values and associated angles
        return np.array([dldt[i, v] for i, v in enumerate(maxind)]), np.array([
            angles[i, v] for i, v in enumerate(maxind)])

    def plume_spreads(self, thresh):
        '''
        Calculate plume spread along characteristic axes
        (https://en.wikipedia.org/wiki/Image_moment).

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of plume spreads along characteristic (eigenvector) axes
        '''
        return np.array([d.plume_spread(thresh) for d in list(self.datasets.values())])

    def plume_spreads_dt(self, thresh):
        '''
        Calculate the derivative of plume spread along characteristic axes
        with respect to time (https://en.wikipedia.org/wiki/Image_moment).

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of the derivative of plume spreads
            along characteristic (eigenvector) axes
        '''
        dsdt = np.gradient(
            np.array([d.plume_spread(thresh) for d in list(self.datasets.values())]),
            self.times, axis=0)
        # Convert nans to zeros, nans occur if there is no plume
        dsdt[np.isnan(dsdt)] = 0
        # Divide by 2 so that this is the effective longitudinal dispersion
        # coefficient (Valochhi, Water Resources Research, 1989)
        dsdt /= 2.
        return dsdt

    def plume_eigenvectors(self, thresh):
        '''
        Calculate plume characteristic axes as eigenvectors

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of plume characteristic axes (eigenvectors)
            where columns contain eigenvectors
        '''
        return np.array([d.plume_eigenvector(thresh) for d in list(self.datasets.values())])

    def plume_angles(self, thresh):
        '''
        Calculate plume moment angles from the positive x-axis in degrees

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Numpy array of plume moment angles from the positive x-axis
        '''
        return np.array([d.plume_angle(thresh) for d in list(self.datasets.values())])

    def plume_stability_animation(self, thresh, animation_writers='ffmpeg', filename=None):
        '''
        Create plume stability animation

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :param animation_writers: Backend to write animations with
        :type animation_writers: string
        :param filename: Name of file to save animation to (e.g., "plume_stability.mp4");
            if None, animation will be displayed
        :type filename: string
        :param return_object: return matplotlib animation object
        :type filename: matplotlib animation object
        '''
        Writer = animation.writers[animation_writers]
        writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)

        ts = self.times
        # Color maps
        fig = plt.figure(figsize=(8, 10))
        # Plot plume area
        ax1 = plt.subplot2grid((6, 4), (0, 0), colspan=4)
        ax1.plot(ts, self.plume_areas(thresh))
        ax1.set_ylabel(r'Area [L$^2$]')
        line1 = ax1.axvline(0, ls='--', c='k')
        # Plot plume area derivative
        ax2 = plt.subplot2grid((6, 4), (1, 0), colspan=4)
        ax2.plot(ts, self.plume_areas_dt(thresh), label='Pressure')
        line2 = ax2.axvline(0, ls='--', c='k')
        ax2.set_ylabel(r'dA/dt [L$^2$/T]')
        # Plot plume centroid derivative
        ax3 = plt.subplot2grid((6, 4), (2, 0), colspan=4)
        #dcdt = self.plume_centroids_dt(thresh)
        dcdt, cangles = self.effective_velocity(thresh)
        ax3.plot(ts, dcdt, label='Speed')
        ax3_2 = ax3.twinx()
        ax3_2.plot(ts, cangles, '.', c='grey', label='angle')
        ax3_2.set_ylim(-180, 180)
        ax3_2.set_ylabel(r'Direction [degrees]')
        line3 = ax3.axvline(0, ls='--', c='k')
        ax3.set_ylabel(r'Mobility [L/T]')
        ax3.legend()
        # Plot plume spread derivative
        ax4 = plt.subplot2grid((6, 4), (3, 0), colspan=4)
        dsdt = self.plume_spreads_dt(thresh)
        dsdt_angles = self.plume_angles(thresh)
        for i in range(len(dsdt_angles[:, 0])):
            if dsdt_angles[i, 0] < 0:
                dsdt_angles[i, 0] += 180.
        ax4.plot(ts, dsdt[:, 0], label=r'Spreading')
        line4 = ax4.axvline(0, ls='--', c='k')
        ax4.set_xlabel('Time')
        ax4.set_ylabel(r'Spreading [L$^2$/T]')
        ax4.legend()
        ax4_2 = ax4.twinx()
        ax4_2.plot(ts, dsdt_angles[:, 0], '.', c='grey', label='angle')
        ax4_2.set_ylim(0., 180)
        ax4_2.set_ylabel(r'Direction [degrees]')
        # Plot colormap of plume
        ax5 = plt.subplot2grid((6, 4), (4, 1), colspan=2, rowspan=2)
        extent = [np.nanmin(self.xm), np.nanmax(self.xm),
                  np.nanmin(self.ym), np.nanmax(self.ym)]
        vmin = 1e20
        vmax = -1.0e20

        for d in list(self.datasets.values()):
            vmin = np.nanmin([vmin, np.nanmin(d.data)])
            vmax = np.nanmax([vmax, np.nanmax(d.data)])
        im = ax5.imshow(np.flipud(self.datasets[self.times[0]].data.T),
                        animated=True, vmin=vmin, vmax=vmax, extent=extent)
        txt = plt.text(0.1, 0.3, 'Time: 0.0', fontsize=14, transform=plt.gcf().transFigure)
        fig.colorbar(im)
        ax5.set_xlabel('Axis1 [L]')
        ax5.set_ylabel('Axis2 [L]')
        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_xlim((0, self.times[-1]))
            ax.axhline(0, c='gray', alpha=0.5)
        def updatefig(i):
            im.set_array(np.flipud(self.datasets[i].data.T))
            txt.set_text('Time: '+str(i))
            line1.set_xdata([i, i])
            line2.set_xdata([i, i])
            line3.set_xdata([i, i])
            line4.set_xdata([i, i])
            return im,
        self.ani = animation.FuncAnimation(fig, updatefig, frames=self.times,
                                      interval=500, blit=True, repeat=False)
        fig.tight_layout()
        if not filename is None:
            ani.save(filename,writer=writer)
        if return_object:
            return ani
        # Put this down here so that both file can be created and ani object
        # can be returned in one call. Otherwise, just show the animation
        if filename is None:
            fig.show()
        else:
            self.ani.save(filename, writer=writer)

    def interpolate(self, pts):
        '''
        Interpolate dataset at given points for all times

        :param pts: 2D Points at which to return interpolated values
        :type pts: array(tuple(float,float))
        :returns: Array of interpolated values
        '''
        return np.array([d.interpolate(pts) for d in list(self.datasets.values())])

class DataSet():
    def __init__(self, data, time, parent):
        """
        Constructor for DataSet object

        :param data: Meshgrid of data
        :type data: Numpy Meshgrid
        :param time: Time associated with data
        :type time: float
        :param parent: Parent object
        :type parent: NRAP-Open-IAM Variable object
         """
        if not data.shape == parent.shape:
            print("Error: Shape of data must be {}".format(parent.shape))
        self.data = ma.masked_invalid(data)
        self.time = time
        # Parent (DataSet)
        self._parent = parent
        # Grandparent (Mesh2D)
        self._gparent = parent._parent
        self.interpolator = None

    @property
    def xs(self):
        return self._gparent.xs

    @property
    def ys(self):
        return self._gparent.ys

    @property
    def xm(self):
        return self._gparent.xm

    @property
    def ym(self):
        return self._gparent.ym

    @property
    def Am(self):
        return self._gparent.Am

    @property
    def shape(self):
        return self._gparent.shape

    def _create_interpolator(self):
        self.interpolator = RegularGridInterpolator((self.xs, self.ys), self.data.transpose())

    def interpolate(self, pts):
        """
        Interpolate data at 2D locations

        :param pts: 2D Points at which to return interpolated values
        :type pts: array(tuple(float,float))
        :returns: Array of interpolated values
        """
        if self.interpolator is None:
            self._create_interpolator()
        return self.interpolator(pts)

    def plume_area(self, thresh):
        '''
        Calculate plume area

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: Plume area
        '''
        # Create a binary mesh grid with ones where vs is greater than thresh, zero otherwise
        vsb = np.zeros_like(self.data)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="invalid value encountered in greater")
            vsb[self.data > thresh] = 1
        # Sum up the binary meshgrid times the meshgrid of areas to get the plume area
        return np.nansum(vsb*self.Am)

    def plume_centroid(self, thresh):
        '''
        Calculate plume centroid (https://en.wikipedia.org/wiki/Image_moment).

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: plume centroid (x,y) value
        '''
        # Copy data
        vs = copy(self.data)
        # Set all value below thresh to zero, ignore annoying warning
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="invalid value encountered in less")
            vs[vs < thresh] = 0
        # Calculate raw moments
        m00 = np.nansum(vs)
        m10 = np.nansum(self.xm*vs)
        m01 = np.nansum(self.ym*vs)
        # Calculate plume centroid
        if m00 == 0: # No plume exists in this case
            xc = yc = np.nan
        else:
            xc = m10/m00
            yc = m01/m00
        return xc, yc

    def plume_covariance(self, thresh):
        '''
        Calculate plume moment covariance (https://en.wikipedia.org/wiki/Image_moment).

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: plume moment covariance or None if undefined (e.g., nonexistent plume)
        '''
        # Copy data
        vs = copy(self.data)
        # Set all value below thresh to zero, ignore annoying warning
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="invalid value encountered in less")
            vs[vs < thresh] = 0
        # Calculate raw moments
        m00 = np.nansum(vs)
        m11 = np.nansum(self.xm*self.ym*vs)
        m20 = np.nansum(self.xm**2*vs)
        m02 = np.nansum(self.ym**2*vs)

        if m00 == 0: # No plume exists in this case
            return None

        # Calculate plume centroid
        xc, yc = self.plume_centroid(thresh)

        # Calculate moment covariances
        cov20 = m20/m00-xc**2
        cov02 = m02/m00-yc**2
        cov11 = m11/m00-xc*yc

        # Change to numpy eigenanalysis so that eigenvectors get calculated as well
        #lm1 = ((cov20+cov02)/2 + np.sqrt(4*cov11**2 + (cov20-cov02)**2)/2)
        #lm2 = ((cov20+cov02)/2 - np.sqrt(4*cov11**2 + (cov20-cov02)**2)/2)
        # Construct covariance matrix
        cov = np.array([[cov20, cov11], [cov11, cov02]])
        return cov

    def plume_spread(self, thresh):
        '''
        Calculate plume spread along characteristic axes

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: plume spread along characteristic (eigenvector) axes
        '''
        # Calculate eigenvalues (w) and eigenvectors (columns of v)
        cov = self.plume_covariance(thresh)
        if cov is None:
            return (np.nan, np.nan)

        w, _ = np.linalg.eig(cov)
        # Return primary eigenvalue first
        if np.abs(w[0]) > np.abs(w[1]):
            return w

        return np.flip(w, axis=0)

    def plume_eigenvector(self, thresh):
        '''
        Calculate plume characteristic axes vectors (eigenvectors)

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: matrix where characteristic axes (eigenvectors) are along the columns
        '''
        # Calculate eigenvalues (w) and eigenvectors (columns of v)
        cov = self.plume_covariance(thresh)
        if cov is None:
            return np.array([[np.nan, np.nan], [np.nan, np.nan]])

        w, v = np.linalg.eig(cov)
        # Return eigenvector associated with primary eigenvalue first
        if np.abs(w[0]) > np.abs(w[1]):
            return v

        return np.flip(v, axis=1)

    def plume_angle(self, thresh):
        '''
        Calculate plume moment angles from the positive x-axis in degrees.
        The primary characteristic axis angle is returned first followed by the secondary

        :param thresh: Cutoff off value for inclusion in plume
        :type thresh: float
        :returns: tuple of angles started with the primary characteristic axis
        '''
        # Calculate eigenvalues (w) and eigenvectors (columns of v)
        cov = self.plume_covariance(thresh)
        if cov is None:
            return np.nan, np.nan

        w, v = np.linalg.eig(cov)
        # Return angle associated with primary eigenvalue first
        if np.abs(w[0]) > np.abs(w[1]):
            return degrees(atan2(v[1, 0], v[0, 0])), degrees(atan2(v[1, 1], v[0, 1]))

        return degrees(atan2(v[1, 1], v[0, 1])), degrees(atan2(v[1, 0], v[0, 0]))

def read_Mesh2D_data(data_file, variable_names, time_points, is_hdf5_file):
    '''
    Read in NRAP-Open-IAM formatted data file into an NRAP-Open-IAM Mesh2D object

    :param data_file: Name of data file
    :type data_file: str
    :param variable_names: Array of names of data to be extracted from data file
    :type variable_names: array-like
    :param time_points: Array of time points at which data is to be extracted
    :type time_points: array-like
    :param is_hdf5_file: Flag indicating whether data file has hdf5 or h5 format
    :type is_hdf5_file: boolean
    :returns: NRAP-Open-IAM Mesh2D object
    '''
    # Check whether data file is of hdf5/h5 format
    if is_hdf5_file:
        d = {}
        with h5py.File(data_file, 'r') as hf:
            all_names = list(hf.keys())
            for nm in ['x', 'y', 'z', 'area']:
                try:
                    d[nm] = hf[nm][()]
                except KeyError:
                    pass

            for nm in variable_names:
                for ind in range(len(time_points)):
                    f_nm = nm + '_' + str(ind + 1)
                    d[f_nm] = hf[f_nm][()]

            # Rewrite all_names to new order for d
            all_names = list(d.keys())
    else:
        d = np.genfromtxt(data_file, names=True, delimiter=',')
        all_names = list(d.dtype.names)

    # Collect unique values of x and y
    xu = np.unique(d['x'])
    yu = np.unique(d['y'])

    # If areas are provided, collect them
    if 'area' in all_names:
        dAc = d['area']
    else:
        # Otherwise, attempt to estimate the area assuming it is a uniform,
        # orthogonal mesh
        # Get unique x and y coordinates, need to round to nearest 10 meters
        # due to stacked mesh offsets
        dxc = np.mean(np.diff(xu))
        dyc = np.mean(np.diff(yu))
        dAc = np.ones_like(d['x'])*dxc*dyc # horizontal area of cells

    # Create meshgrid of zeros with correct dimensions to store cell areas below
    As = np.zeros_like(np.meshgrid(xu, yu)[0])

    # Create Mesh2D object
    M = Mesh2D(xu, yu)

    # Create xu and yu dictionaries to aid in formation of variable meshgrids
    xud = {}
    yud = {}
    for i, x in enumerate(xu):
        xud[str(x)] = i
    for i, y in enumerate(yu):
        yud[str(y)] = i

    # Check to see if data is 3D, and if so, collapse to 2D
    if 'z' in all_names:
        d = collapse_3D_data(d, all_names)

    # Add data to variable objects
    for vnm in variable_names:
        M.add_variable(vnm)
        for ind, t in enumerate(time_points):
            V = M.nans_like_data()
            for x, y, v in zip(d['x'], d['y'], d[vnm+'_'+str(ind+1)]):
                V[xud[str(x)], yud[str(y)]] = v
            M.variables[vnm].add_dataset(V, t)
            del V

    # Add cell areas
    M.Am = M.nans_like_data()
    for x, y, v in zip(d['x'], d['y'], dAc):
        M.Am[xud[str(x)], yud[str(y)]] = v

    return M

def collapse_3D_data(d, all_names):
    """
    Squash 3d data into 2d along z-axis by selecting maximum values.
    """
    # Convert to pandas dataframe
    if isinstance(d, np.ndarray):
        d_p = pd.DataFrame(d, columns=all_names)
    elif isinstance(d, dict):
        d_p = pd.DataFrame.from_dict(d, orient='columns')

    # Get all xy coordinate pairs
    coords = [[x, y] for x, y in list(zip(d['x'], d['y']))]

    # Get unique xy coordinate pairs
    uniq_coords = np.unique(coords, axis=0)

    # Get number of unique coordinates
    coord_num = len(uniq_coords)

    # Pre-allocate pandas dataframe
    c = pd.DataFrame(columns=all_names, index=range(coord_num))

    # For each coordinate pair, get max values in z and write to dataframe
    for i, coord in enumerate(uniq_coords):
        c_temp = d_p[(d_p['x']==coord[0]) & (d_p['y']==coord[1])].max()
        c.iloc[i] = c_temp

    # Remove z-value information
    c = c.drop('z', axis='columns')

    # Convert pandas dataframe to dictionary
    d = c.to_dict('list')

    return d

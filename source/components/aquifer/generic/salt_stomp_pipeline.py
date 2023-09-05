# Convert STOMP plot files into tensorflow datasets
# Tested using python 3.7
import glob
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.base import TransformerMixin
from sklearn.preprocessing import (StandardScaler, MinMaxScaler,
                                   PowerTransformer, QuantileTransformer)
from plot_to_numpy import STOMP_Run

class Run(STOMP_Run):

   # get permeability at all time steps
    def get_x_permeabilities(self):
        x = self.get_variable('X-Direction Intrinsic Permeability',
                              func=lambda x: np.log10(x))
        return list(x.values())

   # get permeability at all time steps
    def get_z_permeabilities(self):
        x = self.get_variable('Z-Direction Intrinsic Permeability',
                              func=lambda x: np.log10(x))
        return list(x.values())

    # get porosity at all time steps
    def get_porosities(self):
        x = self.get_variable('Diffusive Porosity')
        return list(x.values())

    def get_z_coordinates(self):
        x = self.get_centroids('Z')
        return list(x.values())

    def get_x_coordinates(self):
        x = self.get_centroids('X')
        return list(x.values())

    # get injection rate at all time steps
    def get_co2_mass_source_rates(self):
        x = self.get_variable(name='CO2 Mass Source Rate',
                              func=lambda x: np.log1p(x * 360))
        return list(x.values())

    def get_water_mass_source_rates(self):
        x = self.get_variable(name='Water Mass Source Rate',
                              func=lambda x: np.log1p(x * 360))
        return list(x.values())

    def get_salt_mass_source_rates(self):
        x = self.get_variable(name='Salt Mass Source Rate',
                              func=lambda x: np.log1p(x * 360))
        return list(x.values())

    # get cumulative mass injected at all time steps
    def get_co2_mass_source_integrals(self):
        x = self.get_variable_change(name='CO2 Mass Source Integral',
                                     base_time=100.0, func=lambda x: np.log1p(x * 360))
        return list(x.values())

    def get_water_mass_source_integrals(self):
        x = self.get_variable_change(name='Water Mass Source Integral',
                                     base_time=100.0, func=lambda x: np.log1p(x * 360))
        return list(x.values())

    def get_salt_mass_source_integrals(self):
        x = self.get_variable_change(name='Salt Mass Source Integral',
                                     base_time=100.0, func=lambda x: np.log1p(x * 360))
        return list(x.values())

    # get pressure at all time steps
    def get_pressures(self):
        x = self.get_variable('Gas Pressure')
        return list(x.values())

    # get saturation at all time steps
    def get_saturations(self):
        x = self.get_variable('Gas Saturation')
        return list(x.values())

    # get pH at all time steps
    def get_pH(self):
        x = self.get_variable(name='Aqueous h+ Concentration',
                              func=lambda x: -np.log10(x))
        return list(x.values())

    # get TDS at all time steps
    def get_TDS(self):
        salt = self.get_variable('Aqueous Salt Mass Fraction')
        dens = self.get_variable(name='Aqueous Density', func=lambda x: x * 1000)
        x = {key : value * dens[key] for key, value in salt.items() if key in dens}
        return list(x.values())

    # get aqueous CO2 mass fraction
    def get_aq_co2(self):
        x = self.get_variable('Aqueous CO2 Mass Fraction',
                              func=lambda x: np.log10(x + 1.e-30))
        return list(x.values())

    # get aqueous salt mass fraction
    def get_aq_salt(self):
        x = self.get_variable('Aqueous Salt Mass Fraction',
                              func=lambda x: np.log10(x + 1.e-10))
        return list(x.values())

    # get aqueous salt mass fraction
    def get_init_salt(self):
        x = self.get_variable('Aqueous Salt Mass Fraction',
                              func=lambda x: np.log10(x + 1.e-10))
        # for key in x:
        #     print('key',key)
        return [x[100.0] for key in x]

class Images:
# convert CMG GEM simulaton data into images

    def __init__(self, path_to_files, input_transformer=None, output_transformer=None):

        self.path_to_files = path_to_files
        self.input_transformer = input_transformer
        self.output_transformer = output_transformer
        self.image_width = None
        self.image_height = None
        self.image_depth = None

        self.load()
        self.resize()
        self.transform()

    def load(self):
        # open all the run files
        paths = glob.glob(self.path_to_files)
        runs = [Run(path) for path in paths]
        self.image_width = runs[0].nx
        self.image_height = runs[0].ny
        self.image_depth = runs[0].nz

        # read data from files
        time_list = [run.get_times() for run in runs]
        time_series = np.vstack(time_list)

        xperm_list = [run.get_x_permeabilities() for run in runs]
        xperm_array = np.vstack(xperm_list)

        zperm_list = [run.get_z_permeabilities() for run in runs]
        zperm_array = np.vstack(zperm_list)

        xcoord_list = [run.get_x_coordinates() for run in runs]
        xcoord_array = np.vstack(xcoord_list)

        zcoord_list = [run.get_z_coordinates() for run in runs]
        zcoord_array = np.vstack(zcoord_list)

        por_list = [run.get_porosities() for run in runs]
        por_array = np.vstack(por_list)

        # waterrate_list = [run.get_water_mass_source_rates() for run in runs]
        # waterrate_array = np.vstack(waterrate_list)

        watermass_list = [run.get_water_mass_source_integrals() for run in runs]
        watermass_array = np.vstack(watermass_list)

        # co2rate_list = [run.get_co2_mass_source_rates() for run in runs]
        # co2rate_array = np.vstack(co2rate_list)

        co2mass_list = [run.get_co2_mass_source_integrals() for run in runs]
        co2mass_array = np.vstack(co2mass_list)

        # saltrate_list = [run.get_salt_mass_source_rates() for run in runs]
        # saltrate_array = np.vstack(saltrate_list)

        saltmass_list = [run.get_salt_mass_source_integrals() for run in runs]
        saltmass_array = np.vstack(saltmass_list)

        aqco2_list = [run.get_aq_co2() for run in runs]
        aqco2_array = np.vstack(aqco2_list)

        aqsalt_list = [run.get_aq_salt() for run in runs]
        aqsalt_array = np.vstack(aqsalt_list)

        initsalt_list = [run.get_init_salt() for run in runs]
        initsalt_array = np.vstack(initsalt_list)

        # convert values into arrays
        def mapit(values, shape, xcoord=None, ycoord=None, zcoord=None):
            # entire array has same value
            if xcoord is None and ycoord is None and zcoord is None:
                new_array = np.ones(shape)
                for (item, value) in zip(new_array, values):
                    item *= value
            # assign value to single point
            else:
                new_array = np.zeros(shape)
                for (item, value) in zip(new_array, values):
                    item[xcoord-1][ycoord-1][zcoord-1] += value
            return new_array

        # give new arrays as grid
        new_shape = por_array.shape
        # print(new_shape)
        time_array = mapit(time_series, new_shape)

        # combine into input and output images
        # print('x',xcoord_array.shape)
        # print('z',zcoord_array.shape)
        # print('kx',xperm_array.shape)
        # print('kz',zperm_array.shape)
        # print('por',por_array.shape)
        # print('init_salt',initsalt_array.shape)
        # print('co2_mass',co2mass_array.shape)
        # print('watermass',watermass_array.shape)
        # print('saltmass',saltmass_array.shape)
        # print('time',time_array.shape)
        self.input_images = np.stack((
            xcoord_array, zcoord_array, xperm_array, zperm_array, por_array,
            initsalt_array, co2mass_array, watermass_array, saltmass_array, time_array),
            axis=4)
        #self.output_images = np.expand_dims(aqco2_array, axis=4)
        self.output_images = np.expand_dims(aqsalt_array, axis=4)

        # remember min and max values for plotting
        self.input_mins = [xcoord_array.min(), zcoord_array.min(),
                           xperm_array.min(), zperm_array.min(),
                           por_array.min(), initsalt_array.min(),
                           co2mass_array.min(), watermass_array.min(),
                           saltmass_array.min(), time_array.min()]
        self.input_maxs = [xcoord_array.max(), zcoord_array.max(),
                           xperm_array.max(), zperm_array.max(),
                           por_array.max(), initsalt_array.max(),
                           co2mass_array.max(), watermass_array.min(),
                           saltmass_array.max(), time_array.max()]
        #self.output_mins = [aqco2_array.min()]
        #self.output_maxs = [aqco2_array.max()]
        self.output_mins = [aqsalt_array.min()]
        self.output_maxs = [aqsalt_array.max()]
    def transform(self):
        # make a new transformer if one isn't specified
        if not self.input_transformer:
            # self.input_transformer = NDQuantileTransformer(
            #     output_distribution='normal', random_state=0).fit(self.input_resized)
            self.input_transformer = NDMinMaxScaler().fit(self.input_resized)
        if not self.output_transformer:
            # self.output_transformer = NDQuantileTransformer(
            #     output_distribution='normal', random_state=0).fit(self.output_resized)
            self.output_transformer = NDMinMaxScaler().fit(self.output_resized)

        self.input_trans = self.input_transformer.transform(self.input_resized)
        self.output_trans = self.output_transformer.transform(self.output_resized)

    def resize(self):
        # resize images and convert to tensors
        def resize_image(input_image, output_image, height, width):
            input_image = tf.image.resize(input_image, [height, width],
                method=tf.image.ResizeMethod.NEAREST_NEIGHBOR)
            output_image = tf.image.resize(output_image, [height, width],
                method=tf.image.ResizeMethod.NEAREST_NEIGHBOR)
            return input_image, output_image

        def convert_image(input_image, output_image):
            input_image = tf.convert_to_tensor(input_image, dtype=tf.float32)
            output_image = tf.convert_to_tensor(output_image, dtype=tf.float32)
            # input_image, output_image = resize_image(input_image, output_image,
            # self.image_height, self.image_width)
            return input_image, output_image

        self.input_resized, self.output_resized = convert_image(
            self.input_images, self.output_images)

    def compare(self):
        # compare original and transformed data

        def histograms(title, t1, t2):
            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.suptitle(title)
            ax1.hist(tf.reshape(t1, [-1]), bins=20, log=True)
            ax2.hist(tf.reshape(t2, [-1]), bins=20, log=True)

        histograms('perm', self.input_resized[:, :, :, 0], self.input_trans[:, :, :, 0])
        histograms('inj', self.input_resized[:, :, :, 1], self.input_trans[:, :, :, 1])
        histograms('time', self.input_resized[:, :, :, 2], self.input_trans[:, :, :, 2])
        histograms('pres', self.output_resized[:, :, :, 0], self.output_trans[:, :, :, 0])
        histograms('sat', self.output_resized[:, :, :, 1], self.output_trans[:, :, :, 1])
        histograms('prod', self.output_resized[:, :, :, 2], self.output_trans[:, :, :, 2])
        plt.show()

    # plot a single set of input images
    def plot_input(self, i):
        pass

# Apply scikit-learn scalers to N-dimensional arrays by flattening them first

class NDStandardScaler(TransformerMixin):
    # N-dimensional StandardScaler
    def __init__(self, **kwargs):
        self._scaler = StandardScaler(copy=True, **kwargs)
        self._orig_shape = None

    def fit(self, X, **kwargs):
        X = np.array(X)
        # Save the original shape to reshape the flattened X later
        # back to its original shape
        if len(X.shape) > 1:
            self._orig_shape = X.shape[1:]
        X = self._flatten(X)
        self._scaler.fit(X, **kwargs)
        self.data_min_ = self._scaler.data_min_
        self.data_max_ = self._scaler.data_max_
        return self

    def transform(self, X, **kwargs):
        X = np.array(X)
        X = self._flatten(X)
        X = self._scaler.transform(X, **kwargs)
        X = self._reshape(X)
        return X

    def inverse_transform(self, X, **kwargs):
        X = np.array(X)
        X = self._flatten(X)
        X = self._scaler.inverse_transform(X, **kwargs)
        X = self._reshape(X)
        return X

    def _flatten(self, X):
        # Reshape X to <= 2 dimensions
        if len(X.shape) > 2:
            # flatten from right to left
            # n_dims = np.prod(self._orig_shape)
            # X = X.reshape(-1, n_dims)
            # flatten from left to right
            X = X.reshape(-1, X.shape[-1])
        return X

    def _reshape(self, X):
        # Reshape X back to it's original shape
        if len(X.shape) >= 2:
            X = X.reshape(-1, *self._orig_shape)
        return X

class NDMinMaxScaler(NDStandardScaler):
    # N-dimensional MimMaxScaler
    def __init__(self, **kwargs):
        self._scaler = MinMaxScaler(copy=True, **kwargs)
        self._orig_shape = None

class NDMaxAbsScaler(NDStandardScaler):
    # N-dimensional MimMaxScaler
    def __init__(self, **kwargs):
        self._scaler = MaxAbsScaler(copy=True, **kwargs)
        self._orig_shape = None

class NDPowerTransformer(NDStandardScaler):
    # N-dimensional PowerTransformer
    def __init__(self, **kwargs):
        self._scaler = PowerTransformer(copy=True, **kwargs)
        self._orig_shape = None

class NDQuantileTransformer(NDStandardScaler):
    # N-dimensional QuantileTransformer
    def __init__(self, **kwargs):
        self._scaler = QuantileTransformer(copy=True, **kwargs)
        self._orig_shape = None

if __name__ == '__main__':

    # PIPELINE TESTS
    # train = Images('/Volumes/d3a926/nrap/generic_aquifer/runs/48_runs/run1')
    test = Images('/rcfs/projects/nrap/OpenIAM/generic_aquifer_salinity/runs/5000_runs/run1')
    # test = Images('/Volumes/d3a926/nrap/generic_aquifer/runs/5_runs/run1',
    #               train.input_transformer, train.output_transformer)

    # print(train.input_resized.shape)
    print(test.input_resized.shape)

    print(test.input_mins, test.input_maxs)
    print(test.output_mins, test.output_maxs)

    # SAVE FOR LATER - MAY BE NEEDED FOR OTHER DATASETS
    # split images into train, test and validate sets
    # from sklearn.model_selection import train_test_split
    # input_images_train, input_images_test, output_images_train, output_images_test = \
    #     train_test_split(input_images, output_images, test_size=0.30, random_state=42)
    # input_images_valid, input_images_test, output_images_valid, output_images_test = \
    #     train_test_split(input_images_test, output_images_test, test_size=0.30, random_state=42)

    # train_dataset = tf.data.Dataset.from_tensor_slices(
    #     (train.input_trans, train.output_trans))
    # test_dataset = tf.data.Dataset.from_tensor_slices(
    #     (test.input_resized, test.output_resized))

    # print(train_dataset)
    # print(test_dataset)

    # train.compare()

    # SCALER TESTS

    # data = [[[0., 1.], [2., 3.]], [[1., 5.], [2., 9.]]]
    # print(np.array(data))
    # scaler = NDMinMaxScaler()
    # print(scaler.fit_transform(data))

    # X_lognormal = np.random.RandomState(616).lognormal(size=(25, 25, 27))
    # print('X_lognormal')
    # print(X_lognormal)
    # print(X_lognormal.min(), X_lognormal.max())

    # pt_scaler = NDPowerTransformer(method='box-cox', standardize=False).fit(X_lognormal)
    # X_power = pt_scaler.transform(X_lognormal)
    # print(X_power)
    # print(X_power.min(), X_power.max())
    # X_inverse = pt_scaler.inverse_transform(X_power)
    # print(X_inverse)
    # print(X_inverse.min(), X_inverse.max())

    # qt_scaler = NDQuantileTransformer(
    #     output_distribution='normal', random_state=0).fit(X_lognormal)
    # X_quantile = qt_scaler.transform(X_lognormal)
    # print('X_quantile')
    # print(X_quantile)
    # print(X_quantile.min(), X_quantile.max())
    # X_inverse = qt_scaler.inverse_transform(X_quantile)
    # print('X_inverse')
    # print(X_inverse)
    # print(X_inverse.min(), X_inverse.max())

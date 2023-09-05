import numpy as np
import tensorflow as tf
import joblib
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import os
from model_input import generate_input
import warnings

np.random.seed(1)
tf.random.set_seed(2)

current_directory = os.path.dirname(os.path.abspath(__file__))

def model(inputArray):

    input_images = generate_input(inputArray)

    # load correct model and data data transformers based on depth
    depth = int(inputArray[1])
    depths = {'100-499m': (100, 500), '500-899m': (500, 900),
              '900-1299m': (900, 1300), '1300-1699m': (1300, 1700),
              '1700-2099m': (1700, 2100), '2100-2499m': (2100, 2500),
              '2500-2899m': (2500, 2900), '2900-3299m': (2900, 3300),
              '3300-3699m':(3300, 3700), '3700-4099m':(3700, 4100)}
    for interval, depth_lst in depths.items():
        if depth_lst[0] <= depth < depth_lst[1]: # updated to allow for real and integer values
            ddir = current_directory+os.path.sep+str(interval)+os.path.sep
            generator = tf.keras.models.load_model(ddir+'co2_generator.h5', compile=False)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                input_transformer = joblib.load(ddir+'co2_input_transformer.pkl')
                output_transformer = joblib.load(ddir+'co2_output_transformer.pkl')

    # Transform the data
    input_trans = input_transformer.transform(input_images)

    # Convert to tensors
    BATCH_SIZE = 1
    test_dataset = tf.data.Dataset.from_tensor_slices(input_trans)
    test_dataset = test_dataset.batch(BATCH_SIZE)

    for example_input in test_dataset:
        prediction = generator(example_input, training=True)

    transformed_prediction = 10.0**(output_transformer.inverse_transform(
        prediction[0, :, :, :, :]))-1.0e-30
    # return input_images, transformed_prediction
    return transformed_prediction

if __name__ == '__main__':

    # parameters from test run 6 of 20 (2500-2899m)
    thick = 2.188277246944411161e+02
    depth = 2.594924843502177737e+03 # depth is positive in the IAM
    por = 2.141948497788563055e-01
    log_permh = -1.179646354060766456e+01
    log_aniso = 2.751080596962659719e+00
    aquifer_salinity = 4.706932213433123215e-03
    log_co2_rate = 8.269599291104867689e-02
    log_brine_rate = -3.190952070798314111e-01
    reservoir_salinity = 4.101526676315560282e-02

    time = 70 # years
    co2_mass_leaked = 10**log_co2_rate * time * 365.25 * 86400 # kg
    brine_mass_leaked = 10**log_brine_rate * time * 365.25 * 86400 # kg

    inputArray = np.array([thick, depth, por, log_permh, log_aniso, co2_mass_leaked,
                           brine_mass_leaked, aquifer_salinity, reservoir_salinity, time])

    inp = generate_input(inputArray)
    pred = model(inputArray)

    def generate_images(actual_input, actual_prediction):

        min_pred = actual_prediction.min()
        max_pred = actual_prediction.max()

        actual_time = actual_input[0, :, 0, :, 9].max()
        actual_co2_mass = np.expm1(actual_input[0, :, 0, :, 6].max())
        actual_water_mass = np.expm1(actual_input[0, :, 0, :, 7].max())
        actual_salt_mass = np.expm1(actual_input[0, :, 0, :, 8].max())


        input_titles = ['X, m', 'Z, m',
                        'X Permeability, log10 m$^2$', 'Z Permeability, log10 m$^2$',
                        'Porosity', 'Initial Aqueous Salt, mass fraction',
                        'CO$_2$ Mass, log1p kg', 'Water Mass, log1p kg',
                        'Salt Mass, log1p kg', 'Time, days']
        input_colors = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges',
                        'Reds','Greys', 'Purples', 'Blues', 'Greens']
        output_titles = ['Target Aqueous CO$_2$, mass fraction',
                         'Predicted Aqueous CO$_2$, mass fraction',
                         'Aqueous CO$_2$ Difference, mass fraction']
        output_colors = ['Blues', 'Blues', 'RdBu_r']

        crop=50

        if actual_time == 100.0: # plot fixed input

            # static input
            for var in [1, 2, 3, 4]:
                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(13.3, 2.5),
                                       sharex=True, sharey=True)
                col = 0
                ax.set_title(input_titles[var])
                perm = ax.contourf(actual_input[0, :crop, 0, :, 0],
                                   actual_input[0, :crop, 0, :, 1],
                                   actual_input[0, :crop, 0, :, var],
                                   cmap=input_colors[var])
                fig.colorbar(perm, ax=ax)
                plt.tight_layout()
                plt.savefig(pp, format='pdf')
                plt.close()


        else: # plot time-varying input and output

            for var in [5, 6, 7, 8]:
                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(13.3, 2.5),
                                       sharex=True, sharey=True)
                col = 0
                ax.set_title(input_titles[var])
                perm = ax.contourf(actual_input[0, :crop, 0, :, 0],
                                   actual_input[0, :crop, 0, :, 1],
                                   actual_input[0, :crop, 0, :, var],
                                   cmap=input_colors[var])
                fig.colorbar(perm, ax=ax)
                plt.tight_layout()
                plt.savefig(pp, format='pdf')
                plt.close()

            # output
            fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(13.3, 7.5),
                                    sharex=True, sharey=True)
            fig.suptitle(
                'Time: ' + str("{:.1f}".format(actual_time - 100)) +
                ' yr; CO$_2$ Mass: ' + str("{:.2E}".format(actual_co2_mass)) +
                ', kg; Water Mass: ' + str("{:.2E}".format(actual_water_mass)) +
                ', kg; Salt Mass: ' + str("{:.2E}".format(actual_salt_mass)) +
                ', kg', fontsize=16)

            # target
            # col = 0
            # sat = []
            # axs[col].set_title(output_titles[col])
            # sat.append(axs[col].contourf(actual_input[i, :crop, 0, :, 0],
            #                              actual_input[i, :crop, 0, :, 1],
            #                              actual_target[i, :crop, 0, :, 0],
            #                              norm=colors.Normalize(vmin=min_target,
            #                                                    vmax=max_target),
            #                              cmap=output_colors[col]))
            # fig.colorbar(sat[0], ax=axs[col])

            # prediction
            col = 1
            psat = []
            axs[col].set_title(output_titles[col])
            psat.append(axs[col].contourf(actual_input[0, :crop, 0, :, 0],
                                          actual_input[0, :crop, 0, :, 1],
                                          actual_prediction[0, :crop, 0, :, 0],
                                          norm=colors.Normalize(vmin=min_pred,
                                                                vmax=max_pred),
                                          cmap=output_colors[col]))
            fig.colorbar(psat[0], ax=axs[col])

            # # difference
            # col = 2
            # dsat = []
            # axs[col].set_title(output_titles[col])
            # try:
            #     dsat.append(axs[col].contourf(
            #         actual_input[i, :crop, 0, :, 0],
            #         actual_input[i, :crop, 0, :, 1],
            #         actual_difference[i, :crop, 0, :, 0],
            #         norm=colors.TwoSlopeNorm(vcenter=0, vmin=min_difference),
            #         cmap=output_colors[col]))
            # except:
            #     dsat.append(axs[col].contourf(
            #         actual_input[i, :crop, 0, :, 0],
            #         actual_input[i, :crop, 0, :, 1],
            #         actual_difference[i, :crop, 0, :, 0],
            #         norm=colors.TwoSlopeNorm(vcenter=0, vmax=max_difference),
            #         cmap=output_colors[col]))
            # fig.colorbar(dsat[0], ax=axs[col])

            plt.tight_layout()
            plt.savefig(pp, format='pdf')
            plt.close()

    pp = PdfPages('co2_predict.pdf')
    generate_images(inp, pred)
    pp.close()

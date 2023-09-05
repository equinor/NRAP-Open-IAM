import math
import numpy as np
import pandas as pd


# radial grid used in STOMP simulations
def get_vertices(depth, thick, n_vert):
    # grid vertices
    x_vert = np.array([
        0.00, 3.24, 6.73, 10.49, 14.56, 18.94, 23.67, 28.77, 34.27,
        40.20, 46.60, 53.51, 60.96, 69.00, 77.67, 87.02, 97.11,
        107.99, 119.73, 132.40, 146.06, 160.80, 176.70, 193.85, 212.35,
        232.31, 253.84, 277.07, 302.12, 329.15, 358.31, 389.76, 423.69,
        460.29, 499.77, 542.37, 588.31, 637.88, 691.35, 749.03, 811.25,
        878.37, 950.77, 1028.88, 1113.14, 1204.04, 1302.09, 1407.87,
        1521.97, 1645.06, 1777.84, 1921.08, 2075.60, 2242.29, 2422.11,
        2616.08, 2825.33, 3051.06, 3294.57, 3557.25, 3840.61, 4146.29,
        4476.05, 4831.77, 5215.50, 5629.45, 6076.00, 6557.72, 7077.37,
        7637.94, 8242.66, 8895.00, 9598.71, 10357.83, 11176.74, 12060.13,
        13013.09, 14041.10, 15150.06, 16346.35, 17636.85, 19028.97,
        20530.72, 22150.73, 23898.32, 25783.53, 27817.19, 30011.01,
        32377.58, 34930.52, 37684.51, 40655.36, 43860.17, 47317.36,
        51046.79, 55069.92, 59409.86, 64091.57, 69141.96, 74590.06,
        80467.20])

    z_vert = []
    for i in range(n_vert+1):
        z_vert[:0] = [-depth - thick/n_vert*i]
    z_vert = np.array(z_vert)

    return x_vert, z_vert

# convert vertices to centroids
def vert_to_cent(vert):
    return pd.Series(vert).rolling(2).mean().to_numpy()[1:]

def get_dimensions(x_vert, z_vert):
    dz = pd.Series(z_vert).diff().to_numpy()[1:]
    dx = pd.Series(x_vert).diff().to_numpy()[1:]
    dx_mat, dz_mat = np.meshgrid(dx, dz, indexing='ij')

    r2 = np.square(x_vert)
    area = pd.Series(r2).diff().to_numpy()[1:] * math.pi
    area_mat, dz_mat = np.meshgrid(area, dz, indexing='ij')
    vol_mat = area_mat * dz_mat

    return dx_mat, dz_mat, vol_mat

def get_mesh(x_vert, z_vert):
    # grid centroids
    x_cent = vert_to_cent(x_vert)
    z_cent = vert_to_cent(z_vert)

    # create coordinate matrices
    _, x_mat, _, z_mat = np.meshgrid(
        np.array([0]), x_cent, np.array([0]), z_cent, indexing='ij')

    return x_mat, z_mat

def generate_input(inputArray):
    thick = inputArray[0]
    depth = inputArray[1]
    por = inputArray[2]
    log_permh = inputArray[3]
    log_aniso = inputArray[4]
    co2_mass_leaked = inputArray[5]
    brine_mass_leaked = inputArray[6]
    init_aq_salt_mass_frac = inputArray[7]
    leak_aq_salt_mass_frac = inputArray[8]
    time = inputArray[9]

    n_vert = 10 # number of grid cells in the vertical direction
    x_vert, z_vert = get_vertices(depth, thick, n_vert)
    x_mat, z_mat = get_mesh(x_vert, z_vert)

    por_mat = np.full_like(x_mat, por)
    xperm_mat = np.full_like(x_mat, log_permh)
    zperm_mat = np.full_like(x_mat, log_permh-log_aniso)

    well_mat = np.zeros_like(x_mat) # (1, 100, 1, 10)
    well_mat[0, 0, 0, 1:n_vert-1] = 1.0
    co2_mass_mat = well_mat * np.log1p(co2_mass_leaked /(n_vert-2))
    water_mass_leaked = brine_mass_leaked * (1. - leak_aq_salt_mass_frac * 0.93)
    water_mass_mat = well_mat * np.log1p(water_mass_leaked/(n_vert-2))
    salt_mass_leaked = brine_mass_leaked * leak_aq_salt_mass_frac * 0.93
    salt_mass_mat = well_mat * np.log1p(salt_mass_leaked/(n_vert-2))

    time_mat = np.full_like(x_mat, time+100.) # account for 100-year equilibration period

    init_salt_mat = np.ones_like(x_mat) * np.log10(init_aq_salt_mass_frac + 1.e-10)

    input_images = np.stack((
        x_mat, z_mat, xperm_mat, zperm_mat, por_mat, init_salt_mat, co2_mass_mat,
        water_mass_mat, salt_mass_mat, time_mat), axis=4)

    return input_images

if __name__ == "__main__":

    # parameters from test run 1006 of 625 (2100-2499m)
    thick = 1.3822E+02
    depth = 2.2416E+03

    xv, zv = get_vertices(depth, thick, 10)
    dx, dz, vol = get_dimensions(xv, zv)
    print('vol', vol)
    print('dx', dx)
    print('dz', dz)
    print('vol', vol)

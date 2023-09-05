# -*- coding: utf-8 -*-
"""
Created on Sat Mar 03 06:58:10 2018

@author: Seth
"""
import matplotlib.pyplot as plt
import numpy as np


def reservoir_contours(filename, figsiz=(15, 15)):
    raw_data = np.genfromtxt(filename, delimiter=',')
    points = raw_data[:, :2]
    data = raw_data[:, 2:]

    n_time = int(data.shape[1]/2)
    press_fig, p_ax = plt.subplots(4, 4, num='Reservoir Pressure',
                                   figsize=figsiz)
    p_ax = p_ax.flatten()
    sat_fig, sat_ax = plt.subplots(4, 4, num='Reservoir Saturation',
                                   figsize=figsiz)
    sat_ax = sat_ax.flatten()
    fig_num = 1
    i_ax = 0
    for it in range(n_time):
        p_ax[i_ax].tricontourf(points[:, 0], points[:, 1], data[:, it])
        sat_ax[i_ax].tricontourf(points[:, 0],
                                 points[:, 1],
                                 data[:, it+n_time])
        i_ax += 1
        if i_ax >= 16:
            fig_num += 1
            press_fig, p_ax = plt.subplots(4, 4,
                                           num='Reservoir Pressure {}'
                                           .format(fig_num),
                                           figsize=figsiz)
            p_ax = p_ax.flatten()
            sat_fig, sat_ax = plt.subplots(4, 4,
                                           num='Reservoir Saturation {}'
                                           .format(fig_num),
                                           figsize=figsiz)
            sat_ax = sat_ax.flatten()
            i_ax = 0
    return

if __name__ == '__main__':
    filename = 'Reservoir_data_sim01.csv'
    figsiz = (15, 15)
    reservoir_contours(filename, figsiz=figsiz)

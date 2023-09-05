from bisect import bisect

import numpy as np
import matplotlib.pyplot as plt


class Well():
    def __init__(self, x, y, name=''):
        """ Create well with the provided coordinates."""
        self.x = x
        self.y = y
        self.name = name

class Zone():
    def __init__(self, min, max):
        self.min = min
        self.max = max
        self.ticks = []

    def add_ticks(self, min_spacing=1, mul=0):
        if mul > 0:
            self.ticks = [self.min]
            spacing = min_spacing
            while self.ticks[-1] < self.max - min_spacing:
                self.ticks.append(self.ticks[-1] + spacing)
                spacing = spacing * mul
            self.ticks[-1] = self.max
        elif mul < 0:
            self.ticks = [self.max]
            spacing = min_spacing
            while self.ticks[-1] > self.min + min_spacing:
                self.ticks.append(self.ticks[-1] - spacing)
                spacing = spacing * -mul
            self.ticks[-1] = self.min
        elif mul == 0:
            self.ticks = list(np.arange(self.min, self.max, min_spacing))
            self.ticks.append(self.max)

    def build_ticks(self, alist):
        self.ticks.extend(alist)

class Axis():
    def __init__(self, name):
        self.name = name
        self.zones = []

    def add_zone(self, min, max):
        self.zones.append(Zone(min, max))
        return self.zones[-1]

    def get_ticks(self):
        ticks = []
        for zone in self.zones:
            ticks.extend(zone.ticks)
        # remove duplicates
        ticks = list(set(ticks))
        ticks.sort()
        return ticks

class Grid():
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """ Create zone within the provided boundaries."""
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.axis = {}
        self.wells = []

        # set up axes with one big node as default
        xaxis = self.add_axis(name='x')
        yaxis = self.add_axis(name='y')
        zaxis = self.add_axis(name='z')
        xzone = xaxis.add_zone(min=xmin, max=xmax)
        yzone = yaxis.add_zone(min=ymin, max=ymax)
        zzone = zaxis.add_zone(min=zmin, max=zmax)
        xzone.add_ticks(min_spacing=(self.xmax-self.xmin), mul=0)
        yzone.add_ticks(min_spacing=(self.ymax-self.ymin), mul=0)
        zzone.add_ticks(min_spacing=(self.zmax-self.zmin), mul=0)

    def add_axis(self, name):
        self.axis[name] = Axis(name)
        return self.axis[name]

    def build_axis(self, name, alist):
        self.axis[name] = Axis(name)
        zone = self.axis[name].add_zone(min=min(alist), max=max(alist))
        zone.build_ticks([alist])
        return self.axis[name]

    def get_axis(self, name):
        return self.axis[name]

    def add_well_object(self, well):
        self.wells.append(well)

    def add_well(self, x, y, name):
        self.wells.append(Well(x=x, y=y, name=name))

    def get_well_coords(self):
        return [(well.x, well.y) for well in self.wells]

    def get_well_x_coords(self):
        xcoords = [well.x for well in self.wells]
        xcoords.sort()
        return xcoords

    def get_well_y_coords(self):
        ycoords = [well.y for well in self.wells]
        ycoords.sort()
        return ycoords

    def get_well_grid_locations(self):
        coords = self.get_well_coords()
        xticks = self.get_axis('x').get_ticks()
        yticks = self.get_axis('y').get_ticks()
        locs = [(bisect(xticks, coord[0]), bisect(yticks, coord[1])) for coord in coords]
        return locs

    @staticmethod
    def add_mid_end_points(coords, mincoord, maxcoord):
        newcoords = []
        num_coords = len(coords)
        for w in range(num_coords-1):
            newcoords.append((coords[w]+coords[w+1])/2)
        newcoords.append(mincoord)
        newcoords.append(maxcoord)
        newcoords.sort()
        coords.extend(newcoords)
        coords.sort()
        return coords

    def refine_tics_around_wells(self, bounds, axis_name, min_spacing, mul):
        num_bounds = len(bounds)
        for i in range(num_bounds-1):
            # decreasing spacing
            if i%2 == 0:
                # print(axis_name, bounds[i], bounds[i+1], min_spacing, -mul)
                zonex = self.axis[axis_name].add_zone(min=bounds[i], max=bounds[i+1])
                zonex.add_ticks(min_spacing=min_spacing, mul=-mul)
            # increasing spacing
            else:
                # print(axis_name, bounds[i], bounds[i+1], min_spacing, mul)
                zonex = self.axis[axis_name].add_zone(min=bounds[i], max=bounds[i+1])
                zonex.add_ticks(min_spacing=min_spacing, mul=mul)

    def refine_regular_grid(self, num_x, num_y, num_z):
        xaxis = self.get_axis(name='x')
        yaxis = self.get_axis(name='y')
        zaxis = self.get_axis(name='z')
        xzone = xaxis.add_zone(min=self.xmin, max=self.xmax)
        yzone = yaxis.add_zone(min=self.ymin, max=self.ymax)
        zzone = zaxis.add_zone(min=self.zmin, max=self.zmax)
        xzone.add_ticks(min_spacing=(self.xmax-self.xmin)/num_x, mul=0)
        yzone.add_ticks(min_spacing=(self.ymax-self.ymin)/num_y, mul=0)
        zzone.add_ticks(min_spacing=(self.zmax-self.zmin)/num_z, mul=0)

    def refine_grid_around_wells(self, min_x_spacing, min_y_spacing, x_mul, y_mul):
        xcoords = self.get_well_x_coords()
        ycoords = self.get_well_y_coords()
        xbounds = self.add_mid_end_points(xcoords, self.xmin, self.xmax)
        ybounds = self.add_mid_end_points(ycoords, self.ymin, self.ymax)
        self.refine_tics_around_wells(xbounds, 'x', min_x_spacing, x_mul)
        self.refine_tics_around_wells(ybounds, 'y', min_y_spacing, y_mul)

    def get_vertices(self):
        xticks = self.get_axis('x').get_ticks()
        yticks = self.get_axis('y').get_ticks()
        zticks = self.get_axis('z').get_ticks()
        return np.meshgrid(xticks, yticks, zticks, indexing='ij')

    def get_centroids(self):
        xticks = np.array(self.get_axis('x').get_ticks())
        yticks = np.array(self.get_axis('y').get_ticks())
        zticks = np.array(self.get_axis('z').get_ticks())
        xcent = (xticks[1:] + xticks[:-1]) / 2
        ycent = (yticks[1:] + yticks[:-1]) / 2
        zcent = (zticks[1:] + zticks[:-1]) / 2
        return np.meshgrid(xcent, ycent, zcent, indexing='ij')

    # radial distance from well location
    @staticmethod
    def get_radial_distance(x, y, xx, yy):
        rr = np.sqrt((xx-x)**2 + (yy-y)**2)
        return rr

    # plot grid in x-y plane
    def plot_xy(self):
        xaxis = self.get_axis(name='x')
        yaxis = self.get_axis(name='y')
        xloc = xaxis.get_ticks()
        yloc = yaxis.get_ticks()
        _, ax = plt.subplots()
        ax.set_title('nx='+str(len(xloc))+' ny='+str(len(yloc)))
        ax.set_xticks(xloc)
        ax.set_yticks(yloc)
        ax.grid(zorder=0)
        ax.set_xlim(xloc[0], xloc[-1])
        ax.set_ylim(yloc[0], yloc[-1])
        for i, well in enumerate(self.wells):
            ax.scatter(well.x, well.y, zorder=i+2)
            ax.annotate(well.name, (well.x, well.y),
                        xycoords='data', textcoords='data', zorder=1)
        plt.show()

    # print grid in STOMP format
    def print_stomp(self):
        xaxis = self.get_axis(name='x')
        yaxis = self.get_axis(name='y')
        zaxis = self.get_axis(name='z')
        xloc = xaxis.get_ticks()
        yloc = yaxis.get_ticks()
        zloc = zaxis.get_ticks()
        print('Cartesian,')
        print(len(xloc)-1, ',', len(yloc)-1, ',', len(zloc)-1, ',', end='', sep='')

        for i, x in enumerate(xloc):
            if i % 6 == 0:
                print()
            print(x, ',m,', end='', sep='')

        for i, y in enumerate(yloc):
            if i % 6 == 0:
                print()
            print(y, ',m,', end='', sep='')

        for i, z in enumerate(zloc):
            if i % 6 == 0:
                print()
            print(z, ',m,', end='', sep='')
        print()


if __name__ == "__main__":

    # set grid boundaries
    grid = Grid(xmin=0, xmax=1000, ymin=0, ymax=1000, zmin=0, zmax=900)

    # add wells
    grid.add_well(x=739, y=423, name='Well 1')
    grid.add_well(x=584, y=822, name='Well 2')
    grid.add_well(x=287, y=333, name='Well 3')
    grid.add_well(x=519, y=856, name='Well 4')

    grid.refine_grid_around_wells(min_x_spacing=5, min_y_spacing=5, x_mul=1.5, y_mul=1.5)

    grid.plot_xy()

    # regular grid example
    grid2 = Grid(xmin=0, xmax=1000, ymin=0, ymax=1000, zmin=0, zmax=900)
    grid2.refine_regular_grid(num_x=10, num_y=10, num_z=10)
    grid2.print_stomp()

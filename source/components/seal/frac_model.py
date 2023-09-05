#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module contains classes for fracture computations.

Author: Ernest N. Lindner
Date: 05/24/2022

Module Name
    frac_model

Class Functions (26)
    ** Class Pointe **** (6)
        __init__(self, x_value, y_value)
        __repr__(self)  - NOT USED
        replica(self)
        move(self, del_x, del_y)  - NOT USED
        distant(self, end_pt)
        get_coords(self)
    ** Class Segment **** (9)
        __init__(self, start_pt, end_pt)
        replica(self)  - NOT USED
        __repr__(self) - NOT USED
        length(self)
        trend(self)
        params_2(self) - NOT USED
        delist(self)
        get_mid_point(self)
        get_points(self)
    ** Class Rectangle **** (6)
        __init__(self, bottom_left_pt, width, length)
        replica(self)
        area(self)
        side(self, side_number) - NOT USED
        xlimits(self)
        ylimits(self)
    ** Class Fracture **** (5)
        __init__(self, point_1, point_2)
        orient(self) - NOT USED
        extent(self):
        print_fracture(self)
        frac_array(self)

Related Functions (7)
    round_half_up(number, decimals=0)
    create_another_line(stage, original, t_right, left_clip, delta)
    clipping_2(box, clip)
    line_from_center(line_length, theta, center_x, center_y)
    find_point_coordinates(point_1, point_2, extend_a) - NOT USED
    find_determinant(point_1, point_2, point_3, point_4) - NOT USED
    find_intersection(frac_1, frac_2) - NOT USED

Copyright(c) 2019-2022 by Ernest N. Lindner - All Rights Reserved
-------------------------------------------------------------------------------
"""
import math                         # For trigonometry, infinite, hypotenuse
import copy                         # For <deepcopy>

# Constants
MINI_LINE = 0.1                     # Minimum line length (m) > micro-fractures
DISPLAY = 10                        # Round coordinates to decimal
V_SMALL = 1.0E-08                   # Very small number in division by zero
V_SMALL2 = 1.0E-04                  # Very small number in atan


#   ---------------------------------------------------------------------------
class Pointe:
    """A class of geometric points in 2D."""

    def __init__(self, x_value=0.0, y_value=0.0):
        """Create a point in 2D cartesian coordinates.

        Parameters
        ----------
        x_value = (float) coordinate along E-W
        axis (m)
        y_value = (float) coordinate along N-S axis (m)

        Returns
        -------
        None
        """
        self.x_coord = x_value
        self.y_coord = y_value

    def __repr__(self):
        """Return a command string for printing.

        Parameters
        ----------
        N/A

        Returns
        -------
        x_coord, y_coord = (float) x- and y-coordinates of a point
        """
        coordinates = f'({self.x_coord}, {self.y_coord})'
        return coordinates

    def replica(self):
        """Return a copy of a point.

        Parameters
        ----------
        N/A

        Returns
        -------
        None
        """
        new_point = copy.copy(self)
        return new_point

    def move(self, del_x=0.0, del_y=0.0):
        """Move point an increment.

        Parameters
        ----------
        del_x = (float) distance along E-W axis (m)
        del_y = (float) distance along N-S axis (m)

        Notes
        -----
        1. Assumes rectangular coordinate system, N-S (Y-axis), E-W (x-axis)
        """
        self.x_coord += del_x
        self.y_coord += del_y
        # return None

    def distant(self, end_pt):
        """Distance of another point away from current point.

        Parameters
        ----------
        end_pt = (float) another point

        Returns
        -------
        distance = float) distance between points

        Notes
        -----
        1. Python: hypot = sqrt(x*x + y*y).
        """
        delta_x = end_pt.x_coord - self.x_coord
        delta_y = end_pt.y_coord - self.y_coord
        return math.hypot(delta_x, delta_y)

    def get_coords(self):
        """Get X and Y coordinates of a point.

        Parameters
        ----------
        N/A

        Returns
        -------
        x_coord = (float) coordinate along N-S axis (m)
        y_coord = (float) coordinate along E-W axis (m)
        """
        return self.x_coord, self.y_coord


#   ---------------------------------------------------------------------------
class Segment:
    """A class of line segments in 2D."""

    def __init__(self, start_pt, end_pt):
        """Create a line in 2D cartesian coordinates with two points.

        Parameters
        ----------
        start_pt = (pointe class) start point (left)
        end_pt = (pointe class) end point (right)

        Returns
        -------
        None

        Notes
        -----
        1. Convention is that left point is smaller value along X-axis.
        """
        if start_pt.x_coord > end_pt.x_coord:
            self.pt_2 = start_pt
            self.pt_1 = end_pt
        else:
            self.pt_1 = start_pt
            self.pt_2 = end_pt
        # return None

    def replica(self):
        """Return a copy of a line.

        Parameters
        ----------
        N/A

        Returns
        -------
        new_point= (segment class) copy of line
        """
        new_line = copy.deepcopy(self)
        return new_line

    def __repr__(self):
        """Define str for representation of line segment - points.

        Parameters
        ----------
        N/A

        Returns
        -------
        string of line points
        """
        # define format: (pt1.x, pt1.y)--(pt2.x, pt2.y).
        liner = str(f'({self.pt_1.x_coord}, {self.pt_1.y_coord})--'
                    f'({self.pt_2.x_coord}, {self.pt_2.y_coord})')
        return liner

    def length(self):
        """Compute the length of line.

        Parameters
        ----------
        N/A

        Returns
        -------
        distance = distance between two points

        Notes
        -----
        1. Python: hypot = sqrt(x*x + y*y).
        """
        delta_x = self.pt_2.x_coord - self.pt_1.x_coord
        delta_y = self.pt_2.y_coord - self.pt_1.y_coord
        distance = math.hypot(delta_x, delta_y)
        return distance

    def trend(self):
        """Compute strike of line from North.

        Parameters
        ----------
        N/A

        Returns
        -------
        angle = (float) angle from north (deg),
            -> angle must be within < +/-90 deg
        """
        delta_x = self.pt_2.x_coord - self.pt_1.x_coord
        delta_y = self.pt_2.y_coord - self.pt_1.y_coord

        # keep angle positive, above 90 degrees.
        if delta_y < 0.0:
            delta_y = -delta_y
            delta_x = -delta_x

        # prevent crash if delta_y = 0.
        if delta_y < V_SMALL2:
            north_angle = 90.0
        else:
            north_angle = math.degrees(math.atan(delta_x / delta_y))
        return north_angle

    def params_2(self):
        """Define slope and intercept (m, b) for 2D line.

        Parameters
        ----------
        N/A

        Returns
        -------
        incline = (float) slope of line
        beta = (float) line constant

        Notes
        -----
        1. Equation: <y = mx + b>; where m = slope and b = beta.
        """
        delta_x = self.pt_2.x_coord - self.pt_1.x_coord
        delta_y = self.pt_2.y_coord - self.pt_1.y_coord

        if abs(delta_x) < V_SMALL:
            incline = math.inf
            beta = self.pt_1.y_coord
        else:
            incline = delta_y / delta_x
            if abs(incline) < V_SMALL:
                incline = 0.0
            beta = self.pt_1.y_coord - incline * self.pt_1.x_coord
        return incline, beta

    def delist(self):
        """Unlink line parameters for 2D list creation.

        Parameters
        ----------
        N/A

        Returns
        -------
        pt1.x, pt1.y, pt2.x, pt2.y = (float) 2D list of coordinates
        """
        return([(self.pt_1.x_coord, self.pt_1.y_coord),
                (self.pt_2.x_coord, self.pt_2.y_coord)])

    def get_mid_point(self):
        """Get x-, y-coordinates of line's center point.

        Parameters
        ----------
        N/A

        Returns
        -------
        x_new, y_new = (floats) Left corner coordinates for matplotlib
        """
        delta_x = (self.pt_2.x_coord - self.pt_1.x_coord) / 2.0
        delta_y = (self.pt_2.y_coord - self.pt_1.y_coord) / 2.0
        x_new = self.pt_1.x_coord + delta_x
        y_new = self.pt_1.y_coord + delta_y
        return x_new, y_new

    def get_points(self):
        """Get points of 2D segment.

        Parameters
        ----------
        N/A

        Returns
        -------
        pt_1, pt_2 = (pointe class) points defining line
        """
        return self.pt_1, self.pt_2


#   ---------------------------------------------------------------------------
class Rectangle:
    """A class of rectangles in 2D."""

    def __init__(self, bottom_left_pt, width, length):
        """Create a box in 2D cartesian coordinates with two points.

        Parameters
        ----------
        bottom_left_pt = (float) bottom left corner point
        width = (float) x distance (E-W)

        Returns
        -------
        None

        Notes
        -----
        1. Convention is that left point is smaller value along X-axis.
        """
        x_new = bottom_left_pt.x_coord + width
        y_new = bottom_left_pt.y_coord + length

        self.pt_1 = bottom_left_pt.replica()
        self.pt_2 = Pointe(x_new, bottom_left_pt.y_coord)
        self.pt_3 = Pointe(x_new, y_new)
        self.pt_4 = Pointe(bottom_left_pt.x_coord, y_new)

        # return None

    def replica(self):
        """Return a copy of a rectangle.

        Parameters
        ----------
        N/A

        Returns
        -------
        new_box = (rectangle class) new rectangle
        """
        new_box = copy.deepcopy(self)
        return new_box

    def area(self):
        """Compute area of rectangle.

        Parameters
        ----------
        N/A

        Returns
        -------
        box_area = (float) area of rectangle

        Notes
        -----
        Assumes points are numbered counter-clockwise from bottom left,
        starting at bottom left corner => p1.
        """
        x_side = self.pt_2.x_coord - self.pt_1.x_coord
        y_side = self.pt_4.y_coord - self.pt_1.y_coord
        return abs(x_side * y_side)

    def side(self, side_number):
        """Return side of rectangle as line.

        Parameters
        ----------
        side_number = (int) 1,2,3,4 side

        Returns
        -------
        side_line = (segment class) line of selected side

        Notes
        -----
        Assumes sides are numbered counter-clockwise from bottom.
        """
        if side_number == 1:
            side_line = Segment(self.pt_1, self.pt_2)
        elif side_number == 2:
            side_line = Segment(self.pt_2, self.pt_3)
        elif side_number == 3:
            side_line = Segment(self.pt_3, self.pt_4)
        else:
            side_line = Segment(self.pt_4, self.pt_1)
        return side_line

    def xlimits(self):
        """Return x limits of box.

        Parameters
        ----------
        N/A

        Returns
        -------
        x_min, x_max = (floats) minimum/maximum x coordinates

        Notes
        -----
        1. Assumes points are numbered counter-clockwise,
            starting at bottom left corner => pt_1.
        """
        x_min = self.pt_1.x_coord
        x_max = self.pt_2.x_coord
        return x_min, x_max

    def ylimits(self):
        """Return y limits of box.

        Parameters
        ----------
        N/A

        Returns
        -------
        y_min, y_max = (floats) minimum/maximum y coordinates of box

        Notes
        -----
        1. Assumes points are numbered clockwise,
            starting at bottom left corner => pt_1.
        """
        y_min = self.pt_1.y_coord
        y_max = self.pt_4.y_coord
        return y_min, y_max


#   ---------------------------------------------------------------------------
class Fracture:
    """A class of fractures for equivalent permeability."""

    def __init__(self, point_1, point_2):
        """Initialize fracture - in 2D.

        Parameters
        ----------
        profile => Segment
            pt_1 = Pointe - Left
            pt_2 = Pointe - Right
        aperture = fracture aperture (mm)
        strike = fracture strike (deg)
        trans = fracture transmissivity (m^4)
        threshold = threshold pressure (Pa)

        Returns
        -------
        None

        Notes
        -----
        1. fracture is assumed 2-D and smooth.
        2. Assumes the following:
            dip = 0.0           # vertical
            radius = 0.0        # 2D; fracture is line
            filled = None       # open, i.e. no filling
            roughness = 0.0     # smooth
            tortuosity = 0.0    # straight.
        """
        self.profile = Segment(point_1, point_2)
        self.aperture = 0.0
        self.strike = 0.0
        self.trans = 0.0
        self.threshold = 0.0
        # return None

    def orient(self):
        """Get strike of a fracture.

        Parameters
        ----------
        N/A

        Returns
        -------
        result = (float) strike of line
        """
        return self.profile.trend()

    def extent(self):
        """Get length of a fracture.

        Parameters
        ----------
        N/A

        Returns
        -------
        result = (float) length of line
        """
        result = self.profile.length()
        return result

    def print_fracture(self):
        """Print instance properties of a fracture.

        Parameters
        ----------
        N/A

        Returns
        -------
        None
        """
        print("\n Fracture Parameters for instance are:")
        print(f'   Start_point    = {self.profile.pt_1.x_coord}')
        print(f'   End_point      = {self.profile.pt_2.x_coord}')
        print(f'   Aperture       = {self.aperture} mm')
        print(f'   Strike         = {self.strike} deg')
        print(f'   Transmissivity = {self.trans: 0:7.4e} m4')
        print(f'   Threshold      = {self.threshold} Pa')
        # return None

    def frac_array(self):
        """Provide fracture parameters for a list.

        Parameters
        ----------
        N/A

        Returns
        -------
        result = (list) list of instance fracture data
        """
        result = [self.profile.pt_1.x_coord, self.profile.pt_1.y_coord,
                  self.profile.pt_2.x_coord, self.profile.pt_2.y_coord,
                  self.aperture, self.strike, self.trans, self.threshold]

        return result

#   ---------------------------------------------------------------------------
#   RELATED FUNCTIONS


def round_half_up(valu, decimals=0):
    """Round-up a floating point to a defined number of digits.

    Parameters
    ----------
    valu = (float) number to be rounded
    decimals = (int) number of digits per the decimal

    Returns
    -------
    results = (float) rounded number
    """
    multiplier = math.pow(10.0, decimals)
    results = math.floor(valu * multiplier + 0.5) / multiplier

    return results


# noinspection PyTypeChecker,PyTypeChecker
def create_another_line(stage, original, clip_left, clip_right, delta):
    """Define a new 2D line after clipping original line.

    Parameters
    ----------
    stage = (bool) status of clipping effort (True/False)
    original = (line class) original line
    left_clip = (float) t value to left
    right_clip = (float) t value to right
    delta = list with delta x and delta y

    Returns
    -------
    another_line = clipped line
    """
    # Define new line if desired.
    if stage:

        # Correct starting point.
        if clip_left == 0.0:
            # Keep original values if no clipping.
            start_x = original.pt_1.x_coord
            start_y = original.pt_1.y_coord
        else:
            # Clip values.
            start_x = original.pt_1.x_coord + (clip_left * delta[0])
            start_y = original.pt_1.y_coord + (clip_left * delta[1])

        # Correct end point.
        if clip_right == 1.0:
            # keep original values if no clipping.
            end_x = original.pt_2.x_coord
            end_y = original.pt_2.y_coord
        else:
            # Clip values from left.
            end_x = original.pt_1.x_coord + (clip_right * delta[0])
            end_y = original.pt_1.y_coord + (clip_right * delta[1])

        # Round values for clarity - eliminate negative values.
        start_x = round_half_up(start_x, DISPLAY)
        start_y = round_half_up(start_y, DISPLAY)
        end_x = round_half_up(end_x, DISPLAY)
        end_y = round_half_up(end_y, DISPLAY)

        # Create new line.
        start_point = Pointe(start_x, start_y)
        end_point = Pointe(end_x, end_y)
        new_line = Segment(start_point, end_point)
    else:
        new_line = original

    # Check for small lines or points - eliminate.
    if new_line.length() < MINI_LINE:
        stage = False

    return stage, new_line


def clipping_2(box, original):
    """Clip a 2D line to be within a 2D rectangular box.

    Parameters
    ----------
    box = (rectangle class) rectangle with 4 points
    original = (line class) 2D line to be clipped to be inside box

    Returns
    -------
    new_line = (line class) clipped line
    status = (bool) condition of line
        = False => existing line outside of box = trash
        = True => line is inside box (and clipped)

    Notes
    -----
    1) Liang-Barsky algorithm for clipping.
    2) Assumes pt_x2 > pt_x1 of line.
    3) assumes pt_1 of box is lower left, with counter-clockwise
        numbering.
    4) For theory, see URLs:
    https://www.skytopia.com/project/articles/compsci/clipping.html
    https://www.geeksforgeeks.org/liang-barsky-algorithm/
    """
    # Define/translate line variables.
    x_low = original.pt_1.x_coord
    y_low = original.pt_1.y_coord
    delta = [(original.pt_2.x_coord - original.pt_1.x_coord),
             (original.pt_2.y_coord - original.pt_1.y_coord)]

    # Set window/box coordinates.
    win_x = [box.pt_1.x_coord, box.pt_2.x_coord]
    win_y = [box.pt_2.y_coord, box.pt_3.y_coord]

    # Define/translate variables for clipping into lists.
    p_list = [-(delta[0]), delta[0], -(delta[1]), delta[1]]
    q_list = [(x_low - win_x[0]), (win_x[1] - x_low), (y_low - win_y[0]),
              (win_y[1] - y_low)]

    # Clipping code
    #   p[k] = 0                parallel to the clipping boundaries
    #   p[k] = 0 and q[k] < 0   completely outside the boundary
    #   p[k] = 0 and q[k] > 0   inside the parallel clipping boundary
    #   p[k] < 0                line proceeds from outside to inside
    #   p[k] > 0                line proceeds from inside to outside

    # Initialize controls
    tex_0 = 0.0
    tex_1 = 1.0
    status = True
    indx = 0

    # Loop until an error found in clipping or 4 sides.
    while status and indx < 4:

        # Check if line section less than or greater than.
        if p_list[indx] < 0.0:
            r_term = q_list[indx] / p_list[indx]
            if r_term > tex_1:
                # No clip.
                status = False
                break
            if r_term > tex_0:
                # Clip line on left.
                tex_0 = r_term

        elif p_list[indx] > 0.0:
            r_term = q_list[indx] / p_list[indx]
            if r_term < tex_0:
                # No clip
                status = False
                break
            if r_term < tex_1:
                # Clip line on left.
                tex_1 = r_term

        # Check if line section when parallel to a side and outside of box.
        #   * p_list[] == 0
        else:
            if q_list[indx] < 0.0:
                # No intersection.
                status = False
                break

        # increase counter for side.
        indx += 1
    # end loop

    # Define new line.
    status, new_line = create_another_line(status, original,
                                           tex_0, tex_1, delta)

    return status, new_line


# noinspection PyTypeChecker,PyTypeChecker
def line_from_center(line_length, theta, center_x, center_y):
    """Create new line from a center point.

    Parameters
    ----------
    line_length = (float) length of line
    theta = (float) angle from north (radians)
    center_x = (float) X-coordinate of center point
    center_y = (float) Y-coordinate of center point

    Returns
    -------
    new_line = segment instance
    """
    # Define 1/2 delta components.
    new_width = line_length * math.sin(theta) / 2.0
    new_length = line_length * math.cos(theta) / 2.0

    # Define left and right point coordinates of line - from center.
    right_x = center_x + new_width
    right_y = center_y + new_length
    left_x = center_x - new_width
    left_y = center_y - new_length

    # Define new line.
    right_point = Pointe(right_x, right_y)
    left_point = Pointe(left_x, left_y)
    new_line = Segment(left_point, right_point)

    return new_line


def find_point_coordinates(point_1, point_2, extend_a):
    """Find the intersection point using extension format.

    Parameters
    ----------
    point_1 = (pointe class) Pointe #1 - Line #1
    point_2 = (pointe class) Pointe #2 - Line #1
    extend_a = (float) factor for line #1

    Returns
    -------
    intersect = intersection point
    """
    # Define variables.
    pt_x1, pt_y1 = point_1.get_coords()
    pt_x2, pt_y2 = point_2.get_coords()

    # Compute intersection point coordinates.
    x_new = pt_x1 + extend_a * (pt_x2 - pt_x1)
    y_new = pt_y1 + extend_a * (pt_y2 - pt_y1)

    intersect = Pointe(x_new, y_new)

    return intersect


def find_determinant(point_1, point_2, point_3, point_4):
    """Find the determinant and extensions for intersection of segments.

    Parameters
    ----------
    point_1 = (pointe class) Pointe #1 (see seal_model for Class)
    point_2 = (pointe class) Pointe #2
    point_3 = (pointe class) Pointe #3
    point_4 = (pointe class) Pointe #4

    Returns
    -------
    denominator = (float) determinant
    extension_a = (float) numerator of extension (ta)
    extension_b = (float) numerator of extension (tb)

    Notes
    -----
    For theory, see URL:
        1. https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    """
    # Define variables.
    pt_x1, pt_y1 = point_1.get_coords()
    pt_x2, pt_y2 = point_2.get_coords()
    pt_x3, pt_y3 = point_3.get_coords()
    pt_x4, pt_y4 = point_4.get_coords()

    # Compute determinant.
    denominator = ((pt_x4 - pt_x3) * (pt_y1 - pt_y2)
                   - (pt_x1 - pt_x2) * (pt_y4 - pt_y3))

    # Compute extension numerators.
    extension_a = ((pt_y3 - pt_y4) * (pt_x1 - pt_x3) +
                   (pt_x4 - pt_x3) * (pt_y1 - pt_y3))
    extension_b = ((pt_y1 - pt_y2) * (pt_x1 - pt_x3) +
                   (pt_x2 - pt_x1) * (pt_y1 - pt_y3))

    return denominator, extension_a, extension_b


def find_intersection(frac_1, frac_2):
    """Determine if two fracture segments intersect.

    Parameters
    ----------
    frac_1 = (line class) fracture line #1
    frac_2 = (line class) fracture line #2

    Returns
    -------
    intersect = (pointe class) intersection point
    success = (bool) if intersection exists

    Notes
    -----
    1. Error code -999.9 for coordinates.
    """
    # Define basic variables.
    success = False  # default
    point_1, point_2 = frac_1.profile.get_points()
    point_3, point_4 = frac_2.profile.get_points()

    # Compute determinant and numerators for extensions.
    determinant, extend_a, extend_b = \
        find_determinant(point_1, point_2, point_3, point_4)

    # Check for division by zero.
    if abs(determinant) < V_SMALL:
        # Parallel - no intersection!
        intersect = Pointe(-999.9, -999.9)
    else:
        # Get extensions for lines - divide by determinant.
        extend_a /= determinant
        extend_b /= determinant

        # Check if extensions are in-range 0.0 to 1.0.
        # -> otherwise infinite intersection point = no intersection.
        if extend_a < 0.0 or extend_a > 1.0:
            intersect = Pointe(-999.9, -999.9)
        elif extend_b < 0.0 or extend_b > 1.0:
            intersect = Pointe(-999.9, -999.9)

        else:
            # Success - an intersection exists; now compute point on frac #1.
            intersect = find_point_coordinates(point_1, point_2, extend_a)
            success = True

    return success, intersect


#
# -----------------------------------------------------------------------------
# - End of module
# -------1---------2---------3---------4---------5---------6---------7---------8

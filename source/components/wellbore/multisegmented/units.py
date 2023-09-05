"""
Author: Veronika S. Vasylkivska
date: 12/03/2014
"""

def year():
    """ Return a year measured in seconds """
    value = 365*day()
    return value


def month():
    """ Return an average month measured in seconds """
    value = year()/12
    return value


def week():
    """ Return a week measured in seconds """
    value  = 7*day()
    return value


def day():
    """ Return a day measured in seconds """
    value = 24*hour()
    return value


def hour():
    """ Return an hour measured in seconds """
    value = 60*minute()
    return value


def minute():
    """ Return the minute measured in seconds """
    value = 60.0
    return value


def darcy():
    """ Return 1 darcy measured in squared centimeters """
    value = 1.0e+8
    return value


def centimeter():
    """ Return 1 centimeter measured in centimeters """
    value = 1.0
    return value

def micron():
    """ Return 1 micron measured in centimeters """
    value = 1.0e-4
    return value


def meter():
    """ Return 1 meter measured in centimeters """
    value = 100.0
    return value


def cubicmeter():
    """Return 1 cubic meter in cubic centimeters"""
    value = 1000000.0
    return value


def sqrmeter():
    """ Return 1 squared meter in squared centimeters """
    value = 10000.0
    return value

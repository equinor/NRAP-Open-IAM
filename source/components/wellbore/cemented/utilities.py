"""
The module contains the methods for calculations of approximated first and second
derivatives.

Last modified: 12/05/2022
"""
import numpy as np

def first_derivative(f2, f1, dt):
    """ Calculate first derivative. """
    fPrime = (f2-f1)/dt
    return fPrime


def derivative_vector(f, t):
    """ Calculate vector of first derivatives. """
    nElements = np.size(f)
    fPrime = np.zeros(nElements)
    if nElements > 1:
        for ind in range(1, nElements):
            fPrime[ind] = (f[ind]-f[ind-1])/(t[ind]-t[ind-1])
    return fPrime


def second_derivative(f3, f2, f1, dt2, dt1):
    """ Calculate second derivative. """
    fDPrime = (f3-2*f2+f1)/(dt2*dt1)
    return fDPrime


def sec_derivative_vector(f, t):
    """ Calculate vector of second derivatives. """
    nElements = np.size(f)
    fDPrime = np.zeros(nElements)
    if nElements > 2:
        for ind in range(2, nElements):
            fDPrime[ind] = (f[ind]-2*f[ind-1]+f[ind-2])/(
                (t[ind]-t[ind-1])*(t[ind-1]-t[ind-2]))
    return fDPrime

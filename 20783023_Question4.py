"""
Author: Joshua Issa
SID: 20783023
Course: PHYS 375
PS2 - Question 4
"""

import os
import numpy as np
from matplotlib import pyplot as plt

from matplotlib import rc
rc('text', usetex=True)

def part_a(a_path=''):
    """Plots temperature vs s.
    INPUTS:
        a_path - Where to save the plot. If nothing is given then it is displayed.
    OUTPUTS:
        Plots the function.
    """

    Te = 10000
    rho, kappa = 1e-6, 3

    #from lecture 6a : s from 0-1000km --> 0-1e6 meters
    s = np.linspace(0, 1e6, 10000)
    tau = rho * kappa * s

    T = np.power(3/4 * Te**4 * (tau + 2/3), 0.25)

    plt.title(r'$T$ vs $s$')
    plt.ylabel(r'$T$ (K)')
    plt.xlabel(r'$s$ (km)')

    plt.plot(s, T)

    if a_path != '':
        plt.savefig(a_path)
    else: plt.show()

    plt.close('all')

    return T, s

def part_b(T, s, b_path=''):
    """Plots f2 vs s.
    INPUTS:
        T - the temperatures calculated in part A.
        s - the distances from 0 to 1000km from part A.
        b_path - Where to save the plot. If nothing is given then it is displayed.
    OUTPUTS:
        Plots the function.
    """

    def saha(temp):
        rho = 1e-6 #from 4a
        #textbook
        k, me, mp, hbar = 1.381e-23, 9.109e-31, 1.673e-27, 1.055e-34
        eV_to_J = 1.602e-19

        #using the solution to my equation for q3
        return (rho/mp)*np.power(me*k*temp/(2*np.pi * hbar**2), 3/2)*np.exp(-13.6*eV_to_J/(k*temp))

    S_t = saha(T)

    #solution to x^2 + S_t*x - S_t = 0
    f2 = (-S_t - np.sqrt(S_t**2 + 4*S_t))/2 #my f2 is wrong :)

    plt.title(r'$f_2$ vs $s$')
    plt.ylabel(r'$f_2$')
    plt.xlabel(r'$s$ (km)')

    plt.plot(s, f2)

    if b_path != '':
        plt.savefig(b_path)
    else: plt.show()

    plt.close('all')

    """
    At what radius / range of radii (in km) is this fraction largest?
    """

    return f2

def part_c(f2, s, c_path=''):
    """Plots optical depth vs s.
    INPUTS:
        f2 - The frequency of ionized hydrogen (over s) from part B.
        s - the distances from 0 to 1000km from part A.
        b_path - Where to save the plot. If nothing is given then it is displayed.
    OUTPUTS:
        Plots the function.
    """

    rho = 1e-6
    kappa, kbalmer = 3, 3.5e5

    f2 = np.linspace(1, 0, 10000)

    tau_balmer = rho*s*(f2*kbalmer + (1-f2)*kappa) #i think this is wrong
    tau = rho*kappa*s

    plt.title(r'$\tau$ vs $s$')
    plt.ylabel(r'$\tau$')
    plt.xlabel(r'$s$ (km)')

    plt.ylim(0, 1)

    plt.plot(s, tau, label=r'$\tau$')
    plt.plot(s, tau_balmer, label=r'$\tau_{Balmer}$')

    plt.legend()

    if c_path != '':
        plt.savefig(c_path)
    else: plt.show()

    plt.close('all')

    """
    At what radius / range of radii (in km) is this fraction largest?
    """

    return f2


def main():

    here = os.path.dirname(os.path.realpath(__file__))

    a_path = os.path.join(here, 'PS2-Q4a.png')
    T, s = part_a(a_path)

    #getting the wrong graph
    b_path = os.path.join(here, 'PS2-Q4b.png')
    f2 = part_b(T, s, b_path)

    c_path = os.path.join(here, 'PS2-Q4c.png')
    part_c(f2, s, c_path)

if __name__ == '__main__':
    main()

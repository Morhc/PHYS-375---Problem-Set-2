"""
Author: Joshua Issa
SID: 20783023
Course: PHYS 375
PS2 - Question 2
"""

import os
import numpy as np
from matplotlib import pyplot as plt

from matplotlib import rc
rc('text', usetex=True)

def part_a(a_path=''):
    """Plots Blambda vs lambda for T = 5500K.
    INPUTS:
        a_path - Where to save the plot. If nothing is given then it is displayed.
    OUTPUTS:
        Plots the function.
    """

    #Set our constants from the back cover of Ryden
    h, c, k = 6.626e-34, 2.998e8, 1.381e-23

    #Set the temperature
    T = 5500

    #Calculated the maximum of the function by hand
    l_peak = h*c/(5*k*T)

    #Create a series of lambdas from 100nm to 2um
    l = np.linspace(100e-9, 2e-6, 1000)

    #Calculate Blambda
    Bl = 2*h*(c**2)/np.power(l, 5)*np.power(np.exp(h*c/(l*k*T)) - 1, -1)

    plt.title(r'$B_\lambda$ vs $\lambda$')
    plt.xlabel(r'$\lambda$ (m)')
    plt.ylabel(r'$B_\lambda$ (W sr$^{-1}$ m$^{-3}$)')
    plt.ylim(0, np.max(Bl)*1.1)
    plt.xlim(l[0], l[-1])

    plt.plot(l, Bl)

    plt.axvline(l_peak, color='r', ls='--')

    if a_path != '':
        plt.savefig(a_path)
    else: plt.show()

    plt.close('all')

def part_b(b_path=''):
    """Plots Bnu vs nu for T = 5500K.
    INPUTS:
        b_path - Where to save the plot. If nothing is given then it is displayed.
    OUTPUTS:
        Plots the function.
    """

    #Set our constants from the back cover of Ryden
    h, c, k = 6.626e-34, 2.998e8, 1.381e-23

    #Set the temperature
    T = 5500

    #Calculated the maximum of the function by hand
    nu_peak = 3*k*T/h
    corresponding_lambda = c/nu_peak
    print(f'{corresponding_lambda}m')

    #Create a series of nus from 0.15e14 to 30e14 Hz
    nu = np.linspace(0.15e14, 30e14, 1000)

    #Calculate Bnu
    Bnu = (2*h/np.power(c,2)) * np.power(nu, 3) * np.power(np.exp(h*nu/(k*T)) - 1, -1)

    plt.title(r'$B_\nu$ vs $\nu$')
    plt.xlabel(r'$\nu$ (Hz)')
    plt.ylabel(r'$B_\nu$ (W sr$^{-1}$ m$^{-2}$ Hz$^{-1}$)')
    plt.ylim(0, np.max(Bnu)*1.1)
    plt.xlim(nu[0], nu[-1])

    plt.plot(nu, Bnu)

    plt.axvline(nu_peak, color='r', ls='--')

    if b_path != '':
        plt.savefig(b_path)
    else: plt.show()

    plt.close('all')

def part_e(e_path=''):
    """Plots log10(lBl) vs log10(l) for T = 3000K, 5500K, 30000K.
    INPUTS:
        e_path - Where to save the plot. If nothing is given then it is displayed.
    OUTPUTS:
        Plots the function.
    """

    #Set our constants from the back cover of Ryden
    h, c, k = 6.626e-34, 2.998e8, 1.381e-23

    #Create a series of nus from 400 to 800nm
    l = np.linspace(400e-9, 800e-9, 1000)

    #Calculate Bl for each temperature
    Bl = [2*h*(c**2)/np.power(l, 4)*np.power(np.exp(h*c/(l*k*T)) - 1, -1) for T in [3000, 5500, 30000] ]

    plt.title(r'log$_{10}$($B_\lambda$) vs log$_{10}$($\lambda$)')
    plt.xlabel(r'log$_{10}$($\lambda$ (m))')
    plt.ylabel(r'log$_{10}$($B_\lambda$ (W sr$^{-1}$ m$^{-3}$))')

    plt.plot(np.log10(l), np.log10(Bl[0]), label = '$T=3000K$')
    plt.plot(np.log10(l), np.log10(Bl[1]), label = '$T=5500K$')
    plt.plot(np.log10(l), np.log10(Bl[2]), label = '$T=30000K$')

    print(f'The T=3000K star log(Bl) varies by {np.abs(np.log10(Bl[0][-1]) - np.log10(Bl[0][0]))} over the visible wavelengths.')
    print(f'The T=5500K star log(Bl) varies by {np.abs(np.log10(Bl[1][-1]) - np.log10(Bl[1][0]))} over the visible wavelengths.')
    print(f'The T=30000K star log(Bl) varies by {np.abs(np.log10(Bl[2][-1]) - np.log10(Bl[2][0]))} over the visible wavelengths.')

    plt.legend()

    if e_path != '':
        plt.savefig(e_path)
    else: plt.show()

    plt.close('all')

def main():

    here = os.path.dirname(os.path.realpath(__file__))

    a_path = os.path.join(here, 'PS2-Q2a.png')
    part_a(a_path)

    b_path = os.path.join(here, 'PS2-Q2b.png')
    part_b(b_path)

    e_path = os.path.join(here, 'PS2-Q2e.png')
    part_e(e_path)

if __name__ == '__main__':
    main()

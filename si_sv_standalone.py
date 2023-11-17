#!/usr/bin/env python
# coding: utf-8

# This is a python script for coding angle dependent acoustic backscatter
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.special import gamma, factorial, jv
import scipy as sp
import math
import cmath
from sys import exit

# Define important starting variables
c_w = 1485  # m/s speed of sound in water
f_a = 30  # kHz
m_phi = 2 # sand = 1-2; silt = 6.4; clay = 9
# Theta values in degrees
theta_grazing_lower = 70
theta_grazing_upper = 89
theta_grazing_step = .5
# For the integration portion of si, use large upper bound. Step only used if summation used instead
si_integrand_lowerbound = 0  # For si integration
si_integrand_upperbound = 10 ** 3  # For si integration
si_integrand_step = 10 ** -2  # For si integration
si_integrand_vector = np.arange(si_integrand_lowerbound, si_integrand_upperbound, si_integrand_step)  # For si integration

# Define functions to define other important starting variables
def func_rho():
    if -1 <= m_phi < 1:
        return .007797 * m_phi ** 2 - .17057 * m_phi + 2.3139
    elif 1 <= m_phi < 5.3:
        return -.0165406 * m_phi ** 3 + .2290201 * m_phi ** 2 + 1.1069031 * m_phi + 3.0455
    elif 5.3 <= m_phi <= 9:
        return -.0012973 * m_phi + 1.1565
    else:
        raise ValueError('rho: m_phi value not between -1 and 9.')


def func_v():
    if -1 <= m_phi < 1:
        return .002709 * m_phi ** 2 - .056452 * m_phi + 1.2788
    elif 1 <= m_phi < 5.3:
        return -0.0014881 * m_phi ** 3 + 0.0213937 * m_phi ** 2 - 0.1382798 * m_phi + 1.3425
    elif 5.3 <= m_phi <= 9:
        return -0.0024324 * m_phi + 1.0019
    else:
        raise ValueError('v: m_phi value not between -1 and 9.')


def func_gamma():
    if -1 <= m_phi <= 9:
        return 3.25
    else:
        raise ValueError('gamma: m_phi value not between -1 and 9.')


def func_w_2():
    if -1 <= m_phi < 5:
        return .00207 * ((2.03846 - 0.26923 * m_phi) / (1 + .076923 * m_phi)) ** 2
    elif 5 <= m_phi < 9:
        return .0005175
    else:
        raise ValueError('w_2: m_phi value not between -1 and 9.')


def func_kappa_p():
    if -1 <= m_phi < 0:
        return .002709 * m_phi ** 2 - .056452 * m_phi + 1.2788
    elif 0 <= m_phi < 2.6:
        return -0.0014881 * m_phi ** 3 + 0.0213937 * m_phi ** 2 - 0.1382798 * m_phi + 1.3425
    elif 2.6 <= m_phi < 4.5:
        return -0.0024324 * m_phi + 1.0019
    elif 4.5 <= m_phi < 6:
        return -0.0024324 * m_phi + 1.0019
    elif 6 <= m_phi < 9.5:
        return -0.0024324 * m_phi + 1.0019
    elif 9.5 <= m_phi:
        return -0.0024324 * m_phi + 1.0019
    else:
        raise ValueError('kappa_p: m_phi below -1.')


def func_sigma_v():
    if -1 <= m_phi <= 9:
        return 0.004 * alpha_b
    else:
        raise ValueError('sigma_v: m_phi value not between -1 and 9.')


def func_r(tht):
    numerator = rho * v * np.cos(tht) - (1 - v * np.sin(tht) ** 2) ** (1 / 2)
    denominator = rho * c_w * np.cos(tht) + (1 - c_w * np.sin(tht) ** 2) ** (1 / 2)
    return numerator / denominator

# Calculate the remaining starting variables
rho = func_rho()
v = func_v()
gamma_value = func_gamma()  # Spectral exponent... this is a way to measure roughness
w_2 = func_w_2()  # Spectral strength... for 5 < M_phi < 9
kappa_p = func_kappa_p()  # Sediment attenuation constant (dB/m/kHz)
alpha_b = kappa_p * f_a ** 1
sigma_v = func_sigma_v()  # Volume scattering coefficient
alpha = (gamma_value / 2) - 1
sigma_2 = sigma_v / alpha_b
k_a = 2 * np.pi*f_a*(10**3)/c_w/100  # Wave-number in water in cm^-1
c_zeta_2 = (2 * np.pi * w_2 * gamma(2 - alpha) * 2 ** (-2 * alpha)) / (alpha * (1 - alpha) * gamma(1 + alpha))
c_zeta = c_zeta_2 ** .5
# k_c = (gamma_value-2)**(1/(2-gamma_value))/(8*np.pi*w_2*np.cos(tek_a**2)
# zeta = (2*np.pi*w_2*k
r_normal = func_r(0)

# Store the variables in a dictionary
variables = {}
variables['c_w'] = c_w
variables['f_a'] = f_a
variables['k_a'] = k_a
variables['rho'] = rho
variables['v'] = v
variables['gamma_value'] = gamma_value
variables['w_2'] = w_2
variables['kappa_p'] = kappa_p
variables['alpha_b'] = alpha_b
variables['sigma_v'] = sigma_v
variables['alpha'] = alpha
variables['sigma_2'] = sigma_2
variables['r_normal'] = r_normal

# Structure constant


def func_q(tht):  ##This throws some form of error when a 0 is input, it returns a nan
    print(f'func_q tht: {tht}')
    first = np.cos(tht) ** 2
    second = np.sin(tht) ** (-2 * alpha)
    print(f'func_q second: {second}')
    third = c_zeta_2
    fourth = 2 ** (1 - 2 * alpha)
    fifth = k_a ** (2 * (1 - alpha))
    return first * second * third * fourth * fifth


def func_k_c(tht):
    numerator = (gamma_value - 2)
    denominator = (8 * np.pi * w_2 * k_a ** 2 * np.cos(tht) ** .5)
    exponent = (1 / (2 - gamma_value))
    return numerator / denominator ** exponent

    return values


def func_zeta(tht):
    return ((2 * np.pi * w_2 * func_k_c(tht) ** (4 - gamma_value)) / (4 / gamma_value)) ** .5

def func_vf(tht):  # When a 0 is input it returns an infinity
    return (1 - func_r(tht) ** 2) ** 2 * np.cos(tht) ** 2 * (v * np.sin(tht) ** 2) ** -.5


# Transmission term (large-scale interface with small rms slope)
def func_vl(theta):  # right now this is only returning 0s
    coef = 1 / ((np.pi * func_zeta(theta)) * .5)
    integral = integrate.quad(lambda x: func_vf(theta - x) * np.exp(-(x ** 2 / func_zeta(theta) ** 2)), -1*(np.pi / 2 - theta),
                              si_integrand_upperbound)
    return coef * integral[0]


def func_si(theta, lower_bound, upper_bound, step_size):
    # Note alpha is not an input because in this file it is universal variable
    if theta > 0:
        # Define integration by series bounds, samples, and dx
        du = step_size

        # Calculate si
        first = r_normal ** 2
        second = (8*np.pi*np.cos(theta)**2 * np.sin(theta)**2)**-1
        q = func_q(theta)
        # # Using standard integration
        # third = integrate.quad(lambda x: np.exp(-q * x ** (2 * alpha)) * np.i0(x) * x, lower_bound, upper_bound)
        # return first*second*third[0]

        # Using integral of a bessel function
        # third = integrate.quad(lambda u: np.exp(-q * u ** (2 * alpha)) * jv(0, u) * u, lower_bound, upper_bound)
        # return first * second * third[0]

        # Using summation
        sum = 0
        u = lower_bound + step_size
        while u < upper_bound:
            one = jv(0, u)
            two = np.exp(-q*u**(2*alpha))*u
            sum += du*one*two
            u += step_size

        return first*second*sum

    elif theta == 0:
        first = r_normal ** 2
        second = (8 * np.pi * alpha) ** -1
        third = c_zeta ** (-2 / alpha)
        fourth = (2 * k_a ** 2) ** ((alpha - 1) / alpha)
        fifth = gamma(1 / alpha)
        return first * second * third * fourth * fifth

    else:
        raise RuntimeError("Theta is less than zero; cannot compute.")


def func_svl(theta, zeta, r):
    (sigma_2 * func_vl(theta, zeta, r)) / alpha


def s(sis, svls):
    print(f'sis: {sis}')
    print(f'svls: {svls}')
    values = []
    for si, svl in zip(sis, svls):
        values.append(10 * math.log(si + svl))

    return values


def main():
    # 1 to print/plot, 0 to not
    variable_print_flag = True
    svl_plot_flag = True
    si_plot_flag = True
    s_plot_flag = True

    if variable_print_flag is True:  # Line to print variables at the beginning of code
        txt = "Variable List:\n"
        for key, value in variables.items():
            txt +='    ' + key + ": " + str( value) + '\n'

        print(txt)

    # Define Theta range
    theta_grazing_array = np.arange(theta_grazing_lower, theta_grazing_upper, theta_grazing_step)
    theta_incidence_array = 90 - theta_grazing_array
    theta_array = theta_incidence_array*(np.pi/180)
    print(f'theta_array: {theta_array}')
    q_array = np.zeros(len(theta_array))
    for index, theta in enumerate(theta_array):
        q_array[index] = func_q(theta)

    print(f'q_array: {q_array}')


    # Allocate Memory
    svl_array = np.zeros(len(theta_array))
    si_array = np.zeros(len(theta_array))
    s_array = np.zeros(len(theta_array))

    # Calculate si with flag to plot
    for theta_index, theta in enumerate(theta_array):
        si_array[theta_index] = func_si(theta, si_integrand_lowerbound, si_integrand_upperbound, si_integrand_step)

    if si_plot_flag is True:
        x = theta_grazing_array
        y = si_array
        plt.plot(x, y)
        plt.title('Si vs. Grazing Angle')
        plt.xlabel('Grazing Angle (deg)')
        plt.ylabel('Si Value')
        plt.legend()
        plt.grid()
        plt.gcf().autofmt_xdate()  # beautify the x-labels
        plt.show()

    # Calculate svl with flag to plot
    for theta_index, theta in enumerate(theta_array):
        svl_array[theta_index] = sigma_2*func_vl(theta)/(2/(10*np.log10(np.e)))

    if svl_plot_flag is True:
        x = theta_grazing_array
        y = svl_array
        plt.plot(x, y)
        plt.title('Svl vs. Grazing Angle')
        plt.xlabel('Grazing Angle (deg)')
        plt.ylabel('Svl Value')
        plt.legend()
        plt.grid()
        plt.gcf().autofmt_xdate()  # beautify the x-labels
        plt.show()



    # Calculate si with flag to plot
    # Calculate S with a flag to plot

main()




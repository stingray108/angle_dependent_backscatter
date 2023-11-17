# This is a python script for coding angle dependent acoustic backscatter
# This is to test if git commit is working
# A new test
# Here is another test

# 1
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.special import gamma, jv
# import scipy as sp
import math
# import cmath

# Define main constants:
c_w = 1485  # m/s speed of sound in water
f_a = 30  # kHz
k_a = 2 * np.pi*f_a*(10**3)/c_w/100  # Wave-number in water in cm^-1
c_zeta = .1
c_zeta_2 = c_zeta ** 2
r_normal = 1

lowerbound = 0
upperbound = 10**3
step = 10**-1
integrand_vector = np.arange(lowerbound, upperbound, step)


# Define functions:
def func_q(theta, alpha):  # This throws some form of error when a 0 is input, it returns a nan
    first = np.cos(theta)**2
    second = np.sin(theta)**(-2*alpha)
    third = c_zeta_2
    fourth = 2**(1-2*alpha)
    fifth = k_a**(2*(1-alpha))
    return first*second*third*fourth*fifth


def si_integrand(x, q, alpha):
    return np.exp(-q*x**(2*alpha))*jv(0, x)*x


def func_si(theta, q, alpha, lower_bound, upper_bound, step_size):
    if theta > 0:
        # Define integration by series bounds, samples, and dx
        du = step_size

        # Calculate si
        first = r_normal ** 2
        second = (8*np.pi*np.cos(theta)**2 * np.sin(theta)**2)**-1

        # # Using standard integration
        # third = integrate.quad(lambda x: np.exp(-q * x ** (2 * alpha)) * np.i0(x) * x, lower_bound, upper_bound)
        # return first*second*third[0]

        # Using integral of a bessel function
        third = integrate.quad(lambda u: np.exp(-q * u ** (2 * alpha)) * jv(0, u) * u, lower_bound, upper_bound)
        return first * second * third[0]

        # # Using summation
        # sum = 0
        # u = lower_bound + step_size
        # while u < upper_bound:
        #     one = jv(0, u)
        #     two = np.exp(-q*u**(2*alpha))*u
        #     sum += du*one*two
        #     u += step_size
        #
        # return first*second*sum

    elif theta == 0:
        first = r_normal ** 2
        second = (8 * np.pi * alpha) ** -1
        third = c_zeta ** (-2 / alpha)
        fourth = (2 * k_a ** 2) ** ((alpha - 1) / alpha)
        fifth = gamma(1 / alpha)
        return first * second * third * fourth * fifth

    else:
        raise RuntimeError("Theta is less than zero; cannot compute.")


# Define main() function:
def main():
    print(k_a)
    # Define alpha and theta values
    alpha_array = np.array([1, 0.85, 0.75, 0.6, 0.5])

    theta_array_grazing = np.arange(70, 89, .5)
    theta_array_incidence = 90 - theta_array_grazing
    theta_array = theta_array_incidence*(np.pi/180)

    # Allocate memory
    data_points = len(theta_array)
    data_rows = len(alpha_array)
    empty_data_array = np.zeros([data_rows, data_points])
    alpha_q_array = np.copy(empty_data_array)
    si_array = np.copy(empty_data_array)

    # # Iterate through alpha and theta to populate alpha, q data array
    # alpha_index = 0
    # for alpha in alpha_array:
    #     q_index = 0
    #     for theta in theta_array:
    #         alpha_q_array[alpha_index, q_index] = func_q(theta, alpha)
    #         q_index += 1
    #     alpha_index += 1

    #Iterate through alpha and q to populate alpha, Si array
    # Calculate si values
    for alpha_index, alpha in enumerate(alpha_array):
        for theta_index, theta in enumerate(theta_array):
            si_array[alpha_index, theta_index] = \
                func_si(theta, func_q(theta, alpha), alpha, lowerbound, upperbound, step)

    # Convert Si values to dB
    si_array_dB = 10*np.log10(si_array)

    # Plot si values:
    x = theta_array_grazing
    y = 0
    # plt.plot(x, si_array_dB[y, :], label=f'Alpha = {alpha_array[y]}')
    while y < len(alpha_array):
        plt.plot(x, si_array_dB[y, :], label=f'{alpha_array[y]}')
        y += 1
    plt.title('Si Value by Alpha Value')
    plt.xlabel('Grazing Angle (deg)')
    plt.ylabel('Si Value')
    plt.legend()
    plt.grid()
    plt.gcf().autofmt_xdate()  # beautify the x-labels
    plt.show()

    # # Plot integrand sum values:
    # plt.plot(integrand_vector, func_si(theta_array[0], .023, 1, lowerbound, upperbound, step))
    # # while y < len(alpha_array):
    # #     plt.plot(x, alpha_si_array[y, :], label=f'{alpha_array[y]}')
    # #     y += 1
    # plt.title('Integrand sum alpha = 1, q = 3')
    # plt.xlabel('u')
    # plt.ylabel('Cumulative sum for du Value')
    # plt.grid()
    # plt.gcf().autofmt_xdate()  # beautify the x-labels
    # plt.show()


# Run main() function:
main()

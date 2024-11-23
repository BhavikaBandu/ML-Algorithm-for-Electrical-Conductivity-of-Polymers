import numpy as np
import matplotlib.pyplot as plt

# Constants
BOLTZMANN_CONSTANT = 8.617e-05  # Boltzmann constant in eV/K


def calculate_conductivity(T, sigma1, E_a, T1):
    """
    Calculate electrical conductivity using the Arrhenius equation.

    Parameters:
    - T: Temperature at which to calculate the conductivity.
    - sigma1: Electrical conductivity at the reference temperature T1.
    - E_a: Activation energy.
    - T1: Reference temperature.

    Returns:
    - Calculated electrical conductivity at temperature T.
    """
    return sigma1 * np.exp(-E_a / (BOLTZMANN_CONSTANT * T1) * (1 / T - 1 / T1))


def plot_conductivity_temperature(temperature_range, conductivity_values):
    """
    Plot the temperature dependence of electrical conductivity.

    Parameters:
    - temperature_range: Array of temperatures.
    - conductivity_values: Corresponding array of electrical conductivity values.

    Returns:
    - None
    """
    plt.figure(figsize=(8, 6))
    plt.plot(temperature_range, conductivity_values, label='Electrical Conductivity')
    plt.title('Temperature Dependence of Electrical Conductivity')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Electrical Conductivity (S/m)')
    plt.legend()
    plt.show()

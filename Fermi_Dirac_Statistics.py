import numpy as np
import matplotlib.pyplot as plt


def fermi_dirac_distribution(energy, fermi_level, kT):
    return 1 / (1 + np.exp((energy - fermi_level) / kT))


def density_of_states_2d(energy, mass_eff_e, mass_eff_h):
    h_bar = 1.0545718e-34
    pi = np.pi

    dos_e = (mass_eff_e / (pi * h_bar ** 2)) * np.sqrt(2 * mass_eff_e * energy)
    dos_h = (mass_eff_h / (pi * h_bar ** 2)) * np.sqrt(2 * mass_eff_h * energy)

    return dos_e, dos_h


def plot_fermi_dirac_distribution_2d(energy_range, fermi_level_e, fermi_level_h, kT, mass_eff_e, mass_eff_h):
    energy_values = np.linspace(energy_range[0], energy_range[1], 10000)
    dos_e, dos_h = density_of_states_2d(energy_values, mass_eff_e, mass_eff_h)

    peak_e = energy_values[np.argmax(dos_e * fermi_dirac_distribution(energy_values, fermi_level_e, kT))]
    peak_h = energy_values[np.argmax(dos_h * fermi_dirac_distribution(energy_values, fermi_level_h, kT))]
    max_e = np.max(dos_e * fermi_dirac_distribution(energy_values, fermi_level_e, kT))
    max_h = np.max(dos_h * fermi_dirac_distribution(energy_values, fermi_level_h, kT))

    plt.plot([peak_e, peak_e], [0, max_e], color='red', linestyle='dashed', label='Peak')
    plt.plot([peak_h, peak_h], [0, max_h], color='red', linestyle='dashed', label='Peak')
    plt.plot([peak_e, peak_e], [max_e, 0], color='red', linestyle='dashed')
    plt.plot([peak_h, peak_h], [max_h, 0], color='red', linestyle='dashed')

    plt.plot(energy_values, dos_e * fermi_dirac_distribution(energy_values, fermi_level_e, kT), label='Electrons')
    plt.plot(energy_values, dos_h * fermi_dirac_distribution(energy_values, fermi_level_h, kT), label='Holes')
    plt.title('Fermi-Dirac Distribution for 2D Polymeric Sheet')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Carrier Concentration')
    plt.legend()
    plt.grid(True)
    plt.show()

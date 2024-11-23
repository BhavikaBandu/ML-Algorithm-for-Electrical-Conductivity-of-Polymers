from tkinter import *
import numpy as np
import Temperature_Variation
import Fermi_Dirac_Statistics
import Electrical_Conductivity

"""----------------------------------------------------GUI Section--------------------------------------------------"""
def save_input_values():
    with open('input_values.txt', 'w') as file:
        file.write("Fermi-Dirac Statistics:\n")
        file.write(f"Starting Energy (eV): {start_energy.get()}\n")
        file.write(f"Ending Energy (eV): {end_energy.get()}\n")
        file.write(f"Fermi Level (Electrons) (eV): {fermi_level_e.get()}\n")
        file.write(f"Fermi Level (Holes) (eV): {fermi_level_h.get()}\n")
        file.write(f"Temperature (K): {temperature.get()}\n\n")

        file.write("Electrical Conductivity:\n")
        file.write(f"Relative Permittivity: {epsilon_r.get()}\n")
        file.write(f"Electron Mobility (m²/Vs): {mu_n.get()}\n")
        file.write(f"Hole Mobility (m²/Vs): {mu_p.get()}\n")
        file.write(f"Density of States (Conduction Band) (cm^-3): {N_c.get()}\n")
        file.write(f"Density of States (Valence Band) (cm^-3): {N_v.get()}\n")
        file.write(f"Dopant Concentration (cm^-3): {doping_concentration.get()}\n")
        file.write(f"Dopant Type: {dopant_type.get()}\n")
        file.write(f"Dopant Energy Level (eV): {dopant_energy_level.get()}\n")
        file.write(f"Dopant Bandgap (eV): {dopant_bandgap.get()}\n")
        file.write(f"Effective Mass (Dopant) (kg): {effective_mass_dopant.get()}\n\n")

        file.write("Temperature Variation:\n")
        file.write(f"Initial Temperature (K): {T1.get()}\n")
        file.write(f"Electrical Conductivity at T1 (S/m): {sigma1.get()}\n")
        file.write(f"Activation Energy (eV): {E_a.get()}\n")
        file.write(f"Starting Temperature: {StartTemp.get()}\n")
        file.write(f"Ending Temperature: {EndTemp.get()}\n")
        file.write(f"Step Size for the Temperature Range (eV): {TempStep.get()}\n")

        window.destroy()


window = Tk()
window.geometry('1200x700')
window.title('Conductivity Calculator')
window.config(background='#00b3ca')

# Center the window on the screen
window.eval('tk::PlaceWindow . center')

# Welcome Labels
label = Label(window, text='Welcome to Conductivity Calculator!',
              font=('Eras Bold ITC', 30),
              fg='black',
              bg='yellow',
              padx=20,
              pady=20)
label1 = Label(window, text='where you input values and we calculate the conductivity of the polymer for you :)',
               font=('Eras Bold ITC', 20),
               fg='black',
               bg='#e38690')

label.grid(row=0, column=0, columnspan=6, pady=20)
label1.grid(row=1, column=0, columnspan=6, pady=10)

# Group 1: Fermi-Dirac Statistics
Label(window, text='Fermi-Dirac Statistics', font=('Eras Bold ITC', 20)).grid(row=2, column=0, pady=10, columnspan=2,
                                                                              sticky='w')
Label(window, text='Starting Energy (eV):').grid(row=3, column=0, padx=5, pady=5, sticky="e")
start_energy = Entry(window, width=15, font=('Ebrima', 14))
start_energy.grid(row=3, column=1, pady=5, sticky="w")
Label(window, text='Ending Energy (eV):').grid(row=4, column=0, padx=5, pady=5, sticky="e")
end_energy = Entry(window, width=15, font=('Ebrima', 14))
end_energy.grid(row=4, column=1, pady=5, sticky="w")
Label(window, text='Fermi Level (Electrons) (eV):').grid(row=5, column=0, padx=5, pady=5, sticky="e")
fermi_level_e = Entry(window, width=15, font=('Ebrima', 14))
fermi_level_e.grid(row=5, column=1, pady=5, sticky="w")
Label(window, text='Fermi Level (Holes) (eV):').grid(row=6, column=0, padx=5, pady=5, sticky="e")
fermi_level_h = Entry(window, width=15, font=('Ebrima', 14))
fermi_level_h.grid(row=6, column=1, pady=5, sticky="w")
Label(window, text='Temperature (K):').grid(row=7, column=0, padx=5, pady=5, sticky="e")
temperature = Entry(window, width=15, font=('Ebrima', 14))
temperature.grid(row=7, column=1, pady=5, sticky="w")

# Group 2: Electrical Conductivity
Label(window, text='Electrical Conductivity', font=('Eras Bold ITC', 20)).grid(row=2, column=2, pady=10, columnspan=2,
                                                                               sticky='w')
Label(window, text='Relative Permittivity:').grid(row=3, column=2, padx=5, pady=5, sticky="e")
epsilon_r = Entry(window, width=15, font=('Ebrima', 14))
epsilon_r.grid(row=3, column=3, pady=5, sticky="w")
Label(window, text='Electron Mobility (m²/Vs):').grid(row=4, column=2, padx=5, pady=5, sticky="e")
mu_n = Entry(window, width=15, font=('Ebrima', 14))
mu_n.grid(row=4, column=3, pady=5, sticky="w")
Label(window, text='Hole Mobility (m²/Vs):').grid(row=5, column=2, padx=5, pady=5, sticky="e")
mu_p = Entry(window, width=15, font=('Ebrima', 14))
mu_p.grid(row=5, column=3, pady=5, sticky="w")
Label(window, text='Density of States (Conduction Band) (cm^-3):').grid(row=6, column=2, padx=5, pady=5, sticky="e")
N_c = Entry(window, width=15, font=('Ebrima', 14))
N_c.grid(row=6, column=3, pady=5, sticky="w")
Label(window, text='Density of States (Valence Band) (cm^-3):').grid(row=7, column=2, padx=5, pady=5, sticky="e")
N_v = Entry(window, width=15, font=('Ebrima', 14))
N_v.grid(row=7, column=3, pady=5, sticky="w")
Label(window, text='Dopant Concentration (cm^-3):').grid(row=8, column=2, padx=5, pady=5, sticky="e")
doping_concentration = Entry(window, width=15, font=('Ebrima', 14))
doping_concentration.grid(row=8, column=3, pady=5, sticky="w")
Label(window, text='Dopant Type:').grid(row=9, column=2, padx=5, pady=5, sticky="e")
dopant_type = Entry(window, width=15, font=('Ebrima', 14))
dopant_type.grid(row=9, column=3, pady=5, sticky="w")
Label(window, text='Dopant Energy Level (eV):').grid(row=10, column=2, padx=5, pady=5, sticky="e")
dopant_energy_level = Entry(window, width=15, font=('Ebrima', 14))
dopant_energy_level.grid(row=10, column=3, pady=5, sticky="w")
Label(window, text='Dopant Bandgap (eV):').grid(row=11, column=2, padx=5, pady=5, sticky="e")
dopant_bandgap = Entry(window, width=15, font=('Ebrima', 14))
dopant_bandgap.grid(row=11, column=3, pady=5, sticky="w")
Label(window, text='Effective Mass (Dopant) (kg):').grid(row=12, column=2, padx=5, pady=5, sticky="e")
effective_mass_dopant = Entry(window, width=15, font=('Ebrima', 14))
effective_mass_dopant.grid(row=12, column=3, pady=5, sticky="w")

# Group 3: Temperature Variation
Label(window, text='Temperature Variation', font=('Eras Bold ITC', 20)).grid(row=2, column=4, pady=10, columnspan=2,
                                                                             sticky='w')
Label(window, text='Initial Temperature (K):').grid(row=3, column=4, padx=5, pady=5, sticky="e")
T1 = Entry(window, width=15, font=('Ebrima', 14))
T1.grid(row=3, column=5, pady=5, sticky="w")
Label(window, text='Electrical Conductivity at T1 (S/m):').grid(row=4, column=4, padx=5, pady=5, sticky="e")
sigma1 = Entry(window, width=15, font=('Ebrima', 14))
sigma1.grid(row=4, column=5, pady=5, sticky="w")
Label(window, text='Activation Energy (eV):').grid(row=5, column=4, padx=5, pady=5, sticky="e")
E_a = Entry(window, width=15, font=('Ebrima', 14))
E_a.grid(row=5, column=5, pady=5, sticky="w")
Label(window, text='Starting Temperature (K):').grid(row=6, column=4, padx=5, pady=5, sticky="e")
StartTemp = Entry(window, width=15, font=('Ebrima', 14))
StartTemp.grid(row=6, column=5, pady=5, sticky="w")
Label(window, text='Ending Temperature (K):').grid(row=7, column=4, padx=5, pady=5, sticky="e")
EndTemp = Entry(window, width=15, font=('Ebrima', 14))
EndTemp.grid(row=7, column=5, pady=5, sticky="w")
Label(window, text='Step Size Of Temperature Range:').grid(row=8, column=4, padx=5, pady=5, sticky="e")
TempStep = Entry(window, width=15, font=('Ebrima', 14))
TempStep.grid(row=8, column=5, pady=5, sticky="w")


# Save Inputs Button
Button(window, text='Save Inputs', font=('Ebrima', 16), command=save_input_values).grid(row=13, column=3, pady=20)

window.mainloop()
"""-----------------------------------------------------------------------------------------------------------------"""

"""-----------------------------------------------Fermi-Dirac Statistics--------------------------------------------"""
# Read input values from the text file
with open('input_values.txt', 'r') as file:
    lines = file.readlines()

# Extract relevant parameters
start_energy = float(lines[1].split(":")[1].strip())
end_energy = float(lines[2].split(":")[1].strip())
energy_range = (start_energy, end_energy)
fermi_level_e = float(lines[3].split(":")[1].strip())
fermi_level_h = float(lines[4].split(":")[1].strip())
temperature = float(lines[5].split(":")[1].strip())
kT = 8.6173e-5 * temperature
effective_mass_e = 9.1094e-31  # Effective mass of an electron in kg
effective_mass_h = 0.5 * effective_mass_e  # Effective mass of hole in kg

# Call the plotting function with the read parameters
Fermi_Dirac_Statistics.plot_fermi_dirac_distribution_2d(energy_range, fermi_level_e, fermi_level_h, kT, effective_mass_e, effective_mass_h)
"""------------------------------------------------------------------------------------------------------------------"""

"""-------------------------------------------------Temperature Variation-------------------------------------------"""
# Read input values from the text file
with open('input_values.txt', 'r') as file:
    lines = file.readlines()

# Extract relevant parameters
T1 = float(lines[20].split(":")[1].strip())
sigma1 = float(lines[21].split(":")[1].strip())
E_a = float(lines[22].split(":")[1].strip())
start_temp = float(lines[23].split(":")[1].strip())
end_temp = float(lines[24].split(":")[1].strip())
step_temp = float(lines[25].split(":")[1].strip())

# Generate a list of temperatures based on user input
temperature_range = np.arange(start_temp, end_temp + 1, step_temp)

# Calculate conductivity for each temperature using the Arrhenius equation
Temperature_Variation.conductivity_values = [Temperature_Variation.calculate_conductivity(T, sigma1, E_a, T1) for T in temperature_range]

# Plotting
Temperature_Variation.plot_conductivity_temperature(temperature_range, Temperature_Variation.conductivity_values)
"""-----------------------------------------------------------------------------------------------------------------"""

"""-----------------------------------------------Electrical Conductivity--------------------------------------------"""
sigma = Electrical_Conductivity.calculate_electrical_conductivity(Electrical_Conductivity.n, Electrical_Conductivity.p)
Electrical_Conductivity.visualize_results(sigma)
total_conductivity = Electrical_Conductivity.calculate_total_conductivity(sigma)

print("Total Conductivity:", total_conductivity, "S/m")
"""------------------------------------------------------------------------------------------------------------------"""

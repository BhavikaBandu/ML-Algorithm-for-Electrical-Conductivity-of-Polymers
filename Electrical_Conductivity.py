import numpy as np
import matplotlib.pyplot as plt


def einstein_relation(mu):
    Dn = mu * k_B * T / q
    return Dn


def calculate_diffusion_coefficients():
    D_n = einstein_relation(mu_n)
    D_p = einstein_relation(mu_p)
    return D_n, D_p


def calculate_dos_dopant():
    energy_values = np.linspace(0, 1.1, 1000)
    dos_dopant = ((effective_mass_dopant / (np.pi * h_bar ** 2)) *
                  np.sqrt(np.maximum(0, 2 * effective_mass_dopant * (
                          energy_values - (dopant_energy_level - dopant_bandgap)))))
    return dos_dopant


def initialize_simulation_grid():
    n = np.ones((Nx, Ny)) * 1e10 + np.random.rand(Nx, Ny) * 1e8
    p = np.ones((Nx, Ny)) * 1e10 + np.random.rand(Nx, Ny) * 1e8
    phi = np.zeros((Nx, Ny))
    return n, p, phi


def introduce_dopant_effects(n, p):
    if dopant_type == 'donor':
        n += doping_concentration
        average_dos_dopant = np.mean(dos_dopant)
        n += average_dos_dopant * np.exp(-(q * dopant_energy_level) / (k_B * T))
    elif dopant_type == 'acceptor':
        p += doping_concentration
        average_dos_dopant = np.mean(dos_dopant)
        p += average_dos_dopant * np.exp((q * dopant_energy_level) / (k_B * T))
    return n, p


def update_carrier_concentrations(n, p, phi, D_n, D_p):
    n[1:-1, 1:-1] += dt * (D_n / dx ** 2) * (n[2:, 1:-1] - 2 * n[1:-1, 1:-1] + n[:-2, 1:-1]) + \
                     dt * (D_n / dy ** 2) * (n[1:-1, 2:] - 2 * n[1:-1, 1:-1] + n[1:-1, :-2]) - \
                     dt * q * n[1:-1, 1:-1] / (epsilon_r * epsilon_0)

    p[1:-1, 1:-1] += dt * (D_p / dx ** 2) * (p[2:, 1:-1] - 2 * p[1:-1, 1:-1] + p[:-2, 1:-1]) + \
                     dt * (D_p / dy ** 2) * (p[1:-1, 2:] - 2 * p[1:-1, 1:-1] + p[1:-1, :-2]) - \
                     dt * q * p[1:-1, 1:-1] / (epsilon_r * epsilon_0)

    rho = q * (n - p)
    phi[1:-1, 1:-1] += dt * (1 / (epsilon_r * epsilon_0)) * (rho[1:-1, 1:-1] - (n[1:-1, 1:-1] - p[1:-1, 1:-1]))

    return n, p, phi


def calculate_electrical_conductivity(n, p):
    sigma = np.abs(q * (n - p) * (mu_n + mu_p) / (n + p))
    return sigma


def visualize_results(sigma):
    plt.figure(figsize=(12, 8))

    ax1 = plt.subplot(1, 3, 1, projection='3d')
    X, Y = np.meshgrid(np.linspace(0, L, Nx), np.linspace(0, W, Ny))
    surface = ax1.plot_surface(X, Y, sigma * 1e19, cmap='cool', rstride=5, cstride=5, alpha=0.8, linewidth=0.5,
                               antialiased=True)
    ax1.set_xlabel('X-axis (m)')
    ax1.set_ylabel('Y-axis (m)')
    ax1.set_title('3D Surface Plot: Electrical Conductivity')
    plt.colorbar(surface, ax=ax1, label='Electrical Conductivity (S/m)')

    plt.subplot(1, 2, 2)
    image = plt.imshow(sigma * 1e19, extent=(0, L, 0, W), origin='lower', cmap='cool')
    plt.colorbar(image, label='Electrical Conductivity (S/m)')
    plt.xlabel('X-axis (m)')
    plt.ylabel('Y-axis (m)')
    plt.title(f'Image Plot: Electrical Conductivity of 2D Polymer Sheet with {dopant_type} doping')
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)

    plt.show()


def calculate_total_conductivity(sigma):
    total_conductivity = np.mean(sigma) * 10e18
    return total_conductivity


# Constants
q = 1.6e-19
epsilon_0 = 8.85e-12
k_B = 1.38e-23
T = 300
h_bar = 1.0545718e-34

# Simulation parameters
L = 1e-6
W = 1e-6
Nx = 100
Ny = 100
Nt = 1000
dt = 1e-16
dx = L / (Nx - 1)
dy = W / (Ny - 1)

# Material parameters
with open('input_values.txt', 'r') as file:
    lines = file.readlines()

epsilon_r = float(lines[8].split(":")[1].strip())
mu_n = float(lines[9].split(":")[1].strip())
mu_p = float(lines[10].split(":")[1].strip())
N_c = float(lines[11].split(":")[1].strip())
N_v = float(lines[12].split(":")[1].strip())
doping_concentration = float(lines[13].split(":")[1].strip())
dopant_type = lines[14].split(":")[1].strip()
dopant_energy_level = float(lines[15].split(":")[1].strip())
dopant_bandgap = float(lines[16].split(":")[1].strip())
effective_mass_dopant = float(lines[17].split(":")[1].strip())

D_n, D_p = calculate_diffusion_coefficients()
dos_dopant = calculate_dos_dopant()
n, p, phi = initialize_simulation_grid()
n, p = introduce_dopant_effects(n, p)

for t in range(Nt):
    n, p, phi = update_carrier_concentrations(n, p, phi, D_n, D_p)
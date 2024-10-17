import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_simulation_data(filename):
    steps = []
    kinetic_energy = []
    potential_energy = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('Step:'):
                data = line.split()
                steps.append(int(data[1]))
                kinetic_energy.append(float(data[2]))
                potential_energy.append(float(data[3]))
    
    return np.array(steps), np.array(kinetic_energy), np.array(potential_energy)

def plot_energy_conservation(steps, kinetic_energy, potential_energy, method_name):
    total_energy = kinetic_energy + potential_energy
    
    plt.figure(figsize=(12, 8))
    plt.plot(steps, kinetic_energy, label='Kinetic Energy')
    plt.plot(steps, potential_energy, label='Potential Energy')
    plt.plot(steps, total_energy, label='Total Energy')
    plt.xlabel('Simulation Step')
    plt.ylabel('Energy (J)')
    plt.title(f'Energy Conservation with {method_name} method')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f'energy_conservation_{method_name.lower()}.png', dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot energy conservation from simulation data.')
    parser.add_argument('filename', type=str, help='Input filename containing simulation data')
    parser.add_argument('method', type=str, help='Integration method used in the simulation')
    args = parser.parse_args()

    steps, kinetic_energy, potential_energy = read_simulation_data(args.filename)
    plot_energy_conservation(steps, kinetic_energy, potential_energy, args.method)

if __name__ == '__main__':
    main()
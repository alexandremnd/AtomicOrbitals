import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from dataclasses import dataclass

@dataclass
class GaussianPrimitive:
    i: int
    j: int
    k: int
    alpha: float
    x0: float
    y0: float
    z0: float

def get_normalization(gaussian: GaussianPrimitive):



def create_grid(limits=(-5, 5), points=100):
    """Create a 3D grid for evaluating orbitals"""
    x = np.linspace(limits[0], limits[1], points)
    y = np.linspace(limits[0], limits[1], points)
    z = np.linspace(limits[0], limits[1], points)

    # Create meshgrid
    X, Y, Z = np.meshgrid(x, y, z)

    return X, Y, Z

def plot_orbital_density(orbital_densities):
    """Plot the electron density for each orbital"""
    n_orbitals = len(orbital_densities)
    fig = plt.figure(figsize=(5*n_orbitals, 5))

    for i, (orbital_idx, density) in enumerate(orbital_densities.items()):
        # For 3D plotting we typically need to slice or project the data
        # Here we'll take a slice at z=0
        z_mid = density.shape[2] // 2
        density_slice = density[:, :, z_mid]

        ax = fig.add_subplot(1, n_orbitals, i+1)
        im = ax.imshow(density_slice.T, origin='lower',
                       extent=[-5, 5, -5, 5], cmap='viridis')
        plt.colorbar(im, ax=ax)
        ax.set_title(f'Orbital {orbital_idx}')
        ax.set_xlabel('x (Bohr)')
        ax.set_ylabel('y (Bohr)')

    plt.tight_layout()
    plt.savefig('orbital_densities.png', dpi=300)
    plt.show()

def main():
    # Parse the UGBS basis set
    basis_file = './data/BasisSet/ugbs.basis'
    elements_basis: list[GaussianPrimitive] = parse_basis(basis_file)
    grid_points = create_grid(limits=(-10, 10), points=100)

    orbitals_densities = np.zeros((len(elements_basis), grid_points.shape[0], grid_points.shape[1], grid_points.shape[2]))

    for i, orbital in enumerate(elements_basis):
        orbital_densities[i, :, :, :] =


    orbital_densities = calculate_orbital_density(basis_functions, coefficients, grid_points)

    plot_orbital_density(orbital_densities)

if __name__ == "__main__":
    main()
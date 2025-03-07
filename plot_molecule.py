import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from dataclasses import dataclass
from scipy import special

@dataclass
class Atom:
    Z: int
    x: float
    y: float
    z: float

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
    fact_x_exponent = special.gamma(gaussian.i + 1)
    fact_2x_exponent = special.gamma(2 * gaussian.i + 1)

    fact_y_exponent = special.gamma(gaussian.j + 1)
    fact_2y_exponent = special.gamma(2 * gaussian.j + 1)

    fact_z_exponent = special.gamma(gaussian.k + 1)
    fact_2z_exponent = special.gamma(2 * gaussian.k + 1)

    t1 = np.power(2 * gaussian.alpha / np.pi, 3.0 / 4)
    t2 = np.sqrt(np.power(8 * gaussian.alpha, gaussian.i + gaussian.j + gaussian.k))
    t3 = fact_x_exponent * fact_y_exponent * fact_z_exponent
    t4 = fact_2x_exponent * fact_2y_exponent * fact_2z_exponent
    return t1 * t2 * t3 / t4


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
    atoms: list[Atom] = parse_system("data/system.out")
    basis: list[GaussianPrimitive] = parse_basis("data/basis.out")
    coefficient: np.ndarray = parse_coefficients("data/coefficients.out")
    X, Y, Z = create_grid(limits=(-10, 10), points=100)

    orbitals_densities = np.zeros((len(basis), X.shape[0], Y.shape[0], Z.shape[0]))

    for i, orbital in enumerate(basis):
        normalization = get_normalization(orbital)
        orbitals_densities[i, :, :, :] = normalization * (X - orbital.x0)**orbital.i * (Y - orbital.y0)**orbital.j * (Z - orbital.z0)**orbital.k * np.exp(-orbital.alpha * ((X - orbital.x0)**2 + (Y - orbital.y0)**2 + (Z - orbital.z0)**2))

    fig = plt.figure(figsize=(8, 4))



if __name__ == "__main__":
    main()
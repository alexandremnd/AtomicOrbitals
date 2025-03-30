import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from dataclasses import dataclass
from scipy import special
from matplotlib.gridspec import GridSpec

plt.style.use("article.mplstyle")

@dataclass
class Atom:
    Z: int
    x: float
    y: float
    z: float

@dataclass
class Orbital:
    grid_value: np.ndarray

def load_file(file_name: str) -> str:
    """Load a file"""
    content = None
    with open(file_name, 'r') as f:
        content = f.readlines()

    return "".join(content).split("****\n")

def parse_system(atom_description: str) -> list[Atom]:
    """Parse the system file"""
    atoms = []

    atom_description = atom_description.split("\n")
    for line in atom_description[:-1]:
        Z, x, y, z = map(float, line.split())
        atoms.append(Atom(Z, x, y, z))

    return atoms

def get_lim(atoms: list[Atom]):
    x_min = min(atom.x for atom in atoms)
    x_max = max(atom.x for atom in atoms)
    y_min = min(atom.y for atom in atoms)
    y_max = max(atom.y for atom in atoms)
    z_min = min(atom.z for atom in atoms)
    z_max = max(atom.z for atom in atoms)

    lx = x_max - x_min
    ly = y_max - y_min
    lz = z_max - z_min

    # Keep the same aspect ratio as the z axis
    x_min = x_min - (lz - lx) / 2
    x_max = x_max + (lz - lx) / 2
    y_min = y_min - (lz - ly) / 2
    y_max = y_max + (lz - ly) / 2

    return x_min - 1, x_max + 1, y_min - 1, y_max + 1, z_min - 1, z_max + 1

def parse_basis(basis_description: list[str], X, Y, Z) -> list[Orbital]:
    """Parse the basis file"""
    basis = []

    for contracted_decomp in basis_description[:-1]:
        orb_grid = np.zeros_like(X)
        for line in contracted_decomp.split("\n")[:-1]:
            coeff, constant, x, y, z, alpha, i, j, k = map(float, line.split())
            orb_grid += coeff * constant * (X - x)**i * (Y - y)**j * (Z - z)**k * np.exp(-alpha * ((X - x)**2 + (Y - y)**2 + (Z - z)**2))
        basis.append(Orbital(orb_grid))

    return basis

def compute_density(basis: list[Orbital], density_matrix, coefficients, X, Y, Z):
    density = np.zeros_like(basis[0].grid_value)
    n_basis = len(basis)
    points = X.shape[0]

    for i in range(n_basis):
        for j in range(n_basis):
            density += density_matrix[i, j] * basis[i].grid_value * basis[j].grid_value

    return density

def compute_orbital_density(basis: list[Orbital], orbital_index, coefficients, X, Y, Z):
    density = np.zeros_like(basis[0].grid_value)
    n_basis = len(basis)
    points = X.shape[0]

    for i in range(n_basis):
        for j in range(n_basis):
            density += coefficients[i, orbital_index] * coefficients[j, orbital_index] * basis[i].grid_value * basis[j].grid_value

    return density

def create_grid(limits=(-5, 5, -5, 5, -5, 5), points=100):
    """Create a 3D grid for evaluating orbitals"""
    x = np.linspace(limits[0], limits[1], points)
    y = np.linspace(limits[2], limits[3], points)
    z = np.linspace(limits[4], limits[5], points)

    # Create meshgrid
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    return X, Y, Z

def plot(filename, density, X, Y, Z):
    fig = plt.figure(figsize=(5, 10))
    gs = GridSpec(3, 1, figure=fig)
    ax1: plt.Axes = fig.add_subplot(gs[0, 0])
    ax2: plt.Axes = fig.add_subplot(gs[1, 0])
    ax3: plt.Axes = fig.add_subplot(gs[2, 0])

    vmin = density.min()
    vmax = density.max()
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    ax1.set_title("Electronic density ($xy$ plane)")
    ax1.set_xlabel("$x ~ [a_0]$")
    ax1.set_ylabel("$y ~ [a_0]$")
    # cf1 = ax1.contourf(X[:, :, 62], Y[:, :, 62], density[:, :, 62], levels=30, antialiased=True)
    extent = (X[:, :, 50].min(), X[:, :, 50].max(), Y[:, :, 50].min(), Y[:, :, 50].max())
    cf1 = ax1.imshow(density[:, :, 50].T, extent=extent, interpolation="bilinear", norm=norm, origin="lower", aspect="auto")

    ax2.set_title("Electronic density ($xz$ plane)")
    ax2.set_xlabel("$x ~ [a_0]$")
    ax2.set_ylabel("$z ~ [a_0]$")
    # cf2 = ax2.contourf(X[:, 50, :], Z[:, 50, :], density[:, 50, :], levels=30, antialiased=True)
    extent = (X[:, 50, :].min(), X[:, 50, :].max(), Z[:, 50, :].min(), Z[:, 50, :].max())
    cf2 = ax2.imshow(density[:, 50, :].T, extent=extent, interpolation="bilinear", norm=norm, origin="lower", aspect="auto")

    ax3.set_title("Electronic density ($yz$ plane)")
    ax3.set_xlabel("$y ~ [a_0]$")
    ax3.set_ylabel("$z ~ [a_0]$")
    # cf3 = ax3.contourf(Y[50, :, :], Z[50, :, :], density[50, :, :], levels=30, antialiased=True)
    extent = (Y[50, :, :].min(), Y[50, :, :].max(), Z[50, :, :].min(), Z[50, :, :].max())
    cf3 = ax3.imshow(density[50, :, :].T, extent=extent, interpolation="bilinear", norm=norm, origin="lower", aspect="auto")

    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cbar = fig.colorbar(cf1, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Electron Density')
    plt.subplots_adjust(right=0.8, hspace=0.4)

    # plt.tight_layout()
    # plt.savefig(f"figures/{filename}.svg", bbox_inches="tight")
    plt.show()


def main(molecule_name):
    system = load_file(f"data/{molecule_name}/system.out")
    atoms: list[Atom] = parse_system(system[0])
    lim = get_lim(atoms)
    print(f"System limits: {lim}")

    X, Y, Z = create_grid(limits=lim, points=101)
    basis: list[Orbital] = parse_basis(system[1:], X, Y, Z)
    density_matrix = np.loadtxt(f"data/{molecule_name}/density.out")
    coefficients = np.loadtxt(f"data/{molecule_name}/coefficients.out")

    type = input("Enter the type of plot (Total density: 0 | Orbital density: 1): ")
    if type == "0":
        density = compute_density(basis, density_matrix, coefficients, X, Y, Z)
        plot(f"{molecule_name}_density", density, X, Y, Z)
    elif type == "1":
        for i in range(coefficients.shape[1]):
            print(f"Plotting orbital {i}")
            density = compute_orbital_density(basis, i, coefficients, X, Y, Z)
            plot(f"{molecule_name}_orbital_{i}", density, X, Y, Z)


if __name__ == "__main__":
    molecule = input("Enter the molecule name (e.g., h2): ")
    main(molecule)
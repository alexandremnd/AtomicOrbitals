import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

a0 = 5.29177210544e-11 # [m]
a0_angstrom = a0 * 1e10 # [Å]
a0_pm = a0 * 1e12 # [pm]

data = np.loadtxt('data/molecular_optimization.out', delimiter=',')

plt.style.use("article.mplstyle")

minimal_energy = np.min(data[:,1])
minimal_length = data[np.argmin(data[:,1]),0]

plt.plot(data[:,0] * a0_pm, data[:,1], "o")
plt.xlabel('Bond length [pm]')
plt.ylabel('Energy (Hartree)')
plt.title(f"Minimal energy: {minimal_energy*27.2114:.5f} eV at {minimal_length * a0_pm:.5f} pm")
plt.tight_layout()
plt.savefig("molecule.svg")
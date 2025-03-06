import numpy as np
import matplotlib.pyplot as plt

exit()


data = np.loadtxt('data/atom_energies.out')

plt.plot(data[:,0], data[:,1] , 'o-')
plt.xlabel('Number of atoms')
plt.ylabel('Energy (Hartree)')
plt.show()
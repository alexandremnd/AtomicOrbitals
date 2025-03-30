import numpy as np
import matplotlib.pyplot as plt

a0 = 5.29177210544e-11 # [m]
a0_angstrom = a0 * 1e10 # [Å]
a0_pm = a0 * 1e12 # [pm]

data = np.loadtxt('data/h2_6-31g.out', delimiter=',')
data1 = np.loadtxt('data/h2_sto6g.out', delimiter=',')
data2 = np.loadtxt('data/h2_3-21g.out', delimiter=',')
data3 = np.loadtxt('data/h2_rhf.out', delimiter=',')

plt.style.use("article.mplstyle")

minimal_energy = np.min(data[:,1])
minimal_length = data[np.argmin(data[:,1]),0]

plt.plot(data[:,0] * a0_pm, data[:,1], markersize=3, label="6-31G")
plt.plot(data1[:,0] * a0_pm, data1[:,1], markersize=3, label="STO-6G")
plt.plot(data2[:,0] * a0_pm, data2[:,1], markersize=3, label="3-21G")
plt.plot(data3[:,0] * a0_pm, data3[:,1], markersize=3, label="RHF 6-31G")

plt.xlabel('Bond length [pm]')
plt.ylabel('Energy (Hartree)')
plt.title(f"Minimal energy: {minimal_energy:.3f} Hartree at {minimal_length * a0_pm:.2f} pm")

plt.tight_layout()
plt.legend()
plt.savefig("figures/h2_optimization.svg")
plt.show()
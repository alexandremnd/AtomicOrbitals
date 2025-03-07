import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("article.mplstyle")


data1 = np.loadtxt('data/sto6g-314s.out', delimiter=',')
data2 = np.loadtxt('data/ugbs-7360s.out', delimiter=',')
max_Z = data1[-1,0]


df = pd.read_csv("data/experimental_binding_energies.csv", delimiter=",")
neutral_atom = df[(df["Ion Charge"] == 0.0) & (df["At. num"] % 2 == 0) & (df["At. num"] <= max_Z)]

plt.plot(data1[:,0], -data1[:,1], "o", label="Predicted ground state (STO-6G)")
plt.plot(data2[:,0], -data2[:,1], "o", label="Predicted ground state (UGBS)")
plt.errorbar(neutral_atom["At. num"], neutral_atom["Binding Energy (a) (eV)"]/27.2114, yerr=neutral_atom["Uncertainty (b) (eV)"]/27.2114, fmt='-', label="Experimental ground state")
plt.plot(data1[:,0], (neutral_atom["Binding Energy (a) (eV)"]/27.2114 + data1[:,1])/neutral_atom["Binding Energy (a) (eV)"] * 100 * 27.2114, "o", label="STO-6G Error")
plt.plot(data2[:-1,0], (neutral_atom["Binding Energy (a) (eV)"]/27.2114 + data2[:-1,1])/neutral_atom["Binding Energy (a) (eV)"] * 100 * 27.2114, "o", label="UGBS Error")
plt.legend()
plt.xlabel('Atomic number $Z$')
plt.ylabel('Energy (Hartree)')
plt.yscale("log")
plt.savefig("ayo.svg")
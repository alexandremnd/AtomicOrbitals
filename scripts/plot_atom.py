import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("article.mplstyle")

# Load data
data1 = np.loadtxt('data/atom-6-31g.out', delimiter=',')
data3 = np.loadtxt('data/atom-sto6g.out', delimiter=',')
max_Z = data1[-1,0]

# Load experimental data
df = pd.read_csv("data/experimental_binding_energies.csv", delimiter=",")
neutral_atom = df[(df["Ion Charge"] == 0.0) & (df["At. num"] % 2 == 0) & (df["At. num"] <= max_Z)]

# Create figure and primary axis
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot binding energies on the primary y-axis
ax1.plot(data3[:,0], -data3[:,1], "o", label="Predicted (STO-6G)", color='blue')
ax1.plot(data1[:,0], -data1[:,1], "o", label="Predicted (6-31G)", color='red')
ax1.errorbar(neutral_atom["At. num"],
             neutral_atom["Binding Energy (a) (eV)"]/27.2114,
             yerr=neutral_atom["Uncertainty (b) (eV)"]/27.2114,
             fmt='-', label="Experimental", color='black')

# Set up primary axis labels
ax1.set_xlabel('Atomic number $Z$')
ax1.set_ylabel('Binding Energy (Hartree)')
ax1.set_yscale("log")

# Create secondary y-axis for error
ax2 = ax1.twinx()

# Plot error percentages on the secondary y-axis
error1 = ax2.plot(data1[:,0],
                 np.abs(neutral_atom["Binding Energy (a) (eV)"]/27.2114 + data1[:,1])/neutral_atom["Binding Energy (a) (eV)"] * 100,
                 "s", label="6-31G Error", color='red', alpha=0.7)
error3 = ax2.plot(data3[:,0],
                 np.abs(neutral_atom["Binding Energy (a) (eV)"]/27.2114 + data3[:,1])/neutral_atom["Binding Energy (a) (eV)"] * 100,
                 "s", label="STO-6G Error", color='blue', alpha=0.7)

# Set up secondary axis labels
ax2.set_ylabel('Error (%)')

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

# Set title and save
plt.title('Binding Energy of Atoms: Experimental vs. Predicted')
plt.tight_layout()
plt.savefig("figures/binding_energy_comparison.svg", bbox_inches='tight')
plt.show()
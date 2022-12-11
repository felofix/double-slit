import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette('pastel')
sns.set_style('darkgrid')

noslit = np.loadtxt('plotting/datafiles/probmatrix0.txt')
doslit = np.loadtxt('plotting/datafiles/probmatrix2.txt')
t = np.linspace(0, 0.008, len(doslit))
plt.xlabel('t')
plt.ylabel('|1 - $\sum_{i,j}u^{0*}_{ij}u^{0}_{ij}$|')

# Plotting noslit.
plt.plot(t, abs(1 - noslit))
plt.savefig("plotting/figures/noslit_probabilities.pdf")
plt.show()

# Plotting doubleslit.
plt.xlabel('t')
plt.ylabel('|1 - $\sum_{i,j}u^{0*}_{ij}u^{0}_{ij}$|')

plt.plot(t, abs(1 - doslit))
plt.savefig("plotting/figures/doubleslit_probabilities.pdf")
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
colors = sns.color_palette('pastel')
sns.set_style('darkgrid')

singlematrix = np.loadtxt("plotting/datafiles/matrix1.txt")
doublematrix = np.loadtxt("plotting/datafiles/matrix2.txt")
triplematrix = np.loadtxt("plotting/datafiles/matrix3.txt")
detects= singlematrix[-200:-1, int(0.8*199)]/np.sum(singlematrix[int(0.8*199), :])
detectd= doublematrix[-200:-1, int(0.8*199)]/np.sum(doublematrix[int(0.8*199), :])
detectt= triplematrix[-200:-1, int(0.8*199)]/np.sum(triplematrix[int(0.8*199), :])
yvalues = np.linspace(0, 1, len(detectd))

plt.xlabel('Y values')
plt.ylabel('Probability')
plt.plot(yvalues, detects)
plt.savefig("plotting/figures/detector_single.pdf")
plt.show()

plt.xlabel('Y values')
plt.ylabel('Probability')
plt.plot(yvalues, detectd)
plt.savefig("plotting/figures/detector_double.pdf")
plt.show()

plt.xlabel('Y values')
plt.ylabel('Probability')
plt.plot(yvalues, detectt)
plt.savefig("plotting/figures/detector_triple.pdf")
plt.show()

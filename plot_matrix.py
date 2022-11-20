import numpy as np
import matplotlib.pyplot as plt

datareal = np.loadtxt('test_matrix_real.txt')
dataimag = np.loadtxt('test_matrix_imag.txt')
plt.imshow(datareal)
plt.imshow(dataimag)
plt.show()

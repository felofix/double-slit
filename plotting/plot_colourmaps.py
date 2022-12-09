import numpy as np
import matplotlib.pyplot as plt

doublematrix = np.loadtxt("plotting/datafiles/matrix2.txt")
doublerealmatrix = np.loadtxt("plotting/datafiles/realmatrix2.txt")
doubleimagmatrix = np.loadtxt("plotting/datafiles/imagmatrix2.txt")

dt = 2.5e-5
t = np.linspace(dt, 0.002, int(dt/2.5e-5))
#extent = [x_min , x_max, y_min , y_max]

# All.
t0all = doublematrix[0:199, :]
t0001all = doublematrix[int(0.001/dt)*199:(int(0.001/dt)*199 + 199), :]
t002all = doublematrix[-200:-1, :]

# Real.
t0real = doublerealmatrix[0:199, :]
t0001real = doublerealmatrix[int(0.001/dt)*199:(int(0.001/dt)*199 + 199), :]
t002real = doublerealmatrix[-200:-1, :]

# Imag.
t0imag = doubleimagmatrix[0:199, :]
t0001imag = doubleimagmatrix[int(0.001/dt)*199:(int(0.001/dt)*199 + 199), :]
t002imag = doubleimagmatrix[-200:-1, :]

# All
plt.imshow(t0all,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/all_double_slit_at_t=0.pdf')
plt.show()

# All
plt.imshow(t0001all,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/all_double_slit_at_t=0.001.pdf')
plt.show()

# All
plt.imshow(t002all,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/all_double_slit_at_t=0.02.pdf')
plt.show()

# Real
plt.imshow(t0real,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/real_double_slit_at_t=0.pdf')
plt.show()

# Real
plt.imshow(t0001real,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/real_double_slit_at_t=0.001.pdf')
plt.show()

# Real
plt.imshow(t002real,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/real_double_slit_at_t=0.02.pdf')
plt.show()

# Imag
plt.imshow(t0imag,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/imag_double_slit_at_t=0.pdf')
plt.show()

# Imag
plt.imshow(t0001imag,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/imag_double_slit_at_t=0.001.pdf')
plt.show()

# Imag
plt.imshow(t002imag,extent=[0,1,0,1])
plt.colorbar()
plt.savefig('plotting/figures/imag_double_slit_at_t=0.02.pdf')
plt.show()



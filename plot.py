import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('lorenz_out.dat')

plt.plot(data[:,0], data[:,1])
plt.show()

plt.plot(data[:,1], data[:,2])
plt.show()

plt.plot(data[:,0], data[:,2])
plt.show()
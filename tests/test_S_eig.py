import numpy as np
from scipy.linalg import svdvals
import matplotlib.pyplot as plt

S_1p = np.loadtxt('S_test_1p.dat')
S = np.loadtxt('S_test.dat')

sigma = svdvals(S_1p)
print(np.max(sigma),np.min(sigma))
sigma = svdvals(S)
print(np.max(sigma),np.min(sigma))

fig,ax = plt.subplots()
ax.plot(sigma)
ax.set_yscale('log')

plt.show()

import numpy as np
import matplotlib.pyplot as plt

H = np.loadtxt('H_test.dat')
S = np.loadtxt('S_test.dat')

fig_H,ax_H = plt.subplots()
ax_H.spy(H)
ax_H.set_title('H')

fig_S,ax_S = plt.subplots()
ax_S.spy(S)
ax_S.set_title('S')


plt.show()

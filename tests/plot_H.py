import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

H = np.loadtxt('H_test.dat')
S = np.loadtxt('S_test.dat')

H_sp = sp.coo_array(H)
S_sp = sp.coo_array(S)

print(H.shape)
print(f'H nnz: {H_sp.nnz}')
print(f'H_sparsity: {1.0*H_sp.nnz/H.shape[0]**2}')
fig_H,ax_H = plt.subplots()
ax_H.spy(H)
ax_H.set_title('H')

print(f'S nnz: {S_sp.nnz}')
print(f'S_sparsity: {1.0*S_sp.nnz/S.shape[0]**2}')

fig_S,ax_S = plt.subplots()
ax_S.spy(S)
ax_S.set_title('S')

fig_H_mag,ax_H_mag = plt.subplots()
im = ax_H_mag.scatter(H_sp.col,H_sp.row,c=np.log10(H_sp.data),cmap='viridis')
ax_H_mag.set_title('abs(H)')
fig_H_mag.colorbar(im,ax=ax_H_mag)


plt.show()

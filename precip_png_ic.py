from PIL import Image
import numpy as np
im = Image.open("ngkern.png")
pix = im.load()
N, M = im.size
total_points = N * M
Lx, Ly = N, M
alpha=0.2

print("setting initial conditions")
s = np.hstack((np.zeros(total_points), np.full(total_points, alpha)))
for i in range(N):
    for j in range(M):
         s[i + (M - j - 1) * N]= 1.0 - pix[i,j][0]/255

icfile = open('ngkern.ic', 'wb')
np.save(icfile, N-1)
np.save(icfile, M-1)
np.save(icfile, Lx)
np.save(icfile, Ly)
np.save(icfile, s)

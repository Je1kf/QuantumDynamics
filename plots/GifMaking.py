import matplotlib.pyplot as plt
import numpy as np
import subprocess

data = np.loadtxt("../src/stationaryGaussian.txt")
sx = 10
nx = 160
nt = 1000

real = data[:, np.arange(0, round(nx*2 - 2), 2)]
imag = data[:, np.arange(1, round(nx*2 - 1), 2)]

psi = real + 1j*imag

x = np.linspace(-sx/2, sx/2, round(nx))
y = np.linspace(-sx/2, sx/2, round(nx))
x, y = np.meshgrid(x, y)


plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
plt.rcParams["figure.figsize"] = (10, 10)

aux = 0
for k in range(round(nt)):
    if k%10 == 0:
        i = k * nx;
        j = (k+1)*nx - 1;
        ax = plt.pcolor(x, y,
                        abs(psi[round(i):round(j), :])**2,
                        cmap='hot',
                        #vmax=abs(psi).max()**2,
                        vmin=0,
                        )
        plt.colorbar()
        plt.xlabel("x a.u.")
        plt.ylabel("z a.u.")
        plt.tight_layout()
        plt.savefig(f'images/{aux:000002}', dpi=300)
        plt.close()
        aux += 1

subprocess.run(['ffmpeg', '-framerate', '20', '-i', 'images/%02d.png', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', 'stationaryGaussian.mp4'])

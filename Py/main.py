import numpy as np
import matplotlib.pyplot as plt
import os

# Leer los datos
datos = np.loadtxt(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\normal.txt")

g1 = datos[:, 0]
g2 = datos[:, 1]

# Distribución normal teórica
x = np.linspace(-4, 4, 500)
normal = (1 / np.sqrt(2 * np.pi)) * np.exp(-x**2 / 2)

# Histograma de g1
plt.hist(g1, bins=50, density=True)
plt.plot(x, normal)

plt.xlabel("Valor")
plt.ylabel("Densidad")
plt.title("Distribución de g1")

plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/g1.png", dpi=300)
plt.show()
plt.close()

# Histograma de g2
plt.hist(g2, bins=50, density=True)
plt.plot(x, normal)

plt.xlabel("Valor")
plt.ylabel("Densidad")
plt.title("Distribución de g2")

plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/g2.png", dpi=300)
plt.show()
plt.close()

plt.scatter(g1, g2, s=0.1)
plt.xlabel("g1")
plt.ylabel("g2")
plt.title("g2 en función de g1")
plt.grid()
plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/g2_vs_g1.png", dpi=300)
plt.show()
plt.close()
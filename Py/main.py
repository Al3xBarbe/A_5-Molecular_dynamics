import numpy as np
import matplotlib.pyplot as plt
import os

# Plots de generador de numeros aleatoios
if False:

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

# Plots de las soluciones de los metodos (x(h) y p(h))
if False:
    # Nombre del archivo
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\euler_maruyama.txt"

    # Cargar los datos
    # columnas: t, x, p
    t, x, p = np.loadtxt(archivo, unpack=True)

    # Crear las gráficas
    fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    # x(t)
    ax[0].plot(t, x)
    ax[0].set_ylabel("x(t)")
    ax[0].set_title("Posición y momento")
    ax[0].grid(True)

    # v(t)
    ax[1].plot(t, p)
    ax[1].set_xlabel("t")
    ax[1].set_ylabel("p(t)")
    ax[1].grid(True)

    plt.tight_layout()
    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/EM_h0.01_nu10.png", dpi=300)
    plt.show()

    # Nombre del archivo
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\runge_kutta_2.txt"

    # Cargar los datos
    # columnas: t, x, p
    t, x, p = np.loadtxt(archivo, unpack=True)

    # Crear las gráficas
    fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    # x(t)
    ax[0].plot(t, x)
    ax[0].set_ylabel("x(t)")
    ax[0].set_title("Posición y momento")
    ax[0].grid(True)

    # v(t)
    ax[1].plot(t, p)
    ax[1].set_xlabel("t")
    ax[1].set_ylabel("p(t)")
    ax[1].grid(True)

    plt.tight_layout()
    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/RK2_h0.01_nu10.png", dpi=300)
    plt.show()

    # Nombre del archivo
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\verlet_exp_est.txt"

    # Cargar los datos
    # columnas: t, x, p
    t, x, p = np.loadtxt(archivo, unpack=True)

    # Crear las gráficas
    fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    # x(t)
    ax[0].plot(t, x)
    ax[0].set_ylabel("x(t)")
    ax[0].set_title("Posición y momento")
    ax[0].grid(True)

    # v(t)
    ax[1].plot(t, p)
    ax[1].set_xlabel("t")
    ax[1].set_ylabel("p(t)")
    ax[1].grid(True)

    plt.tight_layout()
    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/VerEst_h0.01_nu10.png", dpi=300)
    plt.show()

#Plots de las energias (Ki, V y E), equilibrio termico y equiparticion de energia
if True:

# Nombre del archivo
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\euler_maruyama.txt"
    
    # Cargar los datos
    # columnas: t, x, p
    t, x, p, Ki, V, E = np.loadtxt(archivo, unpack=True)

    beta_inv=1.0

    Ki_media = np.mean(Ki)

    plt.figure()

    plt.plot(t, Ki, label="Energía cinética Ki", linewidth=0.4)
    plt.plot(t, V, label="Energía potencial V", linewidth=0.3)
    plt.plot(t, E, label="Energía total E", linewidth=0.4)

    Ki_teorico=0.5*beta_inv

    plt.axhline(Ki_media, linestyle="--",
            label=f"<Ki> = {Ki_media:.3f}")
    plt.axhline(Ki_teorico, linestyle="-.",
            label=f"Ki teórico = {Ki_teorico:.3f}")
    
    plt.title("Energias")
    plt.xlabel("Tiempo")
    plt.ylabel("Energía")
    plt.legend()
    plt.grid()

    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/EM_Enrg_h0.001_nu0.01.png", dpi=300)
    plt.show()
    plt.close()

# Nombre del archivo
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\runge_kutta_2.txt"
    
    # Cargar los datos
    # columnas: t, x, p
    t, x, p, Ki, V, E = np.loadtxt(archivo, unpack=True)

    beta_inv=1.0

    Ki_media = np.mean(Ki)

    plt.figure()

    plt.plot(t, Ki, label="Energía cinética Ki", linewidth=0.4)
    plt.plot(t, V, label="Energía potencial V", linewidth=0.3)
    plt.plot(t, E, label="Energía total E", linewidth=0.4)

    Ki_teorico=0.5*beta_inv

    plt.axhline(Ki_media, linestyle="--",
            label=f"<Ki> = {Ki_media:.3f}")
    plt.axhline(Ki_teorico, linestyle="-.",
            label=f"Ki teórico = {Ki_teorico:.3f}")

    plt.title("Energias")
    plt.xlabel("Tiempo")
    plt.ylabel("Energía")
    plt.legend()
    plt.grid()

    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/RK2_Enrg_h0.001_nu0.01.png", dpi=300)
    plt.show()
    plt.close()

    # Nombre del archivo
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\verlet_exp_est.txt"
    
    # Cargar los datos
    # columnas: t, x, p
    t, x, p, Ki, V, E = np.loadtxt(archivo, unpack=True)

    beta_inv=1.0

    Ki_media = np.mean(Ki)

    plt.figure()

    plt.plot(t, Ki, label="Energía cinética Ki", linewidth=0.4)
    plt.plot(t, V, label="Energía potencial V", linewidth=0.3)
    plt.plot(t, E, label="Energía total E", linewidth=0.4)

    Ki_teorico=0.5*beta_inv

    plt.axhline(Ki_media, linestyle="--",
            label=f"<Ki> = {Ki_media:.3f}")
    plt.axhline(Ki_teorico, linestyle="-.",
            label=f"Ki teórico = {Ki_teorico:.3f}")

    plt.title("Energias")
    plt.xlabel("Tiempo")
    plt.ylabel("Energía")
    plt.legend()
    plt.grid()

    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/VerEst_Enrg_h0.001_nu0.01.png", dpi=300)
    plt.show()
    plt.close()

    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Plots de generador de numeros aleatoios
if False:

    # Leer los datos
    datos = np.loadtxt(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\normal_f90")

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
    plt.title("Distribución de g1 10k")

    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/g1_10k.png", dpi=300)
    plt.show()
    plt.close()

    # Histograma de g2
    plt.hist(g2, bins=50, density=True)
    plt.plot(x, normal)

    plt.xlabel("Valor")
    plt.ylabel("Densidad")
    plt.title("Distribución de g2 10k")

    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/g2_10k.png", dpi=300)
    plt.show()
    plt.close()

    plt.scatter(g1, g2, s=0.1)
    plt.xlabel("g1")
    plt.ylabel("g2")
    plt.title("g2 en función de g1 10k")
    plt.grid()
    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/g2_vs_g1_10k.png", dpi=300)
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
if False:

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

#Plots del distribución de x(t) y p(t)
if True:

        h = 0.0001
        eta = 10

        # ============================================================
        # GAUSSIANAS TEÓRICAS ESPERADAS
        # ============================================================

        x_teor = np.linspace(-4, 4, 500)
        p_teor = np.linspace(-4, 4, 500)

        gauss_x = (1 / np.sqrt(2 * np.pi)) * np.exp(-x_teor**2 / 2)
        gauss_p = (1 / np.sqrt(2 * np.pi)) * np.exp(-p_teor**2 / 2)

        # ============================================================
        # EULER-MARUYAMA
        # ============================================================

        archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\euler_maruyama.txt"

        t, x, p, Ki, V, E = np.loadtxt(archivo, unpack=True)

        fig, axs = plt.subplots(1, 2, figsize=(12, 5))

        axs[0].hist(x, bins=50, density=True, linewidth=0.5)
        axs[0].plot(x_teor, gauss_x, "--", linewidth=2,
                label="Gaussiana teórica")
        axs[0].set_title("Distribución de x(t)")
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("Densidad")
        axs[0].grid(alpha=0.3)
        axs[0].legend()

        axs[1].hist(p, bins=50, density=True, linewidth=0.5)
        axs[1].plot(p_teor, gauss_p, "--", linewidth=2,
                label="Gaussiana teórica")
        axs[1].set_title("Distribución de p(t)")
        axs[1].set_xlabel("p")
        axs[1].set_ylabel("Densidad")
        axs[1].grid(alpha=0.3)
        axs[1].legend()

        fig.suptitle(
        rf"Euler-Maruyama ($h={h}$, $\eta={eta}$)",
        fontsize=15
        )

        plt.tight_layout()

        plt.savefig(
        rf"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\EM_distr_h{h}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
        )

        plt.show()
        plt.close()


        # ============================================================
        # RUNGE-KUTTA
        # ============================================================

        archivo2 = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\runge_kutta_2.txt"

        print("Leyendo:", archivo2)

        t, x, p, Ki, V, E = np.loadtxt(
        archivo2,
        unpack=True
        )

        fig, axs = plt.subplots(1, 2, figsize=(12, 5))

        axs[0].hist(x, bins=50, density=True, linewidth=0.5)
        axs[0].plot(x_teor, gauss_x, "--", linewidth=2,
                label="Gaussiana teórica")
        axs[0].set_title("Distribución de x(t)")
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("Densidad")
        axs[0].grid(alpha=0.3)
        axs[0].legend()

        axs[1].hist(p, bins=50, density=True, linewidth=0.5)
        axs[1].plot(p_teor, gauss_p, "--", linewidth=2,
                label="Gaussiana teórica")
        axs[1].set_title("Distribución de p(t)")
        axs[1].set_xlabel("p")
        axs[1].set_ylabel("Densidad")
        axs[1].grid(alpha=0.3)
        axs[1].legend()

        fig.suptitle(
        rf"Runge-Kutta 2 ($h={h}$, $\eta={eta}$)",
        fontsize=15
        )

        plt.tight_layout()

        plt.savefig(
        rf"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\RK2_distr_h{h}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
        )

        plt.show()
        plt.close()


        # ============================================================
        # VELOCITY VERLET
        # ============================================================

        archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\verlet_exp_est.txt"

        t, x, p, Ki, V, E = np.loadtxt(
        archivo,
        unpack=True
        )

        fig, axs = plt.subplots(1, 2, figsize=(12, 5))

        axs[0].hist(x, bins=50, density=True, linewidth=0.5)
        axs[0].plot(x_teor, gauss_x, "--", linewidth=2,
                label="Gaussiana teórica")
        axs[0].set_title("Distribución de x(t)")
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("Densidad")
        axs[0].grid(alpha=0.3)
        axs[0].legend()

        axs[1].hist(p, bins=50, density=True, linewidth=0.5)
        axs[1].plot(p_teor, gauss_p, "--", linewidth=2,
                label="Gaussiana teórica")
        axs[1].set_title("Distribución de p(t)")
        axs[1].set_xlabel("p")
        axs[1].set_ylabel("Densidad")
        axs[1].grid(alpha=0.3)
        axs[1].legend()

        fig.suptitle(
        rf"Verlet Explicito ($h={h}$, $\eta={eta}$)",
        fontsize=15
        )

        plt.tight_layout()

        plt.savefig(
        rf"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\VEE_distr_h{h}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
        )

        plt.show()
        plt.close()
        
#Plot termalización de la Energia
if False:

        eta = 10
        h = 0.0001

        archivo_EM = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\termalizacion_EM.txt"

        archivo_RK = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\term_RK.txt"

        archivo_VE = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\term_VE.txt"


        t_EM, x_EM, p_EM, Ki_EM, V_EM = np.loadtxt(
        archivo_EM, unpack=True
        )

        t_RK, x_RK, p_RK, Ki_RK, V_RK = np.loadtxt(
        archivo_RK, unpack=True
        )

        t_VE, x_VE, p_VE, Ki_VE, V_VE = np.loadtxt(
        archivo_VE, unpack=True
        )

        Ki_final_EM = Ki_EM[-1]
        Ki_final_RK = Ki_RK[-1]
        Ki_final_VE = Ki_VE[-1]

        V_final_EM = V_EM[-1]
        V_final_RK = V_RK[-1]
        V_final_VE = V_VE[-1]

        plt.figure(figsize=(8, 5))

        plt.plot(
        t_EM, Ki_EM,
        linewidth=0.6,
        label=r"$\langle K\rangle$"
        )

        plt.plot(
        t_EM, V_EM,
        linewidth=0.6,
        label=r"$\langle V\rangle$"
        )

        plt.axhline(
        0.5,
        linestyle="--",
        linewidth=0.8,
        label=r"$k_BT/2=0.5$"
        )

        plt.title("Euler-Maruyama")
        plt.xlabel("Tiempo")
        plt.ylabel("Energía media")

        plt.legend(
        title=fr"$\eta={eta}$, $h={h}$" + "\n"
                fr"$\langle V\rangle_f={V_final_EM:.4f}$" + "\n"
                fr"$\langle Ki\rangle_f={Ki_final_EM:.4f}$"
        )

        plt.grid(alpha=0.3)

        plt.tight_layout()

        plt.savefig(
                fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\termalizacion_EM_h{h}_eta{eta}.png",
                dpi=300,
                bbox_inches="tight"
                )

        plt.show()
        plt.close()


        # ============================================================
        # RUNGE-KUTTA 2
        # ============================================================

        plt.figure(figsize=(8, 5))

        plt.plot(
        t_RK, Ki_RK,
        linewidth=0.6,
        label=r"$\langle K\rangle$"
        )

        plt.plot(
        t_RK, V_RK,
        linewidth=0.6,
        label=r"$\langle V\rangle$"
        )

        plt.axhline(
        0.5,
        linestyle="--",
        linewidth=0.8,
        label=r"$k_BT/2=0.5$"
        )

        plt.title("Runge-Kutta 2")
        plt.xlabel("Tiempo")
        plt.ylabel("Energía media")

        plt.legend(
        title=fr"$\eta={eta}$, $h={h}$" + "\n"
                fr"$\langle V\rangle_f={V_final_RK:.4f}$" + "\n"
                fr"$\langle Ki\rangle_f={Ki_final_RK:.4f}$"
        )

        plt.grid(alpha=0.3)

        plt.tight_layout()

        plt.savefig(
                fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\termalizacion_RK_h{h}_eta{eta}.png",
                dpi=300,
                bbox_inches="tight"
        )

        plt.show()
        plt.close()


        # ============================================================
        # VERLET
        # ============================================================

        plt.figure(figsize=(8, 5))

        plt.plot(
        t_VE, Ki_VE,
        linewidth=0.6,
        label=r"$\langle K\rangle$"
        )

        plt.plot(
        t_VE, V_VE,
        linewidth=0.6,
        label=r"$\langle V\rangle$"
        )

        plt.axhline(
        0.5,
        linestyle="--",
        linewidth=0.8,
        label=r"$k_BT/2=0.5$"
        )

        plt.title("Verlet")
        plt.xlabel("Tiempo")
        plt.ylabel("Energía media")

        plt.legend(
        title=fr"$\eta={eta}$, $h={h}$" + "\n"
                fr"$\langle V\rangle_f={V_final_VE:.4f}$" + "\n"
                fr"$\langle Ki\rangle_f={Ki_final_VE:.4f}$"
        )

        plt.grid(alpha=0.3)

        plt.tight_layout()

        plt.savefig(
                fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\termalizacion_VE_h{h}_eta{eta}.png",
                dpi=300,
                bbox_inches="tight"
                )

        plt.show()
        plt.close()

#Plot correlación numeros aleatorios
if False:
        # Leer el archivo
        datos = np.loadtxt(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\randcongr")

        # Separar las columnas
        x = datos[:, 0]
        y = datos[:, 1]

        # Representar y en función de x
        plt.scatter(x, y, color="black", s=0.05)
        plt.xlabel("x1")
        plt.ylabel("x2")
        plt.title("Generador_congr 100K")
        plt.grid(True)

        plt.savefig(
                r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\rand_congr_100k.png",
                dpi=300,
                bbox_inches="tight"
                )
        plt.show()
        plt.close()
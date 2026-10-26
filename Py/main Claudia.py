import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.stats import linregress


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
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\verlet_doblepozo.txt"

    # Cargar los datos
    # columnas: t, x, p
    t, x, p = np.loadtxt(archivo, unpack=True)

    # Crear las gráficas
    fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    # x(t)
    ax[0].plot(t, x)
    ax[0].set_ylabel("x(t)")
    ax[0].set_title("Posición y momento (h0.001, eta1.0)")
    ax[0].grid(True)

    # v(t)
    ax[1].plot(t, p)
    ax[1].set_xlabel("t")
    ax[1].set_ylabel("p(t)")
    ax[1].grid(True)

    plt.tight_layout()
    plt.savefig(r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/VerDP_h0.001_nu1.0.png", dpi=300)
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

#distribuciones de x(t) y p(t) doble pozo
if False:

    h = 0.001
    eta = 1.0
    beta_inv = 0.2
    m = 1.0
    sigma =  np.sqrt(beta_inv*m)

    p_teor = np.linspace(-4, 4, 500)

    gauss_p = (1 / (np.sqrt(2 * np.pi)*sigma)) * np.exp(-p_teor**2 / (2*sigma**2))

    # ============================================================
    # VERLET EXPLÍCITO
    # ============================================================

    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\verlet_doblepozo.txt"

    t, x, p, Ki, V, est = np.loadtxt(
        archivo,
        unpack=True
    )

    # ============================================================
    # FIGURA: DISTRIBUCIONES DE x(t) Y p(t)
    # ============================================================

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    # Distribución de x(t)
    axs[0].hist(
        x,
        bins=50,
        density=True,
        linewidth=0.5
    )


    axs[0].set_title("Distribución de x(t)")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("Densidad")
    axs[0].grid(alpha=0.3)
    axs[0].legend()

    # Distribución de p(t)
    axs[1].hist(
        p,
        bins=50,
        density=True,
        linewidth=0.5
    )

    axs[1].plot(
        p_teor,
        gauss_p,
        "--",
        linewidth=2,
        label="Gaussiana teórica"
    )

    axs[1].set_title("Distribución de p(t)")
    axs[1].set_xlabel("p")
    axs[1].set_ylabel("Densidad")
    axs[1].grid(alpha=0.3)
    axs[1].legend()

    # Título general
    fig.suptitle(
        rf"Verlet Explícito ($h={h}$, $\eta={eta}$)",
        fontsize=15
    )

    plt.tight_layout()

    # ============================================================
    # GUARDAR FIGURA
    # ============================================================

    plt.savefig(
        rf"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\VDP_distr_h{h}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close()

#Plots del distribución de x(t) y p(t)
if False:

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

        eta = 1.0
        h = 0.001

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

#Plot equipartición de la Energia O2
if False:
        eta = 1.0
        h = 0.001
        archivo_VE = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\verlet_doblepozo.txt"
        t_VE, x_VE, p_VE, Ki_VE, V_VE = np.loadtxt(
                archivo_VE, unpack=True
                )
        Ki_final_VE = Ki_VE[-1]
        V_final_VE = V_VE[-1]
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
                fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\Equip_VDB_h{h}_eta{eta}.png",
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

#Plots objetivo 2 (h0.001, verlet explicito estocastico)
if False:

    A=2.0
    eta=3.5
    beta_inv = 0.2
    m = 1.0

    sigma =  np.sqrt(beta_inv*m)
    p_teor = np.linspace(-4, 4, 500)
    gauss_p = (1 / (np.sqrt(2 * np.pi)*sigma)) * np.exp(-p_teor**2 / (2*sigma**2))
    
    archivo = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\verlet_doblepozo.txt"
    pozo1 = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\estancia_pozo1.txt"
    pozo2 = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\estancia_pozo2.txt"

    # Cargar los datos
    t, x, p, Ki_m, V_m, Fraccion_m = np.loadtxt(archivo, unpack=True)
    Ki_final_VE = Ki_m[-1]
    V_final_VE = V_m[-1]

    est1 = np.loadtxt(pozo1,ndmin=1)
    est2 = np.loadtxt(pozo2,ndmin=1)

    # Crear las gráficas
    fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    # x(t)
    ax[0].plot(t, x)
    ax[0].set_ylabel("x(t)")
    ax[0].set_title(f"Posición y momento (A{A}, eta{eta})")
    ax[0].grid(True)

    # v(t)
    ax[1].plot(t, p)
    ax[1].set_xlabel("t")
    ax[1].set_ylabel("p(t)")
    ax[1].grid(True)

    plt.tight_layout()
    plt.savefig(rf"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots/VerDP_A{A}_eta{eta}.png", dpi=300)
    plt.show()
    plt.close()

    plt.figure(figsize=(8, 5))
            
    plt.plot(
        t, Ki_m,
        linewidth=0.6,
        label=r"$\langle K\rangle$"
        )
            
    plt.plot(
        t, V_m,
        linewidth=0.6,
        label=r"$\langle V\rangle$"
        )
            
    plt.axhline(
        0.5*beta_inv,
        linestyle="--",
        linewidth=0.8,
        label=fr"$k_BT/2={0.5*beta_inv:.1f}$"
        )
    
    plt.title("Verlet")
    plt.xlabel("Tiempo")
    plt.ylabel("Energía media")
    
    plt.legend(
        title=fr"$\eta={eta}$, $A={A}$" + "\n"
        fr"$\langle V\rangle_f={V_final_VE:.4f}$" + "\n"
        fr"$\langle Ki\rangle_f={Ki_final_VE:.4f}$"
        )
    
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    
    plt.savefig(
        fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\Equip_VDB_A{A}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
        )
    
    plt.show()
    plt.close()

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    # Distribución de x(t)
    axs[0].hist(
        x,
        bins=50,
        density=True,
        linewidth=0.5
    )


    axs[0].set_title("Distribución de x(t)")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("Densidad")
    axs[0].grid(alpha=0.3)

    # Distribución de p(t)
    axs[1].hist(
        p,
        bins=50,
        density=True,
        linewidth=0.5
    )

    axs[1].plot(
        p_teor,
        gauss_p,
        "--",
        linewidth=2,
        label="Gaussiana teórica"
    )

    axs[1].set_title("Distribución de p(t)")
    axs[1].set_xlabel("p")
    axs[1].set_ylabel("Densidad")
    axs[1].grid(alpha=0.3)
    axs[1].legend()

    # Título general
    fig.suptitle(
        rf"Verlet Explícito ($A={A}$, $\eta={eta}$)",
        fontsize=15
    )

    plt.tight_layout()

    # ============================================================
    # GUARDAR FIGURA
    # ============================================================

    plt.savefig(
        rf"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\VDP_distr_A{A}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close()

    plt.figure(figsize=(8, 5))
            
    plt.plot(
        t, Fraccion_m,
        linewidth=0.6,
        )
    
    plt.title(f"Fraccion de tiempo en pozo 2 (A={A}, eta={eta})")
    plt.xlabel("Tiempo")
    plt.ylabel("Fracción de tiempo")
    
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    
    plt.savefig(
        fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\FraccOcup_VDB_A{A}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
        )
    
    plt.show()
    plt.close()

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    # Distribución de x(t)
    axs[0].hist(
        est1,
        bins=50,
        density=True,
        linewidth=0.5
    )


    axs[0].set_title("Distribución pozo 1")
    axs[0].set_xlabel("tiempo de estancia")
    axs[0].set_ylabel("Densidad")
    axs[0].grid(alpha=0.3)

    # Distribución de p(t)
    axs[1].hist(
        est2,
        bins=50,
        density=True,
        linewidth=0.5
    )


    axs[1].set_title("Distribución pozo 2")
    axs[1].set_xlabel("tiempo de estancia")
    axs[1].set_ylabel("Densidad")
    axs[1].grid(alpha=0.3)

    # Título general
    fig.suptitle(
        rf"Tiempos de estancia ($A={A}$, $\eta={eta}$)",
        fontsize=15
    )

    plt.tight_layout()

    # ============================================================
    # GUARDAR FIGURA
    # ============================================================

    plt.savefig(
        rf"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\t_estancia_A{A}_eta{eta}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close()


# Gráficas de termalización: comparación de algoritmos 
if False:

        # Leer tabla
        datos = np.loadtxt("Fortran/Termalizacion_valores.txt", dtype=str)

        # Separar columnas
        algoritmo = datos[:, 0]
        eta = datos[:, 1].astype(float)
        h = datos[:, 2].astype(float)
        Ki = datos[:, 3].astype(float)
        V = datos[:, 4].astype(float)


        def plot_termalizacion(valor_eta, tipo_energia):

                # Seleccionar solo las filas del eta que queremos
                mask_eta = (eta == valor_eta)

                alg_eta = algoritmo[mask_eta]
                h_eta = h[mask_eta]
                Ki_eta = Ki[mask_eta]
                V_eta = V[mask_eta]

                # Separar por algoritmo
                mask_euler = (alg_eta == "Euler")
                mask_rk = (alg_eta == "RK")
                mask_ve = (alg_eta == "VE")

                h_euler = h_eta[mask_euler]
                h_rk = h_eta[mask_rk]
                h_ve = h_eta[mask_ve]

                # Elegir qué energía queremos representar
                if tipo_energia == "cinetica":
                        E_euler = Ki_eta[mask_euler]
                        E_rk = Ki_eta[mask_rk]
                        E_ve = Ki_eta[mask_ve]

                        ylabel = "Energía cinética media"

                elif tipo_energia == "potencial":
                        E_euler = V_eta[mask_euler]
                        E_rk = V_eta[mask_rk]
                        E_ve = V_eta[mask_ve]

                        ylabel = "Energía potencial media"

                else:
                        print("Tipo de energía no válido")
                        return

                # Crear gráfica
                plt.figure()

                plt.plot(h_euler, E_euler, marker="o", label="Euler-Maruyama")
                plt.plot(h_rk, E_rk, marker="o", label="Runge-Kutta 2")
                plt.plot(h_ve, E_ve, marker="o", label="Verlet explícito")

                # Valor teórico esperado por equipartición
                plt.axhline(0.5, linestyle="--", label="Valor teórico = 0.5")

                plt.xscale("log")

                plt.xlabel("h")
                plt.ylabel(ylabel)
                plt.title(f"Termalización: eta = {valor_eta}")

                plt.legend()
                plt.grid()

                plt.show()

#%%
if False:

    A=2
    eta=0.7
    beta_inv = 0.2
    m = 1.0
    
    archivo = "/Users/juansazperez/Documents/A_5-Molecular_dynamics/verlet_doblepozo.txt"

    # Cargar los datos
    t, x, p, Ki_m, V_m, Fraccion_m = np.loadtxt(archivo, unpack=True)
    Ki_final_VE = Ki_m[-1]
    V_final_VE = V_m[-1]

    plt.figure(figsize=(8, 5))
            
    plt.plot(
        t, Ki_m,
        linewidth=0.6,
        label=r"$\langle K\rangle$"
        )
            
    plt.plot(
        t, V_m,
        linewidth=0.6,
        label=r"$\langle V\rangle$"
        )
            
    plt.axhline(
        0.5*beta_inv,
        linestyle="--",
        linewidth=0.8,
        label=fr"$k_BT/2={0.5*beta_inv:.1f}$"
        )
    
    plt.title("Verlet")
    plt.xlabel("Tiempo")
    plt.ylabel("Energía media")
    
    plt.legend(
        title=fr"$\eta={eta}$, $A={A}$" + "\n"
        fr"$\langle V\rangle_f={V_final_VE:.4f}$" + "\n"
        fr"$\langle Ki\rangle_f={Ki_final_VE:.4f}$"
        )
    
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    
    plt.savefig(
    fr"/Users/juansazperez/Documents/A_5-Molecular_dynamics/Py/plots/OBJETIVO_2/Equip_term_DP/amortiguado/A2.0/Equip_VDB_A{A}_eta{eta}.png",
    dpi=300,
    bbox_inches="tight"
)
    
    plt.show()
    plt.close()

#fraccion de ocupacion en funcion de la fuerza
if False:
    ruta = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\OBJETIVO_2\Fracc_ocup(F)\fracc_p2(F)_3.txt"

    datos = np.loadtxt(ruta, comments="#", ndmin=2)

    x = datos[:, 0]
    y = datos[:, 1]

    plt.figure(figsize=(8, 5))
    plt.plot(x, y, "o-", markersize=4)
    plt.xlabel("Fuerza")
    plt.ylabel("Fracción de ocupación del pozo 2")
    plt.title("Dependencia de la fracción de ocupación del pozo 2 con la fuerza")
    plt.grid(True)

    # Leyenda solo con texto (sin símbolo)
    plt.plot([], [], " ", label="Tiempo = 500000\nA = 0.5\n$\\eta$ = 1.0")
    plt.legend(handlelength=0, handletextpad=0)

    plt.tight_layout()
    plt.show()




#Plots Claudia Objetivo 2
if False:
    
    # OBJETIVO 2 - COCIENTE DE OCUPACIONES EN FUNCIÓN DE LA FUERZA

    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path

    # Carpeta donde están los archivos de datos
    carpeta = Path("/Users/juansazperez/Documents/A_5-Molecular_dynamics/Py/plots/OBJETIVO_2/Fracc_ocup(F)")

    # Archivos y parámetros de las tres simulaciones
    simulaciones = [
        ("fracc_p2(F)_1.txt", 0.5, 0.4),
        ("fracc_p2(F)_2.txt", 0.5, 1.0),
        ("fracc_p2(F)_3.txt", 1.0, 0.4),
    ]

    for nombre_archivo, A, eta in simulaciones:

        # Leer datos: primera columna F, segunda P+
        datos = np.loadtxt(carpeta / nombre_archivo, comments="#")

        F = datos[:, 0]
        P = datos[:, 1]

        # Calcular cociente de ocupaciones
        # Evitamos dividir entre cero cuando P = 1
        R = np.divide(
            P,
            1.0 - P,
            out=np.full_like(P, np.nan),
            where=(P < 1.0)
        )

        # Crear una figura independiente
        plt.figure(figsize=(8, 5))

        plt.plot(F, R, "o-", markersize=4)

        plt.xlabel("Fuerza F")
        plt.yscale("log")
        plt.ylabel(r"Cociente de ocupaciones $n_+/n_-$")
        plt.title(
        "Objetivo 2: Cociente de ocupaciones frente a la fuerza\n"
        "Representación en escala semilogarítmica"
        )

        plt.grid(alpha=0.3)
        plt.tight_layout()

        # Guardar en la misma carpeta
        salida = carpeta / f"Cociente_ocup_A{A}_eta{eta}.png"

        plt.savefig(salida, dpi=300, bbox_inches="tight")

        plt.show()
        plt.close()

    # ============================================================
    # OBJETIVO 2 - TIEMPO MEDIO DE ESTANCIA FRENTE A ETA, calculo de los tiempos medios de estancia
    # ============================================================

    import numpy as np
    from pathlib import Path


    carpeta = Path("/Users/juansazperez/Documents/A_5-Molecular_dynamics/Py/plots/OBJETIVO_2/Tiempos_medios_eta")
    print("¿Existe la carpeta?", carpeta.exists())

    for archivo in carpeta.glob("*.txt"):
        print(archivo.name)
    A = 0.5
    etas = [0.1, 0.4, 0.7, 1.0]

    for eta in etas:

        nombre1 = f"estancia_pozo1_A{A:.2f}_eta{eta:.2f}.txt"
        nombre2 = f"estancia_pozo2_A{A:.2f}_eta{eta:.2f}.txt"

        est1 = np.loadtxt(carpeta / nombre1, ndmin=1)
        est2 = np.loadtxt(carpeta / nombre2, ndmin=1)

        print(f"\nA = {A}, eta = {eta}")
        print(f"Pozo 1: {len(est1)} estancias, media = {np.mean(est1):.4f}")
        print(f"Pozo 2: {len(est2)} estancias, media = {np.mean(est2):.4f}")

    # ============================================================
    # OBJETIVO 2 - TIEMPO MEDIO DE ESTANCIA FRENTE A ETA, representacion grafica de los tiempos medios de estancia
    # ============================================================

    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path

    # Carpeta con los archivos de las simulaciones
    carpeta = Path(
        "/Users/juansazperez/Documents/A_5-Molecular_dynamics/"
        "Py/plots/OBJETIVO_2/Tiempos_medios_eta/A2.0"
    )

    # Parametros
    A = 2.0
    etas = [0.1, 0.4, 0.7, 1.0]

    medias1 = []
    medias2 = []

    for eta in etas:

        archivo1 = carpeta / f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
        archivo2 = carpeta / f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"

        est1 = np.loadtxt(archivo1, ndmin=1)
        est2 = np.loadtxt(archivo2, ndmin=1)

        medias1.append(np.mean(est1))
        medias2.append(np.mean(est2))

    # Representacion grafica
    plt.figure(figsize=(8, 5))

    plt.plot(etas, medias1, "o-", label="Pozo 1")
    plt.plot(etas, medias2, "s-", label="Pozo 2")

    plt.xlabel(r"Coeficiente de amortiguamiento $\eta$")
    plt.ylabel(r"Tiempo medio de estancia $\langle \tau \rangle$")

    plt.title(
        "Objetivo 2: Tiempo medio de estancia frente al amortiguamiento"
    )

    plt.legend(title=rf"$A={A}$, $k_BT=0.2$")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Guardar grafica
    salida = carpeta / f"Tiempo_medio_estancia_A{A}.png"

    plt.savefig(salida, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()

    # ============================================================
    # OBJETIVO 2 - ANALISIS ESTADISTICO DE TIEMPOS DE ESTANCIA
    # ============================================================

    import numpy as np
    import pandas as pd
    from pathlib import Path

    # Carpeta principal de nuestros datos
    base = Path(
        "/Users/juansazperez/Documents/A_5-Molecular_dynamics/"
        "Py/plots/OBJETIVO_2/Tiempos_medios_eta"
    )

    # Parametros simulados
    valores_A = [0.5, 1.0, 2.0]
    etas = [0.1, 0.4, 0.7, 1.0]

    resultados = []

    for A in valores_A:

        carpeta = base / f"A{A:.1f}"

        for eta in etas:

            archivo1 = carpeta / f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
            archivo2 = carpeta / f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"

            # Leer tiempos de estancia
            est1 = np.loadtxt(archivo1, ndmin=1)
            est2 = np.loadtxt(archivo2, ndmin=1)

            # Numero de estancias
            N1 = len(est1)
            N2 = len(est2)

            # Tiempos medios
            media1 = np.mean(est1) if N1 > 0 else np.nan
            media2 = np.mean(est2) if N2 > 0 else np.nan

            # Ambos pozos
            todas = np.concatenate((est1, est2))
            Ntotal = len(todas)
            media_total = np.mean(todas) if Ntotal > 0 else np.nan

            resultados.append({
                "A": A,
                "eta": eta,
                "N_pozo1": N1,
                "N_pozo2": N2,
                "N_total": Ntotal,
                "Media_pozo1": media1,
                "Media_pozo2": media2,
                "Media_total": media_total
            })

    # Crear tabla
    tabla = pd.DataFrame(resultados)

    # Mostrar resultados
    print(tabla.round(3).to_string(index=False))

    # Guardar tabla en CSV
    salida = base / "Resumen_tiempos_estancia.csv"
    tabla.to_csv(salida, index=False)

    print("\nTabla guardada en:", salida)

    # ============================================================
    # OBJETIVO 2 - AJUSTES EXPONENCIALES DE TIEMPOS DE ESTANCIA
    # A = 1.0, distintos valores de eta
    # ============================================================

    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path

    base = Path.home() / "Documents/A_5-Molecular_dynamics"

    A = 1.0

    # Valores de eta y cortes elegidos
    cortes_por_eta = {
        0.1: 50.0,
        0.4: 50.0,
        0.7: 50.0,
        1.0: 50.0
    }


    carpeta_datos = (
        base / "Py/plots/OBJETIVO_2/Tiempos_medios_eta/A1.0"
    )

    carpeta_salida = (
        base / "Py/plots/OBJETIVO_2/Distr_testancia/amortiguado/A1.0"
    )
    carpeta_salida.mkdir(parents=True, exist_ok=True)

    for eta, tau_min in cortes_por_eta.items():

        archivo1 = carpeta_datos / f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
        archivo2 = carpeta_datos / f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"

        est1 = np.loadtxt(archivo1, ndmin=1)
        est2 = np.loadtxt(archivo2, ndmin=1)

        fig, axs = plt.subplots(1, 2, figsize=(13, 5))

        for ax, tiempos, nombre in zip(
            axs, [est1, est2], ["Pozo 1", "Pozo 2"]
        ):

            # Histograma normalizado
            densidad, bordes = np.histogram(
                tiempos, bins=50, density=True
            )

            ax.stairs(
                densidad, bordes,
                fill=True,
                alpha=0.4,
                label="Distribución simulada"
            )

            # Seleccion de la cola
            cola = tiempos[tiempos >= tau_min]

            N_total = len(tiempos)
            N_cola = len(cola)

            if N_cola < 2:
                print(f"Datos insuficientes: A={A}, eta={eta}, {nombre}")
                continue

            # Estimacion de maxima verosimilitud
            k = 1.0 / np.mean(cola - tau_min)
            fraccion_cola = N_cola / N_total

            # Curva exponencial ajustada
            x_ajuste = np.linspace(
                tau_min, np.max(tiempos), 500
            )

            y_ajuste = (
                fraccion_cola * k
                * np.exp(-k * (x_ajuste - tau_min))
            )

            # Expresar como C exp(-k*tau)
            C = fraccion_cola * k * np.exp(k * tau_min)

            exponente = int(np.floor(np.log10(C)))
            coeficiente = C / 10**exponente

            ecuacion = (
                rf"$p(\tau)={coeficiente:.3f}\times10^{{{exponente}}}"
                rf"e^{{-{k:.4f}\tau}}$"
            )

            ax.plot(
                x_ajuste, y_ajuste, "--",
                linewidth=2,
                label=ecuacion
            )

            # Indicar corte
            ax.axvline(
                tau_min,
                color="gray",
                linestyle=":",
                linewidth=1.2,
                label=rf"$\tau_{{\min}}={tau_min:g}$"
            )

            # Escala semilogaritmica
            ax.set_yscale("log")

            ax.set_xlabel(r"Tiempo de estancia $\tau$")
            ax.set_ylabel(r"Densidad de probabilidad $p(\tau)$")
            ax.set_title(nombre)
            ax.grid(alpha=0.3, which="both")
            ax.legend(fontsize=9)

            print(
                f"A={A}, eta={eta}, {nombre}: "
                f"N_cola={N_cola}, k={k:.5f}, "
                f"1/k={1/k:.2f}"
            )

        fig.suptitle(
            "Distribución de tiempos de estancia y ajuste de la cola exponencial\n"
            rf"Escala semilogarítmica ($A={A}$, $\eta={eta}$)",
            fontsize=14
        )

        plt.tight_layout()

        salida = carpeta_salida / f"Ajuste_cola_A{A}_eta{eta}.png"

        plt.savefig(salida, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close()

    # ============================================================
    # OBJETIVO 2 - ESTABILIDAD DEL AJUSTE EXPONENCIAL
    # A = 0.5, distintos valores de eta
    # ============================================================

    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    import pandas as pd

    resultados = []

    base = Path.home() / "Documents/A_5-Molecular_dynamics"

    # Parámetros
    A = 0.5
    etas = [0.01, 0.04]
    cortes = [0, 10, 20, 30, 50, 80]


    # Carpeta de datos
    carpeta = base / "Py/plots/OBJETIVO_2/Tiempos_medios_eta/A0.5"

    # Figura: una fila por eta y dos columnas por pozo
    fig, axs = plt.subplots(2, 2, figsize=(12, 9), sharex=True)

    for i, eta in enumerate(etas):

        archivo1 = carpeta / f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
        archivo2 = carpeta / f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"

        est1 = np.loadtxt(archivo1, ndmin=1)
        est2 = np.loadtxt(archivo2, ndmin=1)

        for j, tiempos in enumerate([est1, est2]):

            tasas = []
            cantidades = []

            print(f"\nA={A}, eta={eta}, pozo {j+1}")
            print("Corte     N datos         k         1/k")

            for corte in cortes:

                cola = tiempos[tiempos >= corte]
                N = len(cola)

                if N >= 2 and np.mean(cola - corte) > 0:
                    k = 1.0 / np.mean(cola - corte)
                    tau_car = 1.0 / k
                else:
                    k = np.nan
                    tau_car = np.nan

                tasas.append(k)
                cantidades.append(N)
                resultados.append({
                    "A": A,
                    "eta": eta,
                    "Pozo": j + 1,
                    "tau_min": corte,
                    "N_cola": N,
                    "k": k,
                    "Tiempo_caracteristico": tau_car
                })

                print(f"{corte:5.0f} {N:10d} {k:11.5f} {tau_car:11.3f}")

            ax = axs[i, j]

            ax.plot(cortes, tasas, "o-")

            ax.set_title(rf"Pozo {j+1}, $\eta={eta}$")
            ax.set_xlabel(r"Tiempo mínimo $\tau_{\min}$")
            ax.set_ylabel(r"Tasa de decaimiento $k$")

            # Mostrar todos los valores del eje X
            ax.set_xticks(cortes)
            ax.tick_params(axis="x", labelbottom=True)

            ax.grid(alpha=0.3)

            for corte, k, N in zip(cortes, tasas, cantidades):
                if np.isfinite(k):
                    ax.annotate(
                        f"N={N}",
                        (corte, k),
                        xytext=(0, 8),
                        textcoords="offset points",
                        ha="center",
                        fontsize=8
                    )

    fig.suptitle(
        rf"Estabilidad de la tasa exponencial ($A={A}$)",
        fontsize=15
    )

    # Carpeta de salida
    salida = base / "Py/plots/OBJETIVO_2/Distr_testancia/amortiguado/A0.5/Analisis_cola"
    salida.mkdir(parents=True, exist_ok=True)

    plt.tight_layout()

    plt.savefig(
        salida / "Estabilidad_k_A0.5.png",
        dpi=300,
        bbox_inches="tight"
    )

    tabla = pd.DataFrame(resultados)

    tabla.to_csv(
        salida / "Estabilidad_k_A0.5.csv",
        index=False
    )
    plt.show()
    plt.close()

    # ============================================================
    # OBJETIVO 2 - ESTABILIDAD DEL AJUSTE EXPONENCIAL
    # A = 0.5, distintos valores de eta
    # ============================================================

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pathlib import Path

    # Directorio principal
    base = Path.home() / "Documents/A_5-Molecular_dynamics"

    # Parametros
    A = 0.5
    etas = [0.01, 0.04, 3.5, 5.0]
    cortes = [0, 20, 50, 100, 150, 200, 300]

    # Carpetas
    carpeta_datos = (
        base / "Py/plots/OBJETIVO_2/Tiempos_medios_eta/A0.5"
    )

    carpeta_salida = (
        base / "Py/plots/OBJETIVO_2/Distr_testancia/subamortiguado/A0.5/Analisis_cola"
    )

    carpeta_salida.mkdir(parents=True, exist_ok=True)

    # Lista donde guardaremos los resultados
    resultados = []

    # Una fila por eta y dos columnas por pozo
    fig, axs = plt.subplots(4, 2, figsize=(13, 15), sharex=True)

    for i, eta in enumerate(etas):

        archivo1 = carpeta_datos / f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
        archivo2 = carpeta_datos / f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"

        est1 = np.loadtxt(archivo1, ndmin=1)
        est2 = np.loadtxt(archivo2, ndmin=1)

        for j, tiempos in enumerate([est1, est2]):

            tasas = []
            cantidades = []

            print(f"\nA={A}, eta={eta}, pozo {j+1}")
            print("Corte     N datos        k         1/k")

            for corte in cortes:

                cola = tiempos[tiempos >= corte]
                N = len(cola)

                if N >= 2 and np.mean(cola - corte) > 0:
                    k = 1.0 / np.mean(cola - corte)
                    tau_car = 1.0 / k
                else:
                    k = np.nan
                    tau_car = np.nan

                tasas.append(k)
                cantidades.append(N)

                print(
                    f"{corte:5.0f} {N:10d} "
                    f"{k:10.5f} {tau_car:11.3f}"
                )

                # Guardar resultados en tabla
                resultados.append({
                    "A": A,
                    "eta": eta,
                    "Pozo": j + 1,
                    "tau_min": corte,
                    "N_cola": N,
                    "k": k,
                    "Tiempo_caracteristico": tau_car
                })

            # Representacion
            ax = axs[i, j]

            ax.plot(cortes, tasas, "o-")

            ax.set_title(rf"Pozo {j+1}, $\eta={eta}$")
            ax.set_xlabel(r"Tiempo mínimo $\tau_{{\min}}$")
            ax.set_ylabel(r"Tasa de decaimiento $k$")

            # Numeracion en todos los ejes X
            ax.set_xticks(cortes)
            ax.tick_params(axis="x", labelbottom=True)

            ax.grid(alpha=0.3)

            for corte, k, N in zip(cortes, tasas, cantidades):
                if np.isfinite(k):
                    ax.annotate(
                        f"N={N}",
                        (corte, k),
                        xytext=(0, 8),
                        textcoords="offset points",
                        ha="center",
                        fontsize=8
                    )

    fig.suptitle(
        rf"Estabilidad de la tasa exponencial ($A={A}$)",
        fontsize=15
    )

    plt.tight_layout()

    # Guardar grafica
    salida_png = carpeta_salida / "Estabilidad_k_A0.5_subamortiguado.png"

    plt.savefig(
        salida_png,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close()

    # Guardar tabla de resultados
    tabla = pd.DataFrame(resultados)

    salida_csv = carpeta_salida / "Estabilidad_k_A0.5_subamortiguado.csv"

    tabla.to_csv(salida_csv, index=False)

    print("\nGrafica guardada en:", salida_png)
    print("Tabla guardada en:", salida_csv)

    # ============================================================
    # OBJETIVO 2 - DISTRIBUCIÓN DE TIEMPOS DE ESTANCIA
    # Y AJUSTE DE LA COLA EXPONENCIAL
    # SUBAMORTIGUADO - A = 0.5
    # ============================================================

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from pathlib import Path

    # ------------------------------------------------------------
    # RUTA BASE
    # ------------------------------------------------------------
    base = Path.home() / "Documents/A_5-Molecular_dynamics"

    # ------------------------------------------------------------
    # PARÁMETROS
    # ------------------------------------------------------------
    A = 0.5

    # Solo los valores subamortiguados
    etas = [3.5, 5.0]

    # Cortes elegidos a partir del análisis previo de estabilidad
    cortes_por_eta = {
        3.5: 50.0,
        5.0: 50.0
    }

    # Número de bins del histograma
    bins = 50

    # ------------------------------------------------------------
    # CARPETAS
    # ------------------------------------------------------------
    # OJO: aquí van los txt con los tiempos de estancia
    carpeta_datos = base / "Py/plots/OBJETIVO_2/Tiempos_medios_eta/A0.5"

    # Aquí se guardarán las gráficas y la tabla resumen
    carpeta_salida = (
        base / "Py/plots/OBJETIVO_2/Distr_testancia/sobreamortiguado/A0.5"
    )
    carpeta_salida.mkdir(parents=True, exist_ok=True)

    # Si quieres, puedes guardar el resumen en una subcarpeta aparte
    carpeta_resumen = carpeta_salida / "Analisis_cola"
    carpeta_resumen.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # LISTA PARA GUARDAR RESULTADOS
    # ------------------------------------------------------------
    resultados = []

    # ------------------------------------------------------------
    # BUCLE SOBRE ETAS
    # ------------------------------------------------------------
    for eta in etas:

        tau_min = cortes_por_eta[eta]

        # Nombres de archivos
        archivo1 = carpeta_datos / f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
        archivo2 = carpeta_datos / f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"

        # Cargar datos
        est1 = np.loadtxt(archivo1, ndmin=1)
        est2 = np.loadtxt(archivo2, ndmin=1)

        # Figura con dos subplots: pozo 1 y pozo 2
        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        for j, tiempos in enumerate([est1, est2]):

            ax = axs[j]

            # ----------------------------------------------------
            # HISTOGRAMA
            # ----------------------------------------------------
            densidad, bordes, _ = ax.hist(
                tiempos,
                bins=bins,
                density=True,
                alpha=0.35,
                label="Distribución simulada"
            )

            # ----------------------------------------------------
            # AJUSTE DE LA COLA
            # ----------------------------------------------------
            cola = tiempos[tiempos >= tau_min]
            N_total = len(tiempos)
            N_cola = len(cola)

            if N_cola < 2:
                print(f"No hay suficientes datos en la cola para eta={eta}, pozo {j+1}")
                continue

            # Ajuste exponencial de la cola:
            # p(tau | tau >= tau_min) ~ k exp[-k (tau - tau_min)]
            # Entonces, para la distribución total:
            # p(tau) ~ C exp(-k tau), con
            # C = (N_cola/N_total) * k * exp(k*tau_min)
            k = 1.0 / np.mean(cola - tau_min)
            tau_car = 1.0 / k
            frac_cola = N_cola / N_total
            C = frac_cola * k * np.exp(k * tau_min)

            # Rango para dibujar la recta
            x_fit = np.linspace(tau_min, np.max(tiempos), 400)
            y_fit = C * np.exp(-k * x_fit)

            ax.plot(
                x_fit,
                y_fit,
                "--",
                linewidth=2.5,
                label=rf"$p(\tau) = {C:.3e}\,e^{{-{k:.4f}\tau}}$"
            )

            # Línea vertical indicando el corte
            ax.axvline(
                tau_min,
                linestyle=":",
                linewidth=1.8,
                color="gray",
                label=rf"$\tau_{{\min}} = {tau_min:.0f}$"
            )

            # ----------------------------------------------------
            # ESTÉTICA
            # ----------------------------------------------------
            ax.set_yscale("log")
            ax.set_title(f"Pozo {j+1}", fontsize=18)
            ax.set_xlabel(r"Tiempo de estancia $\tau$", fontsize=15)
            ax.set_ylabel(r"Densidad de probabilidad $p(\tau)$", fontsize=15)
            ax.grid(alpha=0.3, which="both")
            ax.legend(fontsize=12)

            # ----------------------------------------------------
            # GUARDAR RESULTADOS EN TABLA
            # ----------------------------------------------------
            resultados.append({
                "A": A,
                "eta": eta,
                "Pozo": j + 1,
                "tau_min": tau_min,
                "N_total": N_total,
                "N_cola": N_cola,
                "fraccion_cola": frac_cola,
                "k": k,
                "Tiempo_caracteristico_1_sobre_k": tau_car,
                "C_prefactor": C
            })

            # Imprimir por pantalla
            print(f"\nA = {A}, eta = {eta}, pozo {j+1}")
            print(f"N total = {N_total}")
            print(f"N cola  = {N_cola}")
            print(f"tau_min = {tau_min}")
            print(f"k       = {k:.6f}")
            print(f"1/k     = {tau_car:.3f}")
            print(f"C       = {C:.6e}")

        # --------------------------------------------------------
        # TÍTULO GENERAL
        # --------------------------------------------------------
        fig.suptitle(
            rf"Distribución de tiempos de estancia y ajuste de la cola exponencial"
            "\n"
            rf"Escala semilogarítmica ($A={A}$, $\eta={eta}$)",
            fontsize=20
        )

        plt.tight_layout()

        # --------------------------------------------------------
        # GUARDAR FIGURA
        # --------------------------------------------------------
        nombre_fig = carpeta_salida / f"Ajuste_cola_A{A:.1f}_eta{eta:.2f}.png"
        plt.savefig(nombre_fig, dpi=300, bbox_inches="tight")

        plt.show()
        plt.close()

        print(f"\nFigura guardada en: {nombre_fig}")

    # ------------------------------------------------------------
    # GUARDAR TABLA RESUMEN
    # ------------------------------------------------------------
    tabla = pd.DataFrame(resultados)

    salida_csv = carpeta_resumen / f"Resumen_ajustes_subamortiguado_A{A:.1f}.csv"
    tabla.to_csv(salida_csv, index=False)

    print("\nTabla resumen guardada en:")
    print(salida_csv)

    # ============================================================
    # OBJETIVO 2 - ESTABILIDAD DEL AJUSTE EXPONENCIAL
    # A = 1.0 - SUBAMORTIGUADO Y SOBREAMORTIGUADO
    # ============================================================

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from pathlib import Path

    base = Path.home() / "Documents/A_5-Molecular_dynamics"

    A = 1.0

    cortes = [0, 20, 50, 100, 150, 200, 300, 500, 800]

    # Separar los dos grupos
    regimenes = {
        "subamortiguado": [0.01, 0.04],
        "sobreamortiguado": [3.5, 5.0]
    }

    # Carpeta donde están los TXT
    carpeta_datos = (
        base / "Py/plots/OBJETIVO_2/Tiempos_medios_eta/A1.0"
    )

    for regimen, etas in regimenes.items():

        # Carpeta de análisis
        carpeta_salida = (
            base / "Py/plots/OBJETIVO_2/Distr_testancia"
            / regimen / "A1.0" / "Analisis_cola"
        )

        carpeta_salida.mkdir(parents=True, exist_ok=True)

        nombre_png = (
            carpeta_salida /
            f"Estabilidad_k_{regimen}_A{A:.1f}.png"
        )

        nombre_csv = (
            carpeta_salida /
            f"Estabilidad_k_{regimen}_A{A:.1f}.csv"
        )

        # Evitar sobrescribir resultados anteriores
        if nombre_png.exists() or nombre_csv.exists():
            raise FileExistsError(
                f"Ya existen resultados en {carpeta_salida}"
            )

        resultados = []

        # Dos filas (etas) y dos columnas (pozos)
        fig, axs = plt.subplots(
            2, 2, figsize=(13, 9), sharex=True
        )

        for i, eta in enumerate(etas):

            archivo1 = carpeta_datos / (
                f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
            )

            archivo2 = carpeta_datos / (
                f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"
            )

            est1 = np.loadtxt(archivo1, ndmin=1)
            est2 = np.loadtxt(archivo2, ndmin=1)

            for j, tiempos in enumerate([est1, est2]):

                tasas = []
                cantidades = []

                print(f"\nA={A}, eta={eta}, pozo {j+1}")
                print("Corte    N datos      k          1/k")

                for corte in cortes:

                    cola = tiempos[tiempos >= corte]
                    N = len(cola)

                    if N >= 2 and np.mean(cola - corte) > 0:
                        k = 1.0 / np.mean(cola - corte)
                        tau_car = 1.0 / k
                    else:
                        k = np.nan
                        tau_car = np.nan

                    tasas.append(k)
                    cantidades.append(N)

                    resultados.append({
                        "A": A,
                        "eta": eta,
                        "Pozo": j + 1,
                        "tau_min": corte,
                        "N_cola": N,
                        "k": k,
                        "Tiempo_caracteristico": tau_car
                    })

                    print(
                        f"{corte:5.0f} {N:9d} "
                        f"{k:10.6f} {tau_car:11.3f}"
                    )

                # Gráfica de estabilidad
                ax = axs[i, j]

                ax.plot(cortes, tasas, "o-", markersize=4)

                ax.set_title(rf"Pozo {j+1}, $\eta={eta}$")
                ax.set_xlabel(r"Tiempo mínimo $\tau_{\min}$")
                ax.set_ylabel(r"Tasa de decaimiento $k$")

                # Números del eje X en todos los paneles
                ax.set_xticks(cortes)
                ax.tick_params(
                    axis="x", labelbottom=True, labelrotation=35
                )

                ax.grid(alpha=0.3)

                for corte, k, N in zip(
                    cortes, tasas, cantidades
                ):
                    if np.isfinite(k):
                        ax.annotate(
                            f"N={N}",
                            (corte, k),
                            xytext=(0, 8),
                            textcoords="offset points",
                            ha="center",
                            fontsize=7
                        )

        fig.suptitle(
            f"Estabilidad de la tasa exponencial\n"
            rf"Régimen {regimen} ($A={A}$)",
            fontsize=15
        )

        plt.tight_layout()

        # Guardar figura
        plt.savefig(
            nombre_png, dpi=300, bbox_inches="tight"
        )

        plt.show()
        plt.close()

        # Guardar CSV
        tabla = pd.DataFrame(resultados)
        tabla.to_csv(nombre_csv, index=False)

        print(f"\nGráfica guardada: {nombre_png}")
        print(f"Tabla guardada: {nombre_csv}")

    # ============================================================
    # OBJETIVO 2 - DISTRIBUCION DE TIEMPOS DE ESTANCIA
    # Y AJUSTE DE LA COLA EXPONENCIAL
    # SUBAMORTIGUADO - A = 1.0
    # ============================================================

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from pathlib import Path

    # ------------------------------------------------------------
    # PARAMETROS
    # ------------------------------------------------------------

    base = Path.home() / "Documents/A_5-Molecular_dynamics"

    A = 1.0
    etas = [0.01, 0.04]

    # Cortes obtenidos del analisis de estabilidad
    cortes_por_eta = {
        0.01: 300.0,
        0.04: 150.0
    }

    bins = 50

    # ------------------------------------------------------------
    # CARPETAS
    # ------------------------------------------------------------

    # Archivos TXT originales
    carpeta_datos = (
        base / "Py/plots/OBJETIVO_2/Tiempos_medios_eta/A1.0"
    )

    # Graficas definitivas
    carpeta_salida = (
        base / "Py/plots/OBJETIVO_2/Distr_testancia/subamortiguado/A1.0"
    )
    carpeta_salida.mkdir(parents=True, exist_ok=True)

    # Tabla resumen
    carpeta_resumen = carpeta_salida / "Analisis_cola"
    carpeta_resumen.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # COMPROBAR QUE NO VAMOS A SOBRESCRIBIR ARCHIVOS
    # ------------------------------------------------------------

    archivos_salida = [
        carpeta_salida /
        f"Ajuste_cola_subamortiguado_A{A:.1f}_eta{eta:.2f}.png"
        for eta in etas
    ]

    archivo_csv = (
        carpeta_resumen /
        f"Resumen_ajustes_subamortiguado_A{A:.1f}.csv"
    )

    for archivo in archivos_salida + [archivo_csv]:
        if archivo.exists():
            raise FileExistsError(
                f"El archivo ya existe y no se sobrescribira: {archivo}"
            )

    # ------------------------------------------------------------
    # AJUSTES EXPONENCIALES
    # ------------------------------------------------------------

    resultados = []

    for eta, nombre_fig in zip(etas, archivos_salida):

        tau_min = cortes_por_eta[eta]

        archivo1 = carpeta_datos / (
            f"estancia_pozo1_A{A:.1f}_eta{eta:.2f}.txt"
        )

        archivo2 = carpeta_datos / (
            f"estancia_pozo2_A{A:.1f}_eta{eta:.2f}.txt"
        )

        est1 = np.loadtxt(archivo1, ndmin=1)
        est2 = np.loadtxt(archivo2, ndmin=1)

        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        for j, tiempos in enumerate([est1, est2]):

            ax = axs[j]

            # Histograma normalizado de todas las estancias
            ax.hist(
                tiempos,
                bins=bins,
                density=True,
                alpha=0.35,
                label="Distribucion simulada"
            )

            # Seleccionar la cola exponencial
            cola = tiempos[tiempos >= tau_min]

            N_total = len(tiempos)
            N_cola = len(cola)

            if N_cola < 2:
                raise ValueError(
                    f"Datos insuficientes: eta={eta}, pozo {j+1}"
                )

            # Estimacion de maxima verosimilitud
            k = 1.0 / np.mean(cola - tau_min)

            tau_car = 1.0 / k

            fraccion_cola = N_cola / N_total

            # Modelo p(tau) = C exp(-k*tau)
            C = fraccion_cola * k * np.exp(k * tau_min)

            x_fit = np.linspace(
                tau_min, np.max(tiempos), 500
            )

            y_fit = C * np.exp(-k * x_fit)

            # Formato cientifico para la ecuacion
            exponente = int(np.floor(np.log10(C)))
            coeficiente = C / 10**exponente

            ecuacion = (
                rf"$p(\tau)={coeficiente:.3f}\times10^{{{exponente}}}"
                rf"e^{{-{k:.4f}\tau}}$"
            )

            # Dibujar ajuste
            ax.plot(
                x_fit, y_fit, "--",
                linewidth=2.5,
                label=ecuacion
            )

            # Indicar el comienzo de la cola
            ax.axvline(
                tau_min,
                linestyle=":",
                linewidth=1.5,
                color="gray",
                label=rf"$\tau_{{\min}}={tau_min:g}$"
            )

            # Representacion semilogaritmica
            ax.set_yscale("log")

            ax.set_title(f"Pozo {j+1}", fontsize=16)
            ax.set_xlabel(r"Tiempo de estancia $\tau$")
            ax.set_ylabel(r"Densidad de probabilidad $p(\tau)$")

            ax.grid(alpha=0.3, which="both")
            ax.legend(fontsize=10)

            # Guardar resultados numericos
            resultados.append({
                "A": A,
                "eta": eta,
                "Pozo": j + 1,
                "tau_min": tau_min,
                "N_total": N_total,
                "N_cola": N_cola,
                "fraccion_cola": fraccion_cola,
                "k": k,
                "Tiempo_caracteristico": tau_car,
                "C": C
            })

            print(f"\nA={A}, eta={eta}, pozo {j+1}")
            print(f"N total = {N_total}")
            print(f"N cola = {N_cola}")
            print(f"k = {k:.6f}")
            print(f"1/k = {tau_car:.3f}")

        # Titulo general
        fig.suptitle(
            "Distribucion de tiempos de estancia y ajuste de la cola exponencial\n"
            rf"Escala semilogaritmica ($A={A}$, $\eta={eta}$)",
            fontsize=16
        )

        plt.tight_layout()

        # Guardar sin sobrescribir archivos antiguos
        plt.savefig(
            nombre_fig,
            dpi=300,
            bbox_inches="tight"
        )

        plt.show()
        plt.close()

    # ------------------------------------------------------------
    # GUARDAR TABLA RESUMEN
    # ------------------------------------------------------------

    tabla = pd.DataFrame(resultados)

    tabla.to_csv(archivo_csv, index=False)

    print("\nTabla resumen guardada en:")
    print(archivo_csv)

#Term Polimero
if False:
    ruta = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\Termalizacion_polimero.txt"
    N = 4           # número de monómeros
    ke = 100.0       # constante del muelle
    beta_inv = 1.0   # = kBT
    l0=1
    h = 0.001
    frac_zoom = 0.2  # fracción final de la simulación que se amplía


    datos = np.loadtxt(ruta, ndmin=2)
    t  = datos[:, 0]
    Ec = datos[:, 1]
    Ep = datos[:, 2]

    Ec_esp = 1.5 * N * beta_inv
    Ep_esp = 0.5 * (N - 1) * beta_inv
    Ep_esp2 = ke/2 * (3*beta_inv*beta_inv/(ke^2) + l0*l0*beta_inv/ke)/(beta_inv/ke + l0*l0)

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(rf"Termalización polímero ($N={N}$, $k_e={ke:g}$)")

    # Zona del zoom
    i0 = int(len(t) * (1 - frac_zoom))
    tz = t[i0:]

    configs = [
        dict(ax=axs[0], y=Ec, esp=Ec_esp, color="tab:blue",
             lab_esp=rf"$3Nk_BT/2={Ec_esp:.3f}$",
             ylabel="Energía cinética media"),
        dict(ax=axs[1], y=Ep, esp=Ep_esp, color="tab:green",
             lab_esp=rf"$(N-1)k_BT/2={Ep_esp:.3f}, aportacion '1D'$",
             ylabel="Energía potencial media '1D'"),
        dict(ax=axs[1], y=Ep, esp=Ep_esp2, color="tab:red",
             lab_esp=rf"<V>teo ={Ep_esp2:.3f}$, aportacion '3D'",
             ylabel="Energía potencial media '3D'"),
    ]

    for c in configs:
        ax, y, esp, color = c["ax"], c["y"], c["esp"], c["color"]

        ax.plot(t, y, color=color, linewidth=0.8)
        ax.axhline(esp, color="k", linestyle="--", linewidth=0.8,
                   label=c["lab_esp"])
        ax.set_xlabel("Tiempo")
        ax.set_ylabel(c["ylabel"])
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left")

        # Ventana de zoom
        yz = y[i0:]
        axins = ax.inset_axes([0.45, 0.2, 0.5, 0.4])
        axins.plot(tz, yz, color=color, linewidth=0.8)
        axins.axhline(esp, color="k", linestyle="--", linewidth=0.8)
        axins.set_xlim(tz[0], tz[-1])
        ymin = min(yz.min(), esp)
        ymax = max(yz.max(), esp)
        margen = 0.05 * (ymax - ymin) if ymax > ymin else 0.05
        axins.set_ylim(ymin - margen, ymax + margen)
        axins.grid(alpha=0.3)
        axins.tick_params(labelsize=8)
        ax.indicate_inset_zoom(axins, edgecolor="black")

    plt.tight_layout()
    plt.savefig(
        fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\Equip_polimero_N{N}_ke{ke}.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.show()
    plt.close()

#Valor medio Rg^2
if True:
    ruta = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\Radio_giro.txt"
    N = 6             # número de monómeros
    ke = 100.0        # constante del muelle
    h = 0.001
    frac_zoom = 0.5   # fracción final de la simulación que se amplía

    datos = np.loadtxt(ruta, ndmin=2)
    t   = datos[:, 0]
    rg2 = datos[:, 1]

    final = rg2[-1]

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(t, rg2, color="tab:red", linewidth=1.0)

    ax.set_title(rf"Radio de giro del polímero ($N={N}$, $k_e={ke:g}$)")
    ax.set_xlabel("Tiempo")
    ax.set_ylabel(r"$\langle R_G^2\rangle$")
    ax.grid(alpha=0.3)

    # Leyenda de solo texto con el valor final
    ax.plot([], [], " ", label=rf"$\langle R_G^2\rangle_f={final:.3f}$")
    ax.legend(loc="upper right", handlelength=0, handletextpad=0)

    # ---------- Ventana de zoom ----------
    i0 = int(len(t) * (1 - frac_zoom))
    tz, yz = t[i0:], rg2[i0:]

    axins = ax.inset_axes([0.3, 0.3, 0.65, 0.4])   # [x, y, ancho, alto]
    axins.plot(tz, yz, color="tab:red", linewidth=0.8)
    axins.set_xlim(tz[0], tz[-1])
    margen = 0.1 * (yz.max() - yz.min())
    axins.set_ylim(yz.min() - margen, yz.max() + margen)
    axins.locator_params(axis="x", nbins=4)
    axins.grid(alpha=0.3)
    axins.tick_params(labelsize=8)

    ax.indicate_inset_zoom(axins, edgecolor="black")

    plt.tight_layout()
    plt.savefig(
        fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\Rg2_N{N}_ke{ke:g}.png",
        dpi=300, bbox_inches="tight"
    )
    plt.show()
    plt.close()


#plot scatter
if False:
    ruta = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\valores_rg_2.txt"
    ke = 100.0        # constante del muelle

    datos = np.loadtxt(ruta, comments="#", ndmin=2)
    N   = datos[:, 0]
    rg2 = datos[:, 1]

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(N, rg2, color="navy", s=25)

    ax.set_title(rf"Radio de giro frente al tamaño del polímero ($k_e={ke:g}$)")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\langle R_G^2\rangle$")

    plt.tight_layout()
    plt.savefig(
        fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\Rg2_vs_N_puntos_ke{ke:g}.png",
        dpi=300, bbox_inches="tight"
    )
    plt.show()
    plt.close()

#Regresión lineal
if False:
    ruta = r"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Fortran\valores_rg_2.txt"
    ke = 100.0        # constante del muelle

    datos = np.loadtxt(ruta, comments="#", ndmin=2)
    N   = datos[:, 0]
    rg2 = datos[:, 1]

    # Regresión lineal en log-log: ln(Rg2) = alpha*ln(N) + ln(A)
    res = linregress(np.log(N), np.log(rg2))

    alpha, alpha_err = res.slope, res.stderr
    c, c_err = res.intercept, res.intercept_stderr
    r2 = res.rvalue**2

    nu, nu_err = alpha / 2, alpha_err / 2
    A = np.exp(c)
    A_err = A * c_err        # propagación: dA = A * dc

    print("=== RESULTADOS DE LA REGRESIÓN LINEAL (log-log) ===")
    print(f"Exponente alpha (Rg^2 ~ N^alpha) : {alpha:.4f} ± {alpha_err:.4f}")
    print(f"Exponente de escala nu = alpha/2 : {nu:.4f} ± {nu_err:.4f}")
    print(f"Constante A                      : {A:.4f} ± {A_err:.4f}")
    print(f"R^2                              : {r2:.5f}")

    N_fit = np.linspace(N.min(), N.max(), 200)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.scatter(N, rg2, color="navy", s=25, zorder=3, label="Simulación")

    ax.plot(N_fit, A * N_fit**alpha, "k--", linewidth=1.2,
            label=rf"Ajuste: $\langle R_G^2\rangle=({A:.3f}\pm{A_err:.3f})\,N^{{{alpha:.3f}\pm{alpha_err:.3f}}}$" "\n"
                  rf"$\nu={nu:.3f}\pm{nu_err:.3f}$,  $R^2={r2:.5f}$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(rf"Ley de escala del radio de giro ($k_e={ke:g}$)")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$\langle R_G^2\rangle$")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(
        fr"C:\Users\MSI\Desktop\4 FISICA\TeFi_III\A_5\Py\plots\Rg2_vs_N_loglog_ke{ke:g}.png",
        dpi=300, bbox_inches="tight"
    )
    plt.show()
    plt.close()




# %%
# ============================================================
# OBJETIVO 3 - EQUIPARTICION DEL POLIMERO
# ke = 1000, N = 8, 16, 32, 64
# GRAFICAS CON ZOOM Y FORMULAS TEORICAS
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"

# Parametros
valores_N = [8, 16, 32, 64]

ke = 1000
kBT = 1.0
l0 = 1.0

# Carpetas
carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k1000/EQUIP_E"
)

carpeta_salida = (
    base / "Py/plots/OBJETIVO_3/EQUIP_E/k1000"
)

carpeta_salida.mkdir(parents=True, exist_ok=True)

# Comprobar archivos de salida antes de generar figuras
for N in valores_N:
    salida = carpeta_salida / f"Equip_polimero_N{N}_ke{ke}_zoom.png"

    if salida.exists():
        raise FileExistsError(
            f"El archivo ya existe: {salida}"
        )

# ============================================================
# BUCLE SOBRE LOS DISTINTOS TAMANOS
# ============================================================

for N in valores_N:

    archivo = (
        carpeta_datos /
        f"Termalizacion_polimero_N{N}_ke{ke}.txt"
    )

    # Leer datos de Fortran
    t, K_m, V_m = np.loadtxt(archivo, unpack=True)

    # Valores teoricos
    K_teor = 3 * N * kBT / 2

    sigma2 = kBT / ke

    V_enlace = (kBT / 2) * (
        (l0**2 + 3*sigma2) /
        (l0**2 + sigma2)
    )

    V_teor = (N - 1) * V_enlace

    # Imprimir resultados
    print(f"\nN = {N}, ke = {ke}")
    print(f"Tiempo final = {t[-1]:.2f}")
    print(f"<K> simulacion = {K_m[-1]:.5f}")
    print(f"<K> teorico    = {K_teor:.5f}")
    print(f"<V> simulacion = {V_m[-1]:.5f}")
    print(f"<V> teorico    = {V_teor:.5f}")

    # ========================================================
    # CREAR FIGURA
    # ========================================================

    fig, axs = plt.subplots(1, 2, figsize=(15, 6))

    datos = [
        (K_m, K_teor, "Energía cinética media", "tab:blue"),
        (V_m, V_teor, "Energía potencial media", "tab:green")
    ]

    for ax, (energia, valor_teor, titulo, color) in zip(axs, datos):

        # Grafica principal
        ax.plot(
            t, energia,
            color=color,
            linewidth=1
        )

        # Formulas en la leyenda
        if ax is axs[0]:
            etiqueta = (
                rf"$\langle K\rangle = "
                rf"\frac{{3N}}{{2}}k_BT = {K_teor:.3f}$"
            )
        else:
            etiqueta = (
                rf"$\langle V\rangle \simeq "
                rf"\frac{{N-1}}{{2}}k_BT"
                rf" = {V_teor:.3f}$"
            )

        # Linea teorica
        ax.axhline(
            valor_teor,
            color="black",
            linestyle="--",
            linewidth=1,
            label=etiqueta
        )

        ax.set_title(titulo)
        ax.set_xlabel("Tiempo")
        ax.set_ylabel("Energía media")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left")

        # ====================================================
        # ZOOM: ULTIMAS 6000 UNIDADES DE TIEMPO
        # ====================================================

        t_inicio = max(t[0], t[-1] - 6000)
        mascara = t >= t_inicio

        ax_zoom = ax.inset_axes(
            [0.45, 0.12, 0.52, 0.43]
        )

        ax_zoom.plot(
            t, energia,
            color=color,
            linewidth=1
        )

        ax_zoom.axhline(
            valor_teor,
            color="black",
            linestyle="--",
            linewidth=1
        )

        ax_zoom.set_xlim(t_inicio, t[-1])

        valores_zoom = energia[mascara]

        y_min = min(np.min(valores_zoom), valor_teor)
        y_max = max(np.max(valores_zoom), valor_teor)

        margen = max((y_max - y_min) * 0.15, 0.005)

        ax_zoom.set_ylim(
            y_min - margen,
            y_max + margen
        )

        ax_zoom.grid(alpha=0.3)
        ax_zoom.tick_params(
            axis="both",
            labelsize=8
        )

        # Recuadro y lineas de conexion
        ax.indicate_inset_zoom(
            ax_zoom,
            edgecolor="gray"
        )

    # Titulo general
    fig.suptitle(
        rf"Termalización del polímero ($N={N}$, $k_e={ke}$)",
        fontsize=16
    )

    plt.tight_layout()

    # ========================================================
    # GUARDAR FIGURA
    # ========================================================

    salida = (
        carpeta_salida /
        f"Equip_polimero_N{N}_ke{ke}_zoom.png"
    )

    if salida.exists():
        plt.close(fig)
        raise FileExistsError(
            f"No se sobrescribira el archivo: {salida}"
        )

    plt.savefig(
        salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close()

    print("Gráfica guardada en:", salida)

print("\nTodas las gráficas generadas.")

# %%

# %%
# ============================================================
# OBJETIVO 3 - EQUIPARTICION DEL POLIMERO
# ke = 1000
# N = 4, 8, 16, 32, 64
# Inicializacion termica - Sin zoom
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Directorio principal
base = Path.home() / "Documents/A_5-Molecular_dynamics"

# Parametros
valores_N = [4, 8, 16, 32, 64]

ke = 1000
kBT = 1.0
l0 = 1.0

# Carpetas
carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k1000/EQUIP_E"
)

carpeta_salida = (
    base / "Py/plots/OBJETIVO_3/EQUIP_E/k1000"
)

carpeta_salida.mkdir(parents=True, exist_ok=True)

# Configuracion de las graficas
t_inicio_grafica = 5000

# Margen del eje vertical
margen_relativo = 0.35

# ============================================================
# 1. COMPROBAR ARCHIVOS
# ============================================================

for N in valores_N:

    archivo = (
        carpeta_datos /
        f"Termalizacion_polimero_N{N}_ke{ke}_inicio_termico.txt"
    )

    salida = (
        carpeta_salida /
        f"Equip_polimero_N{N}_ke{ke}_inicio_termico_sin_zoom.png"
    )

    if not archivo.is_file():
        raise FileNotFoundError(
            f"No se encuentra el archivo: {archivo}"
        )

    if salida.exists():
        raise FileExistsError(
            f"La grafica ya existe y no se sobrescribira: {salida}"
        )

# ============================================================
# 2. GENERAR LAS CINCO FIGURAS
# ============================================================

for N in valores_N:

    # Cargar datos
    archivo = (
        carpeta_datos /
        f"Termalizacion_polimero_N{N}_ke{ke}_inicio_termico.txt"
    )

    t, K_m, V_m = np.loadtxt(archivo, unpack=True)

    # ========================================================
    # VALORES TEORICOS
    # ========================================================

    # Energia cinetica
    K_teor = 3 * N * kBT / 2

    # Energia potencial, con correccion geometrica 3D
    sigma2 = kBT / ke

    V_enlace = (kBT / 2) * (
        (l0**2 + 3*sigma2) /
        (l0**2 + sigma2)
    )

    V_teor = (N - 1) * V_enlace

    # Mostrar resultados numericos
    print(f"\n{'=' * 45}")
    print(f"N = {N}  |  ke = {ke}")
    print(f"{'=' * 45}")

    print(f"Tiempo de medicion = {t[-1]:.0f}")
    print(f"<K> simulacion = {K_m[-1]:.5f}")
    print(f"<K> teorico    = {K_teor:.5f}")
    print(f"<V> simulacion = {V_m[-1]:.5f}")
    print(f"<V> teorico 3D = {V_teor:.5f}")

    # ========================================================
    # SELECCIONAR INTERVALO
    # ========================================================

    mascara = t >= t_inicio_grafica

    t_graf = t[mascara]
    K_graf = K_m[mascara]
    V_graf = V_m[mascara]

    if len(t_graf) == 0:
        raise ValueError(
            f"No hay datos posteriores a t={t_inicio_grafica} "
            f"para N={N}"
        )

    # ========================================================
    # CREAR FIGURA
    # ========================================================

    fig, axs = plt.subplots(
        1, 2,
        figsize=(15, 5.5)
    )

    datos = [
        (
            K_graf,
            K_teor,
            "Energía cinética media",
            "tab:blue"
        ),
        (
            V_graf,
            V_teor,
            "Energía potencial media",
            "tab:green"
        )
    ]

    for ax, (energia, valor_teor, titulo, color) in zip(
        axs, datos
    ):

        # Curva de simulacion
        ax.plot(
            t_graf,
            energia,
            color=color,
            linewidth=1.2
        )

        # ====================================================
        # LEYENDAS TEORICAS
        # ====================================================

        if ax is axs[0]:

            etiqueta = (
                rf"$\langle K\rangle = "
                rf"\frac{{3N}}{{2}}k_BT"
                rf" = {K_teor:.3f}$"
            )

        else:

            etiqueta = (
                rf"$\langle V\rangle \simeq "
                rf"\frac{{N-1}}{{2}}k_BT$"
                "\n"
                rf"Referencia 3D = {V_teor:.3f}"
            )

        # Linea teorica
        ax.axhline(
            valor_teor,
            color="black",
            linestyle="--",
            linewidth=1,
            label=etiqueta
        )

        # ====================================================
        # FORMATO DE LA GRAFICA
        # ====================================================

        ax.set_title(titulo)
        ax.set_xlabel("Tiempo de medición")
        ax.set_ylabel("Energía media")

        ax.grid(alpha=0.3)
        ax.legend(loc="upper left")

        ax.set_xlim(
            t_graf[0],
            t_graf[-1]
        )

        # ====================================================
        # ESCALA DEL EJE Y CON MARGEN AMPLIO
        # ====================================================

        y_min = min(np.min(energia), valor_teor)
        y_max = max(np.max(energia), valor_teor)

        margen = max(
            margen_relativo * (y_max - y_min),
            0.01 * abs(valor_teor),
            0.03
        )

        ax.set_ylim(
            y_min - margen,
            y_max + margen
        )

    # ========================================================
    # TITULO GENERAL
    # ========================================================

    fig.suptitle(
        rf"Equipartición del polímero ($N={N}$, $k_e={ke}$)",
        fontsize=15
    )

    fig.tight_layout()

    # ========================================================
    # GUARDAR GRAFICA
    # ========================================================

    salida = (
        carpeta_salida /
        f"Equip_polimero_N{N}_ke{ke}_inicio_termico_sin_zoom.png"
    )

    fig.savefig(
        salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print(f"Grafica guardada en: {salida}")

print("\nTodas las graficas generadas correctamente.")



# %%



# %%
# ============================================================
# OBJETIVO 3 - RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS
# ke = 100 - SIMULACIONES CON TIEMPO AMPLIADO
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

valores_N = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]

ke = 100
kBT = 1.0
l0 = 1.0

margen_relativo = 0.30
mostrar_bloques = False

# Carpetas
carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k100/RG"
)

carpeta_salida = (
    base / "Py/plots/OBJETIVO_3/RG/k100"
)

carpeta_salida.mkdir(parents=True, exist_ok=True)

# ============================================================
# VALORES TEORICOS
# ============================================================

sigma2 = kBT / ke

b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

# ============================================================
# GENERAR LAS GRAFICAS
# ============================================================

for N in valores_N:

    archivo = (
        carpeta_datos /
        f"Rg_Ree_N{N}_ke{ke}.txt"
    )

    salida = (
        carpeta_salida /
        f"Rg_Ree_N{N}_ke{ke}_tiempo_ampliado.png"
    )

    if not archivo.is_file():
        print(f"No se encuentra el TXT de N={N}. Se omite.")
        continue

    if salida.exists():
        print(f"La gráfica de N={N} ya existe. Se omite.")
        continue

    # ========================================================
    # 1. CARGAR DATOS
    # ========================================================

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 5:
        raise ValueError(
            f"N={N}: se esperaban 5 columnas, "
            f"pero hay {datos.shape[1]}"
        )

    t = datos[:, 0]

    Rg2_m = datos[:, 1]
    Ree2_m = datos[:, 2]

    Rg2_bloques = datos[:, 3]
    Ree2_bloques = datos[:, 4]

    # ========================================================
    # 2. VALORES TEORICOS
    # ========================================================

    Rg2_teor = (N**2 - 1) * b2 / (6*N)
    Ree2_teor = (N - 1) * b2

    # ========================================================
    # 3. MOSTRAR RESULTADOS
    # ========================================================

    print(f"\n{'=' * 45}")
    print(f"N = {N}, ke = {ke}")
    print(f"Tiempo final de medición = {t[-1]:.0f}")
    print(f"{'=' * 45}")

    print(f"<Rg²> simulación = {Rg2_m[-1]:.6f}")
    print(f"<Rg²> teórico    = {Rg2_teor:.6f}")

    print(f"<Ree²> simulación = {Ree2_m[-1]:.6f}")
    print(f"<Ree²> teórico    = {Ree2_teor:.6f}")

    # ========================================================
    # 4. CREAR FIGURA
    # ========================================================

    fig, axs = plt.subplots(
        1, 2,
        figsize=(15, 5.5)
    )

    series = [
        (
            Rg2_m,
            Rg2_bloques,
            Rg2_teor,
            "Radio de giro al cuadrado",
            r"$\langle R_g^2\rangle_t$",
            "tab:blue",
            rf"Cadena ideal: $\frac{{N^2-1}}{{6N}}b^2"
            rf" = {Rg2_teor:.3f}$"
        ),
        (
            Ree2_m,
            Ree2_bloques,
            Ree2_teor,
            "Distancia entre extremos al cuadrado",
            r"$\langle R_{ee}^2\rangle_t$",
            "tab:green",
            rf"Cadena ideal: $(N-1)b^2"
            rf" = {Ree2_teor:.3f}$"
        )
    ]

    for ax, (
        media, bloques, teor, titulo,
        ylabel, color, etiqueta
    ) in zip(axs, series):

        # Promedio acumulado
        ax.plot(
            t,
            media,
            color=color,
            linewidth=1.5,
            label="Promedio acumulado"
        )

        # Promedios por bloques, opcionales
        if mostrar_bloques:
            ax.plot(
                t,
                bloques,
                color=color,
                alpha=0.35,
                linewidth=1,
                label="Promedio por bloques"
            )

        # Valor teorico
        ax.axhline(
            teor,
            color="black",
            linestyle="--",
            linewidth=1.1,
            label=etiqueta
        )

        ax.set_title(titulo)
        ax.set_xlabel("Tiempo de medición")
        ax.set_ylabel(ylabel)

        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=9)

        # Escala X automatica
        ax.set_xlim(t[0], t[-1])

        # Escala Y con margen
        if mostrar_bloques:
            minimo = min(np.min(media), np.min(bloques), teor)
            maximo = max(np.max(media), np.max(bloques), teor)
        else:
            minimo = min(np.min(media), teor)
            maximo = max(np.max(media), teor)

        margen = max(
            margen_relativo * (maximo - minimo),
            0.05 * teor
        )

        ax.set_ylim(
            max(0, minimo - margen),
            maximo + margen
        )

    # ========================================================
    # 5. TITULO Y GUARDADO
    # ========================================================

    fig.suptitle(
        rf"Evolución del polímero ($N={N}$, $k_e={ke}$)",
        fontsize=15
    )

    fig.tight_layout()

    fig.savefig(
        salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print(f"Gráfica guardada en: {salida}")

print("\nTodas las gráficas disponibles se han generado.")



# %%
import numpy as np
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"
carpeta = base / "Py/plots/OBJETIVO_3/Datos/k100/RG"

for N in [4, 32, 64]:

    archivo = carpeta / f"Rg_Ree_N{N}_ke100.txt"
    datos = np.loadtxt(archivo, ndmin=2)

    t = datos[:, 0]
    rg2 = datos[:, 3]   # Medias por bloques
    ree2 = datos[:, 4]

    print(f"\n======== N = {N} ========")

    # Dividir la medicion en cuatro partes
    indices = np.array_split(np.arange(len(t)), 4)

    for i, idx in enumerate(indices, 1):
        print(
            f"Bloque {i}: "
            f"<Rg²> = {np.mean(rg2[idx]):.4f}  |  "
            f"<Ree²> = {np.mean(ree2[idx]):.4f}"
        )




# %%
# ============================================================
# OBJETIVO 3 - RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS
# N = 64, ke = 100
# Termalizacion: 40000
# Medicion: 200000
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# PARAMETROS
base = Path.home() / "Documents/A_5-Molecular_dynamics"

N = 64
ke = 100
kBT = 1.0
l0 = 1.0

margen_relativo = 0.30

# CARPETAS
carpeta_datos = base / "Py/plots/OBJETIVO_3/Datos/k100/RG"
carpeta_salida = base / "Py/plots/OBJETIVO_3/RG/k100"

carpeta_salida.mkdir(parents=True, exist_ok=True)

# ARCHIVO DE LA NUEVA SIMULACION
archivo = carpeta_datos / "Rg_Ree_N64_ke100_T200000.txt"

# Nombre diferente para conservar las graficas anteriores
salida = carpeta_salida / "Rg_Ree_N64_ke100_T200000.png"

if not archivo.is_file():
    raise FileNotFoundError(f"No se encuentra: {archivo}")

if salida.exists():
    raise FileExistsError(
        f"La grafica ya existe y no se sobrescribira: {salida}"
    )

# ============================================================
# 1. CARGAR LOS DATOS
# ============================================================

datos = np.loadtxt(archivo, ndmin=2)

if datos.shape[1] != 5:
    raise ValueError("El TXT debe contener cinco columnas.")

t = datos[:, 0]

Rg2_m = datos[:, 1]
Ree2_m = datos[:, 2]

Rg2_bloques = datos[:, 3]
Ree2_bloques = datos[:, 4]

print(f"Tiempo final de medicion: {t[-1]:.0f}")

if not np.isclose(t[-1], 200000):
    print("AVISO: el archivo no llega a t = 200000.")

# ============================================================
# 2. VALORES TEORICOS
# ============================================================

sigma2 = kBT / ke

b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

Rg2_teor = (N**2 - 1) * b2 / (6*N)
Ree2_teor = (N - 1) * b2

print("\n===== RESULTADOS FINALES =====")

print(f"<Rg²> simulacion = {Rg2_m[-1]:.6f}")
print(f"<Rg²> teorico    = {Rg2_teor:.6f}")

print(f"<Ree²> simulacion = {Ree2_m[-1]:.6f}")
print(f"<Ree²> teorico    = {Ree2_teor:.6f}")

# ============================================================
# 3. COMPROBAR LOS CUATRO BLOQUES TEMPORALES
# ============================================================

print("\n===== CONVERGENCIA POR BLOQUES =====")

indices = np.array_split(np.arange(len(t)), 4)

for i, idx in enumerate(indices, 1):

    media_rg = np.mean(Rg2_bloques[idx])
    media_ree = np.mean(Ree2_bloques[idx])

    print(
        f"Bloque {i}: "
        f"<Rg²> = {media_rg:.4f}  |  "
        f"<Ree²> = {media_ree:.4f}"
    )

# ============================================================
# 4. CREAR LA FIGURA
# ============================================================

fig, axs = plt.subplots(
    1, 2,
    figsize=(15, 5.5)
)

series = [
    (
        Rg2_m,
        Rg2_teor,
        "Radio de giro al cuadrado",
        r"$\langle R_g^2\rangle_t$",
        "tab:blue",
        rf"Cadena ideal: $\frac{{N^2-1}}{{6N}}b^2"
        rf" = {Rg2_teor:.3f}$"
    ),
    (
        Ree2_m,
        Ree2_teor,
        "Distancia entre extremos al cuadrado",
        r"$\langle R_{ee}^2\rangle_t$",
        "tab:green",
        rf"Cadena ideal: $(N-1)b^2"
        rf" = {Ree2_teor:.3f}$"
    )
]

for ax, (media, teor, titulo, ylabel,
         color, etiqueta) in zip(axs, series):

    # Promedio acumulado
    ax.plot(
        t,
        media,
        color=color,
        linewidth=1.5,
        label="Promedio acumulado"
    )

    # Valor teorico
    ax.axhline(
        teor,
        color="black",
        linestyle="--",
        linewidth=1.1,
        label=etiqueta
    )

    ax.set_title(titulo)
    ax.set_xlabel("Tiempo de medición")
    ax.set_ylabel(ylabel)

    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)

    ax.set_xlim(t[0], t[-1])

    # Margen del eje Y
    minimo = min(np.min(media), teor)
    maximo = max(np.max(media), teor)

    margen = max(
        margen_relativo * (maximo - minimo),
        0.05 * teor
    )

    ax.set_ylim(
        max(0, minimo - margen),
        maximo + margen
    )

# ============================================================
# 5. TITULO Y GUARDADO
# ============================================================

fig.suptitle(
    rf"Evolución del polímero ($N={N}$, $k_e={ke}$)",
    fontsize=15
)

fig.tight_layout()

fig.savefig(
    salida,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close(fig)

print(f"\nGrafica guardada en: {salida}")



# %%
# ============================================================
# OBJETIVO 3 - RADIO DE GIRO
# AJUSTE LINEAL Y LOG-LOG
# ke = 100
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

valores_N = np.array(
    [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
)

ke = 100
kBT = 1.0
l0 = 1.0

# Valor minimo de N incluido en los ajustes
N_ajuste_min = 4

# Carpetas
carpeta_datos = base / "Py/plots/OBJETIVO_3/Datos/k100/RG"
carpeta_salida = base / "Py/plots/OBJETIVO_3/RG/k100"

carpeta_salida.mkdir(parents=True, exist_ok=True)

salida_lineal = (
    carpeta_salida /
    "Rg2_vs_N_ke100_lineal_sin_errores.png"
)

salida_loglog = (
    carpeta_salida /
    "Rg2_vs_N_ke100_loglog_sin_errores.png"
)

# Proteger graficas anteriores
for salida in [salida_lineal, salida_loglog]:
    if salida.exists():
        raise FileExistsError(
            f"El archivo ya existe: {salida}"
        )

# ============================================================
# 2. LEER LOS DATOS DE LAS 13 SIMULACIONES
# ============================================================

Rg2_medios = []

for N in valores_N:

    if N == 64:
        nombre = "Rg_Ree_N64_ke100_T200000.txt"
    else:
        nombre = f"Rg_Ree_N{N}_ke100.txt"

    archivo = carpeta_datos / nombre

    if not archivo.is_file():
        raise FileNotFoundError(
            f"No se encuentra: {archivo}"
        )

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 5:
        raise ValueError(
            f"El archivo N={N} no tiene cinco columnas."
        )

    # Media acumulada final de Rg²
    media = datos[-1, 1]

    if not np.isfinite(media) or media <= 0:
        raise ValueError(
            f"Valor invalido de Rg² para N={N}"
        )

    Rg2_medios.append(media)

    print(
        f"N={N:2d}  |  "
        f"Tiempo={datos[-1,0]:.0f}  |  "
        f"<Rg²>={media:.6f}"
    )

Rg2_medios = np.array(Rg2_medios)

# Puntos que entran en los ajustes
mascara_ajuste = valores_N >= N_ajuste_min

N_fit = valores_N[mascara_ajuste]
Rg2_fit = Rg2_medios[mascara_ajuste]

# ============================================================
# 3. PREDICCION TEORICA DE LA CADENA IDEAL
# ============================================================

sigma2 = kBT / ke

# Longitud cuadratica media del enlace
b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

# Formula teorica completa para N finito
N_suave = np.linspace(4, 64, 500)

Rg2_teor = (
    (N_suave**2 - 1)
    / (6 * N_suave)
    * b2
)

# Pendiente asintotica teorica
m_teor = b2 / 6

# ============================================================
# 4. AJUSTE LINEAL
# <Rg²> = m*N + c
# ============================================================

coef_lineal, cov_lineal = np.polyfit(
    N_fit,
    Rg2_fit,
    1,
    cov=True
)

m, c = coef_lineal

error_m = np.sqrt(cov_lineal[0, 0])
error_c = np.sqrt(cov_lineal[1, 1])

Rg2_pred = m*N_fit + c

SS_res = np.sum((Rg2_fit - Rg2_pred)**2)
SS_tot = np.sum((Rg2_fit - np.mean(Rg2_fit))**2)

R2_lineal = 1 - SS_res/SS_tot

print("\n========== AJUSTE LINEAL ==========")

print(f"Pendiente = {m:.6f} +/- {error_m:.6f}")
print(f"Ordenada  = {c:.6f} +/- {error_c:.6f}")
print(f"R²       = {R2_lineal:.6f}")

print(f"Pendiente teorica = {m_teor:.6f}")

# ============================================================
# 5. GRAFICA DEL AJUSTE LINEAL
# ============================================================

fig, ax = plt.subplots(figsize=(9, 6))

# Puntos de simulacion SIN barras de error
ax.plot(
    valores_N,
    Rg2_medios,
    "o",
    color="tab:blue",
    markersize=6,
    label="Simulación"
)

# Ajuste lineal
ax.plot(
    N_suave,
    m*N_suave+c,
    color="tab:red",
    linewidth=1.7,
    label=(
        rf"Ajuste: $\langle R_g^2\rangle = "
        rf"({m:.4f}\pm{error_m:.4f})N"
        rf"{c:+.4f}$"
        "\n"
        rf"$R^2={R2_lineal:.5f}$"
    )
)

# Prediccion teorica
ax.plot(
    N_suave,
    Rg2_teor,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label=(
        r"Cadena ideal: "
        r"$\langle R_g^2\rangle="
        r"\frac{N^2-1}{6N}b^2$"
    )
)

ax.set_title(
    rf"Radio de giro frente a $N$ ($k_e={ke}$)",
    fontsize=14
)

ax.set_xlabel("Número de monómeros, N")
ax.set_ylabel(r"$\langle R_g^2\rangle$")

ax.grid(alpha=0.3)
ax.legend(loc="upper left", fontsize=9)

ax.set_xlim(2, 67)

fig.tight_layout()

fig.savefig(
    salida_lineal,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close(fig)

# ============================================================
# 6. AJUSTE LOG-LOG
#
# <Rg²> = A*N^(2*nu)
#
# log(<Rg²>) = log(A) + 2*nu*log(N)
# ============================================================

log_N = np.log(N_fit)
log_Rg2 = np.log(Rg2_fit)

coef_log, cov_log = np.polyfit(
    log_N,
    log_Rg2,
    1,
    cov=True
)

alpha, log_A = coef_log

error_alpha = np.sqrt(cov_log[0, 0])

# Exponente de escala
nu = alpha / 2
error_nu = error_alpha / 2

# Prefactor
A = np.exp(log_A)

# Coeficiente de determinacion en escala logaritmica
log_pred = alpha*log_N + log_A

SS_res_log = np.sum((log_Rg2-log_pred)**2)
SS_tot_log = np.sum(
    (log_Rg2-np.mean(log_Rg2))**2
)

R2_log = 1 - SS_res_log/SS_tot_log

print("\n========== AJUSTE LOG-LOG ==========")

print(f"Pendiente 2*nu = {alpha:.6f}")
print(f"Exponente nu   = {nu:.6f} +/- {error_nu:.6f}")
print(f"Prefactor A    = {A:.6f}")
print(f"R² log-log     = {R2_log:.6f}")

print("\nValor teorico esperado: nu = 0.5")

# ============================================================
# 7. GRAFICA LOG-LOG
# ============================================================

fig, ax = plt.subplots(figsize=(9, 6))

# Puntos de simulacion
ax.plot(
    valores_N,
    Rg2_medios,
    "o",
    color="tab:blue",
    markersize=6,
    label="Simulación"
)

# Ajuste de ley de potencia
ax.plot(
    N_suave,
    A*N_suave**alpha,
    color="tab:red",
    linewidth=1.7,
    label=(
        rf"Ajuste: $\langle R_g^2\rangle"
        rf"={A:.4f}N^{{{alpha:.4f}}}$"
        "\n"
        rf"$\nu={nu:.4f}\pm{error_nu:.4f}$"
        "\n"
        rf"$R^2={R2_log:.5f}$"
    )
)

# Prediccion teorica
ax.plot(
    N_suave,
    Rg2_teor,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label="Cadena ideal"
)

# Escalas logaritmicas
ax.set_xscale("log")
ax.set_yscale("log")

ax.set_title(
    rf"Ajuste log-log del radio de giro ($k_e={ke}$)",
    fontsize=14
)

ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$\langle R_g^2\rangle$")

ax.grid(
    True,
    which="both",
    alpha=0.3
)

ax.legend(loc="upper left", fontsize=9)

fig.tight_layout()

fig.savefig(
    salida_loglog,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close(fig)

# ============================================================
# 8. RESUMEN FINAL
# ============================================================

print("\n========== RESUMEN ==========")
print(f"Pendiente lineal = {m:.6f}")
print(f"R² lineal        = {R2_lineal:.6f}")
print(f"Exponente nu     = {nu:.6f} +/- {error_nu:.6f}")
print(f"R² log-log       = {R2_log:.6f}")

print("\nGraficas guardadas:")
print(salida_lineal)
print(salida_loglog)



# %%
# ============================================================
# OBJETIVO 3 - RADIO DE GIRO
# AJUSTE LINEAL Y LOG-LOG
# ke = 100
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

valores_N = np.array(
    [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
)

ke = 100
kBT = 1.0
l0 = 1.0

# Valor minimo de N incluido en los ajustes
N_ajuste_min = 4

# Carpetas
carpeta_datos = base / "Py/plots/OBJETIVO_3/Datos/k100/RG"
carpeta_salida = base / "Py/plots/OBJETIVO_3/RG/k100"

carpeta_salida.mkdir(parents=True, exist_ok=True)

salida_lineal = (
    carpeta_salida /
    "Rg2_vs_N_ke100_lineal_sin_errores.png"
)

salida_loglog = (
    carpeta_salida /
    "Rg2_vs_N_ke100_loglog_sin_errores.png"
)

# Proteger graficas anteriores
for salida in [salida_lineal, salida_loglog]:
    if salida.exists():
        raise FileExistsError(
            f"El archivo ya existe: {salida}"
        )

# ============================================================
# 2. LEER LOS DATOS DE LAS 13 SIMULACIONES
# ============================================================

Rg2_medios = []

for N in valores_N:

    if N == 64:
        nombre = "Rg_Ree_N64_ke100_T200000.txt"
    else:
        nombre = f"Rg_Ree_N{N}_ke100.txt"

    archivo = carpeta_datos / nombre

    if not archivo.is_file():
        raise FileNotFoundError(
            f"No se encuentra: {archivo}"
        )

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 5:
        raise ValueError(
            f"El archivo N={N} no tiene cinco columnas."
        )

    # Media acumulada final de Rg²
    media = datos[-1, 1]

    if not np.isfinite(media) or media <= 0:
        raise ValueError(
            f"Valor invalido de Rg² para N={N}"
        )

    Rg2_medios.append(media)

    print(
        f"N={N:2d}  |  "
        f"Tiempo={datos[-1,0]:.0f}  |  "
        f"<Rg²>={media:.6f}"
    )

Rg2_medios = np.array(Rg2_medios)

# Puntos que entran en los ajustes
mascara_ajuste = valores_N >= N_ajuste_min

N_fit = valores_N[mascara_ajuste]
Rg2_fit = Rg2_medios[mascara_ajuste]

# ============================================================
# 3. PREDICCION TEORICA DE LA CADENA IDEAL
# ============================================================

sigma2 = kBT / ke

# Longitud cuadratica media del enlace
b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

# Formula teorica completa para N finito
N_suave = np.linspace(4, 64, 500)

Rg2_teor = (
    (N_suave**2 - 1)
    / (6 * N_suave)
    * b2
)

# Pendiente asintotica teorica
m_teor = b2 / 6

# ============================================================
# 4. AJUSTE LINEAL
# <Rg²> = m*N + c
# ============================================================

coef_lineal, cov_lineal = np.polyfit(
    N_fit,
    Rg2_fit,
    1,
    cov=True
)

m, c = coef_lineal

error_m = np.sqrt(cov_lineal[0, 0])
error_c = np.sqrt(cov_lineal[1, 1])

Rg2_pred = m*N_fit + c

SS_res = np.sum((Rg2_fit - Rg2_pred)**2)
SS_tot = np.sum((Rg2_fit - np.mean(Rg2_fit))**2)

R2_lineal = 1 - SS_res/SS_tot

print("\n========== AJUSTE LINEAL ==========")

print(f"Pendiente = {m:.6f} +/- {error_m:.6f}")
print(f"Ordenada  = {c:.6f} +/- {error_c:.6f}")
print(f"R²       = {R2_lineal:.6f}")

print(f"Pendiente teorica = {m_teor:.6f}")

# ============================================================
# 5. GRAFICA DEL AJUSTE LINEAL
# ============================================================

fig, ax = plt.subplots(figsize=(9, 6))

# Puntos de simulacion SIN barras de error
ax.plot(
    valores_N,
    Rg2_medios,
    "o",
    color="tab:blue",
    markersize=6,
    label="Simulación"
)

# Ajuste lineal
ax.plot(
    N_suave,
    m*N_suave+c,
    color="tab:red",
    linewidth=1.7,
    label=(
        rf"Ajuste: $\langle R_g^2\rangle = "
        rf"({m:.4f}\pm{error_m:.4f})N"
        rf"{c:+.4f}$"
        "\n"
        rf"$R^2={R2_lineal:.5f}$"
    )
)

# Prediccion teorica
ax.plot(
    N_suave,
    Rg2_teor,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label=(
        r"Cadena ideal: "
        r"$\langle R_g^2\rangle="
        r"\frac{N^2-1}{6N}b^2$"
    )
)

ax.set_title(
    rf"Radio de giro frente a $N$ ($k_e={ke}$)",
    fontsize=14
)

ax.set_xlabel("Número de monómeros, N")
ax.set_ylabel(r"$\langle R_g^2\rangle$")

ax.grid(alpha=0.3)
ax.legend(loc="upper left", fontsize=9)

ax.set_xlim(2, 67)

fig.tight_layout()

fig.savefig(
    salida_lineal,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close(fig)

# ============================================================
# 6. AJUSTE LOG-LOG
#
# <Rg²> = A*N^(2*nu)
#
# log(<Rg²>) = log(A) + 2*nu*log(N)
# ============================================================

log_N = np.log(N_fit)
log_Rg2 = np.log(Rg2_fit)

coef_log, cov_log = np.polyfit(
    log_N,
    log_Rg2,
    1,
    cov=True
)

alpha, log_A = coef_log

error_alpha = np.sqrt(cov_log[0, 0])

# Exponente de escala
nu = alpha / 2
error_nu = error_alpha / 2

# Prefactor
A = np.exp(log_A)

# Coeficiente de determinacion en escala logaritmica
log_pred = alpha*log_N + log_A

SS_res_log = np.sum((log_Rg2-log_pred)**2)
SS_tot_log = np.sum(
    (log_Rg2-np.mean(log_Rg2))**2
)

R2_log = 1 - SS_res_log/SS_tot_log

print("\n========== AJUSTE LOG-LOG ==========")

print(f"Pendiente 2*nu = {alpha:.6f}")
print(f"Exponente nu   = {nu:.6f} +/- {error_nu:.6f}")
print(f"Prefactor A    = {A:.6f}")
print(f"R² log-log     = {R2_log:.6f}")

print("\nValor teorico esperado: nu = 0.5")

# ============================================================
# 7. GRAFICA LOG-LOG
# ============================================================

fig, ax = plt.subplots(figsize=(9, 6))

# Puntos de simulacion
ax.plot(
    valores_N,
    Rg2_medios,
    "o",
    color="tab:blue",
    markersize=6,
    label="Simulación"
)

# Ajuste de ley de potencia
ax.plot(
    N_suave,
    A*N_suave**alpha,
    color="tab:red",
    linewidth=1.7,
    label=(
        rf"Ajuste: $\langle R_g^2\rangle"
        rf"={A:.4f}N^{{{alpha:.4f}}}$"
        "\n"
        rf"$\nu={nu:.4f}\pm{error_nu:.4f}$"
        "\n"
        rf"$R^2={R2_log:.5f}$"
    )
)

# Prediccion teorica
ax.plot(
    N_suave,
    Rg2_teor,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label="Cadena ideal"
)

# Escalas logaritmicas
ax.set_xscale("log")
ax.set_yscale("log")

ax.set_title(
    rf"Ajuste log-log del radio de giro ($k_e={ke}$)",
    fontsize=14
)

ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$\langle R_g^2\rangle$")

ax.grid(
    True,
    which="both",
    alpha=0.3
)

ax.legend(loc="upper left", fontsize=9)

fig.tight_layout()

fig.savefig(
    salida_loglog,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close(fig)

# ============================================================
# 8. RESUMEN FINAL
# ============================================================

print("\n========== RESUMEN ==========")
print(f"Pendiente lineal = {m:.6f}")
print(f"R² lineal        = {R2_lineal:.6f}")
print(f"Exponente nu     = {nu:.6f} +/- {error_nu:.6f}")
print(f"R² log-log       = {R2_log:.6f}")

print("\nGraficas guardadas:")
print(salida_lineal)
print(salida_loglog)


# %%

# %%
# ============================================================
# OBJETIVO 3 - DISTANCIA ENTRE EXTREMOS
# AJUSTE LINEAL Y LOG-LOG
# ke = 100
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. PARAMETROS Y CARPETAS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

valores_N = np.array(
    [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
)

ke = 100
kBT = 1.0
l0 = 1.0

carpeta_datos = base / "Py/plots/OBJETIVO_3/Datos/k100/RG"
carpeta_salida = base / "Py/plots/OBJETIVO_3/RG/k100"

carpeta_salida.mkdir(parents=True, exist_ok=True)

salida_lineal = (
    carpeta_salida / "Ree2_vs_N_ke100_lineal_Claudia.png"
)

salida_loglog = (
    carpeta_salida / "Ree2_vs_N_ke100_loglog_Claudia.png"
)

# ============================================================
# 2. CARGAR LOS DATOS
# ============================================================

Ree2_medios = []
Rg2_medios = []

for N in valores_N:

    if N == 64:
        nombre = "Rg_Ree_N64_ke100_T200000.txt"
        tiempo_esperado = 200000
    else:
        nombre = f"Rg_Ree_N{N}_ke100.txt"
        tiempo_esperado = 100000

    archivo = carpeta_datos / nombre

    if not archivo.is_file():
        raise FileNotFoundError(
            f"No se encuentra: {archivo}"
        )

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 5:
        raise ValueError(
            f"N={N}: el TXT no tiene cinco columnas."
        )

    if not np.isclose(datos[-1, 0], tiempo_esperado):
        raise ValueError(
            f"N={N}: el tiempo final no coincide "
            f"con el esperado."
        )

    # Columna 2: radio de giro al cuadrado
    # Columna 3: distancia entre extremos al cuadrado

    Rg2_medios.append(datos[-1, 1])
    Ree2_medios.append(datos[-1, 2])

    print(
        f"N={N:2d} | "
        f"<Ree²>={datos[-1, 2]:.6f}"
    )

Rg2_medios = np.array(Rg2_medios)
Ree2_medios = np.array(Ree2_medios)

# ============================================================
# 3. VALORES TEORICOS
# ============================================================

sigma2 = kBT / ke

b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

N_suave = np.linspace(4, 64, 500)

Ree2_teor = (N_suave - 1) * b2

print("\nLongitud cuadratica media del enlace:")
print(f"b² = {b2:.6f}")

# ============================================================
# 4. AJUSTE LINEAL
#
# <Ree²> = m*N + c
# Teoria: m = b², c = -b²
# ============================================================

coef, cov = np.polyfit(
    valores_N,
    Ree2_medios,
    1,
    cov=True
)

m, c = coef

error_m = np.sqrt(cov[0, 0])
error_c = np.sqrt(cov[1, 1])

Ree2_ajuste = m*valores_N + c

SS_res = np.sum(
    (Ree2_medios - Ree2_ajuste)**2
)

SS_tot = np.sum(
    (Ree2_medios - np.mean(Ree2_medios))**2
)

R2_lineal = 1 - SS_res/SS_tot

print("\n========== AJUSTE LINEAL ==========")

print(f"Pendiente = {m:.6f} +/- {error_m:.6f}")
print(f"Ordenada  = {c:.6f} +/- {error_c:.6f}")
print(f"R²       = {R2_lineal:.6f}")

print("\nPrediccion teorica:")
print(f"Pendiente = {b2:.6f}")
print(f"Ordenada  = {-b2:.6f}")

# ============================================================
# 5. GRAFICA LINEAL
# ============================================================

fig, ax = plt.subplots(figsize=(9, 6))

# Puntos de simulacion
ax.plot(
    valores_N,
    Ree2_medios,
    "o",
    color="tab:blue",
    markersize=6,
    label="Simulación"
)

# Ajuste
ax.plot(
    N_suave,
    m*N_suave + c,
    color="tab:red",
    linewidth=1.8,
    label=(
        rf"Ajuste: $\langle R_{{ee}}^2\rangle="
        rf"({m:.4f}\pm{error_m:.4f})N{c:+.4f}$"
        "\n"
        rf"$R^2={R2_lineal:.5f}$"
    )
)

# Prediccion teorica
ax.plot(
    N_suave,
    Ree2_teor,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label=(
        r"Cadena ideal: "
        r"$\langle R_{ee}^2\rangle=(N-1)b^2$"
    )
)

ax.set_title(
    rf"Distancia entre extremos al cuadrado ($k_e={ke}$)",
    fontsize=14
)

ax.set_xlabel("Número de monómeros, N")
ax.set_ylabel(r"$\langle R_{ee}^2\rangle$")

ax.grid(alpha=0.3)
ax.legend(loc="upper left", fontsize=9)
ax.set_xlim(2, 67)

fig.tight_layout()

if salida_lineal.exists():
    print("La grafica lineal ya existe; no se sobrescribe.")
else:
    fig.savefig(
        salida_lineal,
        dpi=300,
        bbox_inches="tight"
    )

plt.show()
plt.close(fig)

# ============================================================
# 6. AJUSTE LOG-LOG
#
# <Ree²> = A*N^(2*nu)
# ============================================================

log_N = np.log(valores_N)
log_Ree2 = np.log(Ree2_medios)

coef_log, cov_log = np.polyfit(
    log_N,
    log_Ree2,
    1,
    cov=True
)

alpha, log_A = coef_log

error_alpha = np.sqrt(cov_log[0, 0])

nu = alpha / 2
error_nu = error_alpha / 2

A = np.exp(log_A)

log_ajuste = alpha*log_N + log_A

SS_res_log = np.sum(
    (log_Ree2 - log_ajuste)**2
)

SS_tot_log = np.sum(
    (log_Ree2 - np.mean(log_Ree2))**2
)

R2_log = 1 - SS_res_log/SS_tot_log

# Pendiente aparente de la teoria finita
alpha_teor_finito = np.polyfit(
    log_N,
    np.log((valores_N - 1)*b2),
    1
)[0]

print("\n========== AJUSTE LOG-LOG ==========")

print(f"Pendiente 2*nu = {alpha:.6f}")
print(f"Exponente nu   = {nu:.6f} +/- {error_nu:.6f}")
print(f"Prefactor A    = {A:.6f}")
print(f"R² log-log     = {R2_log:.6f}")

print(
    f"nu aparente de la teoria finita = "
    f"{alpha_teor_finito/2:.6f}"
)

print("nu asintotico teorico = 0.500000")

# ============================================================
# 7. GRAFICA LOG-LOG
# ============================================================

fig, ax = plt.subplots(figsize=(9, 6))

# Puntos
ax.plot(
    valores_N,
    Ree2_medios,
    "o",
    color="tab:blue",
    markersize=6,
    label="Simulación"
)

# Ley de potencia ajustada
ax.plot(
    N_suave,
    A*N_suave**alpha,
    color="tab:red",
    linewidth=1.8,
    label=(
        rf"Ajuste: $\langle R_{{ee}}^2\rangle"
        rf"={A:.4f}N^{{{alpha:.4f}}}$"
        "\n"
        rf"$\nu={nu:.4f}\pm{error_nu:.4f}$"
        "\n"
        rf"$R^2={R2_log:.5f}$"
    )
)

# Teoria exacta para N finito
ax.plot(
    N_suave,
    Ree2_teor,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label="Cadena ideal"
)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_title(
    rf"Ajuste log-log de la distancia entre extremos ($k_e={ke}$)",
    fontsize=14
)

ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$\langle R_{ee}^2\rangle$")

ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="upper left", fontsize=9)

fig.tight_layout()

if salida_loglog.exists():
    print("La grafica log-log ya existe; no se sobrescribe.")
else:
    fig.savefig(
        salida_loglog,
        dpi=300,
        bbox_inches="tight"
    )

plt.show()
plt.close(fig)

# ============================================================
# 8. COMPROBACION DE LA RELACION CON EL RADIO DE GIRO
# ============================================================

print("\n========== RELACION Ree² / Rg² ==========")

for N, ree2, rg2 in zip(
    valores_N,
    Ree2_medios,
    Rg2_medios
):

    razon_sim = ree2 / rg2
    razon_teor = 6*N/(N+1)

    print(
        f"N={N:2d} | "
        f"Simulacion={razon_sim:.4f} | "
        f"Teoria={razon_teor:.4f}"
    )

print("\nGraficas:")
print(salida_lineal)
print(salida_loglog)

# %%

# %%
# ============================================================
# OBJETIVO 3 - RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS
# ke = 1000
# EVOLUCION TEMPORAL PARA TODOS LOS N
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

valores_N = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]

ke = 1000
kBT = 1.0
l0 = 1.0

# Configuracion de las graficas
margen_relativo = 0.30

# Cambiar a True para representar tambien las medias por bloques
mostrar_bloques = False

# Tiempo inicial que se muestra en el eje horizontal
t_inicio_grafica = 0

# ============================================================
# 2. CARPETAS
# ============================================================

carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k1000/RG"
)

carpeta_salida = (
    base / "Py/plots/OBJETIVO_3/RG/k1000"
)

carpeta_salida.mkdir(parents=True, exist_ok=True)

# ============================================================
# 3. VALORES TEORICOS
# ============================================================

# Fluctuaciones de longitud de los enlaces
sigma2 = kBT / ke

# Longitud cuadratica media del enlace (aproximacion 3D)
b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

print(f"Longitud cuadratica media del enlace: b² = {b2:.6f}")

# ============================================================
# 4. RECORRER TODOS LOS TAMANOS
# ============================================================

for N in valores_N:

    archivo = (
        carpeta_datos /
        f"Rg_Ree_N{N}_ke{ke}.txt"
    )

    salida = (
        carpeta_salida /
        f"Rg_Ree_N{N}_ke{ke}_evolucion.png"
    )

    # Comprobar existencia del TXT
    if not archivo.is_file():
        print(f"No se encuentra el archivo de N={N}")
        continue

    # Evitar sobrescribir graficas anteriores
    if salida.exists():
        print(f"La grafica de N={N} ya existe. Se omite.")
        continue

    # ========================================================
    # 5. LEER LOS DATOS
    # ========================================================

    datos = np.loadtxt(
        archivo,
        comments="#",
        ndmin=2
    )

    if datos.shape[1] != 5:
        raise ValueError(
            f"N={N}: se esperaban 5 columnas, "
            f"pero hay {datos.shape[1]}."
        )

    t = datos[:, 0]

    Rg2_m = datos[:, 1]
    Ree2_m = datos[:, 2]

    Rg2_bloques = datos[:, 3]
    Ree2_bloques = datos[:, 4]

    # Tiempos de medicion utilizados en Fortran
    if N == 64:
        tiempo_esperado = 200000
    elif N >= 40:
        tiempo_esperado = 150000
    else:
        tiempo_esperado = 100000

    if not np.isclose(t[-1], tiempo_esperado):
        print(
            f"N={N}: archivo posiblemente incompleto. "
            f"Tiempo final = {t[-1]:.0f}, "
            f"esperado = {tiempo_esperado}."
        )
        continue

    # ========================================================
    # 6. VALORES TEORICOS PARA ESTE N
    # ========================================================

    Rg2_teor = (N**2 - 1) * b2 / (6*N)

    Ree2_teor = (N - 1) * b2

    # ========================================================
    # 7. MOSTRAR LOS RESULTADOS NUMERICOS
    # ========================================================

    print(f"\n{'=' * 48}")
    print(f"N = {N}, ke = {ke}")
    print(f"Tiempo de medicion = {t[-1]:.0f}")
    print(f"{'=' * 48}")

    print(f"<Rg²> simulacion = {Rg2_m[-1]:.6f}")
    print(f"<Rg²> teorico    = {Rg2_teor:.6f}")

    print(f"<Ree²> simulacion = {Ree2_m[-1]:.6f}")
    print(f"<Ree²> teorico    = {Ree2_teor:.6f}")

    # Diferencias relativas
    error_rg = 100 * (Rg2_m[-1] - Rg2_teor) / Rg2_teor

    error_ree = 100 * (Ree2_m[-1] - Ree2_teor) / Ree2_teor

    print(f"Desviacion Rg²  = {error_rg:+.3f} %")
    print(f"Desviacion Ree² = {error_ree:+.3f} %")

    # ========================================================
    # 8. SELECCIONAR EL INTERVALO REPRESENTADO
    # ========================================================

    mascara = t >= t_inicio_grafica

    t_graf = t[mascara]

    # ========================================================
    # 9. CREAR FIGURA
    # ========================================================

    fig, axs = plt.subplots(
        1, 2,
        figsize=(15, 5.5)
    )

    series = [
        (
            Rg2_m[mascara],
            Rg2_bloques[mascara],
            Rg2_teor,
            "Radio de giro al cuadrado",
            r"$\langle R_g^2\rangle_t$",
            "tab:blue",
            (
                r"Cadena ideal: "
                rf"$\frac{{N^2-1}}{{6N}}b^2"
                rf" = {Rg2_teor:.3f}$"
            )
        ),
        (
            Ree2_m[mascara],
            Ree2_bloques[mascara],
            Ree2_teor,
            "Distancia entre extremos al cuadrado",
            r"$\langle R_{ee}^2\rangle_t$",
            "tab:green",
            (
                r"Cadena ideal: "
                rf"$(N-1)b^2 = {Ree2_teor:.3f}$"
            )
        )
    ]

    for ax, (
        media, bloques, teor, titulo,
        ylabel, color, etiqueta
    ) in zip(axs, series):

        # Promedio acumulado
        ax.plot(
            t_graf,
            media,
            color=color,
            linewidth=1.5,
            label="Promedio acumulado"
        )

        # Medias por bloques, si se desea
        if mostrar_bloques:
            ax.plot(
                t_graf,
                bloques,
                color=color,
                alpha=0.30,
                linewidth=1.0,
                label="Promedio por bloques"
            )

        # Valor teorico
        ax.axhline(
            teor,
            color="black",
            linestyle="--",
            linewidth=1.2,
            label=etiqueta
        )

        # Etiquetas
        ax.set_title(titulo)

        ax.set_xlabel("Tiempo de medición")
        ax.set_ylabel(ylabel)

        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=9)

        ax.set_xlim(t_graf[0], t_graf[-1])

        # ====================================================
        # ESCALA VERTICAL CON MARGEN
        # ====================================================

        if mostrar_bloques:
            minimo = min(np.min(media), np.min(bloques), teor)
            maximo = max(np.max(media), np.max(bloques), teor)
        else:
            minimo = min(np.min(media), teor)
            maximo = max(np.max(media), teor)

        margen = max(
            margen_relativo * (maximo - minimo),
            0.05 * teor
        )

        ax.set_ylim(
            max(0, minimo - margen),
            maximo + margen
        )

    # ========================================================
    # 10. TITULO Y GUARDADO
    # ========================================================

    fig.suptitle(
        rf"Evolución del polímero ($N={N}$, $k_e={ke}$)",
        fontsize=15
    )

    fig.tight_layout()

    fig.savefig(
        salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print(f"Grafica guardada en: {salida}")

print("\nTodas las graficas disponibles se han generado.")

# %%



# %%
# ============================================================
# OBJETIVO 3 - COMPROBACION DE CONVERGENCIA
# ke = 1000, N = 32, 48, 64
# RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS
# ============================================================

import numpy as np
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta = base / "Py/plots/OBJETIVO_3/Datos/k1000/RG"

valores_N = [64]

ke = 1000
kBT = 1.0
l0 = 1.0

# Longitud cuadratica media del enlace
sigma2 = kBT / ke

b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

# ============================================================
# ANALIZAR CADA CADENA
# ============================================================

for N in valores_N:

    archivo = carpeta / f"Rg_Ree_N{N}_ke{ke}.txt"

    if not archivo.is_file():
        print(f"No se encuentra el TXT de N={N}")
        continue

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 5:
        raise ValueError(
            f"N={N}: se esperaban cinco columnas."
        )

    t = datos[:, 0]

    # Promedios acumulados
    Rg2_m = datos[:, 1]
    Ree2_m = datos[:, 2]

    # Medias por bloques guardadas por Fortran
    Rg2_bloques = datos[:, 3]
    Ree2_bloques = datos[:, 4]

    # Valores teoricos
    Rg2_teor = (N**2 - 1) * b2 / (6*N)
    Ree2_teor = (N - 1) * b2

    # Dividir la medicion en cuatro partes
    indices = np.array_split(np.arange(len(t)), 4)

    medias_rg = []
    medias_ree = []

    print(f"\n{'=' * 58}")
    print(f"N = {N}, ke = {ke}")
    print(f"Tiempo de medicion = {t[-1]:.0f}")
    print(f"{'=' * 58}")

    for i, idx in enumerate(indices, 1):

        media_rg = np.mean(Rg2_bloques[idx])
        media_ree = np.mean(Ree2_bloques[idx])

        medias_rg.append(media_rg)
        medias_ree.append(media_ree)

        print(
            f"Bloque {i}: "
            f"<Rg²> = {media_rg:.5f} | "
            f"<Ree²> = {media_ree:.5f}"
        )

    # ========================================================
    # COMPARACION ENTRE PRIMER Y ULTIMO BLOQUE
    # ========================================================

    diferencia_rg = (
        100 * (medias_rg[-1] - medias_rg[0])
        / medias_rg[0]
    )

    diferencia_ree = (
        100 * (medias_ree[-1] - medias_ree[0])
        / medias_ree[0]
    )

    print("\nCAMBIO ENTRE BLOQUE 1 Y BLOQUE 4")

    print(f"Rg²:  {diferencia_rg:+.2f} %")
    print(f"Ree²: {diferencia_ree:+.2f} %")

    # ========================================================
    # COMPARACION CON LA TEORIA
    # ========================================================

    print("\nRESULTADOS FINALES")

    print(
        f"<Rg²> = {Rg2_m[-1]:.5f} "
        f"(teoria: {Rg2_teor:.5f})"
    )

    print(
        f"<Ree²> = {Ree2_m[-1]:.5f} "
        f"(teoria: {Ree2_teor:.5f})"
    )

    error_rg = (
        100 * (Rg2_m[-1] - Rg2_teor) / Rg2_teor
    )

    error_ree = (
        100 * (Ree2_m[-1] - Ree2_teor) / Ree2_teor
    )

    print("\nDESVIACION RESPECTO A LA TEORIA")

    print(f"Rg²:  {error_rg:+.2f} %")
    print(f"Ree²: {error_ree:+.2f} %")

print("\nComprobacion terminada.")

# %%




# %%
# ============================================================
# OBJETIVO 3 - AJUSTES DEL TAMANO DEL POLIMERO
# ke = 1000
#
# 1. Rg² frente a N: lineal
# 2. Rg² frente a N: log-log
# 3. Ree² frente a N: lineal
# 4. Ree² frente a N: log-log
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. PARAMETROS Y CARPETAS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

valores_N = np.array(
    [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
)

ke = 1000
kBT = 1.0
l0 = 1.0

carpeta_datos = base / "Py/plots/OBJETIVO_3/Datos/k1000/RG"
carpeta_salida = base / "Py/plots/OBJETIVO_3/RG/k1000"

carpeta_salida.mkdir(parents=True, exist_ok=True)

# Nombres nuevos para proteger otros resultados
nombres_graficas = {
    "rg_lineal": "Rg2_vs_N_ke1000_lineal_Claudia.png",
    "rg_loglog": "Rg2_vs_N_ke1000_loglog_Claudia.png",
    "ree_lineal": "Ree2_vs_N_ke1000_lineal_Claudia.png",
    "ree_loglog": "Ree2_vs_N_ke1000_loglog_Claudia.png"
}

salida_csv = (
    carpeta_salida / "Resumen_Rg_Ree_ke1000_Claudia.csv"
)

# Evitar sobrescribir resultados
salidas = [
    carpeta_salida / nombre
    for nombre in nombres_graficas.values()
] + [salida_csv]

for salida in salidas:
    if salida.exists():
        raise FileExistsError(
            f"El archivo ya existe: {salida}. "
            "Cambia su nombre antes de ejecutar."
        )

# ============================================================
# 2. LEER LAS 13 SIMULACIONES
# ============================================================

Rg2_medios = []
Ree2_medios = []
tiempos = []

for N in valores_N:

    archivo = carpeta_datos / f"Rg_Ree_N{N}_ke{ke}.txt"

    if not archivo.is_file():
        raise FileNotFoundError(
            f"No se encuentra: {archivo}"
        )

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 5:
        raise ValueError(
            f"N={N}: el TXT debe contener cinco columnas."
        )

    # Duracion esperada de la fase de medicion
    if N == 64:
        t_esperado = 200000
    elif N >= 40:
        t_esperado = 150000
    else:
        t_esperado = 100000

    if not np.isclose(
        datos[-1, 0], t_esperado,
        rtol=0, atol=1
    ):
        raise ValueError(
            f"N={N}: tiempo final inesperado."
        )

    rg2 = datos[-1, 1]
    ree2 = datos[-1, 2]

    if not np.isfinite(rg2) or rg2 <= 0:
        raise ValueError(f"Rg² invalido para N={N}")

    if not np.isfinite(ree2) or ree2 <= 0:
        raise ValueError(f"Ree² invalido para N={N}")

    Rg2_medios.append(rg2)
    Ree2_medios.append(ree2)
    tiempos.append(datos[-1, 0])

    print(
        f"N={N:2d} | "
        f"<Rg²>={rg2:.6f} | "
        f"<Ree²>={ree2:.6f}"
    )

Rg2_medios = np.array(Rg2_medios)
Ree2_medios = np.array(Ree2_medios)
tiempos = np.array(tiempos)

# ============================================================
# 3. VALORES TEORICOS
# ============================================================

sigma2 = kBT / ke

# Longitud cuadratica media de los enlaces
b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

N_suave = np.linspace(4, 64, 500)

# Radio de giro
Rg2_teor = (
    (N_suave**2 - 1)
    / (6 * N_suave)
    * b2
)

# Distancia entre extremos
Ree2_teor = (N_suave - 1) * b2

Rg2_teor_puntos = (
    (valores_N**2 - 1)
    / (6 * valores_N)
    * b2
)

Ree2_teor_puntos = (valores_N - 1) * b2

print("\n========== TEORIA ==========")
print(f"b² = {b2:.6f}")
print(f"Pendiente Rg² teorica: {b2/6:.6f}")
print(f"Pendiente Ree² teorica: {b2:.6f}")

# ============================================================
# 4. FUNCIONES DE AJUSTE
# ============================================================

def ajuste_lineal(N, y):

    coef, cov = np.polyfit(
        N, y, 1, cov=True
    )

    m, c = coef

    error_m = np.sqrt(cov[0, 0])
    error_c = np.sqrt(cov[1, 1])

    pred = m*N + c

    R2 = 1 - (
        np.sum((y-pred)**2)
        / np.sum((y-np.mean(y))**2)
    )

    return m, c, error_m, error_c, R2


def ajuste_loglog(N, y):

    log_N = np.log(N)
    log_y = np.log(y)

    coef, cov = np.polyfit(
        log_N, log_y, 1, cov=True
    )

    alpha, log_A = coef

    error_alpha = np.sqrt(cov[0, 0])

    A = np.exp(log_A)

    nu = alpha/2
    error_nu = error_alpha/2

    pred = alpha*log_N + log_A

    R2 = 1 - (
        np.sum((log_y-pred)**2)
        / np.sum((log_y-np.mean(log_y))**2)
    )

    return A, alpha, nu, error_nu, R2


# ============================================================
# 5. AJUSTES DEL RADIO DE GIRO
# ============================================================

m_rg, c_rg, em_rg, ec_rg, R2_rg = ajuste_lineal(
    valores_N, Rg2_medios
)

A_rg, alpha_rg, nu_rg, enu_rg, R2log_rg = ajuste_loglog(
    valores_N, Rg2_medios
)

# ============================================================
# 6. AJUSTES DE LA DISTANCIA ENTRE EXTREMOS
# ============================================================

m_ree, c_ree, em_ree, ec_ree, R2_ree = ajuste_lineal(
    valores_N, Ree2_medios
)

A_ree, alpha_ree, nu_ree, enu_ree, R2log_ree = ajuste_loglog(
    valores_N, Ree2_medios
)

# ============================================================
# 7. FUNCION PARA CREAR LAS GRAFICAS
# ============================================================

def crear_grafica(
    y_datos, y_teor, observable,
    titulo, nombre_archivo,
    tipo, ajuste
):

    fig, ax = plt.subplots(figsize=(9, 6))

    # Puntos simulados, sin barras de error
    ax.plot(
        valores_N,
        y_datos,
        "o",
        color="tab:blue",
        markersize=6,
        label="Simulación"
    )

    if tipo == "lineal":

        m, c, error_m, error_c, R2 = ajuste

        y_fit = m*N_suave + c

        etiqueta = (
            rf"Ajuste: ${observable}="
            rf"({m:.4f}\pm{error_m:.4f})N"
            rf"{c:+.4f}$"
            "\n"
            rf"$R^2={R2:.5f}$"
        )

    else:

        A, alpha, nu, error_nu, R2 = ajuste

        y_fit = A*N_suave**alpha

        etiqueta = (
            rf"Ajuste: ${observable}="
            rf"{A:.4f}N^{{{alpha:.4f}}}$"
            "\n"
            rf"$\nu={nu:.4f}\pm{error_nu:.4f}$"
            "\n"
            rf"$R^2={R2:.5f}$"
        )

    # Ajuste experimental
    ax.plot(
        N_suave,
        y_fit,
        color="tab:red",
        linewidth=1.8,
        label=etiqueta
    )

    # Prediccion teorica
    ax.plot(
        N_suave,
        y_teor,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label="Cadena ideal"
    )

    # Escala logaritmica si corresponde
    if tipo == "loglog":
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
    else:
        ax.set_xlim(2, 67)
        ax.grid(alpha=0.3)

    ax.set_title(titulo, fontsize=14)
    ax.set_xlabel(r"Número de monómeros, $N$")
    ax.set_ylabel(f"${observable}$")

    ax.legend(loc="upper left", fontsize=9)

    fig.tight_layout()

    salida = carpeta_salida / nombre_archivo

    fig.savefig(
        salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print(f"Guardada: {salida.name}")


# ============================================================
# 8. GENERAR LOS CUATRO AJUSTES
# ============================================================

# RADIO DE GIRO - LINEAL
crear_grafica(
    Rg2_medios,
    Rg2_teor,
    r"\langle R_g^2\rangle",
    rf"Radio de giro al cuadrado frente a $N$ ($k_e={ke}$)",
    nombres_graficas["rg_lineal"],
    "lineal",
    (m_rg, c_rg, em_rg, ec_rg, R2_rg)
)

# RADIO DE GIRO - LOG-LOG
crear_grafica(
    Rg2_medios,
    Rg2_teor,
    r"\langle R_g^2\rangle",
    rf"Ajuste log-log del radio de giro ($k_e={ke}$)",
    nombres_graficas["rg_loglog"],
    "loglog",
    (A_rg, alpha_rg, nu_rg, enu_rg, R2log_rg)
)

# DISTANCIA ENTRE EXTREMOS - LINEAL
crear_grafica(
    Ree2_medios,
    Ree2_teor,
    r"\langle R_{ee}^2\rangle",
    rf"Distancia entre extremos al cuadrado ($k_e={ke}$)",
    nombres_graficas["ree_lineal"],
    "lineal",
    (m_ree, c_ree, em_ree, ec_ree, R2_ree)
)

# DISTANCIA ENTRE EXTREMOS - LOG-LOG
crear_grafica(
    Ree2_medios,
    Ree2_teor,
    r"\langle R_{ee}^2\rangle",
    rf"Ajuste log-log de la distancia entre extremos ($k_e={ke}$)",
    nombres_graficas["ree_loglog"],
    "loglog",
    (A_ree, alpha_ree, nu_ree, enu_ree, R2log_ree)
)

# ============================================================
# 9. RESUMEN NUMERICO DE LOS AJUSTES
# ============================================================

print("\n========== RADIO DE GIRO ==========")

print(
    f"Lineal: pendiente = {m_rg:.6f} +/- {em_rg:.6f}"
)
print(f"Ordenada = {c_rg:.6f}")
print(f"R² lineal = {R2_rg:.6f}")

print(
    f"Log-log: nu = {nu_rg:.6f} +/- {enu_rg:.6f}"
)
print(f"R² log-log = {R2log_rg:.6f}")

print("\n========== DISTANCIA ENTRE EXTREMOS ==========")

print(
    f"Lineal: pendiente = {m_ree:.6f} +/- {em_ree:.6f}"
)
print(f"Ordenada = {c_ree:.6f}")
print(f"R² lineal = {R2_ree:.6f}")

print(
    f"Log-log: nu = {nu_ree:.6f} +/- {enu_ree:.6f}"
)
print(f"R² log-log = {R2log_ree:.6f}")

# ============================================================
# 10. GUARDAR TABLA RESUMEN
# ============================================================

tabla = np.column_stack([
    valores_N,
    tiempos,
    Rg2_medios,
    Ree2_medios,
    Rg2_teor_puntos,
    Ree2_teor_puntos
])

np.savetxt(
    salida_csv,
    tabla,
    delimiter=",",
    header=(
        "N,tiempo_medicion,Rg2_media,Ree2_media,"
        "Rg2_teorica,Ree2_teorica"
    ),
    comments="",
    fmt=["%d", "%.0f", "%.8f", "%.8f", "%.8f", "%.8f"]
)

print("\nCuatro graficas generadas correctamente.")
print("Resumen CSV guardado en:", salida_csv)

# %%



# %%
# ============================================================
# OBJETIVO 3 - COMPARACION FINAL
# RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS
# ke = 100 frente a ke = 1000
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

valores_N = np.array(
    [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
)

kBT = 1.0
l0 = 1.0

# Carpeta general donde guardaremos las dos graficas
carpeta_salida = base / "Py/plots/OBJETIVO_3/RG"

carpeta_salida.mkdir(parents=True, exist_ok=True)

# ============================================================
# 2. FUNCION PARA LEER LOS RESULTADOS
# ============================================================

def cargar_resultados(ke):

    carpeta_resumen = (
        base / f"Py/plots/OBJETIVO_3/RG/k{ke}"
    )

    resumen = (
        carpeta_resumen /
        f"Resumen_Rg_Ree_ke{ke}_Claudia.csv"
    )

    if resumen.exists():

        datos = np.genfromtxt(
            resumen,
            delimiter=",",
            names=True,
            encoding="utf-8-sig",
            ndmin=1
        )

        columnas = datos.dtype.names

        necesarias = ["N", "Rg2_media", "Ree2_media"]

        if not all(c in columnas for c in necesarias):
            raise ValueError(
                f"El CSV {resumen.name} no tiene las columnas "
                f"{necesarias}. Tiene: {columnas}"
            )

        N = np.asarray(datos["N"], dtype=int)
        rg2 = np.asarray(datos["Rg2_media"], dtype=float)
        ree2 = np.asarray(datos["Ree2_media"], dtype=float)

    else:

        # Si no existe el CSV, leer los TXT
        carpeta_datos = (
            base / f"Py/plots/OBJETIVO_3/Datos/k{ke}/RG"
        )

        N = valores_N.copy()

        rg2 = []
        ree2 = []

        for n in N:

            if ke == 100 and n == 64:
                nombre = "Rg_Ree_N64_ke100_T200000.txt"
            else:
                nombre = f"Rg_Ree_N{n}_ke{ke}.txt"

            archivo = carpeta_datos / nombre

            if not archivo.exists():
                raise FileNotFoundError(
                    f"No encuentro el resumen de k={ke} "
                    f"ni el TXT: {archivo}. "
                    "Necesitamos conservar al menos el CSV "
                    "con Rg2_media y Ree2_media."
                )

            datos = np.loadtxt(archivo, ndmin=2)

            rg2.append(datos[-1, 1])
            ree2.append(datos[-1, 2])

        rg2 = np.array(rg2)
        ree2 = np.array(ree2)

    # Ordenar por N
    orden = np.argsort(N)

    N = N[orden]
    rg2 = rg2[orden]
    ree2 = ree2[orden]

    if not np.array_equal(N, valores_N):
        raise ValueError(
            f"Los tamaños de k={ke} no coinciden "
            "con los 13 esperados."
        )

    return rg2, ree2


# ============================================================
# 3. CARGAR LAS DOS RIGIDECES
# ============================================================

Rg2_100, Ree2_100 = cargar_resultados(100)

Rg2_1000, Ree2_1000 = cargar_resultados(1000)

# ============================================================
# 4. PREDICCIONES TEORICAS
# ============================================================

def calcular_b2(ke):

    sigma2 = kBT / ke

    return (
        l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
    ) / (l0**2 + sigma2)


b2_100 = calcular_b2(100)
b2_1000 = calcular_b2(1000)

N_suave = np.linspace(4, 64, 500)

# Radio de giro teorico
Rg2_teor_100 = (
    (N_suave**2 - 1)/(6*N_suave)*b2_100
)

Rg2_teor_1000 = (
    (N_suave**2 - 1)/(6*N_suave)*b2_1000
)

# Distancia entre extremos teorica
Ree2_teor_100 = (N_suave - 1)*b2_100
Ree2_teor_1000 = (N_suave - 1)*b2_1000

# ============================================================
# 5. FUNCION PARA CREAR LA GRAFICA COMPARATIVA
# ============================================================

def grafica_comparativa(
    datos_100,
    datos_1000,
    teor_100,
    teor_1000,
    ylabel,
    titulo,
    nombre_salida
):

    # Ajustes lineales independientes
    m100, c100 = np.polyfit(
        valores_N, datos_100, 1
    )

    m1000, c1000 = np.polyfit(
        valores_N, datos_1000, 1
    )

    fig, ax = plt.subplots(figsize=(10, 6.5))

    # --------------------------------------------------------
    # PUNTOS DE SIMULACION
    # --------------------------------------------------------

    ax.plot(
        valores_N,
        datos_100,
        "o",
        color="tab:blue",
        markersize=6,
        label=r"Simulación $k_e=100$"
    )

    ax.plot(
        valores_N,
        datos_1000,
        "o",
        color="tab:orange",
        markersize=6,
        label=r"Simulación $k_e=1000$"
    )

    # --------------------------------------------------------
    # AJUSTES LINEALES
    # --------------------------------------------------------

    ax.plot(
        N_suave,
        m100*N_suave+c100,
        color="tab:blue",
        linewidth=1.8,
        label=(
            rf"Ajuste $k_e=100$: "
            rf"$m={m100:.4f}$"
        )
    )

    ax.plot(
        N_suave,
        m1000*N_suave+c1000,
        color="tab:orange",
        linewidth=1.8,
        label=(
            rf"Ajuste $k_e=1000$: "
            rf"$m={m1000:.4f}$"
        )
    )

    # --------------------------------------------------------
    # PREDICCIONES TEORICAS
    # --------------------------------------------------------

    ax.plot(
        N_suave,
        teor_100,
        color="tab:blue",
        linestyle="--",
        linewidth=1.4,
        alpha=0.8,
        label=r"Teoría $k_e=100$"
    )

    ax.plot(
        N_suave,
        teor_1000,
        color="tab:orange",
        linestyle="--",
        linewidth=1.4,
        alpha=0.8,
        label=r"Teoría $k_e=1000$"
    )

    # --------------------------------------------------------
    # FORMATO
    # --------------------------------------------------------

    ax.set_title(titulo, fontsize=14)

    ax.set_xlabel("Número de monómeros, N")
    ax.set_ylabel(ylabel)

    ax.set_xlim(2, 67)

    ax.grid(alpha=0.3)

    ax.legend(
        loc="upper left",
        fontsize=9,
        ncol=2
    )

    fig.tight_layout()

    salida = carpeta_salida / nombre_salida

    if salida.exists():
        print(
            f"Ya existe {salida.name}. "
            "No se sobrescribe."
        )
    else:
        fig.savefig(
            salida,
            dpi=300,
            bbox_inches="tight"
        )
        print(f"Guardada: {salida}")

    plt.show()
    plt.close(fig)

    print(f"Pendiente ke=100:  {m100:.6f}")
    print(f"Pendiente ke=1000: {m1000:.6f}")


# ============================================================
# 6. COMPARACION DEL RADIO DE GIRO
# ============================================================

grafica_comparativa(
    Rg2_100,
    Rg2_1000,
    Rg2_teor_100,
    Rg2_teor_1000,
    r"$\langle R_g^2\rangle$",
    "Comparación del radio de giro para ambas rigideces",
    "Comparacion_Rg2_ke100_ke1000.png"
)

# ============================================================
# 7. COMPARACION DE DISTANCIA ENTRE EXTREMOS
# ============================================================

grafica_comparativa(
    Ree2_100,
    Ree2_1000,
    Ree2_teor_100,
    Ree2_teor_1000,
    r"$\langle R_{ee}^2\rangle$",
    "Comparación de la distancia entre extremos",
    "Comparacion_Ree2_ke100_ke1000.png"
)

print("\nComparacion final terminada.")

# %%

# %%
# RECUPERAR RESULTADOS DE ke=100 DESDE PYTHON

import numpy as np
from pathlib import Path

variables = ["valores_N", "Rg2_100", "Ree2_100"]

if all(nombre in globals() for nombre in variables):

    base = Path.home() / "Documents/A_5-Molecular_dynamics"

    carpeta = base / "Py/plots/OBJETIVO_3/RG/k100"
    carpeta.mkdir(parents=True, exist_ok=True)

    salida = carpeta / "Resumen_Rg_Ree_ke100_Recuperado.csv"

    tabla = np.column_stack([
        valores_N,
        Rg2_100,
        Ree2_100
    ])

    if salida.exists():
        print("El archivo ya existe. No se sobrescribe.")
    else:
        np.savetxt(
            salida,
            tabla,
            delimiter=",",
            header="N,Rg2_media,Ree2_media",
            comments="",
            fmt=["%d", "%.8f", "%.8f"]
        )

        print("RESULTADOS RECUPERADOS")
        print("Archivo guardado en:", salida)
        print(tabla)

else:
    print("Los resultados ya no estan en la memoria de Python.")
    print("Tendremos que recuperarlos de otra forma.")

# %%

# %%
# RECUPERAR TABLAS DE Rg² Y Ree² PARA ke=100

import numpy as np
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_datos = base / "Py/plots/OBJETIVO_3/Datos/k100/RG"
carpeta_salida = base / "Py/plots/OBJETIVO_3/RG/k100"

carpeta_salida.mkdir(parents=True, exist_ok=True)

valores_N = [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]

ke = 100
sigma2 = 1.0 / ke

# Longitud cuadratica media teorica del enlace
b2 = (1 + 6*sigma2 + 3*sigma2**2) / (1 + sigma2)

resultados = []

for N in valores_N:

    # Para N=64 utilizamos la simulacion larga
    if N == 64:
        archivo = carpeta_datos / "Rg_Ree_N64_ke100_T200000.txt"
    else:
        archivo = carpeta_datos / f"Rg_Ree_N{N}_ke100.txt"

    if not archivo.exists():
        raise FileNotFoundError(f"No encuentro: {archivo}")

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 5:
        raise ValueError(
            f"El archivo de N={N} no tiene las 5 columnas esperadas."
        )

    # Ultimos valores de los promedios acumulados
    Rg2_media = datos[-1, 1]
    Ree2_media = datos[-1, 2]

    # Predicciones teoricas
    Rg2_teorica = (N**2 - 1) * b2 / (6*N)
    Ree2_teorica = (N - 1) * b2

    resultados.append([
        N,
        Rg2_media,
        Ree2_media,
        Rg2_teorica,
        Ree2_teorica
    ])

tabla = np.array(resultados)

# ------------------------------------------------------------
# TABLA 1: DISTANCIA ENTRE EXTREMOS
# ------------------------------------------------------------

salida_ree = (
    carpeta_salida /
    "Ree2_vs_N_ke100_resultados_Claudia.csv"
)

if not salida_ree.exists():

    np.savetxt(
        salida_ree,
        tabla[:, [0, 2, 4]],
        delimiter=",",
        header="N,Ree2_media,Ree2_teorica",
        comments="",
        fmt=["%d", "%.8f", "%.8f"]
    )

    print("Tabla de Ree² guardada.")
else:
    print("La tabla de Ree² ya existe. No se sobrescribe.")

# ------------------------------------------------------------
# TABLA 2: RESUMEN DE AMBAS MAGNITUDES
# ------------------------------------------------------------

salida_resumen = (
    carpeta_salida /
    "Resumen_Rg_Ree_ke100_Claudia.csv"
)

if not salida_resumen.exists():

    np.savetxt(
        salida_resumen,
        tabla,
        delimiter=",",
        header=(
            "N,Rg2_media,Ree2_media,"
            "Rg2_teorica,Ree2_teorica"
        ),
        comments="",
        fmt=["%d", "%.8f", "%.8f", "%.8f", "%.8f"]
    )

    print("Tabla resumen guardada.")
else:
    print("La tabla resumen ya existe. No se sobrescribe.")

# ------------------------------------------------------------
# MOSTRAR LOS RESULTADOS
# ------------------------------------------------------------

print("\n   N        Rg²          Ree²")
print("--------------------------------")

for fila in tabla:
    print(
        f"{int(fila[0]):4d}   "
        f"{fila[1]:10.6f}   "
        f"{fila[2]:10.6f}"
    )

print("\nArchivos guardados en:")
print(carpeta_salida)






# %%
# ============================================================
# OBJETIVO 3 - RADIO DE GIRO SIN CUADRADO
# AJUSTES INDIVIDUALES PARA ke=100 Y ke=1000
#
# Magnitud: Rg_rms = sqrt(<Rg²>)
# Dos graficas por rigidez: normal y log-log
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"
carpeta_rg = base / "Py/plots/OBJETIVO_3/RG"

valores_N = np.array(
    [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
)

kBT = 1.0
l0 = 1.0

# ============================================================
# 1. FUNCION PARA PROCESAR UNA RIGIDEZ
# ============================================================

def analizar_radio_giro(ke):

    carpeta = carpeta_rg / f"k{ke}"

    archivo = (
        carpeta /
        f"Resumen_Rg_Ree_ke{ke}_Claudia.csv"
    )

    if not archivo.exists():
        raise FileNotFoundError(
            f"No encuentro el archivo: {archivo}"
        )

    # --------------------------------------------------------
    # LEER LOS RESULTADOS
    # --------------------------------------------------------

    datos = np.genfromtxt(
        archivo,
        delimiter=",",
        names=True,
        encoding="utf-8-sig",
        ndmin=1
    )

    N = np.asarray(datos["N"], dtype=float)
    Rg2 = np.asarray(datos["Rg2_media"], dtype=float)

    orden = np.argsort(N)
    N = N[orden]
    Rg2 = Rg2[orden]

    if not np.array_equal(N, valores_N):
        raise ValueError(
            f"Los valores de N no son los esperados "
            f"para ke={ke}"
        )

    if np.any(~np.isfinite(Rg2)) or np.any(Rg2 <= 0):
        raise ValueError("Datos de Rg² no validos.")

    # Radio de giro cuadratico medio
    Rg = np.sqrt(Rg2)

    # --------------------------------------------------------
    # PREDICCION TEORICA
    # --------------------------------------------------------

    sigma2 = kBT / ke

    b2 = (
        l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
    ) / (l0**2 + sigma2)

    N_suave = np.linspace(4, 64, 500)

    Rg_teor = np.sqrt(
        (N_suave**2 - 1)*b2/(6*N_suave)
    )

    # --------------------------------------------------------
    # AJUSTE DE POTENCIA
    # Rg = A * N^nu
    # --------------------------------------------------------

    coef, cov = np.polyfit(
        np.log(N),
        np.log(Rg),
        1,
        cov=True
    )

    nu, log_A = coef
    A = np.exp(log_A)

    error_nu = np.sqrt(cov[0, 0])

    Rg_ajuste = A * N_suave**nu

    log_pred = nu*np.log(N) + log_A

    R2 = 1 - (
        np.sum((np.log(Rg)-log_pred)**2)
        / np.sum((np.log(Rg)-np.mean(np.log(Rg)))**2)
    )

    print(f"\n========== ke = {ke} ==========")
    print(f"A = {A:.6f}")
    print(f"nu = {nu:.6f} +/- {error_nu:.6f}")
    print(f"R² log-log = {R2:.6f}")
    print(f"b² teorico = {b2:.6f}")

    # --------------------------------------------------------
    # GRAFICA 1: EJES LINEALES
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        N, Rg,
        "o",
        color="tab:blue",
        markersize=6,
        label="Simulación"
    )

    ax.plot(
        N_suave,
        Rg_ajuste,
        color="tab:red",
        linewidth=1.8,
        label=(
            rf"Ajuste: $R_{{g,\mathrm{{rms}}}}="
            rf"{A:.4f}N^{{{nu:.4f}}}$"
        )
    )

    ax.plot(
        N_suave,
        Rg_teor,
        "--",
        color="black",
        linewidth=1.5,
        label="Cadena ideal"
    )

    ax.set_title(
        rf"Radio de giro frente a $N$ ($k_e={ke}$)"
    )

    ax.set_xlabel("Número de monómeros, N")
    ax.set_ylabel(
        r"$R_{g,\mathrm{rms}}=\sqrt{\langle R_g^2\rangle}$"
    )

    ax.set_xlim(2, 67)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    fig.tight_layout()

    salida1 = (
        carpeta /
        f"Rg_rms_vs_N_ke{ke}_Claudia.png"
    )

    if not salida1.exists():
        fig.savefig(
            salida1,
            dpi=300,
            bbox_inches="tight"
        )
        print("Guardada:", salida1.name)
    else:
        print("Ya existe, no se sobrescribe:", salida1.name)

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # GRAFICA 2: ESCALA LOG-LOG
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        N, Rg,
        "o",
        color="tab:blue",
        markersize=6,
        label="Simulación"
    )

    ax.plot(
        N_suave,
        Rg_ajuste,
        color="tab:red",
        linewidth=1.8,
        label=(
            rf"Ajuste: $R_{{g,\mathrm{{rms}}}}="
            rf"{A:.4f}N^{{{nu:.4f}}}$"
            "\n"
            rf"$\nu={nu:.4f}\pm{error_nu:.4f}$"
            "\n"
            rf"$R^2={R2:.5f}$"
        )
    )

    ax.plot(
        N_suave,
        Rg_teor,
        "--",
        color="black",
        linewidth=1.5,
        label="Cadena ideal"
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_title(
        rf"Ajuste log-log del radio de giro ($k_e={ke}$)"
    )

    ax.set_xlabel("Número de monómeros, N")
    ax.set_ylabel(
        r"$R_{g,\mathrm{rms}}=\sqrt{\langle R_g^2\rangle}$"
    )

    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    fig.tight_layout()

    salida2 = (
        carpeta /
        f"Rg_rms_vs_N_ke{ke}_loglog_Claudia.png"
    )

    if not salida2.exists():
        fig.savefig(
            salida2,
            dpi=300,
            bbox_inches="tight"
        )
        print("Guardada:", salida2.name)
    else:
        print("Ya existe, no se sobrescribe:", salida2.name)

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # GUARDAR TABLA DEL RADIO DE GIRO
    # --------------------------------------------------------

    Rg_teor_puntos = np.sqrt(
        (N**2 - 1)*b2/(6*N)
    )

    tabla = np.column_stack([
        N, Rg, Rg_teor_puntos
    ])

    salida_csv = (
        carpeta /
        f"Rg_rms_vs_N_ke{ke}_resultados_Claudia.csv"
    )

    if not salida_csv.exists():

        np.savetxt(
            salida_csv,
            tabla,
            delimiter=",",
            header="N,Rg_rms,Rg_rms_teorica",
            comments="",
            fmt=["%d", "%.8f", "%.8f"]
        )

        print("Tabla guardada:", salida_csv.name)
    else:
        print("La tabla ya existe. No se sobrescribe.")

    return A, nu, error_nu, R2


# ============================================================
# 2. EJECUTAR LOS AJUSTES INDIVIDUALES
# ============================================================

resultado_100 = analizar_radio_giro(100)
resultado_1000 = analizar_radio_giro(1000)

print("\n========== RESUMEN FINAL ==========")

for ke, resultado in [
    (100, resultado_100),
    (1000, resultado_1000)
]:

    A, nu, error_nu, R2 = resultado

    print(
        f"ke={ke}: "
        f"nu={nu:.5f} +/- {error_nu:.5f}, "
        f"R²={R2:.5f}"
    )

print("\nCuatro graficas individuales terminadas.")


# %%




# %%
# ============================================================
# OBJETIVO 3 - DISTANCIA ENTRE EXTREMOS SIN CUADRADO
# AJUSTES INDIVIDUALES PARA ke=100 Y ke=1000
#
# Magnitud: Ree_rms = sqrt(<Ree²>)
# Dos graficas por rigidez: normal y log-log
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"
carpeta_rg = base / "Py/plots/OBJETIVO_3/RG"

valores_N = np.array(
    [4, 6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 64]
)

kBT = 1.0
l0 = 1.0

# ============================================================
# 1. FUNCION PARA PROCESAR CADA RIGIDEZ
# ============================================================

def analizar_distancia_extremos(ke):

    carpeta = carpeta_rg / f"k{ke}"

    archivo = (
        carpeta /
        f"Resumen_Rg_Ree_ke{ke}_Claudia.csv"
    )

    if not archivo.exists():
        raise FileNotFoundError(
            f"No encuentro el archivo: {archivo}"
        )

    # --------------------------------------------------------
    # LEER RESULTADOS
    # --------------------------------------------------------

    datos = np.genfromtxt(
        archivo,
        delimiter=",",
        names=True,
        encoding="utf-8-sig",
        ndmin=1
    )

    N = np.asarray(datos["N"], dtype=float)
    Ree2 = np.asarray(datos["Ree2_media"], dtype=float)

    orden = np.argsort(N)
    N = N[orden]
    Ree2 = Ree2[orden]

    if not np.array_equal(N, valores_N):
        raise ValueError(
            f"Los tamaños N no coinciden para ke={ke}"
        )

    if np.any(~np.isfinite(Ree2)) or np.any(Ree2 <= 0):
        raise ValueError("Datos de Ree² no validos.")

    # Distancia entre extremos cuadratica media
    Ree = np.sqrt(Ree2)

    # --------------------------------------------------------
    # PREDICCION TEORICA
    # --------------------------------------------------------

    sigma2 = kBT / ke

    b2 = (
        l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
    ) / (l0**2 + sigma2)

    N_suave = np.linspace(4, 64, 500)

    Ree_teor = np.sqrt(
        (N_suave - 1)*b2
    )

    # --------------------------------------------------------
    # AJUSTE DE POTENCIA: Ree = A*N^nu
    # --------------------------------------------------------

    coef, cov = np.polyfit(
        np.log(N),
        np.log(Ree),
        1,
        cov=True
    )

    nu, log_A = coef
    A = np.exp(log_A)

    error_nu = np.sqrt(cov[0, 0])

    Ree_ajuste = A*N_suave**nu

    log_pred = nu*np.log(N) + log_A

    R2 = 1 - (
        np.sum((np.log(Ree)-log_pred)**2)
        / np.sum((np.log(Ree)-np.mean(np.log(Ree)))**2)
    )

    print(f"\n========== ke = {ke} ==========")
    print(f"A = {A:.6f}")
    print(f"nu = {nu:.6f} +/- {error_nu:.6f}")
    print(f"R² log-log = {R2:.6f}")
    print(f"b² teorico = {b2:.6f}")

    # --------------------------------------------------------
    # GRAFICA 1: EJES LINEALES
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        N, Ree,
        "o",
        color="tab:blue",
        markersize=6,
        label="Simulación"
    )

    ax.plot(
        N_suave,
        Ree_ajuste,
        color="tab:red",
        linewidth=1.8,
        label=(
            rf"Ajuste: $R_{{ee,\mathrm{{rms}}}}="
            rf"{A:.4f}N^{{{nu:.4f}}}$"
        )
    )

    ax.plot(
        N_suave,
        Ree_teor,
        "--",
        color="black",
        linewidth=1.5,
        label="Cadena ideal"
    )

    ax.set_title(
        rf"Distancia entre extremos frente a $N$ ($k_e={ke}$)"
    )

    ax.set_xlabel("Número de monómeros, N")
    ax.set_ylabel(
        r"$R_{ee,\mathrm{rms}}=\sqrt{\langle R_{ee}^2\rangle}$"
    )

    ax.set_xlim(2, 67)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    fig.tight_layout()

    salida1 = (
        carpeta /
        f"Ree_rms_vs_N_ke{ke}_Claudia.png"
    )

    if not salida1.exists():
        fig.savefig(
            salida1,
            dpi=300,
            bbox_inches="tight"
        )
        print("Guardada:", salida1.name)
    else:
        print("Ya existe, no se sobrescribe:", salida1.name)

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # GRAFICA 2: ESCALA LOG-LOG
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        N, Ree,
        "o",
        color="tab:blue",
        markersize=6,
        label="Simulación"
    )

    ax.plot(
        N_suave,
        Ree_ajuste,
        color="tab:red",
        linewidth=1.8,
        label=(
            rf"Ajuste: $R_{{ee,\mathrm{{rms}}}}="
            rf"{A:.4f}N^{{{nu:.4f}}}$"
            "\n"
            rf"$\nu={nu:.4f}\pm{error_nu:.4f}$"
            "\n"
            rf"$R^2={R2:.5f}$"
        )
    )

    ax.plot(
        N_suave,
        Ree_teor,
        "--",
        color="black",
        linewidth=1.5,
        label="Cadena ideal"
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_title(
        rf"Ajuste log-log de la distancia entre extremos ($k_e={ke}$)"
    )

    ax.set_xlabel("Número de monómeros, N")
    ax.set_ylabel(
        r"$R_{ee,\mathrm{rms}}=\sqrt{\langle R_{ee}^2\rangle}$"
    )

    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    fig.tight_layout()

    salida2 = (
        carpeta /
        f"Ree_rms_vs_N_ke{ke}_loglog_Claudia.png"
    )

    if not salida2.exists():
        fig.savefig(
            salida2,
            dpi=300,
            bbox_inches="tight"
        )
        print("Guardada:", salida2.name)
    else:
        print("Ya existe, no se sobrescribe:", salida2.name)

    plt.show()
    plt.close(fig)

    # --------------------------------------------------------
    # GUARDAR TABLA DE DISTANCIAS ENTRE EXTREMOS
    # --------------------------------------------------------

    Ree_teor_puntos = np.sqrt(
        (N - 1)*b2
    )

    tabla = np.column_stack([
        N,
        Ree,
        Ree_teor_puntos
    ])

    salida_csv = (
        carpeta /
        f"Ree_rms_vs_N_ke{ke}_resultados_Claudia.csv"
    )

    if not salida_csv.exists():

        np.savetxt(
            salida_csv,
            tabla,
            delimiter=",",
            header="N,Ree_rms,Ree_rms_teorica",
            comments="",
            fmt=["%d", "%.8f", "%.8f"]
        )

        print("Tabla guardada:", salida_csv.name)
    else:
        print("La tabla ya existe. No se sobrescribe.")

    return A, nu, error_nu, R2


# ============================================================
# 2. EJECUTAR LOS AJUSTES INDIVIDUALES
# ============================================================

resultado_100 = analizar_distancia_extremos(100)
resultado_1000 = analizar_distancia_extremos(1000)

print("\n========== RESUMEN FINAL ==========")

for ke, resultado in [
    (100, resultado_100),
    (1000, resultado_1000)
]:

    A, nu, error_nu, R2 = resultado

    print(
        f"ke={ke}: "
        f"nu={nu:.5f} +/- {error_nu:.5f}, "
        f"R²={R2:.5f}"
    )

print("\nCuatro graficas de distancia entre extremos terminadas.")





# %%
# ============================================================
# OBJETIVO 3 - COMPARACION FINAL DE AMBAS RIGIDECES
#
# Radio de giro: sqrt(<Rg²>)
# Distancia entre extremos: sqrt(<Ree²>)
#
# ke = 100 y ke = 1000
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"
carpeta = base / "Py/plots/OBJETIVO_3/RG"

kBT = 1.0
l0 = 1.0

colores = {
    100: "tab:blue",
    1000: "tab:orange"
}

resultados = {}

# ============================================================
# 1. LEER LOS DATOS Y CALCULAR LAS MAGNITUDES
# ============================================================

for ke in [100, 1000]:

    archivo = (
        carpeta / f"k{ke}" /
        f"Resumen_Rg_Ree_ke{ke}_Claudia.csv"
    )

    if not archivo.exists():
        raise FileNotFoundError(
            f"No encuentro el CSV: {archivo}"
        )

    datos = np.genfromtxt(
        archivo,
        delimiter=",",
        names=True,
        encoding="utf-8-sig",
        ndmin=1
    )

    N = np.asarray(datos["N"], dtype=float)
    Rg2 = np.asarray(datos["Rg2_media"], dtype=float)
    Ree2 = np.asarray(datos["Ree2_media"], dtype=float)

    orden = np.argsort(N)

    N = N[orden]
    Rg2 = Rg2[orden]
    Ree2 = Ree2[orden]

    if np.any(Rg2 <= 0) or np.any(Ree2 <= 0):
        raise ValueError(f"Valores invalidos para ke={ke}")

    Rg = np.sqrt(Rg2)
    Ree = np.sqrt(Ree2)

    # Longitud cuadratica media teorica del enlace
    sigma2 = kBT / ke

    b2 = (
        l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
    ) / (l0**2 + sigma2)

    resultados[ke] = {
        "N": N,
        "Rg": Rg,
        "Ree": Ree,
        "b2": b2
    }

# Comprobar que ambas series usan los mismos tamaños
if not np.array_equal(
    resultados[100]["N"],
    resultados[1000]["N"]
):
    raise ValueError("Los valores de N no coinciden.")

N_suave = np.linspace(4, 64, 500)


# ============================================================
# 2. FUNCION PARA GENERAR LAS GRAFICAS
# ============================================================

def comparar_magnitud(magnitud, titulo, ylabel, nombre):

    fig, ax = plt.subplots(figsize=(10, 6.5))

    for ke in [100, 1000]:

        r = resultados[ke]

        N = r["N"]
        y = r[magnitud]
        b2 = r["b2"]
        color = colores[ke]

        # ----------------------------------------------------
        # AJUSTE DE POTENCIA: y = A*N^nu
        # ----------------------------------------------------

        coef, cov = np.polyfit(
            np.log(N),
            np.log(y),
            1,
            cov=True
        )

        nu, log_A = coef
        A = np.exp(log_A)
        error_nu = np.sqrt(cov[0, 0])

        y_ajuste = A * N_suave**nu

        # ----------------------------------------------------
        # PREDICCION TEORICA
        # ----------------------------------------------------

        if magnitud == "Rg":

            y_teorica = np.sqrt(
                (N_suave**2 - 1)*b2/(6*N_suave)
            )

        else:

            y_teorica = np.sqrt(
                (N_suave - 1)*b2
            )

        # ----------------------------------------------------
        # PUNTOS DE SIMULACION
        # ----------------------------------------------------

        ax.plot(
            N,
            y,
            "o",
            color=color,
            markersize=6,
            label=rf"Simulación $k_e={ke}$",
            zorder=3
        )

        # ----------------------------------------------------
        # AJUSTE EXPERIMENTAL
        # ----------------------------------------------------

        ax.plot(
            N_suave,
            y_ajuste,
            color=color,
            linewidth=1.8,
            label=(
                rf"Ajuste $k_e={ke}$: "
                rf"$\nu={nu:.4f}\pm{error_nu:.4f}$"
            )
        )

        # ----------------------------------------------------
        # TEORIA
        # ----------------------------------------------------

        ax.plot(
            N_suave,
            y_teorica,
            "--",
            color=color,
            linewidth=1.5,
            alpha=0.8,
            label=rf"Teoría $k_e={ke}$"
        )

        print(
            f"{magnitud}, ke={ke}: "
            f"A={A:.6f}, "
            f"nu={nu:.6f} +/- {error_nu:.6f}"
        )

    # ========================================================
    # FORMATO DE LA GRAFICA
    # ========================================================

    ax.set_title(titulo, fontsize=14)
    ax.set_xlabel("Número de monómeros, N")
    ax.set_ylabel(ylabel)

    ax.set_xlim(2, 67)
    ax.grid(alpha=0.3)

    ax.legend(
        loc="upper left",
        fontsize=9,
        ncol=2
    )

    fig.tight_layout()

    salida = carpeta / nombre

    if salida.exists():
        print(
            f"El archivo {nombre} ya existe. "
            "No se sobrescribe."
        )
    else:
        fig.savefig(
            salida,
            dpi=300,
            bbox_inches="tight"
        )
        print("Guardada:", salida)

    plt.show()
    plt.close(fig)


# ============================================================
# 3. COMPARACION DEL RADIO DE GIRO
# ============================================================

comparar_magnitud(
    magnitud="Rg",
    titulo="Comparación del radio de giro para ambas rigideces",
    ylabel=r"$R_{g,\mathrm{rms}}=\sqrt{\langle R_g^2\rangle}$",
    nombre="Comparacion_Rg_rms_ke100_ke1000_Claudia.png"
)


# ============================================================
# 4. COMPARACION DE LA DISTANCIA ENTRE EXTREMOS
# ============================================================

comparar_magnitud(
    magnitud="Ree",
    titulo="Comparación de la distancia entre extremos",
    ylabel=r"$R_{ee,\mathrm{rms}}=\sqrt{\langle R_{ee}^2\rangle}$",
    nombre="Comparacion_Ree_rms_ke100_ke1000_Claudia.png"
)

print("\nLas dos comparaciones han terminado.")


# %%

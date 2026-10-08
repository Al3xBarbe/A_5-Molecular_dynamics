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
if True:

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


# --- Gráficas de termalización: comparación de algoritmos ---
import numpy as np
import matplotlib.pyplot as plt

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
        nombre = "cinetica"

    elif tipo_energia == "potencial":
        E_euler = V_eta[mask_euler]
        E_rk = V_eta[mask_rk]
        E_ve = V_eta[mask_ve]

        ylabel = "Energía potencial media"
        nombre = "potencial"

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


# Generar todas las gráficas
plot_termalizacion(0.1, "cinetica")
plot_termalizacion(0.1, "potencial")

plot_termalizacion(1.0, "cinetica")
plot_termalizacion(1.0, "potencial")

plot_termalizacion(10.0, "cinetica")
plot_termalizacion(10.0, "potencial")
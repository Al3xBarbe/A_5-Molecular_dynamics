
# %%
# ============================================================
# OBJETIVO 3 - EQUIPARTICION DE LA ENERGIA
# k_e = 100 | SIMULACIONES ULTRALARGAS
#
# 10 graficas de energia cinetica
# 10 graficas de energia potencial
#
# Curva ROJA: media acumulada de la simulacion
# Linea VERDE: teoria habitual
# Linea MORADA: teoria corregida con jacobiano (solo V)
#
# Incluye la media final simulada como TEXTO en la leyenda.
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

# ============================================================
# 1. RUTAS Y PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k100/EQUIP_E"
)

carpeta_salida = (
    base / "Py/plots/OBJETIVO_3/EQUIP_E/k100"
)

carpeta_salida.mkdir(parents=True, exist_ok=True)

valores_N = [4, 8, 12, 16, 20, 24, 32, 40, 48, 64]

ke = 100.0
kBT = 1.0
l0 = 1.0

# Porcentaje del comienzo de la medicion que no dibujamos.
# No se cambia el calculo de la media acumulada.
RECORTE_INICIAL = 0.001

# Margen alrededor de las curvas y teorias en el eje Y.
MARGEN_Y = 0.10

# ============================================================
# 2. EXPRESIONES TEORICAS
# ============================================================

def energias_teoricas(N):

    K_teo = (3.0 / 2.0) * N * kBT

    V_arm = (N - 1) * kBT / 2.0

    sigma2 = kBT / ke

    V_jac = V_arm * (
        (l0**2 + 3.0 * sigma2) /
        (l0**2 + sigma2)
    )

    return K_teo, V_arm, V_jac


# ============================================================
# 3. AJUSTAR EL EJE Y
# ============================================================

def limites_y(curva, referencias):

    valores = np.concatenate([
        np.asarray(curva, dtype=float),
        np.asarray(referencias, dtype=float)
    ])

    minimo = np.min(valores)
    maximo = np.max(valores)

    amplitud = maximo - minimo

    if amplitud < 1e-12:
        amplitud = max(abs(maximo) * 0.01, 0.01)

    return (
        minimo - MARGEN_Y * amplitud,
        maximo + MARGEN_Y * amplitud
    )


# ============================================================
# 4. DIBUJAR CADA ENERGIA EN UNA IMAGEN INDEPENDIENTE
# ============================================================

def plot_energia(N, tiempo, energia, tipo):

    K_teo, V_arm, V_jac = energias_teoricas(N)

    # Recorte SOLO VISUAL del comienzo de la medicion
    indice_inicio = int(len(tiempo) * RECORTE_INICIAL)

    t = tiempo[indice_inicio:]
    E = energia[indice_inicio:]

    if len(t) < 2:
        raise ValueError(f"N={N}: datos insuficientes")

    media_final = energia[-1]

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # --------------------------------------------------------
    # LINEA ROJA: MEDIA ACUMULADA
    # --------------------------------------------------------

    ax.plot(
        t,
        E,
        color="tab:red",
        linewidth=1.8,
        label="Media acumulada",
        zorder=3
    )

    # --------------------------------------------------------
    # LINEAS TEORICAS
    # --------------------------------------------------------

    if tipo == "K":

        ax.axhline(
            K_teo,
            color="tab:green",
            linestyle="--",
            linewidth=1.8,
            label=(
                rf"$\langle K\rangle_{{\mathrm{{teo}}}}"
                rf"=\frac{{3N}}{{2}}k_BT={K_teo:.4f}$"
            )
        )

        referencias = [K_teo]

        ax.set_title(
            rf"Energía cinética ($N={N}$, $k_e=100$)",
            fontsize=14
        )

        ax.set_ylabel(
            r"$\langle K\rangle_t$",
            fontsize=12
        )

        nombre = (
            f"Energia_cinetica_N{N}_ke100_"
            "SUPERLARGA_Claudia.png"
        )

        texto_media = (
            rf"Media final simulada: "
            rf"$\langle K\rangle={media_final:.4f}$"
        )

    elif tipo == "V":

        ax.axhline(
            V_arm,
            color="tab:green",
            linestyle="--",
            linewidth=1.8,
            label=(
                rf"$\langle V\rangle_{{\mathrm{{arm}}}}"
                rf"=\frac{{N-1}}{{2}}k_BT={V_arm:.4f}$"
            )
        )

        ax.axhline(
            V_jac,
            color="tab:purple",
            linestyle="-.",
            linewidth=1.9,
            label=(
                rf"$\langle V\rangle_{{\mathrm{{Jac}}}}"
                rf"={V_jac:.4f}$"
            )
        )

        referencias = [V_arm, V_jac]

        ax.set_title(
            rf"Energía potencial ($N={N}$, $k_e=100$)",
            fontsize=14
        )

        ax.set_ylabel(
            r"$\langle V\rangle_t$",
            fontsize=12
        )

        nombre = (
            f"Energia_potencial_N{N}_ke100_"
            "SUPERLARGA_Claudia.png"
        )

        texto_media = (
            rf"Media final simulada: "
            rf"$\langle V\rangle={media_final:.4f}$"
        )

    else:
        raise ValueError("Tipo debe ser 'K' o 'V'")

    # --------------------------------------------------------
    # INCLUIR MEDIA FINAL COMO TEXTO EN LEYENDA
    # SIN DIBUJAR NINGUNA LINEA ADICIONAL
    # --------------------------------------------------------

    handle_texto = Line2D(
        [], [],
        linestyle="None",
        color="none",
        label=texto_media
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.append(handle_texto)
    labels.append(texto_media)

    # --------------------------------------------------------
    # FORMATO DE LOS EJES
    # --------------------------------------------------------

    ymin, ymax = limites_y(E, referencias)

    ax.set_ylim(ymin, ymax)
    ax.set_xlim(t[0], t[-1])

    ax.set_xlabel("Tiempo", fontsize=12)

    # Mostrar factor 10^6 automaticamente en el eje X
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((6, 6))
    formatter.set_useOffset(False)

    ax.xaxis.set_major_formatter(formatter)

    ax.grid(
        linestyle=":",
        alpha=0.25
    )

    ax.legend(
        handles,
        labels,
        fontsize=9,
        loc="best",
        framealpha=0.95
    )

    fig.tight_layout()

    # --------------------------------------------------------
    # GUARDAR
    # --------------------------------------------------------

    archivo_salida = carpeta_salida / nombre

    fig.savefig(
        archivo_salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print("Guardada:", archivo_salida.name)


# ============================================================
# 5. LEER LOS DATOS DE TODAS LAS CADENAS
# ============================================================

resumen = []

for N in valores_N:

    archivo = (
        carpeta_datos /
        f"Energia_Rg_Ree_N{N}_ke100_SUPERLARGA.txt"
    )

    archivo_flag = (
        carpeta_datos /
        f"COMPLETADO_N{N}_ke100_SUPERLARGA.flag"
    )

    if not archivo.exists():
        print(f"N={N}: no se encuentra el TXT, se omite.")
        continue

    if not archivo_flag.exists():
        print(
            f"N={N}: no tiene marca de completado; "
            "se omite."
        )
        continue

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 13:
        raise ValueError(
            f"N={N}: se esperaban 13 columnas "
            f"y se han encontrado {datos.shape[1]}"
        )

    if not np.all(np.isfinite(datos)):
        raise ValueError(f"N={N}: hay datos no finitos")

    # --------------------------------------------------------
    # COLUMNAS DEL FORTRAN ULTRALARGA
    # --------------------------------------------------------
    # 0: Tiempo de medicion
    # 1: K instantanea
    # 2: V instantanea
    # 3: Rg
    # 4: Ree
    # 5: Rg2
    # 6: Ree2
    # 7: K acumulada
    # 8: V acumulada
    # 9: Rg acumulada
    # 10: Ree acumulada
    # 11: Rg2 acumulada
    # 12: Ree2 acumulada

    tiempo = datos[:, 0]
    K_media = datos[:, 7]
    V_media = datos[:, 8]

    K_teo, V_arm, V_jac = energias_teoricas(N)

    print("\n" + "=" * 55)
    print(f"N = {N}")
    print("=" * 55)

    print(
        f"<K> final = {K_media[-1]:.6f} "
        f"| teoría = {K_teo:.6f}"
    )

    print(
        f"<V> final = {V_media[-1]:.6f} "
        f"| armónica = {V_arm:.6f} "
        f"| jacobiano = {V_jac:.6f}"
    )

    # Una imagen de energia cinetica
    plot_energia(N, tiempo, K_media, "K")

    # Otra imagen de energia potencial
    plot_energia(N, tiempo, V_media, "V")

    resumen.append([
        N, K_media[-1], K_teo,
        V_media[-1], V_arm, V_jac
    ])


# ============================================================
# 6. GUARDAR RESUMEN NUMERICO
# ============================================================

if resumen:

    archivo_resumen = (
        carpeta_salida /
        "Resumen_equiparticion_ke100_SUPERLARGA_Claudia.csv"
    )

    np.savetxt(
        archivo_resumen,
        np.array(resumen),
        delimiter=",",
        header=(
            "N,K_simulada,K_teorica,"
            "V_simulada,V_armonica,V_jacobiano"
        ),
        comments="",
        fmt=[
            "%d", "%.8f", "%.8f",
            "%.8f", "%.8f", "%.8f"
        ]
    )

    print("\nResumen numérico guardado.")

print("\n======================================")
print("GRÁFICAS DE EQUIPARTICIÓN TERMINADAS")
print("======================================")
print("Carpeta de salida:", carpeta_salida)

# %%

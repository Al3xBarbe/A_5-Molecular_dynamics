
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
print("======================================"
print("Carpeta de salida:", carpeta_salida)

# %%




# %%
# ============================================================
# OBJETIVO 3 - RADIO DE GIRO Y DISTANCIA ENTRE EXTREMOS
# k_e = 100 | SIMULACIONES SUPERLARGAS
#
# Graficas:
#   <Rg^2> frente al tiempo para cada N
#   <Ree^2> frente al tiempo para cada N
#
# ROJO: media acumulada de la simulacion
# VERDE: prediccion teorica
# Media final simulada: solo texto en leyenda
#
# Tambien se guardan las medias finales de:
# Rg, Ree, Rg^2, Ree^2
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

# ============================================================
# 1. PARAMETROS Y RUTAS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k100/EQUIP_E"
)

carpeta_rg = (
    base / "Py/plots/OBJETIVO_3/RG/k100"
)

carpeta_ree = (
    base / "Py/plots/OBJETIVO_3/REE/k100"
)

carpeta_rg.mkdir(parents=True, exist_ok=True)
carpeta_ree.mkdir(parents=True, exist_ok=True)

valores_N = [4, 8, 12, 16, 20, 24, 32, 40, 48, 64]

ke = 100.0
kBT = 1.0
l0 = 1.0

# Recorte VISUAL del comienzo de la medicion.
# Las medias acumuladas siguen usando todos los datos.
RECORTE_INICIAL = 0.001

# Margen vertical respecto a curva y teoria.
MARGEN_Y = 0.10

# ============================================================
# 2. VALORES TEORICOS
# ============================================================

sigma2 = kBT / ke

# Longitud cuadratica media de un enlace
b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (
    l0**2 + sigma2
)

def teorias_radios(N):

    Rg2_teo = (N**2 - 1) / (6*N) * b2

    Ree2_teo = (N - 1) * b2

    return Rg2_teo, Ree2_teo


print(f"Longitud cuadratica media del enlace: b^2 = {b2:.8f}")

# ============================================================
# 3. FUNCION PARA LIMITES DEL EJE Y
# ============================================================

def limites_y_radio(curva, referencia):

    valores = np.append(curva, referencia)

    ymin = np.min(valores)
    ymax = np.max(valores)

    amplitud = ymax - ymin

    if amplitud < 1e-12:
        amplitud = max(abs(ymax)*0.01, 0.01)

    return (
        ymin - MARGEN_Y*amplitud,
        ymax + MARGEN_Y*amplitud
    )


# ============================================================
# 4. FUNCION PARA GRAFICAR
# ============================================================

def graficar_radio(N, tiempo, media, tipo):

    Rg2_teo, Ree2_teo = teorias_radios(N)

    # Recortar solo la parte visible
    i0 = int(RECORTE_INICIAL * len(tiempo))

    t = tiempo[i0:]
    y = media[i0:]

    if len(t) < 2:
        raise ValueError(f"N={N}: datos insuficientes.")

    media_final = media[-1]

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # --------------------------------------------------------
    # MEDIA ACUMULADA: ROJO
    # --------------------------------------------------------

    ax.plot(
        t,
        y,
        color="tab:red",
        linewidth=1.8,
        label="Media acumulada",
        zorder=3
    )

    # --------------------------------------------------------
    # REFERENCIA TEORICA
    # --------------------------------------------------------

    if tipo == "RG":

        teoria = Rg2_teo

        ax.axhline(
            teoria,
            color="tab:green",
            linestyle="--",
            linewidth=1.8,
            label=(
                rf"$\langle R_g^2\rangle_{{\mathrm{{teo}}}}"
                rf"=\frac{{N^2-1}}{{6N}}b^2"
                rf"={teoria:.4f}$"
            )
        )

        texto_final = (
            rf"Media final simulada: "
            rf"$\langle R_g^2\rangle={media_final:.4f}$"
        )

        ax.set_title(
            rf"Radio de giro ($N={N}$, $k_e=100$)",
            fontsize=14
        )

        ax.set_ylabel(r"$\langle R_g^2\rangle_t$", fontsize=12)

        nombre = f"Rg2_tiempo_N{N}_ke100_SUPERLARGA_Claudia.png"
        carpeta = carpeta_rg

    elif tipo == "REE":

        teoria = Ree2_teo

        ax.axhline(
            teoria,
            color="tab:green",
            linestyle="--",
            linewidth=1.8,
            label=(
                rf"$\langle R_{{ee}}^2\rangle_{{\mathrm{{teo}}}}"
                rf"=(N-1)b^2={teoria:.4f}$"
            )
        )

        texto_final = (
            rf"Media final simulada: "
            rf"$\langle R_{{ee}}^2\rangle={media_final:.4f}$"
        )

        ax.set_title(
            rf"Distancia entre extremos ($N={N}$, $k_e=100$)",
            fontsize=14
        )

        ax.set_ylabel(r"$\langle R_{ee}^2\rangle_t$", fontsize=12)

        nombre = f"Ree2_tiempo_N{N}_ke100_SUPERLARGA_Claudia.png"
        carpeta = carpeta_ree

    else:
        raise ValueError("Tipo de radio no valido.")

    # --------------------------------------------------------
    # MEDIA FINAL: SOLO TEXTO EN LA LEYENDA
    # --------------------------------------------------------

    handles, labels = ax.get_legend_handles_labels()

    handles.append(
        Line2D([], [], linestyle="None", color="none")
    )

    labels.append(texto_final)

    # --------------------------------------------------------
    # FORMATO
    # --------------------------------------------------------

    ymin, ymax = limites_y_radio(y, teoria)

    ax.set_ylim(ymin, ymax)
    ax.set_xlim(t[0], t[-1])

    ax.set_xlabel("Tiempo", fontsize=12)

    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((6, 6))
    formatter.set_useOffset(False)

    ax.xaxis.set_major_formatter(formatter)

    ax.grid(
        alpha=0.25,
        linestyle=":"
    )

    ax.legend(
        handles,
        labels,
        fontsize=9,
        loc="best",
        framealpha=0.95
    )

    fig.tight_layout()

    salida = carpeta / nombre

    fig.savefig(
        salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print("Guardada:", salida.name)


# ============================================================
# 5. LEER LOS TXT Y GENERAR LAS 20 GRAFICAS
# ============================================================

resumen_rg = []
resumen_ree = []

for N in valores_N:

    archivo = (
        carpeta_datos /
        f"Energia_Rg_Ree_N{N}_ke100_SUPERLARGA.txt"
    )

    flag = (
        carpeta_datos /
        f"COMPLETADO_N{N}_ke100_SUPERLARGA.flag"
    )

    if not archivo.exists() or not flag.exists():
        print(f"N={N}: faltan archivos completos, se omite.")
        continue

    datos = np.loadtxt(archivo, ndmin=2)

    if datos.shape[1] != 13:
        raise ValueError(
            f"N={N}: se esperaban 13 columnas."
        )

    if not np.all(np.isfinite(datos)):
        raise ValueError(
            f"N={N}: datos no finitos."
        )

    # --------------------------------------------------------
    # COLUMNAS RELEVANTES DEL FORTRAN
    # --------------------------------------------------------

    tiempo = datos[:, 0]

    # Radios instantaneos:
    # 3 Rg
    # 4 Ree
    # 5 Rg^2
    # 6 Ree^2

    # Medias acumuladas:
    # 9  <Rg>
    # 10 <Ree>
    # 11 <Rg^2>
    # 12 <Ree^2>

    Rg_media = datos[-1, 9]
    Ree_media = datos[-1, 10]

    Rg2_media = datos[:, 11]
    Ree2_media = datos[:, 12]

    Rg2_final = Rg2_media[-1]
    Ree2_final = Ree2_media[-1]

    Rg2_teo, Ree2_teo = teorias_radios(N)

    print("\n" + "="*55)
    print(f"N = {N}")
    print("="*55)

    print(f"<Rg>       = {Rg_media:.6f}")
    print(f"<Rg^2>     = {Rg2_final:.6f}")
    print(f"<Rg^2> teo = {Rg2_teo:.6f}")

    print(f"<Ree>       = {Ree_media:.6f}")
    print(f"<Ree^2>     = {Ree2_final:.6f}")
    print(f"<Ree^2> teo = {Ree2_teo:.6f}")

    # --------------------------------------------------------
    # GRAFICAS INDIVIDUALES
    # --------------------------------------------------------

    graficar_radio(N, tiempo, Rg2_media, "RG")
    graficar_radio(N, tiempo, Ree2_media, "REE")

    # --------------------------------------------------------
    # GUARDAR RESULTADOS PARA AJUSTES POSTERIORES
    # --------------------------------------------------------

    resumen_rg.append([
        N,
        Rg_media,
        Rg2_final,
        np.sqrt(Rg2_final),
        Rg2_teo,
        np.sqrt(Rg2_teo)
    ])

    resumen_ree.append([
        N,
        Ree_media,
        Ree2_final,
        np.sqrt(Ree2_final),
        Ree2_teo,
        np.sqrt(Ree2_teo)
    ])


# ============================================================
# 6. GUARDAR TABLAS RESUMEN PARA LOS AJUSTES
# ============================================================

if resumen_rg:

    salida_rg = (
        carpeta_rg /
        "Resumen_Rg_ke100_SUPERLARGA_Claudia.csv"
    )

    np.savetxt(
        salida_rg,
        np.array(resumen_rg),
        delimiter=",",
        header=(
            "N,Rg_media,Rg2_media,Rg_RMS,"
            "Rg2_teorico,Rg_RMS_teorico"
        ),
        comments="",
        fmt=[
            "%d", "%.10f", "%.10f",
            "%.10f", "%.10f", "%.10f"
        ]
    )

    print("\nResumen de radio de giro guardado:")
    print(salida_rg)


if resumen_ree:

    salida_ree = (
        carpeta_ree /
        "Resumen_Ree_ke100_SUPERLARGA_Claudia.csv"
    )

    np.savetxt(
        salida_ree,
        np.array(resumen_ree),
        delimiter=",",
        header=(
            "N,Ree_media,Ree2_media,Ree_RMS,"
            "Ree2_teorico,Ree_RMS_teorico"
        ),
        comments="",
        fmt=[
            "%d", "%.10f", "%.10f",
            "%.10f", "%.10f", "%.10f"
        ]
    )

    print("\nResumen de distancia entre extremos guardado:")
    print(salida_ree)


print("\n==========================================")
print("ANALISIS TEMPORAL DE RADIOS TERMINADO")
print("==========================================")
print("Radio de giro:", carpeta_rg)
print("Distancia entre extremos:", carpeta_ree)


# %%




# %%
# ============================================================
# OBJETIVO 3 - AJUSTES LINEALES DE RADIOS AL CUADRADO
# k_e = 100 | SIMULACIONES SUPERLARGAS
#
# <Rg^2> frente a N
# <Ree^2> frente a N
#
# AZUL: datos simulados
# ROJO: ajuste lineal
# VERDE: expresion teorica
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. CARPETAS Y PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_rg = base / "Py/plots/OBJETIVO_3/RG/k100"
carpeta_ree = base / "Py/plots/OBJETIVO_3/REE/k100"

ke = 100.0
kBT = 1.0
l0 = 1.0

sigma2 = kBT / ke

b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (
    l0**2 + sigma2
)

# ============================================================
# 2. LEER RESUMENES DE LA SIMULACION
# ============================================================

archivo_rg = (
    carpeta_rg /
    "Resumen_Rg_ke100_SUPERLARGA_Claudia.csv"
)

archivo_ree = (
    carpeta_ree /
    "Resumen_Ree_ke100_SUPERLARGA_Claudia.csv"
)

datos_rg = np.genfromtxt(
    archivo_rg, delimiter=",", names=True
)

datos_ree = np.genfromtxt(
    archivo_ree, delimiter=",", names=True
)

N_rg = datos_rg["N"]
Rg2 = datos_rg["Rg2_media"]

N_ree = datos_ree["N"]
Ree2 = datos_ree["Ree2_media"]

# ============================================================
# 3. FUNCION DE AJUSTE LINEAL
# ============================================================

def hacer_ajuste(N, valores, tipo):

    N = np.asarray(N, dtype=float)
    valores = np.asarray(valores, dtype=float)

    if len(N) < 3:
        raise ValueError("Hacen falta al menos 3 tamaños.")

    # Ajuste de la forma y = a*N + c
    parametros, cov = np.polyfit(
        N, valores, 1, cov=True
    )

    a, c = parametros

    errores = np.sqrt(np.diag(cov))
    error_a, error_c = errores

    # Coeficiente de determinacion
    pred = a*N + c

    ss_res = np.sum((valores-pred)**2)
    ss_tot = np.sum((valores-np.mean(valores))**2)

    R2 = 1 - ss_res/ss_tot

    # Eje suave para representar funciones
    N_suave = np.linspace(
        min(N), max(N), 400
    )

    ajuste = a*N_suave + c

    if tipo == "RG":

        teoria = (
            (N_suave**2-1)/(6*N_suave)*b2
        )

        pendiente_teorica = b2/6

        titulo = r"Radio de giro al cuadrado ($k_e=100$)"
        ylabel = r"$\langle R_g^2\rangle$"

        etiqueta_teoria = (
            r"Teoría: $\langle R_g^2\rangle"
            r"=\frac{N^2-1}{6N}b^2$"
        )

        carpeta = carpeta_rg
        nombre = "Ajuste_lineal_Rg2_ke100_SUPERLARGA_Claudia.png"

    elif tipo == "REE":

        teoria = (N_suave-1)*b2

        pendiente_teorica = b2

        titulo = r"Distancia entre extremos al cuadrado ($k_e=100$)"
        ylabel = r"$\langle R_{ee}^2\rangle$"

        etiqueta_teoria = (
            r"Teoría: $\langle R_{ee}^2\rangle"
            r"=(N-1)b^2$"
        )

        carpeta = carpeta_ree
        nombre = "Ajuste_lineal_Ree2_ke100_SUPERLARGA_Claudia.png"

    else:
        raise ValueError("Tipo no reconocido.")

    # ========================================================
    # 4. GRAFICA
    # ========================================================

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.scatter(
        N,
        valores,
        color="tab:blue",
        s=55,
        zorder=4,
        label="Simulación"
    )

    ax.plot(
        N_suave,
        ajuste,
        color="tab:red",
        linewidth=2,
        label=(
            rf"Ajuste: $y=({a:.5f})N"
            rf"{c:+.5f}$"
        )
    )

    ax.plot(
        N_suave,
        teoria,
        color="tab:green",
        linewidth=2,
        linestyle="--",
        label=etiqueta_teoria
    )

    ax.set_title(titulo, fontsize=14)

    ax.set_xlabel(
        "Número de monómeros $N$",
        fontsize=12
    )

    ax.set_ylabel(ylabel, fontsize=12)

    ax.grid(alpha=0.3, linestyle=":")

    ax.legend(fontsize=10, loc="best")

    fig.tight_layout()

    salida = carpeta / nombre

    fig.savefig(
        salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # 5. RESULTADOS NUMERICOS
    # ========================================================

    print("\n" + "="*60)
    print(f"AJUSTE LINEAL: {tipo}")
    print("="*60)

    print(f"Pendiente = {a:.8f} ± {error_a:.8f}")
    print(f"Ordenada  = {c:.8f} ± {error_c:.8f}")

    print(f"R² = {R2:.8f}")

    print(
        f"Pendiente teorica asintotica = "
        f"{pendiente_teorica:.8f}"
    )

    print("Grafica guardada:", salida)

    return a, error_a, c, error_c, R2


# ============================================================
# 6. REALIZAR LOS DOS AJUSTES
# ============================================================

resultado_rg = hacer_ajuste(
    N_rg, Rg2, "RG"
)

resultado_ree = hacer_ajuste(
    N_ree, Ree2, "REE"
)

print("\nAJUSTES LINEALES FINALIZADOS")

# %%





# %%
# RECONSTRUIR LOS CSV DE RG Y REE - ke=100

import numpy as np
import csv
from pathlib import Path

base = Path.home() / "Documents/A_5-Molecular_dynamics"

datos = base / "Py/plots/OBJETIVO_3/Datos/k100/EQUIP_E"

carpeta_rg = base / "Py/plots/OBJETIVO_3/RG/k100"
carpeta_ree = base / "Py/plots/OBJETIVO_3/REE/k100"

carpeta_rg.mkdir(parents=True, exist_ok=True)
carpeta_ree.mkdir(parents=True, exist_ok=True)

N_valores = [4, 8, 12, 16, 20, 24, 32, 40, 48, 64]

ke = 100.0
sigma2 = 1.0 / ke

b2 = (1 + 6*sigma2 + 3*sigma2**2) / (1 + sigma2)

resultados_rg = []
resultados_ree = []

for N in N_valores:

    archivo = datos / f"Energia_Rg_Ree_N{N}_ke100_SUPERLARGA.txt"

    if not archivo.exists():
        raise FileNotFoundError(
            f"No encuentro los datos originales de N={N}: {archivo}"
        )

    datos_N = np.loadtxt(archivo, comments="#", ndmin=2)

    # Columnas Python (numeradas desde cero):
    # 9 = Rg acumulado
    # 10 = Ree acumulado
    # 11 = Rg² acumulado
    # 12 = Ree² acumulado

    ultima = datos_N[-1]

    Rg_media = ultima[9]
    Ree_media = ultima[10]

    Rg2_media = ultima[11]
    Ree2_media = ultima[12]

    Rg_RMS = np.sqrt(Rg2_media)
    Ree_RMS = np.sqrt(Ree2_media)

    Rg2_teo = (N**2 - 1)/(6*N) * b2
    Ree2_teo = (N - 1) * b2

    resultados_rg.append([
        N,
        Rg_media,
        Rg2_media,
        Rg_RMS,
        Rg2_teo,
        np.sqrt(Rg2_teo)
    ])

    resultados_ree.append([
        N,
        Ree_media,
        Ree2_media,
        Ree_RMS,
        Ree2_teo,
        np.sqrt(Ree2_teo)
    ])

# Guardar los CSV en las carpetas k100 originales
# para que los encuentre el script de ajustes.

archivo_rg = (
    carpeta_rg / "Resumen_Rg_ke100_SUPERLARGA_Claudia.csv"
)

archivo_ree = (
    carpeta_ree / "Resumen_Ree_ke100_SUPERLARGA_Claudia.csv"
)

with open(archivo_rg, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "N", "Rg_media", "Rg2_media", "Rg_RMS",
        "Rg2_teorico", "Rg_RMS_teorico"
    ])

    writer.writerows(resultados_rg)

with open(archivo_ree, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "N", "Ree_media", "Ree2_media", "Ree_RMS",
        "Ree2_teorico", "Ree_RMS_teorico"
    ])

    writer.writerows(resultados_ree)

print("CSV reconstruidos correctamente:")
print(archivo_rg)
print(archivo_ree)











# %%
# ============================================================
# OBJETIVO 3 - AJUSTES DE RADIOS SIN CUADRAR
# ke = 100 | SIMULACIONES SUPERLARGAS
#
# COMPARACION CON LA MISMA REFERENCIA TEORICA RMS
#
# Para cada magnitud:
#   - Ajuste potencial R = A*N^nu
#   - Representacion log-log
#
# AZUL: datos de simulacion
# ROJO: ajuste potencial
# VERDE: prediccion teorica RMS
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import csv

# ============================================================
# 1. CARPETAS Y PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_rg = base / "Py/plots/OBJETIVO_3/RG/k100"
carpeta_ree = base / "Py/plots/OBJETIVO_3/REE/k100"

salida_rg = carpeta_rg / "normal"
salida_ree = carpeta_ree / "normal"

salida_rg.mkdir(parents=True, exist_ok=True)
salida_ree.mkdir(parents=True, exist_ok=True)

ke = 100.0
kBT = 1.0
l0 = 1.0

sigma2 = kBT / ke

b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

b = np.sqrt(b2)

# ============================================================
# 2. LEER LOS CSV
# ============================================================

archivo_rg = (
    carpeta_rg / "Resumen_Rg_ke100_SUPERLARGA_Claudia.csv"
)

archivo_ree = (
    carpeta_ree / "Resumen_Ree_ke100_SUPERLARGA_Claudia.csv"
)

# Buscar tambien en "cuadrado" si los moviste
if not archivo_rg.exists():
    archivo_rg = (
        carpeta_rg / "cuadrado" /
        "Resumen_Rg_ke100_SUPERLARGA_Claudia.csv"
    )

if not archivo_ree.exists():
    archivo_ree = (
        carpeta_ree / "cuadrado" /
        "Resumen_Ree_ke100_SUPERLARGA_Claudia.csv"
    )

if not archivo_rg.exists():
    raise FileNotFoundError(archivo_rg)

if not archivo_ree.exists():
    raise FileNotFoundError(archivo_ree)

datos_rg = np.genfromtxt(
    archivo_rg, delimiter=",", names=True,
    encoding="utf-8", ndmin=1
)

datos_ree = np.genfromtxt(
    archivo_ree, delimiter=",", names=True,
    encoding="utf-8", ndmin=1
)

N_rg = np.asarray(datos_rg["N"], dtype=float)
N_ree = np.asarray(datos_ree["N"], dtype=float)

Rg_RMS = np.asarray(datos_rg["Rg_RMS"], dtype=float)
Rg_medio = np.asarray(datos_rg["Rg_media"], dtype=float)

Ree_RMS = np.asarray(datos_ree["Ree_RMS"], dtype=float)
Ree_medio = np.asarray(datos_ree["Ree_media"], dtype=float)

# ============================================================
# 3. REFERENCIAS TEORICAS
# ============================================================

def teoria_rg(N):
    return b * np.sqrt((N**2 - 1)/(6*N))


def teoria_ree(N):
    return b * np.sqrt(N - 1)


# ============================================================
# 4. AJUSTE POTENCIAL POR LOGARITMOS
# ============================================================

def ajustar_potencia(N, R):

    N = np.asarray(N, dtype=float)
    R = np.asarray(R, dtype=float)

    x = np.log(N)
    y = np.log(R)

    parametros, cov = np.polyfit(
        x, y, 1, cov=True
    )

    nu, lnA = parametros

    A = np.exp(lnA)

    error_nu = np.sqrt(cov[0, 0])
    error_A = A * np.sqrt(cov[1, 1])

    pred_log = nu*x + lnA

    ss_res = np.sum((y - pred_log)**2)
    ss_tot = np.sum((y - np.mean(y))**2)

    R2 = 1 - ss_res/ss_tot

    return A, error_A, nu, error_nu, R2


# ============================================================
# 5. GRAFICAS NORMAL Y LOG-LOG
# ============================================================

def graficar_ajustes(
    N, R, nombre, simbolo, tipo,
    carpeta_salida, es_rms
):

    orden = np.argsort(N)

    N = np.asarray(N, dtype=float)[orden]
    R = np.asarray(R, dtype=float)[orden]

    A, error_A, nu, error_nu, R2 = ajustar_potencia(N, R)

    N_suave = np.linspace(
        N.min(), N.max(), 500
    )

    ajuste = A * N_suave**nu

    if tipo == "RG":

        teoria = teoria_rg(N_suave)

        etiqueta_teoria = (
            r"$\sqrt{\langle R_g^2\rangle}_{\mathrm{teo}}"
            r"=b\sqrt{\frac{N^2-1}{6N}}$"
        )

    elif tipo == "REE":

        teoria = teoria_ree(N_suave)

        etiqueta_teoria = (
            r"$\sqrt{\langle R_{ee}^2\rangle}_{\mathrm{teo}}"
            r"=b\sqrt{N-1}$"
        )

    else:
        raise ValueError("Tipo de radio incorrecto")

    # Indicar la diferencia entre promedio directo y RMS
    if not es_rms:
        etiqueta_teoria += " (referencia RMS)"

    etiqueta_ajuste = (
        rf"Ajuste: $R=({A:.4f})N^{{{nu:.4f}}}$"
    )

    # ========================================================
    # GRAFICA 1 - EJES LINEALES
    # ========================================================

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.scatter(
        N, R,
        color="tab:blue",
        s=60,
        zorder=4,
        label=simbolo
    )

    ax.plot(
        N_suave,
        ajuste,
        color="tab:red",
        linewidth=2,
        label=etiqueta_ajuste
    )

    ax.plot(
        N_suave,
        teoria,
        color="tab:green",
        linestyle="--",
        linewidth=1.9,
        label=etiqueta_teoria
    )

    ax.set_title(
        rf"{nombre} ($k_e=100$)",
        fontsize=14
    )

    ax.set_xlabel("Número de monómeros $N$", fontsize=12)
    ax.set_ylabel(simbolo, fontsize=12)

    ax.grid(alpha=0.3, linestyle=":")

    ax.legend(fontsize=9, loc="best")

    fig.tight_layout()

    nombre_normal = (
        f"{tipo}_{'RMS' if es_rms else 'MEDIA'}"
        "_ajuste_potencial_ke100_Claudia.png"
    )

    fig.savefig(
        carpeta_salida / nombre_normal,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # GRAFICA 2 - EJES LOG-LOG
    # ========================================================

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.scatter(
        N, R,
        color="tab:blue",
        s=60,
        zorder=4,
        label=simbolo
    )

    ax.plot(
        N_suave,
        ajuste,
        color="tab:red",
        linewidth=2,
        label=(
            rf"Ajuste: $\ln R=\ln({A:.4f})"
            rf"+({nu:.4f}\pm{error_nu:.4f})\ln N$"
        )
    )

    ax.plot(
        N_suave,
        teoria,
        color="tab:green",
        linestyle="--",
        linewidth=1.9,
        label=etiqueta_teoria
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_title(
        rf"{nombre} - log-log ($k_e=100$)",
        fontsize=14
    )

    ax.set_xlabel(r"$N$", fontsize=12)
    ax.set_ylabel(simbolo, fontsize=12)

    ax.grid(
        which="both",
        alpha=0.3,
        linestyle=":"
    )

    ax.legend(fontsize=9, loc="best")

    fig.tight_layout()

    nombre_loglog = (
        f"{tipo}_{'RMS' if es_rms else 'MEDIA'}"
        "_ajuste_loglog_ke100_Claudia.png"
    )

    fig.savefig(
        carpeta_salida / nombre_loglog,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print("\n" + "="*55)
    print(nombre)
    print("="*55)

    print(f"A  = {A:.6f} +/- {error_A:.6f}")
    print(f"nu = {nu:.6f} +/- {error_nu:.6f}")
    print(f"R2 log-log = {R2:.8f}")

    return [nombre, A, error_A, nu, error_nu, R2]


# ============================================================
# 6. EJECUTAR LOS CUATRO AJUSTES
# ============================================================

resultados = []

resultados.append(
    graficar_ajustes(
        N_rg,
        Rg_RMS,
        "Raíz del promedio cuadrático del radio de giro",
        r"$\sqrt{\langle R_g^2\rangle}$",
        "RG",
        salida_rg,
        es_rms=True
    )
)

resultados.append(
    graficar_ajustes(
        N_rg,
        Rg_medio,
        "Promedio directo del radio de giro",
        r"$\langle R_g\rangle$",
        "RG",
        salida_rg,
        es_rms=False
    )
)

resultados.append(
    graficar_ajustes(
        N_ree,
        Ree_RMS,
        "Raíz del promedio cuadrático entre extremos",
        r"$\sqrt{\langle R_{ee}^2\rangle}$",
        "REE",
        salida_ree,
        es_rms=True
    )
)

resultados.append(
    graficar_ajustes(
        N_ree,
        Ree_medio,
        "Promedio directo de la distancia entre extremos",
        r"$\langle R_{ee}\rangle$",
        "REE",
        salida_ree,
        es_rms=False
    )
)

# ============================================================
# 7. GUARDAR RESULTADOS DE LOS AJUSTES
# ============================================================

archivo_resultados = (
    base / "Py/plots/OBJETIVO_3/"
    "Ajustes_radios_normales_ke100_Claudia.csv"
)

with open(archivo_resultados, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "Magnitud", "A", "error_A",
        "nu", "error_nu", "R2_log"
    ])

    writer.writerows(resultados)

print("\nAJUSTES TERMINADOS")
print("Radio de giro:", salida_rg)
print("Distancia entre extremos:", salida_ree)




















# %%
# ============================================================
# OBJETIVO 3 - EQUIPARTICION DE ENERGIAS
# ke = 1000 | SIMULACIONES ULTRALARGAS
#
# Para cada N:
#   1. Energia cinetica media acumulada
#   2. Energia potencial media acumulada
#
# ROJO: media acumulada de la simulacion
# VERDE: prediccion teorica
# MORADO: correccion por jacobiano (solo V)
#
# Las leyendas incluyen expresiones y valores numericos.
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. PARAMETROS Y RUTAS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k1000/EQUIP_E"
)

carpeta_salida = (
    base / "Py/plots/OBJETIVO_3/EQUIP_E/k1000"
)

carpeta_salida.mkdir(parents=True, exist_ok=True)

N_valores = [4, 8, 12, 16, 20, 24, 32, 40, 48, 64]

ke = 1000.0
kBT = 1.0
l0 = 1.0

sigma2 = kBT / ke

# Mostrar el 85 % final de la evolucion temporal
RECORTE_INICIAL = 0.002

# ============================================================
# 2. PREDICCIONES TEORICAS
# ============================================================

def energia_cinetica_teorica(N):
    return 1.5 * N * kBT


def energia_potencial_armonica(N):
    return 0.5 * (N - 1) * kBT


def energia_potencial_jacobiano(N):
    return (
        0.5 * (N - 1) * kBT
        * (l0**2 + 3*sigma2)
        / (l0**2 + sigma2)
    )


# ============================================================
# 3. FUNCION PARA REPRESENTAR CADA ENERGIA
# ============================================================

def graficar_energia(
    N, tiempo, media_acumulada,
    tipo, teorias, archivo_salida
):

    tiempo = np.asarray(tiempo, dtype=float)
    media_acumulada = np.asarray(media_acumulada, dtype=float)

    # Seleccionar el tramo temporal representado
    indice_inicio = int(
        RECORTE_INICIAL * len(tiempo)
    )

    t = tiempo[indice_inicio:]
    y = media_acumulada[indice_inicio:]

    valor_final = media_acumulada[-1]

    fig, ax = plt.subplots(figsize=(10, 6))

    # --------------------------------------------------------
    # MEDIA ACUMULADA DE LA SIMULACION
    # --------------------------------------------------------

    ax.plot(
        t,
        y,
        color="tab:red",
        linewidth=1.8,
        label="Media acumulada (simulación)"
    )

    # --------------------------------------------------------
    # PREDICCIONES TEORICAS
    # --------------------------------------------------------

    for valor, etiqueta, color, estilo in teorias:

        ax.axhline(
            y=valor,
            color=color,
            linestyle=estilo,
            linewidth=1.8,
            label=etiqueta
        )

    # --------------------------------------------------------
    # VALOR FINAL SIMULADO (SOLO EN LA LEYENDA)
    # --------------------------------------------------------

    ax.plot(
        [],
        [],
        linestyle="None",
        label=(
            rf"Media final simulada: ${valor_final:.6f}$"
        )
    )

    # --------------------------------------------------------
    # TITULO Y EJES
    # --------------------------------------------------------

    ax.set_title(
        rf"Energía {tipo} media ($N={N}$, $k_e=1000$)",
        fontsize=14
    )

    ax.set_xlabel("Tiempo", fontsize=12)

    if tipo == "cinética":
        ax.set_ylabel(
            r"$\langle K\rangle$",
            fontsize=12
        )
    else:
        ax.set_ylabel(
            r"$\langle V\rangle$",
            fontsize=12
        )

    # Notacion cientifica en el eje temporal
    ax.ticklabel_format(
        axis="x",
        style="sci",
        scilimits=(6, 6),
        useMathText=True
    )

    # --------------------------------------------------------
    # LIMITES VERTICALES
    # --------------------------------------------------------

    # Incluir todos los valores simulados visibles
    # y las predicciones teoricas.

    valores = np.concatenate([
        y,
        np.array([teoria[0] for teoria in teorias])
    ])

    minimo = np.min(valores)
    maximo = np.max(valores)

    margen = max(
        0.15 * (maximo - minimo),
        0.002 * max(abs(maximo), 1.0)
    )

    ax.set_ylim(
        minimo - margen,
        maximo + margen
    )

    # --------------------------------------------------------
    # FORMATO Y GUARDADO
    # --------------------------------------------------------

    ax.grid(
        alpha=0.3,
        linestyle=":"
    )

    ax.legend(
        loc="best",
        fontsize=9
    )

    fig.tight_layout()

    fig.savefig(
        archivo_salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)


# ============================================================
# 4. LEER DATOS Y GENERAR LAS GRAFICAS
# ============================================================

resumen = []

for N in N_valores:

    archivo = (
        carpeta_datos /
        f"Energia_Rg_Ree_N{N}_ke1000_ULTRALARGA.txt"
    )

    if not archivo.exists():
        raise FileNotFoundError(
            f"No se encuentra el archivo para N={N}:\n{archivo}"
        )

    datos = np.loadtxt(
        archivo,
        comments="#",
        ndmin=2
    )

    if datos.shape[1] < 13:
        raise ValueError(
            f"El archivo N={N} tiene "
            f"{datos.shape[1]} columnas, "
            f"pero esperamos al menos 13."
        )

    # --------------------------------------------------------
    # COLUMNAS DEL ARCHIVO
    # --------------------------------------------------------
    #
    # 0  = tiempo de medicion
    # 1  = K instantanea
    # 2  = V instantanea
    # 3  = Rg
    # 4  = Ree
    # 5  = Rg²
    # 6  = Ree²
    # 7  = K media acumulada
    # 8  = V media acumulada
    # 9  = Rg media acumulada
    # 10 = Ree media acumulada
    # 11 = Rg² media acumulada
    # 12 = Ree² media acumulada
    # --------------------------------------------------------

    tiempo = datos[:, 0]

    K_acum = datos[:, 7]
    V_acum = datos[:, 8]

    # Predicciones teoricas
    K_teo = energia_cinetica_teorica(N)

    V_arm = energia_potencial_armonica(N)
    V_jac = energia_potencial_jacobiano(N)

    # ========================================================
    # GRAFICA DE ENERGIA CINETICA
    # ========================================================

    graficar_energia(
        N,
        tiempo,
        K_acum,
        tipo="cinética",
        teorias=[
            (
                K_teo,
                (
                    r"Teoría: "
                    r"$\langle K\rangle="
                    r"\frac{3N}{2}k_BT"
                    rf"={K_teo:.6f}$"
                ),
                "tab:green",
                "--"
            )
        ],
        archivo_salida=(
            carpeta_salida /
            f"Energia_cinetica_N{N}_ke1000_ULTRALARGA_Claudia.png"
        )
    )

    # ========================================================
    # GRAFICA DE ENERGIA POTENCIAL
    # ========================================================

    graficar_energia(
        N,
        tiempo,
        V_acum,
        tipo="potencial",
        teorias=[
            (
                V_arm,
                (
                    r"Armónica: "
                    r"$\langle V\rangle="
                    r"\frac{N-1}{2}k_BT"
                    rf"={V_arm:.6f}$"
                ),
                "tab:green",
                "--"
            ),
            (
                V_jac,
                rf"$\langle V\rangle_{{\mathrm{{Jac}}}} = {V_jac:.6f}$",
                "tab:purple",
                "-."
            )
        ],
        archivo_salida=(
            carpeta_salida /
            f"Energia_potencial_N{N}_ke1000_ULTRALARGA_Claudia.png"
        )
    )

    # ========================================================
    # GUARDAR RESULTADOS NUMERICOS
    # ========================================================

    resumen.append([
        N,
        K_acum[-1],
        K_teo,
        V_acum[-1],
        V_arm,
        V_jac
    ])

    print(
        f"N={N:2d} | "
        f"K_sim={K_acum[-1]:.6f} | "
        f"K_teo={K_teo:.6f} | "
        f"V_sim={V_acum[-1]:.6f} | "
        f"V_arm={V_arm:.6f} | "
        f"V_jac={V_jac:.6f}"
    )


# ============================================================
# 5. GUARDAR CSV RESUMEN
# ============================================================

resumen = np.asarray(resumen)

archivo_resumen = (
    carpeta_salida /
    "Resumen_equiparticion_ke1000_ULTRALARGA_Claudia.csv"
)

np.savetxt(
    archivo_resumen,
    resumen,
    delimiter=",",
    header=(
        "N,K_simulada,K_teorica,"
        "V_simulada,V_armonica,V_jacobiano"
    ),
    comments="",
    fmt="%.10f"
)

print("\n" + "="*55)
print("EQUIPARTICION DE ENERGIAS - ke=1000")
print("="*55)

print(f"Graficas generadas: {2 * len(N_valores)}")
print(f"Carpeta de salida: {carpeta_salida}")
print(f"Resumen numerico: {archivo_resumen}")












# %%
# ============================================================
# OBJETIVO 3 - RADIOS AL CUADRADO EN FUNCION DEL TIEMPO
# ke = 1000 | SIMULACIONES ULTRALARGAS
#
# RG:  <Rg²>  frente al tiempo
# REE: <Ree²> frente al tiempo
#
# ROJO: media acumulada de la simulacion
# VERDE: prediccion teorica
#
# Cada radio se guarda en su propia carpeta "cuadrado".
# Se generan tambien dos CSV resumen.
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import csv

# ============================================================
# 1. RUTAS Y PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_datos = (
    base / "Py/plots/OBJETIVO_3/Datos/k1000/EQUIP_E"
)

carpeta_rg = (
    base / "Py/plots/OBJETIVO_3/RG/k1000/cuadrado"
)

carpeta_ree = (
    base / "Py/plots/OBJETIVO_3/REE/k1000/cuadrado"
)

carpeta_rg.mkdir(parents=True, exist_ok=True)
carpeta_ree.mkdir(parents=True, exist_ok=True)

N_valores = [4, 8, 12, 16, 20, 24, 32, 40, 48, 64]

ke = 1000.0
kBT = 1.0
l0 = 1.0

sigma2 = kBT / ke

# Longitud cuadratica media del enlace
b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

# Mostrar el 85 % final de la evolucion temporal
RECORTE_INICIAL = 0.001

print(f"b² teorico = {b2:.8f}")

# ============================================================
# 2. PREDICCIONES TEORICAS
# ============================================================

def Rg2_teorico(N):
    return ((N**2 - 1) / (6*N)) * b2


def Ree2_teorico(N):
    return (N - 1) * b2


# ============================================================
# 3. FUNCION PARA REPRESENTAR LOS RADIOS
# ============================================================

def graficar_radio(
    N, tiempo, media_acumulada,
    valor_teorico, tipo, carpeta_salida
):

    tiempo = np.asarray(tiempo, dtype=float)
    media_acumulada = np.asarray(
        media_acumulada, dtype=float
    )

    indice_inicio = int(
        RECORTE_INICIAL * len(tiempo)
    )

    t = tiempo[indice_inicio:]
    y = media_acumulada[indice_inicio:]

    valor_final = media_acumulada[-1]

    fig, ax = plt.subplots(figsize=(10, 6))

    # --------------------------------------------------------
    # RADIO DE GIRO
    # --------------------------------------------------------

    if tipo == "RG":

        titulo = (
            rf"Radio de giro al cuadrado "
            rf"($N={N}$, $k_e=1000$)"
        )

        ylabel = r"$\langle R_g^2\rangle$"

        etiqueta_simulacion = (
            r"Simulación: $\langle R_g^2\rangle$"
        )

        etiqueta_teoria = (
            r"Teoría: "
            r"$\langle R_g^2\rangle="
            r"\frac{N^2-1}{6N}b^2"
            rf"={valor_teorico:.6f}$"
        )

        nombre_archivo = (
            f"Rg2_tiempo_N{N}_ke1000_ULTRALARGA_Claudia.png"
        )

    # --------------------------------------------------------
    # DISTANCIA ENTRE EXTREMOS
    # --------------------------------------------------------

    elif tipo == "REE":

        titulo = (
            rf"Distancia entre extremos al cuadrado "
            rf"($N={N}$, $k_e=1000$)"
        )

        ylabel = r"$\langle R_{ee}^2\rangle$"

        etiqueta_simulacion = (
            r"Simulación: $\langle R_{ee}^2\rangle$"
        )

        etiqueta_teoria = (
            r"Teoría: "
            r"$\langle R_{ee}^2\rangle="
            r"(N-1)b^2"
            rf"={valor_teorico:.6f}$"
        )

        nombre_archivo = (
            f"Ree2_tiempo_N{N}_ke1000_ULTRALARGA_Claudia.png"
        )

    else:
        raise ValueError("El tipo debe ser RG o REE")

    # --------------------------------------------------------
    # CURVA ROJA - SIMULACION
    # --------------------------------------------------------

    ax.plot(
        t,
        y,
        color="tab:red",
        linewidth=1.8,
        label=etiqueta_simulacion
    )

    # --------------------------------------------------------
    # LINEA VERDE - TEORIA
    # --------------------------------------------------------

    ax.axhline(
        valor_teorico,
        color="tab:green",
        linestyle="--",
        linewidth=1.8,
        label=etiqueta_teoria
    )

    # --------------------------------------------------------
    # MEDIA FINAL SOLO EN LA LEYENDA
    # --------------------------------------------------------

    ax.plot(
        [],
        [],
        linestyle="None",
        label=rf"Media final simulada: ${valor_final:.6f}$"
    )

    # --------------------------------------------------------
    # FORMATO
    # --------------------------------------------------------

    ax.set_title(titulo, fontsize=14)

    ax.set_xlabel("Tiempo", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)

    ax.ticklabel_format(
        axis="x",
        style="sci",
        scilimits=(6, 6),
        useMathText=True
    )

    minimo = min(np.min(y), valor_teorico)
    maximo = max(np.max(y), valor_teorico)

    margen = max(
        0.15 * (maximo - minimo),
        0.002 * max(abs(maximo), 1.0)
    )

    ax.set_ylim(
        minimo - margen,
        maximo + margen
    )

    ax.grid(alpha=0.3, linestyle=":")
    ax.legend(loc="best", fontsize=10)

    fig.tight_layout()

    # --------------------------------------------------------
    # GUARDAR GRAFICA
    # --------------------------------------------------------

    archivo_salida = carpeta_salida / nombre_archivo

    fig.savefig(
        archivo_salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)


# ============================================================
# 4. LEER LAS SIMULACIONES Y GENERAR LAS GRAFICAS
# ============================================================

resumen_rg = []
resumen_ree = []

for N in N_valores:

    archivo = (
        carpeta_datos /
        f"Energia_Rg_Ree_N{N}_ke1000_ULTRALARGA.txt"
    )

    if not archivo.exists():
        raise FileNotFoundError(
            f"No encuentro los datos originales de N={N}:\n{archivo}"
        )

    datos = np.loadtxt(
        archivo,
        comments="#",
        ndmin=2
    )

    if datos.shape[1] < 13:
        raise ValueError(
            f"El archivo de N={N} tiene "
            f"{datos.shape[1]} columnas. Se esperaban 13."
        )

    # --------------------------------------------------------
    # COLUMNAS DEL ARCHIVO
    # --------------------------------------------------------
    #
    # 0  = Tiempo
    # 3  = Rg instantaneo
    # 4  = Ree instantaneo
    # 5  = Rg² instantaneo
    # 6  = Ree² instantaneo
    # 9  = <Rg> acumulado
    # 10 = <Ree> acumulado
    # 11 = <Rg²> acumulado
    # 12 = <Ree²> acumulado
    # --------------------------------------------------------

    tiempo = datos[:, 0]

    Rg_acum = datos[:, 9]
    Ree_acum = datos[:, 10]

    Rg2_acum = datos[:, 11]
    Ree2_acum = datos[:, 12]

    # Medias finales simuladas
    Rg_media = Rg_acum[-1]
    Ree_media = Ree_acum[-1]

    Rg2_media = Rg2_acum[-1]
    Ree2_media = Ree2_acum[-1]

    # Raices cuadraticas medias
    Rg_RMS = np.sqrt(Rg2_media)
    Ree_RMS = np.sqrt(Ree2_media)

    # Predicciones teoricas
    Rg2_teo = Rg2_teorico(N)
    Ree2_teo = Ree2_teorico(N)

    # --------------------------------------------------------
    # GRAFICA DEL RADIO DE GIRO
    # --------------------------------------------------------

    graficar_radio(
        N,
        tiempo,
        Rg2_acum,
        Rg2_teo,
        "RG",
        carpeta_rg
    )

    # --------------------------------------------------------
    # GRAFICA DE DISTANCIA ENTRE EXTREMOS
    # --------------------------------------------------------

    graficar_radio(
        N,
        tiempo,
        Ree2_acum,
        Ree2_teo,
        "REE",
        carpeta_ree
    )

    # --------------------------------------------------------
    # RESUMENES PARA LOS AJUSTES EN FUNCION DE N
    # --------------------------------------------------------

    resumen_rg.append([
        N,
        Rg_media,
        Rg2_media,
        Rg_RMS,
        Rg2_teo,
        np.sqrt(Rg2_teo)
    ])

    resumen_ree.append([
        N,
        Ree_media,
        Ree2_media,
        Ree_RMS,
        Ree2_teo,
        np.sqrt(Ree2_teo)
    ])

    print(
        f"N={N:2d} | "
        f"Rg² sim={Rg2_media:.6f}, teo={Rg2_teo:.6f} | "
        f"Ree² sim={Ree2_media:.6f}, teo={Ree2_teo:.6f}"
    )


# ============================================================
# 5. GUARDAR CSV DEL RADIO DE GIRO
# ============================================================

archivo_rg = (
    carpeta_rg /
    "Resumen_Rg_ke1000_ULTRALARGA_Claudia.csv"
)

with open(archivo_rg, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "N",
        "Rg_media",
        "Rg2_media",
        "Rg_RMS",
        "Rg2_teorico",
        "Rg_RMS_teorico"
    ])

    writer.writerows(resumen_rg)


# ============================================================
# 6. GUARDAR CSV DE DISTANCIA ENTRE EXTREMOS
# ============================================================

archivo_ree = (
    carpeta_ree /
    "Resumen_Ree_ke1000_ULTRALARGA_Claudia.csv"
)

with open(archivo_ree, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "N",
        "Ree_media",
        "Ree2_media",
        "Ree_RMS",
        "Ree2_teorico",
        "Ree_RMS_teorico"
    ])

    writer.writerows(resumen_ree)


print("\n" + "="*60)
print("RADIOS AL CUADRADO - ke=1000 FINALIZADOS")
print("="*60)

print(f"Graficas generadas: {2 * len(N_valores)}")

print("\nRadio de giro:")
print(carpeta_rg)
print("Resumen:", archivo_rg)

print("\nDistancia entre extremos:")
print(carpeta_ree)
print("Resumen:", archivo_ree)






# %%
# ============================================================
# OBJETIVO 3 - AJUSTES LINEALES DE LOS RADIOS AL CUADRADO
# ke = 1000 | SIMULACIONES ULTRALARGAS
#
# Ajustes:
#   <Rg²>  = a_g*N  + c_g
#   <Ree²> = a_ee*N + c_ee
#
# AZUL: datos de simulacion
# ROJO: ajuste lineal
# VERDE: prediccion teorica
#
# Cada grafica se guarda en su carpeta "cuadrado".
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. RUTAS Y PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_rg = (
    base / "Py/plots/OBJETIVO_3/RG/k1000/cuadrado"
)

carpeta_ree = (
    base / "Py/plots/OBJETIVO_3/REE/k1000/cuadrado"
)

carpeta_rg.mkdir(parents=True, exist_ok=True)
carpeta_ree.mkdir(parents=True, exist_ok=True)

archivo_rg = (
    carpeta_rg /
    "Resumen_Rg_ke1000_ULTRALARGA_Claudia.csv"
)

archivo_ree = (
    carpeta_ree /
    "Resumen_Ree_ke1000_ULTRALARGA_Claudia.csv"
)

for archivo in [archivo_rg, archivo_ree]:
    if not archivo.exists():
        raise FileNotFoundError(
            f"No se encuentra el CSV:\n{archivo}"
        )

ke = 1000.0
kBT = 1.0
l0 = 1.0

sigma2 = kBT / ke

b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

# ============================================================
# 2. LEER LOS RESULTADOS
# ============================================================

datos_rg = np.genfromtxt(
    archivo_rg,
    delimiter=",",
    names=True,
    encoding="utf-8",
    ndmin=1
)

datos_ree = np.genfromtxt(
    archivo_ree,
    delimiter=",",
    names=True,
    encoding="utf-8",
    ndmin=1
)

N_rg = np.asarray(datos_rg["N"], dtype=float)
Rg2 = np.asarray(datos_rg["Rg2_media"], dtype=float)

N_ree = np.asarray(datos_ree["N"], dtype=float)
Ree2 = np.asarray(datos_ree["Ree2_media"], dtype=float)

# ============================================================
# 3. PREDICCIONES TEORICAS
# ============================================================

def teoria_rg2(N):
    return (N**2 - 1)/(6*N) * b2


def teoria_ree2(N):
    return (N - 1) * b2


# ============================================================
# 4. FUNCION PARA EL AJUSTE LINEAL
# ============================================================

def ajuste_lineal(N, valores):

    parametros, cov = np.polyfit(
        N,
        valores,
        1,
        cov=True
    )

    a, c = parametros

    error_a = np.sqrt(cov[0, 0])
    error_c = np.sqrt(cov[1, 1])

    prediccion = a*N + c

    ss_res = np.sum(
        (valores - prediccion)**2
    )

    ss_tot = np.sum(
        (valores - np.mean(valores))**2
    )

    R2 = 1 - ss_res/ss_tot

    return a, error_a, c, error_c, R2


# ============================================================
# 5. FUNCION PARA REPRESENTAR LOS AJUSTES
# ============================================================

def graficar_ajuste(
    N, valores, tipo,
    teoria, carpeta_salida
):

    orden = np.argsort(N)

    N = N[orden]
    valores = valores[orden]

    a, error_a, c, error_c, R2 = ajuste_lineal(
        N, valores
    )

    N_suave = np.linspace(
        N.min(),
        N.max(),
        500
    )

    valores_ajuste = a*N_suave + c
    valores_teoricos = teoria(N_suave)

    # --------------------------------------------------------
    # RADIO DE GIRO
    # --------------------------------------------------------

    if tipo == "RG":

        titulo = (
            r"Ajuste lineal del radio de giro al cuadrado "
            r"($k_e=1000$)"
        )

        ylabel = r"$\langle R_g^2\rangle$"

        etiqueta_datos = (
            r"Simulación: $\langle R_g^2\rangle$"
        )

        etiqueta_teoria = (
            r"Teoría: "
            r"$\langle R_g^2\rangle="
            r"\frac{N^2-1}{6N}b^2$"
        )

        archivo_salida = (
            carpeta_salida /
            "Ajuste_lineal_Rg2_ke1000_ULTRALARGA_Claudia.png"
        )

    # --------------------------------------------------------
    # DISTANCIA ENTRE EXTREMOS
    # --------------------------------------------------------

    elif tipo == "REE":

        titulo = (
            r"Ajuste lineal de la distancia entre extremos "
            r"al cuadrado ($k_e=1000$)"
        )

        ylabel = r"$\langle R_{ee}^2\rangle$"

        etiqueta_datos = (
            r"Simulación: $\langle R_{ee}^2\rangle$"
        )

        etiqueta_teoria = (
            r"Teoría: "
            r"$\langle R_{ee}^2\rangle=(N-1)b^2$"
        )

        archivo_salida = (
            carpeta_salida /
            "Ajuste_lineal_Ree2_ke1000_ULTRALARGA_Claudia.png"
        )

    else:
        raise ValueError("Tipo incorrecto: debe ser RG o REE")

    # ========================================================
    # GRAFICA
    # ========================================================

    fig, ax = plt.subplots(figsize=(9, 6))

    # Datos simulados
    ax.scatter(
        N,
        valores,
        s=65,
        color="tab:blue",
        zorder=4,
        label=etiqueta_datos
    )

    # Ajuste lineal
    ax.plot(
        N_suave,
        valores_ajuste,
        color="tab:red",
        linewidth=2,
        label=(
            rf"Ajuste: $y=({a:.5f})N"
            rf"{c:+.5f}$"
        )
    )

    # Prediccion teorica
    ax.plot(
        N_suave,
        valores_teoricos,
        color="tab:green",
        linestyle="--",
        linewidth=1.8,
        label=etiqueta_teoria
    )

    # Formato
    ax.set_title(titulo, fontsize=14)

    ax.set_xlabel(
        "Número de monómeros $N$",
        fontsize=12
    )

    ax.set_ylabel(ylabel, fontsize=12)

    ax.grid(alpha=0.3, linestyle=":")

    ax.legend(fontsize=10, loc="best")

    fig.tight_layout()

    fig.savefig(
        archivo_salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # RESULTADOS NUMERICOS
    # ========================================================

    print("\n" + "="*55)
    print(titulo)
    print("="*55)

    print(f"Pendiente: a = {a:.8f} ± {error_a:.8f}")
    print(f"Ordenada: c = {c:.8f} ± {error_c:.8f}")
    print(f"R² del ajuste = {R2:.8f}")

    if tipo == "RG":

        print(
            f"Pendiente teorica asintotica "
            f"b²/6 = {b2/6:.8f}"
        )

    else:

        print(
            f"Pendiente teorica exacta "
            f"b² = {b2:.8f}"
        )

        print(
            f"Ordenada teorica exacta "
            f"-b² = {-b2:.8f}"
        )

    return [
        tipo,
        a,
        error_a,
        c,
        error_c,
        R2
    ]


# ============================================================
# 6. EJECUTAR AJUSTES
# ============================================================

resultado_rg = graficar_ajuste(
    N_rg,
    Rg2,
    "RG",
    teoria_rg2,
    carpeta_rg
)

resultado_ree = graficar_ajuste(
    N_ree,
    Ree2,
    "REE",
    teoria_ree2,
    carpeta_ree
)

# ============================================================
# 7. GUARDAR RESUMEN DE LOS AJUSTES
# ============================================================

archivo_resumen = (
    base / "Py/plots/OBJETIVO_3/"
    "Resumen_ajustes_lineales_ke1000_ULTRALARGA_Claudia.csv"
)

import csv

with open(archivo_resumen, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "Radio",
        "Pendiente_a",
        "Error_a",
        "Ordenada_c",
        "Error_c",
        "R2"
    ])

    writer.writerow(resultado_rg)
    writer.writerow(resultado_ree)

print("\n" + "="*55)
print("AJUSTES LINEALES FINALIZADOS")
print("="*55)

print("Grafica Rg²:", carpeta_rg)
print("Grafica Ree²:", carpeta_ree)
print("Resumen:", archivo_resumen)
# %%


















# %%
# ============================================================
# OBJETIVO 3 - AJUSTES DE RADIOS SIN CUADRAR
# ke = 1000 | SIMULACIONES ULTRALARGAS
#
# MAGNITUDES:
#   1. sqrt(<Rg²>)
#   2. <Rg>
#   3. sqrt(<Ree²>)
#   4. <Ree>
#
# AJUSTES:
#   R = A*N^nu
#   ln(R) = ln(A) + nu*ln(N)
#
# AZUL: simulacion
# ROJO: ajuste potencial
# VERDE: prediccion teorica RMS
#
# SALIDAS:
#   RG/k1000/normal
#   REE/k1000/normal
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import csv

# ============================================================
# 1. RUTAS Y PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

carpeta_rg = base / "Py/plots/OBJETIVO_3/RG/k1000"
carpeta_ree = base / "Py/plots/OBJETIVO_3/REE/k1000"

salida_rg = carpeta_rg / "normal"
salida_ree = carpeta_ree / "normal"

salida_rg.mkdir(parents=True, exist_ok=True)
salida_ree.mkdir(parents=True, exist_ok=True)

ke = 1000.0
kBT = 1.0
l0 = 1.0

sigma2 = kBT / ke

# Longitud cuadratica media del enlace
b2 = (
    l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
) / (l0**2 + sigma2)

b = np.sqrt(b2)

print(f"b² teorico = {b2:.8f}")
print(f"b teorico  = {b:.8f}")

# ============================================================
# 2. LEER LOS CSV DE LOS RADIOS
# ============================================================

archivo_rg = (
    carpeta_rg / "cuadrado" /
    "Resumen_Rg_ke1000_ULTRALARGA_Claudia.csv"
)

archivo_ree = (
    carpeta_ree / "cuadrado" /
    "Resumen_Ree_ke1000_ULTRALARGA_Claudia.csv"
)

# Permitir tambien que los CSV esten en la carpeta k1000
if not archivo_rg.exists():
    archivo_rg = (
        carpeta_rg /
        "Resumen_Rg_ke1000_ULTRALARGA_Claudia.csv"
    )

if not archivo_ree.exists():
    archivo_ree = (
        carpeta_ree /
        "Resumen_Ree_ke1000_ULTRALARGA_Claudia.csv"
    )

if not archivo_rg.exists():
    raise FileNotFoundError(
        f"No encuentro el CSV de Rg:\n{archivo_rg}"
    )

if not archivo_ree.exists():
    raise FileNotFoundError(
        f"No encuentro el CSV de Ree:\n{archivo_ree}"
    )

datos_rg = np.genfromtxt(
    archivo_rg,
    delimiter=",",
    names=True,
    encoding="utf-8",
    ndmin=1
)

datos_ree = np.genfromtxt(
    archivo_ree,
    delimiter=",",
    names=True,
    encoding="utf-8",
    ndmin=1
)

# ============================================================
# 3. EXTRAER LAS CUATRO MAGNITUDES
# ============================================================

N_rg = np.asarray(datos_rg["N"], dtype=float)
N_ree = np.asarray(datos_ree["N"], dtype=float)

Rg_RMS = np.asarray(
    datos_rg["Rg_RMS"], dtype=float
)

Rg_medio = np.asarray(
    datos_rg["Rg_media"], dtype=float
)

Ree_RMS = np.asarray(
    datos_ree["Ree_RMS"], dtype=float
)

Ree_medio = np.asarray(
    datos_ree["Ree_media"], dtype=float
)

# ============================================================
# 4. PREDICCIONES TEORICAS
# ============================================================

def teoria_rg_rms(N):

    N = np.asarray(N, dtype=float)

    return b * np.sqrt(
        (N**2 - 1)/(6*N)
    )


def teoria_ree_rms(N):

    N = np.asarray(N, dtype=float)

    return b * np.sqrt(N - 1)


# ============================================================
# 5. AJUSTE POTENCIAL
# ============================================================

def ajustar_potencia(N, R):

    N = np.asarray(N, dtype=float)
    R = np.asarray(R, dtype=float)

    if np.any(N <= 0) or np.any(R <= 0):
        raise ValueError(
            "Los valores de N y R deben ser positivos."
        )

    # Linealizacion:
    # ln(R) = ln(A) + nu*ln(N)

    x = np.log(N)
    y = np.log(R)

    parametros, cov = np.polyfit(
        x, y, 1, cov=True
    )

    nu = parametros[0]
    lnA = parametros[1]

    A = np.exp(lnA)

    error_nu = np.sqrt(cov[0, 0])
    error_lnA = np.sqrt(cov[1, 1])

    error_A = A * error_lnA

    # Coeficiente de determinacion en espacio logaritmico
    y_ajuste = nu*x + lnA

    ss_res = np.sum((y - y_ajuste)**2)
    ss_tot = np.sum((y - np.mean(y))**2)

    R2_log = 1 - ss_res/ss_tot

    return A, error_A, nu, error_nu, R2_log


# ============================================================
# 6. FUNCION PARA GENERAR LOS DOS AJUSTES
# ============================================================

def graficar_ajustes(
    N, R, nombre, simbolo, tipo,
    carpeta_salida, es_rms
):

    N = np.asarray(N, dtype=float)
    R = np.asarray(R, dtype=float)

    # Ordenar los puntos por tamaño de cadena
    orden = np.argsort(N)

    N = N[orden]
    R = R[orden]

    # Ajuste potencial
    A, error_A, nu, error_nu, R2 = ajustar_potencia(
        N, R
    )

    N_suave = np.linspace(
        N.min(),
        N.max(),
        500
    )

    R_ajuste = A * N_suave**nu

    # --------------------------------------------------------
    # REFERENCIA TEORICA SEGUN EL RADIO
    # --------------------------------------------------------

    if tipo == "RG":

        R_teorico = teoria_rg_rms(N_suave)

        etiqueta_teoria = (
            r"$\sqrt{\langle R_g^2\rangle}_{\mathrm{teo}}"
            r"=b\sqrt{\frac{N^2-1}{6N}}$"
        )

    elif tipo == "REE":

        R_teorico = teoria_ree_rms(N_suave)

        etiqueta_teoria = (
            r"$\sqrt{\langle R_{ee}^2\rangle}_{\mathrm{teo}}"
            r"=b\sqrt{N-1}$"
        )

    else:
        raise ValueError(
            "El tipo debe ser 'RG' o 'REE'."
        )

    # En el promedio directo la linea verde es una
    # referencia RMS, no su propia prediccion exacta.
    if not es_rms:
        etiqueta_teoria += " (referencia RMS)"

    # Etiqueta del ajuste potencial
    etiqueta_ajuste = (
        rf"Ajuste: $R=({A:.4f})N^{{{nu:.4f}}}$"
    )

    # ========================================================
    # GRAFICA 1 - ESCALA NORMAL
    # ========================================================

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.scatter(
        N,
        R,
        color="tab:blue",
        s=60,
        zorder=4,
        label=simbolo
    )

    ax.plot(
        N_suave,
        R_ajuste,
        color="tab:red",
        linewidth=2,
        label=etiqueta_ajuste
    )

    ax.plot(
        N_suave,
        R_teorico,
        color="tab:green",
        linestyle="--",
        linewidth=1.8,
        label=etiqueta_teoria
    )

    ax.set_title(
        rf"{nombre} ($k_e=1000$)",
        fontsize=14
    )

    ax.set_xlabel(
        "Número de monómeros $N$",
        fontsize=12
    )

    ax.set_ylabel(
        simbolo,
        fontsize=12
    )

    ax.grid(
        alpha=0.3,
        linestyle=":"
    )

    ax.legend(
        loc="best",
        fontsize=9
    )

    fig.tight_layout()

    nombre_archivo = (
        f"{tipo}_{'RMS' if es_rms else 'MEDIA'}"
        "_ajuste_potencial_ke1000_ULTRALARGA_Claudia.png"
    )

    fig.savefig(
        carpeta_salida / nombre_archivo,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # GRAFICA 2 - ESCALA LOG-LOG
    # ========================================================

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.scatter(
        N,
        R,
        color="tab:blue",
        s=60,
        zorder=4,
        label=simbolo
    )

    ax.plot(
        N_suave,
        R_ajuste,
        color="tab:red",
        linewidth=2,
        label=(
            rf"Ajuste: $\ln R=\ln({A:.4f})"
            rf"+({nu:.4f}\pm{error_nu:.4f})\ln N$"
        )
    )

    ax.plot(
        N_suave,
        R_teorico,
        color="tab:green",
        linestyle="--",
        linewidth=1.8,
        label=etiqueta_teoria
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_title(
        rf"{nombre} - escala log-log ($k_e=1000$)",
        fontsize=14
    )

    ax.set_xlabel(
        r"$N$",
        fontsize=12
    )

    ax.set_ylabel(
        simbolo,
        fontsize=12
    )

    ax.grid(
        which="both",
        alpha=0.3,
        linestyle=":"
    )

    ax.legend(
        loc="best",
        fontsize=9
    )

    fig.tight_layout()

    nombre_archivo_log = (
        f"{tipo}_{'RMS' if es_rms else 'MEDIA'}"
        "_ajuste_loglog_ke1000_ULTRALARGA_Claudia.png"
    )

    fig.savefig(
        carpeta_salida / nombre_archivo_log,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    # ========================================================
    # RESULTADOS NUMERICOS
    # ========================================================

    print("\n" + "="*60)
    print(nombre)
    print("="*60)

    print(f"A = {A:.8f} ± {error_A:.8f}")
    print(f"nu = {nu:.8f} ± {error_nu:.8f}")
    print(f"R² log-log = {R2:.8f}")
    print("nu teorico asintotico = 0.50000000")

    return [
        nombre,
        A,
        error_A,
        nu,
        error_nu,
        R2
    ]


# ============================================================
# 7. AJUSTES DEL RADIO DE GIRO
# ============================================================

resultados = []

# Raiz del promedio del cuadrado
resultados.append(
    graficar_ajustes(
        N_rg,
        Rg_RMS,
        "Raíz del promedio cuadrático del radio de giro",
        r"$\sqrt{\langle R_g^2\rangle}$",
        "RG",
        salida_rg,
        es_rms=True
    )
)

# Promedio directo
resultados.append(
    graficar_ajustes(
        N_rg,
        Rg_medio,
        "Promedio directo del radio de giro",
        r"$\langle R_g\rangle$",
        "RG",
        salida_rg,
        es_rms=False
    )
)


# ============================================================
# 8. AJUSTES DE LA DISTANCIA ENTRE EXTREMOS
# ============================================================

# Raiz del promedio del cuadrado
resultados.append(
    graficar_ajustes(
        N_ree,
        Ree_RMS,
        "Raíz del promedio cuadrático entre extremos",
        r"$\sqrt{\langle R_{ee}^2\rangle}$",
        "REE",
        salida_ree,
        es_rms=True
    )
)

# Promedio directo
resultados.append(
    graficar_ajustes(
        N_ree,
        Ree_medio,
        "Promedio directo de la distancia entre extremos",
        r"$\langle R_{ee}\rangle$",
        "REE",
        salida_ree,
        es_rms=False
    )
)


# ============================================================
# 9. GUARDAR RESUMEN DE LOS CUATRO AJUSTES
# ============================================================

archivo_resultados = (
    base / "Py/plots/OBJETIVO_3/"
    "Ajustes_radios_normales_ke1000_ULTRALARGA_Claudia.csv"
)

with open(archivo_resultados, "w", newline="") as f:

    writer = csv.writer(f)

    writer.writerow([
        "Magnitud",
        "A",
        "error_A",
        "nu",
        "error_nu",
        "R2_log"
    ])

    writer.writerows(resultados)

print("\n" + "="*60)
print("AJUSTES DE RADIOS SIN CUADRAR - ke=1000")
print("="*60)

print(f"Gráficas generadas: {2 * len(resultados)}")

print("\nRadio de giro:")
print(salida_rg)

print("\nDistancia entre extremos:")
print(salida_ree)

print("\nResumen de ajustes:")
print(archivo_resultados)

# %%
















# %%
# ============================================================
# OBJETIVO 3 - COMPARACION DE RADIOS ENTRE DISTINTOS ke
#
# Modelos:
#   ke = 100
#   ke = 1000
#   Enlaces rigidos (ke -> infinito)
#
# Comparaciones:
#   1. <Rg²> frente a N
#   2. sqrt(<Rg²>) frente a N
#   3. <Rg> frente a N
#   4. <Ree²> frente a N
#   5. sqrt(<Ree²>) frente a N
#   6. <Ree> frente a N
#
# Salidas:
#   OBJETIVO_3/RG/
#   OBJETIVO_3/REE/
#
# No se modifican los datos originales.
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. CARPETAS Y PARAMETROS
# ============================================================

base = Path.home() / "Documents/A_5-Molecular_dynamics"

objetivo = base / "Py/plots/OBJETIVO_3"

salida_rg = objetivo / "RG"
salida_ree = objetivo / "REE"

salida_rg.mkdir(parents=True, exist_ok=True)
salida_ree.mkdir(parents=True, exist_ok=True)

kBT = 1.0
l0 = 1.0

# ============================================================
# 2. BUSQUEDA DE CSV
# ============================================================

def buscar_csv(tipo, ke):

    etiqueta = (
        "SUPERLARGA" if ke == 100
        else "ULTRALARGA"
    )

    nombre = (
        f"Resumen_{tipo}_ke{ke}_"
        f"{etiqueta}_Claudia.csv"
    )

    raiz = objetivo / tipo / f"k{ke}"

    candidatos = [
        raiz / "cuadrado" / nombre,
        raiz / nombre,
        raiz / "normal" / nombre
    ]

    for archivo in candidatos:
        if archivo.exists():
            print("Archivo encontrado:", archivo)
            return archivo

    raise FileNotFoundError(
        f"No se encuentra el CSV de {tipo}, ke={ke}.\n"
        + "\n".join(str(p) for p in candidatos)
    )


def leer_csv(tipo, ke):

    archivo = buscar_csv(tipo, ke)

    datos = np.genfromtxt(
        archivo,
        delimiter=",",
        names=True,
        encoding="utf-8",
        ndmin=1
    )

    return {
        "N": np.asarray(datos["N"], dtype=float),
        "media": np.asarray(
            datos[f"{'Rg' if tipo == 'RG' else 'Ree'}_media"],
            dtype=float
        ),
        "media2": np.asarray(
            datos[f"{'Rg' if tipo == 'RG' else 'Ree'}2_media"],
            dtype=float
        )
    }


# ============================================================
# 3. CARGAR LAS CUATRO SERIES
# ============================================================

rg100 = leer_csv("RG", 100)
rg1000 = leer_csv("RG", 1000)

ree100 = leer_csv("REE", 100)
ree1000 = leer_csv("REE", 1000)

# ============================================================
# 4. LONGITUD CUADRATICA MEDIA DEL ENLACE
# ============================================================

def b2_elastico(ke):

    sigma2 = kBT / ke

    return (
        l0**4 + 6*l0**2*sigma2 + 3*sigma2**2
    ) / (l0**2 + sigma2)


b2_100 = b2_elastico(100)
b2_1000 = b2_elastico(1000)
b2_rigido = l0**2

print("\nLONGITUD CUADRATICA MEDIA DEL ENLACE")
print(f"ke=100:  b²={b2_100:.8f}")
print(f"ke=1000: b²={b2_1000:.8f}")
print(f"Rigido:  b²={b2_rigido:.8f}")


# ============================================================
# 5. EXPRESIONES TEORICAS
# ============================================================

def Rg2_teoria(N, b2):

    return ((N**2 - 1)/(6*N))*b2


def Ree2_teoria(N, b2):

    return (N - 1)*b2


def teoria_radio(N, b2, tipo, magnitud):

    if tipo == "RG":
        radio2 = Rg2_teoria(N, b2)

    elif tipo == "REE":
        radio2 = Ree2_teoria(N, b2)

    else:
        raise ValueError("Radio no reconocido")

    if magnitud == "cuadrado":
        return radio2

    if magnitud == "raiz":
        return np.sqrt(radio2)

    raise ValueError("Magnitud no reconocida")


# ============================================================
# 6. FUNCION PARA COMPARAR MODELOS
# ============================================================

def graficar_comparacion(
    tipo, magnitud, datos100, datos1000,
    carpeta_salida
):

    # --------------------------------------------------------
    # DATOS DE LAS SIMULACIONES
    # --------------------------------------------------------

    N100 = datos100["N"]
    N1000 = datos1000["N"]

    if magnitud == "cuadrado":

        y100 = datos100["media2"]
        y1000 = datos1000["media2"]

    elif magnitud == "raiz":

        y100 = np.sqrt(datos100["media2"])
        y1000 = np.sqrt(datos1000["media2"])

    elif magnitud == "media":

        y100 = datos100["media"]
        y1000 = datos1000["media"]

    else:
        raise ValueError("Magnitud no valida")

    # Ordenar los tamaños
    orden100 = np.argsort(N100)
    orden1000 = np.argsort(N1000)

    N100 = N100[orden100]
    y100 = y100[orden100]

    N1000 = N1000[orden1000]
    y1000 = y1000[orden1000]

    N_min = min(N100.min(), N1000.min())
    N_max = max(N100.max(), N1000.max())

    N_suave = np.linspace(N_min, N_max, 600)

    # --------------------------------------------------------
    # NOMBRES Y NOTACION
    # --------------------------------------------------------

    if tipo == "RG":

        radio = "Rg"
        titulo_radio = "radio de giro"

        if magnitud == "cuadrado":
            ylabel = r"$\langle R_g^2\rangle$"
            titulo = "Radio de giro cuadrático medio"

        elif magnitud == "raiz":
            ylabel = r"$\sqrt{\langle R_g^2\rangle}$"
            titulo = "Raíz cuadrática media del radio de giro"

        else:
            ylabel = r"$\langle R_g\rangle$"
            titulo = "Radio de giro medio"

        formula = (
            r"$\frac{N^2-1}{6N}$"
            if magnitud == "cuadrado"
            else r"$\sqrt{\frac{N^2-1}{6N}}$"
        )

    elif tipo == "REE":

        radio = "Ree"
        titulo_radio = "distancia entre extremos"

        if magnitud == "cuadrado":
            ylabel = r"$\langle R_{ee}^2\rangle$"
            titulo = "Distancia cuadrática media entre extremos"

        elif magnitud == "raiz":
            ylabel = r"$\sqrt{\langle R_{ee}^2\rangle}$"
            titulo = "Raíz cuadrática media entre extremos"

        else:
            ylabel = r"$\langle R_{ee}\rangle$"
            titulo = "Distancia media entre extremos"

        formula = (
            r"$(N-1)$"
            if magnitud == "cuadrado"
            else r"$\sqrt{N-1}$"
        )

    else:
        raise ValueError("Tipo no valido")

    # ========================================================
    # GRAFICA
    # ========================================================

    fig, ax = plt.subplots(figsize=(10, 6.5))


    # --------------------------------------------------------
    # SIMULACIONES: SOLO PUNTOS, SIN LINEAS
    # --------------------------------------------------------

    ax.scatter(
        N100,
        y100,
        color="tab:blue",
        marker="o",
        s=60,
        zorder=4,
        label=r"Simulación $k_e=100$"
    )

    ax.scatter(
        N1000,
        y1000,
        color="tab:red",
        marker="s",
        s=60,
        zorder=4,
        label=r"Simulación $k_e=1000$"
    )

    # --------------------------------------------------------
    # PREDICCIONES TEORICAS
    # --------------------------------------------------------

    # Las formulas RMS y cuadraticas se aplican
    # unicamente a las magnitudes correspondientes.
    #
    # Para los promedios directos no superponemos
    # predicciones RMS, porque no son la misma magnitud.

    if magnitud != "media":

        teo100 = teoria_radio(
            N_suave, b2_100, tipo, magnitud
        )

        teo1000 = teoria_radio(
            N_suave, b2_1000, tipo, magnitud
        )

        teo_rigido = teoria_radio(
            N_suave, b2_rigido, tipo, magnitud
        )

        ax.plot(
            N_suave,
            teo100,
            color="tab:blue",
            linestyle="--",
            linewidth=1.8,
            alpha=0.85,
            label=r"Teoría elástica $k_e=100$"
        )

        ax.plot(
            N_suave,
            teo1000,
            color="tab:red",
            linestyle="--",
            linewidth=1.8,
            alpha=0.85,
            label=r"Teoría elástica $k_e=1000$"
        )

        ax.plot(
            N_suave,
            teo_rigido,
            color="tab:green",
            linestyle="-.",
            linewidth=2.3,
            label="Modelo rígido: " + formula
        )

    # --------------------------------------------------------
    # FORMATO
    # --------------------------------------------------------

    ax.set_title(
        titulo + r" ($k_e=100,\ 1000$)",
        fontsize=14
    )

    ax.set_xlabel(
        "Número de monómeros $N$",
        fontsize=12
    )

    ax.set_ylabel(ylabel, fontsize=12)

    ax.grid(alpha=0.3, linestyle=":")

    ax.legend(
        loc="best",
        fontsize=9
    )

    fig.tight_layout()

    # --------------------------------------------------------
    # GUARDAR DIRECTAMENTE EN RG O REE
    # --------------------------------------------------------

    nombre_archivo = (
        f"Comparacion_{radio}_{magnitud}_"
        "ke100_ke1000_rigido_Claudia.png"
    )

    ruta_salida = carpeta_salida / nombre_archivo

    fig.savefig(
        ruta_salida,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig)

    print("Guardada:", ruta_salida)


# ============================================================
# 7. COMPARACIONES DEL RADIO DE GIRO
# ============================================================

graficar_comparacion(
    "RG", "cuadrado",
    rg100, rg1000, salida_rg
)

graficar_comparacion(
    "RG", "raiz",
    rg100, rg1000, salida_rg
)

graficar_comparacion(
    "RG", "media",
    rg100, rg1000, salida_rg
)


# ============================================================
# 8. COMPARACIONES DE DISTANCIA ENTRE EXTREMOS
# ============================================================

graficar_comparacion(
    "REE", "cuadrado",
    ree100, ree1000, salida_ree
)

graficar_comparacion(
    "REE", "raiz",
    ree100, ree1000, salida_ree
)

graficar_comparacion(
    "REE", "media",
    ree100, ree1000, salida_ree
)


# ============================================================
# 9. FIN
# ============================================================

print("\n" + "="*60)
print("COMPARACIONES ENTRE MODELOS FINALIZADAS")
print("="*60)

print("Se han generado 6 graficas.")
print("Radio de giro:", salida_rg)
print("Distancia entre extremos:", salida_ree)

# %%

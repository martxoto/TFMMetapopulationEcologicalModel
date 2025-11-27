import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Configuración Estética ===
# Usamos un estilo limpio para publicaciones académicas
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 10, 'figure.titlesize': 16})

# === 1. Cargar Datos ===
filename = 'results.txt'
print(f"Leyendo datos de {filename}...")

try:
    # Leemos el archivo generado por C++
    # Asumimos columnas: K, Robustez, SurvPlants, SurvInsects, Service, GiniP, GiniI
    df = pd.read_csv(filename, sep='\s+', comment='#', header=None,
                     names=['K', 'Robustez', 'SurvPlants', 'SurvInsects', 'Service', 'GiniP', 'GiniI'])
except FileNotFoundError:
    print(f"ERROR: No se encuentra '{filename}'. Ejecuta primero la simulación en C++.")
    exit()

# === 2. Procesar Datos ===
# Normalizamos el Eje X: 0.0 = Inicio, 1.0 = Todas las plantas eliminadas
max_k = df['K'].max()
df['Fraccion_Eliminada'] = df['K'] / max_k

# === 3. Crear Panel de Gráficas ===
# Creamos una figura con 3 sub-gráficas (columnas)
fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

# --- GRÁFICA A: ROBUSTEZ Y SUPERVIVENCIA ---
ax1 = axes[0]
# Línea principal: Robustez Total
ax1.plot(df['Fraccion_Eliminada'], df['Robustez'], color='navy', linewidth=2.5, label='Robustez Total')
ax1.fill_between(df['Fraccion_Eliminada'], df['Robustez'], color='navy', alpha=0.1)

# Eje secundario para conteos absolutos
ax1b = ax1.twinx()
ax1b.plot(df['Fraccion_Eliminada'], df['SurvPlants'], color='forestgreen', linestyle='--', label='Plantas Vivas')
ax1b.plot(df['Fraccion_Eliminada'], df['SurvInsects'], color='darkorange', linestyle='-.', label='Insectos Vivos')

# Etiquetas y Estilo
ax1.set_title("Colapso del Ecosistema")
ax1.set_xlabel("Fracción de Plantas Eliminadas")
ax1.set_ylabel("Fracción de Robustez (0-1)", color='navy')
ax1b.set_ylabel("Número de Especies", color='gray')
ax1.set_ylim(-0.05, 1.05)

# Leyenda combinada
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1b.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower left', fontsize=9)


# --- GRÁFICA B: SERVICIO DE POLINIZACIÓN ---
ax2 = axes[1]
# Asumiendo que el servicio es la biomasa total
ax2.plot(df['Fraccion_Eliminada'], df['Service'], color='darkgoldenrod', linewidth=2)
ax2.fill_between(df['Fraccion_Eliminada'], df['Service'], color='gold', alpha=0.3)

ax2.set_title("Servicio de Polinización (Función)")
ax2.set_xlabel("Fracción de Plantas Eliminadas")
ax2.set_ylabel("Biomasa Total de Polinizadores")
ax2.grid(True, linestyle=':', alpha=0.6)


# --- GRÁFICA C: DIVERSIDAD (Gini-Simpson) ---
ax3 = axes[2]
ax3.plot(df['Fraccion_Eliminada'], df['GiniP'], color='forestgreen', linewidth=2, label='Diversidad Plantas')
ax3.plot(df['Fraccion_Eliminada'], df['GiniI'], color='darkorange', linewidth=2, label='Diversidad Insectos')

ax3.set_title("Pérdida de Diversidad")
ax3.set_xlabel("Fracción de Plantas Eliminadas")
ax3.set_ylabel("Índice Gini-Simpson (1-D)")
ax3.set_ylim(-0.05, 1.05)
ax3.legend(loc='lower left')
ax3.grid(True, linestyle=':', alpha=0.6)


# === 4. Guardar y Mostrar ===
output_name = 'panel_resultados_tfm.png'
plt.savefig(output_name, dpi=300)
print(f"Gráfica guardada exitosamente como '{output_name}'")

plt.show()

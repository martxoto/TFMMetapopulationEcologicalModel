import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns # Opcional, para estilos más bonitos

# Configuración de estilo (opcional, hace que las gráficas se vean más académicas)
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 12})

def plot_robustness(filename='results.txt'):
    print(f"Leyendo datos de {filename}...")
    
    try:
        # Leemos el archivo. 
        # Asumimos que la primera línea empieza por # y es el header, o la saltamos.
        # Usamos sep='\s+' para que entienda espacios múltiples como separador.
        # Nombres de columnas basados en tu código C++: K, Robustez, SurvPlants, SurvInsects
        df = pd.read_csv(filename, sep='\s+', comment='#', header=None, 
                         names=['K', 'Robustez', 'SurvPlants', 'SurvInsects'])
    except FileNotFoundError:
        print(f"Error: No se encuentra el archivo {filename}.")
        print("Asegúrate de haber ejecutado la simulación en C++ primero.")
        return

    # Normalizamos el eje X (Fracción de plantas eliminadas)
    # K=0 es 0% eliminado. K=Total es 100% eliminado.
    max_k = df['K'].max()
    df['Fraction_Removed'] = df['K'] / max_k

    # --- CREAR LA GRÁFICA ---
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # 1. Eje Y Principal (Izquierda): Robustez (Fracción Total)
    color = 'tab:blue'
    ax1.set_xlabel('Fracción de Plantas Eliminadas (Rareza: Rara $\\to$ Común)')
    ax1.set_ylabel('Robustez del Ecosistema\n(Fracción de especies vivas)', color=color, fontsize=12)
    
    # Línea de Robustez
    line1 = ax1.plot(df['Fraction_Removed'], df['Robustez'], color=color, linewidth=2.5, label='Robustez Total')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim(-0.05, 1.05) # Márgenes para que se vea bien el 0 y el 1

    # Área sombreada bajo la curva (R50 o AUC visual)
    ax1.fill_between(df['Fraction_Removed'], df['Robustez'], color=color, alpha=0.1)

    # 2. Eje Y Secundario (Derecha): Número absoluto de especies
    ax2 = ax1.twinx()  
    
    # Línea de Plantas
    color_plants = 'tab:green'
    line2 = ax2.plot(df['Fraction_Removed'], df['SurvPlants'], color=color_plants, linestyle='--', label='Plantas Vivas')
    
    # Línea de Insectos
    color_insects = 'tab:orange'
    line3 = ax2.plot(df['Fraction_Removed'], df['SurvInsects'], color=color_insects, linestyle='-.', label='Insectos Vivos')
    
    ax2.set_ylabel('Número de Especies Supervivientes', color='gray', fontsize=12)
    ax2.tick_params(axis='y', labelcolor='gray')

    # --- DETALLES FINALES ---
    
    # Título informativo
    plt.title(f'Curva de Tolerancia a la Extinción\n(Total Pasos K={max_k})', fontsize=14, fontweight='bold')
    
    # Leyenda combinada
    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='lower left', frameon=True, fancybox=True, framealpha=0.9)

    plt.tight_layout()
    
    # Guardar y Mostrar
    output_img = 'grafica_robustez.png'
    plt.savefig(output_img, dpi=300)
    print(f"Gráfica guardada como '{output_img}'")
    plt.show()

if __name__ == "__main__":
    plot_robustness()
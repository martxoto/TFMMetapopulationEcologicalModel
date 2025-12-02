import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import time

# === CONFIGURACIÓN ===
EXECUTABLE = './exp1.exe'  # Asegúrate de que este es el nombre de tu exe
STEP_H = "0.01"                    # El paso de integración que usas siempre
METRIC_FILE = "robustnessD.txt" # El archivo donde guardaste la R final

# Valores de D a probar (Escala logarítmica suele ser mejor para barridos amplios)
# Probamos desde aislamiento (0) hasta conectividad extrema (100)
D_values = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0]

results_R = []

print(f"Iniciando barrido de D para {len(D_values)} valores...")
print("-" * 40)
print(f"{'D':<10} | {'R (Robustez)':<15}")
print("-" * 40)

for D in D_values:
    # 1. Preparar la entrada para el programa C++
    # El programa espera: h (enter) D (enter)
    input_str = f"{STEP_H}\n{D}\n"
    
    # 2. Ejecutar C++ y pasarle los datos
    try:
        # subprocess.run lanza el programa, le pasa el input y espera que acabe
        process = subprocess.run(
            [EXECUTABLE], 
            input=input_str, 
            text=True, 
            capture_output=True, # Capturamos lo que imprime por pantalla para que no ensucie
            check=True
        )
        
        # 3. Leer el resultado
        if os.path.exists(METRIC_FILE):
            with open(METRIC_FILE, 'r') as f:
                r_val = float(f.read().strip())
                results_R.append(r_val)
                print(f"{D:<10} | {r_val:.5f}")
        else:
            print(f"Error: No se generó {METRIC_FILE} para D={D}")
            results_R.append(0)
            
    except subprocess.CalledProcessError as e:
        print(f"Error ejecutando C++ para D={D}: {e}")
        results_R.append(0)

# === GRAFICAR ===
sns.set_theme(style="whitegrid")
plt.figure(figsize=(10, 6))

# Pintamos línea y puntos
plt.plot(D_values, results_R, 'o-', color='navy', linewidth=2, markersize=8)

# Destacar el máximo
max_R = max(results_R)
best_D = D_values[results_R.index(max_R)]
plt.plot(best_D, max_R, 'ro', markersize=12, label=f'Óptimo (D={best_D}, R={max_R:.3f})')

plt.title("Sensibilidad a la Dispersión", fontsize=14, fontweight='bold')
plt.xlabel("Tasa de Dispersión ($D$)", fontsize=12)
plt.ylabel("Robustez del Ecosistema ($R$)", fontsize=12)
plt.xscale('symlog', linthresh=0.1) # Escala especial: logarítmica pero permite el 0
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()

plt.tight_layout()
plt.savefig('optimizacion_D.png', dpi=300)
print("-" * 40)
print(f"Gráfica guardada como 'optimizacion_D.png'. Mejor D = {best_D}")
plt.show()
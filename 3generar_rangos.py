import pandas as pd
import matplotlib.pyplot as plt
from Bio import AlignIO
import numpy as np
import os
import time

# --- CONFIGURACIÓN GIT-FRIENDLY ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)

# NOTA: Cambia "p" por "c" si estás trabajando con Clases, igual que en el script 2
NIVEL_TAXONOMICO = "p" 

# Rutas Dinámicas (Coinciden con lo que generó el Script 2)
ARCHIVO_ALINEADO = os.path.join(BASE_DIR, "data", f"consensos_{NIVEL_TAXONOMICO}_alineados.aln")
CARPETA_RESULTADOS = os.path.join(BASE_DIR, "resultados")
EXCEL_SALIDA = os.path.join(CARPETA_RESULTADOS, f"reporte_barcodes_{NIVEL_TAXONOMICO}.xlsx")
GRAFICO_SALIDA = os.path.join(CARPETA_RESULTADOS, f"grafico_unicidad_{NIVEL_TAXONOMICO}.png")

# Parámetros Científicos
TAMANO_VENTANA = 100    # Tamaño del barcode (bp)
MINIMO_DIFERENCIAS = 25 # Umbral de seguridad (letras distintas)

def adn_to_int_matrix(alignment):
    """
    Convierte el alineamiento a una matriz numérica para NumPy.
    A=1, C=2, G=3, T=4, Gap/N=0
    """
    mapeo = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0, 'N': 0, 'a':1, 'c':2, 'g':3, 't':4}
    n_seqs = len(alignment)
    largo = alignment.get_alignment_length()
    
    matriz = np.zeros((n_seqs, largo), dtype=np.int8)
    
    for i, record in enumerate(alignment):
        matriz[i, :] = [mapeo.get(c, 0) for c in str(record.seq)]
        
    return matriz, [r.id for r in alignment], [str(r.seq) for r in alignment]

def buscar_ventanas_vectorizadas(alignment, ventana, umbral_diff):
    print("⏳ Digitalizando secuencias (Matriz NumPy)...")
    matriz_gen, nombres, seqs_txt = adn_to_int_matrix(alignment)
    n_seqs, largo_gen = matriz_gen.shape
    
    hallazgos = []
    # Matriz para gráfico: (N_Phyla, Posiciones)
    perfil_unicidad = np.zeros((n_seqs, largo_gen - ventana), dtype=np.int16) 

    print(f"🚀 Analizando {largo_gen - ventana} ventanas de {ventana} bp...")
    print(f"   Matriz de procesamiento: {n_seqs}x{n_seqs} comparaciones simultáneas.")
    start_time = time.time()

    # --- SLIDING WINDOW (VECTORIZADO) ---
    for i in range(largo_gen - ventana):
        
        # 1. Extraer Ventana: Matriz (N, Ventana)
        slice_ventana = matriz_gen[:, i : i + ventana]
        
        # 2. Filtro de Calidad (Gaps)
        conteo_gaps = np.sum(slice_ventana == 0, axis=1)
        es_valido = conteo_gaps <= (ventana * 0.5)
        
        if not np.any(es_valido): continue
            
        # 3. MÁGICA DE NUMPY: Broadcasting 3D
        # Comparamos todas las filas contra todas las filas instantáneamente
        matriz_A = slice_ventana[:, np.newaxis, :] # (N, 1, W)
        matriz_B = slice_ventana[np.newaxis, :, :] # (1, N, W)
        
        # Matriz de distancias (N, N)
        matriz_distancias = np.sum(matriz_A != matriz_B, axis=2)
        
        # Ajuste: Distancia conmigo mismo es Infinito (para que no sea el mínimo)
        np.fill_diagonal(matriz_distancias, 9999)
        
        # 4. Encontrar el vecino más cercano (Mínimo por fila)
        vecino_mas_cercano = np.min(matriz_distancias, axis=1)
        
        # 5. Guardar datos para gráfico
        perfil_unicidad[:, i] = np.where(es_valido, vecino_mas_cercano, 0)
        
        # 6. Detectar Hallazgos
        indices_ganadores = np.where(es_valido & (vecino_mas_cercano >= umbral_diff))[0]
        
        for idx in indices_ganadores:
            dist = vecino_mas_cercano[idx]
            hallazgos.append({
                'Phylum': nombres[idx],
                'Inicio': i + 1,
                'Fin': i + ventana,
                'Diferencias_Minimas': int(dist),
                'Secuencia_Rango': seqs_txt[idx][i : i + ventana]
            })

        if i % 200 == 0:
            print(f"   ... {i}/{largo_gen} bp", end='\r')

    print(f"\n✨ Completado en {round(time.time() - start_time, 2)} segundos.")
    return hallazgos, perfil_unicidad, nombres

# --- INICIO DEL PROGRAMA ---
if __name__ == "__main__":
    print(f"--- 03. BUSCADOR DE BARCODES ({NIVEL_TAXONOMICO}) ---")

    # Crear carpeta de resultados si no existe
    if not os.path.exists(CARPETA_RESULTADOS):
        os.makedirs(CARPETA_RESULTADOS)

    if not os.path.exists(ARCHIVO_ALINEADO):
        print(f"❌ Error: No encuentro '{ARCHIVO_ALINEADO}'")
        print("   Ejecuta el paso 2 (Generar Consensos) primero.")
        exit()

    try:
        print("📂 Cargando alineamiento...")
        align = AlignIO.read(ARCHIVO_ALINEADO, "clustal")
        
        datos, matriz_grafico, nombres = buscar_ventanas_vectorizadas(align, TAMANO_VENTANA, MINIMO_DIFERENCIAS)
        
        # 1. GENERAR EXCEL
        if datos:
            df = pd.DataFrame(datos)
            df = df.sort_values(by=['Diferencias_Minimas', 'Phylum'], ascending=[False, True])
            
            # Eliminamos duplicados adyacentes (opcional, limpieza visual)
            # df = df.drop_duplicates(subset=['Phylum', 'Diferencias_Minimas'])
            
            df.to_excel(EXCEL_SALIDA, index=False)
            print(f"\n✅ REPORTE GENERADO: {EXCEL_SALIDA}")
            print(f"   Hallazgos totales: {len(datos)}")
            if not df.empty:
                print(f"   🥇 Mejor hallazgo: {df.iloc[0]['Phylum']} ({df.iloc[0]['Diferencias_Minimas']} diferencias)")
        else:
            print("\n⚠️ No se encontraron rangos diagnósticos con estos parámetros.")

        # 2. GENERAR GRÁFICO
        print("📊 Generando gráfico de unicidad...")
        plt.figure(figsize=(15, 8))
        
        # Graficar solo el Top 10 para no saturar
        maximos = np.max(matriz_grafico, axis=1)
        indices_top = np.argsort(maximos)[-10:] 
        
        for idx in indices_top:
            plt.plot(matriz_grafico[idx], label=nombres[idx], alpha=0.8, linewidth=1.5)
            
        plt.title(f"Perfil de Unicidad (Ventana {TAMANO_VENTANA}bp) - Top 10 Grupos")
        plt.xlabel("Posición de Inicio")
        plt.ylabel("Diferencias mínimas vs. el resto")
        plt.axhline(y=MINIMO_DIFERENCIAS, color='r', linestyle='--', label='Umbral')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        plt.savefig(GRAFICO_SALIDA)
        print(f"✅ Gráfico guardado en: {GRAFICO_SALIDA}")
        # plt.show() # Descomenta si quieres verlo al momento

    except Exception as e:
        print(f"\n❌ Ocurrió un error inesperado: {e}")
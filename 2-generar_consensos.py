import os
import random
import subprocess
from collections import Counter
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- CONFIGURACIÓN DE RUTAS RELATIVAS (GIT FRIENDLY) ---
# 1. Detectamos dónde está este script ahora mismo
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# 2. Subimos un nivel para llegar a la carpeta raíz 'problema_1'
BASE_DIR = os.path.dirname(SCRIPT_DIR)

# 3. Definimos las rutas usando la base dinámica
# NOTA: Si en el futuro usas Clases, cambia "secuencias_por_phylum" por "secuencias_por_c"
CARPETA_INPUT = os.path.join(BASE_DIR, "secuencias_por_phylum")
MUSCLE_EXE = os.path.join(BASE_DIR, "muscle.exe")
ARCHIVO_SALIDA = os.path.join(BASE_DIR, "todos_los_consensos_robustos.fasta")

# --- PARÁMETROS CIENTÍFICOS ---
MUESTRA_PARA_CONSENSO = 50 
UMBRAL_CONSENSO = 0.6  # 60% de acuerdo necesario para evitar 'N'

def obtener_consenso_frecuencial(alignment, umbral):
    """
    Calcula el consenso basándose en la distribución de frecuencias.
    Si ninguna base supera el umbral, coloca una 'N'.
    """
    consenso_str = ""
    largo = alignment.get_alignment_length()
    
    for i in range(largo):
        columna = alignment[:, i]
        
        # Filtramos gaps '-' y Ns existentes
        bases_reales = [b for b in columna if b not in ['-', 'N']]
        
        if not bases_reales:
            consenso_str += "-" # Si todo era gap, mantenemos gap
            continue

        total_bases = len(bases_reales)
        conteo = Counter(bases_reales)
        
        # Obtenemos la base más común
        base_top, cantidad = conteo.most_common(1)[0]
        
        # Calculamos frecuencia
        frecuencia = cantidad / total_bases
        
        if frecuencia >= umbral:
            consenso_str += base_top
        else:
            consenso_str += "N" # Ambigüedad/Variabilidad
        
    return Seq(consenso_str)

# --- INICIO DEL PROGRAMA ---
if __name__ == "__main__":
    print(f"--- GENERADOR DE CONSENSOS (Frecuencia > {int(UMBRAL_CONSENSO*100)}%) ---")
    
    # Validaciones de seguridad para Git
    if not os.path.exists(MUSCLE_EXE):
        print(f"❌ ERROR: No encuentro 'muscle.exe' en la carpeta raíz:\n   {BASE_DIR}")
        print("   Por favor, descarga muscle y ponlo ahí.")
        exit()
        
    if not os.path.exists(CARPETA_INPUT):
        print(f"❌ ERROR: No encuentro la carpeta de datos:\n   {CARPETA_INPUT}")
        print("   Ejecuta el script 1_separar_secuencias.py primero.")
        exit()

    lista_consensos = []
    archivos = [f for f in os.listdir(CARPETA_INPUT) if f.endswith(".fasta")]

    print(f"Procesando {len(archivos)} grupos taxonómicos...")

    for archivo in archivos:
        grupo_name = archivo.replace(".fasta", "")
        ruta_input = os.path.join(CARPETA_INPUT, archivo)
        
        # Archivos temporales (se sobrescriben, no importa el nombre)
        temp_fasta = "temp_consenso.fasta"
        temp_aln = "temp_consenso.aln"
        
        try:
            seqs = list(SeqIO.parse(ruta_input, "fasta"))
            
            # Caso 1: Muy pocos datos
            if len(seqs) < 3:
                print(f"⚠️ {grupo_name}: Pocos datos (<3). Consenso directo.")
                registro = SeqRecord(seqs[0].seq, id=grupo_name, description="Representante_Directo")
                lista_consensos.append(registro)
                continue
            
            # Caso 2: Submuestreo y Alineamiento
            seleccion = random.sample(seqs, min(len(seqs), MUESTRA_PARA_CONSENSO))
            SeqIO.write(seleccion, temp_fasta, "fasta")
            
            # Llamada a MUSCLE
            cmd = [MUSCLE_EXE, "-in", temp_fasta, "-out", temp_aln, "-clw", "-maxiters", "1", "-diags", "-quiet"]
            subprocess.run(cmd, check=True)
            
            # Cálculo de Consenso Frecuencial
            alineamiento = AlignIO.read(temp_aln, "clustal")
            consenso_seq = obtener_consenso_frecuencial(alineamiento, UMBRAL_CONSENSO)
            
            registro = SeqRecord(consenso_seq, id=grupo_name, description="Consenso_Frecuencial")
            lista_consensos.append(registro)
            
            print(f"✅ {grupo_name}: OK")
            
        except Exception as e:
            print(f"❌ Error en {grupo_name}: {e}")

    # Guardar Resultado Final
    if lista_consensos:
        SeqIO.write(lista_consensos, ARCHIVO_SALIDA, "fasta")
        print(f"\n🎉 ¡PROCESO TERMINADO!")
        print(f"Archivo generado: {ARCHIVO_SALIDA}")
        print(f"Total consensos: {len(lista_consensos)}")
        
        # Limpieza de archivos temporales
        if os.path.exists("temp_consenso.fasta"): os.remove("temp_consenso.fasta")
        if os.path.exists("temp_consenso.aln"): os.remove("temp_consenso.aln")
    else:
        print("\n⚠️ No se generaron consensos.")
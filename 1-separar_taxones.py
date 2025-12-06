import os
import re
from Bio import SeqIO

# --- CONFIGURACIÓN (LO ÚNICO QUE CAMBIAS) ---

# ¿Qué nivel taxonómico quieres separar?
# "p__" = Phylum
# "c__" = Clase (Útil para descomponer Arthropoda)
# "o__" = Orden
TAG_TAXONOMICO = "p__" 

# Nombres de tus archivos (deben estar en la misma carpeta que este script)
ARCHIVO_TAXONES = "taxones.fasta"
ARCHIVO_SECUENCIAS = "secuencias.fasta"

# Carpeta de salida (se crea sola)
# Ejemplo: si usas p__, la carpeta será "secuencias_por_p"
CARPETA_SALIDA = f"secuencias_por_{TAG_TAXONOMICO.replace('__','')}"

# --- FUNCIÓN DE LIMPIEZA PROFUNDA ---
def limpiar_nombre_taxon(nombre_sucio):
    """
    Recibe algo como: 'phylum_class_Arthropoda_123'
    Devuelve: 'Arthropoda'
    Esta función sirve para Phylum, Clase, Orden, etc.
    """
    nombre = nombre_sucio
    
    # 1. Lista de prefijos "basura" que suelen aparecer en bases de datos
    # (El orden importa: del más largo al más corto)
    prefijos_a_borrar = [
        "phylum_class_order_family_",
        "phylum_class_order_",
        "phylum_class_",
        "phylum_",
        "class_order_family_",
        "class_order_",
        "class_",
        "order_family_",
        "order_",
        "family_",
        "genus_"
    ]
    
    # Borramos prefijos recursivamente
    for prefijo in prefijos_a_borrar:
        if nombre.startswith(prefijo):
            nombre = nombre.replace(prefijo, "")
    
    # 2. Quitamos números finales (ej: _55528)
    nombre = re.sub(r'_\d+$', '', nombre)
    
    # 3. Quitamos caracteres prohibidos en nombres de archivo de Windows/Linux
    # (Esto evita errores al crear el archivo .fasta)
    nombre = re.sub(r'[\\/*?:"<>|]', "", nombre)
    
    # 4. Reemplazamos espacios por guiones bajos
    nombre = nombre.strip().replace(" ", "_")
    
    return nombre

def extraer_taxon(linea_taxonomia, tag_buscado):
    """
    Busca el tag (ej 'c__') y extrae el nombre limpio.
    """
    # Dividimos por punto y coma (jerarquía estándar)
    partes = linea_taxonomia.split(';')
    taxon_raw = "Unknown"
    
    for parte in partes:
        parte = parte.strip()
        if tag_buscado in parte:
            # Encontramos el nivel buscado
            inicio = parte.find(tag_buscado)
            taxon_raw = parte[inicio:]
            break
            
    if taxon_raw == "Unknown":
        return None

    # Quitamos el tag de búsqueda (ej: "p__")
    nombre_sin_tag = taxon_raw.replace(tag_buscado, "")
    
    # Aplicamos la limpieza profunda
    nombre_final = limpiar_nombre_taxon(nombre_sin_tag)
    
    # Filtros de calidad de la base de datos
    if not nombre_final: return None
    if nombre_final.lower() in ["unassigned", "possible_chimera", "uncultured"]:
        return None
        
    return nombre_final

# --- INICIO DEL PROGRAMA ---
if __name__ == "__main__":
    # Usamos rutas absolutas basadas en la ubicación del script (Mejor práctica para Git)
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    ruta_taxones = os.path.join(BASE_DIR, ARCHIVO_TAXONES)
    ruta_secuencias = os.path.join(BASE_DIR, ARCHIVO_SECUENCIAS)
    ruta_salida = os.path.join(BASE_DIR, CARPETA_SALIDA)

    if not os.path.exists(ruta_salida):
        os.makedirs(ruta_salida)

    print(f"--- SEPARADOR DE SECUENCIAS (Nivel: {TAG_TAXONOMICO}) ---")
    print(f"Leyendo desde: {BASE_DIR}")
    
    # 1. CARGAR MAPEO
    mapa_id_taxon = {}
    print("Cargando taxonomía...")
    
    try:
        with open(ruta_taxones, 'r', encoding='utf-8', errors='ignore') as f:
            for linea in f:
                linea = linea.strip()
                if not linea: continue
                
                partes = linea.split(None, 1)
                if len(partes) >= 2:
                    id_seq = partes[0].strip()
                    tax_str = partes[1]
                    
                    nombre_limpio = extraer_taxon(tax_str, TAG_TAXONOMICO)
                    if nombre_limpio:
                        mapa_id_taxon[id_seq] = nombre_limpio
                        
        print(f"✅ Mapeo completado: {len(mapa_id_taxon)} IDs encontrados para {TAG_TAXONOMICO}.")
        
    except FileNotFoundError:
        print(f"❌ Error: No encuentro '{ARCHIVO_TAXONES}'.")
        exit()

    # 2. PROCESAR SECUENCIAS
    print(f"Separando archivos en: {CARPETA_SALIDA}")
    
    if not os.path.exists(ruta_secuencias):
        print(f"❌ Error: No encuentro '{ARCHIVO_SECUENCIAS}'.")
        exit()

    count = 0
    archivos_generados = set()
    
    input_seq_iterator = SeqIO.parse(ruta_secuencias, "fasta")
    
    for record in input_seq_iterator:
        grupo = mapa_id_taxon.get(record.id)
        
        if grupo:
            archivo_destino = os.path.join(ruta_salida, f"{grupo}.fasta")
            
            with open(archivo_destino, "a") as f_out:
                SeqIO.write(record, f_out, "fasta")
            
            archivos_generados.add(grupo)
            count += 1
        else:
            # Opcional: Guardar los no identificados
            # with open(os.path.join(ruta_salida, "Sin_Clasificar.fasta"), "a") as f_out:
            #    SeqIO.write(record, f_out, "fasta")
            pass
            
        if count % 100000 == 0:
            print(f"   ... procesadas {count} secuencias", end='\r')

    print(f"\n\n¡LISTO! 🎉")
    print(f"Total secuencias organizadas: {count}")
    print(f"Total grupos (archivos) creados: {len(archivos_generados)}")
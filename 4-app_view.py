import streamlit as st
import pandas as pd
import os

# --- CONFIGURACIÓN DE PÁGINA ---
st.set_page_config(page_title="Bio-Barcodes Explorer", page_icon="🧬", layout="wide")

# --- RUTAS RELATIVAS (GIT FRIENDLY) ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR) # Subimos a 'problema_1'
RESULTADOS_DIR = os.path.join(BASE_DIR, "resultados")

st.title("🧬 Explorador de Barcodes Genéticos")
st.markdown("Herramienta para visualizar regiones diagnósticas generadas por el pipeline bioinformático.")

# --- BARRA LATERAL: CARGA DE DATOS ---
st.sidebar.header("📂 Cargar Datos")

# 1. Buscar archivos automáticamente
archivos_disponibles = []
if os.path.exists(RESULTADOS_DIR):
    archivos_disponibles = [f for f in os.listdir(RESULTADOS_DIR) if f.endswith(".xlsx") or f.endswith(".csv")]

archivo_seleccionado = None

# Opción A: Seleccionar de la carpeta
if archivos_disponibles:
    archivo_seleccionado = st.sidebar.selectbox(
        "Reportes encontrados en 'resultados/':", 
        archivos_disponibles,
        index=0
    )
    ruta_archivo = os.path.join(RESULTADOS_DIR, archivo_seleccionado)
else:
    st.sidebar.warning("⚠️ No se encontraron reportes en la carpeta 'resultados'.")

# Opción B: Subir manualmente
uploaded_file = st.sidebar.file_uploader("O sube tu archivo manualmente:", type=["xlsx", "csv"])

if uploaded_file:
    df = pd.read_excel(uploaded_file) if uploaded_file.name.endswith('.xlsx') else pd.read_csv(uploaded_file)
    st.sidebar.success("Archivo manual cargado.")
elif archivo_seleccionado:
    try:
        df = pd.read_excel(ruta_archivo) if ruta_archivo.endswith('.xlsx') else pd.read_csv(ruta_archivo)
        st.sidebar.success(f"Cargado: {archivo_seleccionado}")
    except Exception as e:
        st.error(f"Error leyendo el archivo: {e}")
        st.stop()
else:
    st.info("👈 Por favor, selecciona o sube un archivo de reporte para comenzar.")
    st.stop()

# --- FUNCIONES DE AYUDA ---
def calcular_calidad(secuencia):
    """Calcula % de bases puras (A,C,G,T) vs ruido (N, -)"""
    s = str(secuencia).upper()
    if len(s) == 0: return 0
    ruido = s.count('N') + s.count('-')
    return 100 * (1 - (ruido / len(s)))

def limpiar_blast(secuencia):
    """Quita guiones y N para NCBI BLAST"""
    return str(secuencia).upper().replace("-", "").replace("N", "").replace(" ", "")

# --- LÓGICA PRINCIPAL ---
if df is not None:
    # Asegurar columnas necesarias
    cols_req = ['Phylum', 'Diferencias_Minimas', 'Secuencia_Rango']
    if not all(c in df.columns for c in cols_req):
        st.error(f"El archivo debe tener las columnas: {cols_req}")
        st.stop()

    # Filtros laterales
    grupos = sorted(df['Phylum'].astype(str).unique())
    grupo_elegido = st.sidebar.selectbox("Filtrar por Grupo (Phylum/Clase):", grupos)
    
    # Filtrar datos del grupo
    datos = df[df['Phylum'] == grupo_elegido].copy()
    
    # Calcular Calidad
    datos['Calidad'] = datos['Secuencia_Rango'].apply(calcular_calidad)
    
    # Ordenar: Prioridad a Diferencias altas Y Calidad alta
    # (A veces es mejor 20 diferencias con 100% calidad que 25 diferencias llena de Ns)
    datos = datos.sort_values(by=['Diferencias_Minimas', 'Calidad'], ascending=[False, False])
    
    # --- MOSTRAR RESULTADOS ---
    st.divider()
    
    if not datos.empty:
        mejor = datos.iloc[0]
        
        # Encabezado del Ganador
        col_main, col_score = st.columns([3, 1])
        
        with col_main:
            st.subheader(f"🏆 Mejor Candidato para: *{grupo_elegido}*")
            
            # Secuencia Limpia
            seq_blast = limpiar_blast(mejor['Secuencia_Rango'])
            st.code(seq_blast, language="text")
            
            # Advertencias
            if len(seq_blast) < 15:
                st.warning("⚠️ La secuencia es muy corta después de limpiar las 'N'. No recomendada para BLAST.")
            if mejor['Calidad'] < 80:
                st.warning(f"⚠️ Calidad baja ({int(mejor['Calidad'])}%). Contiene muchas 'N' (ambigüedades).")
                st.caption("Esto suele ocurrir en grupos taxonómicos muy diversos donde el consenso no fue claro.")
            else:
                st.caption("✅ Copia esta secuencia y pégala en NCBI Nucleotide BLAST.")

        with col_score:
            st.metric("Diferencias", f"{mejor['Diferencias_Minimas']} bp", help="Letras distintas vs. el grupo más parecido")
            st.metric("Calidad", f"{int(mejor['Calidad'])}%", help="% de secuencia sin 'N' o huecos")
            st.metric("Rango", f"{mejor.get('Inicio', '?')} - {mejor.get('Fin', '?')}")

        # --- TABLA DETALLADA ---
        st.subheader("📋 Otros candidatos")
        
        # Formatear tabla
        tabla_show = datos[['Inicio', 'Fin', 'Diferencias_Minimas', 'Calidad', 'Secuencia_Rango']].head(50)
        st.dataframe(
            tabla_show.style.background_gradient(subset=['Diferencias_Minimas'], cmap='Greens'),
            use_container_width=True,
            hide_index=True
        )
        
    else:
        st.warning("No hay datos para este grupo.")

# --- FOOTER ---
st.markdown("---")
st.caption(f"Leyendo datos desde: `{BASE_DIR}`")
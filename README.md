Pipeline de Descubrimiento de Biomarcadores Genéticos (CO1) 🧬

Este proyecto implementa un flujo de trabajo bioinformático completo para identificar regiones diagnósticas ("Barcodes") en el gen mitocondrial Citocromo Oxidasa I (CO1) a partir de grandes volúmenes de datos taxonómicos.

🚀 Descripción del Proyecto

El objetivo es analizar millones de secuencias de ADN para encontrar regiones específicas que permitan diferenciar un grupo taxonómico (Phylum, Clase, Orden) de todos los demás.

El pipeline transforma datos crudos en una herramienta de diagnóstico visual, pasando por etapas de limpieza, alineamiento, consenso estadístico y minería de datos vectorizada.

🛠️ Tecnologías y Librerías

El proyecto está construido íntegramente en Python 3 utilizando las siguientes librerías científicas:

| Librería | Función Principal |
| Biopython (Bio) | Manejo de archivos FASTA, lectura de secuencias y parsing. |
| NumPy (numpy) | Cálculos matriciales de alto rendimiento y vectorización (aceleración x100). |
| Pandas (pandas) | Manipulación de tablas de datos y generación de reportes Excel. |
| Matplotlib (matplotlib) | Visualización de datos, perfiles de unicidad y mapas de calor. |
| Streamlit (streamlit) | Creación de la interfaz web interactiva para usuarios finales. |
| MUSCLE | Algoritmo externo de alineamiento múltiple de secuencias (MSA). |

⚙️ Estructura del Pipeline

Los scripts están numerados para ejecutarse secuencialmente. Cada uno prepara el terreno para el siguiente.

1️⃣ 01_separar_taxones.py (Data Wrangling)

Objetivo: Organizar el caos inicial.

Método: Lee los archivos gigantes taxones.fasta y secuencias.fasta. Separa las secuencias en archivos individuales según su nivel taxonómico (Phylum, Clase, etc.).

Limpieza: Filtra organismos no identificados, quimeras o entradas con taxonomía incompleta para evitar ruido.

2️⃣ 02_generar_consensos.py (Reducción Estadística)

Objetivo: Simplificar la complejidad biológica.

Método:

Toma una submuestra estadística aleatoria (n=150) de cada grupo.

Alinea internamente con MUSCLE.

Calcula una Secuencia Consenso basada en frecuencia.

Criterio de Calidad: Si una posición no tiene >60% de acuerdo, se marca como N (ambigüedad).

Finalmente, alinea todos los consensos entre sí para permitir la comparación.

3️⃣ 03_buscar_barcodes.py (Minería de Datos)

Objetivo: Encontrar la "aguja en el pajar".

Método:

Convierte el ADN a matrices numéricas (int8) usando NumPy.

Aplica una técnica de Sliding Window (Ventana Deslizante) de 100bp.

Calcula la Distancia de Hamming de "Todos contra Todos" simultáneamente (Broadcasting).

Identifica ventanas donde un grupo es significativamente distinto a su vecino más cercano.

4️⃣ 04_app_visualizador.py (Interfaz de Usuario)

Objetivo: Hacer los datos accesibles.

Características:

Panel interactivo web.

Filtra resultados por grupo taxonómico.

Muestra la secuencia limpia y lista para validación en NCBI BLAST.

Evalúa la calidad de la secuencia (porcentaje de nucleótidos vs ambigüedades).

📦 Instalación y Requisitos

Asegúrate de tener Python 3.8+ instalado.

Instala las dependencias necesarias:

pip install biopython numpy pandas matplotlib openpyxl streamlit



Descarga el ejecutable de MUSCLE (v3.8 recomendado) y colócalo en la carpeta raíz del proyecto.

▶️ Ejecución

Sigue el orden numérico de los scripts dentro de la carpeta scripts/:

# Paso 1: Organizar datos
python scripts/01_separar_taxones.py

# Paso 2: Generar consensos y alinear
python scripts/02_generar_consensos.py

# Paso 3: Buscar rangos diagnósticos
python scripts/03_buscar_barcodes.py

# Paso 4: Abrir la App Web
streamlit run scripts/04_app_visualizador.py



📊 Resultados Esperados

El pipeline generará en la carpeta resultados/:

Reportes Excel (.xlsx): Tablas detalladas con los mejores candidatos a biomarcadores, ordenados por diferenciación y calidad.

Gráficos de Unicidad (.png): Visualización de los picos de variabilidad a lo largo del gen.

Validación

Los barcodes generados han sido validados in-silico mediante matrices de distancia y externamente mediante NCBI Nucleotide BLAST, confirmando la especificidad taxonómica de las regiones seleccionadas.

Desarrollado para Investigathón Bioinformática 2025.

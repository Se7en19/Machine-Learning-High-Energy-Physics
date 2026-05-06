# Machine Learning High Energy Physics — Bar Strip Detector

Proyecto de simulación Geant4 y clasificación machine learning para separación de μ⁺/π⁺ con detector de barras centelladoras BC404.

## Estructura del proyecto

```
Machine-Learning-High-Energy-Physics/
├── Bar Strip Detector/
│   ├── simulation_mixed/          # Simulación Geant4 (haz mixto μ⁺/π⁺)
│   ├── img/                       # Plots generados
│   └── plot_all.py               # Script unificado de visualización
│
├── Classifier/                     # Clasificadores ML
├── Pion and muon simulation/       # Simulaciones individuales
├── Muon_dEdx_Simulation/           # Estudios dE/dx de muones
│
└── .scratch/                       # Issues y notas locales
```

## Pipeline end-to-end

### 1. Simulación Geant4

```bash
cd "Bar Strip Detector/simulation_mixed/build"
. /home/relka/Documentos/geant4-install/bin/geant4.sh
./sim ../barrido_continuo.mac
```

### 2. Generar plots

```bash
cd "Bar Strip Detector"
python plot_all.py --mixed "../Classifier/data/mixed/output_run*.root" --out img/
```

### 3. Clasificador

Abrir Jupyter Notebook:
```bash
cd "Classifier"
jupyter notebook muon_pion_classifier.ipynb
```

## Física del detector

- **Absorbedor**: G4_Fe 70×70×70 cm (4.17 λ_I)
- **Detector activo**: 2 capas de 20 barras BC404 (1 m × 5 cm × 1 cm)
  - Capa 1 (z = 100.5 cm): barras a lo largo de X, miden Y
  - Capa 2 (z = 103.5 cm): barras a lo largo de Y, miden X
- **Fuente**: punto fijo (0, 0, −2 m), mezcla 50/50 μ⁺/π⁺
- **Barrido**: 80 runs log, 50 MeV/c a 10 GeV/c, 2000 eventos/run

## Definición de un hit

Un evento se considera detectado cuando la partícula primaria (TrackID = 1) cumple esta secuencia:
1. Atraviesa el bloque de hierro sin ser absorbida
2. Entra en una barra de Capa 1 y deposita energía (dE/dx ≥ 0.05 MeV/mm)
3. Sale de Capa 1 y recorre los 2 cm de gap libre entre capas
4. Entra en una barra de Capa 2 y deposita energía (dE/dx ≥ 0.05 MeV/mm)
5. Sale del otro lado del centellador

El umbral de dE/dx = 0.05 MeV/mm elimina pasos espurios sin afectar la señal real (MIP ≈ 0.17 MeV/mm en BC404).

## Regímenes de eficiencia μ⁺

| Régimen | p₀ (MeV/c) | Comportamiento | ε |
|---------|-----------|----------------|---|
| I | < 500 | Muón no atraviesa Fe | ≈ 0% |
| II | 500–700 | Transición sigmoidal | Subida |
| III | > 700 | Meseta | ~89 ± 1% |

## Dependencias

- Geant4 11.x
- ROOT
- CMake ≥ 3.5
- Python: numpy, pandas, matplotlib, uproot, scikit-learn, xgboost, jupyter

## Contacto

Este proyecto forma parte del trabajo de investigación en Machine Learning aplicado a Física de Altas Energías.

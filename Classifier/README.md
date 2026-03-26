# Clasificador μ⁺ vs π⁺ — XGBoost

Clasifica muones y piones a partir de lo que registra el centellador BC404. El problema es que en el plateau relativista (p > 1 GeV/c) las dos partículas depositan energía de forma muy similar, así que el clasificador tiene que apoyarse en varias variables a la vez.

---

## Estructura

```
Classifier/
├── data/
│   ├── muon/    # output_run18.root ... output_run79.root  (62 archivos)
│   └── pion/    # output_run19.root ... output_run79.root  (61 archivos)
└── muon_pion_classifier.ipynb
```

Los archivos ROOT empiezan en run 18 (muones) y run 19 (piones) porque por debajo de ~170-183 MeV/c las partículas no atraviesan el absorbedor de 5 cm de hierro.

Cada archivo tiene un TTree llamado `Hits` con estos campos:

| Columna | Descripción | Unidades |
|---|---|---|
| `fEvent` | ID del evento | entero |
| `fX`, `fY`, `fZ` | Centro del detector (constante: 0, 0, 1100 mm) | mm |
| `fEdep` | Energía depositada en el paso | MeV |
| `fdEdx` | dE/dx = fEdep / longitud del paso | MeV/mm |
| `Ekin` | Energía cinética al inicio del paso | MeV |
| `TOF` | Tiempo de vuelo global | ns |
| `TrackLength` | Longitud total de traza acumulada | mm |
| `ScatteringAng` | Ángulo de dispersión en el paso | rad |
| `Momentum` | Módulo del momento | MeV/c |

Nota: `fX`, `fY`, `fZ` son el centro del volumen del detector, no la posición del paso. Son constantes para todos los hits. No sirven como features.

---

## Entorno

```bash
conda activate ML_HE_Physics
# numpy >= 2.0, pandas, uproot 5.x, xgboost, scikit-learn
```

---

## Pipeline del notebook

### 1. Carga de datos y curva de eficiencia

Se cargan todos los archivos ROOT por partícula. Para cada run se cuenta cuántos eventos únicos produjeron al menos un hit:

```python
epsilon(p) = N_eventos_detectados / 1000
```

Eso da la curva de eficiencia de detección vs momento. Para muones, epsilon sube de 0 a ~1 alrededor de 170 MeV/c y se mantiene alta. Para piones, el umbral es ~183 MeV/c y la eficiencia en el plateau es un poco menor (~85-90%) porque algunos piones interaccionan inelásticamente en el hierro.

### 2. Feature engineering

Cada evento puede tener varios pasos (steps) dentro del centellador. Se agregan en una sola fila por evento:

| Feature | Descripción |
|---|---|
| `n_steps` | Número de pasos registrados en el centellador |
| `edep_total`, `edep_mean`, `edep_max`, `edep_std` | Depósito de energía total y sus estadísticos |
| `dedx_mean`, `dedx_max`, `dedx_std` | dE/dx medio, máximo y desviación |
| `ekin_entry`, `ekin_exit`, `ekin_loss` | Energía cinética al entrar y salir del centellador |
| `tof_entry`, `tof_exit`, `tof_range` | Tiempo de vuelo al entrar, salir y rango temporal |
| `track_entry`, `track_exit`, `track_in_scint` | Longitud de traza acumulada al entrar y salir |
| `scat_mean`, `scat_max`, `scat_sum` | Ángulo de dispersión dentro del centellador |
| `momentum_entry`, `momentum_exit`, `momentum_loss` | Momento al entrar y salir |

Las columnas de posición (`fX`, `fY`, `fZ`) se excluyen porque son constantes.

### 3. Balance de clases

Se toman hasta 50 000 eventos por clase:

```python
N = min(50_000, len(events_mu), len(events_pi))
```

El dataset final tiene N muones y N piones mezclados y barajados.

### 4. Entrenamiento

```python
XGBClassifier(
    n_estimators=400,
    max_depth=5,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    eval_metric='logloss',
    early_stopping_rounds=20
)
```

Split 80/20 estratificado por clase. El modelo para en cuanto el log-loss en validación deja de bajar.

### 5. Evaluación

- Curva ROC y AUC
- Matriz de confusión
- Reporte de clasificación (precision, recall, F1)
- Curva de aprendizaje (log-loss vs round)
- Importancia de features por ganancia

### 6. Análisis por rango de momento

Se mapea cada evento a su momento de beam usando `run_id → MOMENTA_GeV`, y se calculan tasa de error y ROC-AUC por bin de momento:

| Rango | Por qué es interesante |
|---|---|
| 0.5-1 GeV/c | Zona de transición, eficiencia todavía sube |
| 1-2 GeV/c | Antes del plateau, mayor separación en dE/dx |
| 2-4 GeV/c | Plateau inicial |
| 4-7 GeV/c | Plateau estable, separación solo por TOF y scattering |
| 7-10 GeV/c | Alto momento, poca diferencia entre partículas |

---

## Por qué es difícil este problema

En el plateau (p > 1 GeV/c), a igual momento, μ⁺ y π⁺ depositan prácticamente la misma energía por unidad de longitud. La curva de Bethe-Bloch en función de bγ es universal para partículas cargadas pesadas, y a esos momenta ambas tienen bγ > 3. Las diferencias que quedan son:

- **TOF**: el muón es más ligero, llega antes al centellador a igual momento.
- **TrackLength**: recorre una trayectoria distinta antes de entrar al centellador.
- **ScatteringAngle**: el pión puede tener dispersión hadrónica residual.

Debajo de 500 MeV/c, la diferencia de masa (m_π/m_μ = 1.32) produce una separación visible en dE/dx porque las dos partículas están en puntos distintos de la curva de Bethe-Bloch. Ahí el clasificador lo tiene más fácil.

---

## Referencias

- Simulaciones: `../Pion and muon simulation/simulation_mu/` y `../simulation_pi/`
- Plots Bethe-Bloch: `../Pion and muon simulation/plot_bethe_bloch.py`
- Imagenes: `../Pion and muon simulation/img/`

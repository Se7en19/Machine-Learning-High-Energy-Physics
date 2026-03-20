# Clasificador muón vs pión — XGBoost

Clasifica muones (μ+) vs piones (π+) a partir de los hits registrados en las simulaciones Geant4.

---

## Estructura del proyecto

```
Classifier/
├── data/
│   ├── muon/          # 10 archivos output_run0.root ... output_run9.root
│   └── pion/          # 10 archivos output_run0.root ... output_run9.root
└── muon_pion_classifier.ipynb
```

Cada archivo ROOT tiene un TTree llamado `Hits` con estas columnas:

| Columna | Descripción |
|---|---|
| `fEvent` | ID del evento |
| `fX`, `fY`, `fZ` | Posición del hit (mm) |
| `fEdep` | Energía depositada (MeV) |
| `fdEdx` | Pérdida de energía por longitud (MeV/mm) |
| `Ekin` | Energía cinética (MeV) |
| `TOF` | Tiempo de vuelo (ns) |
| `TrackLength` | Longitud de trayectoria (mm) |
| `ScatteringAng` | Ángulo de scattering múltiple (rad) |
| `Momentum` | Momento (MeV/c) |

Los 10 runs por partícula cubren de 1.0 GeV a ~11 GeV en escala logarítmica, 1 000 eventos cada uno, para un total de 10 000 eventos por clase.

---

## Entorno

```bash
conda activate ML_HE_Physics
# numpy 2.2.6 + pandas 2.2.3 + uproot 5.x + xgboost + scikit-learn + shap + optuna
```

> numpy tiene que ser >= 2.0 para que pandas no rompa al crear índices string.
> Si aparece `Cannot convert numpy.ndarray`, ejecutar: `pip install "numpy>=2.0"`

---

## Plan del notebook (`muon_pion_classifier.ipynb`)

### Paso 1 — Carga de datos ✅

Cargar todos los archivos ROOT y concatenar en un DataFrame de hits:

```python
branches = ['fEvent', 'fX', 'fY', 'fZ', 'fEdep', 'fdEdx',
            'Ekin', 'TOF', 'TrackLength', 'ScatteringAng', 'Momentum']

def load_hits(pattern):
    files = sorted(glob(pattern))
    dfs = []
    for run_id, path in enumerate(files):
        with uproot.open(path) as f:
            tree = f['Hits']
            data = {b: tree[b].array(library='np').astype(np.float64) for b in branches}
            arr = np.column_stack([data[b] for b in branches])
            df = pd.DataFrame(arr, columns=branches)
            df['fEvent'] = df['fEvent'].astype(int)
            df['run_id'] = run_id
            df['event_uid'] = df['run_id'].astype(str) + '_' + df['fEvent'].astype(str)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

hits_mu = load_hits('data/muon/output_run*.root')
hits_pi = load_hits('data/pion/output_run*.root')
```

---

### Paso 2 — Feature engineering ✅

Agregar hits → 1 fila por evento usando `groupby('event_uid')`:

| Grupo | Features |
|---|---|
| Conteo | `n_hits`, `n_unique_cells` |
| Energía dep. | `edep_sum`, `edep_max`, `edep_std` |
| dE/dx | `dedx_mean`, `dedx_max`, `dedx_std` |
| Ekin | `ekin_first`, `ekin_last`, `ekin_loss` |
| TOF | `tof_first`, `tof_last`, `tof_range` |
| Track length | `track_first`, `track_last`, `track_mean` |
| Scattering | `scat_mean`, `scat_max`, `scat_std` |
| Geometría | `radial_spread`, `z_span` |

Resultado: DataFrame `df` con ~22 features + columna `label` (1=muón, 0=pión).

---

### Paso 3 — Train/test split ✅

```python
from sklearn.model_selection import train_test_split

FEATURES = [c for c in df.columns if c not in ('event_uid', 'label')]
X = df[FEATURES]
y = df['label']

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)
```

---

### Paso 4 — Entrenamiento XGBoost ✅

```python
from xgboost import XGBClassifier

model = XGBClassifier(
    n_estimators=300,
    max_depth=4,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    eval_metric='logloss',
    early_stopping_rounds=20,
    random_state=42
)

model.fit(
    X_train, y_train,
    eval_set=[(X_test, y_test)],
    verbose=50
)
```

---

### Paso 5 — Evaluación ✅

Métricas a reportar:
- ROC-AUC (métrica principal)
- Confusion matrix
- Classification report (precision, recall, F1)
- Feature importance con SHAP
- Curvas de aprendizaje

```python
from sklearn.metrics import roc_auc_score, confusion_matrix, classification_report
import shap

y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:, 1]

print('ROC-AUC:', roc_auc_score(y_test, y_prob))
print(classification_report(y_test, y_pred, target_names=['pion', 'muon']))

# SHAP
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values, X_test)
```

---

### Paso 6 — Análisis por energía ✅

Ver dónde falla el clasificador separando por bins de energía:

```python
# Añadir Ekin media del evento al DataFrame antes del split
df['ekin_mean_event'] = ...  # ekin promedio de los hits del evento

# Bins: [10, 50, 100, 300, 1000] MeV
# Para cada bin: calcular ROC-AUC y accuracy
# Esperado: peor separación en la región MIP (βγ ~ 3-4) donde μ y π son idénticos
```

---

### Paso 7 — Optimización de hiperparámetros ⬜ PENDIENTE

Usar Optuna para buscar los mejores hiperparámetros:

```python
import optuna

def objective(trial):
    params = {
        'n_estimators': trial.suggest_int('n_estimators', 100, 500),
        'max_depth': trial.suggest_int('max_depth', 3, 8),
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
        'subsample': trial.suggest_float('subsample', 0.6, 1.0),
        'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
    }
    model = XGBClassifier(**params, eval_metric='logloss', random_state=42)
    model.fit(X_train, y_train)
    return roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])

study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=50)
print('Best params:', study.best_params)
```

---

### Paso 8 — Clase `MuonPionClassifier` ⬜ PENDIENTE (entregable semana 13)

Todo el pipeline en una clase:

```python
class MuonPionClassifier:
    def __init__(self, **xgb_params): ...
    def load_data(self, muon_pattern, pion_pattern): ...
    def engineer_features(self, hits_df, label): ...
    def train(self, X_train, y_train): ...
    def evaluate(self, X_test, y_test): ...
    def predict(self, hits_df): ...
    def save(self, path): ...
    def load(self, path): ...
```

---

## Resultados

### Métricas globales

| Métrica | Valor |
|---|---|
| ROC-AUC | 1.000 |
| Eventos de prueba | 4 000 (2 000 muones, 2 000 piones) |
| Clasificaciones correctas | 3 999 / 4 000 |
| Tasa de error global | 0.025 % |

Un solo evento fue mal clasificado: un pión predicho como muón. Ningún muón fue clasificado como pión.

### Matriz de confusión

|  | Predicho: pión (0) | Predicho: muón (1) |
|---|---|---|
| Real: pión (0) | 1 999 | 1 |
| Real: muón (1) | 0 | 2 000 |

### Curva ROC y curva de aprendizaje

AUC = 1.00. El log-loss desciende desde ~0.65 en el primer round hasta prácticamente cero, con convergencia cerca del round 250. No hay señal de sobreajuste.

### Importancia de características (gain)

| Feature | Gain |
|---|---|
| `dedx_std` | 1 802.4 |
| `dedx_mean` | 851.1 |
| `dedx_max` | 631.9 |
| `tof_range` | 81.2 |
| `tof_first` | 8.6 |
| `n_unique_cells` | 5.7 |
| `n_hits` | 5.0 |
| `radial_spread` | 4.5 |
| `scat_max` | 3.0 |
| `edep_sum` | 2.2 |
| `track_mean` | 2.1 |
| `edep_std` | 1.6 |
| `edep_max` | 1.5 |
| `scat_std` | 1.4 |
| `track_first` | 1.1 |

Las tres variables de dE/dx concentran el 84 % de la ganancia total. El clasificador discrimina principalmente por la variabilidad del dE/dx paso a paso: el pión inicia cascadas hadrónicas que producen fluctuaciones grandes en los depósitos, el muón deja una traza de ionización comparativamente uniforme. `tof_range` entra como cuarta variable porque los secundarios de la cascada llegan al detector con tiempos retrasados respecto al primario.

### Error por rango de energía cinética

| Rango | N eventos | Error (%) | ROC-AUC |
|---|---|---|---|
| 1.0–1.5 GeV | 825 | 0.00 | 1.0 |
| 1.5–2.5 GeV | 753 | 0.00 | 1.0 |
| 2.5–4.5 GeV | 789 | 0.00 | 1.0 |
| 4.5–7.5 GeV | 824 | 0.12 | 1.0 |
| 7.5–11 GeV | 808 | 0.00 | 1.0 |

El único rango con errores es 4.5–7.5 GeV, con una tasa de 0.12 %. A esas energías la probabilidad de interacción hadrónica del pión es intermedia: algunos eventos no desarrollan una cascada completa y su perfil de depósito se acerca al del muón. El AUC de 1.0 en todos los rangos confirma que el modelo ordena correctamente los scores de probabilidad en cualquier región del espectro, aunque en ese intervalo cometa algunas asignaciones de clase en el umbral de decisión.

---

## Notas de física

- En la región MIP (mínimo de ionización, βγ ~ 3-4, Ekin ~ 300 MeV para piones), dE/dx es idéntico para μ y π. Aquí el clasificador tiene que apoyarse en `ScatteringAng` y la varianza de `fEdep`.
- Los piones sufren interacciones hadrónicas que producen fluctuaciones grandes en `fEdep` y `ScatteringAng`. Son las features más útiles a alta energía.
- A baja energía (< 100 MeV), la diferencia de masa (μ: 105.7 MeV, π: 139.6 MeV) hace la separación más sencilla por cinemática.

---

## Referencias

- Proyecto: https://www.fcfm.buap.mx/mrodriguez/ProyectoPPML.html (Fase 3, semanas 10-13)
- Simulaciones Geant4: `../Pion and muon simulation/`

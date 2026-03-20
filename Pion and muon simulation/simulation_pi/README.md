# Simulación Geant4 — Pión (π⁺)

Simulación de piones positivos en un detector de hierro segmentado. Produce los datasets ROOT que usa el clasificador μ⁺/π⁺.

---

## Descripción general

El cañón dispara piones (π⁺) con energía cinética variable en la dirección +z. Cada partícula atraviesa una cuadrícula de 100 celdas de hierro (10 × 10). A diferencia del muón, el pión puede iniciar cascadas hadrónicas al interactuar con los núcleos de hierro, generando decenas o cientos de secundarios por evento, cada uno con sus propios pasos registrados como hits.

El barrido cubre 11 puntos entre 1.0 y 2.6 GeV en escala logarítmica, 1 000 eventos por punto, lo que da 11 000 eventos en total.

---

## Geometría del mundo

| Parámetro | Valor |
|---|---|
| Material | G4_AIR |
| Semi-longitudes (x, y, z) | 0.6 m, 0.6 m, 8.0 m |
| Extensión en z | −8.0 m a +8.0 m |
| Posición del cañón | (0, 0, 0) |
| Dirección del haz | +z |

El mundo mide 8 m en semi-longitud z para contener sin recortes el detector, que llega hasta z = 7.05 m.

---

## Geometría del detector

| Parámetro | Valor |
|---|---|
| Material | G4_Fe (ρ = 7.874 g/cm³) |
| Celdas | 100 en total, cuadrícula 10 × 10 en x-y |
| Dimensiones de cada celda | 10 cm × 10 cm × 7 m |
| Centro en z | 3.55 m |
| Rango en z | 0.05 m a 7.05 m |
| Cobertura transversal | −0.5 m a +0.5 m en x e y |

Centro de cada celda en el plano transversal:

```
x_i = -0.5 m + (i + 0.5) × 0.1 m   (i = 0 … 9)
y_j = -0.5 m + (j + 0.5) × 0.1 m   (j = 0 … 9)
```

Índice de copia: `j + i×10`.

---

## Lista de física

```
G4EmStandardPhysics      — ionización, bremsstrahlung, procesos EM estándar
G4OpticalPhysics         — fotones ópticos
G4HadronPhysicsFTFP_BERT — hadrones (BERT < 3-6 GeV, FTFP por encima)
```

FTFP_BERT es la lista estándar para física hadrónica de altas energías en Geant4. Para el π⁺ es especialmente relevante porque el pión puede interactuar inelásticamente con los núcleos de hierro, produciendo secundarios (π⁰, protones de retroceso, neutrones) que a su vez dejan hits en el detector. El umbral de transición BERT → FTFP está entre 3 y 6 GeV para piones, aunque en este barrido las energías no superan los 2.6 GeV, por lo que el modelo BERT maneja todas las interacciones hadrónicas.

---

## Fuente primaria

| Parámetro | Valor |
|---|---|
| Partícula | π⁺ (masa = 139.57 MeV/c²) |
| Posición | (0, 0, 0) |
| Dirección | (0, 0, 1) |
| Energía | variable, definida por barrido_continuo.mac |
| Partículas por evento | 1 |

---

## Barrido en energía

El archivo `barrido_continuo.mac` contiene 11 runs:

| Parámetro | Valor |
|---|---|
| Energía mínima | 1.0000 GeV |
| Energía máxima | 2.6102 GeV |
| Número de runs | 11 |
| Espaciado | logarítmico uniforme |
| Eventos por run | 1 000 |
| Total de eventos | 11 000 |

Los archivos de salida son `output_run0.root` a `output_run10.root`.

---

## Definición de un hit

Un hit se registra cada vez que un paso de una partícula ocurre dentro de una celda sensible. El detector implementa `G4VSensitiveDetector::ProcessHits`, que Geant4 llama automáticamente para cada paso cuyo volumen de pre-step sea el volumen sensible.

No hay umbral de energía mínima. En eventos de pión con cascada hadrónica, una sola celda puede acumular cientos de hits de secundarios: `n_unique_cells` es bajo pero `n_hits` puede ser muy alto.

Columnas del árbol `Hits`:

| Columna ROOT | Descripción | Fuente en Geant4 |
|---|---|---|
| `fEvent` | ID del evento dentro del run | `G4Event::GetEventID()` |
| `fX`, `fY`, `fZ` | Centro geométrico de la celda (mm) | `G4VPhysicalVolume::GetTranslation()` |
| `fEdep` | Energía depositada en el paso (MeV) | `G4Step::GetTotalEnergyDeposit()` |
| `fdEdx` | Energía por unidad de longitud (MeV/mm) | `fEdep / G4Step::GetStepLength()` |
| `Ekin` | Energía cinética al inicio del paso (MeV) | `G4StepPoint::GetKineticEnergy()` |
| `TOF` | Tiempo global al inicio del paso (ns) | `G4StepPoint::GetGlobalTime()` |
| `TrackLength` | Longitud total de traza hasta ese paso (mm) | `G4Track::GetTrackLength()` |
| `ScatteringAng` | Ángulo entre dirección pre-step y post-step (rad) | `dirPre.angle(dirPost)` |
| `Momentum` | Módulo del momento al inicio del paso (MeV/c) | `G4StepPoint::GetMomentum().mag()` |

`fX`, `fY`, `fZ` son el centro geométrico de la celda, no la posición exacta del paso. Todas las filas de una misma celda tienen los mismos valores de estas tres columnas.

---

## Salida

Cada run genera un archivo `output_runN.root` con un `TTree` llamado `Hits`. El archivo se abre en `BeginOfRunAction` y se cierra en `EndOfRunAction` usando `G4AnalysisManager`.

---

## Compilación y ejecución

```bash
cd build
cmake ..
make -j4

./sim ../barrido_continuo.mac
```

---

## Física del pión en hierro

El π⁺ es un hadrón (quark u, antiquark d̄). Cuando entra en el hierro, puede interactuar mediante la fuerza fuerte con los núcleos de Fe, iniciando una cascada hadrónica. Esto produce tres efectos observables que lo distinguen del muón:

La cascada convierte buena parte de la energía cinética en secundarios que se frenan dentro del detector. El depósito de energía total por evento (`edep_sum`) es mucho mayor que en el muón a igual energía de entrada.

Los secundarios se dispersan lateralmente fuera del eje del haz. El `radial_spread` y el `ScatteringAng` promedio son mayores que los del muón, que viaja en línea casi recta.

Las interacciones hadrónicas tienen fluctuaciones intrínsecamente grandes. La distribución de `fEdep` paso a paso es más ancha que la de ionización pura, con colas más pesadas.

En la región MIP (βγ ≈ 3–4, unos 400–600 MeV para el π⁺), el dE/dx del pión por ionización es prácticamente igual al del muón. A 1–2.6 GeV el pión está por encima del MIP y las cascadas hadrónicas son el mecanismo dominante de separación entre las dos partículas.

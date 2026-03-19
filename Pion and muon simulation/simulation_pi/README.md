# Simulación Geant4 — Pión (π⁺)

Simulación de piones positivos atravesando un detector de hierro segmentado. Genera los datasets ROOT que alimentan el clasificador μ⁺/π⁺.

---

## Descripción general

Un cañón de partículas dispara piones (π⁺) con energía cinética variable a lo largo del eje z. Las partículas atraviesan un arreglo de 10 000 celdas de hierro (100 × 100) y cada interacción dentro de una celda sensible queda registrada como un *hit* en un árbol ROOT.

El barrido en energía cubre de 1.0 GeV a 10.0 GeV en 60 puntos distribuidos en escala logarítmica, con 1 000 eventos por punto. Cada run produce un archivo ROOT independiente.

---

## Geometría del mundo

| Parámetro | Valor |
|---|---|
| Volumen world | caja de aire (G4_AIR), 1 m × 1 m × 10 m |
| Material world | G4_AIR |
| Origen del cañón | (0, 0, 0) mm |
| Dirección del haz | +z |

El mundo tiene semi-longitudes (0.5, 0.5, 5.0) m en (x, y, z).

---

## Geometría del detector

El detector es un arreglo de 100 × 100 celdas de hierro dispuesto perpendicularmente al haz.

| Parámetro | Valor |
|---|---|
| Material | G4_Fe (hierro, ρ = 7.874 g/cm³) |
| Dimensiones de cada celda | 1 cm × 1 cm × 7 m |
| Número de celdas | 10 000 (100 × 100 en x-y) |
| Centro del arreglo en z | 3.55 m |
| Rango en z | 0.05 m → 7.05 m |
| Cobertura transversal | −0.5 m → +0.5 m en x e y |

Cada celda ocupa una posición única en la cuadrícula. Su centro en el plano transversal es:

```
x_i = -0.5 m + (i + 0.5) × 0.01 m   (i = 0 … 99)
y_j = -0.5 m + (j + 0.5) × 0.01 m   (j = 0 … 99)
```

El índice de copia de cada celda es `j + i×100`.

---

## Lista de física

```
G4EmStandardPhysics      — procesos electromagnéticos estándar
G4OpticalPhysics         — fotones ópticos
G4HadronPhysicsFTFP_BERT — hadrones: FTFP por encima de 3–6 GeV, BERT por debajo
```

FTFP_BERT es la lista recomendada para física hadrónica de altas energías. Para el π⁺ esto es especialmente relevante: el pión puede iniciar cascadas hadrónicas a través de interacciones nucleares inelásticas (absorción, producción de secundarios), lo que genera depósitos de energía muy distintos a los del muón. El umbral de transición BERT→FTFP está entre 3 y 6 GeV para piones.

---

## Fuente primaria (cañón)

| Parámetro | Valor |
|---|---|
| Partícula | π⁺ (pión positivo, masa = 139.57 MeV/c²) |
| Posición | (0, 0, 0) |
| Dirección | (0, 0, 1) — paralela al eje z |
| Energía | variable (ver barrido) |
| Número de partículas por evento | 1 |

El valor por defecto en el código es 100 GeV de momento; en el barrido de producción la energía la sobreescribe `barrido_continuo.mac`.

---

## Barrido en energía

El archivo `barrido_continuo.mac` define 60 runs en escala logarítmica:

| Parámetro | Valor |
|---|---|
| Energía mínima | 1.0000 GeV |
| Energía máxima | 10.000 GeV |
| Número de puntos | 60 |
| Espaciado | logarítmico uniforme |
| Eventos por run | 1 000 |
| Total de eventos generados | 60 000 |

Cada run produce un archivo `output_runN.root` (N = 0 … 59).

---

## Definición de un hit

Un hit se registra cuando un paso (*step*) de una partícula ocurre dentro de una celda sensible. El detector implementa `G4VSensitiveDetector::ProcessHits`, que se invoca automáticamente por Geant4 para cada paso cuyo volumen de pre-step es el volumen sensible.

No existe umbral de energía mínima: cualquier paso dentro de una celda, independientemente de `fEdep`, genera una fila en el árbol. En la práctica, los piones producen cascadas hadrónicas que generan muchos secundarios (π⁰, protones de retroceso, neutrones), cada uno de los cuales contribuye con sus propios hits si penetra una celda sensible.

Por cada hit se registran las siguientes cantidades:

| Columna ROOT | Descripción | Fuente en Geant4 |
|---|---|---|
| `fEvent` | ID del evento dentro del run | `G4Event::GetEventID()` |
| `fX`, `fY`, `fZ` | Posición del centro de la celda detectora (mm) | `G4VPhysicalVolume::GetTranslation()` |
| `fEdep` | Energía total depositada en el paso (MeV) | `G4Step::GetTotalEnergyDeposit()` |
| `fdEdx` | Energía depositada dividida por la longitud del paso (MeV/mm) | `fEdep / G4Step::GetStepLength()` |
| `Ekin` | Energía cinética de la partícula al inicio del paso (MeV) | `G4StepPoint::GetKineticEnergy()` (pre-step) |
| `TOF` | Tiempo global al inicio del paso (ns) | `G4StepPoint::GetGlobalTime()` (pre-step) |
| `TrackLength` | Longitud total recorrida por la traza hasta ese paso (mm) | `G4Track::GetTrackLength()` |
| `ScatteringAng` | Ángulo entre la dirección pre-step y post-step (rad) | `dirPre.angle(dirPost)` |
| `Momentum` | Módulo del momento al inicio del paso (MeV/c) | `G4StepPoint::GetMomentum().mag()` |

> **Nota sobre fX, fY, fZ**: estas columnas no corresponden a la posición exacta del paso sino al centro geométrico de la celda que lo contiene. Todas las filas de una misma celda tienen los mismos valores de fX, fY, fZ. En el caso del pión, una sola celda puede concentrar cientos de hits de secundarios, lo que resulta en `n_unique_cells` bajo pero `n_hits` muy alto.

---

## Salida

Cada run genera un archivo `output_runN.root` en el directorio de ejecución con un único `TTree` llamado `Hits`. El archivo se abre, escribe y cierra en `BeginOfRunAction` / `EndOfRunAction` usando `G4AnalysisManager`.

---

## Compilación y ejecución

```bash
cd build
cmake ..
make -j4

# Barrido completo en energía (60 runs × 1000 eventos)
./sim ../barrido_continuo.mac
```

---

## Notas de física relevantes

El pión cargado es un hadrón (quark u y antiquark d̄). A energías de GeV interactúa con los núcleos de hierro mediante la fuerza fuerte, produciendo cascadas hadrónicas con características muy distintas a las del muón:

- **Mayor deposición de energía**: la interacción hadrónica convierte gran parte de la energía cinética en secundarios que se detienen dentro del detector.
- **Dispersión lateral mayor**: los secundarios de la cascada se alejan del eje del haz, aumentando `radial_spread` y `ScatteringAng`.
- **Fluctuaciones en fEdep**: las interacciones hadrónicas tienen varianza intrínsecamente mayor que la ionización electromagnética, generando distribuciones de Landau más anchas.

En la región MIP (βγ ≈ 3–4, aproximadamente 400–600 MeV para el π⁺), el dE/dx del pión por ionización es prácticamente idéntico al del muón. A energías de GeV, el pión ya está bien por encima del MIP y la diferencia con el muón se vuelve más pronunciada a través de las cascadas hadrónicas, que el μ⁺ no puede iniciar.

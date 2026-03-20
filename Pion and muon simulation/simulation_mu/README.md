# Simulación Geant4 — Muón (μ⁺)

Simulación de muones positivos en un detector de hierro segmentado. Produce los datasets ROOT que usa el clasificador μ⁺/π⁺.

---

## Descripción general

El cañón dispara muones (μ⁺) con energía cinética variable en la dirección +z. Cada partícula atraviesa una cuadrícula de 100 celdas de hierro (10 × 10), y cada paso dentro de una celda sensible genera una fila en el árbol ROOT.

El barrido cubre 10 puntos entre 1.0 y 10.0 GeV en escala logarítmica, 1 000 eventos por punto, lo que da 10 000 eventos en total.

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

Para el μ⁺, los procesos que importan son ionización y bremsstrahlung. Las interacciones hadrónicas no aplican al muón directamente, aunque sí pueden afectar a los secundarios que produce en el hierro.

---

## Fuente primaria

| Parámetro | Valor |
|---|---|
| Partícula | μ⁺ (masa = 105.66 MeV/c²) |
| Posición | (0, 0, 0) |
| Dirección | (0, 0, 1) |
| Energía | variable, definida por barrido_continuo.mac |
| Partículas por evento | 1 |

---

## Barrido en energía

El archivo `barrido_continuo.mac` contiene 10 runs:

| Parámetro | Valor |
|---|---|
| Energía mínima | 1.0000 GeV |
| Energía máxima | 10.000 GeV |
| Número de runs | 10 |
| Espaciado | logarítmico uniforme |
| Eventos por run | 1 000 |
| Total de eventos | 10 000 |

Los archivos de salida son `output_run0.root` a `output_run9.root`.

---

## Definición de un hit

Un hit se registra cada vez que un paso de una partícula ocurre dentro de una celda sensible. El detector implementa `G4VSensitiveDetector::ProcessHits`, que Geant4 llama automáticamente para cada paso cuyo volumen de pre-step sea el volumen sensible.

No hay umbral de energía mínima: cualquier paso dentro de una celda genera una fila en el árbol, independientemente de `fEdep`.

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

## Física del muón en hierro

A energías de GeV el muón está en la región relativista de Bethe-Bloch. Pierde energía principalmente por ionización, con una tasa cercana al mínimo (MIP). En 7 m de hierro, los muones de menor energía del barrido (cercanos a 1 GeV) pueden detenerse dentro del detector, dejando una traza completa. A partir de unos pocos GeV el muón atraviesa el volumen sin frenarse del todo.

El muón no inicia cascadas hadrónicas. Eso es lo que lo distingue del pión en este experimento: su traza es limpia, estrecha en el plano transversal, y con depósitos de energía relativamente uniformes a lo largo del recorrido.

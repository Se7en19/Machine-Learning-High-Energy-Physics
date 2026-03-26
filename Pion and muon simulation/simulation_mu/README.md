# Simulación Geant4 — Muón (μ⁺)

Dispara muones positivos contra un absorbedor de hierro delgado y registra lo que llega al centellador. Los archivos ROOT que genera son la mitad del dataset que usa el clasificador.

---

## Qué hace esta simulación

Un cañón de partículas lanza μ⁺ desde z = -2 m en dirección +z. Las partículas atraviesan 5 cm de hierro, cruzan 1 m de vacío, y si llegan al centellador BC404, se registra el evento. Solo se guarda la traza primaria (TrackID = 1), no los secundarios.

El barrido cubre 80 momenta distintos entre 50 MeV/c y 10 GeV/c en escala logarítmica. A menos de ~150 MeV/c el muón se detiene en el hierro y no produce hits, así que los primeros archivos con datos reales aparecen alrededor del run 18.

---

## Geometría

```
[μ⁺ gun, z=-2m]  →→→  [Fe 5cm, z=0-5cm]  →→→  [1m vacío]  →→→  [BC404, z=105-115cm]
```

### Mundo

| Parámetro | Valor |
|---|---|
| Material | G4_Galactic (vacío) |
| Semi-longitudes | 6 m × 6 m × 5 m |

### Absorbedor de hierro

| Parámetro | Valor |
|---|---|
| Material | G4_Fe (rho = 7.874 g/cm³) |
| Dimensiones | 5 cm × 5 cm × 5 cm |
| Centro en z | 2.5 cm |

5 cm es 0.30 longitudes de interacción nuclear (lambda_I ~ 16.77 cm en Fe). Suficiente para absorber piones de bajo momento pero no para detener muones relativistas.

### Centellador (detector activo)

| Parámetro | Valor |
|---|---|
| Material | G4_PLASTIC_SC_VINYLTOLUENE (BC404, rho = 1.032 g/cm³) |
| Dimensiones | 10 m × 10 m × 10 cm |
| Centro en z | 110 cm |
| Cara delantera | z = 105 cm |
| Cara trasera | z = 115 cm |

El centellador es intencionalmente grande en x-y para capturar cualquier dispersión transversal.

---

## Lista de física

```
FTFP_BERT
```

Para el muón los procesos relevantes son ionización y bremsstrahlung. La física hadrónica no le aplica directamente, pero está activa por coherencia con la simulación de piones.

---

## Fuente primaria

| Parámetro | Valor |
|---|---|
| Partícula | μ⁺ (masa = 105.66 MeV/c²) |
| Posición | (0, 0, -2 m) |
| Dirección | +z |
| Momento | variable, definido por barrido_continuo.mac |
| Eventos por run | 1 000 |

El momento se fija con `/gun/momentumAmp X GeV`. Geant4 calcula internamente la energía cinética como E_kin = sqrt(p² + m²) - m.

---

## Barrido en momento

80 runs logarítmicamente espaciados:

| Parámetro | Valor |
|---|---|
| Momento mínimo | 50 MeV/c (run 0) |
| Momento máximo | 10 000 MeV/c (run 79) |
| Runs totales | 80 |
| Eventos por run | 1 000 |
| Eventos generados | 80 000 |

Los runs 0 a 17 (p < ~170 MeV/c) no producen hits porque el muón no alcanza el centellador. Los primeros archivos con datos son output_run18.root en adelante.

El rango cubre desde la región de alto dE/dx (parte descendente de Bethe-Bloch, bγ < 2) hasta el plateau de Fermi (p > 500 MeV/c).

---

## Definición de un hit

Un hit es un paso de la partícula primaria dentro del volumen del centellador BC404. Se filtra explícitamente por `TrackID == 1` para ignorar secundarios.

Columnas del árbol `Hits`:

| Columna | Descripción | Unidades |
|---|---|---|
| `fEvent` | ID del evento dentro del run | entero |
| `fX`, `fY`, `fZ` | Centro del volumen detector (constante: 0, 0, 1100 mm) | mm |
| `fEdep` | Energía depositada en el paso | MeV |
| `fdEdx` | fEdep / longitud del paso | MeV/mm |
| `Ekin` | Energía cinética al inicio del paso | MeV |
| `TOF` | Tiempo global al inicio del paso | ns |
| `TrackLength` | Longitud total acumulada de la traza | mm |
| `ScatteringAng` | Ángulo entre dirección pre-step y post-step | rad |
| `Momentum` | Módulo del momento al inicio del paso | MeV/c |

Nota: `fX`, `fY`, `fZ` son el centro del volumen físico del detector, no la posición exacta del paso. Todas las filas tienen los mismos valores (0, 0, 1100).

---

## Salida

Los archivos se guardan directamente en `../../../Classifier/data/muon/output_runN.root`. Cada uno contiene un TTree llamado `Hits`.

---

## Compilación y ejecución

```bash
cd build
cmake ..
make -j4
./sim ../barrido_continuo.mac
```

---

## Curva de eficiencia

La fracción de eventos que produce al menos un hit en el centellador es:

```
epsilon(p) = N_detectados / 1000
```

Para muones, epsilon sube de 0 en los runs de bajo momento (partícula frenada en el hierro) a ~1 en cuanto el momento supera ~170 MeV/c. A partir de ahí se mantiene cerca de 1 en todo el plateau.

---

## Física del muón

El muón es un leptón, no interacciona mediante la fuerza fuerte. Su pérdida de energía en el centellador sigue la distribución de Landau con MPV en torno a 0.17 MeV/mm para partículas relativistas en BC404. No genera cascadas hadrónicas. Eso hace que su traza sea limpia y su dE/dx predecible, lo que es tanto una ventaja para el clasificador como el motivo por el que resulta difícil de distinguir del pión en el plateau relativista.

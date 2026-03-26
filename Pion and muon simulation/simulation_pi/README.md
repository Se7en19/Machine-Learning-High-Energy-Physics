# Simulación Geant4 — Pión (π⁺)

Igual que la simulación de muones, pero con piones. La diferencia importante es que el pión puede iniciar cascadas hadrónicas al pasar por el hierro, y eso complica todo.

---

## Qué hace esta simulación

Un cañón dispara π⁺ desde z = -2 m. Las partículas atraviesan 5 cm de hierro, cruzan 1 m de vacío, y si el pión llega al centellador BC404, se registra el evento. Solo se guarda la traza primaria (TrackID = 1).

A momenta bajos, los piones se frenan en el hierro con más facilidad que los muones porque son más pesados y también pueden interactuar inelásticamente con los núcleos de Fe. Los primeros archivos con datos aparecen alrededor del run 19 (p ~ 183 MeV/c).

---

## Geometría

```
[π⁺ gun, z=-2m]  →→→  [Fe 5cm, z=0-5cm]  →→→  [1m vacío]  →→→  [BC404, z=105-115cm]
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

Con 5 cm de Fe se absorben la mayoría de los piones por debajo de ~300 MeV/c, ya sea por frenado o por interacción hadrónica. Eso deja la región de bajo momento sin datos, que es justo la región donde muones y piones se separan mejor cinématicamente.

### Centellador (detector activo)

| Parámetro | Valor |
|---|---|
| Material | G4_PLASTIC_SC_VINYLTOLUENE (BC404, rho = 1.032 g/cm³) |
| Dimensiones | 10 m × 10 m × 10 cm |
| Centro en z | 110 cm |
| Cara delantera | z = 105 cm |
| Cara trasera | z = 115 cm |

---

## Lista de física

```
FTFP_BERT
```

FTFP_BERT es la lista estándar para física hadrónica en Geant4. Para el pión es importante porque puede producir interacciones inelásticas con los núcleos de hierro, aunque en esta simulación esos secundarios no se registran (se filtra TrackID = 1). El efecto visible es que algunos piones no llegan al centellador porque se absorben en el hierro.

---

## Fuente primaria

| Parámetro | Valor |
|---|---|
| Partícula | π⁺ (masa = 139.57 MeV/c²) |
| Posición | (0, 0, -2 m) |
| Dirección | +z |
| Momento | variable, definido por barrido_continuo.mac |
| Eventos por run | 1 000 |

---

## Barrido en momento

80 runs logarítmicamente espaciados, idéntico al barrido de muones:

| Parámetro | Valor |
|---|---|
| Momento mínimo | 50 MeV/c (run 0) |
| Momento máximo | 10 000 MeV/c (run 79) |
| Runs totales | 80 |
| Eventos por run | 1 000 |

Los runs 0 a 18 (p < ~183 MeV/c) no producen hits. El primer archivo con datos es output_run19.root.

El umbral es ligeramente más alto que para muones (~170 MeV/c) porque el pión es más pesado y tiene mayor probabilidad de interacción hadrónica en el hierro.

---

## Definición de un hit

Un paso de la partícula primaria (TrackID = 1) dentro del volumen BC404.

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

---

## Salida

Los archivos se guardan en `../../../Classifier/data/pion/output_runN.root`. Mismo formato que la simulación de muones.

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

```
epsilon(p) = N_detectados / 1000
```

La eficiencia sube de 0 a ~1 entre 150 y 300 MeV/c. En la región del plateau (p > 500 MeV/c) se estabiliza pero no llega al 100% porque algunos piones siguen interaccionando inelásticamente en el hierro y no alcanzan el centellador. Eso explica por qué hay menos hits por run en piones (~2400) que en muones (~3400).

---

## Física del pión

El pión es un hadrón (par quark-antiquark). Tiene dos modos de interacción en el hierro que no tiene el muón:

**Interacción fuerte con núcleos de Fe.** Produce cascadas hadrónicas que pueden absorber al pión antes de llegar al centellador. A bajo momento esto domina y es la principal fuente de ineficiencia.

**Mayor masa.** A igual momento, el pión tiene menor bγ que el muón (m_π/m_μ = 1.32). Eso significa que a la misma p en GeV/c, el pión deposita más energía por unidad de longitud. Esta diferencia de dE/dx vs p es el discriminador principal que usa el clasificador.

En el plateau de Fermi (p > 1 GeV/c), las dos curvas convergen y la separación se vuelve difícil usando solo dE/dx. Ahí entran el TOF, la longitud de traza y el ángulo de scattering como variables complementarias.

# Muon dE/dx Simulation (Geant4 B2a Modificado)

## 1. Objetivo científico

Este repositorio implementa una simulación de pérdida de energía por ionización (`dE/dx`) para muones (`mu+`, `mu-`) en cámaras gaseosas de Xenón, basada en el ejemplo `B2a` de Geant4 y modificada para:

- generar un espectro de momento de muones en un rango relativista;
- registrar depósitos de energía paso a paso en volúmenes sensibles;
- reconstruir `dE/dx` por cámara y por evento primario;
- producir un `TTree` ROOT para análisis comparativo con curvas tipo Bethe-Bloch.

El caso de uso principal es validación física de la dependencia de `dE/dx` respecto a `p/q`, `beta` y `beta*gamma`.

## 2. Modelo físico implementado

### 2.1 Procesos de física

La lista de física usada es:

- `FTFP_BERT`
- `G4StepLimiterPhysics` (para respetar límite máximo de paso en la región tracker)

La ionización de partículas cargadas en materia proviene de los procesos EM de Geant4 incluidos en `FTFP_BERT`.

### 2.2 Variable de interés: dE/dx

En cada cámara y para cada evento, se acumula:

- `Edep_total = sum_i Edep_i`
- `L_total = sum_i stepLength_i`

y se define:

- `dE/dx = Edep_total / L_total`

Con las unidades internas de Geant4 usadas en este código:

- `Edep` en MeV
- `stepLength` en mm
- `dE/dx` en MeV/mm

### 2.3 Rigidez

Se guarda la rigidez como:

- `rigidity = |p| / q`

con signo determinado por la carga (`mu+` positivo, `mu-` negativo). En la implementación actual se toma el momento del primer hit primario del evento.

## 3. Geometría y materiales

La geometría se define en [`src/DetectorConstruction.cc`](src/DetectorConstruction.cc).

### 3.1 Materiales

- Mundo y volumen tracker: `G4_AIR`
- Blanco (target): `G4_Pb`
- Cámaras sensibles: `G4_Xe`

### 3.2 Dimensiones (configuración por defecto)

- Número de cámaras: `5`
- Espaciado entre centros de cámara: `80 cm`
- Espesor de cámara (eje z): `20 cm` (semi-espesor `10 cm`)
- Tracker cilíndrico: semi-longitud `240 cm`, radio máximo `240 cm`
- Blanco de Pb: cilindro de radio `2.5 cm`, longitud `5 cm`
- Mundo: cubo de lado `588 cm`

Posiciones de centros de cámara (z): `-160, -80, 0, 80, 160 cm`.

Radii de cámaras (crecientes por diseño): `24, 78, 132, 186, 240 cm`.

## 4. Generación primaria de eventos

Definida en [`src/PrimaryGeneratorAction.cc`](src/PrimaryGeneratorAction.cc).

Por evento:

- partícula: `mu+` o `mu-` con probabilidad 50/50;
- dirección: fija a `+z`;
- vértice: cerca del borde de entrada del mundo (`z = -Z_world + 1 um`);
- momento `p`: muestreo log-uniforme en `[200, 3000] MeV/c`;
- energía cinética: calculada a partir de `E = sqrt(p^2 + m_mu^2)` con `m_mu = 105.6583755 MeV/c^2`.

## 5. Detección sensible y reconstrucción por evento

### 5.1 Hits (paso a paso)

En [`src/TrackerSD.cc`](src/TrackerSD.cc), para cada paso en cámaras sensibles:

- se ignoran pasos con `Edep == 0` o `stepLength == 0`;
- se registra `trackID`, `chamberNb`, `Edep`, posición, longitud de paso, momento pre-step y carga.

### 5.2 Reducción a observables por cámara

En [`src/EventAction.cc`](src/EventAction.cc):

- solo se usa la traza primaria (`trackID == 1`);
- se acumula `Edep` y longitud por cámara;
- se calcula `dE/dx` por cámara;
- se escribe una fila de ntuple por cada cámara con longitud acumulada positiva.

## 6. Salida de datos ROOT

Configurada en [`src/RunAction.cc`](src/RunAction.cc).

Archivo de salida:

- `B2a_output.root`

`TTree`: `Hits`

Columnas:

1. `eventID` (`int`)
2. `trackID` (`int`, en esta reducción siempre `1`)
3. `chamberNb` (`int`, `0..4`)
4. `dedx` (`double`, MeV/mm)
5. `rigidity` (`double`, MeV/c con signo de carga)
6. `ekin` (`double`)  
   Nota: actualmente esta columna se rellena con el mismo valor que `rigidity` (no con energía cinética real).
7. `particle` (`string`)  
   Nota: en el estado actual no se asigna explícitamente en `TrackerSD`, por lo que puede quedar vacío.

## 7. Postprocesado y validación física

El macro ROOT [`analysis/plot_dedx.C`](analysis/plot_dedx.C) genera:

- `dE/dx vs p/q`
- `dE/dx vs beta`
- `dE/dx vs beta*gamma`
- hits por cámara
- distribución tipo Landau por cámara
- mapa espacial/fallback por cámara

Incluye una función analítica de Bethe-Bloch con parámetros efectivos de Xenón gaseoso para superposición cualitativa con la simulación.

## 8. Flujo de ejecución recomendado

### 8.1 Requisitos

- Geant4 11.x con módulo de análisis ROOT habilitado (`g4analysis`)
- ROOT (para ejecutar `analysis/plot_dedx.C`)
- CMake >= 3.16
- compilador C++17 compatible

### 8.2 Compilar

```bash
cmake -S . -B build
cmake --build build -j
```

### 8.3 Ejecutar simulación en batch

```bash
cd build
./exampleB2a ../run_dedx.mac
```

### 8.4 Ejecutar análisis ROOT

```bash
cd build
root -l -q '../analysis/plot_dedx.C'
```

## 9. Comandos de UI disponibles

Vía `DetectorMessenger`:

- `/B2/det/setTargetMaterial <G4_MATERIAL>`
- `/B2/det/setChamberMaterial <G4_MATERIAL>`
- `/B2/det/stepMax <valor> <unidad>`

## 10. Limitaciones técnicas actuales (importantes)

1. `CMakeLists.txt` hereda una lista de scripts del ejemplo B2a original (`exampleB2a.out`, `exampleB2.in`, `gui.mac`, `run1.mac`, `run2.mac`, `init_vis.mac`, `vis.mac`) que no están presentes en este repositorio.  
   En el estado actual, `cmake -S . -B build` falla hasta ajustar esa sección.

2. `run_dedx.mac` define energías de `/gun/energy`, pero el generador primario personalizado vuelve a fijar energía y partícula en cada evento.  
   Resultado: el barrido de energías del macro no controla la cinemática efectiva mientras `PrimaryGeneratorAction` conserve la lógica actual.

3. La columna `ekin` del ntuple no contiene energía cinética real en la implementación actual.

4. El branch `particle` puede quedar vacío si no se propaga explícitamente el nombre de partícula al hit.

## 11. Interpretación física esperada

Si la estadística es suficiente y se corrigen/entienden los puntos anteriores, se espera observar:

- separación por signo en `dE/dx vs p/q` (`mu-` a la izquierda, `mu+` a la derecha);
- tendencia `1/beta^2` a bajos `beta*gamma`;
- región MIP alrededor de `beta*gamma ~ 3-4`;
- crecimiento relativista suave a altos `beta*gamma`;
- distribuciones de depósito con cola tipo Landau.

## 12. Referencias de código

- Ejecutable y physics list: [`exampleB2a.cc`](exampleB2a.cc)
- Geometría y materiales: [`src/DetectorConstruction.cc`](src/DetectorConstruction.cc)
- Generación primaria: [`src/PrimaryGeneratorAction.cc`](src/PrimaryGeneratorAction.cc)
- Detector sensible: [`src/TrackerSD.cc`](src/TrackerSD.cc)
- Reducción por evento: [`src/EventAction.cc`](src/EventAction.cc)
- Escritura ROOT: [`src/RunAction.cc`](src/RunAction.cc)
- Análisis y gráficas: [`analysis/plot_dedx.C`](analysis/plot_dedx.C)

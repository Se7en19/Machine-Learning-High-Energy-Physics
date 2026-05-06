# Simulación μ⁺ — Bar Strip Detector

Geant4 para μ⁺ en el detector de barras centelladoreas con absorbedor de 70 cm de hierro.

## Compilar y ejecutar

```bash
cd build
cmake ..
make -j4
./sim ../barrido_continuo.mac
```

## Geometría

- Gun: μ⁺ desde z = −2 m, dirección +z
- Absorbedor: Fe 70×70×70 cm, z = 0–70 cm (≈4.17 λ_I)
- Capa 1: 20 barras BC404 (1m×5cm×1cm) a lo largo de X, z = 100.5 cm
- Capa 2: 20 barras BC404 (1m×5cm×1cm) a lo largo de Y, z = 103.5 cm

## Salida

`../../../Classifier/data/muon_bars/output_runN.root`  (N = 0–79)

TTree `Hits`: fX, fY, fZ, fEdep, fdEdx, Ekin, TOF, TrackLength, ScatteringAng, Momentum, fEvent, layerID, barID

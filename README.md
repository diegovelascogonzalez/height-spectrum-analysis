# Height Spectrum (v4) — PN-only membrane midplane fluctuations

Compute the **membrane midplane height fluctuation spectrum** ⟨|h(q)|²⟩ from an MD trajectory that contains **only P and N atoms** (e.g., `PN.dcd` + `PN.pdb`).  
This version bins modes by **Fourier index rings** (k-index radius) *per frame*, averages within each ring, and then computes **time-averaged spectra** and **error bars**.

---

## Features

- Computes midplane height spectrum using leaflet FFTs:
  - FFT top and bottom leaflet height fields independently
  - Midplane Fourier mode: **h(q) = (h_top(q) + h_bottom(q))/2**
- Grids atoms to a 2D lattice and fills empty cells via **weighted neighbor averaging**
- Optional detrending: remove per-leaflet spatial mean (**REMOVE_MEAN** toggle)
- Two error estimates for the binned spectrum:
  - **Flyvbjerg–Petersen (blocking)** error estimate
  - **Fixed M-block** error estimate (default `M = 10`)
- Writes both:
  - **Unbinned** spectrum (legacy/continuity output)
  - **Binned** spectrum (recommended)

---

## Requirements

### Software
- A C++17 compiler (e.g., `g++ >= 7`, `clang++ >= 7`)
- [Chemfiles](https://chemfiles.org/) (trajectory I/O)
- [FFTW3](https://www.fftw.org/) (2D FFTs)

### Libraries
You must be able to link:
- `-lchemfiles`
- `-lfftw3`
- `-lm`

---

## Build

From the folder containing `heightspectrum_v4.cc`:

```bash
g++ heightspectrum_v4.cc -o heightspectrum -std=c++17 -lchemfiles -lfftw3 -lm

# Height Spectrum v4  
Midplane Fluctuation Analysis for PN-Only Membrane Trajectories

C++ program to compute the membrane midplane height fluctuation spectrum  
⟨|h(q)|²⟩ from molecular dynamics trajectories containing only phosphate (P) and nitrogen (N) atoms.

---

## 1. What does this code do?

This program:

- Reads `PN.dcd` and `PN.pdb` using Chemfiles
- Constructs 2D height fields for top and bottom leaflets
- Performs 2D FFTs using FFTW3
- Computes midplane Fourier modes:

  h(q) = (h_top(q) + h_bottom(q)) / 2

- Bins modes by Fourier index rings
- Computes time-averaged spectra
- Computes statistically robust error bars using:
  - Flyvbjerg–Petersen blocking
  - Fixed M-block error estimation
- Outputs publication-ready spectrum files

This implementation is designed for large membrane systems and PN-only trajectories.

---

## 2. Features

- Leaflet-resolved FFT
- Midplane spectrum construction
- Ring-based k-index binning
- Optional spatial mean removal (detrending)
- Flyvbjerg–Petersen error analysis
- Fixed M-block error estimate
- Handles large systems efficiently
- Skips invalid frames safely

---

## 3. Requirements

### Compiler

- C++17 compatible compiler (e.g., g++ ≥ 7)

### Libraries

- Chemfiles
- FFTW3

Link flags required:

    -lchemfiles -lfftw3 -lm

---

## 4. Build Instructions

From the directory containing `heightspectrum_v4.cc`:

```bash
g++ heightspectrum_v4.cc -o heightspectrum -std=c++17 -lchemfiles -lfftw3 -lm
```

If libraries are not found at runtime:

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

Adjust according to your system or module environment.

---

## 5. Usage

This version does NOT use command-line arguments.  
All configuration is done directly in the source file.

### Step 1 — Edit configuration block

Modify these variables near the top of `heightspectrum_v4.cc`:

```cpp
const std::string folder          = "/path/to/data/";
const std::string trajectory_file = folder + "PN.dcd";
const std::string pdb_file        = folder + "PN.pdb";

const double L_grid   = 16.44830556;
const bool   centered = true;
```

### Step 2 — Run

```bash
./heightspectrum
```

---

## 6. Inputs

Required files:

- PN.dcd — trajectory file
- PN.pdb — topology file

Both must be readable by Chemfiles.

### IMPORTANT: Leaflet Assignment

Atoms must be ordered such that:

- First half of atoms → top leaflet
- Second half of atoms → bottom leaflet

Leaflet assignment is index-based.

---

## 7. Outputs

All output files are written to the working directory.

### A.dat

Mean projected membrane area.

Column 1:
- ⟨Lx · Ly⟩

---

### height_spectrum_avg.dat

Unbinned time-averaged spectrum.

Column 1:
- q (1/length)

Column 2:
- ⟨|h(q)|²⟩ (length⁴)

---

### height_spectrum_avg_binned.dat  (Recommended)

Ring-binned spectrum (includes q = 0).

Column 1:
- q

Column 2:
- ⟨|h(q)|²⟩

---

### height_spectrum_avg_error_binned.dat

Flyvbjerg–Petersen blocking error (q = 0 excluded).

Column 1:
- q

Column 2:
- ⟨|h(q)|²⟩

Column 3:
- σ_FP(q)

---

### height_spectrum_avg_error_fix_blocking_binned.dat

Fixed M-block error estimate (default M = 10).

Column 1:
- q

Column 2:
- ⟨|h(q)|²⟩

Column 3:
- σ_fixed(q)

---

## 8. Minimal Example

Compile:

```bash
g++ heightspectrum_v4.cc -o heightspectrum -std=c++17 -lchemfiles -lfftw3 -lm
```

Run:

```bash
./heightspectrum
```

Plot example (gnuplot):

```bash
gnuplot -e "set logscale xy; plot 'height_spectrum_avg_error_binned.dat' u 1:2:3 w yerrorbars"
```

---

## 9. Units

If input coordinates are in length units L:

- q has units 1/L
- ⟨|h(q)|²⟩ has units L⁴

---

## 10. Error Analysis

### Flyvbjerg–Petersen (1989)

Blocking procedure accounting for time correlation in MD trajectories.

### Fixed M-Block

Divides trajectory into M equal blocks (default M = 10).

---

## 11. Reproducibility Notes

Results depend on:

- Grid spacing (L_grid)
- Whether coordinates are centered
- Detrending toggle
- Number of frames
- Blocking parameters
- Box dimensions per frame

The code:

- Skips frames with invalid box dimensions
- Handles uneven block sizes in fixed blocking
- Uses internal convergence safeguards for FP blocking

---

## 12. Post-Processing

Recommended file for bending modulus fitting:

    height_spectrum_avg_error_binned.dat

Typical workflow:

- Select low-q regime
- Fit q⁴ ⟨|h(q)|²⟩ scaling
- Extract bending modulus κ

---

## 13. Citation

If you use this code, please cite:

Flyvbjerg, H. & Petersen, H. G.  
Error estimates on averages of correlated data.  
J. Chem. Phys. 91, 461 (1989).

Also cite:

- FFTW3
- Chemfiles

Example BibTeX:

```bibtex
@article{FlyvbjergPetersen1989,
  title = {Error estimates on averages of correlated data},
  journal = {J. Chem. Phys.},
  year = {1989},
  volume = {91},
  pages = {461--466}
}
```

---

## 14. Code Availability Statement (for manuscripts)

Height-spectrum analysis was performed using an in-house C++ implementation based on Chemfiles and FFTW3 (this repository).

---

## 15. Acknowledgments

This work used resources of the Oak Ridge Leadership Computing Facility at Oak Ridge National Laboratory, supported by the U.S. Department of Energy Office of Science under Contract No. DE-AC05-00OR22725. Access to the Frontier supercomputer was provided through an approved allocation.

---

## 16. License

Add a LICENSE file before public release (MIT or BSD-3-Clause recommended).

---

## 17. Author

Diego L. Velasco  
University of Delaware  
Lyman Lab

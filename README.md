# EVPA Rotation Analysis

**A Python package for detecting and analyzing Electric Vector Position Angle (EVPA) rotations in polarimetric monitoring data of Active Galactic Nuclei (AGN).**

---

## Table of Contents

1. [Features](#features)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Data Format](#data-format)
5. [Methodology](#methodology)
6. [Output](#output)
7. [Examples](#examples)
8. [Testing](#testing)
9. [Citation & Contact](#citation--contact)
10. [References](#references)

---

## 1. Features

* **Data Loading & Preprocessing**

  * Automatic filtering: PD/σ\_PD ≥ 3 \[1]
  * MJD conversion & gap segmentation

* **Angle Ambiguity Resolution**

  * Error-weighted 180° wrap optimization \[2]

* **Bayesian Blocks Segmentation**

  * `astropy.stats.bayesian_blocks` with measurement uncertainties \[3]

* **Rotation Detection**

  * Local extrema pairing (max↔min)
  * Statistical tests: Student’s t-test \[4], binomial test \[5]

* **Visualization**

  * Publication-ready EVPA & Bayesian Blocks plots

* **CLI Support**

  * Batch or single-source processing
  * Customizable thresholds (amplitude, p-values)

---

## 2. Installation

**Via pip:**

```bash
pip install -e .
```

**Via conda:**

```bash
conda env create -f environment.yml
conda activate evpa-rotation
```

**Development mode:**

```bash
git clone https://github.com/glykanastasia/blazar-evpa-rotation-detector.git
cd evpa-rotation-analysis
pip install -e .[dev,docs,plotting]
```

---

## 3. Quick Start

### 3.1 Python API

```python
from evpa_rotation import RotationAnalyzer

analyzer = RotationAnalyzer(data_file="monitoring_data.csv")
results = analyzer.analyze_all_sources(
    p0=0.001, diff_threshold=0, 
    t_test_threshold=0.05, binom_threshold=0.0625
)
print(analyzer.get_summary_statistics())
```

### 3.2 Command Line

```bash
# Full dataset
evpa-analyze --data monitoring_data.csv --output results/

# Specific sources & stricter criteria
evpa-analyze --data monitoring_data.csv \
  --sources "3C 454.3" "CTA 102" \
  --min-amp 30 --t-threshold 0.01 --binom-threshold 0.01 \
  --show-plots
```

---

## 4. Data Format

* **Data Source:** Monitoring data compiled by Blinov et al. (2021) \[6]

| Column          | Description                   |
| --------------- | ----------------------------- |
| `J2000_name`    | Source identifier             |
| `Julian_date`   | Observation date (Julian Day) |
| `EVPA[deg]`     | EVPA in degrees               |
| `err_EVPA[deg]` | EVPA uncertainty (degrees)    |
| `PD[%]`         | Polarization degree (%)       |
| `err_PD[%]`     | Polarization uncertainty (%)  |

\-------------------------------------|
\| `J2000_name`     | Source identifier                   |
\| `Julian_date`    | Observation date (Julian Day)       |
\| `EVPA[deg]`      | EVPA in degrees                     |
\| `err_EVPA[deg]`  | EVPA uncertainty (degrees)          |
\| `PD[%]`          | Polarization degree (%)             |
\| `err_PD[%]`      | Polarization uncertainty (%)        |

---

## 5. Methodology

1. **Preprocessing**

   * Convert to MJD, filter PD/σ\_PD ≥ 3 \[1], segment gaps (30 d).

2. **Angle Adjustment**

   * Minimize error-weighted Δψ between points, resolve 180° wraps \[2].

3. **Bayesian Blocks**

   * Identify change points with heteroscedastic errors \[3].

4. **Rotation Detection**

   * Pair extrema → compute Δψ, Δt.
   * **t-test:** H₀: EVPA distribution centered on 180° \[4].
   * **Binomial test:** random orientation around 180° \[5].

5. **Criteria**

   * p\_t ≤ 0.05, p\_b ≤ 0.0625, Δψ ≥ 0° (customizable)

---

## 6. Output

* **CSV**: `rotation_events.csv` (detailed event table)
* **Plots**: PNG per event (time series & blocks)
* **Summary**: aggregated statistics

Columns: Source, Period, t\_start, t\_end, ΔMJD, AMPLITUDE, AMPLITUDE\_error, t-test p, Binom test p, n\_points, extrema\_type

---

## 7. Examples

### Single Source

```python
analyzer = RotationAnalyzer(data_file="data.csv")
rotations = analyzer.analyze_source("3C 454.3", show_plots=True)
```

### Custom Parameters

```python
results = analyzer.analyze_all_sources(
  p0=1e-4, diff_threshold=45,
  t_test_threshold=1e-2, binom_threshold=1e-2
)
```

### Plotting

```python
from evpa_rotation import plot_source_overview

plot_source_overview(analyzer.load_data("3C 454.3"), "3C 454.3")
```

---

## 8. Testing

```bash
pytest tests/             # Run tests
pytest --cov=evpa_rotation  # Coverage report
```

---

## 9. Citation & Contact

**Citation:**

> Glykopoulou Anastasiia (2025). *EVPA Rotation Analysis: blazar-evpa-rotation-detector*. GitHub repository. Commit d4dbd2b (Aug 7, 2025). [https://github.com/glykanastasia/blazar-evpa-rotation-detector](https://github.com/glykanastasia/blazar-evpa-rotation-detector)

**Contact & Support:**

* **Developer:** Glykopoulou Anastasiia (Bachelor Student)
  Department of Physics, University of Crete
  Institute of Astrophysics, Foundation for Research and Technology – Hellas (FORTH)
  Email: [glykopoulouanastasi@gmail.com](mailto:glykopoulouanastasi@gmail.com)
* **GitHub Issues:** [https://github.com/glykanastasia/blazar-evpa-rotation-detector/issues](https://github.com/glykanastasia/blazar-evpa-rotation-detector/issues)

---

## 10. References

1. Angel, J. R. P. & Stockman, H. S. (1980). *Optical and Infrared Polarization of Active Extragalactic Objects*. Annual Review of Astronomy and Astrophysics, 18, 321–361. [https://doi.org/10.1146/annurev.aa.18.090180.001541](https://doi.org/10.1146/annurev.aa.18.090180.001541)
2. Hovatta, T., et al. (2016). *Discrete correlation functions and EVPA swings in blazars*. Astronomy & Astrophysics, 596, A78. [https://doi.org/10.1051/0004-6361/201628967](https://doi.org/10.1051/0004-6361/201628967)
3. Scargle, J. D., et al. (2013). *Bayesian Block Representations of Poisson Data*. The Astrophysical Journal, 764, 167. [https://doi.org/10.1088/0004-637X/764/2/167](https://doi.org/10.1088/0004-637X/764/2/167)
4. Blinov, D., et al. (2015). *RoboPol: polarization rotations in blazars*. Monthly Notices of the Royal Astronomical Society, 453, 1669–1683. [https://doi.org/10.1093/mnras/stv1698](https://doi.org/10.1093/mnras/stv1698)
5. Myserlis, I., et al. (2018). *Polarimetric monitoring and EVPA rotations in AGN*. Astronomy & Astrophysics, 614, A123. [https://doi.org/10.1051/0004-6361/201732039](https://doi.org/10.1051/0004-6361/201732039)
6. Blinov, D., et al. (2021). *Long-term optical polarization monitoring of blazars with RoboPol*. Monthly Notices of the Royal Astronomical Society, **XXX**, YYY–ZZZ. [https://doi.org/10.1093/mnras/abcd1234](https://doi.org/10.1093/mnras/abcd1234)

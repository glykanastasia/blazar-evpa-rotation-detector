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
  * Automatic filtering: PD/σ_PD ≥ 3 [1]
  * MJD conversion & gap segmentation

* **Angle Ambiguity Resolution**
  * Error-weighted 180° wrap optimization [2]

* **Bayesian Blocks Segmentation**
  * `astropy.stats.bayesian_blocks` with measurement uncertainties [3]

* **Rotation Detection**
  * Local extrema pairing (max↔min)
  * Statistical tests: Student's t-test [4], binomial test [5]

* **Visualization**
  * Publication-ready EVPA & Bayesian Blocks plots

* **CLI Support**
  * Batch or single-source processing
  * Customizable thresholds (amplitude, p-values)

---

## 2. Installation

### Option 1: Download ZIP (Quick Start - No Git Required)

**Perfect for quick testing or if you don't have Git installed:**

1. **Download:** Click the green "Code" button on GitHub → "Download ZIP"
2. **Extract:** Unzip the downloaded file
3. **Setup environment:**
   ```bash
   cd blazar-evpa-rotation-detector-main
   conda env create -f environment.yml
   conda activate evpa-rotation
   ```
4. **Run analysis:**
   ```bash
   python analyze.py
   ```

⚠️ **Important:** You must be in the extracted directory to run the code.

---

### Option 2: Via Git Clone + pip

```bash
git clone https://github.com/glykanastasia/blazar-evpa-rotation-detector.git
cd blazar-evpa-rotation-detector
pip install -e .
```

---

### Option 3: Via conda (from cloned repo)

```bash
git clone https://github.com/glykanastasia/blazar-evpa-rotation-detector.git
cd blazar-evpa-rotation-detector
conda env create -f environment.yml
conda activate evpa-rotation
```

---

### Option 4: Development mode (for contributors)

```bash
git clone https://github.com/glykanastasia/blazar-evpa-rotation-detector.git
cd blazar-evpa-rotation-detector
pip install -e .[dev,docs,plotting]
```

---

## 3. Quick Start

### 3.1 Python API

```python
from evpa_rotation import RotationAnalyzer

analyzer = RotationAnalyzer(data_file="data/monitoring_data.csv")
results = analyzer.analyze_all_sources(
    p0=0.001, diff_threshold=90,
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

* **Data Source:** Monitoring data compiled by Blinov et al. (2021) [6]

| Column          | Description                   |
|-----------------|-------------------------------|
| `J2000_name`    | Source identifier             |
| `Julian_date`   | Observation date (Julian Day) |
| `EVPA[deg]`     | EVPA in degrees               |
| `err_EVPA[deg]` | EVPA uncertainty (degrees)    |
| `PD[%]`         | Polarization degree (%)       |
| `err_PD[%]`     | Polarization uncertainty (%)  |

---

## 5. Methodology

1. **Preprocessing**
   * Convert to MJD, filter PD/σ_PD ≥ 3 [1], segment gaps (30 d).

2. **Angle Adjustment**
   * Minimize error-weighted Δψ between points, resolve 180° wraps [2].

3. **Bayesian Blocks**
   * Identify change points with heteroscedastic errors [3].

4. **Rotation Detection**
   * Pair extrema → compute Δψ, Δt.
   * **t-test:** H₀: EVPA distribution centered on 180° [4].
   * **Binomial test:** random orientation around 180° [5].

5. **Criteria**
   * p_t ≤ 0.05, p_b ≤ 0.0625, Δψ ≥ 0° (customizable)

---

## 6. Output

* **CSV**: `rotation_events.csv` (detailed event table)
* **Plots**: PNG per event (time series & blocks)
* **Summary**: aggregated statistics

Columns: Source, Period, t_start, t_end, ΔMJD, AMPLITUDE, AMPLITUDE_error, t-test p, Binom test p, n_points, extrema_type

---

## 7. Examples

### Single Source

```python
analyzer = RotationAnalyzer(data_file="data.csv")
rotations = analyzer.analyze_source("3C 454.3", show_plots=True)
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

plot_source_overview(analyzer.load_data("3C 454.3"), "3C 454.3")
```

---

## 8. Testing

```bash
pytest tests/             # Run tests
pytest --cov=evpa_rotation  # Coverage report
```

---

## 9. Citation & Contact

### How to Cite This Software

**If you use this software in your research, please cite it!** This helps support continued development and allows others to discover the tool.

#### BibTeX Entry

```bibtex
@software{glykopoulou2025evpa,
  author = {Glykopoulou, Anastasia},
  title = {EVPA Rotation Analysis: A Python package for detecting EVPA rotations in AGN},
  year = {2025},
  url = {https://github.com/glykanastasia/blazar-evpa-rotation-detector},
  version = {1.0.0}
}
```

#### Text Citation

> Glykopoulou, A. (2025). *EVPA Rotation Analysis: blazar-evpa-rotation-detector* (Version 1.0.0) [Computer software]. https://github.com/glykanastasia/blazar-evpa-rotation-detector

#### APA Style

> Glykopoulou, A. (2025). *EVPA Rotation Analysis: blazar-evpa-rotation-detector* (Version 1.0.0) [Computer software]. GitHub. https://github.com/glykanastasia/blazar-evpa-rotation-detector

---

### Acknowledgments

If this software contributed to your published work, we kindly request that you:

1. **Cite the software** using one of the formats above
2. **Mention it in your acknowledgments section**, for example:

> "This research made use of the EVPA Rotation Analysis package (Glykopoulou, 2025), available at https://github.com/glykanastasia/blazar-evpa-rotation-detector"

3. **Reference the underlying methodology** from:
   - Blinov et al. (2015) for rotation detection criteria (MNRAS, 453, 1669)
   - Blinov et al. (2021) for the monitoring data compilation (MNRAS, 501, 3715)

---

### Contact & Support

* **Developer:** Glykopoulou Anastasia (Bachelor Student)  
  Department of Physics, University of Crete  
  Institute of Astrophysics, Foundation for Research and Technology – Hellas (FORTH)  
  Email: [glykopoulouanastasia@gmail.com](mailto:glykopoulouanastasia@gmail.com)

* **Report Issues:** [GitHub Issues](https://github.com/glykanastasia/blazar-evpa-rotation-detector/issues)

* **Contributions:** Pull requests are welcome! Please see our contributing guidelines.

---

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 10. References



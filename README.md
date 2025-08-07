# EVPA Rotation Analysis

A comprehensive Python package for detecting and analyzing Electric Vector Position Angle (EVPA) rotations in polarimetric monitoring data of astronomical sources, particularly blazars and other Active Galactic Nuclei (AGN).

## Features

* **Data Loading & Preprocessing**
  Load CSV monitoring data with automatic quality filtering (PD/σ\_PD ≥ 3) \[1].

* **Angle Adjustment**
  Resolve the 180° EVPA ambiguity via an error-weighted optimization that minimizes sequential angular jumps \[2].

* **Bayesian Blocks Analysis**
  Segment time series into statistically significant intervals using the Bayesian blocks algorithm with heteroscedastic uncertainties \[3].

* **Rotation Detection**
  Identify rotation events by pairing adjacent extrema (max→min or min→max) and applying Student’s t-test and binomial test to assess significance \[4, 5].

* **Publication-Quality Plots**
  Generate high-resolution visualizations of EVPA, Bayesian blocks, and detected rotations.

* **Command-Line Interface**
  Batch-process entire datasets or individual sources with customizable amplitude and p-value thresholds.

## Installation

### Using pip

```bash
pip install -e .
```

### Using conda

```bash
conda env create -f environment.yml
conda activate evpa-rotation
```

### Development Installation

```bash
git clone https://github.com/yourusername/evpa-rotation-analysis.git
cd evpa-rotation-analysis
pip install -e .[dev,docs,plotting]
```

## Quick Start

### Python API

```python
from evpa_rotation import RotationAnalyzer

# Initialize analyzer
analyzer = RotationAnalyzer(data_file="monitoring_data.csv")

# Analyze all sources
results = analyzer.analyze_all_sources(
    p0=0.001,                    # Bayesian blocks prior [3]
    diff_threshold=0,            # Minimum amplitude (°)
    t_test_threshold=0.05,       # t-test p-value threshold [4]
    binom_threshold=0.0625       # Binomial test p-value threshold [5]
)

# Get summary statistics
summary = analyzer.get_summary_statistics()
print(summary)
```

### Command Line Interface

```bash
# Analyze all sources in data file
evpa-analyze --data monitoring_data.csv --output results/

# Analyze specific sources with custom thresholds
evpa-analyze --data monitoring_data.csv \
             --sources "3C 454.3" "CTA 102" \
             --min-amp 30 \
             --t-threshold 0.01 \
             --binom-threshold 0.01 \
             --show-plots
```

## Data Format

The package expects CSV files with the following columns:

* `J2000_name`: Source name
* `Julian_date`: Julian date of observation
* `EVPA[deg]`: Electric Vector Position Angle (degrees)
* `err_EVPA[deg]`: EVPA measurement error (degrees)
* `PD[%]`: Polarization degree (%)
* `err_PD[%]`: Polarization degree error (%)

## Methodology

1. **Data Preprocessing**

   * Convert Julian Date to Modified Julian Date (MJD).
   * Filter observations with PD/σ\_PD ≥ 3 \[1].
   * Segment by observation gaps (default: 30 days).

2. **EVPA Angle Adjustment**

   * Optimize 180° wraps to minimize error-weighted angular differences between consecutive points \[2].

3. **Bayesian Blocks Analysis**

   * Use `astropy.stats.bayesian_blocks` incorporating measurement uncertainties \[3].
   * Identify local maxima and minima within each block.

4. **Rotation Detection**

   * Pair adjacent extrema (max→min or min→max).
   * Compute rotation amplitude Δψ and duration Δt (days).
   * Apply Student’s t-test (H₀: distribution centered on 180°) \[4].
   * Apply binomial test for random orientation around 180° \[5].

5. **Statistical Criteria**
   Default thresholds (customizable):

   * t-test p ≤ 0.05
   * Binomial test p ≤ 0.0625
   * Minimum amplitude ≥ 0°

## Output

The analysis produces:

1. **CSV Results**: `rotation_events.csv` with detected rotation events
2. **Individual Plots**: PNG files for each rotation event
3. **Summary Statistics**: Overview of all detected rotations

## Results Columns

* `Source`
* `Period`
* `t_start`, `t_end` (MJD)
* `ΔMJD` (duration in days)
* `AMPLITUDE` (degrees)
* `AMPLITUDE_error`
* `t-test p`
* `Binom test p`
* `n_points`
* `extrema_type` (max→min or min→max)

## Examples

### Analyze Single Source

```python
from evpa_rotation import RotationAnalyzer

analyzer = RotationAnalyzer(data_file="data.csv")
rotations = analyzer.analyze_source("3C 454.3", show_plots=True)
```

### Custom Analysis Parameters

```python
# More stringent criteria
results = analyzer.analyze_all_sources(
    p0=0.0001,              # Higher sensitivity [3]
    diff_threshold=45,      # Require ≥45° amplitude
    t_test_threshold=0.01,  # More stringent t-test
    binom_threshold=0.01    # More stringent binomial test
)
```

### Plotting Functions

```python
from evpa_rotation import plot_bayesian_blocks_with_extrema, plot_source_overview

# Plot Bayesian blocks analysis
plot_bayesian_blocks_with_extrema(data, "3C 454.3")

# Plot source overview
plot_source_overview(data, "3C 454.3")
```

## Testing

Run the test suite:

```bash
pytest tests/
```

With coverage:

```bash
pytest --cov=evpa_rotation tests/
```

## Citation

If you use this package in published research, please cite:

> **YourName**, A., et al. (2025). *EVPA Rotation Analysis: A Python Toolbox for Polarimetric Monitoring of AGN*. Journal of Astronomical Software, **XX**(X), 1–12.

## Contact

For questions, feature requests, or bug reports, please contact:

* **Email**: [your.email@institution.edu](mailto:your.email@institution.edu)
* **GitHub Issues**: [https://github.com/yourusername/evpa-rotation-analysis/issues](https://github.com/yourusername/evpa-rotation-analysis/issues)

## References

1. Angel, J. R. P. & Stockman, H. S. (1980). Optical and Infrared Polarization of Active Extragalactic Objects. *Annual Review of Astronomy and Astrophysics*, 18, 321–361.
2. Hovatta, T., et al. (2016). Discrete correlation functions and EVPA swings in blazars. *Astronomy & Astrophysics*, 596, A78.
3. Scargle, J. D., et al. (2013). Bayesian Block Representations of Poisson Data. *The Astrophysical Journal*, 764, 167.
4. Blinov, D., et al. (2015). RoboPol: polarization rotations in blazars. *Monthly Notices of the Royal Astronomical Society*, 453, 1669–1683.
5. Myserlis, I., et al. (2018). Polarimetric monitoring and EVPA rotations in AGN. *Astronomy & Astrophysics*, 614, A123.

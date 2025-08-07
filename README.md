# EVPA Rotation Analysis

A comprehensive Python package for detecting and analyzing Electric Vector Position Angle (EVPA) rotations in polarimetric monitoring data of astronomical sources, particularly blazars and other Active Galactic Nuclei (AGN).

## Features

- **Data Loading & Preprocessing**: Load CSV monitoring data with automatic quality filtering
- **Angle Adjustment**: Resolve 180° EVPA ambiguity using error-weighted optimization
- **Bayesian Blocks Analysis**: Segment time series into statistically significant intervals
- **Rotation Detection**: Identify rotation events using statistical tests (t-test, binomial test)
- **Publication-Quality Plots**: Generate publication-ready visualizations
- **Command-Line Interface**: Easy-to-use CLI for batch processing

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
    p0=0.001,                    # Bayesian blocks prior
    diff_threshold=0,            # Minimum amplitude
    t_test_threshold=0.05,       # t-test p-value threshold
    binom_threshold=0.0625       # Binomial test p-value threshold
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
             --show-plots
```

## Data Format

The package expects CSV files with the following columns:

- `J2000_name`: Source name
- `Julian_date`: Julian date of observation
- `EVPA[deg]`: Electric Vector Position Angle in degrees
- `err_EVPA[deg]`: EVPA measurement error
- `PD[%]`: Polarization degree in percent
- `err_PD[%]`: Polarization degree error

## Methodology

### 1. Data Preprocessing
- Convert Julian dates to Modified Julian Date (MJD)
- Filter observations with significant polarization detection (PD/σ_PD ≥ 3)
- Segment data by observation gaps (default: 30 days)

### 2. EVPA Angle Adjustment
- Resolve 180° ambiguity using error-weighted optimization
- Minimize angular differences between consecutive measurements
- Account for measurement uncertainties

### 3. Bayesian Blocks Analysis
- Segment time series into statistically optimal intervals
- Use `astropy.stats.bayesian_blocks` with measurement errors
- Identify local maxima and minima in segmented data

### 4. Rotation Detection
- Analyze adjacent extrema pairs (max→min or min→max)
- Apply statistical tests:
  - **t-test**: Test if EVPA distribution differs from 180°
  - **Binomial test**: Test if EVPA values are randomly distributed around 180°
- Calculate rotation amplitude and duration

### 5. Statistical Criteria
Default thresholds (can be customized):
- t-test p-value ≤ 0.05
- Binomial test p-value ≤ 0.0625
- Minimum amplitude ≥ 0° (customizable)

## Output

The analysis produces:

1. **CSV Results**: `rotation_events.csv` with detected rotation events
2. **Individual Plots**: PNG files for each rotation event
3. **Summary Statistics**: Overview of all detected rotations

### Results Columns

- `Source`: Source name
- `Period`: Observation period index
- `t_start`, `t_end`: Start and end times (MJD)
- `ΔMJD`: Duration in days
- `AMPLITUDE`: EVPA rotation amplitude in degrees
- `AMPLITUDE_error`: Amplitude uncertainty
- `t-test p`: t-test p-value
- `Binom test p`: Binomial test p-value
- `n_points`: Number of data points in rotation
- `extrema_type`: Type of rotation (max→min or min→max)

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
    p0=0.0001,              # More sensitive Bayesian blocks
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

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)

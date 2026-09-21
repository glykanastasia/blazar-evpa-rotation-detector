# Automated Detection of Polarization Angle Rotations in Blazars

[![DOI](https://img.shields.io/badge/DOI-10.1051%2F0004--6361%2F202558360-blue)](https://doi.org/10.1051/0004-6361/202558360)
[![A&A](https://img.shields.io/badge/A%26A-710%2C%20A260-orange)](https://doi.org/10.1051/0004-6361/202558360)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.8-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A Python pipeline for the automated detection and statistical analysis of **Electric Vector Position Angle (EVPA) rotations** in blazars, using long-term polarimetric monitoring data.

This repository contains the analysis pipeline developed and used in:

> **Glykopoulou, A., Liodakis, I., & Blinov, D. (2026)**
> *Automating the detection of polarization angle rotations in blazars: Re-analysis of RoboPol data reveals 27 new rotations*
> **Astronomy & Astrophysics, 710, A260**
> https://doi.org/10.1051/0004-6361/202558360

The pipeline combines error-weighted EVPA ambiguity correction, Bayesian Blocks segmentation, local-extrema identification, and statistical validation into a reproducible procedure for identifying EVPA rotation events.

---

## Table of Contents

- [Overview](#overview)
- [Dataset](#dataset)
- [Method](#method)
- [Rotation Parameters](#rotation-parameters)
- [Results Reproduced by the Pipeline](#results-reproduced-by-the-pipeline)
- [Comparison with Previous RoboPol Catalogues](#comparison-with-previous-robopol-catalogues)
- [Fermi-LAT Comparison](#fermi-lat-comparison)
- [Installation](#installation)
- [Package Structure](#package-structure)
- [Usage](#usage)
- [Input Data](#input-data)
- [Output](#output)
- [Testing](#testing)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [References](#references)
- [Acknowledgements](#acknowledgements)
- [License](#license)
- [Author](#author)

---

## Overview

EVPA rotations are large, systematic changes in the polarization position angle observed in blazars. Because the EVPA is defined modulo 180°, apparent discontinuities can arise from the angle representation itself and must be corrected before searching for rotations.

The pipeline:

- preprocesses polarimetric monitoring data;
- applies an error-weighted correction of the 180° EVPA ambiguity;
- divides the data into contiguous observing segments;
- applies Bayesian Blocks to the EVPA time series;
- identifies local maxima and minima;
- constructs candidate rotations from adjacent extrema;
- applies statistical significance tests;
- calculates rotation amplitude, duration, and angular velocity;
- produces rotation catalogues and diagnostic figures.

The goal is to reduce subjective choices in identifying EVPA rotations and to provide a reproducible framework for large polarimetric monitoring datasets.

---

## Dataset

The analysis uses the **RoboPol optical polarimetric monitoring dataset** from the 2013–2017 observing campaign. It includes sources from the main and control samples, as well as additional individually observed sources.

The data supporting the published results are openly available through the **Harvard Dataverse**. <!-- TODO: add Dataverse DOI/link -->

---

## Method

### 1. Preprocessing

Monitoring data are read from CSV files containing the polarimetric measurements of each source. Julian Dates are converted to Modified Julian Dates (MJD), and observations are sorted chronologically.

Only measurements with significant polarization are retained:

$$
\frac{PD}{\sigma_{PD}} \geq 3
$$

### 2. EVPA ambiguity correction

Because EVPA is defined modulo 180°, consecutive observations can show artificial jumps of about 180°.

The pipeline applies an **error-weighted ambiguity correction**. Each new observation is compared with the previously adjusted value, taking the combined measurement uncertainties into account:

- if the angular difference is small enough, the measured value is kept;
- otherwise, integer multiples of 180° are added, and the equivalent angle giving the smallest uncertainty-weighted difference is selected.

Gaps longer than **30 days** are treated as segment boundaries, because reliable ambiguity correction across such gaps is not possible.

### 3. Time segmentation

Each light curve is split into contiguous observing segments. A segment is kept only if:

- it contains **more than four data points**;
- no two consecutive observations are separated by more than **30 days**.

Accepted segments are analysed independently.

### 4. Bayesian Blocks

The adjusted EVPA time series is segmented with the **Bayesian Blocks** algorithm (Scargle et al. 2013), implemented in `astropy.stats`. Since EVPA values are continuous measurements with Gaussian uncertainties, the pipeline uses:

```python
fitness="measures"
```

The published analysis adopts a false-alarm probability of:

```text
p0 = 0.001
```

### 5. Identification of local extrema

Local maxima and minima are found in the Bayesian Blocks representation using:

```python
scipy.signal.argrelextrema
```

The first and last blocks are checked explicitly so that boundary extrema are not missed.

### 6. Candidate rotation events

Adjacent extrema form a candidate rotation when their types alternate:

```text
maximum → minimum
minimum → maximum
```

The start and end of each candidate are set by the corresponding Bayesian Block edges.

Unlike the earlier RoboPol procedure, this method does **not impose a smoothness constraint** on the EVPA curve. Candidates are identified from significant extrema pairs, regardless of the behaviour of the curve between them.

### 7. Rotation selection criteria

A candidate is retained only if it meets all of the following:

| Criterion | Requirement |
|---|---|
| Minimum amplitude | $\Delta\theta \geq 90^\circ$ |
| Internal structure | ≥ 3 Bayesian Blocks after sub-segmentation |
| Number of measurements | ≥ 4 data points in the interval |
| One-sample Student's t-test | p-value ≤ 0.05 |
| Binomial test | p-value ≤ 0.0625 |

Both statistical tests are applied to the EVPA measurements within the candidate interval.

---

## Rotation Parameters

For each accepted event the pipeline computes:

**Amplitude**

$$
\Delta\theta = \theta_{\max} - \theta_{\min}
$$

with measurement uncertainties propagated to the amplitude uncertainty.

**Duration**

$$
\Delta t = t_{\mathrm{end}} - t_{\mathrm{start}}
$$

measured from the Bayesian Block boundaries, with uncertainty derived from the block-boundary uncertainties.

**Angular velocity**

$$
\dot{\theta} = \frac{\Delta\theta}{\Delta t}
$$

with uncertainty propagated from the amplitude and duration uncertainties.

---

## Results Reproduced by the Pipeline

Applied to the RoboPol 2013–2017 dataset, the pipeline identifies:

| Quantity | Value |
|---|---|
| EVPA rotations | 48 |
| Unique sources | 25 |
| Previously unreported rotations | 27 |
| Rotations overlapping previous RoboPol catalogues | 21 |

Sources with multiple detected events include:

| Source | Rotations |
|---|---|
| RBPLJ2232+1143 | 6 |
| RBPLJ1751+0939 | 4 |
| RBPLJ1800+7828 | 4 |
| RBPLJ2253+1608 | 4 |

Properties of the detected rotations:

| Property | Value |
|---|---|
| Amplitude range | 90.8° – 359.7° |
| Duration range | 7.0 – 111.3 days |
| Mean amplitude | 159.7° ± 66.5° |
| Mean duration | 44.4 ± 23.8 days |
| Mean rate | 5.01° day⁻¹ |

The largest amplitude (359.7°) is found in RBPLJ1751+0939 and the smallest (90.8°) in RBPLJ1635+3808.

Of the **27 previously unreported rotations**:

- **16** occur in seasons already covered by previous RoboPol studies;
- **11** occur in the final 2016–2017 season, extending the rotation search into the last RoboPol observing season.

---

## Comparison with Previous RoboPol Catalogues

The detections were cross-matched with the rotations reported by Blinov et al. (2015, 2016a,b, 2018). On average, the Bayesian Blocks events are:

| Property | Compared with previous catalogues |
|---|---|
| Amplitude | ~10% larger |
| Duration | ~2× longer |
| Angular velocity | ~⅔ slower |

The paper interprets these as methodological effects of adaptive Bayesian Blocks segmentation compared with the earlier manual/event-segmentation approach.

---

## Fermi-LAT Comparison

The paper also examines the link between EVPA rotations and contemporaneous gamma-ray activity. Seven-day binned Fermi-LAT light curves were used to compute the peak-to-mean gamma-ray energy-flux ratio during each rotation interval.

All 48 events were analysed with Spearman rank correlations and a Monte Carlo analysis of **50,000 realizations**.

| Rotation property | Spearman ρ | Median p-value |
|---|---|---|
| Duration | 0.601 ± 0.049 | 6.4 × 10⁻⁶ |
| Maximum EVPA amplitude | 0.146 ± 0.016 | 0.322 |

The duration correlation stays significant (p < 0.05) in all Monte Carlo realizations; the amplitude correlation does not.

---

## Installation

### Requirements

- Python ≥ 3.8
- NumPy
- pandas
- SciPy
- Matplotlib
- Astropy
- tqdm

Additional development and plotting dependencies are listed in `environment.yml`.

### Clone the repository

```bash
git clone https://github.com/glykanastasia/blazar-evpa-rotation-detector.git
cd blazar-evpa-rotation-detector
```

### Create and activate the Conda environment

```bash
conda env create -f environment.yml
conda activate evpa-rotation
```

### Install the package

```bash
pip install -e .
```

---

## Package Structure

```text
blazar-evpa-rotation-detector/
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── setup.py
├── src/
│   └── evpa_rotation/
│       ├── __init__.py
│       ├── angles.py
│       ├── bayesian_blocks.py
│       ├── loader.py
│       ├── main.py
│       └── plotter.py
└── tests/
    └── test_angles.py
```

---

## Usage

### Python API

The main analysis class is `RotationAnalyzer`:

```python
from evpa_rotation import RotationAnalyzer

analyzer = RotationAnalyzer(
    data_file="monitoring_data.csv",
    output_dir="rotation_analysis",
)

results = analyzer.analyze_all_sources(
    p0=0.001,
    t_test_threshold=0.05,
    binom_threshold=0.0625,
)

summary = analyzer.get_summary_statistics()
print(summary)
```

Analyse specific sources only:

```python
results = analyzer.analyze_all_sources(
    sources=["RBPLJ2232+1143", "RBPLJ1751+0939"],
    p0=0.001,
)
```

Load data on its own:

```python
from evpa_rotation import load_data

data = load_data("monitoring_data.csv")
print(data.head())
```

### Command line

After installation, the `evpa-analyze` command is available:

```bash
evpa-analyze \
    --data monitoring_data.csv \
    --output rotation_analysis \
    --p0 0.001
```

Selected sources:

```bash
evpa-analyze \
    --data monitoring_data.csv \
    --output rotation_analysis \
    --sources RBPLJ2232+1143 RBPLJ1751+0939
```

All options:

```bash
evpa-analyze --help
```

---

## Input Data

The pipeline expects a CSV table with at least the following columns (or equivalent):

| Column | Description |
|---|---|
| `J2000_name` | Source identifier |
| `MJD` | Modified Julian Date |
| `PD[%]` | Polarization degree |
| `err_PD[%]` | Polarization-degree uncertainty |
| `EVPA[deg]` | Electric Vector Position Angle |
| `err_EVPA[deg]` | EVPA uncertainty |

The exact preprocessing used for the published RoboPol analysis is described in Section 2 of the paper.

---

## Output

The analysis produces a rotation catalogue and diagnostic plots:

```text
rotation_analysis/
├── rotation_events.csv
└── plots/
```

`rotation_events.csv` includes:

- Source
- Amplitude and amplitude uncertainty
- Start MJD and end MJD
- Duration and duration uncertainty
- t-test p-value
- Binomial-test p-value
- Rotation rate and rotation-rate uncertainty

The plots show the EVPA measurements, the Bayesian Blocks representation, and the detected rotation intervals.

---

## Testing

Run the unit tests:

```bash
pytest
```

With coverage:

```bash
pytest --cov=evpa_rotation
```

---

## Reproducibility

This repository supports the reproducibility of the analysis in the associated publication. It contains the scripts used for:

- EVPA rotation detection;
- Bayesian Blocks segmentation;
- figure generation.

For the full scientific methodology and the published rotation catalogue, see the paper and its Appendix A.

---

## Citation

If you use this software in your research, please cite:

```bibtex
@article{Glykopoulou2026EVPA,
  author  = {Glykopoulou, Anastasia and Liodakis, Ioannis and Blinov, Dmitry},
  title   = {Automating the detection of polarization angle rotations in blazars:
             Re-analysis of RoboPol data reveals 27 new rotations},
  journal = {Astronomy \& Astrophysics},
  year    = {2026},
  volume  = {710},
  pages   = {A260},
  doi     = {10.1051/0004-6361/202558360}
}
```

A [`CITATION.cff`](CITATION.cff) file is also provided for software citation.

---

## References

<details>
<summary>Show references</summary>

- Abdo, A. A., Ackermann, M., Agudo, I., et al. 2010, *ApJ*, 716, 30
- Abdollahi, S., Ajello, M., Baldini, L., et al. 2023, *ApJS*, 265, 31
- Angelakis, E., Hovatta, T., Blinov, D., et al. 2016, *MNRAS*, 463, 3365
- Astropy Collaboration (Price-Whelan, A. M., et al.) 2022, *ApJ*, 935, 167
- Blandford, R., Meier, D., & Readhead, A. 2019, *ARA&A*, 57, 467
- Blinov, D. & Pavlidou, V. 2019, *Galaxies*, 7, 46
- Blinov, D., Pavlidou, V., Papadakis, I., et al. 2015, *MNRAS*, 453, 1669
- Blinov, D., Pavlidou, V., Papadakis, I., et al. 2016a, *MNRAS*, 462, 1775
- Blinov, D., Pavlidou, V., Papadakis, I. E., et al. 2016b, *MNRAS*, 457, 2252
- Blinov, D., Pavlidou, V., Papadakis, I., et al. 2018, *MNRAS*, 474, 1296
- Blinov, D., Kiehlmann, S., Pavlidou, V., et al. 2020, *MNRAS*, 501, 3715
- Britzen, S., Fendt, C., Witzel, G., et al. 2018, *MNRAS*, 478, 3199
- Di Gesu, L., Marshall, H. L., Ehlert, S. R., et al. 2023, *Nature Astronomy*, 7, 1245
- Hovatta, T. & Lindfors, E. 2019, *New Astronomy Reviews*, 87, 101541
- Kiehlmann, S., Blinov, D., Pearson, T. J., & Liodakis, I. 2017, *MNRAS*, 472, 3589
- King, O. G., Blinov, D., Ramaprakash, A. N., et al. 2014, *MNRAS*, 442, 1706
- Maksym, W. P., Liodakis, I., Saade, M. L., et al. 2025, *ApJ*, 986, 230
- Marscher, A. P. 2014, *ApJ*, 780, 87
- Marscher, A. P., Jorstad, S. G., Larionov, V. M., et al. 2008, *Nature*, 452, 966
- Otero-Santos, J., Acosta-Pulido, J. A., Becerra González, J., et al. 2023, *MNRAS*, 523, 4504
- Panopoulou, G., Tassis, K., Blinov, D., et al. 2015, *MNRAS*, 452, 715
- Pavlidou, V., Angelakis, E., Myserlis, I., et al. 2014, *MNRAS*, 442, 1693
- Raiteri, C. M., Villata, M., Acosta-Pulido, J. A., et al. 2017, *Nature*, 552, 374
- Savchenko, S. S., Morozova, D. A., Jorstad, S. G., et al. 2024, *Astrophysical Bulletin*, 79, 186
- Scargle, J. D., Norris, J. P., Jackson, B., & Chiang, J. 2013, *ApJ*, 764, 167

</details>

---

## Acknowledgements

This work uses the RoboPol optical polarimetric monitoring data.

A.G. and I.L. were funded by the European Union ERC-2022-STG BOOTES project (grant agreement No. 101076343). D.B. acknowledges support from the European Research Council under grant agreement No. 101040021.

---

## License

This project is distributed under the **MIT License**. See [`LICENSE`](LICENSE) for the full text.

---

## Author

**Anastasia Glykopoulou**

Department of Physics, University of Crete, GR-70013 Heraklion, Greece
Institute of Astrophysics, Foundation for Research and Technology–Hellas (FORTH), GR-70013 Heraklion, Greece

Email: glykopoulouanastasia@gmail.com
Repository: https://github.com/glykanastasia/blazar-evpa-rotation-detector

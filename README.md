# Automated Detection of Polarization Angle Rotations in Blazars

A Python pipeline for the automated detection and statistical analysis of **Electric Vector Position Angle (EVPA) rotations** in blazars using long-term polarimetric monitoring data.

This repository contains the analysis pipeline developed and used in:

> **Glykopoulou, A., Liodakis, I., & Blinov, D. (2026)**
> *Automating the detection of polarization angle rotations in blazars*
> *Re-analysis of RoboPol data reveals 27 new rotations*
> **Astronomy & Astrophysics, 710, A260**
> https://doi.org/10.1051/0004-6361/202558360

The pipeline combines error-weighted EVPA ambiguity correction, Bayesian Blocks segmentation, local-extrema identification, and statistical validation to provide a reproducible procedure for identifying EVPA rotation events.

---

## Overview

EVPA rotations are large and systematic changes in the polarization position angle observed in blazars. Because EVPA is defined modulo \(180^\circ\), apparent discontinuities can arise from the angle representation itself and must be corrected before searching for rotations.

The pipeline developed in this work provides an automated analysis framework that:

* preprocesses polarimetric monitoring data;
* applies an error-weighted correction of the \(180^\circ\) EVPA ambiguity;
* divides the data into contiguous observing segments;
* applies Bayesian Blocks to the EVPA time series;
* identifies local maxima and minima;
* constructs candidate rotations from adjacent extrema;
* applies statistical significance tests;
* calculates rotation amplitude, duration, and angular velocity;
* produces rotation catalogues and diagnostic figures.

The goal is to reduce subjective choices in the identification of EVPA rotations and provide a reproducible framework for large polarimetric monitoring datasets.

---

## Dataset

The analysis presented in the paper uses the **RoboPol optical polarimetric monitoring dataset**, covering the 2013–2017 observing campaign.

The data include the sources observed by RoboPol during the monitoring programme, including sources from the main and control samples as well as additional individually observed sources.

The data supporting the published results are openly available through the **Harvard Dataverse**.

---

## Method

### 1. Preprocessing

The monitoring data are read from CSV files containing the polarimetric measurements of each source.

Julian Dates are converted to Modified Julian Dates (MJD), and observations are sorted chronologically.

To retain sufficiently significant polarization measurements, the analysis applies:

$$
\frac{PD}{\sigma_{PD}} \geq 3.
$$

Only measurements satisfying this criterion are retained for the subsequent analysis.

---

### 2. EVPA ambiguity correction

EVPA measurements are defined modulo \(180^\circ\). Therefore, consecutive observations can show artificial jumps of approximately \(180^\circ\).

The pipeline applies an **error-weighted ambiguity correction**. For each new observation, the measured EVPA is compared with the previously adjusted value while accounting for the combined measurement uncertainties.

If the angular difference is sufficiently small, the measured value is retained. Otherwise, integer multiples of \(180^\circ\) are added and the equivalent angle producing the smallest uncertainty-weighted difference is selected.

Observing gaps larger than **30 days** are treated as natural segment boundaries because reliable ambiguity correction across such gaps is not possible.

---

### 3. Time segmentation

Each source light curve is divided into contiguous observing segments.

A segment is retained only when:

* it contains **more than four data points**;
* no consecutive observations are separated by more than **30 days**.

These accepted segments are analysed independently.

---

### 4. Bayesian Blocks

The adjusted EVPA time series is segmented using the **Bayesian Blocks** algorithm of Scargle et al. (2013), implemented through `astropy.stats`.

Because EVPA measurements are treated as continuous measurements with Gaussian uncertainties, the pipeline uses:

```python
fitness="measures"
```

The analysis presented in the paper adopts:

```text
p0 = 0.001
```

for the Bayesian Blocks false-alarm probability parameter.

---

### 5. Identification of local extrema

Local maxima and minima are identified from the Bayesian-Blocks representation of the EVPA time series.

The implementation uses:

```python
scipy.signal.argrelextrema
```

and explicitly checks the first and last Bayesian Blocks so that boundary extrema are not missed.

---

### 6. Candidate rotation events

Adjacent extrema are considered candidate rotations when their types alternate:

```text
maximum → minimum
minimum → maximum
```

The temporal boundaries of a candidate event are defined by the corresponding Bayesian-Block edges.

Unlike the earlier RoboPol identification procedure, the method does **not impose a smoothness constraint** on the EVPA curve. Candidate rotations are identified from statistically significant extrema pairs irrespective of the derivative behaviour of the intervening light curve.

---

### 7. Rotation selection criteria

A candidate event is retained only when it satisfies the following conditions.

#### Minimum amplitude

The EVPA amplitude must satisfy

$$
\Delta\theta \geq 90^\circ.
$$

This criterion is used to select large EVPA variations.

#### Internal Bayesian-Blocks structure

The candidate interval must contain at least **three Bayesian Blocks** after sub-segmentation.

#### Minimum number of measurements

At least **four data points** must fall within the candidate interval.

#### Statistical validation

Each candidate is tested using:

1. a one-sample Student's t-test;
2. a binomial test.

The statistical tests are applied to the EVPA measurements within the candidate interval.

The adopted thresholds in the analysis are:

```text
t-test p-value      ≤ 0.05
binomial p-value    ≤ 0.0625
```

---

## Rotation parameters

For each accepted event, the pipeline calculates the following quantities.

### Amplitude

The rotation amplitude is determined from the difference between the maximum and minimum EVPA values:

$$
\Delta\theta =
\theta_{\max}-\theta_{\min}.
$$

Measurement uncertainties are propagated to obtain the amplitude uncertainty.

### Duration

The rotation duration is measured from the Bayesian-Blocks boundaries:

$$
\Delta t = t_{\rm end}-t_{\rm start}.
$$

The uncertainty is derived from the uncertainties associated with the block boundaries.

### Angular velocity

The rotation rate is calculated as

$$
\dot{\theta} =
\frac{\Delta\theta}{\Delta t}.
$$

The uncertainty is propagated from the uncertainties in amplitude and duration.

---

## Results reproduced by the pipeline

Application of the pipeline to the RoboPol 2013–2017 dataset identified:

```text
48 EVPA rotations
25 unique sources
27 previously unreported rotations
21 rotations overlapping previous RoboPol catalogues
```

Several sources contain multiple detected events, including:

```text
RBPLJ2232+1143       6 rotations
RBPLJ1751+0939       4 rotations
RBPLJ1800+7828       4 rotations
RBPLJ2253+1608       4 rotations
```

The detected rotations span:

```text
Amplitude:       90.8° – 359.7°
Duration:        7.0 – 111.3 days
Mean amplitude:  159.7° ± 66.5°
Mean duration:   44.4 ± 23.8 days
Mean rate:       5.01° day⁻¹
```

The largest amplitude detected is \(359.7^\circ\) for RBPLJ1751+0939, while the smallest is \(90.8^\circ\) for RBPLJ1635+3808.

The pipeline identified **27 previously unreported rotations**:

```text
16 events in seasons already covered by previous RoboPol studies
11 events in the final 2016–2017 season
```

The latter extends the rotation search into the final RoboPol observing season.

---

## Comparison with previous RoboPol catalogues

The published study cross-matched the Bayesian-Blocks detections with rotations reported by Blinov et al. (2015, 2016a,b, 2018).

The comparison found systematic differences between the two methodologies. On average, the Bayesian-Blocks events were:

```text
~10% larger in amplitude
~2 times longer in duration
~2/3 slower in angular velocity
```

These differences were interpreted in the paper as methodological effects associated with adaptive Bayesian-Blocks segmentation compared with the previous manual/event-segmentation approach.

---

## Fermi–LAT comparison

The paper also investigates the relationship between EVPA rotations and contemporaneous gamma-ray activity.

Seven-day binned Fermi–LAT light curves were used to calculate the peak-to-mean gamma-ray energy-flux ratio during each EVPA rotation interval.

The analysis considered all 48 detected events and used Spearman rank correlations together with a Monte Carlo analysis of **50,000 realizations**.

For rotation duration:

$$
\rho = 0.601 \pm 0.049,
\qquad
p_{\rm med}=6.4\times10^{-6}.
$$

For maximum EVPA amplitude:

$$
\rho = 0.146 \pm 0.016,
\qquad
p_{\rm med}=0.322.
$$

The paper reports that the duration correlation remains significant in all Monte Carlo realizations at \(p<0.05\), while the amplitude correlation does not.

---

## Installation

### Requirements

The package requires:

```text
Python >= 3.8
NumPy
pandas
SciPy
Matplotlib
Astropy
tqdm
```

Additional development and plotting dependencies are provided in `environment.yml`.

### Clone the repository

```bash
git clone https://github.com/glykanastasia/blazar-evpa-rotation-detector.git
cd blazar-evpa-rotation-detector
```

### Create the Conda environment

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate evpa-rotation
```

### Install the package

```bash
pip install -e .
```

---

## Package structure

```text
blazar-evpa-rotation-detector/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── setup.py
│
├── src/
│   └── evpa_rotation/
│       ├── __init__.py
│       ├── angles.py
│       ├── bayesian_blocks.py
│       ├── loader.py
│       ├── main.py
│       └── plotter.py
│
└── tests/
    └── test_angles.py
```

---

## Python API

The main analysis class is:

```python
from evpa_rotation import RotationAnalyzer
```

A basic analysis can be run with:

```python
from evpa_rotation import RotationAnalyzer

analyzer = RotationAnalyzer(
    data_file="monitoring_data.csv",
    output_dir="rotation_analysis"
)

results = analyzer.analyze_all_sources(
    p0=0.001,
    t_test_threshold=0.05,
    binom_threshold=0.0625
)

summary = analyzer.get_summary_statistics()

print(summary)
```

A specific set of sources can be analysed with:

```python
results = analyzer.analyze_all_sources(
    sources=[
        "RBPLJ2232+1143",
        "RBPLJ1751+0939"
    ],
    p0=0.001
)
```

Data can also be loaded independently:

```python
from evpa_rotation import load_data

data = load_data("monitoring_data.csv")

print(data.head())
```

---

## Command-line usage

After installation, the package provides the `evpa-analyze` command.

Example:

```bash
evpa-analyze \
    --data monitoring_data.csv \
    --output rotation_analysis \
    --p0 0.001
```

Selected sources can be analysed with:

```bash
evpa-analyze \
    --data monitoring_data.csv \
    --output rotation_analysis \
    --sources RBPLJ2232+1143 RBPLJ1751+0939
```

See all available options with:

```bash
evpa-analyze --help
```

---

## Input data

The analysis pipeline expects a CSV table containing, at minimum, information equivalent to:

| Column          | Description                     |
| --------------- | ------------------------------- |
| `J2000_name`    | Source identifier               |
| `MJD`           | Modified Julian Date            |
| `PD[%]`         | Polarization degree             |
| `err_PD[%]`     | Polarization-degree uncertainty |
| `EVPA[deg]`     | Electric Vector Position Angle  |
| `err_EVPA[deg]` | EVPA uncertainty                |

The exact preprocessing requirements used for the published RoboPol analysis are described in Section 2 of the paper.

---

## Output

The analysis produces a rotation catalogue and diagnostic plots.

A typical output directory is:

```text
rotation_analysis/
├── rotation_events.csv
└── plots/
```

The rotation catalogue contains quantities including:

```text
Source
Amplitude
Amplitude uncertainty
Start MJD
End MJD
Duration
Duration uncertainty
t-test p-value
Binomial-test p-value
Rotation rate
Rotation-rate uncertainty
```

The pipeline also generates plots showing the EVPA measurements, Bayesian-Blocks representation, and detected rotation intervals.

---

## Testing

The repository contains unit tests under:

```text
tests/
```

Run the tests with:

```bash
pytest
```

For coverage:

```bash
pytest --cov=evpa_rotation
```

---

## Reproducibility

The software repository is intended to support reproducibility of the analysis presented in the associated publication.

The paper explicitly makes the analysis pipeline available through this repository, including the scripts used for:

* EVPA rotation detection;
* Bayesian-Blocks segmentation;
* figure generation.

For the complete description of the scientific methodology and the published rotation catalogue, refer to the paper and its Appendix A.

---

## Citation

Please cite the following paper when using this software in research:

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

The repository also contains a `CITATION.cff` file for software citation.

---

## References

Abdo, A. A., Ackermann, M., Agudo, I., et al. 2010, *ApJ*, 716, 30

Abdollahi, S., Ajello, M., Baldini, L., et al. 2023, *ApJS*, 265, 31

Angelakis, E., Hovatta, T., Blinov, D., et al. 2016, *MNRAS*, 463, 3365

Astropy Collaboration (Price-Whelan, A. M., et al.) 2022, *ApJ*, 935, 167

Blandford, R., Meier, D., & Readhead, A. 2019, *ARA&A*, 57, 467

Blinov, D. & Pavlidou, V. 2019, *Galaxies*, 7, 46

Blinov, D., Pavlidou, V., Papadakis, I., et al. 2015, *MNRAS*, 453, 1669

Blinov, D., Pavlidou, V., Papadakis, I., et al. 2016a, *MNRAS*, 462, 1775

Blinov, D., Pavlidou, V., Papadakis, I. E., et al. 2016b, *MNRAS*, 457, 2252

Blinov, D., Pavlidou, V., Papadakis, I., et al. 2018, *MNRAS*, 474, 1296

Blinov, D., Kiehlmann, S., Pavlidou, V., et al. 2020, *MNRAS*, 501, 3715

Britzen, S., Fendt, C., Witzel, G., et al. 2018, *MNRAS*, 478, 3199

Di Gesu, L., Marshall, H. L., Ehlert, S. R., et al. 2023, *Nature Astronomy*, 7, 1245

Hovatta, T. & Lindfors, E. 2019, *New Astronomy Reviews*, 87, 101541

Kiehlmann, S., Blinov, D., Pearson, T. J., & Liodakis, I. 2017, *MNRAS*, 472, 3589

King, O. G., Blinov, D., Ramaprakash, A. N., et al. 2014, *MNRAS*, 442, 1706

Maksym, W. P., Liodakis, I., Saade, M. L., et al. 2025, *ApJ*, 986, 230

Marscher, A. P. 2014, *ApJ*, 780, 87

Marscher, A. P., Jorstad, S. G., Larionov, V. M., et al. 2008, *Nature*, 452, 966

Otero-Santos, J., Acosta-Pulido, J. A., Becerra González, J., et al. 2023, *MNRAS*, 523, 4504

Panopoulou, G., Tassis, K., Blinov, D., et al. 2015, *MNRAS*, 452, 715

Pavlidou, V., Angelakis, E., Myserlis, I., et al. 2014, *MNRAS*, 442, 1693

Raiteri, C. M., Villata, M., Acosta-Pulido, J. A., et al. 2017, *Nature*, 552, 374

Savchenko, S. S., Morozova, D. A., Jorstad, S. G., et al. 2024, *Astrophysical Bulletin*, 79, 186

Scargle, J. D., Norris, J. P., Jackson, B., & Chiang, J. 2013, *ApJ*, 764, 167

---

## Acknowledgements

This work uses the RoboPol optical polarimetric monitoring data.

A.G. and I.L. were funded by the European Union ERC-2022-STG BOOTES project (grant agreement No. 101076343). D.B. acknowledges support from the European Research Council under grant agreement No. 101040021.

---

## License

This project is distributed under the **MIT License**.

See [`LICENSE`](LICENSE) for the full license text.

---

## Author

**Anastasia Glykopoulou**

Department of Physics
University of Crete
GR-70013 Heraklion, Greece

Institute of Astrophysics
Foundation for Research and Technology–Hellas (FORTH)
GR-70013 Heraklion, Greece

Email: `glykopoulouanastasia@gmail.com`

Repository:

https://github.com/glykanastasia/blazar-evpa-rotation-detector

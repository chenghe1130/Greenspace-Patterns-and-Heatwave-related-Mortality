
# Small Accessible Urban Parks Mitigate Heat-related Mortality

[![License](https://img.shields.io/badge/License-Research-blue.svg)]()
[![R Version](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)]()
[![DOI](https://img.shields.io/badge/DOI-pending-orange)]()

## Overview

This repository contains the statistical analysis code for our paper:

> **Small accessible urban parks mitigate heat-related mortality**  
> *He, C., Yin, P., Chen, R., et al.*

**Key Finding**: Urban green space configuration—not just total coverage—plays a crucial role in reducing heatwave-related mortality. Cities with more dispersed, accessible green spaces and regular boundaries show significantly lower mortality risks during heatwaves.

### What this code does

This analysis examines how different green space landscape patterns modify the relationship between heatwaves and mortality across 265 Chinese cities (2013-2019). We use:

- **Stage 1**: Distributed lag non-linear models (DLNM) to quantify city-specific heatwave-mortality associations
- **Stage 2**: Random-effects meta-analysis to assess how five landscape metrics modify these associations
- **Key metrics**: PLAND (coverage), SPLIT (dispersion), LSI (boundary regularity), PD (patch density), LPI (largest patch size)

---

## Requirements

### Software
- R (≥ 4.0.0)
- Required packages: `mgcv`, `dlnm` (≥2.4.7), `mvmeta` (≥1.0.3), `splines`, `dplyr`

```r
install.packages(c("mgcv", "dlnm", "mvmeta", "splines", "dplyr"))
```

### Data
You will need two datasets:

1. **Daily mortality and meteorological data** (`ccches_heatwave.csv`)
   - Columns: `code_city`, `death`, `temp`, `rh`, `dow`, `time`, `season`, heatwave indicators
   
2. **City-level green space metrics** (`city_metadata_greenspace.csv`)
   - Columns: `code_city`, `PLAND`, `PD`, `LPI`, `LSI`, `SPLIT` (categorized as 1/2/3 for Low/Medium/High)

---

## Advanced Usage

### Sensitivity Analyses

```r
# Test different lag periods
results_lag14 <- analyze_heatwave_mortality(data, lag_days = 14)

# Test different degrees of freedom
results_sens <- analyze_heatwave_mortality(data, df_time = 6, df_humidity = 6)
```

### Stratified Analyses

```r
# Age-specific analysis
data_young <- data_main %>% mutate(death = death_5_64)
results_young <- analyze_heatwave_mortality(data_young)

data_old <- data_main %>% mutate(death = death_65plus)
results_old <- analyze_heatwave_mortality(data_old)

# Climate zone-specific analysis
data_temperate <- subset(data, code_city %in% temperate_cities)
results_temperate <- analyze_heatwave_mortality(data_temperate)
```

---

## Results

### Output File
`heatwave_greenspace_results_all.csv` contains:

| Column | Description |
|--------|-------------|
| `Heatwave_Definition` | 92.5%_2day / 92.5%_3day / 95%_2day / 95%_3day |
| `Landscape_Metric` | PLAND / PD / LPI / LSI / SPLIT |
| `Category` | Low / Medium / High |
| `RR` | Relative risk (point estimate) |
| `Lower_CI` | Lower 95% confidence interval |
| `Upper_CI` | Upper 95% confidence interval |
| `P_value` | Significance of difference vs. Low category |

### Interpretation
- **RR > 1**: Increased mortality during heatwaves
- **P < 0.05**: Significant difference between categories
- Compare High vs Low categories for each metric

---

## Methods

### Statistical Models

**Stage 1 - City-specific GAM**:
```r
death ~ crossbasis(heatwave, lag=10) + factor(dow) + ns(rh, 4) + ns(time, 12)
```
- Family: quasi-Poisson
- Crossbasis: B-spline (var, 4 df) + natural spline (lag, 4 df, log-knots)
- Controls: day of week, humidity (4 df), time trend (4 df/year)

**Stage 2 - Meta-analysis**:
```r
mvmeta(coef ~ 1, vcov, method = "reml")
```
- Random-effects model
- REML estimation
- Stratified by landscape metric categories

---

## Citation

If you use this code, please cite:

```bibtex
@article{he2024greenspace,
  title={Small accessible urban parks mitigate heat-related mortality},
  author={He, Cheng and Yin, Peng and Chen, Renjie and others},
  journal={[Journal Name]},
  year={2026},
  doi={[DOI]}
}
```

---

## Contact

**Lead Author**:
- Cheng He: chenghe@hsph.harvard.edu

For code-related questions, please open an issue on GitHub.

---

## License

This code is provided for research and academic purposes. For commercial use, please contact the authors.

---

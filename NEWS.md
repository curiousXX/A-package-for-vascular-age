# VasAge 0.1.0

## Initial Release

### Features
- First public release of the `VasAge` R package.
- Implements vascular age estimation using KDM (Klemera–Doubal Method) and PCA (Principal Component Analysis).
- Core functions include:
  - `Vas_calcu()`: calculates individual vascular age using internal standardization parameters.
- `Vas_train()`: allows users to train their own vascular age models.
- `data_process()`: assists with data cleaning, log transformation, and missing value handling.

### Data
- Includes an example dataset `test_data` (from a Chinese cohort with over 4500 participants), which can be loaded using `data("test_data")` for demonstration and testing purposes.
- Internal standardization parameters are stored in `sysdata.rda` and used exclusively within package functions.

### Documentation
- Comprehensive function documentation with reproducible examples is provided.

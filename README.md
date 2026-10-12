# AMDP — Absolute Magnetic Data Processing

Python software package to process ADS/DIM absolute magnetic measurement data
and calculate spot baseline values.

Authors: László Merényi, Ádám Domján — Mingeo Ltd.

## Requirements

- Python 3.10+ (the code uses `match` statements)
- Dependencies listed in [`requirements.txt`](requirements.txt):

  ```bash
  python -m pip install -r requirements.txt
  ```

Optional dependencies, needed only for some input formats:

| Package | Needed for |
|---|---|
| `geomagpy` (MagPy) | `DIM-data-files-format: MagPy-autodetect` (also requires `magpy_wrapper.py`, see below) |
| `cdflib` | reading CDF data files |

## Features
- Can process non standard geomagnetic absolute measurements:
   - Handles any number of measurement points/samples (ADS case) per absolute measurement
   - Azimuth data can be read from separate data files (ADS case)
   - Different methods/solvers to get spot baselines: residual, Gauss-Newton, Levenberg–Marquardt algorithm
- Command line usage with YAML configuration files: for most use scenarios, calling Python functions and modifying the Python source code is not required
- Can run calculations with several parameter sets at once, for easy and fast comparison of results.


## Configuration files

The program is driven by a set of YAML configuration files:

| File | Purpose |
|---|---|
| `amdp-run-parameters.yaml` | Run settings: site/pillar/instrument definition, data-sets selection, calculation parameters, calculation variants |
| `piers-data.yaml` | Observatories/test sites pier data |
| `observatory-datasets.yaml` | Observatory magnetic data set (data file locations) definitions |
| `ui-preferences.yaml` | Display/output preferences (plot options, decimal precision, optional data cache directory, etc.) |

Example versions of all four files are in [`sample-config-files/`](sample-config-files/).
Copy them to your own configuration directory and adapt them.

`amdp-run-parameters.yaml`, `piers-data.yaml` and `observatory-datasets.yaml`
contain machine- and project-specific settings and are excluded from version
control via `.gitignore`. `ui-preferences.yaml` is tracked in the repository
as a shared default.

Notes on `amdp-run-parameters.yaml`:

- `calculation-type` must be `spot-baselines` (the only method supported in this version).
- Supported values of `DIM-data-files-format`: `ADS-raw-dat-files` (default),
  `Mingeo-manual-XLSX`, `MagPy-autodetect`.
  `MagPy-autodetect` needs the `magpy_wrapper.py` module, which is not yet
  included in this repository; this format is therefore not available yet.
- Data file paths starting with `./` are resolved relative to the directory of
  `amdp-run-parameters.yaml`, not to the current directory.

## Usage

Run from the project directory, optionally giving the directory that contains the
configuration files:

```bash
python amdp.py                      # configuration files in the current directory
python amdp.py path/to/config-dir
```

The program reads amdp-run-parameters.yaml from the given configuration directory, or from the current directory if none is given. 
The other three files are also looked up there first; if one is not found, it is taken from the current directory instead.

Individual configuration files can be overridden on the command line:

```bash
python amdp.py path/to/config-dir \
  --run-parameters path/to/amdp-run-parameters.yaml \
  --piers-data path/to/piers-data.yaml \
  --observatory-datasets path/to/observatory-datasets.yaml \
  --ui-preferences path/to/ui-preferences.yaml
```

To run the program from a Python IDE (e.g. Spyder), set `project_path` in
[`amdp-runner.py`](amdp-runner.py) and run that file.

## Project structure

Main program:

- `amdp.py` — main entry point
- `amdp-runner.py` — helper to start `amdp.py` from a Python IDE

DIM/ADS measurement processing:

- `dim_data_process.py` — DIM/ADS measurement processing, spot baseline calculation
- `dim_measurement_loader.py` — loading DIM/ADS measurement data
- `dim_instr.py` — software model of the ADS theodolite
- `dim_azimuth.py` — azimuth / true North calculations
- `dim_utils.py` — DIM/ADS-specific variables and helper functions
- `amdp_base_formulas.py` — core absolute measurement formulas

Input data:

- `dat_file_reader.py` — reading time-stamped records from ASCII DAT (incl. IAGA-2002), SC and CDF files
- `sc_file_reader.py` — reading maglin SC compressed binary files
- `absmeas_xlsx_import.py` — importing absolute measurements from XLSX tables
- `magpy_wrapper.py` — interface to the MagPy library (not yet included in this repository)

Geomagnetic data and utilities:

- `geomag_utils.py`, `geomag_utils_preferences.py` — geomagnetic utilities and preferences
- `geomag_obs_dataset.py` — observatory data set handling
- `didd_utils.py` — DIDD-related utilities

Supporting modules:

- `angle_utils.py`, `trigonometric_utils.py`, `model_utils.py` — angle, trigonometric and model-fitting utilities
- `plot_utils.py` — plotting helpers
- `cli_utils.py` — console output helpers
- `misc_utils.py` — miscellaneous shared utilities (incl. Excel date/time conversions)

## License

This project is licensed under the MIT License — see [`LICENSE`](LICENSE) for details.

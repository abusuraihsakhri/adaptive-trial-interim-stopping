# Adaptive Trial Interim Stopping

A Python library, command-line tool, and browser calculator for group-sequential efficacy boundaries and interim-monitoring calculations in clinical trials.

The statistical core uses the canonical joint-normal/Brownian-motion model to calibrate two-sided O'Brien-Fleming and Pocock efficacy boundaries for the requested alpha level, number of looks, and information fractions. It also provides Lan-DeMets alpha-spending functions, conditional-power calculations, and simple sample-size re-estimation utilities.

## Features

- Two-sided O'Brien-Fleming and Pocock group-sequential efficacy boundaries.
- Arbitrary strictly increasing information fractions.
- Lan-DeMets O'Brien-Fleming-type and Pocock-type alpha spending.
- Conditional power from an interim Z statistic under an explicit future-effect assumption.
- Variance-ratio blinded and effect-ratio unblinded sample-size re-estimation utilities.
- CSV batch processing from the command line.
- Static browser calculator that runs locally and does not require a backend.
- Optional FastAPI service for programmatic use.

## Browser application

The browser application is in `web/`. All calculations run in the browser; entered values are not sent to a server by the application. The only browser storage used is `localStorage` for the light/dark theme preference.

## Installation

Python 3.10 or later is required.

```bash
git clone https://github.com/abusuraihsakhri/adaptive-trial-interim-stopping.git
cd adaptive-trial-interim-stopping
python -m pip install -e .
```

For development tests:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

For the optional API server:

```bash
python -m pip install -e ".[server]"
adaptive-trial serve
```

## Command-line examples

Calculate five equally spaced O'Brien-Fleming boundaries:

```bash
adaptive-trial boundaries --method obf --alpha 0.05 --looks 5
```

Use custom information fractions:

```bash
adaptive-trial boundaries --method pocock --alpha 0.05 --looks 4 --fractions 0.2,0.45,0.7,1
```

Calculate a Lan-DeMets spending schedule:

```bash
adaptive-trial spending --method obf --alpha 0.05 --fractions 0.25,0.5,0.75,1
```

Calculate conditional power from an interim Z statistic:

```bash
adaptive-trial futility --interim-z 1.4 --information-fraction 0.5 --threshold 0.1
```

Batch-process the included design examples:

```bash
adaptive-trial batch -i sample.csv -o results.csv
```

Use `--json` with the `boundaries`, `spending`, `futility`, or `ssr` commands for machine-readable output.

## Statistical notes

Boundary calibration controls the requested two-sided Type I error across the supplied interim looks under the canonical group-sequential model. O'Brien-Fleming boundaries are represented by a constant Brownian-scale boundary, producing stringent early Z thresholds that decrease toward the final analysis. Pocock boundaries use a constant Z threshold across looks. Lan-DeMets spending is implemented separately and is useful when the actual information fractions differ from the originally planned timing.

The conditional-power and sample-size re-estimation functions are compact analytical utilities, not a complete adaptive-design package. They do not replace a protocol-specific simulation study, statistical analysis plan, independent data monitoring process, or regulatory review. Trial adaptations should be prospectively specified and their operating characteristics evaluated before use in a confirmatory study.

## Technology and compatibility

The core library uses only the Python standard library. The optional API requires FastAPI and Uvicorn. The browser application uses plain HTML, CSS, and JavaScript with no external runtime assets and is intended for current versions of Chrome, Edge, Firefox, and Safari.

## License

MIT License. See `LICENSE`.

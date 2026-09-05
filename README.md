# Adaptive Trial Interim Stopping & Boundary Decision Support

A Python biostatistics library and CLI tool for adaptive clinical trial design, group sequential interim monitoring, and sample size re-estimation. Implements O'Brien-Fleming and Pocock stopping boundaries, Lan-DeMets alpha spending functions, conditional power estimation for futility termination, and blinded/unblinded sample size adjustments under FDA Adaptive Clinical Trial Guidelines.

Requires Python standard library only (zero external runtime dependencies).

---

## Features

- **Group Sequential Stopping Boundaries:**
  - **O'Brien-Fleming (OBF):** Conservative early efficacy boundaries preserving overall type I error rate $\alpha$.
  - **Pocock Boundaries:** Uniform critical z-value boundaries across all interim looks.
- **Lan-DeMets Alpha Spending:** Flexible spending approaches ($\alpha^*(t)$) accommodating irregular look intervals and varying information fractions.
- **Futility & Conditional Power Monitoring:** Computes conditional power given interim observed effect sizes and remaining information fraction to evaluate early futility stopping rules.
- **Sample Size Re-estimation (SSR):**
  - **Blinded SSR:** Adjusts total sample size based on pooled variance inflation.
  - **Unblinded SSR:** Recalibrates target enrollment based on observed interim effect sizes.
- **Multi-Agent Adaptive Trial Coordinator:** Evaluates multi-parameter clinical trial telemetry alerts across risk boundaries.
- **Tabular Batch Processing:** Batch evaluation of trial interim looks and task records via CSV.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies. `pytest` is optional for running the unit tests.

```bash
git clone https://github.com/abusuraihsakhri/adaptive-trial-interim-stopping.git
cd adaptive-trial-interim-stopping
```

## Environment Variables

| Variable | Description | Required |
|----------|-------------|----------|
| `AUDIT_SECRET_KEY` | Secret key for HMAC-SHA256 audit trail signing. If not set, an ephemeral key is generated at runtime and a warning is emitted. | No (but recommended for production) |

---

## CLI Usage

### 1. Single Task / Look Evaluation
Run interim audit on trial parameters:
```bash
python -m adaptive_trial.cli audit --task-id LOOK-01 --target ARM-B --primary 29.4 --secondary 15.1
```
Output as JSON:
```bash
python -m adaptive_trial.cli audit --task-id LOOK-01 --target ARM-B --primary 29.4 --secondary 15.1 --json
```

### 2. Batch CSV Processing
Process trial interim roster from CSV:
```bash
python -m adaptive_trial.cli batch -i sample.csv -o results.csv
```

### 3. Supervisory Query
Query configuration and trial guidelines:
```bash
python -m adaptive_trial.cli chat "What standard is applied for stopping boundaries?"
```

---

## Python API Quickstart

```python
from adaptive_boundaries import (
    OBFlemingBoundary,
    PocockBoundary,
    LanDeMetsSpending,
    FutilityAssessor,
    SampleSizeReestimator,
)

# 1. Compute 5-look O'Brien-Fleming boundaries (alpha = 0.05)
obf = OBFlemingBoundary(total_alpha=0.05, num_looks=5).compute_boundaries()
for look in obf:
    print(f"Look {look.look_number}: Z-bound = {look.z_bound}, alpha-spent = {look.alpha_spent}")

# 2. Assess futility with conditional power at 50% information fraction
futility = FutilityAssessor().assess(
    observed_effect=0.35,
    target_power=0.80,
    information_fraction=0.50,
    futility_threshold=0.10,
)
print(f"Futile: {futility.is_futile} | Conditional Power: {futility.conditional_power:.2%}")

# 3. Blinded sample size re-estimation
ssr = SampleSizeReestimator().blinded_reestimate(
    initial_n=200,
    current_variance=1.4,
    expected_variance=1.0,
)
print(f"Revised N: {ssr.revised_n} (Inflation: {ssr.inflation_factor})")
```

---

## Security Features

- **PHI Outbound Guard:** Detects and blocks protected health information (PHI) including MRNs, SSNs, email addresses, phone numbers, and patient names before any outbound transmission.
- **HMAC-SHA256 Audit Trail:** Cryptographically signed, tamper-evident audit log chain for all task evaluations. Set `AUDIT_SECRET_KEY` environment variable for production deployments.
- **Input Validation:** All boundary calculation and re-estimation functions validate parameters and raise descriptive errors for invalid inputs.

---

## Running Tests

Run the test suite using standard `pytest`:

```bash
pytest -v
```

The test suite includes:
- `tests/test_adaptive_trial.py` — Core adaptive trial functionality, CLI, and boundary calculations
- `tests/test_enrichment.py` — Enrichment engine feature tests
- `tests/test_security_and_validation.py` — Security (PHI guard, audit trail) and input validation tests


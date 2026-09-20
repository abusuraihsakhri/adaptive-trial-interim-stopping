#!/usr/bin/env python3
"""Statistical utilities for group-sequential interim monitoring."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
import math
from statistics import NormalDist
from typing import List, Optional, Sequence, Tuple

_NORMAL = NormalDist()
_SQRT_2PI = math.sqrt(2.0 * math.pi)


class BoundaryMethod(str, Enum):
    OBRIEN_FLEMING = "O'Brien-Fleming"
    POCOCK = "Pocock"
    LAN_DEMETS = "Lan-DeMets"


@dataclass(frozen=True)
class StoppingBoundary:
    method: str
    look_number: int
    information_fraction: float
    z_bound: float
    alpha_spent: float
    cumulative_alpha: float
    stopping_rule: str = "efficacy"


@dataclass(frozen=True)
class SpendingFunctionResult:
    function_name: str
    alpha_at_looks: List[float]
    cumulative_alpha: List[float]
    total_alpha_spent: float


@dataclass(frozen=True)
class FutilityAnalysis:
    is_futile: bool
    conditional_power: float
    predicted_success_probability: float
    recommendation: str
    monitoring_threshold: float
    interim_z: float
    final_critical_z: float


@dataclass(frozen=True)
class SampleSizeReestimate:
    original_n: int
    revised_n: int
    inflation_factor: float
    reason: str


def _validate_alpha(alpha: float, name: str = "alpha") -> None:
    if not math.isfinite(alpha) or not 0.0 < alpha < 1.0:
        raise ValueError(f"{name} must be between 0 and 1, got {alpha}")


def _validate_information_fractions(values: Sequence[float]) -> Tuple[float, ...]:
    if not values:
        raise ValueError("information_fractions must not be empty")
    fractions = tuple(float(value) for value in values)
    if any(not math.isfinite(t) or not 0.0 < t <= 1.0 for t in fractions):
        raise ValueError("information_fractions must be in (0, 1]")
    if any(right <= left for left, right in zip(fractions, fractions[1:])):
        raise ValueError("information_fractions must be strictly increasing")
    return fractions


def _normal_pdf(x: float, sd: float) -> float:
    z = x / sd
    return math.exp(-0.5 * z * z) / (_SQRT_2PI * sd)


def _brownian_limit(constant: float, t: float, method: str) -> float:
    if method == "obf":
        return constant
    if method == "pocock":
        return constant * math.sqrt(t)
    raise ValueError(f"Unknown boundary method: {method}")


def _trapezoid(values: Sequence[float], step: float) -> float:
    return step * (0.5 * values[0] + sum(values[1:-1]) + 0.5 * values[-1])


def _survival_path(
    constant: float,
    information_fractions: Tuple[float, ...],
    method: str,
    grid_points: int,
) -> List[float]:
    """P(no efficacy-boundary crossing through each look) under H0."""
    first_t = information_fractions[0]
    limit = _brownian_limit(constant, first_t, method)
    step = 2.0 * limit / (grid_points - 1)
    xs = [-limit + i * step for i in range(grid_points)]
    density = [_normal_pdf(x, math.sqrt(first_t)) for x in xs]
    survival = [_trapezoid(density, step)]
    previous_t = first_t

    for t in information_fractions[1:]:
        previous_xs = xs
        previous_density = density
        previous_step = step
        transition_sd = math.sqrt(t - previous_t)

        limit = _brownian_limit(constant, t, method)
        step = 2.0 * limit / (grid_points - 1)
        xs = [-limit + i * step for i in range(grid_points)]
        density = []

        for x in xs:
            total = 0.0
            for j, (previous_x, previous_value) in enumerate(zip(previous_xs, previous_density)):
                weight = 0.5 if j in (0, grid_points - 1) else 1.0
                total += weight * previous_value * _normal_pdf(x - previous_x, transition_sd)
            density.append(total * previous_step)

        survival.append(_trapezoid(density, step))
        previous_t = t

    return survival


@lru_cache(maxsize=128)
def _calibrate_constant(
    alpha: float,
    information_fractions: Tuple[float, ...],
    method: str,
    grid_points: int = 161,
) -> float:
    """Calibrate a symmetric boundary so total crossing probability is alpha."""
    _validate_alpha(alpha)
    low = 0.1
    high = max(4.0, _NORMAL.inv_cdf(1.0 - alpha / 2.0) + 1.5)

    def crossing_probability(constant: float) -> float:
        return 1.0 - _survival_path(
            constant, information_fractions, method, grid_points
        )[-1]

    while crossing_probability(high) > alpha:
        high *= 1.25
        if high > 12.0:
            raise RuntimeError("Could not bracket the group-sequential boundary")

    for _ in range(30):
        midpoint = (low + high) / 2.0
        if crossing_probability(midpoint) > alpha:
            low = midpoint
        else:
            high = midpoint
    return (low + high) / 2.0


class _CanonicalBoundary:
    method_key: str
    method_label: str

    def __init__(
        self,
        total_alpha: float = 0.05,
        num_looks: int = 5,
        information_fractions: Optional[Sequence[float]] = None,
        grid_points: int = 161,
    ):
        _validate_alpha(total_alpha, "total_alpha")
        if num_looks < 1:
            raise ValueError(f"num_looks must be >= 1, got {num_looks}")
        if grid_points < 81 or grid_points % 2 == 0:
            raise ValueError("grid_points must be an odd integer >= 81")

        if information_fractions is None:
            fractions = tuple(i / num_looks for i in range(1, num_looks + 1))
        else:
            fractions = _validate_information_fractions(information_fractions)
            if not math.isclose(fractions[-1], 1.0, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError("group-sequential information_fractions must end at 1.0")
            num_looks = len(fractions)

        self.total_alpha = float(total_alpha)
        self.num_looks = int(num_looks)
        self.information_fractions = fractions
        self.grid_points = int(grid_points)

    def compute_boundaries(self) -> List[StoppingBoundary]:
        constant = _calibrate_constant(
            self.total_alpha,
            self.information_fractions,
            self.method_key,
            self.grid_points,
        )
        survival = _survival_path(
            constant,
            self.information_fractions,
            self.method_key,
            self.grid_points,
        )
        cumulative = [max(0.0, min(1.0, 1.0 - value)) for value in survival]
        incremental = [cumulative[0]] + [
            max(0.0, cumulative[i] - cumulative[i - 1])
            for i in range(1, len(cumulative))
        ]

        rows = []
        for index, t in enumerate(self.information_fractions):
            z_boundary = (
                constant / math.sqrt(t)
                if self.method_key == "obf"
                else constant
            )
            rows.append(
                StoppingBoundary(
                    method=self.method_label,
                    look_number=index + 1,
                    information_fraction=round(t, 6),
                    z_bound=round(z_boundary, 4),
                    alpha_spent=round(incremental[index], 7),
                    cumulative_alpha=round(cumulative[index], 7),
                )
            )
        return rows


class OBFlemingBoundary(_CanonicalBoundary):
    """Two-sided O'Brien-Fleming efficacy boundaries under the canonical model."""

    method_key = "obf"
    method_label = "O'Brien-Fleming"


class PocockBoundary(_CanonicalBoundary):
    """Two-sided Pocock efficacy boundaries under the canonical model."""

    method_key = "pocock"
    method_label = "Pocock"


class LanDeMetsSpending:
    """Lan-DeMets cumulative alpha-spending functions."""

    @staticmethod
    def obf_spending(t: float, alpha: float = 0.05) -> float:
        _validate_alpha(alpha)
        if t <= 0.0:
            return 0.0
        if t > 1.0:
            raise ValueError("information fraction t must be <= 1")
        z_alpha = _NORMAL.inv_cdf(1.0 - alpha / 2.0)
        return 2.0 * (1.0 - _NORMAL.cdf(z_alpha / math.sqrt(t)))

    @staticmethod
    def pocock_spending(t: float, alpha: float = 0.05) -> float:
        _validate_alpha(alpha)
        if t <= 0.0:
            return 0.0
        if t > 1.0:
            raise ValueError("information fraction t must be <= 1")
        return alpha * math.log(1.0 + (math.e - 1.0) * t)

    def compute_looks(
        self,
        information_fractions: Sequence[float],
        alpha: float = 0.05,
        method: str = "OBF",
    ) -> SpendingFunctionResult:
        _validate_alpha(alpha)
        fractions = _validate_information_fractions(information_fractions)
        method_key = method.strip().upper().replace("_", "-")
        if method_key in {"OBF", "O'BRIEN-FLEMING", "OBRIEN-FLEMING"}:
            spending_function = self.obf_spending
            label = "O'Brien-Fleming"
        elif method_key == "POCOCK":
            spending_function = self.pocock_spending
            label = "Pocock"
        else:
            raise ValueError("method must be 'OBF' or 'Pocock'")

        cumulative = [spending_function(t, alpha) for t in fractions]
        incremental = [cumulative[0]] + [
            max(0.0, cumulative[i] - cumulative[i - 1])
            for i in range(1, len(cumulative))
        ]
        return SpendingFunctionResult(
            function_name=f"Lan-DeMets {label}-type",
            alpha_at_looks=[round(value, 8) for value in incremental],
            cumulative_alpha=[round(value, 8) for value in cumulative],
            total_alpha_spent=round(cumulative[-1], 8),
        )


class FutilityAssessor:
    """Conditional power for a non-binding futility rule."""

    @staticmethod
    def assess_from_z(
        interim_z: float,
        information_fraction: float,
        final_critical_z: float = 1.959963984540054,
        assumed_final_mean_z: Optional[float] = None,
        futility_threshold: float = 0.10,
    ) -> FutilityAnalysis:
        if not math.isfinite(interim_z):
            raise ValueError("interim_z must be finite")
        if not 0.0 < information_fraction <= 1.0:
            raise ValueError("information_fraction must be in (0, 1]")
        if not math.isfinite(final_critical_z) or final_critical_z <= 0.0:
            raise ValueError("final_critical_z must be positive and finite")
        if not 0.0 <= futility_threshold <= 1.0:
            raise ValueError("futility_threshold must be in [0, 1]")
        if assumed_final_mean_z is not None and not math.isfinite(assumed_final_mean_z):
            raise ValueError("assumed_final_mean_z must be finite")

        t = information_fraction
        if t == 1.0:
            conditional_power = 1.0 if interim_z >= final_critical_z else 0.0
        else:
            drift = (
                interim_z / math.sqrt(t)
                if assumed_final_mean_z is None
                else assumed_final_mean_z
            )
            conditional_mean = math.sqrt(t) * interim_z + drift * (1.0 - t)
            conditional_sd = math.sqrt(1.0 - t)
            conditional_power = 1.0 - _NORMAL.cdf(
                (final_critical_z - conditional_mean) / conditional_sd
            )

        conditional_power = max(0.0, min(1.0, conditional_power))
        is_futile = conditional_power < futility_threshold
        comparison = "below" if is_futile else "at or above"
        recommendation = (
            f"Conditional power {conditional_power:.1%} is {comparison} the "
            f"non-binding futility threshold {futility_threshold:.1%}."
        )
        return FutilityAnalysis(
            is_futile=is_futile,
            conditional_power=round(conditional_power, 6),
            predicted_success_probability=round(conditional_power, 6),
            recommendation=recommendation,
            monitoring_threshold=futility_threshold,
            interim_z=round(interim_z, 6),
            final_critical_z=round(final_critical_z, 6),
        )

    def assess(
        self,
        observed_effect: float,
        information_fraction: float = 0.5,
        variance: float = 1.0,
        futility_threshold: float = 0.10,
        planned_information: float = 100.0,
        alpha: float = 0.05,
        assumed_effect: Optional[float] = None,
    ) -> FutilityAnalysis:
        """Compatibility interface using an effect estimate and information scale."""
        if not math.isfinite(observed_effect):
            raise ValueError("observed_effect must be finite")
        if not 0.0 < information_fraction <= 1.0:
            raise ValueError("information_fraction must be in (0, 1]")
        if not math.isfinite(variance) or variance <= 0.0:
            raise ValueError("variance must be positive and finite")
        if not math.isfinite(planned_information) or planned_information <= 0.0:
            raise ValueError("planned_information must be positive and finite")
        _validate_alpha(alpha)
        if assumed_effect is not None and not math.isfinite(assumed_effect):
            raise ValueError("assumed_effect must be finite")

        interim_information = planned_information * information_fraction
        interim_z = observed_effect * math.sqrt(interim_information / variance)
        future_effect = observed_effect if assumed_effect is None else assumed_effect
        assumed_final_mean_z = future_effect * math.sqrt(planned_information / variance)
        return self.assess_from_z(
            interim_z=interim_z,
            information_fraction=information_fraction,
            final_critical_z=_NORMAL.inv_cdf(1.0 - alpha / 2.0),
            assumed_final_mean_z=assumed_final_mean_z,
            futility_threshold=futility_threshold,
        )


class SampleSizeReestimator:
    """Simple variance- and effect-ratio sample-size re-estimation utilities."""

    def blinded_reestimate(
        self,
        initial_n: int,
        current_variance: float,
        expected_variance: float,
    ) -> SampleSizeReestimate:
        if initial_n <= 0:
            raise ValueError(f"initial_n must be positive, got {initial_n}")
        if not math.isfinite(current_variance) or current_variance <= 0.0:
            raise ValueError("current_variance must be positive and finite")
        if not math.isfinite(expected_variance) or expected_variance <= 0.0:
            raise ValueError("expected_variance must be positive and finite")

        inflation = current_variance / expected_variance
        revised = max(initial_n, int(math.ceil(initial_n * inflation)))
        return SampleSizeReestimate(
            original_n=initial_n,
            revised_n=revised,
            inflation_factor=round(inflation, 6),
            reason="Variance-ratio blinded re-estimation",
        )

    def unblinded_reestimate(
        self,
        initial_n: int,
        observed_effect: float,
        expected_effect: float,
    ) -> SampleSizeReestimate:
        if initial_n <= 0:
            raise ValueError(f"initial_n must be positive, got {initial_n}")
        if not math.isfinite(observed_effect) or observed_effect <= 0.0:
            raise ValueError("observed_effect must be positive and finite")
        if not math.isfinite(expected_effect) or expected_effect <= 0.0:
            raise ValueError("expected_effect must be positive and finite")

        inflation = (expected_effect / observed_effect) ** 2
        revised = max(initial_n, int(math.ceil(initial_n * inflation)))
        return SampleSizeReestimate(
            original_n=initial_n,
            revised_n=revised,
            inflation_factor=round(inflation, 6),
            reason="Squared effect-ratio unblinded re-estimation",
        )

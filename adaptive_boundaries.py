#!/usr/bin/env python3
"""
Adaptive Clinical Trial: O'Brien-Fleming & Pocock Interim Stopping Boundaries,
Spending Functions, Sample Size Re-estimation, and Futility Assessment.
"""
from dataclasses import dataclass, field
from typing import List, Dict, Optional
import math
from enum import Enum


class BoundaryMethod(str, Enum):
    OBRIEN_FLEMING = "O'Brien-Fleming"
    POCOCK = "Pocock"
    LAN_DEMETS = "Lan-DeMets"


@dataclass
class InterimLook:
    look_number: int
    information_fraction: float  # 0-1
    sample_size_n: int
    observed_effect: float
    observed_se: float
    alpha_spent: float


@dataclass
class StoppingBoundary:
    method: str
    look_number: int
    z_bound: float
    alpha_spent: float
    cumulative_alpha: float
    stopping_rule: str  # efficacy, futility, none


@dataclass
class SpendingFunctionResult:
    function_name: str
    alpha_at_looks: List[float]
    total_alpha_spent: float


@dataclass
class FutilityAnalysis:
    is_futile: bool
    conditional_power: float
    predicted_success_probability: float
    recommendation: str
    monitoring_threshold: float


@dataclass
class SampleSizeReestimate:
    original_n: int
    revised_n: int
    inflation_factor: float
    reason: str


class OBFlemingBoundary:
    """O'Brien-Fleming alpha spending function for group sequential designs."""

    def __init__(self, total_alpha: float = 0.05, num_looks: int = 5):
        if not (0 < total_alpha < 1):
            raise ValueError(f"total_alpha must be between 0 and 1, got {total_alpha}")
        if num_looks < 1:
            raise ValueError(f"num_looks must be >= 1, got {num_looks}")
        self.total_alpha = total_alpha
        self.num_looks = num_looks

    def compute_boundaries(self) -> List[StoppingBoundary]:
        boundaries = []
        cumulative = 0.0
        for i in range(1, self.num_looks + 1):
            information_fraction = i / self.num_looks
            z_bound = self._obf_z(information_fraction)
            alpha_at_look = self._alpha_spent_at_look(i)
            cumulative += alpha_at_look

            stopping = "none"
            if i == self.num_looks:
                stopping = "efficacy" if z_bound < 2.0 else "none"
            elif alpha_at_look > 0.001:
                stopping = "efficacy"

            boundaries.append(StoppingBoundary(
                method="O'Brien-Fleming",
                look_number=i,
                z_bound=round(z_bound, 4),
                alpha_spent=round(alpha_at_look, 6),
                cumulative_alpha=round(cumulative, 6),
                stopping_rule=stopping,
            ))
        return boundaries

    def _obf_z(self, t: float) -> float:
        if t <= 0:
            return float('inf')
        z0 = 2.7965  # standard OBF critical value at final
        return z0 / math.sqrt(t)

    def _alpha_spent_at_look(self, look: int) -> float:
        t_prev = (look - 1) / self.num_looks
        t_curr = look / self.num_looks
        alpha_prev = 2 * (1 - self._norm_cdf(self._obf_z(max(t_prev, 0.001))))
        alpha_curr = 2 * (1 - self._norm_cdf(self._obf_z(t_curr)))
        return max(0, self.total_alpha * (alpha_curr / alpha_prev)) if alpha_prev > 0 else 0

    @staticmethod
    def _norm_cdf(x: float) -> float:
        return 0.5 * (1 + math.erf(x / math.sqrt(2)))


class PocockBoundary:
    """Pocock alpha spending function with constant boundary values."""

    def __init__(self, total_alpha: float = 0.05, num_looks: int = 5):
        if not (0 < total_alpha < 1):
            raise ValueError(f"total_alpha must be between 0 and 1, got {total_alpha}")
        if num_looks < 1:
            raise ValueError(f"num_looks must be >= 1, got {num_looks}")
        self.total_alpha = total_alpha
        self.num_looks = num_looks

    def compute_boundaries(self) -> List[StoppingBoundary]:
        boundaries = []
        z0 = 2.414  # Pocock critical value (approximate for 5 looks)
        cumulative = 0.0

        for i in range(1, self.num_looks + 1):
            alpha_at_look = self.total_alpha / self.num_looks
            cumulative += alpha_at_look
            stopping = "efficacy" if i < self.num_looks else "efficacy"

            boundaries.append(StoppingBoundary(
                method="Pocock",
                look_number=i,
                z_bound=round(z0, 4),
                alpha_spent=round(alpha_at_look, 6),
                cumulative_alpha=round(cumulative, 6),
                stopping_rule=stopping,
            ))
        return boundaries


class LanDeMetsSpending:
    """Lan-DeMets alpha spending functions (alpha-spending approach)."""

    @staticmethod
    def obf_spending(t: float, alpha: float = 0.05) -> float:
        if t <= 0:
            return 0.0
        return alpha * (2 - 2 * math.erf(math.sqrt(1.0 / t) * 1.96 / math.sqrt(2)))

    @staticmethod
    def pocock_spending(t: float, alpha: float = 0.05) -> float:
        if t <= 0:
            return 0.0
        return alpha * math.log(1 + (math.e - 1) * t)

    def compute_looks(self, information_fractions: List[float], alpha: float = 0.05,
                       method: str = "OBF") -> SpendingFunctionResult:
        if not (0 < alpha < 1):
            raise ValueError(f"alpha must be between 0 and 1, got {alpha}")
        if not information_fractions:
            raise ValueError("information_fractions must not be empty")
        if any(not (0 < t <= 1) for t in information_fractions):
            raise ValueError("information_fractions must be in (0, 1]")
        spending_fn = self.obf_spending if method == "OBF" else self.pocock_spending
        alpha_at_looks = []
        prev = 0.0
        for t in information_fractions:
            spent = spending_fn(t, alpha)
            alpha_at_looks.append(round(spent - prev, 6))
            prev = spent

        return SpendingFunctionResult(
            function_name=f"Lan-DeMets {method}",
            alpha_at_looks=alpha_at_looks,
            total_alpha_spent=round(prev, 6),
        )


class FutilityAssessor:
    """Conditional power and futility monitoring for adaptive trials."""

    def assess(self, observed_effect: float, target_power: float = 0.8,
               information_fraction: float = 0.5, variance: float = 1.0,
               futility_threshold: float = 0.1) -> FutilityAnalysis:
        if not (0 <= information_fraction <= 1):
            raise ValueError(f"information_fraction must be in [0, 1], got {information_fraction}")
        if not (0 <= target_power <= 1):
            raise ValueError(f"target_power must be in [0, 1], got {target_power}")
        if not (0 <= futility_threshold <= 1):
            raise ValueError(f"futility_threshold must be in [0, 1], got {futility_threshold}")
        remaining_info = 1.0 - information_fraction
        if remaining_info <= 0 or variance <= 0:
            return FutilityAnalysis(False, 0, 0, "Insufficient information", futility_threshold)

        se_remaining = math.sqrt(variance / (remaining_info * 100))
        z_remaining = observed_effect / se_remaining if se_remaining > 0 else 0
        conditional_power = self._norm_cdf(z_remaining - 1.96) if z_remaining > 0 else 0

        is_futile = conditional_power < futility_threshold

        if is_futile:
            rec = f"Conditional power {conditional_power:.1%} < threshold {futility_threshold:.1%}. Recommend futility stop."
        elif conditional_power < 0.3:
            rec = f"Conditional power {conditional_power:.1%} is low. Continue monitoring closely."
        else:
            rec = f"Conditional power {conditional_power:.1%} is acceptable. Continue trial."

        return FutilityAnalysis(
            is_futile=is_futile,
            conditional_power=round(conditional_power, 4),
            predicted_success_probability=round(conditional_power, 4),
            recommendation=rec,
            monitoring_threshold=futility_threshold,
        )

    @staticmethod
    def _norm_cdf(x: float) -> float:
        return 0.5 * (1 + math.erf(x / math.sqrt(2)))


class SampleSizeReestimator:
    """Blinded and unblinded sample size re-estimation for adaptive designs."""

    def blinded_reestimate(self, initial_n: int, current_variance: float,
                            expected_variance: float, target_power: float = 0.8) -> SampleSizeReestimate:
        if initial_n <= 0:
            raise ValueError(f"initial_n must be positive, got {initial_n}")
        if current_variance <= 0 or expected_variance <= 0:
            return SampleSizeReestimate(initial_n, initial_n, 1.0, "Invalid variance")

        inflation = current_variance / expected_variance
        revised = int(math.ceil(initial_n * inflation))
        revised = max(revised, initial_n)  # never reduce below original

        return SampleSizeReestimate(
            original_n=initial_n,
            revised_n=revised,
            inflation_factor=round(inflation, 3),
            reason=f"Blinded variance ratio: {inflation:.2f}",
        )

    def unblinded_reestimate(self, initial_n: int, observed_effect: float,
                              expected_effect: float, alpha: float = 0.05,
                              target_power: float = 0.8) -> SampleSizeReestimate:
        if initial_n <= 0:
            raise ValueError(f"initial_n must be positive, got {initial_n}")
        if observed_effect <= 0 or expected_effect <= 0:
            return SampleSizeReestimate(initial_n, initial_n, 1.0, "Invalid effect sizes")

        ratio = (expected_effect / observed_effect) ** 2
        revised = int(math.ceil(initial_n * ratio))
        revised = max(revised, initial_n)

        return SampleSizeReestimate(
            original_n=initial_n,
            revised_n=revised,
            inflation_factor=round(ratio, 3),
            reason=f"Effect size ratio: {ratio:.2f}",
        )

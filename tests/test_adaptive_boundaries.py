import pytest

from adaptive_boundaries import (
    FutilityAssessor,
    LanDeMetsSpending,
    OBFlemingBoundary,
    PocockBoundary,
    SampleSizeReestimator,
)


def test_obrien_fleming_five_look_reference_values():
    rows = OBFlemingBoundary(total_alpha=0.05, num_looks=5).compute_boundaries()
    assert len(rows) == 5
    assert rows[0].z_bound == pytest.approx(4.56, abs=0.03)
    assert rows[-1].z_bound == pytest.approx(2.04, abs=0.02)
    assert rows[-1].cumulative_alpha == pytest.approx(0.05, abs=2e-4)
    assert all(a.z_bound > b.z_bound for a, b in zip(rows, rows[1:]))


def test_pocock_five_look_reference_value():
    rows = PocockBoundary(total_alpha=0.05, num_looks=5).compute_boundaries()
    assert len(rows) == 5
    assert rows[0].z_bound == pytest.approx(2.413, abs=0.02)
    assert all(row.z_bound == rows[0].z_bound for row in rows)
    assert rows[-1].cumulative_alpha == pytest.approx(0.05, abs=2e-4)


def test_single_look_matches_standard_normal_critical_value():
    row = OBFlemingBoundary(total_alpha=0.05, num_looks=1).compute_boundaries()[0]
    assert row.z_bound == pytest.approx(1.96, abs=0.002)
    assert row.cumulative_alpha == pytest.approx(0.05, abs=2e-4)


def test_custom_information_fractions_are_supported():
    rows = OBFlemingBoundary(information_fractions=[0.25, 0.5, 1.0]).compute_boundaries()
    assert [row.information_fraction for row in rows] == [0.25, 0.5, 1.0]
    assert rows[-1].cumulative_alpha == pytest.approx(0.05, abs=2e-4)


def test_boundary_validation():
    with pytest.raises(ValueError):
        OBFlemingBoundary(total_alpha=0)
    with pytest.raises(ValueError):
        PocockBoundary(num_looks=0)
    with pytest.raises(ValueError):
        OBFlemingBoundary(information_fractions=[0.5, 0.4, 1.0])
    with pytest.raises(ValueError):
        OBFlemingBoundary(information_fractions=[0.25, 0.5, 0.75])


def test_lan_demets_obf_spending_reaches_alpha_at_final_look():
    result = LanDeMetsSpending().compute_looks([0.2, 0.4, 0.6, 0.8, 1.0], 0.05, "OBF")
    assert result.total_alpha_spent == pytest.approx(0.05, abs=1e-8)
    assert result.cumulative_alpha[-1] == pytest.approx(0.05, abs=1e-8)
    assert all(x >= 0 for x in result.alpha_at_looks)


def test_lan_demets_pocock_spending_reaches_alpha_at_final_look():
    result = LanDeMetsSpending().compute_looks([0.25, 0.5, 0.75, 1.0], 0.025, "pocock")
    assert result.total_alpha_spent == pytest.approx(0.025, abs=1e-8)


def test_lan_demets_rejects_unknown_method_and_unsorted_fractions():
    calc = LanDeMetsSpending()
    with pytest.raises(ValueError):
        calc.compute_looks([0.5, 1.0], method="other")
    with pytest.raises(ValueError):
        calc.compute_looks([0.5, 0.5, 1.0])


def test_conditional_power_observed_trend():
    result = FutilityAssessor.assess_from_z(2.0, 0.5, futility_threshold=0.1)
    assert result.conditional_power == pytest.approx(0.8903, abs=0.003)
    assert result.is_futile is False


def test_conditional_power_final_look_is_deterministic():
    assert FutilityAssessor.assess_from_z(2.1, 1.0).conditional_power == 1.0
    assert FutilityAssessor.assess_from_z(1.0, 1.0).conditional_power == 0.0


def test_conditional_power_validation():
    with pytest.raises(ValueError):
        FutilityAssessor.assess_from_z(1.0, 0.0)
    with pytest.raises(ValueError):
        FutilityAssessor.assess_from_z(1.0, 0.5, futility_threshold=1.1)


def test_sample_size_reestimation():
    calc = SampleSizeReestimator()
    assert calc.blinded_reestimate(100, 1.5, 1.0).revised_n == 150
    assert calc.unblinded_reestimate(100, 0.4, 0.5).revised_n == 157


def test_sample_size_reestimation_never_reduces_n():
    calc = SampleSizeReestimator()
    assert calc.blinded_reestimate(100, 0.8, 1.0).revised_n == 100
    assert calc.unblinded_reestimate(100, 0.6, 0.5).revised_n == 100

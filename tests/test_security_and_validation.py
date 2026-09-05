"""
Security and input validation tests for adaptive-trial-interim-stopping.
"""
import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
import warnings
from adaptive_boundaries import (
    OBFlemingBoundary,
    PocockBoundary,
    LanDeMetsSpending,
    FutilityAssessor,
    SampleSizeReestimator,
)
from agents.base import PHIGuard, SecurityException, AuditTrail


class TestOBFlemingBoundaryValidation:
    def test_valid_construction(self):
        obf = OBFlemingBoundary(total_alpha=0.05, num_looks=5)
        assert obf.total_alpha == 0.05
        assert obf.num_looks == 5

    def test_invalid_alpha_zero(self):
        with pytest.raises(ValueError, match="total_alpha must be between 0 and 1"):
            OBFlemingBoundary(total_alpha=0.0, num_looks=5)

    def test_invalid_alpha_negative(self):
        with pytest.raises(ValueError, match="total_alpha must be between 0 and 1"):
            OBFlemingBoundary(total_alpha=-0.05, num_looks=5)

    def test_invalid_alpha_too_large(self):
        with pytest.raises(ValueError, match="total_alpha must be between 0 and 1"):
            OBFlemingBoundary(total_alpha=1.0, num_looks=5)

    def test_invalid_num_looks_zero(self):
        with pytest.raises(ValueError, match="num_looks must be >= 1"):
            OBFlemingBoundary(total_alpha=0.05, num_looks=0)

    def test_invalid_num_looks_negative(self):
        with pytest.raises(ValueError, match="num_looks must be >= 1"):
            OBFlemingBoundary(total_alpha=0.05, num_looks=-3)

    def test_single_look(self):
        obf = OBFlemingBoundary(total_alpha=0.05, num_looks=1)
        boundaries = obf.compute_boundaries()
        assert len(boundaries) == 1


class TestPocockBoundaryValidation:
    def test_valid_construction(self):
        poc = PocockBoundary(total_alpha=0.05, num_looks=5)
        assert poc.total_alpha == 0.05
        assert poc.num_looks == 5

    def test_invalid_alpha(self):
        with pytest.raises(ValueError, match="total_alpha must be between 0 and 1"):
            PocockBoundary(total_alpha=0.0, num_looks=5)

    def test_invalid_num_looks(self):
        with pytest.raises(ValueError, match="num_looks must be >= 1"):
            PocockBoundary(total_alpha=0.05, num_looks=0)


class TestLanDeMetsSpendingValidation:
    def test_valid_computation(self):
        ldm = LanDeMetsSpending()
        result = ldm.compute_looks([0.2, 0.4, 0.6, 0.8, 1.0])
        assert len(result.alpha_at_looks) == 5

    def test_invalid_alpha(self):
        ldm = LanDeMetsSpending()
        with pytest.raises(ValueError, match="alpha must be between 0 and 1"):
            ldm.compute_looks([0.5], alpha=0.0)

    def test_empty_fractions(self):
        ldm = LanDeMetsSpending()
        with pytest.raises(ValueError, match="information_fractions must not be empty"):
            ldm.compute_looks([])

    def test_invalid_fraction_zero(self):
        ldm = LanDeMetsSpending()
        with pytest.raises(ValueError, match="information_fractions must be in"):
            ldm.compute_looks([0.0, 0.5])

    def test_invalid_fraction_negative(self):
        ldm = LanDeMetsSpending()
        with pytest.raises(ValueError, match="information_fractions must be in"):
            ldm.compute_looks([-0.1, 0.5])

    def test_invalid_fraction_too_large(self):
        ldm = LanDeMetsSpending()
        with pytest.raises(ValueError, match="information_fractions must be in"):
            ldm.compute_looks([0.5, 1.1])


class TestFutilityAssessorValidation:
    def test_valid_assessment(self):
        fa = FutilityAssessor()
        result = fa.assess(observed_effect=0.5, information_fraction=0.5)
        assert result.conditional_power >= 0

    def test_invalid_information_fraction_negative(self):
        fa = FutilityAssessor()
        with pytest.raises(ValueError, match="information_fraction must be in"):
            fa.assess(observed_effect=0.5, information_fraction=-0.1)

    def test_invalid_information_fraction_too_large(self):
        fa = FutilityAssessor()
        with pytest.raises(ValueError, match="information_fraction must be in"):
            fa.assess(observed_effect=0.5, information_fraction=1.5)

    def test_invalid_target_power(self):
        fa = FutilityAssessor()
        with pytest.raises(ValueError, match="target_power must be in"):
            fa.assess(observed_effect=0.5, target_power=1.5)

    def test_invalid_futility_threshold(self):
        fa = FutilityAssessor()
        with pytest.raises(ValueError, match="futility_threshold must be in"):
            fa.assess(observed_effect=0.5, futility_threshold=-0.1)


class TestSampleSizeReestimatorValidation:
    def test_blinded_valid(self):
        ssr = SampleSizeReestimator()
        result = ssr.blinded_reestimate(initial_n=100, current_variance=1.5, expected_variance=1.0)
        assert result.revised_n == 150

    def test_blinded_invalid_initial_n(self):
        ssr = SampleSizeReestimator()
        with pytest.raises(ValueError, match="initial_n must be positive"):
            ssr.blinded_reestimate(initial_n=0, current_variance=1.5, expected_variance=1.0)

    def test_blinded_negative_initial_n(self):
        ssr = SampleSizeReestimator()
        with pytest.raises(ValueError, match="initial_n must be positive"):
            ssr.blinded_reestimate(initial_n=-10, current_variance=1.5, expected_variance=1.0)

    def test_unblinded_invalid_initial_n(self):
        ssr = SampleSizeReestimator()
        with pytest.raises(ValueError, match="initial_n must be positive"):
            ssr.unblinded_reestimate(initial_n=0, observed_effect=0.5, expected_effect=0.8)


class TestPHIGuard:
    def test_clean_text_passes(self):
        PHIGuard.assert_no_phi("Normal clinical trial data with no PHI")

    def test_mrn_detected(self):
        with pytest.raises(SecurityException, match="PHI Outbound Guard Violation"):
            PHIGuard.assert_no_phi("Patient MRN-12345678")

    def test_ssn_detected(self):
        with pytest.raises(SecurityException, match="PHI Outbound Guard Violation"):
            PHIGuard.assert_no_phi("SSN: 123-45-6789")

    def test_email_detected(self):
        with pytest.raises(SecurityException, match="PHI Outbound Guard Violation"):
            PHIGuard.assert_no_phi("Contact patient@example.com")

    def test_phone_detected(self):
        with pytest.raises(SecurityException, match="PHI Outbound Guard Violation"):
            PHIGuard.assert_no_phi("Call 555-123-4567")

    def test_empty_text_passes(self):
        PHIGuard.assert_no_phi("")

    def test_none_text_passes(self):
        PHIGuard.assert_no_phi(None)

    def test_redact_phi(self):
        redacted = PHIGuard.redact_phi("Patient MRN-12345678 and email test@example.com")
        assert "MRN-12345678" not in redacted
        assert "test@example.com" not in redacted
        assert "[REDACTED_IDENTIFIER]" in redacted


class TestAuditTrail:
    def test_audit_trail_with_explicit_key(self):
        trail = AuditTrail(secret_key="test-key-123")
        entry = trail.log("actor1", "tier1", "TEST_EVENT", {"data": "value"})
        assert "current_hash" in entry
        assert trail.verify_integrity() is True

    def test_audit_trail_generates_key_without_env(self):
        # Ensure no env var is set
        old_key = os.environ.pop("AUDIT_SECRET_KEY", None)
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                trail = AuditTrail()
                # Should have issued a warning about missing key
                assert any("AUDIT_SECRET_KEY" in str(warning.message) for warning in w)
        finally:
            if old_key is not None:
                os.environ["AUDIT_SECRET_KEY"] = old_key

    def test_audit_trail_chain_integrity(self):
        trail = AuditTrail(secret_key="chain-test-key")
        trail.log("actor1", "tier1", "EVENT_1", {"step": 1})
        trail.log("actor2", "tier2", "EVENT_2", {"step": 2})
        trail.log("actor3", "tier3", "EVENT_3", {"step": 3})
        assert trail.verify_integrity() is True
        assert len(trail.get_trail()) == 3

    def test_audit_trail_tamper_detection(self):
        trail = AuditTrail(secret_key="tamper-test-key")
        trail.log("actor1", "tier1", "EVENT_1", {"step": 1})
        trail.log("actor2", "tier2", "EVENT_2", {"step": 2})
        # Tamper with the first entry
        trail.logs[0]["current_hash"] = "TAMPERED_HASH"
        assert trail.verify_integrity() is False

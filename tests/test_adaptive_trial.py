import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from adaptive_trial.models import FrontierPayload, ExecutionStatus
from adaptive_trial.engine import FrontierDomainEngine
from adaptive_trial.agents import AlphaSpendingBoundaryAgent, ConditionalPowerEstimatorAgent, FutilityStoppingEvaluatorAgent, AdaptiveTrialCoordinator
from adaptive_trial.cli import main


def test_sub_agents():
    a1 = AlphaSpendingBoundaryAgent()
    p1 = FrontierPayload("T1", "KEY-01", primary_metric=35.0, secondary_metric=4.0, status_descriptor="NOMINAL")
    alerts1 = a1.audit(p1)
    assert len(alerts1) == 1
    assert alerts1[0].status == ExecutionStatus.ELEVATED_RISK

    a2 = ConditionalPowerEstimatorAgent()
    p2 = FrontierPayload("T2", "KEY-02", primary_metric=10.0, secondary_metric=15.0, status_descriptor="NOMINAL", is_critical_flag=True)
    alerts2 = a2.audit(p2)
    assert len(alerts2) == 1
    assert alerts2[0].status == ExecutionStatus.CRITICAL_INTERVENTION

    a3 = FutilityStoppingEvaluatorAgent()
    p3 = FrontierPayload("T3", "KEY-03", primary_metric=10.0, secondary_metric=4.0, status_descriptor="DISCORDANT_ANOMALY")
    alerts3 = a3.audit(p3)
    assert len(alerts3) == 1


def test_coordinator():
    coord = AdaptiveTrialCoordinator()
    p_nominal = FrontierPayload("T4", "KEY-04", primary_metric=12.0, secondary_metric=4.0, status_descriptor="NOMINAL")
    dossier = coord.process(p_nominal)
    assert dossier["overall_status"] == ExecutionStatus.NOMINAL.value
    assert dossier["total_alerts"] == 0

    ans = coord.query_supervisory_chat("What standard is applied?")
    assert "FDA Adaptive Clinical Trial Guidelines" in ans or "specifications" in ans


def test_cli():
    assert main(["audit", "--task-id", "CLI-01"]) == 0
    assert main(["chat", "What", "is", "the", "system", "status?"]) == 0
    assert main(["audit", "--task-id", "CLI-02", "--json"]) == 0


def test_cli_batch_sample_csv(tmp_path):
    sample_csv = Path(__file__).resolve().parent.parent / "sample.csv"
    out_file = tmp_path / "out.csv"
    res = main(["batch", "-i", str(sample_csv), "-o", str(out_file)])
    assert res == 0
    assert out_file.exists()


def test_adaptive_boundaries_direct():
    from adaptive_boundaries import (
        OBFlemingBoundary,
        PocockBoundary,
        LanDeMetsSpending,
        FutilityAssessor,
        SampleSizeReestimator,
    )
    obf = OBFlemingBoundary(total_alpha=0.05, num_looks=5).compute_boundaries()
    assert len(obf) == 5
    assert obf[0].z_bound > obf[-1].z_bound

    poc = PocockBoundary(total_alpha=0.05, num_looks=5).compute_boundaries()
    assert len(poc) == 5
    assert poc[0].z_bound == poc[-1].z_bound

    ldm = LanDeMetsSpending().compute_looks([0.2, 0.4, 0.6, 0.8, 1.0])
    assert len(ldm.alpha_at_looks) == 5
    assert ldm.total_alpha_spent > 0

    fut = FutilityAssessor().assess(observed_effect=0.5, information_fraction=0.5)
    assert not fut.is_futile
    assert fut.conditional_power > 0.5

    ssr = SampleSizeReestimator().blinded_reestimate(initial_n=100, current_variance=1.5, expected_variance=1.0)
    assert ssr.revised_n == 150


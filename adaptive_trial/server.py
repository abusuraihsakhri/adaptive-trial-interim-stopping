"""Optional FastAPI service exposing the statistical calculators."""
from __future__ import annotations

from dataclasses import asdict
from typing import Literal

from adaptive_boundaries import FutilityAssessor, LanDeMetsSpending, OBFlemingBoundary, PocockBoundary, SampleSizeReestimator


def create_app():
    try:
        from fastapi import FastAPI
        from pydantic import BaseModel, Field
    except ImportError as exc:
        raise RuntimeError("Install server dependencies with: pip install -e '.[server]'") from exc

    app = FastAPI(
        title="Adaptive Trial Interim Stopping",
        description="Group-sequential efficacy boundaries, alpha spending, conditional power, and sample-size re-estimation.",
        version="3.0.0",
    )

    class BoundaryRequest(BaseModel):
        method: Literal["obf", "pocock"] = "obf"
        alpha: float = Field(0.05, gt=0, lt=1)
        looks: int = Field(5, ge=1, le=20)

    class SpendingRequest(BaseModel):
        method: Literal["obf", "pocock"] = "obf"
        alpha: float = Field(0.05, gt=0, lt=1)
        information_fractions: list[float]

    class FutilityRequest(BaseModel):
        interim_z: float
        information_fraction: float = Field(..., gt=0, le=1)
        final_critical_z: float = Field(1.959963984540054, gt=0)
        assumed_final_mean_z: float | None = None
        futility_threshold: float = Field(0.10, ge=0, le=1)

    class BlindedSSRRequest(BaseModel):
        initial_n: int = Field(..., gt=0)
        current_variance: float = Field(..., gt=0)
        expected_variance: float = Field(..., gt=0)

    @app.get("/health")
    def health():
        return {"status": "ok", "version": "3.0.0"}

    @app.post("/api/boundaries")
    def boundaries(req: BoundaryRequest):
        cls = OBFlemingBoundary if req.method == "obf" else PocockBoundary
        return [asdict(x) for x in cls(req.alpha, req.looks).compute_boundaries()]

    @app.post("/api/spending")
    def spending(req: SpendingRequest):
        return asdict(LanDeMetsSpending().compute_looks(req.information_fractions, req.alpha, req.method))

    @app.post("/api/futility")
    def futility(req: FutilityRequest):
        return asdict(FutilityAssessor.assess_from_z(
            req.interim_z,
            req.information_fraction,
            req.final_critical_z,
            req.assumed_final_mean_z,
            req.futility_threshold,
        ))

    @app.post("/api/ssr/blinded")
    def blinded_ssr(req: BlindedSSRRequest):
        return asdict(SampleSizeReestimator().blinded_reestimate(
            req.initial_n,
            req.current_variance,
            req.expected_variance,
        ))

    return app

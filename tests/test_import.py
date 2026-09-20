import adaptive_trial


def test_public_imports():
    assert adaptive_trial.__version__ == "3.0.0"
    assert adaptive_trial.OBFlemingBoundary().compute_boundaries()

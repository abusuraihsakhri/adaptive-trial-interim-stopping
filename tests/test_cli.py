import csv
import json

from adaptive_trial.cli import main


def test_boundaries_json(capsys):
    assert main(["boundaries", "--method", "obf", "--looks", "3", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert len(payload) == 3
    assert payload[-1]["cumulative_alpha"] == 0.05


def test_spending_json(capsys):
    assert main(["spending", "--method", "pocock", "--fractions", "0.25,0.5,0.75,1", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["total_alpha_spent"] == 0.05


def test_futility_json(capsys):
    assert main(["futility", "--interim-z", "2", "--information-fraction", "0.5", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert 0 <= payload["conditional_power"] <= 1


def test_batch(tmp_path):
    src = tmp_path / "input.csv"
    out = tmp_path / "output.csv"
    src.write_text("method,alpha,looks\nobf,0.05,3\npocock,0.05,2\n", encoding="utf-8")
    assert main(["batch", "-i", str(src), "-o", str(out)]) == 0
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert len(rows) == 5
    assert {row["method"] for row in rows} == {"O'Brien-Fleming", "Pocock"}

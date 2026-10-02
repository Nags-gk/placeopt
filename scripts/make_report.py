"""Bundle results/*.json into docs/data/results.json for the dashboard."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
out = {k: json.loads((ROOT / f"results/{k}.json").read_text())
       for k in ("flows", "exact", "viz") if (ROOT / f"results/{k}.json").exists()}
(ROOT / "docs/data").mkdir(parents=True, exist_ok=True)
(ROOT / "docs/data/results.json").write_text(json.dumps(out, separators=(",", ":")))
print("wrote docs/data/results.json", {k: len(v) for k, v in out.items()})

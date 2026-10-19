"""Export everything the dashboard shows into one JSON file."""
import json, sys, warnings
from pathlib import Path
from Bio import SeqIO
sys.path.insert(0, str(Path(__file__).parent))
import codonrescore as cr, evaluate as ev
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
bm = json.load(open(ROOT/"results/benchmark.json")); an = json.load(open(ROOT/"results/analysis.json"))
out = {"benchmark": bm, "analysis": an, "records": {}}
for fn, label in ev.RECORDS:
    rec = SeqIO.read(ROOT/"data"/fn, "gb"); seq = str(rec.seq).upper()
    truth, neutral = ev.load_truth(rec)
    M, model, tau, info = ev.methods_for(seq, with_single=len(seq) >= 20000)
    M.pop("_aux")
    tk = {ev.stop_key(t) for t in truth}
    R = {"length": len(seq), "truth": [[t["begin"], t["end"], t["strand"], t["product"]] for t in truth],
         "neutral": [[n["begin"], n["end"], n["strand"]] for n in neutral], "methods": {}}
    for name, pred in M.items():
        keys = {ev.stop_key(p) for p in pred}
        calls = []
        for p in pred:
            if ev.stop_key(p) in tk: c = "TP"
            elif any(n["strand"] == p["strand"] and cr.overlap(p, n) >= 0.5*p["length"] for n in neutral): c = "N"
            else: c = "FP"
            calls.append([p["begin"], p["end"], p["strand"], c])
        R["methods"][name] = {"calls": calls, "missed": [i for i, t in enumerate(truth) if ev.stop_key(t) not in keys]}
    out["records"][label] = R
json.dump(out, open(ROOT/"dashboard/data.json", "w"), separators=(",", ":"))
print(Path(ROOT/"dashboard/data.json").stat().st_size)

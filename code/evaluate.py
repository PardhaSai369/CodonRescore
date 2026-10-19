"""Run all methods on the annotated records and score against GenBank CDS annotation."""
import gzip, json, sys, warnings
from pathlib import Path
import numpy as np
import pyrodigal
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
import codonrescore as cr

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
RECORDS = [("NC_000932.gb", "A. thaliana chloroplast"), ("NC_005816.gb", "Y. pestis plasmid pPCP1"),
           ("NC_001422.gbk", "phiX174 phage")]


def load_truth(rec):
    truth, neutral = [], []
    for f in rec.features:
        if f.type != "CDS":
            continue
        parts = f.location.parts
        pseudo = "pseudo" in f.qualifiers or "pseudogene" in f.qualifiers
        fuzzy = any(not isinstance(p.start, int.__mro__[0]) and False for p in parts)
        if len(parts) > 1 or pseudo:
            for p in parts:
                neutral.append(dict(begin=int(p.start) + 1, end=int(p.end), strand=int(p.strand)))
            continue
        p = parts[0]
        truth.append(dict(begin=int(p.start) + 1, end=int(p.end), strand=int(p.strand),
                          product=f.qualifiers.get("product", [""])[0]))
    for t in truth:
        t["length"] = t["end"] - t["begin"] + 1
    return truth, neutral


def stop_key(g):
    return (g["strand"], g["end"] if g["strand"] == 1 else g["begin"])


def pyro_calls(seq, meta):
    gf = pyrodigal.GeneFinder(meta=meta)
    if not meta:
        gf.train(seq.encode())
    out = []
    for g in gf.find_genes(seq.encode()):
        out.append(dict(begin=g.begin, end=g.end, strand=g.strand, length=g.end - g.begin + 1,
                        pscore=float(g.score), start=g.start_type))
    return out


def score(pred, truth, neutral):
    tkeys = {stop_key(t): t for t in truth}
    tp_stop = tp_exact = 0
    matched = set()
    fp = neu = 0
    for p in pred:
        k = stop_key(p)
        if k in tkeys:
            tp_stop += 1
            matched.add(k)
            t = tkeys[k]
            if t["begin"] == p["begin"] and t["end"] == p["end"]:
                tp_exact += 1
        else:
            if any(n["strand"] == p["strand"] and cr.overlap(p, n) >= 0.5 * p["length"] for n in neutral):
                neu += 1
            else:
                fp += 1
    n_eval = tp_stop + fp
    rec = len(matched) / len(truth)
    prec = tp_stop / n_eval if n_eval else 0.0
    f1 = 2 * rec * prec / (rec + prec) if rec + prec else 0.0
    exact_rec = tp_exact / len(truth)
    return dict(n_pred=len(pred), TP=tp_stop, FP=fp, neutral=neu, FN=len(truth) - len(matched), recall=rec,
                precision=prec, f1=f1, exact_recall=exact_rec), matched


def methods_for(seq, with_single=True):
    orfs = cr.find_orfs(seq, 90)
    M = {}
    M["LongORF>=300"] = [o for o in cr.find_orfs(seq, 300)]
    if with_single:
        M["Pyrodigal-single"] = pyro_calls(seq, False)
    pm = pyro_calls(seq, True)
    M["Pyrodigal-meta"] = pm
    try:
        model, tau, info = cr.train_codonrescore(seq, orfs)
    except ValueError as e:
        M["_aux"] = {}
        return M, None, None, dict(n_seeds=0, n_shadow=0, youden=float("nan"), seed_scores=np.array([]), shadow_scores=np.array([]), failed=str(e))
    M["CodonLLR-only"] = cr.call_codon_only(seq, model, tau, orfs)
    for base in (["Pyrodigal-single"] if with_single else []) + ["Pyrodigal-meta"]:
        final, pruned, rescued = cr.call_hybrid(seq, M[base], model, tau, orfs)
        M["CodonRescore(" + base.split("-")[1] + ")"] = final
        M.setdefault("_aux", {})[base] = dict(pruned=pruned, rescued=rescued)
    return M, model, tau, info


def main():
    results, per_gene = {}, {}
    allpred_matched = {}
    aux_all = {}
    for fn, label in RECORDS:
        rec = SeqIO.read(ROOT / "data" / fn, "gb")
        seq = str(rec.seq).upper()
        truth, neutral = load_truth(rec)
        M, model, tau, info = methods_for(seq, with_single=len(seq) >= 20000)
        aux = M.pop("_aux")
        results[label] = dict(length=len(seq), n_truth=len(truth), n_neutral_features=len(neutral),
                              tau=(tau if tau is not None else float("nan")), n_seeds=info["n_seeds"], trained=("failed" not in info), n_shadow=info["n_shadow"], youden=info["youden"], methods={})
        per_gene[label] = {}
        for name, pred in M.items():
            s, matched = score(pred, truth, neutral)
            results[label]["methods"][name] = s
            per_gene[label][name] = [stop_key(t) in matched for t in truth]
        results[label]["aux"] = {b: dict(pruned=len(a["pruned"]), rescued=len(a["rescued"])) for b, a in aux.items()}
        # which pruned/rescued were right?
        tk = {stop_key(t) for t in truth}
        results[label]["aux_detail"] = {b: dict(pruned_true=sum(stop_key(g) in tk for g in a["pruned"]),
                                                rescued_true=sum(stop_key(g) in tk for g in a["rescued"]))
                                        for b, a in aux.items()}
        np.savez(ROOT / "results" / f"scores_{fn.split('.')[0]}.npz", seed=info["seed_scores"], shadow=info["shadow_scores"], tau=(tau if tau is not None else np.nan))
    json.dump(results, open(ROOT / "results" / "benchmark.json", "w"), indent=1)
    json.dump(per_gene, open(ROOT / "results" / "per_gene_hits.json", "w"))
    print(json.dumps({k: {m: {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in s.items()}
                          for m, s in v["methods"].items()} for k, v in results.items()}, indent=1))
    for k, v in results.items():
        print(k, {a: b for a, b in v.items() if a not in ("methods",)})


if __name__ == "__main__":
    main()

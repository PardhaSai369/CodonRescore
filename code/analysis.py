"""Sensitivity, bootstrap CIs, C. diphtheriae agreement run, and figures."""
import gzip, json, sys, warnings
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
import codonrescore as cr
import evaluate as ev

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
rng = np.random.default_rng(2026)
C = {"base": "#6B7280", "pyro": "#0072B2", "ours": "#D55E00", "alt": "#009E73", "yel": "#E69F00"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 200,
                     "font.family": "DejaVu Sans"})

rec = SeqIO.read(ROOT / "data/NC_000932.gb", "gb")
seq = str(rec.seq).upper()
truth, neutral = ev.load_truth(rec)
orfs = cr.find_orfs(seq, 90)
model, tau, info = cr.train_codonrescore(seq, orfs)
single = ev.pyro_calls(seq, False)
meta = ev.pyro_calls(seq, True)
out = {}


# ---------- bootstrap CI for F1 (resample truth genes for recall, predictions for precision)
def boot_f1(pred, B=2000):
    tk = {ev.stop_key(t) for t in truth}
    hit_t = np.array([ev.stop_key(t) in {ev.stop_key(p) for p in pred} for t in truth])
    cls = []  # 1 TP, 0 FP, None neutral
    for p in pred:
        if ev.stop_key(p) in tk:
            cls.append(1)
        elif any(n["strand"] == p["strand"] and cr.overlap(p, n) >= 0.5 * p["length"] for n in neutral):
            continue
        else:
            cls.append(0)
    cls = np.array(cls)
    f = []
    for _ in range(B):
        r = hit_t[rng.integers(0, len(hit_t), len(hit_t))].mean()
        p_ = cls[rng.integers(0, len(cls), len(cls))].mean()
        f.append(2 * r * p_ / (r + p_) if r + p_ else 0)
    return np.percentile(f, [2.5, 97.5]).tolist()


sets = {"Pyrodigal-single": single, "Pyrodigal-meta": meta}
for b, P in (("single", single), ("meta", meta)):
    fin, pr, re_ = cr.call_hybrid(seq, P, model, tau, orfs)
    sets["CodonRescore(%s)" % b] = fin
sets["CodonLLR-only"] = cr.call_codon_only(seq, model, tau, orfs)
sets["LongORF>=300"] = cr.find_orfs(seq, 300)
out["bootstrap_f1_ci"] = {k: boot_f1(v) for k, v in sets.items()}

# ---------- sensitivity: tau offset and prune length (post-hoc, NOT used for tuning)
offs = np.round(np.arange(-0.3, 0.41, 0.1), 2)
sens = {"single": [], "meta": []}
for b, P in (("single", single), ("meta", meta)):
    for d in offs:
        fin, _, _ = cr.call_hybrid(seq, P, model, tau + d, orfs)
        s, _ = ev.score(fin, truth, neutral)
        sens[b].append(s)
out["sens_tau"] = {b: [dict(offset=float(d), f1=s["f1"], precision=s["precision"], recall=s["recall"]) for d, s in zip(offs, v)] for b, v in sens.items()}
pl = [150, 300, 450, 600]
out["sens_prune"] = {}
for b, P in (("single", single), ("meta", meta)):
    out["sens_prune"][b] = []
    for L in pl:
        fin, _, _ = cr.call_hybrid(seq, P, model, tau, orfs, prune_len=L)
        s, _ = ev.score(fin, truth, neutral)
        out["sens_prune"][b].append(dict(prune_len=L, f1=s["f1"], precision=s["precision"], recall=s["recall"]))

# ---------- length-stratified recall on truth genes
bins = [(0, 300), (300, 600), (600, 1200), (1200, 10**6)]
strat = {}
for name in ["Pyrodigal-single", "Pyrodigal-meta", "CodonRescore(single)", "CodonRescore(meta)"]:
    keys = {ev.stop_key(p) for p in sets[name]}
    strat[name] = []
    for lo, hi in bins:
        T = [t for t in truth if lo <= t["length"] < hi]
        strat[name].append([sum(ev.stop_key(t) in keys for t in T), len(T)])
out["length_strat"] = strat
out["length_bins"] = bins

# ---------- false-positive anatomy for Pyrodigal-single on chloroplast
tk = {ev.stop_key(t) for t in truth}
fps = []
for p in single:
    if ev.stop_key(p) in tk:
        continue
    if any(n["strand"] == p["strand"] and cr.overlap(p, n) >= 0.5 * p["length"] for n in neutral):
        continue
    fps.append(p)
annot = [f for f in rec.features if f.type in ("tRNA", "rRNA")]
out["single_fp"] = dict(n=len(fps), median_len=float(np.median([p["length"] for p in fps])),
                        in_rna=sum(any(cr.overlap(p, dict(begin=int(a.location.start) + 1, end=int(a.location.end))) > 0 for a in annot) for p in fps),
                        overlap_true_gene_otherframe=sum(any(cr.overlap(p, t) > 0.5 * p["length"] for t in truth) for p in fps))
out["single_fp_len_lt300"] = sum(p["length"] < 300 for p in fps)

# ---------- C. diphtheriae (no truth): agreement
with gzip.open(ROOT / "data/NCTC11397.fna.gz", "rt") as fh:
    cseq = "".join(l.strip() for l in fh if not l.startswith(">")).upper()
c_orfs = cr.find_orfs(cseq, 90)
cm, ct, ci = cr.train_codonrescore(cseq, c_orfs)
c_single = ev.pyro_calls(cseq, False)
c_meta = ev.pyro_calls(cseq, True)
res = {}
for b, P in (("single", c_single), ("meta", c_meta)):
    fin, pr, re_ = cr.call_hybrid(cseq, P, cm, ct, c_orfs)
    res[b] = dict(n_pyro=len(P), n_final=len(fin), pruned=len(pr), rescued=len(re_),
                  pruned_median_len=float(np.median([g["length"] for g in pr])) if pr else None,
                  rescued_median_len=float(np.median([g["length"] for g in re_])) if re_ else None)
ks, km = {ev.stop_key(g) for g in c_single}, {ev.stop_key(g) for g in c_meta}
res["agree_single_meta_stop"] = len(ks & km) / len(ks | km)
res["GC"] = (cseq.count("G") + cseq.count("C")) / len(cseq)
res["length"] = len(cseq)
res["n_seeds"] = ci["n_seeds"]; res["tau"] = ct; res["youden"] = ci["youden"]
res["coding_density_single"] = sum(g["length"] for g in c_single) / len(cseq)
out["cdiph"] = res
np.savez(ROOT / "results/scores_NCTC11397.npz", seed=ci["seed_scores"], shadow=ci["shadow_scores"], tau=ct)

# chloroplast descriptors
out["chloro"] = dict(length=len(seq), GC=(seq.count("G") + seq.count("C")) / len(seq), n_orfs90=len(orfs), n_seeds=info["n_seeds"],
                     n_shadow=info["n_shadow"], tau=tau, youden=info["youden"], seed_min_len=info["seed_min_len"],
                     n_truth=len(truth), n_neutral=len(neutral))
json.dump(out, open(ROOT / "results/analysis.json", "w"), indent=1, default=float)

# ================= figures
bm = json.load(open(ROOT / "results/benchmark.json"))
# Fig 1: F1 / precision / recall chloroplast
names = ["LongORF>=300", "Pyrodigal-single", "Pyrodigal-meta", "CodonLLR-only", "CodonRescore(single)", "CodonRescore(meta)"]
cols = [C["base"], C["pyro"], C["pyro"], C["alt"], C["ours"], C["ours"]]
m = bm["A. thaliana chloroplast"]["methods"]
fig, ax = plt.subplots(figsize=(3.4, 2.5))
x = np.arange(len(names)); w = 0.27
for i, (key, lab, hatch) in enumerate((("precision", "Precision", None), ("recall", "Recall", "///"), ("f1", "F1", "...."))):
    ax.bar(x + (i - 1) * w, [m[n][key] for n in names], w, color=cols, hatch=hatch, edgecolor="white", linewidth=0.4, label=lab)
ax.set_xticks(x); ax.set_xticklabels(["LongORF", "Pyro-\nsingle", "Pyro-\nmeta", "Codon\nLLR", "CR\n(single)", "CR\n(meta)"], fontsize=6.5)
ax.set_ylim(0, 1.05); ax.set_ylabel("Score")
h = [plt.Rectangle((0, 0), 1, 1, fc="#bbb", hatch=ht, ec="white") for ht in (None, "///", "....")]
ax.legend(h, ["Precision", "Recall", "F1"], fontsize=6.5, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.16))
fig.tight_layout(); fig.savefig(ROOT / "figs/fig_f1_bars.pdf"); fig.savefig(ROOT / "figs/fig_f1_bars.png"); plt.close()

# Fig 2: LLR distributions
fig, ax = plt.subplots(figsize=(3.4, 2.2))
bins_ = np.linspace(-0.6, 1.0, 40)
ax.hist(info["shadow_scores"], bins=bins_, alpha=0.75, color=C["base"], label="Shadow ORFs (n=%d)" % len(info["shadow_scores"]))
ax.hist(info["seed_scores"], bins=bins_, alpha=0.75, color=C["ours"], label="Seed ORFs (n=%d)" % len(info["seed_scores"]))
ax.axvline(tau, color="k", lw=1, ls="--"); ax.text(tau + 0.02, ax.get_ylim()[1] * 0.85, r"$\tau$=%.3f" % tau, fontsize=7)
ax.set_xlabel("Per-codon log-likelihood ratio (nats)"); ax.set_ylabel("ORFs"); ax.legend(fontsize=6.5, frameon=False)
fig.tight_layout(); fig.savefig(ROOT / "figs/fig_llr.pdf"); fig.savefig(ROOT / "figs/fig_llr.png"); plt.close()

# Fig 3: tau sensitivity
fig, ax = plt.subplots(figsize=(3.4, 2.2))
for b, col in (("single", C["pyro"]), ("meta", C["ours"])):
    ax.plot(offs, [s["f1"] for s in sens[b]], marker="o", ms=3, lw=1.4, color=col, label="CodonRescore(%s)" % b)
    ax.axhline(bm["A. thaliana chloroplast"]["methods"]["Pyrodigal-" + b]["f1"], color=col, lw=0.8, ls=":")
ax.set_xlabel(r"Threshold offset $\Delta\tau$"); ax.set_ylabel("F1"); ax.legend(fontsize=6.5, frameon=False)
ax.text(0.02, 0.5, "dotted = unmodified Pyrodigal\n(blue: single, orange: meta)", fontsize=6)
fig.tight_layout(); fig.savefig(ROOT / "figs/fig_sens.pdf"); fig.savefig(ROOT / "figs/fig_sens.png"); plt.close()

# Fig 4: length-stratified recall
fig, ax = plt.subplots(figsize=(3.4, 2.2))
labs = ["<300", "300-599", "600-1199", ">=1200"]
x = np.arange(4); w = 0.2
for i, (n, col, ht) in enumerate((("Pyrodigal-single", C["pyro"], "///"), ("Pyrodigal-meta", C["pyro"], None), ("CodonRescore(single)", C["ours"], "///"), ("CodonRescore(meta)", C["ours"], None))):
    ax.bar(x + (i - 1.5) * w, [a / b if b else 0 for a, b in strat[n]], w, color=col, hatch=ht, edgecolor="white", linewidth=0.4, label=n)
ax.set_xticks(x); ax.set_xticklabels([f"{l}\n(n={strat['Pyrodigal-meta'][i][1]})" for i, l in enumerate(labs)], fontsize=7)
ax.set_xlabel("True CDS length (nt)"); ax.set_ylabel("Recall"); ax.set_ylim(0, 1.0); ax.legend(fontsize=5.8, frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(); fig.savefig(ROOT / "figs/fig_lenrecall.pdf"); fig.savefig(ROOT / "figs/fig_lenrecall.png"); plt.close()
print(json.dumps(out, indent=1, default=float))

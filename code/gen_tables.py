"""Generate LaTeX tables (shared by both manuscript versions) from results JSON."""
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
bm = json.load(open(ROOT / "results/benchmark.json"))
an = json.load(open(ROOT / "results/analysis.json"))
T = ROOT / "paper"
f3 = lambda x: "%.3f" % x
ORDER = ["LongORF>=300", "Pyrodigal-single", "Pyrodigal-meta", "CodonLLR-only", "CodonRescore(single)", "CodonRescore(meta)"]
tex = lambda s: s.replace(">=", r"$\geq$")

# main chloroplast table with CI
m = bm["A. thaliana chloroplast"]["methods"]
rows = []
for n in ORDER:
    s = m[n]; lo, hi = an["bootstrap_f1_ci"][n]
    rows.append(f"{tex(n)} & {s['n_pred']} & {s['TP']} & {s['FP']} & {s['neutral']} & {f3(s['precision'])} & {f3(s['recall'])} & {f3(s['exact_recall'])} & {f3(s['f1'])} [{lo:.2f}--{hi:.2f}] \\\\")
(T / "tab_chloro.tex").write_text(
    "\\begin{tabular}{lrrrrrrrl}\n\\toprule\nMethod & Pred. & TP & FP & Neut. & Prec. & Rec. & Exact rec. & F1 [95\\% CI] \\\\\n\\midrule\n"
    + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

# all genomes (F1, P, R)
labs = list(bm.keys())
head = " & ".join(f"\\multicolumn{{3}}{{c}}{{{l}}}" for l in ["Chloroplast", "pPCP1 plasmid", "phiX174"])
rows = []
for n in ["LongORF>=300", "Pyrodigal-single", "Pyrodigal-meta", "CodonLLR-only", "CodonRescore(single)", "CodonRescore(meta)"]:
    cells = []
    for l in labs:
        if n in bm[l]["methods"]:
            s = bm[l]["methods"][n]
            cells += [f3(s["precision"]), f3(s["recall"]), f3(s["f1"])]
        else:
            cells += ["--"] * 3
    rows.append(tex(n) + " & " + " & ".join(cells) + " \\\\")
(T / "tab_all.tex").write_text(
    "\\begin{tabular}{lccccccccc}\n\\toprule\nMethod & " + head + " \\\\\n"
    "& P & R & F1 & P & R & F1 & P & R & F1 \\\\\n\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

# dataset
ch = an["chloro"]
gc = {"NC_005816": 45.3, "NC_001422": 44.8}
ds = [
    ("NC\\_000932.1", "\\textit{A. thaliana} chloroplast", "154,478", f"{100*ch['GC']:.1f}", "85 / 70 / 34", "Truth + training"),
    ("NC\\_005816.1", "\\textit{Y. pestis} plasmid pPCP1", "9,609", "45.3", "10 / 10 / 0", "Truth"),
    ("NC\\_001422.1", "Coliphage phiX174", "5,386", "44.8", "11 / 8 / 6", "Truth"),
    ("NZ\\_LN831026.1", "\\textit{C. diphtheriae} NCTC11397", "2,463,666", f"{100*an['cdiph']['GC']:.1f}", "not used", "Agreement only"),
]
(T / "tab_data.tex").write_text(
    "\\begin{tabular}{llrrcl}\n\\toprule\nAccession & Organism / record & Length (bp) & GC (\\%) & CDS: all / eval. / neutral$^\\dagger$ & Role \\\\\n\\midrule\n"
    + "\n".join(" & ".join(r) + " \\\\" for r in ds) + "\n\\bottomrule\n\\end{tabular}\n")

# sensitivity table
rows = []
for b in ("single", "meta"):
    for x in an["sens_tau"][b]:
        if x["offset"] in (-0.1, 0.0, 0.1, 0.2):
            rows.append(f"CR({b}) & {x['offset']:+.1f} & {f3(x['precision'])} & {f3(x['recall'])} & {f3(x['f1'])} \\\\")
(T / "tab_sens.tex").write_text("\\begin{tabular}{lrrrr}\n\\toprule\nVariant & $\\Delta\\tau$ & Prec. & Rec. & F1 \\\\\n\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

# C. diphtheriae
c = an["cdiph"]
rows = [
    f"Pyrodigal-single calls & {c['single']['n_pyro']:,} & -- \\\\",
    f"Pyrodigal-meta calls & {c['meta']['n_pyro']:,} & -- \\\\",
    f"Stop-codon concordance single vs.\\ meta (Jaccard) & {c['agree_single_meta_stop']:.3f} & -- \\\\",
    f"CodonRescore(single): pruned / rescued & {c['single']['pruned']} / {c['single']['rescued']} & median {c['single']['pruned_median_len']:.0f} / {c['single']['rescued_median_len']:.0f} nt \\\\",
    f"CodonRescore(meta): pruned / rescued & {c['meta']['pruned']} / {c['meta']['rescued']} & median {c['meta']['pruned_median_len']:.0f} / {c['meta']['rescued_median_len']:.0f} nt \\\\",
    f"Seed ORFs / Youden $J$ / $\\tau$ & {c['n_seeds']:,} / {c['youden']:.3f} / {c['tau']:.4f} & -- \\\\",
]
(T / "tab_cdiph.tex").write_text("\\begin{tabular}{lcc}\n\\toprule\nQuantity & Value & Note \\\\\n\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")
print("ok")

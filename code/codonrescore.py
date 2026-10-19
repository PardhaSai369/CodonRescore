"""CodonRescore: self-trained codon-model rescoring of ab initio gene calls.

Pure numpy implementation. No annotation is ever used for training; GenBank
annotation is used only by evaluate.py to score predictions.
"""
import numpy as np

COMP = str.maketrans("ACGT", "TGCA")
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}


def revcomp(s):
    return s.translate(COMP)[::-1]


def encode(s):
    a = np.frombuffer(s.encode(), dtype=np.uint8)
    out = np.full(a.shape, 4, dtype=np.int8)
    for ch, v in CODE.items():
        out[a == ord(ch)] = v
    return out


def _codon_ids(enc):
    """codon id at every position i (covering i..i+2); -1 when ambiguous."""
    n = len(enc)
    if n < 3:
        return np.zeros(0, dtype=np.int16)
    a, b, c = enc[:-2].astype(np.int16), enc[1:-1].astype(np.int16), enc[2:].astype(np.int16)
    ids = a * 16 + b * 4 + c
    ids[(a > 3) | (b > 3) | (c > 3)] = -1
    return ids


def _id(codon):
    return CODE[codon[0]] * 16 + CODE[codon[1]] * 4 + CODE[codon[2]]


STOP_IDS = np.array(sorted(_id(c) for c in STOPS))
START_IDS = np.array(sorted(_id(c) for c in STARTS))


def find_orfs(seq, min_len=90):
    """Longest-ORF-per-stop candidates on both strands (linear, complete only).

    Returns list of dicts with 1-based inclusive coords on the forward strand:
    begin, end, strand (+1/-1), length (nt, incl. stop), start codon.
    """
    L = len(seq)
    res = []
    for strand, s in ((1, seq), (-1, revcomp(seq))):
        ids = _codon_ids(encode(s))
        is_stop = np.isin(ids, STOP_IDS)
        is_start = np.isin(ids, START_IDS)
        for f in range(3):
            idx = np.arange(f, len(ids), 3)
            stops = idx[is_stop[idx]]
            starts = idx[is_start[idx]]
            prev = f - 3  # position of previous stop codon (virtual)
            for st in stops:
                lo = prev + 3
                j = np.searchsorted(starts, lo)
                if j < len(starts) and starts[j] < st:
                    b = int(starts[j])
                    ln = int(st) + 3 - b
                    if ln >= min_len:
                        if strand == 1:
                            res.append(dict(begin=b + 1, end=int(st) + 3, strand=1, length=ln, start=s[b:b + 3]))
                        else:
                            res.append(dict(begin=L - (int(st) + 3) + 1, end=L - b, strand=-1, length=ln, start=s[b:b + 3]))
                prev = int(st)
    return res


def orf_seq(seq, o):
    sub = seq[o["begin"] - 1:o["end"]]
    return sub if o["strand"] == 1 else revcomp(sub)


def overlap(a, b):
    return max(0, min(a["end"], b["end"]) - max(a["begin"], b["begin"]) + 1)


# ----------------------------------------------------------------- model
class PeriodicMarkov:
    """3-periodic order-2 Markov chain (codon-position specific) vs homogeneous order-2 background."""

    def __init__(self, alpha=1.0):
        self.alpha = alpha

    @staticmethod
    def _counts(seqs, periodic):
        n = np.zeros((3 if periodic else 1, 16, 4))
        for s in seqs:
            e = encode(s)
            if len(e) < 4:
                continue
            ctx = e[:-2].astype(int) * 4 + e[1:-1]
            nxt = e[2:].astype(int)
            ok = (e[:-2] < 4) & (e[1:-1] < 4) & (e[2:] < 4)
            pos = (np.arange(2, len(e)) % 3) if periodic else np.zeros(len(e) - 2, dtype=int)
            np.add.at(n, (pos[ok], ctx[ok], nxt[ok]), 1)
        return n

    def fit(self, coding_seqs, background_seqs):
        c = self._counts(coding_seqs, True) + self.alpha
        b = self._counts(background_seqs, False) + self.alpha
        self.logc = np.log(c / c.sum(2, keepdims=True))
        self.logb = np.log(b / b.sum(2, keepdims=True))
        return self

    def llr_per_codon(self, s):
        e = encode(s)
        if len(e) < 6:
            return 0.0
        ctx = e[:-2].astype(int) * 4 + e[1:-1]
        nxt = e[2:].astype(int)
        ok = (e[:-2] < 4) & (e[1:-1] < 4) & (e[2:] < 4)
        pos = np.arange(2, len(e)) % 3
        d = self.logc[pos[ok], ctx[ok], nxt[ok]] - self.logb[0, ctx[ok], nxt[ok]]
        return float(d.sum() / (len(e) / 3.0))


# ------------------------------------------------------- self-training
def pick_seeds(orfs, min_len=450, max_ov_frac=0.3):
    """Greedy: longest ORFs first; reject if overlapping an accepted seed by > max_ov_frac of itself."""
    seeds = []
    for o in sorted((o for o in orfs if o["length"] >= min_len), key=lambda o: -o["length"]):
        if all(overlap(o, s) <= max_ov_frac * o["length"] for s in seeds):
            seeds.append(o)
    return seeds


def shadow_orfs(orfs, seeds, min_frac=0.5):
    """ORFs (>=90 nt) on other frames/strand that overlap a seed by >= min_frac of themselves."""
    seedset = {(s["begin"], s["end"], s["strand"]) for s in seeds}
    out = []
    for o in orfs:
        if (o["begin"], o["end"], o["strand"]) in seedset:
            continue
        for s in seeds:
            if overlap(o, s) >= min_frac * o["length"]:
                out.append(o)
                break
    return out


def train_codonrescore(seq, orfs=None, seed_min_len=450):
    """Returns (model, tau, info). tau from 2-fold cross-fitted Youden threshold, seeds vs shadow ORFs."""
    if orfs is None:
        orfs = find_orfs(seq)
    seeds = []
    for ml in (seed_min_len, 360, 300):
        seeds = pick_seeds(orfs, ml)
        if len(seeds) >= 10:
            break
    seed_min_len = ml
    if len(seeds) < 10:
        raise ValueError("too few seed ORFs (%d)" % len(seeds))
    bg = [seq, revcomp(seq)]
    folds = [seeds[0::2], seeds[1::2]]
    sc_seed, sc_shadow = [], []
    for k in range(2):
        train = [orf_seq(seq, o) for o in folds[1 - k]]
        m = PeriodicMarkov().fit(train, bg)
        held = folds[k]
        sc_seed += [m.llr_per_codon(orf_seq(seq, o)) for o in held]
        sh = shadow_orfs(orfs, held)
        sc_shadow += [m.llr_per_codon(orf_seq(seq, o)) for o in sh]
    sc_seed, sc_shadow = np.array(sc_seed), np.array(sc_shadow)
    cand = np.unique(np.concatenate([sc_seed, sc_shadow]))
    best, tau = -9, 0.0
    for t in cand:
        j = (sc_seed >= t).mean() + (sc_shadow < t).mean() - 1
        if j > best:
            best, tau = j, float(t)
    model = PeriodicMarkov().fit([orf_seq(seq, o) for o in seeds], bg)
    info = dict(seed_min_len=seed_min_len, n_seeds=len(seeds), n_shadow=len(sc_shadow), youden=float(best), tau=tau,
                seed_scores=sc_seed, shadow_scores=sc_shadow)
    return model, tau, info


# ------------------------------------------------- calling procedures
def _resolve(cands, max_ov=60, opp_frac=0.6):
    """Greedy conflict resolution by total log-likelihood ratio (score*codons)."""
    kept = []
    for o in sorted(cands, key=lambda o: -o["llr_total"]):
        ok = True
        for k in kept:
            ov = overlap(o, k)
            if ov == 0:
                continue
            if o["strand"] == k["strand"] and ov > max_ov:
                ok = False
                break
            if o["strand"] != k["strand"] and ov > opp_frac * min(o["length"], k["length"]):
                ok = False
                break
        if ok:
            kept.append(o)
    return sorted(kept, key=lambda o: o["begin"])


def annotate(seq, genes, model):
    for g in genes:
        g["llr"] = model.llr_per_codon(orf_seq(seq, g))
        g["llr_total"] = g["llr"] * g["length"] / 3.0
    return genes


def call_codon_only(seq, model, tau, orfs=None, min_len=90):
    orfs = find_orfs(seq, min_len) if orfs is None else orfs
    cands = annotate(seq, [dict(o) for o in orfs], model)
    cands = [c for c in cands if c["llr"] >= tau]
    return _resolve(cands)


def call_hybrid(seq, pyro_genes, model, tau, orfs=None, prune_len=300, rescue_min_len=120):
    """CodonRescore: prune weak short Pyrodigal calls, rescue strong uncalled ORFs."""
    P = annotate(seq, [dict(g) for g in pyro_genes], model)
    kept = [g for g in P if not (g["llr"] < tau and g["length"] < prune_len)]
    pruned = [g for g in P if g not in kept]
    orfs = find_orfs(seq, rescue_min_len) if orfs is None else orfs
    cands = annotate(seq, [dict(o) for o in orfs if o["length"] >= rescue_min_len], model)
    rescued = []
    for c in sorted(cands, key=lambda c: -c["llr_total"]):
        if c["llr"] < tau:
            continue
        if any(overlap(c, k) > 60 or (c["strand"] != k["strand"] and overlap(c, k) > 0.6 * min(c["length"], k["length"]))
               for k in kept + rescued):
            continue
        rescued.append(c)
    final = sorted(kept + rescued, key=lambda g: g["begin"])
    return final, pruned, rescued

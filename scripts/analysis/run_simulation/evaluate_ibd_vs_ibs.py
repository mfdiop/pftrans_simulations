#!/usr/bin/env python3

"""
evaluate_ibd_vs_ibs.py

Inputs expected per scenario (folder):
    - simulated_samples.vcf[.gz]            (VCF for the scenario, same sample ordering as pedigree)
    - pedigree_simplified.csv               (columns: parent_individual, child_individual, parent_time, child_time)
    - sample_map.json                       (maps sample_name -> individual_id, optional coords)
    - ibd_segments.csv                      (columns: sample1, sample2, start, end, length_bp)  -- produced by IBD caller
    - config.json                           (optional, to read recombination rate used)

Outputs:
    - PR-AUC and PR curves for IBD and IBS
    - PR-AUC(IBD) - PR-AUC(IBS) vs recomb_rate plot
    - max k (generations) with recall@precision>=0.8
    - recall vs geographic distance bins
    - scenario_summary.csv (one row per scenario with key metrics)

 Syntax:
    python evaluate_ibd_vs_ibs.py --scenarios scenario1_dir scenario2_dir ... --outdir output_dir  
"""

import os, json, math, argparse
from collections import defaultdict
import numpy as np # pyright: ignore[reportMissingImports]
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import allel  # scikit-allel
import tskit
import networkx as nx
from sklearn.metrics import precision_recall_curve, auc, average_precision_score

# -------------------------
# Utility parsers & helpers
# -------------------------
def load_sample_map(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)

def load_ibd_segments(path):
    """
    Generic loader for IBD callers. Expects columns: sample1, sample2, start, end, length_bp
    """
    df = pd.read_csv(path)
    # normalize sample names
    df['sample1'] = df['sample1'].astype(str)
    df['sample2'] = df['sample2'].astype(str)
    if 'length_bp' not in df.columns:
        if 'end' in df.columns and 'start' in df.columns:
            df['length_bp'] = df['end'] - df['start']
        elif 'length' in df.columns:
            df['length_bp'] = df['length']
        else:
            raise ValueError("IBD segments need length (bp) or start/end columns")
    return df

def ibd_total_length_per_pair(ibd_df):
    """
    Sum total IBD length per unordered pair -> returns dict {(s1,s2): total_bp}
    Unordered pair uses tuple(sorted((a,b))).
    """
    d = defaultdict(float)
    for _, r in ibd_df.iterrows():
        a, b = str(r['sample1']), str(r['sample2'])
        key = tuple(sorted((a,b)))
        d[key] += float(r['length_bp'])
    return d

def compute_ibs_matrix_from_vcf(vcf_path, samples=None):
    """
    Compute pairwise IBS distance (proportion mismatching alleles) using scikit-allel.
    Returns dict {(a,b): ibs_distance} where smaller = more similar.
    """
    callset = allel.read_vcf(vcf_path, fields=['samples','calldata/GT','variants/POS'])
    gt = callset['calldata/GT']  # shape (n_variants, n_samples, ploidy)
    samples_list = list(callset['samples'])
    if samples is not None:
        # restrict to specified sample list (by name)
        idx = [samples_list.index(s) for s in samples]
        gt = gt[:, idx, :]
        samples_list = [samples_list[i] for i in idx]

    # For haploid data, GT will have last dim=1; collapse to 2D array (variants x samples)
    if gt.shape[2] == 1:
        geno = gt[:, :, 0].astype('int8')
    else:
        # for diploid, convert to allele counts or flatten haplotypes; here we compute mismatch proportion on major allele
        # simplest: convert to max allele per sample per site (not ideal); assuming haploid simulation this won't be used
        geno = gt[:,:,0]
    n_variants = geno.shape[0]
    n_samples = geno.shape[1]
    ibs = {}
    # vectorized mismatch counts using numpy / bit operations if alleles encoded 0/1
    for i in range(n_samples):
        ai = geno[:, i]
        for j in range(i+1, n_samples):
            aj = geno[:, j]
            # ignore sites with missing (-1 or 255)
            mask = (ai >= 0) & (aj >= 0)
            if mask.sum() == 0:
                dist = 1.0
            else:
                mismatches = (ai[mask] != aj[mask]).sum()
                dist = mismatches / mask.sum()
            s1, s2 = samples_list[i], samples_list[j]
            ibs[(s1, s2)] = dist
    return ibs

def pairs_from_sample_list(samples):
    pairs = []
    n = len(samples)
    for i in range(n):
        for j in range(i+1, n):
            pairs.append((samples[i], samples[j]))
    return pairs

# -------------------------
# Build ground truth pair sets
# -------------------------
def build_ground_truth_pairs(pedigree_csv, sample_map=None):
    """
    Input: pedigree CSV with parent_individual, child_individual columns.
    sample_map: maps sample_name -> individual_id (or inverse)
    Return:
      GT_pairs: set of frozenset({sampleA, sampleB}) for those pairs where both samples correspond to child/parent
      Also returns a mapping sample_name->individual_id used.
    """
    ped = pd.read_csv(pedigree_csv)
    # we expect ped uses individual ids, map to sample ids via sample_map if provided
    # sample_map can be: { "sample_name": individual_id } OR { "individual_id": "sample_name" }
    if sample_map is None:
        # assume sample names are individual ids as strings
        # produce GT on individual ids as strings
        gt = set()
        for _, r in ped.iterrows():
            a = str(int(r['parent_individual']))
            b = str(int(r['child_individual']))
            gt.add(frozenset((a, b)))
        return gt, None

    # determine mapping direction
    sm = sample_map
    # convert numpy types to native
    sm = {str(k): (int(v) if isinstance(v, (int, np.integer)) else v) for k,v in sm.items()}
    # if map is sample_name->individual_id, invert to individual->sample
    # Heuristic: keys look like sample names (strings with letters) and values ints -> sample->ind
    sample_to_ind = {}
    ind_to_sample = {}
    # allow both possibilities by checking types
    first_val = next(iter(sm.values()))
    if isinstance(first_val, int):
        # map is sample->ind
        for s, ind in sm.items():
            sample_to_ind[s] = int(ind)
            ind_to_sample[str(int(ind))] = s
    else:
        # map is ind->sample
        for ind, s in sm.items():
            ind_to_sample[str(int(ind))] = s
            sample_to_ind[s] = int(ind)

    gt = set()
    sample_list = set(ind_to_sample.values())
    for _, r in ped.iterrows():
        p = str(int(r['parent_individual']))
        c = str(int(r['child_individual']))
        s_p = ind_to_sample.get(p, None)
        s_c = ind_to_sample.get(c, None)
        if s_p is None or s_c is None:
            continue
        gt.add(frozenset((s_p, s_c)))
    return gt, ind_to_sample

# -------------------------
# Compute PR-AUC given scores and ground truth labels
# -------------------------
def compute_pr_auc_for_pairs(pair_list, score_dict, GT_set):
    """
    pair_list: list of (s1,s2) tuples (ordered)
    score_dict: dict with unordered pair keys (sorted) -> score (higher = more evidence of link)
    GT_set: set of frozenset pairs (ground truth, unordered)
    Returns: precision, recall, thresholds, pr_auc
    """
    y_true = []
    y_score = []
    for (a,b) in pair_list:
        key = tuple(sorted((a,b)))
        y_score.append(score_dict.get(key, 0.0))
        y_true.append(1 if frozenset((a,b)) in GT_set else 0)
    y_true = np.array(y_true)
    y_score = np.array(y_score)
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    pr_auc = auc(recall, precision)
    return precision, recall, thresholds, pr_auc, y_true, y_score

# -------------------------
# Convert IBD length -> score
# -------------------------
def ibd_score_from_segments(ibd_length_dict, pair_list):
    """
    We use total length (bp) as the score (larger = more evidence).
    """
    scores = {}
    for (a,b) in pair_list:
        scores[tuple(sorted((a,b)))] = ibd_length_dict.get(tuple(sorted((a,b))), 0.0)
    return scores

# -------------------------
# Evaluate scenario
# -------------------------
def evaluate_scenario(scenario_dir, vcf_path, pedigree_csv, ibd_segments_csv, sample_map_json=None, recomb_rate=None, outdir=None):
    os.makedirs(outdir, exist_ok=True)
    # 1) load sample_map
    sample_map = load_sample_map(sample_map_json) if sample_map_json else None

    # 2) load IBD segments
    ibd_df = load_ibd_segments(ibd_segments_csv)
    ibd_length_dict = ibd_total_length_per_pair(ibd_df)

    # 3) compute IBS distances from VCF (lower=more similar). Convert to similarity score = (1 - distance)
    ibs_dict = compute_ibs_matrix_from_vcf(vcf_path)
    ibs_score = {tuple(sorted(k)): 1.0 - v for k,v in ibs_dict.items()}  # higher = more similar

    # 4) build pair list (all unordered pairs present in either ibd or ibs)
    samples = sorted(list(set([s for pair in ibs_dict.keys() for s in pair])))
    pair_list = pairs_from_sample_list(samples)

    # 5) build ground truth set
    GT_set, ind_to_sample = build_ground_truth_pairs(pedigree_csv, sample_map)

    # 6) IBD PR/AUC
    ibd_scores = ibd_score_from_segments(ibd_length_dict, pair_list)
    p_ibd, r_ibd, thr_ibd, pr_auc_ibd, y_true, y_score_ibd = compute_pr_auc_for_pairs(pair_list, ibd_scores, GT_set)

    # 7) IBS PR/AUC
    # ibs_score already prepared
    p_ibs, r_ibs, thr_ibs, pr_auc_ibs, _, y_score_ibs = compute_pr_auc_for_pairs(pair_list, ibs_score, GT_set)

    # Save PR curves and AUC
    summary = {
        "scenario_dir": scenario_dir,
        "recomb_rate": recomb_rate,
        "pr_auc_ibd": float(pr_auc_ibd),
        "pr_auc_ibs": float(pr_auc_ibs),
        "pr_auc_diff": float(pr_auc_ibd - pr_auc_ibs),
        "n_pairs": len(pair_list),
        "n_gt_pairs": int(sum(y_true))
    }
    # Save curves
    pd.DataFrame({"precision": p_ibd, "recall": r_ibd}).to_csv(os.path.join(outdir, "pr_ibd.csv"), index=False)
    pd.DataFrame({"precision": p_ibs, "recall": r_ibs}).to_csv(os.path.join(outdir, "pr_ibs.csv"), index=False)

    # Plot PR curves
    plt.figure(figsize=(6,5))
    plt.plot(r_ibd, p_ibd, label=f'IBD PR-AUC={pr_auc_ibd:.3f}')
    plt.plot(r_ibs, p_ibs, label=f'IBS PR-AUC={pr_auc_ibs:.3f}')
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.legend()
    plt.title(f"PR curves (scenario: {os.path.basename(scenario_dir)})")
    plt.savefig(os.path.join(outdir, "pr_curves.png"), bbox_inches="tight")
    plt.close()

    # 8) For each k compute recall@precision>=0.8
    # Build pedigree graph (directed parent->child) and compute shortest path lengths (generations)
    ped = pd.read_csv(pedigree_csv)
    # Build undirected graph for generation distances: nodes = individual ids or sample names (use sample names)
    G = nx.Graph()
    # map individuals to sample names if sample_map provided
    if sample_map:
        # build ind->sample mapping; if sample_map is sample->ind invert
        sm = sample_map
        first_val = next(iter(sm.values()))
        if isinstance(first_val, int):
            ind_to_sample = {str(v): k for k, v in sm.items()}
        else:
            ind_to_sample = {str(k): v for k,v in sm.items()}
    else:
        ind_to_sample = {}

    for _, r in ped.iterrows():
        p = str(int(r['parent_individual']))
        c = str(int(r['child_individual']))
        s_p = ind_to_sample.get(p, p)  # fallback to using individual id string
        s_c = ind_to_sample.get(c, c)
        G.add_edge(s_p, s_c)

    # consider all ks up to some max
    max_k = 20
    best_k = -1
    recall_at_prec_thresh = {}
    for k in range(1, max_k+1):
        # Build GT_k: pairs with shortest path <= k
        GTk = set()
        nodes_list = list(G.nodes())
        for i in range(len(nodes_list)):
            for j in range(i+1, len(nodes_list)):
                a = nodes_list[i]; b = nodes_list[j]
                try:
                    d = nx.shortest_path_length(G, a, b)
                except nx.NetworkXNoPath:
                    continue
                if d <= k:
                    GTk.add(frozenset((a, b)))
        # compute PR curve for IBD scores but treat y_true_k accordingly
        _, _, _, pr_auc_k, y_true_k, _ = compute_pr_auc_for_pairs(pair_list, ibd_scores, GTk)
        # Compute precision at desired recall/threshold: find threshold achieving precision>=0.8 and report recall
        precision, recall, thresholds = precision_recall_curve([1 if frozenset((a,b)) in GTk else 0 for (a,b) in pair_list],
                                                              [ibd_scores[tuple(sorted((a,b)))] for (a,b) in pair_list])
        # find max recall where precision >= 0.8
        pr_arr = np.array(precision); rec_arr = np.array(recall)
        idxs = np.where(pr_arr >= 0.8)[0]
        rec_at_prec = float(rec_arr[idxs].max()) if len(idxs) > 0 else 0.0
        recall_at_prec_thresh[k] = rec_at_prec
        if rec_at_prec >= 0.8:
            best_k = k

    # 9) Spatial bins: if sample_map includes coords (x,y) compute distance bins and recall per bin
    recall_vs_dist = None
    if sample_map and 'coords' in sample_map.get(next(iter(sample_map)), {}):
        # sample_map assumed: sample_name -> {"ind": id, "coords": [x,y]}
        coords = {s: tuple(sample_map[s]['coords']) for s in sample_map}
        # compute pairwise distances
        pair_dists = []
        pairs = pair_list
        for (a,b) in pairs:
            xa, ya = coords[a]; xb, yb = coords[b]
            dist = math.hypot(xa-xb, ya-yb)
            pair_dists.append(dist)
        # define bins
        bins = np.linspace(0, max(pair_dists)+1e-6, 6)
        bin_idx = np.digitize(pair_dists, bins)
        # choose ibd threshold achieving precision >= 0.8 (global), find threshold
        precision, recall, thresholds = precision_recall_curve([1 if frozenset((a,b)) in GT_set else 0 for (a,b) in pair_list],
                                                              [ibd_scores[tuple(sorted((a,b)))] for (a,b) in pair_list])
        idxs = np.where(np.array(precision) >= 0.8)[0]
        chosen_thr = thresholds[idxs[-1]] if len(idxs)>0 else np.percentile(list(ibd_scores.values()), 90)
        # compute recall per bin
        bin_recalls = {}
        for b in range(1, len(bins)+1):
            idxs_b = [i for i,v in enumerate(pair_list) if bin_idx[i]==b]
            if len(idxs_b)==0:
                bin_recalls[b] = np.nan
                continue
            # selected positives by threshold
            selected = [1 if ibd_scores[tuple(sorted(pair_list[i]))] >= chosen_thr else 0 for i in idxs_b]
            truth = [1 if frozenset(pair_list[i]) in GT_set else 0 for i in idxs_b]
            tp = sum([1 for s,t in zip(selected, truth) if s==1 and t==1])
            fn = sum(truth) - tp
            recall = tp / (tp + fn) if (tp+fn)>0 else np.nan
            bin_recalls[b] = recall
        recall_vs_dist = {'bins': bins.tolist(), 'recalls': bin_recalls}

    # Save summary CSV
    scenario_summary = {
        'pr_auc_ibd': pr_auc_ibd, 'pr_auc_ibs': pr_auc_ibs, 'pr_auc_diff': pr_auc_ibd - pr_auc_ibs,
        'best_k_recall80': best_k, 'recall_per_k': recall_at_prec_thresh
    }
    with open(os.path.join(outdir, 'scenario_summary.json'), 'w') as fh:
        json.dump(scenario_summary, fh, indent=2)

    return scenario_summary

# -------------------------
# CLI
# -------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenarios", nargs="+", help="List of scenario directories (each must contain VCF, pedigree CSV and IBD segments file)", required=True)
    parser.add_argument("--vcf_name", default="simulated_samples.vcf", help="VCF filename inside scenario dir")
    parser.add_argument("--pedigree_name", default="pedigree_simplified.csv", help="pedigree filename inside scenario dir")
    parser.add_argument("--ibd_name", default="ibd_segments.csv", help="IBD segments filename inside scenario dir")
    parser.add_argument("--sample_map_name", default="sample_map.json", help="sample map filename inside scenario dir")
    parser.add_argument("--outdir", default="eval_out", help="Output directory for evaluation")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    summary_rows = []
    for scen in args.scenarios:
        vcf_path = os.path.join(scen, args.vcf_name)
        ped_path = os.path.join(scen, args.pedigree_name)
        ibd_path = os.path.join(scen, args.ibd_name)
        sample_map = os.path.join(scen, args.sample_map_name)
        outdir = os.path.join(args.outdir, os.path.basename(scen))

        os.makedirs(outdir, exist_ok=True)

        # attempt to read recomb rate from a config.json in the scenario dir
        config_path = os.path.join(scen, "config.json")
        recomb_rate = None

        if os.path.exists(config_path):
            conf = json.load(open(config_path))
            recomb_rate = conf.get("recombination_rate", None)

        summary = evaluate_scenario(scen, vcf_path, ped_path, ibd_path, sample_map_json=sample_map, recomb_rate=recomb_rate, outdir=outdir)
        summary_rows.append({'scenario': scen, 'recomb_rate': recomb_rate, **summary})

    # Save combined summary table
    pd.DataFrame(summary_rows).to_csv(os.path.join(args.outdir, "all_scenarios_summary_table.csv"), index=False)
    
    # Plot PR-AUC diff vs recomb_rate if available
    df = pd.DataFrame(summary_rows)
    if 'recomb_rate' in df.columns and df['recomb_rate'].notnull().all():
        df['recomb_rate'] = df['recomb_rate'].astype(float)
        plt.figure(figsize=(6,4))
        sns.lineplot(x='recomb_rate', y='pr_auc_diff', data=df, marker='o')
        plt.xscale('log')
        plt.xlabel('Recombination rate (per bp per gen)')
        plt.ylabel('PR-AUC(IBD) - PR-AUC(IBS)')
        plt.title('IBD vs IBS performance across recombination rates')
        plt.savefig(os.path.join(args.outdir, "pr_auc_diff_vs_recomb.png"), bbox_inches='tight')
        plt.close()

    print("Evaluation complete. See output:", args.outdir)

if __name__ == "__main__":
    main()

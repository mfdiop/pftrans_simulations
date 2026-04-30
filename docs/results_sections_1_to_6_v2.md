# Results — Sections 1–6

### Revised following strict results-section style constraints

------------------------------------------------------------------------

## Section 1

### Transmission-linked parasite pairs retain only partial genealogical separation, establishing a biological ceiling on detectable signal

Figure 1 shows the distribution of true pairwise coalescence times (TMRCA) for directly transmitted and non-directly transmitted parasite pairs under two simulation conditions: baseline parameters (left panel) and a recombination sweep (right panel). In both panels, the x-axis represents generational depth to the most recent common ancestor (G, from 0 to 50), and distributions are shown separately for direct pairs (red) and non-direct pairs (grey).

Under baseline conditions, the distribution of direct pairs was heavily concentrated at shallow genealogical depths, with the majority of mass below G = 5 and a rapid decline toward zero at greater depths. The distribution of non-direct pairs was broadly dispersed across the full range of G values shown, with a flat to gradually declining profile. The two distributions overlapped over the G = 3–15 range, where both classes contributed non-trivial density. At G values above approximately 20, the direct pair distribution approached zero while non-direct pairs retained detectable density.

Under the recombination sweep condition (Figure 1, right panel), the distribution of direct pairs was qualitatively similar to baseline. The non-direct pair distribution showed increased mass at intermediate genealogical depths relative to baseline, with a broader and more uniform profile. The overlap region between the two distributions was wider under the sweep than at baseline, extending further toward shallower genealogical depths.

In both panels, a subset of non-direct pairs coalesced at depths overlapping with the bulk of the direct pair distribution. This overlap was present at baseline and was not eliminated by the recombination sweep.

> **Figure 1.** TMRCA distributions for directly and non-directly transmitted parasite pairs. Left: baseline simulation parameters. Right: recombination sweep condition. Direct pairs (red); non-direct pairs (grey). Overlap region visible between G = 3 and G = 15 in both panels; the overlap region is wider under the recombination sweep.

Supplementary Figure S1 shows pairwise IBD fraction distributions at baseline for both true and hmmIBD-inferred values across direct and non-direct pair classes, providing the sequence-level counterpart to the genealogical distributions in Figure 1.

------------------------------------------------------------------------

## Section 2

### Recombination progressively erodes the genomic signal of shared ancestry, widening the gap between available and recoverable transmission information

Figure 2 summarises how recombination rate affects both the availability and recoverability of IBD signal across the five recombination rates examined (r = 10\^-9, 10\^-8, 6.67 × 10\^-7, and 10\^-6 per base pair per generation). The figure contains four panels, each addressing a distinct aspect of signal availability or inference fidelity.

Panel A shows the distribution of true IBD proportions for direct (red) and non-direct (grey) pairs at each recombination rate. At r = 10\^-9, both distributions showed elevated IBD proportions relative to higher rates, with direct pairs displaying a wide spread toward high IBD values and non-direct pairs showing a broad background distribution. As recombination rate increased, both distributions shifted toward lower IBD values. The separation between the direct and non-direct distributions was most pronounced at r = 6.67 × 10\^-7, where the direct pair distribution retained a rightward tail absent in the non-direct distribution. At r = 10\^-6, both distributions converged toward low IBD values, with reduced separation between classes.

Panel B shows the corresponding inferred IBD distributions (hmmIBD output) at the same recombination rates. At r = 10\^-9, inferred IBD values for non-direct pairs were elevated relative to true IBD values, producing a rightward displacement of the non-direct distribution toward the direct pair range. At r = 10\^-6, inferred IBD values for direct pairs were reduced relative to true IBD, producing a leftward shift of the direct pair distribution toward the non-direct range. At intermediate rates, the inferred distributions more closely tracked the true distributions shown in Panel A.

Panel C shows the divergence between mean true and mean inferred IBD as a function of recombination rate, separately for direct and non-direct pairs. The divergence for non-direct pairs was positive at r = 10\^-9 (inferred exceeding true) and decreased toward zero at intermediate rates before reversing sign at r = 10\^-6. The divergence for direct pairs was negative at high recombination rates (true exceeding inferred) and smallest in magnitude at r ≈ 6.67 × 10\^-7. The two divergence curves crossed near the intermediate recombination rate.

Panel D shows the Pearson correlation between true and inferred IBD across recombination rates. The correlation was lowest at r = 10\^-9 (r ≈ 0.53, 95% CI 0.49–0.58), increased to a peak of approximately 0.95 at r = 10\^-7, then declined to approximately 0.65 at r = 10\^-6. The profile was non-monotonic across the full range of recombination rates tested.

> **Figure 2.** Effect of recombination rate on the availability and recoverability of transmission signal. **(A)** True IBD proportion distributions for direct (red) and non-direct (grey) pairs across recombination rates. **(B)** Inferred IBD proportions (hmmIBD) at the same rates. **(C)** Divergence between mean true and mean inferred IBD as a function of recombination rate. **(D)** Pearson correlation between true and inferred IBD across recombination rates. Non-monotonic correlation profile peaks at r ≈ 10\^-7 and declines at both extremes.

Supplementary Figure S2 shows the effect of transmission bottleneck size (Ne = 1,000–20,000) on AUPRC across methods and recombination conditions. AUPRC differences attributable to bottleneck size ranged from 0.011 to 0.014 across the levels tested. Supplementary Table S1 lists the full simulation parameter space, replicate counts, and random seeds for all scenarios in this section.

------------------------------------------------------------------------

## Section 3 ★

### Transmission detection undergoes a threshold transition governed by recombination rate, independent of inference approach

Figure 3 characterises how detection performance varied across recombination rates and inference approaches using two complementary metrics: F1 score as a function of recombination rate (Panel A) and precision-recall curves under contrasting simulation scenarios (Panel B).

Panel A shows F1 score on the y-axis against recombination rate (log scale, x-axis) for IBD (red), IBS (blue), and phylogenetic cophenetic distance (Phylo, green), under baseline parameters. All three approaches produced F1 scores in the range of 0.02–0.07 at recombination rates of 10\^-9 and 10\^-8, with no approach exceeding F1 = 0.10 at these rates under any parameter combination examined in the full factorial analysis. Between r = 10\^-8 and r = 6.67 × 10\^-7, F1 scores increased sharply for IBD and IBS, with IBD reaching F1 ≈ 0.54 and IBS reaching F1 ≈ 0.48 at r = 6.67 × 10\^-7. Phylo showed a smaller and more variable increase over the same range. The transition from near-zero F1 to values above 0.40 was confined to approximately one order of magnitude of recombination rate for IBD and IBS. The location of this transition was consistent across the three migration rates and three sample sizes included in the full factorial analysis. At r = 10\^-6, IBD and IBS F1 scores were comparable to or slightly above values at r = 6.67 × 10\^-7; Phylo showed heterogeneous responses at the highest recombination rate.

Expressed as AUPRC, no method exceeded 0.62 at r = 10\^-9 under any demographic or surveillance configuration tested (Supplementary Figure S4; Supplementary Table S2). At r = 10\^-6 and n = 100, best AUPRC across methods reached 0.79–0.82. At n ≥ 200, no parameter combination produced AUPRC \> 0.80.

Panel B shows precision-recall curves for all three approaches under three scenarios: baseline (solid lines), recombination sweep (dashed lines), and the full factorial condition aggregated across parameter combinations (shaded envelopes). Each curve plots precision on the y-axis against recall on the x-axis. Under baseline, IBD reached the highest precision values at low recall levels, with the curve positioned above IBS, which in turn was above Phylo across most of the recall range. Under the recombination sweep, precision values at any fixed recall level were lower for all three approaches relative to baseline. The reduction in precision was most pronounced for Phylo and smallest for IBD. Under the full factorial scenario, the outer envelope of the precision-recall space — the maximum precision achievable at each recall level across all approaches — contracted relative to the baseline and sweep conditions, with the contraction most visible at low-to-intermediate recall values.

> **Figure 3.** **(A)** F1 score as a function of recombination rate for IBD (red), IBS (blue), and Phylo (green) under baseline parameters. Near-zero F1 at r ≤ 10\^-8 for all approaches; sharp increase between 10\^-8 and 6.67 × 10\^-7 for IBD and IBS. **(B)** Precision–recall curves under baseline (solid), recombination sweep (dashed), and full factorial (shaded envelopes) conditions. Precision values at fixed recall levels decrease from baseline to sweep to factorial for all approaches; the outer achievable envelope contracts progressively across scenarios.

Supplementary Figure S4 shows the fraction of simulation replicates achieving AUPRC \> 0.70 across the sample size × recombination rate parameter space for each method. Supplementary Table S2 provides complete AUPRC values with 95% bootstrap confidence intervals for all method, scenario, recombination rate, sample size, and migration combinations.

------------------------------------------------------------------------

## Section 4

### The detection threshold reflects coalescent geometry: genealogical overlap between transmission classes propagates into sequence-based representations

Figure 4 shows pairwise cophenetic distance distributions for directly and non-directly transmitted parasite pairs under baseline (left panel) and recombination sweep (right panel) simulation conditions. Cophenetic distances were derived from the maximum-likelihood phylogenetic tree reconstructed by IQ-TREE for each simulation replicate.

Under baseline conditions, the cophenetic distance distributions for direct (red) and non-direct (grey) pairs showed overlapping profiles across most of the distance range. The direct pair distribution was concentrated at shorter cophenetic distances, with a peak at lower values and a right tail extending to intermediate distances. The non-direct pair distribution spanned a broader range, with mass at both short and long distances. Compared to the TMRCA distributions in Figure 1 (left panel), the overlap between direct and non-direct classes was more extensive in cophenetic distance space: the direct pair distribution extended further into the range occupied by non-direct pairs, and fewer distance values were exclusive to one class.

Under the recombination sweep (Figure 4, right panel), the cophenetic distance distributions for both pair classes shifted and converged relative to baseline. The direct pair distribution broadened, with reduced mass at the shortest distances and increased mass at intermediate values. The non-direct pair distribution showed a leftward shift, with more mass concentrated at shorter cophenetic distances. The net effect was a reduced separation between the two distributions under the sweep compared to baseline. The convergence between classes was more pronounced in cophenetic distance space than in the corresponding TMRCA comparison between the two panels of Figure 1.

> **Figure 4.** Cophenetic distance distributions for directly and non-directly transmitted pairs. Left: baseline parameters. Right: recombination sweep. Direct pairs (red); non-direct pairs (grey). Overlap between classes is more extensive in cophenetic distance space than in the TMRCA distributions shown in Figure 1; convergence between classes is greater under the sweep than at baseline.

Supplementary Figure S6 provides mean AUPRC values with 95% confidence intervals as continuous line plots across the full recombination rate range for all three methods under each migration condition, showing the full performance gradient for all approaches including the phylogenetic method characterised here.

------------------------------------------------------------------------

## Section 5

### Population structure modulates but cannot override the recombination-governed detectability ceiling

Figure 5 summarises how migration rate (m = 0.001, 0.01, and 0.05) affected AUPRC across recombination rates and inference approaches. Panels A through D show mean AUPRC (± 95% CI) as a function of recombination rate for IBD (A), IBS (B), Phylo (C), and best AUPRC across methods (D), each with three lines corresponding to the three migration rates. Panel E shows the AUPRC gradient across the full recombination rate × migration space as a summary characterisation.

In Panels A through D, AUPRC increased with recombination rate for all three approaches and at all three migration rates. The transition from low to higher AUPRC values occurred at the same recombination rate stratum regardless of migration condition, with the rank ordering of recombination rates preserved across all migration values. Within any given recombination stratum, differences in AUPRC attributable to migration rate were present but small. The maximum within-stratum difference between the lowest and highest migration conditions did not exceed 0.024 AUPRC units across all method and recombination rate combinations. The direction of the migration effect was consistent: higher migration rates were associated with marginally higher AUPRC values in all panels, most visibly at recombination rates above r = 6.67 × 10\^-7.

Below r = 10\^-8, all three migration conditions produced similarly low AUPRC values for all three approaches, with overlapping confidence intervals across migration conditions at this recombination stratum. Above r = 6.67 × 10\^-7, the three migration lines were more separated but remained within a narrow band, with confidence intervals overlapping in most panels.

Panel E shows the AUPRC gradient as a two-dimensional summary across recombination rate (x-axis) and migration rate (y-axis). The dominant gradient in Panel E was oriented along the recombination axis. Variation along the migration axis at any fixed recombination rate was visually narrow relative to the variation along the recombination axis across the full range tested.

> **Figure 5.** Migration rate effects on AUPRC across recombination rates and methods. **(A–D)** Mean AUPRC ± 95% CI for IBD (A), IBS (B), Phylo (C), and best across methods (D) at migration rates m = 0.001 (red), m = 0.01 (orange), and m = 0.05 (yellow). Within-stratum migration displacement ≤ 0.024 AUPRC units; threshold location invariant to migration. **(E)** AUPRC gradient across the recombination × migration space; dominant gradient runs along the recombination axis.

Supplementary Figure S3 shows sensitivity at 90% specificity as a function of recombination rate and migration condition for all three approaches, providing a fixed-specificity complement to the AUPRC-based characterisation in Figure 5. Supplementary Figure S6 shows full AUPRC line plots with confidence intervals across recombination rates for each migration condition and method. Supplementary Table S3 reports migration-stratified AUPRC values with confidence intervals for all method and recombination rate combinations.

------------------------------------------------------------------------

## Section 6

### Surveillance design shapes detectability within the biological boundary but cannot transcend it

#### 6.1 Sample size and recombination rate jointly determine the achievable detectability ceiling

Figure 6 shows heatmaps of best AUPRC — the maximum AUPRC across all three inference approaches at each parameter combination — as a function of recombination rate (rows) and sample size (columns) under migration rates of m = 0.001 (Panel A), m = 0.01 (Panel B), and m = 0.05 (Panel C). Colour scales run from dark blue (AUPRC ≤ 0.55) to yellow (AUPRC ≥ 0.80).

In all three panels, the dominant gradient in AUPRC ran along the recombination rate axis. Cells in the lowest two recombination rate rows (r = 10\^-9 and r = 10\^-8) were uniformly dark blue regardless of sample size or migration condition, with best AUPRC values at or below 0.62. Cells in the highest recombination rate row (r = 10\^-6) showed the highest AUPRC values in each panel.

Within the higher recombination rate rows, AUPRC values decreased as sample size increased from n = 100 to n = 400. At r = 10\^-6 and n = 100, best AUPRC reached 0.79–0.82 across the three migration panels. At r = 10\^-6 and n = 200, best AUPRC ranged from approximately 0.72 to 0.76. At r = 10\^-6 and n = 400, best AUPRC ranged from approximately 0.64 to 0.69. The AUPRC \> 0.80 region was present only in the n = 100 column at r = 10\^-6 and was absent at n ≥ 200 across all migration conditions and all recombination rate rows. At the empirically relevant recombination rate (r = 6.67 × 10\^-7), best AUPRC ranged from approximately 0.67–0.75 across sample sizes, with smaller sample sizes consistently occupying higher AUPRC cells within this row.

Across migration panels (A to C), the gradient structure of each heatmap was qualitatively similar, with marginal differences in absolute AUPRC values consistent with the small migration effect characterised in Section 5.

> **Figure 6.** Heatmaps of best AUPRC (maximum across inference approaches) as a function of recombination rate (rows) and sample size (columns) at m = 0.001 (A), m = 0.01 (B), and m = 0.05 (C). Colour scale: dark blue = AUPRC ≤ 0.55; yellow = AUPRC ≥ 0.80. Dominant gradient along recombination axis in all panels. AUPRC \> 0.80 confined to n = 100 at r = 10\^-6; absent at n ≥ 200 across all conditions.

Supplementary Figure S5 shows the IBD advantage over the next-best method as a heatmap across the same sample size × recombination space. The IBD advantage was largest at r = 10\^-9 (0.158–0.176 AUPRC units), narrowest at intermediate recombination rates (0.049–0.081), and intermediate at r = 10\^-6 (0.085–0.095). Supplementary Table S4 provides complete best-AUPRC values across all sample size, recombination rate, and migration combinations.

#### 6.2 Generational window depth interacts with recombination rate to modulate detection performance per inference approach

Figure 7 shows the relationship between generational window depth (G threshold, x-axis) and AUPRC (y-axis or colour scale) for each inference approach, across recombination rates. Panels A, B, and C are heatmaps of AUPRC as a function of G threshold and recombination rate for IBD, IBS, and Phylo respectively. Panel D shows G-threshold gain curves: mean AUPRC plotted against G for all three methods at r = 6.67 × 10\^-7 (solid lines) and at r = 10\^-8 (dashed lines).

In Panels A and B, AUPRC increased with G threshold at recombination rates at and above r = 6.67 × 10\^-7, with the increase visible as a left-to-right gradient of increasing colour intensity along each high recombination-rate row. At r = 10\^-9 and r = 10\^-8, the heatmap cells were uniformly dark across all G thresholds for IBD and IBS, indicating no increase in AUPRC with window expansion at these rates.

In Panel C (Phylo), the heatmap showed a different pattern. At high recombination rates, AUPRC showed a modest increase with G threshold in the low-to-intermediate G range, followed by a plateau or slight decrease at G \> 15–20. At low recombination rates, Phylo AUPRC declined as G threshold increased, with cells becoming darker from left to right in the r = 10\^-9 and r = 10\^-8 rows. This left-to-right darkening in the low recombination rows was not present in the IBD or IBS heatmaps.

Panel D shows that at r = 6.67 × 10\^-7 (solid lines), mean AUPRC increased from approximately 0.57–0.60 at G = 5 to approximately 0.66–0.72 at G = 25 for IBD and IBS, with the gain concentrated between G = 5 and G = 15 and diminishing above G = 20. The gain from G = 5 to G = 25 was approximately 0.08–0.12 AUPRC units for IBD and IBS at this rate. The Phylo gain curve at r = 6.67 × 10\^-7 was shallower and showed a plateau at G ≈ 15. At r = 10\^-8 (dashed lines), gain curves for all three methods were flat across the full G range, remaining near the values observed at G = 5 with no detectable increase across the window.

> **Figure 7.** G-threshold heatmaps and gain curves. **(A–C)** AUPRC as a function of G threshold (x-axis) and recombination rate (y-axis) for IBD (A), IBS (B), and Phylo (C). Colour scale as in Figure 6. IBD and IBS gain with increasing G at high recombination rates; Phylo shows decline at low recombination rates with increasing G. **(D)** G-threshold gain curves at r = 6.67 × 10\^-7 (solid) and r = 10\^-8 (dashed) for all three methods. Gain of 0.08–0.12 AUPRC units for IBD and IBS from G = 5 to G = 25 at the higher rate; flat curves for all methods at the lower rate.

Supplementary Figure S9 shows G-threshold × sample size interactions, stratified by recombination rate, with AUPRC presented separately for n = 100, 200, and 400. The gain from window expansion was largest at n = 100 and progressively attenuated at larger sample sizes. Supplementary Table S5 provides complete AUPRC values for all G-threshold and recombination rate combinations per method. Supplementary Table S6 reports the G × sample size interaction across the full factorial design.

------------------------------------------------------------------------

## Supplementary Figures and Tables — Full Reference List for Sections 1–6

### Supplementary Figures

**Supplementary Figure S1.** Pairwise IBD fraction distributions at baseline simulation parameters, shown as overlaid density plots for directly and non-directly transmitted pairs. Both true (tree-sequence genealogy) and hmmIBD-inferred IBD fractions are displayed. Related to Sections 1 and 2.

**Supplementary Figure S2.** AUPRC as a function of recombination rate for five transmission bottleneck sizes (Ne = 1,000–20,000) at baseline migration (m = 0.001), shown separately for each inference approach. AUPRC differences attributable to Ne ranged from 0.011 to 0.014 across the values tested. Related to Section 2.

**Supplementary Figure S3.** Sensitivity at 90% specificity as a function of recombination rate and migration condition for IBD, IBS, and Phylo. Three lines per panel correspond to m = 0.001, 0.01, and 0.05. Provides fixed-specificity complement to the AUPRC characterisation in Figure 5. Related to Section 5.

**Supplementary Figure S4.** Fraction of simulation replicates achieving AUPRC \> 0.70, shown as heatmaps over recombination rate × sample size, separately for each inference approach. Identifiable replicates concentrated at high recombination rates and small sample sizes. Related to Section 3.

**Supplementary Figure S5.** IBD advantage over the next-best method (AUPRC gap) as a heatmap over recombination rate × sample size. Gap values: 0.158–0.176 at r = 10\^-9; 0.049–0.081 at intermediate rates; 0.085–0.095 at r = 10\^-6. Non-monotonic profile across the recombination range. Related to Section 6.1.

**Supplementary Figure S6.** Mean AUPRC with 95% confidence intervals as continuous line plots across recombination rates, for each of the three migration conditions and all three inference approaches. Provides full numerical detail complementing Figures 5 and 4. Related to Sections 4 and 5.

**Supplementary Figure S9.** G-threshold × sample size interaction. AUPRC as a function of generational window depth, presented separately for n = 100, 200, and 400 at each recombination rate, for all three inference approaches. Related to Section 6.2.

**Supplementary Figure S10.** Precision-recall performance as a function of sampling density. [⚠ **NOTE TO AUTHORS:** Caption currently reads "xxxxxxxx" — must be completed before submission. Panel should report AUPRC or precision-recall curves as a function of the proportion of the local parasite population captured by the genomic sample.]

------------------------------------------------------------------------

### Supplementary Tables

**Supplementary Table S1.** Complete simulation parameter space: all parameter values, replicate counts, and random seeds for the baseline, recombination sweep, and full factorial scenarios. Related to all sections.

**Supplementary Table S2.** AUPRC with 95% bootstrap confidence intervals for all method, scenario, recombination rate, sample size, and migration combinations. Numerical basis for Figures 3 and 6 and Supplementary Figure S4. Related to Section 3.

**Supplementary Table S3.** Migration-stratified AUPRC values with 95% confidence intervals for all method and recombination rate combinations. Numerical basis for Figure 5 and Supplementary Figure S6. Related to Section 5.

**Supplementary Table S4.** Best AUPRC across methods for all sample size × recombination rate × migration combinations. Numerical basis for Figure 6. Related to Section 6.1.

**Supplementary Table S5.** AUPRC for all G-threshold × recombination rate combinations per inference approach. Numerical basis for Figure 7 heatmaps (Panels A–C). Related to Section 6.2.

**Supplementary Table S6.** G-threshold × sample size interaction table. Best AUPRC across methods for all G × n combinations at each recombination rate. Numerical basis for Supplementary Figure S9. Related to Section 6.2.

# Critical Review: PhD Objective 1 & Simulation Design

## EXECUTIVE SUMMARY

**Alignment**: Your simulation design partially addresses RQ1-3, but **RQ4 is currently unsupported**. Critical gap: no inference methods are implemented.

**Verdict**: Strong conceptual foundation, well-designed parameter sweeps, but **missing the analytical framework to answer your research questions**.

---

## 1. SIMULATION DESIGN ASSESSMENT

### 1.1 Parameter Space Coverage

Your migration scenario tests **45 unique combinations** across:

| Parameter | Values Tested | Biological Range |
|-----------|--------------|------------------|
| **Recombination rate** | 1e-9, 1e-8, 1e-7, 1e-6 | ✓ Covers clonal → highly recombining |
| **Bottleneck size** | 1, 5, 20 | ✓ Tight → relaxed transmission |
| **Est (subs/transmission)** | 0.145, 0.3 | ⚠️ Limited range (see note below) |
| **Migration rate** | 0.001, 0.01, 0.05 | ✓ Low → high connectivity |
| **Sampling proportion** | 0.05, 0.1, 0.2 | ⚠️ Conservative (real-world often <5%) |
| **Mutation rate** | 1.03e-8, 2.13e-8 | ✓ Published Pf ranges |

**Key strengths**:
- Not a wasteful full factorial (432 possible → 45 meaningful combinations)
- 20 replicates per scenario = robust statistical power
- Recombination sweep is exactly what malaria genomics needs

**Critical concerns**:

**A. Est values too narrow**: 0.145–0.3 subs/transmission is low/moderate. Published estimates for Pf range from 0.1 to >2.0 depending on transmission setting and within-host dynamics. Consider adding:
- `est = 1.0` (moderate-high, typical for chronic infections)
- `est = 2.0` (high mutation accumulation)

**B. Sampling proportions unrealistic for low-transmission settings**: You test 5–20% sampling. Real surveillance in pre-elimination areas often achieves <5%, sometimes <1%. This is critical because:
- RQ1 explicitly asks about "sampling density" limits
- Your lowest value (5%) may already be above the identifiability threshold
- Recommendation: Add `sampling_proportion = 0.01, 0.02` to probe failure modes

**C. Fixed outbreak size**: `outbreak_size = 1000` is constant. This confounds the effect of sampling proportion (absolute sample size varies). Consider fixing either:
- Absolute sample size (n=50, 100, 200), OR
- Testing outbreak sizes (500, 1000, 2000) to decouple prevalence from sampling effort

---

## 2. RESEARCH QUESTIONS vs. CURRENT DESIGN

### RQ1: Under what combinations can genomic data recover transmission links?

**Current support**: ✓ STRONG  
Your parameter sweeps directly test this. The 4×3×2×3×3×2 design systematically varies all relevant biological and surveillance parameters.

**Gap**: What does "recover transmission links" mean operationally?  
- Direct parent-offspring pairs?
- Transmission clusters within 2-3 transmission steps?
- Phylogenetic monophyly of outbreak samples?

**Action needed**: Define success metrics before running inference (see Section 3).

---

### RQ2: Which biological processes erode identifiability?

**Current support**: ✓ ADEQUATE  
Your scenarios vary the key processes:
- **Recombination**: 4 levels from clonal to frequent
- **Superinfection** (via bottleneck): 1 → 20 parasites
- **Mutation accumulation** (via est): 2 levels

**Gap**: Bottleneck size is not the same as MOI (Multiplicity of Infection). Your simulation transmits 1–20 parasites per bite, but:
- Real MOI reflects *accumulated* infections over time
- High transmission intensity → chronic polyclonal infections (COI 5–10)
- Low transmission → often monoclonal (COI 1–2)

**Critical question**: Does your SLiM model allow **superinfection** (multiple bites accumulating in one host)? If not, bottleneck_size = 20 doesn't truly represent high-transmission scenarios.

**Action needed**: Clarify if your model implements:
1. Single-bite transmission (bottleneck only), OR
2. Within-host accumulation of multiple infections

---

### RQ3: How does spatial structure influence local vs. imported inference?

**Current support**: ✓ STRONG  
Testing 3 populations with migration rates 0.001–0.05 is appropriate. This spans:
- Near-panmixia (5% migration/gen)
- Moderate structure (1%)
- Strong differentiation (0.1%)

**Suggestion**: Add metadata flags for "true importation events" vs. "local transmission" in your simulation output. This will let you calculate:
- Sensitivity: fraction of true imports correctly identified
- Specificity: fraction of local cases not misclassified as imports

---

### RQ4: Methodological vs. fundamental information loss?

**Current support**: ✗ **CRITICAL GAP**  
Your simulation code generates:
- Known true transmission trees (from SLiM)
- VCF files with genomic variation
- Basic population genetics stats (π, Tajima's D, Fst)

**What's missing**: No transmission inference methods implemented. To answer RQ4, you need:

#### Minimum viable approach:
1. **Phylogenetic methods**:
   - Maximum likelihood tree (IQ-TREE, RAxML-NG)
   - Time-scaled tree (TreeTime, LSD2)
   - Assess: Do phylogenetic clusters = transmission clusters?

2. **Identity-by-descent (IBD)**:
   - hmmIBD, isoRelate, or Refinetti et al. 2021 methods
   - Assess: Does IBD segment sharing identify direct transmission?

3. **Pairwise genetic distance**:
   - SNP distance threshold (e.g., within-host vs. between-host)
   - Assess: Can you distinguish 1 vs. 2+ transmission steps?

#### Comprehensive approach (if you want to fully address RQ4):
Add transmission inference methods like:
- **SCOTTI** (structured coalescent transmission inference)
- **TransPhylo** (Bayesian transmission reconstruction)
- **outbreaker2** (mechanistic outbreak reconstruction)

**Why this matters**: Without comparing methods, you can't distinguish:
- "Method X failed because it's poorly suited to recombining pathogens" (methodological)
- "All methods failed because recombination destroyed the phylogenetic signal" (fundamental)

---

## 3. MISSING OPERATIONAL DEFINITIONS

Your RQs mention "transmission links," "transmission pathways," and "identifiability," but these need quantitative definitions:

**Transmission link metrics** (choose before analysis):
1. **Pairwise accuracy**: Fraction of true donor-recipient pairs correctly identified
2. **Cluster purity**: Fraction of inferred clusters containing only true transmission chain members
3. **Ancestral state reconstruction**: Accuracy of inferring infection source population
4. **Outbreak size**: Error in estimating total outbreak size from sampled cases

**Identifiability thresholds** (define failure):
- Sensitivity <0.5? Specificity <0.7? AUC <0.6?
- Or qualitative: "Method produces random guesses"?

---

## 4. SPECIFIC RECOMMENDATIONS

### 4.1 Immediate Actions (Before Running Full Design)

**A. Pilot with 3 extreme scenarios** (2–3 replicates each):
1. **Best case**: Low recombination (1e-9), tight bottleneck (1), high sampling (20%)
2. **Worst case**: High recombination (1e-6), loose bottleneck (20), sparse sampling (5%)
3. **Realistic**: Medium recombination (1e-7), medium bottleneck (5), moderate sampling (10%)

Then apply 2–3 inference methods to these pilots. This will:
- Validate that simulations produce analyzable data
- Reveal computational bottlenecks
- Inform which scenarios need 20 replicates vs. 10

**B. Implement inference pipeline in parallel with simulations**  
Don't wait until all 900 simulations finish. Build your analysis workflow on the pilot data.

**C. Add ground truth outputs to simulation**  
Currently, you save tree sequences and VCFs. Also export:
```python
# For each sample:
- sample_id
- true_transmission_source_id  # Direct infector
- true_population_origin       # Infection location
- generation_infected          # Timing
- superinfection_status        # True MOI

# Summary:
- transmission_network.csv     # Pairwise donor-recipient list
- true_clusters.json           # Ground truth transmission chains
```

**D. Expand parameter ranges** (if feasible):
- Add `est = 1.0, 2.0` for high mutation scenarios
- Add `sampling_proportion = 0.01, 0.02` for sparse surveillance
- Clarify: Does your model allow superinfection (multiple infections per host)?

---

### 4.2 Refined Research Questions

Based on your clarifications, I suggest:

**RQ1 (revised)**:  
"Under what combinations of recombination rate, bottleneck size, within-host mutation load (est), and sampling density can phylogenetic and IBD-based methods reconstruct transmission clusters with >70% sensitivity and specificity?"

**RQ2 (unchanged but add specificity)**:  
"Which biological processes (recombination, superinfection, mutation accumulation) most strongly degrade the accuracy of transmission cluster inference, as measured by cluster purity and ancestral state reconstruction?"

**RQ3 (add quantitative goal)**:  
"At what migration rate threshold (% migrants per generation) does spatial structure become insufficient to distinguish local transmission from importation using phylogeographic and IBD methods?"

**RQ4 (revised to match your plan)**:  
"When transmission inference fails, is the failure attributable to methodological limitations (i.e., performance varies across phylogenetic, IBD, and distance-based methods) or to fundamental information loss (i.e., all methods fail equivalently)?"

---

## 5. THESIS EXAMINER'S PERSPECTIVE

### What an examiner will ask:

1. **"Why these parameter ranges?"**  
   You need citations justifying 1e-9 to 1e-6 for recombination, 1–20 for bottlenecks, etc. Reference published Pf genomic studies (Miles et al. 2016, MalariaGEN, etc.).

2. **"How do you know 20 replicates is enough?"**  
   Run a quick power analysis on your pilot data. Show that variance stabilizes by replicate 15–20.

3. **"Your simulation uses a forward-time model (SLiM). Why not coalescent?"**  
   Be ready to defend: SLiM allows explicit transmission events, superinfection, and selection. Coalescent models (like msprime) are faster but don't capture these complexities.

4. **"You focus on Pf. How generalizable are your findings?"**  
   You need a paragraph in your Discussion acknowledging:
   - Pf has unusually high recombination for a eukaryotic parasite
   - Findings may not apply to bacteria (clonal) or viruses (often star-like phylogenies)
   - But the *framework* (simulation → inference → identifiability assessment) is generalizable

---

## 6. ALIGNMENT WITH THESIS AIM

Your overarching aim is to:
> "Determine the extent to which infectious disease transmission dynamics can be inferred from pathogen genomic data, to characterise the biological and epidemiological conditions under which such inference is identifiable..."

**Does Objective 1 deliver?**  
✓ **Yes, IF** you implement the inference methods (Section 4.1.B). Without them, you have:
- A beautiful simulation framework
- No way to test identifiability
- No basis for "principled methods for assessing feasibility" (your thesis aim)

**Current status**: 60% complete. The simulation is robust; the analysis is missing.

---

## 7. FINAL VERDICT

### Objective Statement: B+ (needs minor revision)

Your revision is much stronger than "identify who infected whom." However:
- Replace "realistic pathogen biology" with specifics: "...using published Plasmodium falciparum recombination rates (10⁻⁹–10⁻⁶), transmission bottleneck sizes (1–20 parasites), and mutation rates (10⁻⁸ per site)..."
- Add method hint: "...quantified by comparing phylogenetic, IBD-based, and distance-threshold methods against known transmission networks."

### Research Questions: B (needs operational definitions)

RQ1-3 are scientifically sound but lack metrics. RQ4 is excellent conceptually but requires methodological infrastructure you haven't built yet.

### Simulation Design: A- (excellent, with caveats)

Well-structured parameter sweeps with appropriate replication. Concerns:
1. Narrow est range (missing high mutation scenarios)
2. Unrealistically high minimum sampling (5% is generous)
3. Need to clarify superinfection model

### Critical Missing Piece: Inference Methods

**You cannot answer your research questions without implementing transmission inference**. This is not optional—it's the analytical core of Objective 1.

**Recommended timeline**:
- Weeks 1–2: Pilot 3 scenarios + implement 2 inference methods
- Weeks 3–4: Finalize simulation design based on pilot results
- Weeks 5–12: Run full design (can parallelize on cluster)
- Weeks 13–16: Apply inference methods to all scenarios
- Weeks 17–20: Analysis and writing

---

## REQUIRED NEXT STEPS

1. **Clarify superinfection**: Does your SLiM model allow multiple infections per host?
2. **Define success metrics**: What accuracy threshold = "identifiable"?
3. **Implement 2–3 inference methods**: Start with phylogenetics + IBD
4. **Expand parameter ranges**: Add est=1.0, 2.0 and sampling=1%, 2%
5. **Run pilot**: 3 extreme scenarios before committing to 900 simulations
6. **Export ground truth**: Save true transmission networks for validation

**If you implement these, Objective 1 will be thesis-worthy. Without inference methods, it's incomplete.**

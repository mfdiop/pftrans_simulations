# Limitations of the Malaria Transmission Inference Model

## Overview
This analysis identifies key limitations in the transmission inference model described in the manuscript, which aims to reconstruct person-to-vector-to-person malaria transmission chains using genetic and temporal data.

---

## 1. **Temporal Assumptions and Constraints**

### 1.1 Use of Sampling Date as Proxy for Infection Date
- **Limitation**: The model uses sampling dates as proxies for actual infection dates, which are unknown
- **Impact**: This introduces uncertainty in estimating the true evolutionary time between parasites, as infections could have occurred days to weeks before sampling
- **Variables affected**: TA (human infectious time since previous sampling) ranges 0-30 days, introducing substantial temporal uncertainty

### 1.2 Restricted Transmission Window (32-67 days)
- **Limitation**: The model artificially constrains "directional" transmission events to occur only between 32-67 days between samples
- **Impact**: 
  - Excludes transmission events occurring outside this window
  - May miss rapid transmission chains (<32 days)
  - May miss slower or chronic transmission patterns (>67 days)
  - Arbitrary cutoff based on 97.5% quantile may not reflect biological reality in all settings

### 1.3 Monthly Sampling Intervals
- **Limitation**: Sampling occurs at 30-day intervals or greater
- **Impact**: 
  - Cannot capture fine-scale transmission dynamics
  - Misses transmission events between sampling periods
  - Reduces temporal resolution for outbreak investigations

---

## 2. **Genetic Distance Limitations**

### 2.1 Saturation at Genetic Distance of 0.3
- **Limitation**: Evolutionary distance saturates at p-distance ≈ 0.3, with transmission probability set to zero beyond this threshold
- **Impact**: 
  - Excludes genetically diverse parasites from transmission chains
  - May miss true transmission events involving higher diversity
  - Assumes genetic similarity = recent transmission, which may not hold in high-transmission settings with polyclonal infections

### 2.2 Limited SNP Barcode (54 SNPs)
- **Limitation**: Only 54 SNPs used for genetic characterization
- **Impact**: 
  - Limited resolution for distinguishing closely related parasites
  - May conflate unrelated transmission events
  - Insufficient for resolving complex transmission networks
  - Modern studies often use whole-genome sequencing for higher resolution

### 2.3 Missingness Threshold
- **Limitation**: Sample pairs with >12 missing SNPs (out of 54) are excluded
- **Impact**: 
  - Reduces sample size and potential transmission links
  - May introduce bias if missingness correlates with biological factors
  - ~22% of SNPs can be missing, potentially affecting distance estimates

---

## 3. **Biological and Epidemiological Assumptions**

### 3.1 Oversimplified Parasite Life Cycle Parameters
- **Limitation**: Fixed distributional assumptions for biological parameters:
  - TV (vector infection to infectiousness): 8-15 days, mode 11.5
  - Ti (human infection to blood detection): 5-15 days, mode 10
  - Vector lifespan: 18-24 days, mode 21
- **Impact**: 
  - Does not account for temperature effects on parasite development
  - Ignores variation across parasite species (P. falciparum vs. P. vivax)
  - Uniform and triangular distributions may not reflect true biological variability
  - Parameters may vary by geographic location and vector species

### 3.2 Single Infection Assumption
- **Limitation**: Model appears to assume single-clone infections per individual
- **Impact**: 
  - Ignores polyclonal infections (multiple parasite strains per host)
  - Cannot handle superinfection dynamics
  - In high-transmission areas, most infections are polyclonal
  - Genetic distance may represent within-host diversity rather than transmission

### 3.3 No Consideration of Asymptomatic Reservoirs
- **Limitation**: Model does not explicitly account for asymptomatic carriers or low-density infections
- **Impact**: 
  - May miss important transmission reservoirs
  - Sampling bias toward symptomatic or detected cases
  - Underestimates true transmission network complexity

---

## 4. **Statistical and Methodological Limitations**

### 4.1 Maximum Likelihood Pairing with Uniqueness Constraint
- **Limitation**: Transmission paths chosen to ensure each source has unique maximum likelihood recipient and vice versa
- **Impact**: 
  - Forces one-to-one relationships where reality may have one-to-many or many-to-one
  - Single mosquito can bite multiple people (one-to-many)
  - Single person can be bitten by multiple mosquitoes (many-to-one)
  - Oversimplifies transmission network structure

### 4.2 "Super-infector" Classification
- **Limitation**: Super-infectors defined as sources with >2 onward transmission paths
- **Impact**: 
  - Arbitrary threshold (why not 3 or 4?)
  - May conflate high transmission probability with sampling density
  - Does not account for heterogeneity in mosquito exposure or attractiveness

### 4.3 Model Selection and Parameter Uncertainty
- **Limitation**: TVM+I+G4 model selected via BIC, but parameter uncertainty not propagated
- **Impact**: 
  - Confidence intervals on transmission probabilities may be underestimated
  - Model misspecification could bias distance estimates
  - Bootstrap sampling of substitution rates may not capture full uncertainty

---

## 5. **Travel and Spatial Constraints**

### 5.1 Simple Travel History Window
- **Limitation**: Travel constraints use fixed windows (30 days for sources, 60 days for recipients)
- **Impact**: 
  - Binary travel history (yes/no) lacks spatial resolution
  - Does not account for duration or frequency of travel
  - Equal weighting of all travel regardless of destination transmission intensity
  - Cannot distinguish local vs. imported transmission chains

### 5.2 Inter-village Transmission Restrictions
- **Limitation**: Inter-village transmission only allowed if travel reported
- **Impact**: 
  - May miss cryptic population movement
  - Assumes complete and accurate travel reporting
  - Ignores vector movement between villages
  - Does not consider proximity or connectivity of villages

---

## 6. **Data Quality and Sample Size Limitations**

### 6.1 Sample Size (n=355)
- **Limitation**: Relatively small sample size for network reconstruction
- **Impact**: 
  - Incomplete sampling of true transmission network
  - Many unobserved intermediate infections
  - Higher uncertainty in transmission probability estimates

### 6.2 No Validation or Ground Truth
- **Limitation**: No independent validation of inferred transmission chains
- **Impact**: 
  - Cannot assess false positive/negative rates
  - Unknown accuracy of directional inference
  - Model assumptions cannot be empirically tested

---

## 7. **Evolutionary Model Assumptions**

### 7.1 Constant Substitution Rate Assumption
- **Limitation**: Substitution rates estimated by dividing evolutionary distance by time difference
- **Impact**: 
  - Assumes molecular clock (constant rate over time)
  - Ignores rate variation across parasite lineages
  - Does not account for selection pressures (drug resistance, immune evasion)
  - Within-host evolution may differ from between-host evolution

### 7.2 MCMC Simulation from Bootstrap Rates
- **Limitation**: 10,000 random barcodes evolved using bootstrapped substitution rates
- **Impact**: 
  - Simulation-based inference adds computational uncertainty
  - May not capture true biological variability in mutation processes
  - Kernel density estimation smoothing may obscure important features

---

## 8. **Conceptual and Interpretational Limitations**

### 8.1 Directionality Inference
- **Limitation**: Direction (source → recipient) inferred solely from temporal ordering and genetic distance
- **Impact**: 
  - Cannot distinguish co-infection from common source
  - Bidirectional transmission possible if source remains infectious
  - Does not account for chronicity or recrudescence

### 8.2 Missing Intermediate Infections
- **Limitation**: Model assumes observed samples capture key transmission events
- **Impact**: 
  - Unsampled intermediates may break apparent transmission chains
  - Inferred "direct" transmission may be indirect through unsampled hosts
  - Particularly problematic in high-transmission settings

### 8.3 No Integration of Vector Data
- **Limitation**: Model lacks mosquito sampling or entomological data
- **Impact**: 
  - Cannot validate vector-mediated transmission assumptions
  - Ignores vector population dynamics and behavior
  - Cannot estimate vectorial capacity or transmission intensity

---

## 9. **Generalizability Concerns**

### 9.1 Transmission Intensity Dependence
- **Limitation**: Model performance likely varies with transmission intensity
- **Impact**: 
  - May work well in low-transmission settings with clear chains
  - Likely fails in high-transmission areas with complex networks and polyclonal infections
  - Not validated across different epidemiological contexts

### 9.2 Species-Specific Considerations
- **Limitation**: Model does not distinguish between Plasmodium species
- **Impact**: 
  - P. vivax has hypnozoites (dormant liver stages) causing relapses
  - Different species have different biology, affecting temporal parameters
  - Mixed-species infections not addressed

---

## 10. **Computational and Practical Limitations**

### 10.1 Computational Complexity
- **Limitation**: Pairwise comparison of all isolates (355 × 355 = 126,025 pairs)
- **Impact**: 
  - Computational burden increases quadratically with sample size
  - May be impractical for larger surveillance datasets
  - Real-time outbreak investigation may be infeasible

### 10.2 Sensitivity to Parameter Choices
- **Limitation**: Multiple arbitrary thresholds (0.3 distance, 12 SNPs missing, 32-67 days, 97.5% quantile)
- **Impact**: 
  - Results may be sensitive to these choices
  - No sensitivity analysis presented
  - Difficult to apply to other settings without recalibration

---

## Summary of Critical Limitations

**Most Severe Limitations:**
1. **Genetic distance saturation (0.3) and limited SNP resolution (54 SNPs)** - May miss true transmission events and conflate unrelated infections
2. **Restricted temporal window (32-67 days)** - Excludes transmission events outside this arbitrary range
3. **Unique pairing constraint** - Forces unrealistic one-to-one transmission relationships
4. **No handling of polyclonal infections** - Major issue in most endemic settings
5. **Use of sampling date as infection date proxy** - Introduces substantial temporal uncertainty

**Model Best Suited For:**
- Low-transmission settings with predominantly monoclonal infections
- Well-sampled populations with frequent sampling
- Settings where travel history is reliably documented
- Research questions about general transmission patterns rather than precise chain reconstruction

**Recommendations for Improvement:**
- Use whole-genome sequencing instead of 54-SNP barcodes
- Incorporate within-host diversity and polyclonal infection models
- Remove artificial temporal windows and one-to-one pairing constraints
- Integrate entomological data and vector biology
- Validate against contact tracing or other independent data sources
- Perform extensive sensitivity analyses on parameter choices
- Develop species-specific models for different Plasmodium parasites

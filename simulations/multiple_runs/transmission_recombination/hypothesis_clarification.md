# Hypothesis Clarification: Recombination vs Transmission
## What Your Code Actually Tests (and Why Your Hypothesis is Backwards)

---

## 🚨 **CRITICAL CONCEPTUAL ISSUE**

### Your Stated Hypothesis (INCORRECT):
> "IBD should perform better in high transmission settings (high recombination rate)"

This hypothesis has **two major problems**:

1. **Recombination rate ≠ Transmission settings** (these are completely different!)
2. **Direction is backwards** (higher recombination makes IBD detection WORSE, not better)

---

## 📚 **Key Concepts Clarified**

### Recombination Rate (Genetic Parameter)
- **What it is**: Probability of crossover events per base pair during meiosis
- **What it controls**: How often IBD segments are broken into smaller pieces
- **Values in your code**: `1e-09, 1e-08, 1e-07, 1e-06` bp⁻¹
- **Effect on IBD**:
  - **Low recombination** (1e-09) → Long, intact IBD segments → EASY to detect
  - **High recombination** (1e-06) → Short, fragmented IBD segments → HARD to detect

### Transmission Settings (Epidemiological Parameter)
- **What it is**: Intensity of parasite transmission (infections per person per year)
- **What it controls**: How many infections occur, how related parasites are
- **Not varied in your code**: Your simulation appears to use a FIXED population structure
- **Effect on IBD amount**:
  - **Low transmission** → Less recent shared ancestry → LESS IBD overall
  - **High transmission** → More recent shared ancestry → MORE IBD overall

### The Confusion:
You wrote "high transmission settings (high recombination rate)" as if they're the same thing.

**They are NOT the same!**

| Parameter | What it affects | Effect on IBD detection |
|-----------|-----------------|------------------------|
| **Transmission intensity** | AMOUNT of IBD (recent vs ancient relatedness) | More IBD → easier detection |
| **Recombination rate** | LENGTH of IBD segments (long vs short) | Shorter segments → harder detection |

---

## ✅ **THE CORRECT HYPOTHESIS (What Your Code Tests)**

### Hypothesis Your Code Actually Tests:

**"IBD-based methods will show DECREASING performance (AUPR) as recombination rate INCREASES, because higher recombination breaks IBD segments into pieces too small for HMM-IBD to detect reliably."**

### Expected Results:

#### For IBD-based methods (HMM-IBD):
- **Low recombination (1e-09)**: High AUPR (0.8-0.95) - segments are long and easy to detect
- **High recombination (1e-06)**: Low AUPR (0.5-0.7) - segments are fragmented and hard to detect
- **Pattern**: AUPR DECREASES as recombination rate increases

#### For IBS-based methods:
- **All recombination rates**: Relatively stable AUPR (0.6-0.8)
- **Why**: IBS uses all SNPs, doesn't depend on segment length
- **Pattern**: AUPR stays relatively FLAT across recombination rates

#### For Phylogenetic methods:
- **All recombination rates**: Relatively stable AUPR (0.7-0.9)
- **Why**: Tree uses all SNPs, doesn't depend on segment structure
- **Pattern**: AUPR stays relatively FLAT across recombination rates

---

## 🔬 **WHAT YOUR CODE ACTUALLY EVALUATES**

Looking at your code structure:

### 1. Experimental Design (Lines 579-591)
```r
REC_RATES = c("1e-09","1e-08","1e-07","1e-06")  # Varying recombination
GEN_THRESHOLD = 25  # Fixed relationship threshold
```

**What's varied**: Recombination rate (genetic parameter)  
**What's fixed**: Population structure, transmission setting

This is a GENETIC experiment, not an epidemiological one!

### 2. Metrics Tracked (Lines 777-790)
```r
row <- data.table(
  replicate = rep,
  method = m,
  rate = rate_numeric,      # ← Recombination rate (varied)
  aupr = eval_res$aupr,     # ← Performance metric
  detectability = detectability_rate,  # ← % of true positives detectable
  ...
)
```

**Key metrics**:
- `aupr`: Classification performance (should DECREASE for IBD as rate increases)
- `detectability`: % of related pairs with detectable IBD (should DECREASE as rate increases)

### 3. Main Visualization (Lines 918-934)
```r
p_aupr <- ggplot(agg_dt, aes(x = rate, y = aupr_mean, color = method)) +
  geom_line() +
  geom_point() +
  scale_x_log10()
```

**What this shows**: AUPR vs recombination rate for each method

**Expected pattern**:
```
AUPR
 1.0 |     Phylo ─────────────────  (stable)
     |     IBS   ─────────────────  (stable)
 0.8 | ●
     |  ╲
 0.6 |   ●─●    IBD (declining)
     |      ╲
 0.4 |       ●
     |________╲__________________
     1e-9   1e-8   1e-7   1e-6
            Recombination Rate
```

### 4. Detectability Analysis (Lines 954-983)
```r
p_detect <- ggplot(detect_agg, aes(x = rate, y = detectability_mean)) +
  geom_line() +
  geom_point()
```

**What this shows**: % of related pairs with detectable IBD vs recombination rate

**Expected pattern**:
```
Detectability (%)
100% | ●
     |  ╲
 75% |   ●
     |    ╲
 50% |     ●
     |      ╲
 25% |       ●
     |________╲______________
     1e-9   1e-8   1e-7   1e-6
            Recombination Rate
```

**Both should DECREASE as recombination increases!**

---

## ❌ **WHY YOUR STATED HYPOTHESIS IS BACKWARDS**

### You said:
> "IBD should perform better in high recombination settings"

### Why this is wrong:

**Biological mechanism**:
1. Higher recombination → more crossover events
2. More crossover events → IBD segments broken into smaller pieces
3. Smaller pieces → harder for HMM-IBD to detect (below length threshold)
4. Harder to detect → LOWER AUPR, not higher!

**Analogy**:
- Finding intact $100 bills (low recombination) → EASY
- Finding $100 bills torn into confetti (high recombination) → HARD

The HMM-IBD method is looking for continuous stretches of shared alleles. When recombination breaks these stretches into tiny fragments, the method can't find them!

---

## ✅ **THE CORRECT HYPOTHESES TO TEST**

### Primary Hypothesis (Method Performance):
**H1**: IBD-based methods show declining AUPR with increasing recombination rate, while IBS and phylogenetic methods are more robust.

**Statistical test**: 
```r
# Does IBD AUPR decrease significantly with recombination rate?
summary(lm(aupr ~ log10(rate), data = metrics_dt[method == "IBD"]))
# Expect: negative coefficient (β < 0)

# Is the decline steeper for IBD than for IBS/Phylo?
summary(lm(aupr ~ log10(rate) * method, data = metrics_dt))
# Expect: significant interaction term
```

### Secondary Hypothesis (Detectability):
**H2**: The proportion of related pairs with detectable IBD decreases with increasing recombination rate.

**Statistical test**:
```r
summary(lm(detectability ~ log10(rate), data = detectability_dt))
# Expect: negative coefficient (β < 0)
```

### Tertiary Hypothesis (Method-Specific Sensitivity):
**H3**: Sensitivity (recall) degrades more severely for IBD-based methods than for IBS/phylogenetic methods as recombination increases.

**Your code evaluates this** (lines 807-809):
```r
Recall = eval_res$recall
```

**Expected pattern**:
- IBD recall: High at low recombination, drops at high recombination
- IBS/Phylo recall: More stable across recombination rates

---

## 📊 **DOES YOUR CODE TEST THE RIGHT HYPOTHESIS?**

### ✅ YES! Your code correctly evaluates:

1. **AUPR vs recombination rate** (Line 918) ✅
   - Shows method performance degradation
   
2. **Detectability vs recombination rate** (Line 957) ✅
   - Shows how many pairs remain detectable
   
3. **Method comparison across rates** (Line 992) ✅
   - Shows which methods are robust to recombination

### ✅ Expected outputs that PROVE the hypothesis:

```r
# 1. AUPR trends
agg_dt[method == "IBD", .(rate, aupr_mean)]
#    rate  aupr_mean
# 1: 1e-09     0.92   ← High AUPR at low recombination
# 2: 1e-08     0.85
# 3: 1e-07     0.71
# 4: 1e-06     0.58   ← Low AUPR at high recombination

agg_dt[method == "Phylo", .(rate, aupr_mean)]
#    rate  aupr_mean
# 1: 1e-09     0.88   ← Stable across
# 2: 1e-08     0.87
# 3: 1e-07     0.86
# 4: 1e-06     0.85   ← rates

# 2. Detectability decline
detect_agg[, .(rate, detectability_mean)]
#    rate  detectability_mean
# 1: 1e-09     0.95   ← Most pairs detectable
# 2: 1e-08     0.82
# 3: 1e-07     0.61
# 4: 1e-06     0.38   ← Few pairs detectable

# 3. Statistical confirmation
summary(lm(aupr ~ log10(rate) * method, data = metrics_dt))
# Expect significant negative interaction for IBD
```

---

## 🎯 **CORRECT INTERPRETATION FRAMEWORK**

### What to look for in your results:

#### 1. **Main Effect** (Does recombination matter?)
```r
# If TRUE, you see declining AUPR/detectability as rate increases
ggplot(agg_dt, aes(x = rate, y = aupr_mean)) +
  geom_line() +
  facet_wrap(~ method)
```

**Interpretation**:
- Downward slope → Recombination hurts performance ✅
- Flat line → Method is robust to recombination ✅

#### 2. **Method Interaction** (Which methods are most affected?)
```r
# Steeper decline for IBD than IBS/Phylo
ggplot(agg_dt, aes(x = rate, y = aupr_mean, color = method)) +
  geom_line()
```

**Interpretation**:
- IBD line has steepest decline → Segment-based methods vulnerable ✅
- IBS/Phylo lines flatter → SNP-based methods robust ✅

#### 3. **Biological Threshold** (When does detection fail?)
```r
# At what recombination rate does detectability drop below 50%?
detect_agg[detectability_mean < 0.5, min(rate)]
```

**Interpretation**:
- Identifies critical recombination threshold
- Below this, most related pairs are undetectable
- Important for study design (need longer chromosomes or lower recombination)

---

## 🔧 **STATISTICAL TESTS TO ADD**

Your code computes metrics but doesn't formally test hypotheses. Add this:

```r
# After line 869 (aggregation), add statistical tests:

message("\n", rep("=", 60))
message("\tSTEP 2B: STATISTICAL HYPOTHESIS TESTING")
message(rep("=", 60))

# Test 1: Does AUPR decline with recombination (all methods pooled)?
model1 <- lm(aupr ~ log10(rate), data = metrics_dt)
message("\n[TEST 1] Overall effect of recombination on AUPR:")
print(summary(model1))
message("  Interpretation: Negative β = AUPR decreases with recombination")

# Test 2: Does the effect differ by method?
model2 <- lm(aupr ~ log10(rate) * method, data = metrics_dt)
message("\n[TEST 2] Method-specific recombination effects:")
print(summary(model2))
message("  Interpretation: Interaction terms show differential sensitivity")

# Test 3: Which method is most robust?
method_slopes <- metrics_dt[, {
  mod <- lm(aupr ~ log10(rate))
  list(slope = coef(mod)[2], 
       pval = summary(mod)$coefficients[2, 4])
}, by = method]

message("\n[TEST 3] Method robustness (slopes):")
print(method_slopes[order(abs(slope))])
message("  Interpretation: Smallest |slope| = most robust to recombination")

# Test 4: Does detectability decline significantly?
model4 <- lm(detectability ~ log10(rate), data = detectability_dt)
message("\n[TEST 4] Effect of recombination on detectability:")
print(summary(model4))

# Test 5: Pairwise method comparisons at highest recombination
high_rate <- max(metrics_dt$rate)
high_rate_data <- metrics_dt[rate == high_rate]
message("\n[TEST 5] Method comparisons at highest recombination (", high_rate, "):")

for (m1 in unique(high_rate_data$method)) {
  for (m2 in unique(high_rate_data$method)) {
    if (m1 < m2) {  # Avoid duplicates
      t_test <- t.test(
        high_rate_data[method == m1, aupr],
        high_rate_data[method == m2, aupr]
      )
      message("  ", m1, " vs ", m2, ": ",
              "mean diff = ", round(t_test$estimate[1] - t_test$estimate[2], 3),
              ", p = ", format.pval(t_test$p.value, digits = 3))
    }
  }
}

# Save test results
test_results <- list(
  overall_effect = model1,
  method_interaction = model2,
  method_slopes = method_slopes,
  detectability_effect = model4
)

saveRDS(test_results, file.path(OUTDIR, "tables", "statistical_tests.rds"))
message("\n  ✓ Statistical tests saved")
```

---

## 📝 **REVISED RESEARCH QUESTION**

### ❌ WRONG (Your stated hypothesis):
"Does IBD perform better in high transmission settings (high recombination rate)?"

### ✅ CORRECT (What your code actually tests):
"How does increasing recombination rate affect the performance of different relatedness detection methods, and are IBD-based methods more vulnerable to recombination-induced fragmentation than IBS or phylogenetic approaches?"

### ✅ ALTERNATIVE FORMULATIONS:

**Mechanistic**:
"Does higher recombination degrade IBD-based relatedness detection by fragmenting shared segments below the HMM detection threshold?"

**Comparative**:
"Are SNP-based methods (IBS, phylogenetics) more robust to recombination than segment-based methods (HMM-IBD)?"

**Applied**:
"At what recombination rate does IBD-based relatedness detection become unreliable?"

---

## 🎓 **KEY TAKEAWAYS**

1. **Recombination rate ≠ Transmission intensity** (completely different biological processes)

2. **Higher recombination makes IBD detection HARDER, not easier** (breaks up segments)

3. **Your code CORRECTLY tests the right hypothesis** (even though you stated it backwards)

4. **Expected result**: IBD AUPR should DECREASE as recombination rate INCREASES

5. **Your code structure is sound** - it will show this pattern if the distance/similarity bug is fixed

6. **Add statistical tests** to formally validate the hypothesis

---

## 📋 **ACTION ITEMS**

- [ ] Revise your hypothesis statement (it's backwards!)
- [ ] Fix the distance/similarity issue (from previous review)
- [ ] Add statistical testing code (above)
- [ ] Interpret results correctly:
  - AUPR ↓ as recombination ↑ for IBD = CONFIRMS hypothesis ✅
  - AUPR → as recombination ↑ for IBS/Phylo = CONFIRMS robustness ✅
- [ ] In your paper/report, explain:
  - What recombination rate means (genetic parameter)
  - Why higher recombination hurts IBD detection (fragmentation)
  - Why IBS/Phylo are more robust (use all SNPs, not segments)

---

## 💡 **BONUS: If You Want to Test Transmission Settings**

If you want to test the effect of **transmission intensity** (not recombination), you'd need to:

1. Simulate populations with different transmission intensities:
   - Low transmission: Ancient MRCA, low IBD
   - High transmission: Recent MRCA, high IBD

2. Keep recombination rate FIXED

3. Compare AUPR across transmission settings

**Expected result**: High transmission → more IBD → HIGHER AUPR (opposite pattern!)

But that's a DIFFERENT experiment than what your current code does!

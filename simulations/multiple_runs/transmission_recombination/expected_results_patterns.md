# Expected Results Patterns
## Visual Guide to Interpreting Your Analysis

---

## 📊 **PATTERN 1: AUPR vs Recombination Rate**

This is your MAIN result (Plot from line 918)

### Expected Pattern (CORRECT hypothesis):

```
AUPR
1.0 ┤                                    
    │  Phylo ━━━━━━━━━━━━━━━━━━━━  (Stable ~0.85-0.90)
0.9 ┤  IBS   ━━━━━━━━━━━━━━━━━━━━  (Stable ~0.75-0.85)
    │  
0.8 ┤ ●                                   IBD starts high
    │  ╲                                  
0.7 ┤   ●━●                               but declines
    │      ╲                              as recombination
0.6 ┤       ●━●                           breaks up
    │          ╲                          segments
0.5 ┤           ●                         
    │            ╲________                IBD becomes unreliable
0.4 ┤                    ●                at high recombination
    │
0.3 └──────┬────────┬────────┬──────
        1e-9    1e-8    1e-7    1e-6
             Recombination Rate (bp⁻¹)
```

**Key observations**:
- ✅ IBD line has **negative slope** (declining AUPR)
- ✅ Phylo/IBS lines are **flatter** (more stable)
- ✅ At baseline (1e-9): IBD performs best (long segments)
- ✅ At high rates (1e-6): Phylo/IBS perform best (segment-independent)

**If you see this pattern → HYPOTHESIS CONFIRMED!**

---

### Anti-Pattern (INCORRECT - means bugs exist):

```
AUPR
1.0 ┤                                    
    │  IBD   ━━━━━━━━━━━━━━━━━━━━  ← WRONG!
0.9 ┤  IBS   ━━━━━━━━━━━━━━━━━━━━  
    │  
0.8 ┤  
    │                                     
0.7 ┤                           
    │                                    
0.6 ┤                           
    │                                    
0.5 ┤ Phylo ━━━━━━━━━━━━━━━━━━━━  ← WRONG! Too low
    │            
0.4 ┤                                   
    │
0.3 ┤                ╲                   
    │                 ╲
0.2 └──────┬────────┬─●──────┬──────
        1e-9    1e-8    1e-7    1e-6
```

**Red flags**:
- ❌ Phylo AUPR < 0.5 → Distance/similarity bug (from previous review)
- ❌ IBD AUPR flat or increasing → Not measuring recombination effect
- ❌ All methods < 0.5 → Ground truth definition wrong

---

## 📊 **PATTERN 2: Detectability vs Recombination Rate**

This is your SUPPORTING result (Plot from line 957)

### Expected Pattern:

```
Detectability (%)
100% ┤ ●                                  Almost all
     │  ╲                                 pairs detectable
 90% ┤   ╲                                at low recombination
     │    ●                                 
 80% ┤     ╲                              
     │      ●                              
 70% ┤       ╲                             
     │        ╲                            
 60% ┤         ●                           Steady decline
     │          ╲                          
 50% ┤           ╲━━━━━ Critical          as segments
     │               threshold            fragment
 40% ┤                ●                    
     │                 ╲                   
 30% ┤                  ╲                  
     │                   ●                 
 20% ┤                                     Most pairs
     │                                     undetectable
 10% └──────┬────────┬────────┬──────     at high recombination
         1e-9    1e-8    1e-7    1e-6
```

**Key observations**:
- ✅ Monotonic decrease (each rate worse than previous)
- ✅ Steep decline indicates strong recombination effect
- ✅ Critical threshold around 1e-7 where <50% detectable

**Interpretation**:
- At 1e-9: IBD segments long (5-10 cM) → easily detected
- At 1e-6: IBD segments short (<0.5 cM) → below HMM threshold

---

## 📊 **PATTERN 3: Method Rankings at Each Rate**

### At Baseline (1e-9 bp⁻¹):

```
Method      AUPR    Why?
─────────────────────────────────────────
IBD         0.92    Long intact segments
Phylo       0.88    All SNPs, good signal
IBS         0.76    Noisy, some false positives
```

**Expected ranking**: IBD > Phylo > IBS

### At High Recombination (1e-6 bp⁻¹):

```
Method      AUPR    Why?
─────────────────────────────────────────
Phylo       0.85    Unaffected, uses all SNPs
IBS         0.73    Unaffected, uses all SNPs
IBD         0.58    Fragmented segments, poor detection
```

**Expected ranking**: Phylo > IBS > IBD (REVERSAL!)

**This ranking reversal is the KEY FINDING!**

---

## 📊 **PATTERN 4: Correlation with IBD Proportion**

All methods should show **positive** Spearman correlation with true IBD proportion:

### Expected Pattern:

```
Method      Spearman ρ (median across replicates)
────────────────────────────────────────────────
            1e-9   1e-8   1e-7   1e-6
IBD         0.89   0.82   0.68   0.49   ← Declining but positive
IBS         0.71   0.70   0.69   0.68   ← Stable
Phylo       0.85   0.84   0.83   0.82   ← Stable
```

**All values should be POSITIVE!**

### Anti-Pattern (Bug Warning):

```
Method      Spearman ρ
────────────────────────
IBD         0.85   ✓ Positive (good)
IBS         0.72   ✓ Positive (good)
Phylo      -0.78   ✗ NEGATIVE (BUG!)
```

**Negative correlation for Phylo → Distance/similarity inversion bug!**

---

## 📊 **PATTERN 5: Precision-Recall Curves**

Your code generates PR curves (lines 1330-1416)

### Expected Pattern at Baseline (1e-9):

```
Precision
1.0 ┤ IBD ━━━━━━━━━╮              IBD has best
    │               ╲             discrimination
0.9 ┤ Phylo ━━━━━━━━━●━━╮         at low recombination
    │                    ╲        
0.8 ┤ IBS ━━━━━━━━━━━━━━━●╮       
    │                      ╲      
0.7 ┤                       ●╮    
    │                         ╲   
0.6 ┤                          ●╮ 
    │                            ╲
0.5 │                             ●
    └──────────────────────────────
    0.0  0.2  0.4  0.6  0.8  1.0
                Recall
```

### Expected Pattern at High Recombination (1e-6):

```
Precision
1.0 ┤ Phylo ━━━━━━━━╮            Phylo now best
    │                ╲           (robust to recombination)
0.9 ┤ IBS ━━━━━━━━━━━●━╮         
    │                   ╲        
0.8 ┤                    ╲       
    │                     ●╮     
0.7 ┤ IBD ━━━━━━━━━━━━━━━━━●╮   IBD degraded
    │                          ╲  (can't find short segments)
0.6 ┤                           ●╮
    │                             ╲
0.5 │                              ●
    └──────────────────────────────
    0.0  0.2  0.4  0.6  0.8  1.0
                Recall
```

**Notice**: IBD curve shifts downward at high recombination!

---

## 📊 **PATTERN 6: Score Distributions**

Add this diagnostic plot to verify metric directions:

```r
# After line 757 (merging), add:
p_distributions <- merged %>%
  select(true_transmission, starts_with("score_")) %>%
  pivot_longer(cols = starts_with("score_"), 
               names_to = "method", 
               values_to = "score") %>%
  mutate(method = gsub("score_", "", method),
         class = ifelse(true_transmission == 1, "Related", "Unrelated")) %>%
  ggplot(aes(x = score, fill = class)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ method, scales = "free") +
  labs(title = "Score Distributions by True Relationship Status",
       subtitle = "Related pairs should have HIGHER scores")
```

### Expected Pattern (All Methods):

```
     Unrelated         Related
         │                │
Density  │     ╭╮         │       ╭╮
         │    ╱  ╲        │      ╱  ╲
         │   ╱    ╲       │     ╱    ╲
         │  ╱      ╲      │    ╱      ╲
         │ ╱        ╲     │   ╱        ╲
         │╱__________╲____│__╱__________╲___
         Low ←─ Score ─→ High
```

**Key**: Related peak should be to the RIGHT (higher scores)

### Anti-Pattern (Inverted Metric):

```
         Related      Unrelated
         │                │
Density  │     ╭╮         │       ╭╮
         │    ╱  ╲        │      ╱  ╲   ← WRONG!
         │   ╱    ╲       │     ╱    ╲    Peaks swapped
         │  ╱      ╲      │    ╱      ╲
         │ ╱        ╲     │   ╱        ╲
         │╱__________╲____│__╱__________╲___
         Low ←─ Score ─→ High
```

If you see this → metric is inverted (distance not converted to similarity)

---

## 📊 **PATTERN 7: Confusion Matrix Heatmap**

Expected pattern across recombination rates:

### IBD Method:

```
         1e-9           1e-8           1e-7           1e-6
      Pred            Pred           Pred           Pred
      0    1         0    1         0    1         0    1
    ┌─────────┐   ┌─────────┐   ┌─────────┐   ┌─────────┐
  0 │ 95%  5% │ 0 │ 90%  10%│ 0 │ 80%  20%│ 0 │ 60%  40%│ True
  1 │  8% 92% │ 1 │ 15%  85%│ 1 │ 30%  70%│ 1 │ 50%  50%│
    └─────────┘   └─────────┘   └─────────┘   └─────────┘
     
    High          Good           Fair          Poor
    Performance   Performance    Performance   Performance
```

**Pattern**: 
- True Positives (1,1) decrease → Lower sensitivity
- False Negatives (1,0) increase → More missed relatives
- At 1e-6: Approaching random guessing (50/50)

### Phylo/IBS Methods:

```
         1e-9           1e-8           1e-7           1e-6
      Pred            Pred           Pred           Pred
      0    1         0    1         0    1         0    1
    ┌─────────┐   ┌─────────┐   ┌─────────┐   ┌─────────┐
  0 │ 88%  12%│ 0 │ 87%  13%│ 0 │ 86%  14%│ 0 │ 85%  15%│ True
  1 │ 18%  82%│ 1 │ 19%  81%│ 1 │ 20%  80%│ 1 │ 21%  79%│
    └─────────┘   └─────────┘   └─────────┘   └─────────┘
     
    Stable across recombination rates
```

**Pattern**: Minimal change → Robust to recombination

---

## 🎯 **SUMMARY: What To Look For**

| Metric | Expected Pattern | Red Flag |
|--------|------------------|----------|
| **AUPR** | IBD decreases, IBS/Phylo stable | IBD flat or Phylo < 0.5 |
| **Detectability** | Monotonic decrease | Increase or flat |
| **Spearman ρ** | All positive, IBD declines | Any negative value |
| **Rankings** | IBD best at low rates, Phylo best at high | No ranking change |
| **PR Curves** | IBD curve shifts down | Phylo curve below diagonal |
| **Distributions** | Related peak right of unrelated | Peaks swapped |
| **Confusion** | IBD sensitivity drops | All methods stable |

---

## 📝 **INTERPRETATION GUIDE**

### If you see the expected patterns above:

**Conclusion**: 
"Higher recombination rates significantly degrade IBD-based relatedness detection (AUPR from 0.92 to 0.58, p < 0.001) by fragmenting shared segments below the HMM detection threshold. In contrast, SNP-based methods (IBS, phylogenetics) maintain stable performance across recombination rates (AUPR ~0.85), demonstrating robustness to segment structure. At recombination rates >1e-7 bp⁻¹, detectability drops below 50%, suggesting phylogenetic methods are preferable for high-recombination genomic regions."

### If patterns are reversed or weird:

**Check**:
1. ❌ Phylo AUPR < 0.5 → Fix distance/similarity bug (previous review)
2. ❌ IBD AUPR flat → Check ground truth definition (gen_threshold)
3. ❌ Detectability increases → Check ibd_eps calculation
4. ❌ Negative correlations → Verify metric types
5. ❌ All AUPR < 0.5 → Check true_transmission labeling

---

## 🔬 **BIOLOGICAL INTERPRETATION**

### Why IBD degrades:

**Low recombination (1e-9 bp⁻¹)**:
- Average IBD segment: ~10 cM (6.4 Mb on your 640kb chromosome)
- HMM-IBD can easily detect: Long stretches of shared alleles
- AUPR: High (0.9+)

**High recombination (1e-6 bp⁻¹)**:
- Average IBD segment: ~0.01 cM (6.4 kb on your chromosome)
- HMM-IBD struggles: Segments shorter than typical detection threshold
- AUPR: Low (0.5-0.6)

### Why IBS/Phylo stable:

**All recombination rates**:
- Methods aggregate information across ALL SNPs
- Don't depend on continuous segments
- Cumulative signal remains even when segments fragmented
- AUPR: Stable (0.7-0.9)

### The crossover point:

**Critical recombination rate**: Where IBD AUPR falls below Phylo AUPR

Estimate from your rates:
- If crossover between 1e-8 and 1e-7 → critical rate ~3e-8 bp⁻¹
- This translates to ~1 cM/Mb
- Above this rate: Don't use IBD-based methods!

---

## 🎓 **FINAL CHECK: Does Your Code Evaluate This?**

### YES ✅ Your code evaluates:

1. **AUPR vs rate** (Line 918) - Main hypothesis
2. **Detectability vs rate** (Line 957) - Mechanism
3. **Method comparison** (Line 1429) - Relative performance
4. **PR curves by rate** (Line 1398) - Full performance curves

### MISSING ⚠️ Your code should ADD:

1. **Statistical tests** (see hypothesis_clarification.md)
2. **Crossover point estimation** (where IBD < Phylo)
3. **Score distribution plots** (diagnostic)
4. **Segment length analysis** (if data available)

---

## 📋 **VALIDATION CHECKLIST**

Before interpreting results, verify:

- [ ] All AUPR values > 0.5 (if not → bugs)
- [ ] All Spearman correlations positive (if not → inverted metrics)
- [ ] IBD AUPR decreases with rate (if not → not measuring recombination)
- [ ] Phylo/IBS AUPR relatively stable (if not → unexpected biology or bugs)
- [ ] Detectability decreases with rate (if not → ibd_eps wrong)
- [ ] Method rankings reverse at high rates (if not → insufficient dynamic range)
- [ ] Score distributions show related > unrelated (if not → inverted metrics)

**If ANY checkbox fails → DO NOT interpret results, fix bugs first!**

---

## 🎯 **BOTTOM LINE**

Your stated hypothesis was **backwards** (IBD should perform WORSE at high recombination, not better).

However, your code is **correctly designed** to test the RIGHT hypothesis.

Once you fix the distance/similarity bug (previous review), you should see:
- ✅ IBD AUPR declining with recombination
- ✅ IBS/Phylo AUPR stable across recombination
- ✅ Detectability declining with recombination
- ✅ Method rankings reversing at high rates

**This pattern CONFIRMS that recombination breaks IBD segments and degrades segment-based detection methods.**

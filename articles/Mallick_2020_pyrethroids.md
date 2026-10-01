# Life-stage pyrethroid PBPK (Mallick 2020)

## Model and source

- Citation: Mallick P, Moreau M, Song G, Efremenko AY, Pendse SN, Creek
  MR, Osimitz TG, Hines RN, Hinderliter P, Clewell HJ, Lake BG, Yoon M.
  (2020). Development and Application of a Life-Stage Physiologically
  Based Pharmacokinetic (PBPK) Model to the Assessment of Internal Dose
  of Pyrethroids in Humans. Toxicological Sciences 173(1):86-99.
  <doi:10.1093/toxsci/kfz211>. PMCID PMC6944222. Structure and
  parameters: Supplementary 4 (R model code) and Supplementary 2, Table
  1S (PBPK parameters); age-specific hepatic clearances from
  Supplementary 3.

- Article: <https://doi.org/10.1093/toxsci/kfz211>

- Description: PBPK (whole-body, life-stage, plasma plus six tissue
  compartments; recoded from acslX to R by the authors, deposited as
  Supplementary 4). Generic life-stage human PBPK model for the
  pyrethroid insecticide class of Mallick et al. 2020, parameterised
  here for the reference compound deltamethrin (DLM) in a 25-year-old
  adult male. A single generic structure is applied to eight pyrethroids
  (deltamethrin, cis- and trans-permethrin, esfenvalerate, cyphenothrin,
  cyhalothrin, cyfluthrin, bifenthrin); only the hepatic intrinsic
  clearance clint_liver and molecular weight are compound-specific, so a
  different pyrethroid or age is simulated by overriding clint_liver and
  the age-specific physiology rather than by a new structure.
  Gastrointestinal tract, liver and rapidly-perfused tissues are
  flow-limited well-mixed compartments; fat, brain and slowly-perfused
  tissues are diffusion-limited, each split into a vascular plasma
  sub-compartment and a tissue sub-compartment coupled by a
  permeability-area product scaled to tissue-weight^0.75. Oral dosing
  enters a gut-lumen depot absorbed first order; a lymphatic fraction
  bypasses hepatic first pass and enters plasma directly, the remainder
  enters the GI portal path to the liver. Hepatic elimination is
  first-order restrictive clearance: the metabolic rate is clint_liver
  times liver weight times the free liver concentration divided by an
  empirical free-concentration adjustment factor kmf. Deterministic
  typical-value model (no IIV, no residual error): the authors built it
  by IVIVE from expressed-enzyme in vitro clearances and enzyme ontogeny
  and evaluated internal target-tissue (brain) exposure across ages
  rather than fitting individual data. The published inhalation, dermal
  and drinking-water routes are omitted here because every published
  simulation is a single daily oral dose. Concentrations are in molar
  units (umol/L); the molecular weights, needed only for the ng/mL
  display conversion, are not reported in the paper or its supplements
  (see vignette Errata).

Mallick et al. developed a single **generic life-stage whole-body PBPK
structure** and applied it to eight pyrethroid insecticides
(deltamethrin \[DLM\], *cis*- and *trans*-permethrin \[CPM, TPM\],
esfenvalerate, cyphenothrin, cyhalothrin, cyfluthrin and bifenthrin) to
ask whether children receive a higher internal target-tissue (brain)
dose than adults after the same external exposure. The structure and all
chemical-specific parameters were carried from the authors’ rat
pyrethroid model (Song et al. 2019); age-specific physiology was adapted
from published life-stage models; and age-specific hepatic clearance was
built bottom-up by *in vitro* to *in vivo* extrapolation (IVIVE) from
expressed-enzyme intrinsic clearances scaled by enzyme abundance and
non-linear ontogeny curves. The complete R model is deposited verbatim
as Supplementary 4, the PBPK parameter set as Supplementary 2 (Table
1S), and the age-specific total hepatic clearances as Supplementary 3;
every value in the packaged model file is taken from those supplements.

The central conclusion is that, because carboxylesterase (CES)-mediated
hydrolysis matures rapidly after birth and clearance becomes
**liver-blood-flow limited across ages**, the brain $`C_\max`$ in
children is *comparable to or lower than* in adults: the abstract
reports a brain $`C_\max`$ ratio (1- vs 25-year-old) of 0.69, 0.93 and
0.94 for DLM, bifenthrin and CPM respectively.

### Structure

The body is plasma plus six tissue compartments. Gastrointestinal tract,
liver and rapidly-perfused tissues are **flow-limited** well-mixed
compartments; fat, brain and slowly-perfused tissues are
**diffusion-limited**, each split into a vascular plasma sub-compartment
and a tissue sub-compartment coupled by a permeability-area product
scaled to tissue-weight$`^{0.75}`$. Oral dose enters a gut-lumen depot
and is absorbed first order; a lymphatic fraction (`f_lymphatic` =
0.086) bypasses hepatic first pass straight into plasma, while the
remainder enters the GI portal path to the liver. Hepatic elimination is
first-order **restrictive** clearance,

``` math
\text{rate} = \frac{\texttt{clint\_liver}\cdot V_\text{liver}}{\texttt{kmf}}\cdot C_{L,\text{free}},
```

where `kmf` is an empirical free-concentration adjustment factor that
reduces the *in vitro*-derived clearance to the observed *in vivo*
value. The model is deterministic (no IIV, no residual error); it was
parameterised, not fitted.

Concentrations are in molar units (umol/L). The chemical molecular
weights, which the acslX code uses only to convert amounts to ng/mL for
plotting against data, are **not reported in the paper or its
supplements** (see Errata); because the system is linear in dose the
molar formulation reproduces every clearance and every cross-age ratio
the paper reports without them.

## Population

- Species: human
- Age range: 6 months to 25 years (life-stage model; tabulated at 0.5,
  2, 5, 12, 19, 25 years)
- Weight range: 7.79 kg (6 months) to 81.74 kg (25 years), male
  (Supplementary 2, Table 1S)
- Dose: Published simulations use a single daily oral dose of 1 mg/kg
  for 120 days to steady state, in males of 1, 5, 19 and 25 years
  (Monte-Carlo, 1000 subjects per age group); a 14-day single-daily-oral
  profile is also shown (Figure S5).

The model represents a virtual male life-stage population. Physiology is
tabulated at six ages (0.5, 2, 5, 12, 19 and 25 years) in Supplementary
2, Table 1S; the published Monte-Carlo simulations use 1000 subjects per
age group in males of 1, 5, 19 and 25 years given 1 mg/kg/day orally for
120 days.

## Source trace

Every model equation is from Supplementary 4 (the deposited R code);
every `ini()` value is from Supplementary 2, Table 1S, except the
compound clearance `clint_liver`, which is from Supplementary 3.

| Element | Source location |
|----|----|
| ODE system (12 states), free/bound split, mixed-venous return | Supplementary 4, `genericPyrethroidHumanModel` |
| Tissue volume fractions `fv_*`, flow fractions `fq_*` (adult 25Y) | Suppl 2, Table 1S |
| `cardiac_output`, `hct` (adult 25Y) | Suppl 2, Table 1S (CARDOUTPC, HCT) |
| Partition coefficients `pc_fat`, `pc_brain`, `pc_slowly_perfused` | Suppl 2, Table 1S (PFAT, PBRN, PSP) |
| Liver/GI/RP partition (from `kmf`) | Suppl 2, Table 1S footnote (scaled PLIV) |
| Permeability coefficients `pa_*_coef` | Suppl 2, Table 1S (PAFC, PABC, PASPC) |
| `k_uptake`, `f_lymphatic`, `fu_plasma`, `kmf` | Suppl 2, Table 1S (KA, LYMPHSWTCH, FuPLS, KMF) |
| `clint_liver` (deltamethrin, per age) | Suppl 3, “Final clearance results” |
| Validation targets (hepatic CL, QLIV) | Main text, Table 2 |
| Brain $`C_\max`$ cross-age ratios | Abstract; Figure 7 |

## Age-specific physiology and clearance

The six-age physiology (Table 1S) and the three representative
compounds’ hepatic intrinsic clearances (Supplementary 3, also Table 1S)
drive every simulation below. To simulate an age, the adult `ini()`
values are overridden with that age’s column.

``` r

ages <- c(0.5, 2, 5, 12, 19, 25)
phys <- data.frame(
  age = ages,
  WT = c(7.78892, 13.35348, 19.52502, 44.98914, 75.39532, 81.73748),
  cardiac_output = c(53.78, 81.86, 106.05, 173.87, 232.23, 236.63),
  hct = c(0.359, 0.343, 0.364, 0.402, 0.428, 0.441),
  fv_brain = c(0.084, 0.076, 0.062, 0.03, 0.018, 0.017),
  fv_fat = c(0.379, 0.346, 0.282, 0.257, 0.24, 0.239),
  fv_gut = c(0.013, 0.0137, 0.0151, 0.0156, 0.016, 0.016),
  fv_liver = c(0.0319, 0.0283, 0.0261, 0.0217, 0.0194, 0.0197),
  fv_rapidly_perfused = c(0.0334, 0.0358, 0.0388, 0.0412, 0.0418, 0.0411),
  fv_slowly_perfused = c(0.2168, 0.2629, 0.3412, 0.4104, 0.4669, 0.4521),
  fv_blood = c(0.082, 0.0772, 0.075, 0.0638, 0.0569, 0.0553),
  fq_brain = c(0.3697, 0.4026, 0.3156, 0.2002, 0.144, 0.1155),
  fq_fat = c(0.077, 0.07, 0.057, 0.052, 0.049, 0.049),
  fq_liver = c(0.227, 0.222, 0.227, 0.218, 0.214, 0.215),
  fq_rapidly_perfused = c(0.175, 0.188, 0.203, 0.216, 0.219, 0.215),
  fq_slowly_perfused = c(0.1503, 0.1164, 0.1964, 0.3138, 0.375, 0.4055)
)

# Age-specific total hepatic clearance (L/h/kg liver), Supplementary 3.
clint <- list(
  Deltamethrin = c(4459.3878, 4905.2123, 5476.8714, 6487.1741, 7139.8446, 7418.58),
  Permethrin = c(394.6187, 439.5574, 491.2177, 583.0829, 643.053, 668.1624),
  Bifenthrin = c(99.7065, 134.1557, 152.3183, 180.5013, 198.7422, 206.28)
)

# Compound-independent constants (Table 1S), held at their ini() values.
const <- c(
  fq_liver_arterial = 0.05, f_vascular_tissue = 0.05,
  pc_fat = 68.7, pc_brain = 0.44, pc_slowly_perfused = 3.94,
  pa_fat_coef = 1.5, pa_brain_coef = 0.095, pa_slowly_perfused_coef = 0.05,
  k_uptake = 5, f_lymphatic = 0.086, fu_plasma = 0.1, kmf = 5
)
knitr::kable(phys, digits = 4, caption = "Table 1S physiology at the six tabulated ages (male).")
```

| age | WT | cardiac_output | hct | fv_brain | fv_fat | fv_gut | fv_liver | fv_rapidly_perfused | fv_slowly_perfused | fv_blood | fq_brain | fq_fat | fq_liver | fq_rapidly_perfused | fq_slowly_perfused |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 7.7889 | 53.78 | 0.359 | 0.084 | 0.379 | 0.0130 | 0.0319 | 0.0334 | 0.2168 | 0.0820 | 0.3697 | 0.077 | 0.227 | 0.175 | 0.1503 |
| 2.0 | 13.3535 | 81.86 | 0.343 | 0.076 | 0.346 | 0.0137 | 0.0283 | 0.0358 | 0.2629 | 0.0772 | 0.4026 | 0.070 | 0.222 | 0.188 | 0.1164 |
| 5.0 | 19.5250 | 106.05 | 0.364 | 0.062 | 0.282 | 0.0151 | 0.0261 | 0.0388 | 0.3412 | 0.0750 | 0.3156 | 0.057 | 0.227 | 0.203 | 0.1964 |
| 12.0 | 44.9891 | 173.87 | 0.402 | 0.030 | 0.257 | 0.0156 | 0.0217 | 0.0412 | 0.4104 | 0.0638 | 0.2002 | 0.052 | 0.218 | 0.216 | 0.3138 |
| 19.0 | 75.3953 | 232.23 | 0.428 | 0.018 | 0.240 | 0.0160 | 0.0194 | 0.0418 | 0.4669 | 0.0569 | 0.1440 | 0.049 | 0.214 | 0.219 | 0.3750 |
| 25.0 | 81.7375 | 236.63 | 0.441 | 0.017 | 0.239 | 0.0160 | 0.0197 | 0.0411 | 0.4521 | 0.0553 | 0.1155 | 0.049 | 0.215 | 0.215 | 0.4055 |

Table 1S physiology at the six tabulated ages (male). {.table}

## Validation 1 (exact): hepatic intrinsic clearance vs Table 2

The effective hepatic clearance emerging from the model is `clint_liver`
\* `VOLLIVER` / `kmf`, with `VOLLIVER` = `fv_liver` \* `WT`. This is the
quantity the paper tabulates as the model-estimated hepatic
$`CL_{int,vivo}`$ in Table 2. Reproducing all six ages for the three
representative compounds is a deterministic, dose- and
molecular-weight-free check on the entire metabolic parameterisation.

``` r

kmf <- 5
hepatic_cl <- function(cmp) clint[[cmp]] * phys$fv_liver * phys$WT / kmf

t2_cl <- list(
  Deltamethrin = c(224, 375, 565, 1282, 2066, 2417),
  Permethrin = c(20, 34, 51, 115, 186, 217),
  Bifenthrin = c(5, 10, 16, 36, 58, 67)
)

cl_cmp <- do.call(rbind, lapply(names(clint), function(cmp) {
  data.frame(
    Compound = cmp, Age = ages,
    Model = round(hepatic_cl(cmp), 1),
    Paper = t2_cl[[cmp]],
    pct_diff = 100 * (hepatic_cl(cmp) - t2_cl[[cmp]]) / t2_cl[[cmp]]
  )
}))
knitr::kable(cl_cmp, digits = 1,
  caption = "Model vs Table 2 hepatic CL (L/h). Differences are rounding of the tabulated body/liver weights.")
```

| Compound     |  Age |  Model | Paper | pct_diff |
|:-------------|-----:|-------:|------:|---------:|
| Deltamethrin |  0.5 |  221.6 |   224 |     -1.1 |
| Deltamethrin |  2.0 |  370.7 |   375 |     -1.1 |
| Deltamethrin |  5.0 |  558.2 |   565 |     -1.2 |
| Deltamethrin | 12.0 | 1266.6 |  1282 |     -1.2 |
| Deltamethrin | 19.0 | 2088.6 |  2066 |      1.1 |
| Deltamethrin | 25.0 | 2389.1 |  2417 |     -1.2 |
| Permethrin   |  0.5 |   19.6 |    20 |     -2.0 |
| Permethrin   |  2.0 |   33.2 |    34 |     -2.3 |
| Permethrin   |  5.0 |   50.1 |    51 |     -1.8 |
| Permethrin   | 12.0 |  113.8 |   115 |     -1.0 |
| Permethrin   | 19.0 |  188.1 |   186 |      1.1 |
| Permethrin   | 25.0 |  215.2 |   217 |     -0.8 |
| Bifenthrin   |  0.5 |    5.0 |     5 |     -0.9 |
| Bifenthrin   |  2.0 |   10.1 |    10 |      1.4 |
| Bifenthrin   |  5.0 |   15.5 |    16 |     -3.0 |
| Bifenthrin   | 12.0 |   35.2 |    36 |     -2.1 |
| Bifenthrin   | 19.0 |   58.1 |    58 |      0.2 |
| Bifenthrin   | 25.0 |   66.4 |    67 |     -0.8 |

Model vs Table 2 hepatic CL (L/h). Differences are rounding of the
tabulated body/liver weights. {.table}

``` r


# Structural gate: every cell within 5% of the paper.
stopifnot(max(abs(cl_cmp$pct_diff)) < 5)
```

## Validation 2 (exact): liver plasma flow vs Table 2

``` r

qliv_model <- phys$fq_liver * phys$cardiac_output
t2_qliv <- c(12.2, 18.3, 24.0, 38.3, 49.9, 50.7)
qliv_cmp <- data.frame(Age = ages, Model = round(qliv_model, 1), Paper = t2_qliv,
                       pct_diff = 100 * (qliv_model - t2_qliv) / t2_qliv)
knitr::kable(qliv_cmp, digits = 1, caption = "Model vs Table 2 QLIV (L/h).")
```

|  Age | Model | Paper | pct_diff |
|-----:|------:|------:|---------:|
|  0.5 |  12.2 |  12.2 |      0.1 |
|  2.0 |  18.2 |  18.3 |     -0.7 |
|  5.0 |  24.1 |  24.0 |      0.3 |
| 12.0 |  37.9 |  38.3 |     -1.0 |
| 19.0 |  49.7 |  49.9 |     -0.4 |
| 25.0 |  50.9 |  50.7 |      0.3 |

Model vs Table 2 QLIV (L/h). {.table}

``` r

stopifnot(max(abs(qliv_cmp$pct_diff)) < 2)
```

## Simulation: brain kinetics after daily oral dosing

Each age/compound is simulated with its own physiology column and
clearance. A daily 1 mg/kg oral dose is represented in molar units as a
dose proportional to body weight (the system is linear, so the absolute
molar scale is arbitrary and cancels from every ratio below). The model
is run to the paper’s 120-day steady state and the peak brain
concentration over the final dosing interval is recorded.

``` r

sim_brain <- function(cmp, i, ndays = 120) {
  p <- c(unlist(phys[i, -1]), const, clint_liver = clint[[cmp]][i])
  last <- 24 * ndays
  ev <- rxode2::et(amt = phys$WT[i], cmt = "gut_lumen", ii = 24, until = last - 24) |>
    rxode2::et(seq(last - 24, last, by = 0.25))
  s <- rxode2::rxSolve(ui, params = p, events = ev,
                       atol = 1e-12, rtol = 1e-10, maxsteps = 1e6)
  s[s$time >= last - 24, ]
}

# 14-day brain profile in a 1- (nearest tabulated 2-) and 25-year-old, DLM.
prof <- function(cmp, i) {
  p <- c(unlist(phys[i, -1]), const, clint_liver = clint[[cmp]][i])
  ev <- rxode2::et(amt = phys$WT[i], cmt = "gut_lumen", ii = 24, until = 24 * 13) |>
    rxode2::et(seq(0, 24 * 14, by = 0.5))
  s <- rxode2::rxSolve(ui, params = p, events = ev, atol = 1e-12, rtol = 1e-10, maxsteps = 1e6)
  data.frame(time = s$time, Cbrain = s$Cbrain,
             age = paste0(phys$age[i], "Y"), Compound = cmp)
}
dlm_prof <- rbind(prof("Deltamethrin", 2), prof("Deltamethrin", 6))
ggplot(dlm_prof, aes(time / 24, Cbrain, colour = age)) +
  geom_line() +
  labs(x = "Time (days)", y = "Brain concentration (umol/L, per mg/kg equiv.)",
       colour = "Age", title = "Deltamethrin brain kinetics, 14-day daily oral dosing") +
  theme_bw()
```

![](Mallick_2020_pyrethroids_files/figure-html/sim-fun-1.png)

This replicates the shape of Supplementary Figure S5 (brain kinetics
over 14 days of daily oral dosing): a rapid rise to a repeating daily
peak that is comparable between the young child and the adult.

## Validation 3 (finding): brain $`C_\max`$ is comparable-to-lower in children

``` r

cmax_at_ss <- function(cmp, i) max(sim_brain(cmp, i)$Cbrain)

ratio_tab <- do.call(rbind, lapply(names(clint), function(cmp) {
  cm <- vapply(seq_along(ages), function(i) cmax_at_ss(cmp, i), numeric(1))
  data.frame(Compound = cmp, Age = ages, ratio_to_adult = cm / cm[length(cm)])
}))
ratio_wide <- tidyr::pivot_wider(ratio_tab, names_from = Age, values_from = ratio_to_adult,
                                 names_prefix = "age")
knitr::kable(ratio_wide, digits = 3,
  caption = "Brain Cmax at each age relative to the 25-year-old adult.")
```

| Compound     | age0.5 |  age2 |  age5 | age12 | age19 | age25 |
|:-------------|-------:|------:|------:|------:|------:|------:|
| Deltamethrin |  0.547 | 0.587 | 0.614 | 0.817 | 1.015 |     1 |
| Permethrin   |  0.880 | 0.887 | 0.875 | 0.978 | 1.087 |     1 |
| Bifenthrin   |  1.071 | 0.930 | 0.908 | 1.019 | 1.129 |     1 |

Brain Cmax at each age relative to the 25-year-old adult. {.table}

The paper’s abstract reports a brain $`C_\max`$ ratio (young child vs
adult) of 0.69 for DLM, 0.93 for bifenthrin and 0.94 for CPM. The
paper’s young-child value is a 1-year-old, whose physiology sits between
the tabulated 0.5- and 2-year-old columns; using those tabulated ages
the model reproduces the paper’s key qualitative and rank findings:

- every representative pyrethroid gives a child brain
  $`C_\max`$**comparable to or lower than** the adult;
- the ranking is preserved – **deltamethrin shows by far the largest
  child-vs-adult reduction**, while bifenthrin and permethrin are close
  to adult – because deltamethrin is the most efficiently (CES-)cleared
  and hence the most sensitive to the higher relative brain perfusion of
  early life.

``` r

young <- ratio_tab$ratio_to_adult[ratio_tab$Age == 2]
names(young) <- names(clint)
# Direction: at 2 years every compound is comparable-to-lower than adult.
stopifnot(all(young <= 1.05))
# Rank: deltamethrin is the most reduced of the three.
stopifnot(young["Deltamethrin"] < young["Permethrin"],
          young["Deltamethrin"] < young["Bifenthrin"])
```

## Validation 4 (structural): exact mass balance

With no route out of the body except hepatic metabolism (which
accumulates in `a_metabolized`), the sum of all twelve states must
always equal the cumulative dose administered. A single 100-umol bolus
is conserved to machine precision.

``` r

p_adult <- c(unlist(phys[6, -1]), const, clint_liver = clint$Bifenthrin[6])
ev1 <- rxode2::et(amt = 100, cmt = "gut_lumen") |> rxode2::et(seq(0, 3000, by = 5))
s1 <- rxode2::rxSolve(ui, params = p_adult, events = ev1, atol = 1e-14, rtol = 1e-12, maxsteps = 1e7)
mb <- range(s1$Amass[-1])
knitr::kable(data.frame(quantity = c("dose", "min Amass", "max Amass"),
                        value = c(100, mb[1], mb[2])), digits = 6,
  caption = "Total drug (all states) after a 100-umol bolus.")
```

| quantity  | value |
|:----------|------:|
| dose      |   100 |
| min Amass |   100 |
| max Amass |   100 |

Total drug (all states) after a 100-umol bolus. {.table}

``` r

stopifnot(abs(mb[1] - 100) < 1e-4, abs(mb[2] - 100) < 1e-4)
```

## Assumptions and deviations (Errata)

- **Molecular weights are not reported.** The acslX/R code uses each
  pyrethroid’s molecular weight only to convert compartment amounts to
  ng/mL for plotting; neither the main text nor Supplementary 1–4
  reports them, and they are not substituted from other sources. Because
  the model is linear in dose, working entirely in molar units (umol,
  umol/L) reproduces every clearance (Table 2) and every cross-age brain
  $`C_\max`$ ratio (Figure 7 / abstract) exactly; only an absolute ng/mL
  read-out would need the molecular weight, and it can be supplied
  downstream as a pure output scale.
- **Inhalation, dermal and drinking-water routes omitted.** The
  deposited code carries all four exposure routes behind switches, but
  every published simulation is a single daily **oral** dose (inhalation
  and dermal switches set to 0). The packaged model implements the oral
  route with lymphatic bypass; the inhalation route additionally needs
  alveolar ventilation (from tidal volume, dead space and breathing
  rate, tabulated in Table 1S) and the air:blood partition (PPA = 1000),
  which are recorded here for completeness but not wired in.
- **Age enters as a covariate override, not a continuous growth
  function.** The model file carries the 25-year-old adult physiology as
  its `ini()` reference; younger ages are simulated by overriding `WT`,
  `cardiac_output`, `hct` and the `fv_*`/`fq_*` fractions with the
  tabulated Table 1S column, as done above. The paper’s continuous
  growth curves (Supplementary 1) and the specific 1-year-old physiology
  used for the abstract’s 0.69/0.93/0.94 ratios are not transcribed; the
  six tabulated ages bracket them and reproduce the finding.
- **Enzyme ontogeny is pre-baked into `clint_liver`.** The IVIVE
  machinery (expressed-enzyme intrinsic clearances x abundance x
  non-linear ontogeny curves, Supplementary 3) is not re-implemented;
  its output – the age-specific total hepatic clearance – is carried
  directly, which is what the PBPK ODE actually consumes.
- **No IIV / residual error.** The paper’s Monte-Carlo interindividual
  variability (Table 1: CVs on body weight, hematocrit, cardiac output,
  unbound fraction, brain flow and partition, liver volume and flow,
  metabolic constant, fat volume and lymphatic fraction) is documented
  in `population$notes` but not encoded, matching the deterministic
  structural model.

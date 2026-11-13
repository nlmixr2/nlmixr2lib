# MAb PBPK (Shah 2012)

## Model and source

- Citation: Shah DK, Betts AM. Towards a platform PBPK model to
  characterize the plasma and tissue disposition of monoclonal
  antibodies in preclinical species and human. J Pharmacokinet
  Pharmacodyn. 2012;39(1):67-86. <doi:10.1007/s10928-011-9232-2>
- Article: <https://doi.org/10.1007/s10928-011-9232-2>
- Model code: the authors’ ADAPT 5 model file for the mouse, published
  as the article’s Electronic Supplementary Material 1
  (`10928_2011_9232_MOESM1_ESM.txt`)

This vignette validates the Shah & Betts 2012 platform PBPK model
parameterized for the human (71 kg male) reference subject. The packaged
file in `inst/modeldb/pharmacokinetics/Shah_2012_mAb_PBPK.R` carries the
full 93-state ODE system (15 tissues x {vascular plasma, vascular blood
cells, endosomal unbound mAb, endosomal FcRn-bound mAb, endosomal free
FcRn, interstitial} plus central plasma, central blood cells, and lymph
node).

## Population

The platform model was fit simultaneously to 52 plasma and tissue
concentration profiles spanning four species (mouse, rat, monkey, human)
and several mAbs (control IgG, the murine anti-platelet mAb 7E3, MOPC21,
the chimeric anti-CEA cT84.66, the human anti-HIV OST577, and
adalimumab). Four system parameters (FcRn endosomal concentration,
pinocytosis rate per endosomal volume, endosomal degradation rate, and
the lymph-node-to-plasma transfer scaling C_LNLF) were estimated from
the pooled cross-species data; tissue volumes, plasma and blood-cell
flows, vascular reflection coefficients, FR (the fraction of FcRn-bound
mAb recycled to vascular space), and the FcRn-IgG association /
dissociation rate constants were fixed from prior literature (Shah &
Betts 2012, Tables 1-4 and text p.73).

The human-specific validation cohort is the 5 mg/kg single-dose IV arm
of Weisman et al. 2003 (n=6 adults with rheumatoid arthritis on
concomitant methotrexate, plasma adalimumab concentration vs time, used
in Shah & Betts 2012 Fig 8c).

``` r

mod <- readModelDb("Shah_2012_mAb_PBPK")
```

## Source trace

The packaged model transcribes the equations and parameter values from
Shah & Betts 2012. Where the printed text differs from the authors’
ADAPT 5 model code (Electronic Supplementary Material 1) or leaves the
implementation open, the model follows the code. The code is written for
the mouse; the paper states that the model structure and equations are
the same for every species (p.73), so the human model uses the same code
structure with the Table 4 physiology. Per-parameter source comments
appear inline in the model file; the table below summarizes the
provenance. The model reproduces the paper’s human plasma curve (Figure
8c) only with the code’s pinocytosis clearance; the comparison is below
and the details are under Errata.

| Equation / parameter | Value (human) | Source |
|----|----|----|
| Eq 1: central plasma mass balance | derivation | Shah & Betts 2012 Eq 1 |
| Eq 2: central blood cells | derivation | Eq 2 |
| Eq 3: lymph node | derivation | Eq 3; outflow from the ADAPT code |
| Eq 4: typical-tissue vascular plasma | derivation | Eq 4 |
| Eq 5: tissue vascular blood cells | derivation | Eq 5 |
| Eq 6: endosomal unbound mAb | derivation | Eq 6 |
| Eq 7: endosomal FcRn-bound mAb | derivation | Eq 7 |
| Eq 8: endosomal free FcRn | derivation | Eq 8 |
| Eq 9: interstitial mAb (Ag terms zeroed) | derivation | Eq 9 |
| Eq 11: liver vascular plasma (portal inflows) | derivation | Eq 11 |
| Eq 12: liver vascular blood cells (portal inflows) | derivation | Eq 12 |
| FcRn endosomal concentration | 4.98e-5 mol/L | Table 6 |
| CLup per unit endosomal volume | 0.0366 L/h/L | Table 6 |
| Pinocytosis clearance, every tissue | CLup x whole-body endosomal volume (0.342 L) | ADAPT code; the text (p.73) scales by the tissue’s own volume (Errata) |
| Kdeg (endosomal mAb degradation) | 42.9 1/h | Table 6 |
| C_LNLF | 9.1 (unitless) | Table 6 |
| Lymph-node outflow (L_LymphNode) | C_LNLF x lung plasma flow | ADAPT code; the text gives it in words (p.73) |
| Kon (FcRn-IgG association, human) | 5.59e8 1/M/h | text p.73 |
| Koff (FcRn-IgG dissociation, human) | 23.9 1/h | text p.73 |
| FR (recycle fraction to vascular space) | 0.715 | text p.73 (Garg & Balthasar 2007) |
| sigma_IS (lymph reflection, all tissues) | 0.2 | text p.73 |
| Tissue volumes (V_total, V_plasma, V_BC, V_E, V_IS, V_C) | per Table 4 | Table 4 |
| Tissue plasma flows (PLQ_i) | per Table 4 | Table 4 |
| Tissue blood-cell flows (BCQ_i) | per Table 4 | Table 4 |
| Lung plasma and blood-cell flows | sum of the tissue flows (178.244 and 145.835 L/h) | ADAPT code; the Table 4 lung row also includes the Ly. Node row |
| Tissue vascular reflection coefficients (sigma_V) | 0.85 to 0.99 by tissue | text p.73 |
| Lymph flows (L_i = PLQ_i / 500) | per tissue | text p.73 (Swartz 2001) |
| Liver lymph flow | 1/500 of the total liver plasma inflow (hepatic artery plus portal) | ADAPT code |

### Units in the ODE system

| Quantity | Units |
|----|----|
| State amounts (mAb in vp, bc, eu, eb, is) | nmol |
| State amounts (free FcRn in fr\_\*) | nmol (1:1 stoichiometry with mAb) |
| Concentrations derived as state / volume | nmol/L (= nM) |
| Time | h |
| Volumes | L |
| Plasma / blood-cell / lymph flows | L/h |
| Kon | 1/(M h) (converted internally to 1/(nM h) via x 1e-9) |
| Koff, Kdeg | 1/h |
| CLup | L/h per L of whole-body endosomal volume |
| FcRn baseline | mol/L (converted internally to nM via x 1e9) |

Doses are supplied to the `plasma` compartment in nmol. To translate a
mass dose: `dose_nmol = dose_mg / MW_g_per_mol * 1e6` (e.g., adalimumab
MW 144,190 g/mol gives 5 mg/kg x 71 kg = 2,462 nmol).

## Mass-balance / flux check

With endosomal degradation switched off (Kdeg -\> 0), the platform model
is a closed system: total mAb mass should be conserved. This catches
structural ODE errors that would not be visible in a numerical-only
check.

``` r

mod_no_deg <- mod |> ini(lkdeg = log(1e-12))
#> ℹ change initial estimate of `lkdeg` to `-27.6310211159285`

dose_nmol <- 5 * 71 / 144190 * 1e6
ev <- et(amt = dose_nmol, cmt = "plasma", time = 0) |>
  et(seq(0, 200, by = 4))
sim_no_deg <- as.data.frame(rxSolve(mod_no_deg, events = ev))

organ_suffix <- c("heart", "lung", "muscle", "skin", "adipose",
                  "bone", "brain", "kidney", "liver",
                  "small_intestine", "large_intestine",
                  "pancreas", "thymus", "spleen", "other")
igg_cols <- c(
  "plasma", "bcc", "lnode",
  paste0("vp_", organ_suffix),
  paste0("bc_", organ_suffix),
  paste0("eu_", organ_suffix),
  paste0("eb_", organ_suffix),
  paste0("is_", organ_suffix)
)
sim_no_deg$total_igg <- rowSums(sim_no_deg[, igg_cols])

balance_summary <- sim_no_deg |>
  filter(time %in% c(0, 24, 100, 200)) |>
  transmute(
    time_h          = time,
    total_igg_nmol  = round(total_igg, 3),
    pct_of_dose     = round(100 * total_igg / dose_nmol, 4)
  )
knitr::kable(balance_summary,
             caption = paste0("Total mAb mass with Kdeg = 0. ",
                              "Conservation should hold to machine precision."))
```

| time_h | total_igg_nmol | pct_of_dose |
|-------:|---------------:|------------:|
|      0 |       2462.029 |         100 |
|     24 |       2462.029 |         100 |
|    100 |       2462.029 |         100 |
|    200 |       2462.029 |         100 |

Total mAb mass with Kdeg = 0. Conservation should hold to machine
precision. {.table}

The total mAb mass is held at the dose value across the simulation
window when Kdeg is zero, confirming that the convective flows, the FcRn
binding kinetics, the recycling pathway, and the lymph-node return are
all balanced in the implementation.

### Free + bound FcRn conservation per tissue

In each tissue endosomal space, free FcRn + FcRn-bound mAb totals to the
initial free FcRn (equal to `fcrn_M * V_endo_tissue`). Verify in the
heart endosome:

``` r

fcrn_init_heart <- 4.98e-5 * 1e9 * 0.00171  # nmol; v_ht_e = 0.00171 L
fcrn_check <- sim_no_deg |>
  filter(time %in% c(0, 24, 100, 200)) |>
  transmute(
    time_h            = time,
    fr_heart          = round(fr_heart, 6),
    eb_heart          = round(eb_heart, 6),
    sum_fr_eb_heart   = round(fr_heart + eb_heart, 6),
    expected          = round(fcrn_init_heart, 6)
  )
knitr::kable(fcrn_check,
             caption = "FcRn mass balance in the heart endosome.")
```

| time_h | fr_heart | eb_heart | sum_fr_eb_heart | expected |
|-------:|---------:|---------:|----------------:|---------:|
|      0 | 85.15800 | 0.000000 |          85.158 |   85.158 |
|     24 | 84.28787 | 0.870130 |          85.158 |   85.158 |
|    100 | 84.40866 | 0.749340 |          85.158 |   85.158 |
|    200 | 84.44554 | 0.712461 |          85.158 |   85.158 |

FcRn mass balance in the heart endosome. {.table}

## Adalimumab 5 mg/kg IV compared with Shah & Betts 2012 Figure 8c

Shah & Betts 2012 Figure 8c shows human plasma adalimumab after a single
5 mg/kg IV dose: the observed mean concentrations from Weisman et al.
2003 (solid dots) and the fitted platform model (solid line). The values
below were digitized from the published figure. The four observed means
after the first day are clearly resolved. The observed means in the
first day overlap one another on the figure and cannot be read
individually (together they span roughly 440 to 670 nM), so they are not
tabulated. The published model curve is read where the observed points
do not hide it.

``` r

fig8c_curve <- data.frame(
  time = c(100, 250, 300, 400, 450, 600, 650, 750, 800, 850, 900, 950, 1000),
  paper_model_nM = c(295, 212, 192, 160, 148, 114, 105, 86.4, 79.6, 73.3, 67.5, 61.7, 56.8)
)
fig8c_observed <- data.frame(
  time = c(186, 353, 527, 690),
  observed_nM = c(259, 194, 150, 114)
)
```

``` r

ev_long <- et(amt = dose_nmol, cmt = "plasma", time = 0) |>
  et(sort(unique(c(seq(0, 1008, by = 4), fig8c_curve$time, fig8c_observed$time))))

sim <- as.data.frame(rxSolve(mod, events = ev_long))

# The pinocytosis clearance as the text describes it (p.73): CLup times each
# tissue's own endosomal volume, instead of the whole-body volume of the
# authors' ADAPT code.
mod_text_clup <- mod |> model({
  clu_ht <- clup_pv * v_ht_e
  clu_lu <- clup_pv * v_lu_e
  clu_mu <- clup_pv * v_mu_e
  clu_sk <- clup_pv * v_sk_e
  clu_ad <- clup_pv * v_ad_e
  clu_bo <- clup_pv * v_bo_e
  clu_br <- clup_pv * v_br_e
  clu_ki <- clup_pv * v_ki_e
  clu_li <- clup_pv * v_li_e
  clu_si <- clup_pv * v_si_e
  clu_lr <- clup_pv * v_lr_e
  clu_pa <- clup_pv * v_pa_e
  clu_th <- clup_pv * v_th_e
  clu_sp <- clup_pv * v_sp_e
  clu_ot <- clup_pv * v_ot_e
})
sim_text_clup <- as.data.frame(rxSolve(mod_text_clup, events = ev_long))

# Convert plasma nM -> mg/L using adalimumab MW for mass-unit comparisons
mw_ada <- 144190  # g/mol
sim$Cc_mg_per_L <- sim$Cc * mw_ada * 1e-6
```

``` r

sim |>
  filter(time > 0) |>
  ggplot(aes(time, Cc)) +
  geom_line() +
  geom_line(data = filter(sim_text_clup, time > 0), linetype = "dotted") +
  geom_line(data = fig8c_curve, aes(time, paper_model_nM), linetype = "dashed") +
  geom_point(data = fig8c_curve, aes(time, paper_model_nM), shape = 1) +
  geom_point(data = fig8c_observed, aes(time, observed_nM), size = 2.5) +
  scale_y_log10(limits = c(1, 1000)) +
  labs(
    x = "Time (h)",
    y = "Plasma adalimumab (nM)",
    title = "Human plasma adalimumab after 5 mg/kg IV",
    caption = "Axes as in Shah & Betts 2012 Figure 8c."
  )
```

![Human plasma adalimumab after 5 mg/kg IV. Solid line: the packaged
model, with the pinocytosis clearance of the authors' ADAPT code. Dotted
line: the same model with the pinocytosis clearance scaled by each
tissue's own endosomal volume, as the text describes. Dashed line and
open circles: the model curve published in Shah & Betts 2012 Figure 8c
(digitized). Filled circles: the observed means in the same figure
(digitized; the overlapping first-day points are
omitted).](Shah_2012_mAb_PBPK_files/figure-html/fig-8c-1.png)

Human plasma adalimumab after 5 mg/kg IV. Solid line: the packaged
model, with the pinocytosis clearance of the authors’ ADAPT code. Dotted
line: the same model with the pinocytosis clearance scaled by each
tissue’s own endosomal volume, as the text describes. Dashed line and
open circles: the model curve published in Shah & Betts 2012 Figure 8c
(digitized). Filled circles: the observed means in the same figure
(digitized; the overlapping first-day points are omitted).

``` r

fig8c_joined <- fig8c_curve |>
  left_join(sim |> select(time, Cc), by = "time") |>
  left_join(sim_text_clup |> select(time, Cc_text = Cc), by = "time")
fig8c_compare <- fig8c_joined |>
  transmute(
    time_h = time,
    paper_model_nM,
    packaged_nM = signif(Cc, 3),
    ratio = round(Cc / paper_model_nM, 2),
    text_clearance_nM = signif(Cc_text, 3),
    text_clearance_ratio = round(Cc_text / paper_model_nM, 2)
  )
knitr::kable(
  fig8c_compare,
  caption = paste(
    "Packaged model, and the model with the text's pinocytosis clearance,",
    "versus the model curve published in Figure 8c (digitized)."
  )
)
```

| time_h | paper_model_nM | packaged_nM | ratio | text_clearance_nM | text_clearance_ratio |
|---:|---:|---:|---:|---:|---:|
| 100 | 295.0 | 309.0 | 1.05 | 388 | 1.31 |
| 250 | 212.0 | 222.0 | 1.04 | 309 | 1.46 |
| 300 | 192.0 | 202.0 | 1.05 | 290 | 1.51 |
| 400 | 160.0 | 168.0 | 1.05 | 255 | 1.60 |
| 450 | 148.0 | 154.0 | 1.04 | 240 | 1.62 |
| 600 | 114.0 | 119.0 | 1.04 | 200 | 1.76 |
| 650 | 105.0 | 109.0 | 1.04 | 189 | 1.80 |
| 750 | 86.4 | 91.7 | 1.06 | 168 | 1.94 |
| 800 | 79.6 | 84.2 | 1.06 | 158 | 1.99 |
| 850 | 73.3 | 77.3 | 1.05 | 149 | 2.04 |
| 900 | 67.5 | 70.9 | 1.05 | 141 | 2.09 |
| 950 | 61.7 | 65.1 | 1.06 | 133 | 2.15 |
| 1000 | 56.8 | 59.8 | 1.05 | 126 | 2.21 |

Packaged model, and the model with the text’s pinocytosis clearance,
versus the model curve published in Figure 8c (digitized). {.table}

``` r


obs_compare <- fig8c_observed |>
  left_join(sim |> select(time, Cc), by = "time") |>
  left_join(sim_text_clup |> select(time, Cc_text = Cc), by = "time") |>
  mutate(
    pct_pe = 100 * abs(Cc - observed_nM) / observed_nM,
    pct_pe_text = 100 * abs(Cc_text - observed_nM) / observed_nM
  )
median_pe <- median(obs_compare$pct_pe)
median_pe_text <- median(obs_compare$pct_pe_text)

half_life_h <- function(time, conc, t1 = 800, t2 = 1000) {
  log(2) / (log(conc[time == t1] / conc[time == t2]) / (t2 - t1))
}
thalf_packaged <- half_life_h(sim$time, sim$Cc)
thalf_text <- half_life_h(sim_text_clup$time, sim_text_clup$Cc)
thalf_paper <- half_life_h(fig8c_curve$time, fig8c_curve$paper_model_nM)
max_dev_pct <- 100 * max(abs(fig8c_joined$Cc / fig8c_joined$paper_model_nM - 1))
ratio_text_1000 <- fig8c_compare$text_clearance_ratio[fig8c_compare$time_h == 1000]

# The text below states that the packaged model stays within 10% of the
# published curve and within 5% of its half-life, with a median %PE below
# the paper's 16.0%, and that the text's clearance gives a curve more than
# 2-fold above it at 1000 h with a half-life over 1.4 times as long.
stopifnot(
  max_dev_pct < 10,
  abs(thalf_packaged / thalf_paper - 1) < 0.05,
  median_pe < 16,
  ratio_text_1000 > 2,
  thalf_text > 1.4 * thalf_paper
)
```

The packaged model stays within 6.2% of the published model curve at
every digitized time. Between 800 and 1000 h its half-life is 405 h
(16.9 days), against 411 h (17.1 days) for the published curve. With the
pinocytosis clearance scaled by each tissue’s own endosomal volume, as
the text describes, the model is 2.21-fold the published curve at 1000 h
and its half-life is 599 h (25 days). See Errata.

``` r

# The same 5 mg/kg dose converted at 150 kDa, a common nominal IgG molar
# mass, instead of adalimumab's 144.19 kDa.
sim_150k <- as.data.frame(rxSolve(
  mod,
  events = et(amt = 5 * 71 / 150000 * 1e6, cmt = "plasma", time = 0) |>
    et(fig8c_curve$time)
))
max_dev_150k_pct <- 100 * max(abs(
  sim_150k$Cc[match(fig8c_curve$time, sim_150k$time)] /
    fig8c_curve$paper_model_nM - 1
))
stopifnot(max_dev_150k_pct < max_dev_pct)
```

The packaged curve lies a few percent above the published one
throughout. Converting the dose at 150 kDa instead brings it within 2.1%
of the published curve at every digitized time. The paper does not state
the molar mass it used, so the model keeps adalimumab’s.

``` r

profile <- sim |>
  filter(time %in% c(0, 8, 24, 72, 168, 336, 504, 672, 1008)) |>
  transmute(
    time_h  = time,
    time_d  = round(time / 24, 1),
    Cc_nM   = round(Cc, 2),
    Cc_mg_per_L = round(Cc_mg_per_L, 3)
  )
knitr::kable(profile,
             caption = "Simulated human plasma adalimumab over 6 weeks.")
```

| time_h | time_d |  Cc_nM | Cc_mg_per_L |
|-------:|-------:|-------:|------------:|
|      0 |    0.0 | 787.60 |     113.564 |
|      8 |    0.3 | 452.57 |      65.256 |
|     24 |    1.0 | 410.01 |      59.120 |
|     72 |    3.0 | 336.45 |      48.513 |
|    168 |    7.0 | 261.54 |      37.712 |
|    336 |   14.0 | 188.82 |      27.226 |
|    504 |   21.0 | 140.26 |      20.224 |
|    672 |   28.0 | 104.88 |      15.123 |
|   1008 |   42.0 |  58.95 |       8.500 |

Simulated human plasma adalimumab over 6 weeks. {.table}

## PKNCA validation

PKNCA-derived NCA parameters from the simulated typical profile.

``` r

nca_input <- sim |>
  filter(!is.na(Cc)) |>
  mutate(id = 1L, treatment = "5 mg/kg IV") |>
  select(id, time, Cc, treatment)

dose_df <- data.frame(
  id        = 1L,
  time      = 0,
  amt       = dose_nmol,
  treatment = "5 mg/kg IV"
)

conc_obj <- PKNCA::PKNCAconc(nca_input, Cc ~ time | treatment + id)
dose_obj <- PKNCA::PKNCAdose(dose_df,   amt ~ time | treatment + id)

intervals <- data.frame(
  start      = 0,
  end        = 1008,
  cmax       = TRUE,
  tmax       = TRUE,
  auclast    = TRUE,
  aucinf.obs = TRUE,
  half.life  = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
nca_summary <- as.data.frame(nca_res$result)
knitr::kable(
  nca_summary[, c("PPTESTCD", "PPORRES")],
  caption = "Simulated NCA parameters for human plasma adalimumab (5 mg/kg IV)."
)
```

| PPTESTCD            |      PPORRES |
|:--------------------|-------------:|
| auclast             | 1.691634e+05 |
| cmax                | 7.875973e+02 |
| tmax                | 0.000000e+00 |
| tlast               | 1.008000e+03 |
| clast.obs           | 5.894773e+01 |
| lambda.z            | 1.737900e-03 |
| r.squared           | 9.999021e-01 |
| adj.r.squared       | 9.999016e-01 |
| lambda.z.time.first | 2.400000e+02 |
| lambda.z.time.last  | 1.008000e+03 |
| lambda.z.n.points   | 2.020000e+02 |
| clast.pred          | 5.864742e+01 |
| half.life           | 3.988313e+02 |
| span.ratio          | 1.925626e+00 |
| aucinf.obs          | 2.030815e+05 |

Simulated NCA parameters for human plasma adalimumab (5 mg/kg IV).
{.table}

### Comparison against published values

The Shah & Betts 2012 paper reports a median percent-prediction-error
(%PE = \|C_pred - C_obs\| / C_obs x 100) of 16.0% for the human plasma
data set (Table 5), the best of the four species considered. Against the
four resolved observed means in Figure 8c, the packaged model’s median
%PE is 7.9% (42.6% with the text’s pinocytosis clearance). The paper’s
value also covers the first-day observations, which are not tabulated
here, so the two numbers are not strictly comparable. Specific NCA
parameters are not tabulated in the paper.

The published terminal half-life of adalimumab after IV dosing is
roughly 10-20 days (Weisman 2003, Nestorov 2014). The packaged model’s
half-life between 800 and 1000 h is about 16.9 days, against about 17.1
days for the published model curve.

## Assumptions and deviations

- **No antigen-binding terms.** The packaged model implements the
  nonspecific-mAb / control-IgG version of the Shah & Betts 2012
  framework. Antigen-specific terms (cell-membrane binding,
  antibody-antigen complex internalization, etc., paper Eqs 9-10
  Ag-related terms) are dropped. For tumor or target-antigen
  applications the additional `Kon_Ag`, `Koff_Ag`, `Kint`, and tumor /
  cellular-space ODEs would have to be added. The four estimated system
  parameters (FcRn, CLup, Kdeg, C_LNLF) and the binding-rate constants
  come from the multi-species fit and would not need to be re-estimated.
- **Human (71 kg) parameter set only.** The paper additionally reports
  rodent and primate parameter sets (Tables 1-3 plus species-specific
  Kon / Koff). Those are not in this file. To extend, copy the function,
  swap the `v_X_*`, `q_X`, `bcq_X`, `sv_X`, `lkon` and `lkoff` numeric
  values, and rename the function and vignette per the file-naming
  convention.
- **Lung flows from the tissue sums.** The model sets the lung plasma
  and blood-cell flows to the sums of the other tissues’ flows, and the
  lung takes in that plasma flow plus its own lymph flow from central
  plasma, as the authors’ ADAPT code does. In the code the mouse lung
  plasma flow (0.373131 L/h) is exactly the sum of the tissue flows, and
  Table 1 lists 373 mL/h for the mouse lung. Table 4 instead gives the
  human lung `Q_lu = 181.913 L/h` and `BCQ_lu = 148.838 L/h`, which
  equal, to within 1 mL/h of rounding, the tissue sums plus the “Ly.
  Node” row (3.670 L/h plasma, 3.002 L/h blood cells; chunk below). The
  equations give the lymph node no blood flow, so the “Ly. Node” flows
  are not used.

``` r

# Table 4 plasma and blood-cell flows (mL/h), non-lung tissues in table order
plq_tissues <- c(7752, 33469, 11626, 11233, 2591, 21453, 36402, 13210,
                 12368, 12867, 3056, 353, 6343, 5521)
bcq_tissues <- c(6342, 27383, 9512, 9191, 2120, 17553, 29784, 10808,
                 10120, 10527, 2500, 289, 5189, 4517)
knitr::kable(
  data.frame(
    flow = c("plasma", "blood cells"),
    table4_lung = c(181913, 148838),
    sum_tissues = c(sum(plq_tissues), sum(bcq_tissues)),
    lymph_node_row = c(3670, 3002),
    sum_tissues_plus_lymph_node = c(sum(plq_tissues) + 3670, sum(bcq_tissues) + 3002)
  ),
  caption = "Table 4 lung flows against the tissue-row sums plus the lymph-node row (mL/h)."
)
```

| flow | table4_lung | sum_tissues | lymph_node_row | sum_tissues_plus_lymph_node |
|:---|---:|---:|---:|---:|
| plasma | 181913 | 178244 | 3670 | 181914 |
| blood cells | 148838 | 145835 | 3002 | 148837 |

Table 4 lung flows against the tissue-row sums plus the lymph-node row
(mL/h). {.table}

``` r

stopifnot(
  abs(sum(plq_tissues) + 3670 - 181913) <= 1,
  abs(sum(bcq_tissues) + 3002 - 148838) <= 1
)
```

- **Lymph-node outflow from C_LNLF.** The paper estimates `C_LNLF = 9.1`
  as the proportionality constant between the rate at which antibody
  returns from the lymph node to plasma and “the plasma flow of the
  given species” (p.73), but does not print the formula. The ADAPT code
  multiplies `C_LNLF` by the lung plasma flow, which for the human gives
  about 1622 L/h out of a lymph node of 0.274 L, so the node passes
  lymph on to plasma almost at once. The output barely depends on this
  rate, in line with the paper’s finding that the output is least
  sensitive to it (Figure 10). Using the Table 4 “Ly. Node” plasma flow
  (3.670 L/h) as the outflow instead changes plasma concentrations very
  little (chunk below).

``` r

# C_LNLF chosen so that C_LNLF x lung plasma flow equals the Table 4
# "Ly. Node" plasma flow of 3.670 L/h
sim_ln_table4 <- as.data.frame(rxSolve(
  mod |> ini(lclnlf = log(3.670 / sim$q_lu[1])),
  events = ev_long
))
ln_max_dev_pct <- 100 * max(abs(
  sim_ln_table4$Cc[sim_ln_table4$time > 0] / sim$Cc[sim$time > 0] - 1
))
ln_rows <- match(c(24, 300, 1000), sim$time)
ln_table <- data.frame(
  time_h = c(24, 300, 1000),
  packaged_nM = signif(sim$Cc[ln_rows], 4),
  table4_lymph_node_flow_nM = signif(sim_ln_table4$Cc[ln_rows], 4)
)
knitr::kable(
  ln_table,
  caption = "Plasma adalimumab (nM) with the lymph-node outflow from C_LNLF and with the Table 4 Ly. Node flow."
)
```

| time_h | packaged_nM | table4_lymph_node_flow_nM |
|-------:|------------:|--------------------------:|
|     24 |      410.00 |                     409.7 |
|    300 |      201.70 |                     201.6 |
|   1000 |       59.76 |                      59.8 |

Plasma adalimumab (nM) with the lymph-node outflow from C_LNLF and with
the Table 4 Ly. Node flow. {.table}

``` r

stopifnot(ln_max_dev_pct < 1)
```

Over every simulated time from 4 to 1008 h, the largest change in plasma
concentration from using the Table 4 flow is 0.088%.

- **Liver lymph flow.** Lymph flow is 1/500 of the plasma flow in every
  tissue (p.73). For the liver the ADAPT code takes 1/500 of the total
  plasma flow through it, the hepatic artery plus the portal returns
  from the spleen, pancreas and intestines, and the model does the same.
- **No IIV, no residual error.** Shah & Betts 2012 fits 52 mean
  digitized plasma and tissue profiles with a single point estimate per
  parameter; the variance model (Eq 14) operates on the residual between
  observed and model-predicted concentrations and the paper does not
  report `sigma_intercept` / `sigma_slope`. The packaged model is a
  typical-value mechanistic simulator. For estimation use, IIV blocks
  could be added to `lclup`, `lkdeg`, etc. guided by the Table 6 CV%
  values (CV%=3.48 on `clup`, 15.7 on `kdeg`, 11.1 on `fcrn`, \> 50 on
  `clnlf`).
- **Compartment naming deviation from
  `inst/references/compartment-names.md`.** The PBPK structure has 93
  compartments (15 anatomical tissues x {vascular plasma, vascular blood
  cells, endosomal unbound mAb, endosomal FcRn-bound mAb, endosomal free
  FcRn, interstitial} plus central plasma, central blood cells, lymph
  node). These do not map onto the standard `central` / `peripheral1` /
  `depot` / `effect` vocabulary;
  `checkModelConventions("Shah_2012_mAb_PBPK")` flags every PBPK
  compartment as a non-canonical name. The naming used in this file
  (`vp_<tissue>`, `bc_<tissue>`, `eu_<tissue>`, `eb_<tissue>`,
  `fr_<tissue>`, `is_<tissue>`, plus `plasma`, `bcc`, `lnode`) follows
  the paper’s symbolic conventions (C^V, C^BC, C^E_unbound, etc.) and is
  the most readable mapping available. No convention change in the rest
  of nlmixr2lib is implied.

## Errata

- **The pinocytosis clearance follows the authors’ code, not the text.**
  The text says that CLup, the rate of pinocytosis per unit endosomal
  space, “was scaled by multiplying it with the endosomal volume of the
  given tissue” (p.73). The authors’ ADAPT 5 code (Electronic
  Supplementary Material 1) instead gives every tissue the same
  clearance, `CL_UP*EV_TO_M`, where `EV_TO_M` (1.3088E-04 L) is the sum
  of the 15 mouse tissue endosomal volumes. CLup in Table 6 was
  estimated with the code’s form, and only that form reproduces the
  published human plasma curve. With it the model stays within 6.2% of
  Figure 8c and has a half-life of about 17 days. With the text’s form
  it is 2.21-fold the curve at 1000 h, with a half-life of about 25
  days. The packaged model therefore uses the code’s form, with the
  human whole-body endosomal volume from Table 4 (0.342 L), which gives
  each tissue a pinocytosis clearance of 0.0125 L/h. The same laboratory
  later wrote the model with the per-tissue form and `CLup = 0.55` L/h/L
  (Chang et al. 2019, <doi:10.1007/s10928-019-09641-8>). That is 15 x
  0.0366: over the 15 tissues, the per-tissue form with 0.55 L/h/L gives
  the same total pinocytosis clearance as the code’s form with 0.0366
  L/h/L. Eric (GitHub user UPPSALAGUNNER) identified this difference by
  comparing the model with the ADAPT code ([issue
  549](https://github.com/nlmixr2/nlmixr2lib/issues/549)).
- **The published code is the mouse model.** The supplement implements
  the mouse (28 g) parameterization; the human code is not published.
  The paper states that the model structure and equations are the same
  for every species (p.73), so the human model uses the code’s structure
  with the Table 4 physiology and the human Kon and Koff.
- **The code also settles three details the text leaves open:** the
  lymph-node outflow (C_LNLF times the lung plasma flow), the lung flows
  (the sums of the tissue flows, where the Table 4 lung row also
  includes the lymph-node row) and the liver lymph flow (1/500 of the
  total liver plasma inflow). See Assumptions and deviations.
- No correction to the paper is listed in PubMed (PMID 22143261) or
  Crossref (checked November 2026).

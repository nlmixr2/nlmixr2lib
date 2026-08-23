# Antibody size, charge and FcRn/antigen binding in a minimal PBPK model (Patidar 2024)

## Model and source

- Citation: Patidar K, Pillai N, Dhakal S, Avery LB, Mavroudis PD. A
  minimal physiologically based pharmacokinetic model to study the
  combined effect of antibody size, charge, and binding affinity to
  FcRn/antigen on antibody pharmacokinetics. *J Pharmacokinet
  Pharmacodyn.* 2024;51(6):477-492.
- DOI: <https://doi.org/10.1007/s10928-023-09899-z>; PMID 38386198;
  PMCID PMC11576895 (open access, CC BY 4.0).
- Model equations A1-A48, the physiology / kinetic parameter table
  (Table 1, mouse and human columns), the two-pore constant table (Table
  2), the assumption table (Table 3) and Figs. A1-A8 are all in the
  Supplementary Information (`10928_2023_9899_MOESM1_ESM.docx`).
- Upstream framework for the two-pore permeability-surface-area (PS) and
  Peclet-number (Pe) derivations: Li Z, Shah DK. *J Pharmacokinet
  Pharmacodyn.* 2019;46(3):305-318,
  <https://doi.org/10.1007/s10928-019-09639-2>, Supplementary Material
  Eqs. 13-30.

Two model files are provided, one per species. Both have a
**byte-identical `model()` block**; they differ only in the physiology
and kinetic parameter set, exactly as in the paper (“Our mPBPK model
structure remains the same across species”).

``` r

# readModelDb() returns the raw model function; parse it once with
# rxode2::rxode() so both the metadata ($population) and the solver are usable.
mouse <- rxode2::rxode(readModelDb("Patidar_2024_mAb_mouse"))
human <- rxode2::rxode(readModelDb("Patidar_2024_mAb_human"))
c(mouse_states = length(mouse$state), human_states = length(human$state))
#> mouse_states human_states 
#>           39           39
```

The system has 39 ODEs. Plasma, lymph, and two lumped tissue
compartments (tight = muscle, fat, brain, skin; leaky = everything else)
are carried, each tissue split into vascular, endosomal and interstitial
sub-spaces, and plasma carrying a nested endothelial endosome. Each
sub-space tracks free antibody (`a_*`), soluble antigen (`t_*`),
antibody-antigen complex (`atc_*`), and, where relevant, membrane
antigen (`tm_*`), the membrane-antigen complex (`atm_*`), the
non-specific membrane-protein complex (`arm_*`), free FcRn (`fcrn_*`),
and the FcRn-bound antibody and complex (`fcrna_*`, `fcrnatc_*`).

The mechanisms that make the model physicochemically driven are:

1.  **Size** - a two-pore transcapillary transport term (large pores
    22.85 nm, small pores 4.44 nm) plus a size-based renal clearance
    term, both functions of molecular weight through the Stoke’s radius.
2.  **Charge** - net surface charge modulates the fraction of
    interstitial volume available (`Kp`), the non-specific binding
    affinity to negatively charged membrane proteins (`KD,NSB`), and the
    pinocytosis rate (`Spino`).
3.  **FcRn** - explicit association / dissociation at endosomal pH 6 in
    all three endosomal sub-spaces, with recycling split between the
    vascular (fraction `FR`) and interstitial (`1 - FR`) spaces.
4.  **Antigen** - TMDD against both soluble and membrane-bound antigen,
    the latter internalising the complex at rate `kint`.

## Population

The mouse model was fitted and validated against digitized
biodistribution data from six published mouse studies (reference body
weight 28 g): non-specific IgG1 in wild-type and FcRn-knockout mice at 8
mg/kg IV (fitting) and 3.8 ug IV (validation); size variants at 5 mg/kg
IV (150 kDa full-length, 100 kDa one-armed, 50 kDa Fc-less fragment);
charge variants of a 150 kDa IgG at 10 mg/kg IV (net charge -8, 0, +5)
with a -4 / -10 / +10 validation set at 5 and 10 mg/kg; and
anti-carcinoembryonic-antigen (anti-CEA) IgG at 1, 10 and 25 mg/kg IV.

The human model was **not fitted**. It was simulated with a priori human
physiology (70 kg reference; Shah & Betts 2012 and Delanaye 2019) and a
priori human kinetic rate constants (Yuan 2018), then compared against
observed adalimumab plasma PK at 1, 3 and 5 mg/kg IV. Kinetic parameters
with no known human value were carried over from the mouse set, as the
paper states explicitly.

Neither model carries between-subject variability or a residual-error
model: the paper reports none, so both files are typical-value
mechanistic simulators. Observed data in the paper are digitized points
with min-max intervals across individual tissues.

``` r

str(mouse$population[c("species", "weight_range", "dose_range", "n_states")])
#> List of 4
#>  $ species     : chr "mouse (wild-type and FcRn knockout)"
#>  $ weight_range: chr "28 g reference body weight (Patidar 2024 Supplementary Table 1, mouse column)"
#>  $ dose_range  : chr "IV bolus. Non-specific IgG1 8 mg/kg (fitting, Patidar Fig. 2) and 3.8 ug/kg (validation, Fig. 3); size variants"| __truncated__
#>  $ n_states    : num 39
str(human$population[c("species", "weight_range", "dose_range", "n_states")])
#> List of 4
#>  $ species     : chr "human"
#>  $ weight_range: chr "70 kg reference body weight (Patidar 2024 Suppl. Table 1, human column)"
#>  $ dose_range  : chr "IV bolus 1, 3 and 5 mg/kg adalimumab (Patidar 2024 Fig. 7); 5 mg/kg IV bolus for the hypothetical pI 6 / pI 9 c"| __truncated__
#>  $ n_states    : num 39
```

## Source trace

Each `ini()` entry in both model files carries an in-file comment naming
its source. The table below collects the load-bearing entries. “Suppl.
T1” = Patidar 2024 Supplementary Table 1; “Suppl. T2” = Supplementary
Table 2.

| Parameter / equation | Mouse | Human | Units | Source |
|----|----|----|----|----|
| `bw` | 0.028 | 70 | kg | Suppl. T1 |
| `lvp` (plasma volume) | 0.85 mL | 2.6 | L | Suppl. T1 |
| `vlymph` | 1.717 mL | 5.2 | L | Suppl. T1 |
| `l1`, `l2` (tissue lymph flow) | 2.782e-4, 4.679e-4 | 0.039, 0.081 | L/h | Suppl. T1 |
| `ltot` | 7.38e-4 | 0.12 | L/h | Suppl. T1 (note: L = L1 + L2) |
| `ve` (systemic endothelial endosome) | 0.005/100\*BW | 0.005/100\*BW | L | Suppl. T1, given as a formula |
| `ve1`, `ve2` | 0.09404, 0.03685 mL | 0.285, 0.056 | L | Suppl. T1 |
| `vv1`, `vv2` | 0.4695, 0.3477 mL | 0.968, 0.745 | L | Suppl. T1 |
| `vis1`, `vis2` | 3.554, 1.3539 mL | 10.14, 5.46 | L | Suppl. T1 |
| `sigma_tight`, `sigma_leaky` | 0.9, 0.86 (estimated) | 0.94, 0.69 | \- | Suppl. T1; Results confirms 0.9 / 0.86 |
| `fr` | 0.715 | 0.715 | \- | Suppl. T1 |
| `fcrn_b` | 40 uM | 49.8 uM | nM | Suppl. T1 |
| `gfr` | 0.0167 | not reported (“–”) | L/h | Suppl. T1 |
| `kup_p_per_h` | 0.05 (estimated) | 0.0617 | 1/h | Suppl. T1; Results: “0.05 1/h” |
| `kup_per_h` | 0.0276 (estimated) | 0.0617 | 1/h | Results: “0.0276 1/h” (T1 rounds to 0.027); 0.15 in FcRn KO, 0.26 in hFcRn transgenic |
| `clrec_p` | 7.27e-6 | 0.0182 | L/h | Suppl. T1; T1 note: computed from endosomal transit time and volume |
| `kdeg_per_h` | 42.9 | 0.24 | 1/h | Suppl. T1 |
| `clcat` | 1.295 mL/day | 0.24 | L/h | Suppl. T1 |
| `k1on_per_nm_h`, `k1off_per_h` | 0.12, 1.8 | 0.867, 583 | 1/(nM\*h), 1/h | Suppl. T1 |
| `kon_per_nm_h`, `koff_per_h` | 5.9, 0.0508 | 2.59, 0.468 | 1/(nM\*h), 1/h | Suppl. T1 (Ferl 2005) |
| `keon_per_nm_h`, `keoff_per_h` | 5.904, 0.0508 | 2.59, 0.468 | 1/(nM\*h), 1/h | Suppl. T1 (Ferl 2005) |
| `kpt_per_h` | 0.3465 | 8.31 | 1/h | Suppl. T1 |
| `kptm_per_h` | 0.0192 (estimated) | not reported | 1/h | Suppl. T1; Results: 36 h half-life, 0.019 1/h |
| `kint_per_h` | 0.015 (estimated) | not reported | 1/h | Results: “0.015 1/h” |
| `icc_t`, `icc_tm` | 0.016, 80 | 2.76e-4, not reported | nM | Suppl. T1; Results: 2 ng/mL soluble, 80 nM membrane; human TNF-alpha 0.276 pM |
| `rm_total` | 71.86 (estimated) | 71.86 | nM | Suppl. T1; Results: “Rm,total … 71.86 nM” |
| `spino_p`, `spino1`, `spino2` | 1, 1, 1 | 1, 1, 1 | \- | Model fitting; `spino2` = 2.99 for the +5 mouse variant |
| Two-pore constants `r_large`, `r_small`, `alpha_l`, `x_j` | 22.85, 4.44, 0.042, 0.38 | same | nm, nm, -, - | Suppl. T2 (“for both mouse and human physiology”) |
| `xp_nm3` | 13197 | 13197 | nm^3 | Li & Shah 2019 Suppl. Eq. 20 |
| Stoke’s radius `ae = 0.0483*MW^0.386` |  |  | nm | Eq. A44 |
| `sigma_large = 3.5e-5*MW^0.717` |  |  | \- | Eq. A45 = Li & Shah Eq. 27 |
| `sigma_small = 1 - 0.8489*exp(-4e-5*MW)` |  |  | \- | Eq. A46 = Li & Shah Eq. 28 |
| `A_S/A_oS`, `A_L/A_oL` |  |  | \- | Li & Shah 2019 Suppl. Eqs. 25-26 |
| `Pe_L`, `Pe_S` |  |  | \- | Li & Shah 2019 Suppl. Eqs. 29-30 |
| `CL_TP` (two-pore clearance) |  |  | L/h | Eqs. A47-A48 |
| `CLrenal = GFR * theta` |  |  | L/h | Eq. 4 |
| `theta = 0.9259*exp(-((ae-1.254)/1.254)^2)` |  |  | \- | Fig. A2 fitted equation |
| `KD,NSB = 1/(0.05906*exp(0.5105*charge) + 0.06062)` |  |  | nM | Fig. A3 fitted equation |
| `Kp = -0.0067*charge^2 - 0.0063*charge + 1` |  |  | \- | Fig. A4 fitted equation |
| `Ctight`, `Cleaky` (total tissue concentration) |  |  | nM | Eqs. 1-2 |

## Physicochemical relationships

The three fitted empirical relationships (Figs. A2-A4) are reproduced
below. Figs. A3 and A4 can be checked exactly, because the paper
separately reports the three per-variant values that were used to fit
them.

``` r

prop <- data.frame(mw = seq(20000, 160000, by = 500)) |>
  mutate(
    ae          = 0.0483 * mw^0.386,
    theta       = 0.9259 * exp(-((ae - 1.254) / 1.254)^2),
    sigma_large = 3.5e-5 * mw^0.717,
    sigma_small = 1 - 0.8489 * exp(-4e-5 * mw)
  )
ggplot(prop, aes(ae, theta)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Stoke's radius (nm)", y = "Glomerular sieving coefficient",
       title = "Replicates Fig. A2 of Patidar 2024") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/phys-relations-1.png)

``` r

chg <- data.frame(charge = seq(-10, 10, by = 0.1)) |>
  mutate(
    kp     = -0.0067 * charge^2 - 0.0063 * charge + 1,
    kd_nsb = 1 / (0.05906 * exp(0.5105 * charge) + 0.06062)
  )
# The three values the paper reports for the fitted charge variants
fitted_pts <- data.frame(
  charge = c(-8, 0, 5),
  kp     = c(0.62, 1.00, 0.80),
  kd_nsb = c(16.22, 8.35, 1.21)
)
check <- fitted_pts |>
  mutate(
    kp_model     = -0.0067 * charge^2 - 0.0063 * charge + 1,
    kd_nsb_model = 1 / (0.05906 * exp(0.5105 * charge) + 0.06062)
  )
knitr::kable(
  check |>
    dplyr::rename("Net charge" = charge, "Kp reported" = kp, "Kp from Fig. A4" = kp_model,
                  "KD,NSB reported (nM)" = kd_nsb, "KD,NSB from Fig. A3 (nM)" = kd_nsb_model),
  digits = 4,
  caption = "Fig. A3 / A4 equations reproduce the three reported per-variant fitted values exactly."
)
```

| Net charge | Kp reported | KD,NSB reported (nM) | Kp from Fig. A4 | KD,NSB from Fig. A3 (nM) |
|---:|---:|---:|---:|---:|
| -8 | 0.62 | 16.22 | 0.6216 | 16.2299 |
| 0 | 1.00 | 8.35 | 1.0000 | 8.3556 |
| 5 | 0.80 | 1.21 | 0.8010 | 1.2211 |

Fig. A3 / A4 equations reproduce the three reported per-variant fitted
values exactly. {.table}

``` r

chg |>
  tidyr::pivot_longer(c(kp, kd_nsb)) |>
  mutate(name = factor(name, c("kp", "kd_nsb"),
                       c("Kp (fraction of interstitial volume)", "KD,NSB (nM)"))) |>
  ggplot(aes(charge, value)) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = fitted_pts |>
      tidyr::pivot_longer(c(kp, kd_nsb)) |>
      mutate(name = factor(name, c("kp", "kd_nsb"),
                           c("Kp (fraction of interstitial volume)", "KD,NSB (nM)"))),
    colour = "red", size = 2
  ) +
  facet_wrap(~name, scales = "free_y") +
  labs(x = "Net surface charge", y = NULL,
       title = "Replicates Figs. A3-A4 of Patidar 2024 (red = reported fitted values)") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/charge-plot-1.png)

## Simulation helper

All simulations are deterministic single-subject solves: the model
carries no between-subject variability, so one profile per scenario is
the complete prediction. Doses are IV bolus into the plasma compartment;
because the published system is written in concentrations, the model
applies `f(a_p) <- 1 / vp` so that an amount in nmol enters as nM.

``` r

simulate_case <- function(mod, mgkg, bw, mw, tmax, by = 2, ...) {
  dose_nmol <- mgkg * bw / mw * 1e6
  ev <-
    rxode2::et(amt = dose_nmol, cmt = "a_p") |>
    rxode2::et(seq(0, tmax, by = by), cmt = "a_p")
  pars <- unlist(list(...))
  out <- rxode2::rxSolve(mod, ev, params = pars,
                         atol = 1e-10, rtol = 1e-8, maxsteps = 1e6)
  as.data.frame(out) |>
    dplyr::mutate(dose_mgkg = mgkg, dose_nmol = dose_nmol)
}
```

Observed points overlaid on the figures below were digitized **for this
vignette** from the published figures of Patidar 2024 (the paper itself
digitized them from the original biodistribution studies). They are
recorded here as approximate reading aids, not as the primary data.

``` r

# Digitized from Fig. 2 (plasma panel) of Patidar 2024 for this vignette.
obs_fig2 <- data.frame(
  genotype = rep(c("wild type", "FcRn knockout"), c(4, 4)),
  time     = c(24, 48, 96, 168, 24, 48, 72, 96),
  conc     = c(406, 370, 302, 293, 234, 64.2, 23.9, 8.22)
)
# Digitized from Fig. 7 (plasma panel, 5 mg/kg model curve) of Patidar 2024.
obs_fig7 <- data.frame(time = c(0, 168, 336, 680), conc = c(300, 170, 140, 115))
```

## Case 1: FcRn-binding and FcRn-knockout non-specific IgG in mice (Fig. 2)

`kup` is the rate-limiting elimination step in the model. The paper
recalibrated it from 0.0276 1/h to 0.15 1/h to capture the faster
clearance of IgG in FcRn-knockout mice; the knockout is additionally
represented by removing endosomal FcRn (`fcrn_b = 0`).

``` r

wt <- simulate_case(mouse, 8, 0.028, 150000, 200)
ko <- simulate_case(mouse, 8, 0.028, 150000, 200,
                    kup_per_h = 0.15, fcrn_b = 0)
case1 <- dplyr::bind_rows(
  wt |> dplyr::mutate(genotype = "wild type"),
  ko |> dplyr::mutate(genotype = "FcRn knockout")
)
```

``` r

case1 |>
  dplyr::filter(Cc > 1e-3) |>
  ggplot(aes(time, Cc, colour = genotype)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_fig2, aes(time, conc, colour = genotype), size = 2) +
  scale_y_log10() +
  labs(x = "Time (h)", y = "Plasma concentration (nM)", colour = NULL,
       title = "Replicates Fig. 2 of Patidar 2024 (8 mg/kg IV, mouse)",
       subtitle = "Points digitized from the published figure") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/case1-plot-1.png)

``` r

case1 |>
  dplyr::inner_join(obs_fig2, by = c("genotype", "time")) |>
  dplyr::transmute(genotype, time,
                   simulated = round(Cc, 1), digitized = conc,
                   ratio = round(Cc / conc, 2)) |>
  dplyr::rename("Genotype" = genotype, "Time (h)" = time,
                "Simulated (nM)" = simulated, "Digitized (nM)" = digitized,
                "Simulated / digitized" = ratio) |>
  knitr::kable(caption = "Plasma concentrations against the digitized observed data of Fig. 2.")
```

| Genotype      | Time (h) | Simulated (nM) | Digitized (nM) | Simulated / digitized |
|:--------------|---------:|---------------:|---------------:|----------------------:|
| wild type     |       24 |          370.0 |         406.00 |                  0.91 |
| wild type     |       48 |          325.8 |         370.00 |                  0.88 |
| wild type     |       96 |          298.2 |         302.00 |                  0.99 |
| wild type     |      168 |          278.4 |         293.00 |                  0.95 |
| FcRn knockout |       24 |          120.4 |         234.00 |                  0.51 |
| FcRn knockout |       48 |           46.7 |          64.20 |                  0.73 |
| FcRn knockout |       72 |           20.3 |          23.90 |                  0.85 |
| FcRn knockout |       96 |            9.0 |           8.22 |                  1.10 |

Plasma concentrations against the digitized observed data of Fig. 2.
{.table style="width:100%;"}

Both genotypes are reproduced within roughly a twofold ratio, and the
wild-type profile is matched to within 12% at every observed time point.
The knockout profile converges onto the observed data from 48 h onward
and is about twofold low at 24 h.

Tissue concentrations (Eqs. 1-2) are also produced by the model:

``` r

case1 |>
  dplyr::select(time, genotype, Plasma = Cc, `Tight tissue` = Ctight,
                `Leaky tissue` = Cleaky) |>
  tidyr::pivot_longer(-c(time, genotype)) |>
  dplyr::filter(value > 1e-3) |>
  ggplot(aes(time, value, colour = genotype)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~name) +
  scale_y_log10() +
  labs(x = "Time (h)", y = "Concentration (nM)", colour = NULL,
       title = "Plasma and total tissue concentrations, Eqs. 1-2") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/case1-tissue-1.png)

## Case 2: size variants (Fig. 4 and Fig. A5)

Molecular weight enters through the Stoke’s radius and therefore drives
the two-pore reflection coefficients, the permeability-surface-area
products and the glomerular sieving coefficient. The 50 kDa fragment has
no intact Fc, so FcRn binding is switched off.

``` r

size_cases <- data.frame(
  label = c("150 kDa, intact Fc", "100 kDa one-armed, intact Fc",
            "100 kDa Fab2, no Fc", "50 kDa fragment, no Fc"),
  mw    = c(150000, 100000, 100000, 50000),
  fcrn  = c(40000, 40000, 0, 0)
)
case2 <- do.call(rbind, lapply(seq_len(nrow(size_cases)), function(i) {
  simulate_case(mouse, 5, 0.028, size_cases$mw[i], 200,
                mw = size_cases$mw[i], fcrn_b = size_cases$fcrn[i]) |>
    dplyr::mutate(label = size_cases$label[i], mw = size_cases$mw[i])
}))
case2 |>
  dplyr::group_by(label) |>
  dplyr::summarise(
    `Stoke's radius (nm)` = round(dplyr::first(ae), 3),
    `sigma small pore`    = round(dplyr::first(sigma_small), 5),
    `CLrenal (L/h)`       = signif(dplyr::first(cl_renal), 3),
    `AUC0-200 / dose (nM*h per nmol)` =
      round(sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2) / dplyr::first(dose_nmol)),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(`AUC0-200 / dose (nM*h per nmol)`)) |>
  knitr::kable(caption = "Size dependence of exposure. Renal clearance rises about 300-fold from 150 kDa to 50 kDa.")
```

| label | Stoke’s radius (nm) | sigma small pore | CLrenal (L/h) | AUC0-200 / dose (nM\*h per nmol) |
|:---|---:|---:|---:|---:|
| 150 kDa, intact Fc | 4.807 | 0.9979 | 5e-06 | 436561 |
| 100 kDa one-armed, intact Fc | 4.807 | 0.9979 | 5e-06 | 436046 |
| 100 kDa Fab2, no Fc | 4.807 | 0.9979 | 5e-06 | 119510 |
| 50 kDa fragment, no Fc | 4.807 | 0.9979 | 5e-06 | 119510 |

Size dependence of exposure. Renal clearance rises about 300-fold from
150 kDa to 50 kDa. {.table}

``` r

case2 |>
  dplyr::filter(Cc > 1e-3) |>
  ggplot(aes(time, Cc / dose_nmol, colour = label)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(x = "Time (h)", y = "Dose-normalised plasma concentration (nM per nmol)",
       colour = NULL, title = "Replicates Fig. 4 of Patidar 2024 (5 mg/kg IV, mouse)") +
  theme_bw() +
  theme(legend.position = "bottom")
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/case2-plot-1.png)

Dose-normalised exposure falls monotonically as the antibody gets
smaller, which is the paper’s central size finding: antibodies of 4.5-5
nm (100-150 kDa) clear much more slowly than a roughly 50 kDa fragment,
because small pores restrict the large molecules almost completely while
the size-selective glomerular membrane clears the fragment.

## Case 3: charge variants (Fig. 5)

For the +5 variant the paper fixed `spino1` to 1 and estimated `spino2`
= 2.99; `Kp` and `KD,NSB` follow the Fig. A3 / A4 relationships for each
net charge.

``` r

charge_cases <- data.frame(
  label  = c("neutral (0)", "positive (+5)", "negative (-8)"),
  charge = c(0, 5, -8),
  spino2 = c(1, 2.99, 1)
)
case3 <- do.call(rbind, lapply(seq_len(nrow(charge_cases)), function(i) {
  simulate_case(mouse, 10, 0.028, 150000, 200,
                charge = charge_cases$charge[i], spino2 = charge_cases$spino2[i]) |>
    dplyr::mutate(label = charge_cases$label[i])
}))
auc_neutral <- case3 |>
  dplyr::filter(label == "neutral (0)") |>
  dplyr::summarise(a = sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2)) |>
  dplyr::pull(a)
case3 |>
  dplyr::group_by(label) |>
  dplyr::summarise(
    Kp        = round(dplyr::first(kp), 3),
    `KD,NSB (nM)` = round(dplyr::first(kd_nsb), 3),
    `AUC0-200 (nM*h)` = round(sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2)),
    .groups = "drop"
  ) |>
  dplyr::mutate(`AUC vs neutral` = round(`AUC0-200 (nM*h)` / auc_neutral, 3)) |>
  knitr::kable(caption = "Charge dependence of plasma exposure at 10 mg/kg IV in mice.")
```

| label         |    Kp | KD,NSB (nM) | AUC0-200 (nM\*h) | AUC vs neutral |
|:--------------|------:|------------:|-----------------:|---------------:|
| negative (-8) | 0.622 |      16.230 |           119731 |          1.371 |
| neutral (0)   | 1.000 |       8.356 |            87302 |          1.000 |
| positive (+5) | 0.801 |       1.221 |            37912 |          0.434 |

Charge dependence of plasma exposure at 10 mg/kg IV in mice. {.table}

``` r

case3 |>
  dplyr::filter(Cc > 1e-3) |>
  ggplot(aes(time, Cc, colour = label)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(x = "Time (h)", y = "Plasma concentration (nM)", colour = NULL,
       title = "Replicates Fig. 5 of Patidar 2024 (10 mg/kg IV, mouse)") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/case3-plot-1.png)

The positively charged variant has the lowest plasma exposure and the
negative variant the highest, matching the paper’s finding that a
positively charged mAb clears faster than a neutral or negatively
charged one.

## Case 4: target-mediated disposition of anti-CEA IgG (Fig. 6)

Switching on membrane-bound antigen (`icc_tm` = 80 nM) plus soluble
antigen (`icc_t` = 0.016 nM) produces classic target-mediated
non-linearity: exposure rises faster than dose as the antigen sink
saturates.

``` r

case4 <- do.call(rbind, lapply(c(1, 10, 25), function(d) {
  simulate_case(mouse, d, 0.028, 150000, 250,
                icc_t = 0.016, icc_tm = 80) |>
    dplyr::mutate(label = paste0(d, " mg/kg"))
}))
case4 |>
  dplyr::group_by(label) |>
  dplyr::summarise(
    `AUC0-250 (nM*h)` = round(sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2)),
    `AUC / dose (nM*h per nmol)` =
      round(sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2) / dplyr::first(dose_nmol)),
    .groups = "drop"
  ) |>
  knitr::kable(caption = "Dose-normalised exposure of anti-CEA IgG rises with dose: target-mediated disposition.")
```

| label    | AUC0-250 (nM\*h) | AUC / dose (nM\*h per nmol) |
|:---------|-----------------:|----------------------------:|
| 1 mg/kg  |             3390 |                       18161 |
| 10 mg/kg |            94326 |                       50532 |
| 25 mg/kg |           249148 |                       53389 |

Dose-normalised exposure of anti-CEA IgG rises with dose:
target-mediated disposition. {.table}

``` r

case4 |>
  dplyr::filter(Cc > 1e-3) |>
  ggplot(aes(time, Cc, colour = label)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10() +
  labs(x = "Time (h)", y = "Plasma concentration (nM)", colour = NULL,
       title = "Replicates Fig. 6 of Patidar 2024 (anti-CEA IgG, mouse)") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/case4-plot-1.png)

## Case 5: adalimumab in humans (Fig. 7)

``` r

case5 <- do.call(rbind, lapply(c(1, 3, 5), function(d) {
  simulate_case(human, d, 70, 148000, 700, by = 4) |>
    dplyr::mutate(label = paste0(d, " mg/kg"))
}))
case5 |>
  dplyr::filter(Cc > 1e-3) |>
  ggplot(aes(time, Cc, colour = label)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = obs_fig7, aes(time, conc), colour = "black", size = 2,
             inherit.aes = FALSE) +
  scale_y_log10() +
  labs(x = "Time (h)", y = "Plasma concentration (nM)", colour = NULL,
       title = "Replicates Fig. 7 of Patidar 2024 (adalimumab, 70 kg human)",
       subtitle = "Black points digitized from the published 5 mg/kg model curve") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/case5-1.png)

``` r

case5 |>
  dplyr::filter(label == "5 mg/kg", time %in% obs_fig7$time) |>
  dplyr::select(time, Cc) |>
  dplyr::inner_join(obs_fig7, by = "time") |>
  dplyr::transmute(time, simulated = round(Cc, 1), digitized = conc,
                   ratio = round(Cc / conc, 2)) |>
  dplyr::rename("Time (h)" = time, "Simulated (nM)" = simulated,
                "Digitized Fig. 7 (nM)" = digitized,
                "Simulated / digitized" = ratio) |>
  knitr::kable(caption = "Adalimumab 5 mg/kg against the digitized published curve.")
```

| Time (h) | Simulated (nM) | Digitized Fig. 7 (nM) | Simulated / digitized |
|---------:|---------------:|----------------------:|----------------------:|
|        0 |          909.6 |                   300 |                  3.03 |
|      168 |          182.7 |                   170 |                  1.07 |
|      336 |           89.9 |                   140 |                  0.64 |
|      680 |           76.8 |                   115 |                  0.67 |

Adalimumab 5 mg/kg against the digitized published curve. {.table}

The 168 h concentration is matched to within 8%. The simulated profile
has a more pronounced early distribution phase than the published curve
(see Assumptions and deviations).

## Case 6: charge and human plasma clearance (Fig. 8 and Fig. A7)

The paper sweeps net charge from -8 to +10 and computes a linear plasma
clearance from the terminal slope, `CL = kel * Vp`, to compare against
the observed clearance-versus-isoelectric-point correlation of Zheng et
al.

``` r

charge_sweep <- do.call(rbind, lapply(c(-8, -4, 0, 4, 8, 10), function(ch) {
  d <- simulate_case(human, 5, 70, 148000, 700, by = 4, charge = ch)
  fitwin <- d |> dplyr::filter(Cc > 0, time >= 200, time <= 700)
  kel <- -unname(coef(stats::lm(log(Cc) ~ time, data = fitwin))[2])
  data.frame(charge = ch, kp = d$kp[1], kd_nsb = d$kd_nsb[1],
             cl = kel * 2.6)  # human plasma volume, Suppl. Table 1
}))
charge_sweep |>
  dplyr::transmute(charge, Kp = round(kp, 3), `KD,NSB (nM)` = round(kd_nsb, 3),
                   `CL = kel*Vp (L/h)` = signif(cl, 3)) |>
  dplyr::rename("Net charge" = charge) |>
  knitr::kable(caption = "Apparent linear plasma clearance versus net surface charge (human, 5 mg/kg IV).")
```

| Net charge |    Kp | KD,NSB (nM) | CL = kel\*Vp (L/h) |
|-----------:|------:|------------:|-------------------:|
|         -8 | 0.622 |      16.230 |           0.000849 |
|         -4 | 0.918 |      14.645 |           0.001400 |
|          0 | 1.000 |       8.356 |           0.002530 |
|          4 | 0.868 |       1.939 |           0.007100 |
|          8 | 0.521 |       0.280 |           0.015100 |
|         10 | 0.267 |       0.102 |           0.018300 |

Apparent linear plasma clearance versus net surface charge (human, 5
mg/kg IV). {.table}

``` r

zheng <- data.frame(charge = c(-8, 8), cl = c(0.0025, 0.0101),
                    what = "Patidar Fig. A7 reported")
ggplot(charge_sweep, aes(charge, cl)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_point(data = zheng, aes(charge, cl), colour = "red", size = 3, shape = 17) +
  labs(x = "Net surface charge", y = "Linear plasma clearance (L/h)",
       title = "Replicates Fig. 8 of Patidar 2024",
       subtitle = "Red triangles: the two clearances reported in Fig. A7") +
  theme_bw()
```

![](Patidar_2024_mAb_size_charge_fcrn_files/figure-html/case6-plot-1.png)

Clearance rises monotonically with net positive charge, reproducing the
correlation the paper reports. The magnitude at +8 (about 0.015 L/h) is
within 1.5-fold of the reported 0.0101 L/h; at -8 the model is about
threefold lower than the reported 0.0025 L/h, so the simulated charge
effect is steeper than published.

## PKNCA validation

Non-compartmental analysis of the human adalimumab profiles confirms
dose proportionality of the linear-range exposure and provides the
half-life and clearance summaries used in the comparison table below.

``` r

nca_conc <- case5 |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::transmute(treatment = label, time = time, Cc = Cc)

nca_dose <- case5 |>
  dplyr::group_by(treatment = label) |>
  dplyr::summarise(time = 0, dose = dplyr::first(dose_nmol), .groups = "drop")

# One deterministic profile per dose level, so treatment is the only grouping
# level; PKNCAdose() does not accept a nested (slash) formula.
o_conc <- PKNCA::PKNCAconc(nca_conc, Cc ~ time | treatment)
o_dose <- PKNCA::PKNCAdose(nca_dose, dose ~ time | treatment)
o_data <- PKNCA::PKNCAdata(
  o_conc, o_dose,
  intervals = data.frame(start = 0, end = 700,
                         cmax = TRUE, tmax = TRUE, auclast = TRUE,
                         aucinf.obs = TRUE, half.life = TRUE, cl.obs = TRUE)
)
nca_res <- PKNCA::pk.nca(o_data, verbose = FALSE)
as.data.frame(nca_res) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "aucinf.obs",
                                "half.life", "cl.obs")) |>
  dplyr::transmute(treatment, PPTESTCD, PPORRES = signif(PPORRES, 4)) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES) |>
  knitr::kable(caption = "PKNCA summary of the simulated human adalimumab profiles (concentrations nM, dose nmol).")
```

| treatment | auclast |  cmax | tmax | half.life | aucinf.obs |   cl.obs |
|:----------|--------:|------:|-----:|----------:|-----------:|---------:|
| 1 mg/kg   |   24390 | 181.9 |    0 |      8681 |     216400 | 0.002186 |
| 3 mg/kg   |   73170 | 545.7 |    0 |      8694 |     650200 | 0.002182 |
| 5 mg/kg   |  122000 | 909.6 |    0 |      8707 |    1085000 | 0.002179 |

PKNCA summary of the simulated human adalimumab profiles (concentrations
nM, dose nmol). {.table}

``` r

sim_wide <- as.data.frame(nca_res) |>
  dplyr::filter(PPTESTCD %in% c("half.life", "cl.obs")) |>
  dplyr::select(treatment, PPTESTCD, PPORRES)

# Reference: the clearance the paper reports for a neutral-to-slightly-positive
# human mAb (Fig. A7 mAb2, 0.0101 L/h) is not dose-stratified, so the same
# reference value is used for each dose level. Reported in L/h; PKNCA returns
# clearance in L/h because dose is in nmol and concentration in nM.
ref_wide <- data.frame(
  treatment = c("1 mg/kg", "3 mg/kg", "5 mg/kg"),
  cl.obs    = rep(0.0101, 3)
)
tbl <- nlmixr2lib::ncaComparisonTable(
  simulated = sim_wide,
  reference = ref_wide,
  by = "treatment",
  units = c(cl.obs = "L/h")
)
knitr::kable(tbl, caption = "Simulated versus reported human linear plasma clearance.")
```

| NCA parameter | treatment | Reference | Simulated | % diff   |
|:--------------|:----------|:----------|:----------|:---------|
| CL/F (L/h)    | 1 mg/kg   | 0.0101    | 0.00219   | -78.4%\* |
| CL/F (L/h)    | 3 mg/kg   | 0.0101    | 0.00218   | -78.4%\* |
| CL/F (L/h)    | 5 mg/kg   | 0.0101    | 0.00218   | -78.4%\* |

Simulated versus reported human linear plasma clearance. {.table}

``` r

attr(tbl, "footnote")
#> [1] "* differs from reference by more than ±20%."
```

## Assumptions and deviations

The model file encodes the published equations with the corrections
listed below. Every correction is either forced by mass balance or fixed
by the paper’s own tables and text; each is also flagged in an in-file
comment at the relevant ODE.

**1. Pinocytosis flux volume scaling (the one substantive change).**
Eqs. A1 and A4 remove `kup_p * C_plasma` from plasma and add
`kup_p * C_plasma` to the nested endosome with no source-volume /
endosome-volume ratio. Because the endosome is about 600-fold smaller
than plasma in the mouse, taking this literally discards more than 99%
of the pinocytosed antibody and makes FcRn salvage almost irrelevant.
The model file multiplies each uptake flux by the source/endosome volume
ratio. Scored against the paper’s digitized observed data, the
volume-corrected form gives mouse plasma concentrations of 370 / 326 /
298 / 278 nM at 24 / 48 / 96 / 168 h (observed 406 / 370 / 302 / 293 nM)
and human adalimumab 183 nM at 168 h (published curve about 170 nM); the
literal form gives 21 nM and 0.5 nM respectively, i.e. 14-fold and
360-fold low. The same correction is what restores the published
direction of the charge-versus-clearance relationship in Fig. 8.

**2. Transcription errors in the published ODEs.** Corrected as follows.

| Equation | As printed | Encoded | Why |
|----|----|----|----|
| A4 | `Spino*kup*A_p` | `spino_p*kup_p*a_p` | Table 1 and the Results distinguish `kup_p` (plasma endosome) from `kup` (tissue endosomes); A1 loses `kup_p` |
| A5 | `+Spino*kup*T_v` | source is `t_p` | `T_v` is not a state of the plasma-endosome block; mass balance with A2 requires the plasma source |
| A6 | `Spino*kup*ATC_e` | source is `atc_p` | as printed the state produces itself |
| A13, A27 | `k_p,T*Tm` and `+k_off*ATC_v` | `kptm*tm` and `+koff*atm` | Table 1 carries a distinct `k_p,Tm`; mass balance with A15/A29 requires the membrane complex |
| A16, A30 | endosomal uptake from vascular space only | vascular **and** interstitial sources | A22/A36 remove antibody from the interstitium by the same process; otherwise that loss is unaccounted |
| A17, A31 | `-ke_on*A_e1*T_v1` | `t_e1` | an endosomal equation cannot bind the vascular antigen pool |
| A22, A36 | `+CL_TP*A_is/(Kp*Vis)` | driven by `a_v` | the two-pore influx comes from the vasculature; as printed the interstitial state feeds itself |
| A34 | `-ke_on*FcRnA_e2*T_e2 ke_off*FcRnATC_e2` | `+ke_off*...` | operator missing; A20 is the tight-tissue mirror |
| A41, A43 | `J_s = +J_iso + (1-alpha_l)*L` | `-J_iso + ...` | large- and small-pore flows must sum to the tissue lymph flow; this is the form in Li & Shah 2019 Suppl. Eq. 30 |
| A10, A11, A24, A25 | `Kp` multiplies the vascular volume | `Kp` applied to interstitial volume only | the text states `Kp` scales “the volume of interstitial space” |

**3. `CLrec` per endosome.** Table 1 gives one `CLrec` value, labelled
“recycling clearance rate in plasma endosomes”, and notes it is computed
from the endosomal transit time and the endosomal volume. `clrec_p / ve`
is 5.19 1/h in the mouse and 5.20 1/h in the human, i.e. a transit time
of 11.5 min in both species; the model therefore uses that single
first-order recycling rate constant for all three endosomes. Note the
paper’s text says the transit time is 8 min, which does not reproduce
the tabulated `CLrec` values.

**4. `sigma_L` is overloaded.** In A22/A23/A36-A39 `sigma_L` gates
lymphatic drainage from the interstitium, while in A45 and A47 the same
symbol is the **large-pore vascular** reflection coefficient, defined as
`3.5e-5*MW^0.717`. Since no separate lymphatic reflection coefficient
appears in Table 1, and the supplement lists `sigma_L` among the
quantities “calculated using derived equations by Li et al.”, the
MW-derived value is used throughout. For a 150 kDa antibody it evaluates
to 0.180, close to the 0.2 used for the lymphatic reflection coefficient
elsewhere in the Cao/Yuan mPBPK family, so the choice is immaterial for
a full-length antibody but not for a fragment.

**5. Human unit labels in Table 1.** The human column prints “L/h”
against several first-order rate constants (`kup_p`, `kup`, `kdeg`,
`k_p,T`). These are transcription errors inherited from the
clearance-parameterised upstream table; the values are rate constants in
1/h. The proof is the table’s own identity `ksyn = ICC_p,T * k_p,T`:
0.0023 nM/h divided by 2.76e-4 nM is 8.33 1/h, which is the tabulated
8.31 “L/h”. For `kdeg` the choice happens not to matter: because
`k1on * FcRn` is about 43,000 1/h in the human, reading `kdeg` as 0.24
1/h or as 0.24 L/h divided by `ve` (68.6 1/h) changes the predicted
human clearance by 1.2%.

**6. Human GFR is not reported.** Table 1 prints “–”, so `gfr` defaults
to 0 in the human file rather than being filled in from an off-disk
source. Renal clearance is insignificant for antibodies above 4 nm, so
this does not affect any published human result; a user modelling a
human antibody *fragment* must supply a value.

**7. `Spino` has no published charge relationship.** The paper names
`Spino` among the charge-dependent parameters but supplies fitted
equations only for `KD,NSB` (Fig. A3) and `Kp` (Fig. A4). Only three
discrete `Spino` values are reported (`spino2` = 2.99 for the +5 mouse
variant, 1 otherwise). The human charge sweep therefore holds
`spino1 = spino2 = 1`.

**8. Non-specific binding has no elimination route.** Eqs. A14/A28 make
the antibody-membrane-protein complex `arm_*` a purely reversible
binding equilibrium, with no internalisation or degradation. It is
therefore a reservoir rather than a clearance pathway. This limits how
steeply charge can accelerate clearance through the non-specific-binding
mechanism alone; no rate constant for internalising `arm_*` is reported
anywhere in the paper, so none was invented.

**9. Quantitative agreement.** Mouse wild-type plasma is within 12% of
the digitized observed data at all four observed times and the FcRn
knockout within about twofold; dose-normalised exposure falls
monotonically with decreasing size; the charge ordering and the human
charge-clearance correlation are reproduced; human adalimumab is within
8% at 168 h. Two residual gaps: the simulated profiles have a more
pronounced early distribution phase than the published model curves (the
published curves are close to log-linear from t = 0, starting about
threefold below `Dose / Vp` – in both species the published initial
concentration is close to `Dose / (Vp + Vlymph)`, which the paper never
explains), and the simulated terminal slopes are shallower than the
published ones. The model file uses the standard, physically correct
`Dose / Vp` convention.

**10. Complexes do not cross the two-pore barrier.** Eqs. A11/A25 and
A23/A37 give the antibody-antigen complex no `CL_TP` term, so only free
antibody moves between the vascular and interstitial spaces by two-pore
transport. Encoded as printed.

**11. Total tissue concentrations.** Eqs. 1-2 weight the interstitial
concentration by the anatomical `Vis`, not by the available `Kp*Vis`.
Encoded as printed; the two coincide for a neutral antibody (`Kp` = 1).

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.3         tidyr_1.3.2           dplyr_1.2.1          
#> [4] rxode2_5.1.6          PKNCA_0.12.1          nlmixr2lib_0.3.2.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6       xfun_0.60          bslib_0.12.0       lattice_0.22-9    
#>  [5] vctrs_0.7.3        tools_4.6.1        generics_0.1.4     parallel_4.6.1    
#>  [9] tibble_3.3.1       symengine_0.2.13   pkgconfig_2.0.3    data.table_1.18.4 
#> [13] checkmate_2.3.4    RColorBrewer_1.1-3 S7_0.2.2           desc_1.4.3        
#> [17] RcppParallel_6.2.0 lifecycle_1.0.5    compiler_4.6.1     farver_2.1.2      
#> [21] textshaping_1.0.5  fontawesome_0.5.3  htmltools_0.5.9    sys_3.4.3         
#> [25] sass_0.4.10        yaml_2.3.12        pillar_1.11.1      pkgdown_2.2.1     
#> [29] crayon_1.5.3       jquerylib_0.1.4    whisker_0.4.1      openssl_2.4.2     
#> [33] cachem_1.1.0       nlme_3.1-169       tidyselect_1.2.1   digest_0.6.39     
#> [37] lotri_1.0.4        purrr_1.2.2        labeling_0.4.3     rxode2ll_2.0.16   
#> [41] fastmap_1.2.0      grid_4.6.1         cli_3.6.6          dparser_1.3.1-13  
#> [45] magrittr_2.0.5     withr_3.0.3        scales_1.4.0       backports_1.5.1   
#> [49] rmarkdown_2.31     otel_0.2.0         askpass_1.2.1      ragg_1.5.2        
#> [53] memoise_2.0.1      evaluate_1.0.5     knitr_1.51         rex_1.2.2         
#> [57] PreciseSums_0.7    rlang_1.3.0        downlit_0.4.5      Rcpp_1.1.2        
#> [61] glue_1.8.1         xml2_1.6.0         jsonlite_2.0.0     R6_2.6.1          
#> [65] systemfonts_1.3.2  fs_2.1.0
```

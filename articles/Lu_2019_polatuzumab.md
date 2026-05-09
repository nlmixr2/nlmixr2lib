# Polatuzumab vedotin acMMAE + MMAE integrated population PK model (Lu 2019)

``` r

library(nlmixr2lib)
library(PKNCA)
#> 
#> Attaching package: 'PKNCA'
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(rxode2)
#> rxode2 5.0.2 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tidyr)
library(ggplot2)
```

## Model and source

- Citation: Lu D, Lu T, Gibiansky L, Li X, Li C, Agarwal P, Shemesh CS,
  Shi R, Dere RC, Hirata J, Miles D, Chanu P, Girish S, Jin JY.
  *Integrated Two-Analyte Population Pharmacokinetic Model of
  Polatuzumab Vedotin in Patients With Non-Hodgkin Lymphoma.* CPT
  Pharmacometrics Syst Pharmacol. 2020;9(1):48-59.
- Article: <https://doi.org/10.1002/psp4.12482> (PMID 31749251)
- Supplement: NONMEM control stream + equations (open access; published
  alongside the article).

The model couples two submodels (paper Figure 1):

1.  **acMMAE submodel** – antibody-conjugated MMAE, two-compartment,
    with three parallel elimination pathways from the central
    compartment:

    - `CL_NS` (nonspecific linear clearance) declining slowly with time
      on repeated dosing via a sigmoidal Hill function
      `CL_NS(t) = CL_INF * (1 + CL_INF,EMAX * T50^gamma / (T50^gamma + t^gamma))`,
    - `CL_t` (rapidly-declining linear clearance)
      `CL_t(t) = CL_T0 * exp(-kdes * t)`, and
    - `CL_MM` (Michaelis-Menten clearance)
      `CL_MM = Vmax * V1 / (KM + A1/V1)`.

2.  **MMAE submodel** – unconjugated MMAE, two-compartment, with
    apparent parameters (true fraction of formation cannot be estimated,
    so `CL_MMAE` and `V_MMAE` absorb the formation-fraction scalar). The
    central MMAE compartment is driven by the three acMMAE elimination
    pathways with relative conversion fractions modulated by a
    time-decaying multiplier on FRAC_NS (Lu 2019 Eq. 2):

        KINPUT = FRAC_NS * (CL_NS + FRAC_CLT * CL_t + FRAC_MM * CL_MM) / V1
        FRAC_NS(t) = FRAC_0 * (1 + FRAC_T * exp(-alpha * t))

    plus parallel linear (`CL_MMAE`) and Michaelis-Menten
    (`Vmax_MMAE / KSS`) elimination from the MMAE central compartment.

Modeling is in MMAE-equivalent micrograms (a pola dose of 1.8 mg/kg in
an 80 kg patient corresponds to
`1.8 mg/kg * 80 kg * 1000 ug/mg * 3.65 * 718 / 145001 = 2603 ug` of
acMMAE into the central acMMAE compartment, where 3.65 is the mean
drug-to-antibody ratio, 718 g/mol is the MMAE molecular weight, and
145001 g/mol is the antibody-drug conjugate molecular weight; Methods
section).

## Population

The integrated final model was fit to pooled data from four clinical
studies (Lu 2019 Table S1) in 460 adult patients with non-Hodgkin
lymphoma (NHL):

- **DCS4968g** (NCT01290549, Phase I/Ib) – 95 patients with relapsed/
  refractory (R/R) B-cell NHL or chronic lymphocytic leukaemia,
  single-agent pola at 0.1, 0.25, 0.5, 1.0, 1.8, or 2.4 mg/kg q3w; or
  pola 2.4 mg/kg q3w + rituximab. Used the most-intensive PK sampling
  and contributed \>30 plasma samples / patient.
- **GO27834 / ROMULUS** (NCT01691898, Phase Ib/II) – 142 R/R B-cell NHL
  (DLBCL or FL) patients on pola 1.8 or 2.4 mg/kg q3w + rituximab or
  obinutuzumab.
- **GO29365** (NCT02257567, Phase Ib/II) – 106 R/R FL or DLBCL patients
  on pola 1.8 mg/kg q3w (DLBCL) or q4w (FL) + bendamustine + rituximab
  or obinutuzumab.
- **GO29044** (NCT01992653, Phase Ib/II) – 45 previously untreated DLBCL
  patients on pola 1.0-2.4 mg/kg q3w + rituximab or obinutuzumab +
  cyclophosphamide + doxorubicin.

The pooled dataset contributed 4,215 acMMAE and 4,194 unconjugated MMAE
concentration-time pairs after exclusions (post-first-dose BLQ values:
4.5% of acMMAE and 9.2% of MMAE; observations more than 1,008 hours = 6
weeks since the most recent dose were excluded as not clinically
relevant). Sensitivity analyses were performed at body weight 5th-95th
percentiles of 48.7 to 118 kg (Lu 2019 Figure 3 annotation).

The reference subject for the model’s typical-value parameters is a 75
kg male, R/R, non-Asian, ECOG \>= 1, single-agent (no
rituximab/obinutuzumab), normal hepatic function, with serum albumin 35
g/L, baseline B-cell count of 1 cell/uL, and SPD tumor burden of 5,000
mm^2.

``` r

mod_meta <- rxode2::rxode(readModelDb("Lu_2019_polatuzumab"))
#> ℹ parameter labels from comments will be replaced by 'label()'
str(mod_meta$meta$population, max.level = 1)
#> List of 13
#>  $ n_subjects       : int 460
#>  $ n_studies        : int 4
#>  $ n_observations   : chr "4215 acMMAE + 4194 unconjugated MMAE concentration-time pairs (Lu 2019 Results section first paragraph)"
#>  $ age_range        : chr "Adults with NHL (Lu 2019 does not tabulate age in the main paper; the four constituent studies enrolled adult p"| __truncated__
#>  $ weight_range     : chr "5th-95th percentile 48.7-118 kg per Lu 2019 Figure 3 sensitivity-analysis annotation; reference 75 kg used for "| __truncated__
#>  $ sex_female_pct   : num NA
#>  $ race_ethnicity   : chr "Tested as Asian vs non-Asian indicator only in the final model (RACE_ASIAN); a multiplicative effect on Vc (e_a"| __truncated__
#>  $ disease_state    : chr "Relapsed/refractory or previously untreated B-cell non-Hodgkin lymphoma (NHL): diffuse large B-cell lymphoma (D"| __truncated__
#>  $ dose_range       : chr "Pola 0.1-2.4 mg/kg IV every 3 weeks (Q3W) as monotherapy or in combination with rituximab, obinutuzumab, bendam"| __truncated__
#>  $ regions          : chr "Multi-regional (the four constituent studies were global Phase I/Ib/II trials NCT01290549, NCT01691898, NCT0225"| __truncated__
#>  $ studies          : chr "DCS4968g (NCT01290549, Phase I/Ib single-agent and Pola+R run-in), GO27834 / ROMULUS (NCT01691898, Phase Ib/II "| __truncated__
#>  $ reference_subject: chr "75 kg, ALB 35 g/L, TUMSZ 5000 mm^2 SPD, B-cell 1 cell/uL (so max(1, BLBCELL) = 1), male, R/R, non-Asian, normal"| __truncated__
#>  $ notes            : chr "Population characteristics drawn from Lu 2019 main text Results, Figure 3 sensitivity-analysis annotations, and"| __truncated__
```

## Source trace

Every [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html)
entry in `inst/modeldb/specificDrugs/Lu_2019_polatuzumab.R` carries an
in-file comment pointing to Lu 2019 Table 1, Table 2, or Table S3. The
table below collates them.

| Parameter (nlmixr2lib) | Source quantity | Value | Source location |
|----|----|----|----|
| `lkdes` | kdes | 0.0046 1/hour | Table 1, theta1 |
| `lcl_time` | CL_TIME (initial) | 0.00623 L/hour | Table 1, theta2 |
| `lcl` | CL_SS | 0.0344 L/hour | Table 1, theta3 |
| `lvc` | V1 | 3.15 L | Table 1, theta4 |
| `lvp` | V2 | 3.98 L | Table 1, theta5 |
| `lq` | Q | 0.0145 L/hour | Table 1, theta6 |
| `lvmax` | Vmax (acMMAE) | 0.0203 ng/mL/hr | Table 1, theta7 |
| `lkm_ac` | KM (acMMAE) | 0.604 ng/mL | Table 1, theta8 |
| `clss_emax` | CL_SS,EMAX | 0.223 | Table 1, theta9 |
| `lt50_mo` | T50 (months) | 3.53 months | Table 1, theta10 |
| `gamma_ns` | gamma | 2.27 | Table 1, theta11 |
| `lvc_mmae` | V_MMAE (apparent) | 82.2 L | Table 1, theta12 |
| `lcl_mmae` | CL_MMAE (apparent) | 1.89 L/hour | Table 1, theta13 |
| `lq_mmae` | Q_MMAE (apparent) | 36.3 L/hour | Table 1, theta14 |
| `lvp_mmae` | V2_MMAE (apparent) | 200 L | Table 1, theta15 |
| `lvmax_mmae` | Vmax (MMAE) | 0.0307 ng/mL/hr | Table 1, theta16 |
| `lkss_mmae` | KSS (MMAE) | 0.581 ng/mL | Table 1, theta17 |
| `lfrac_clt` | FRAC_CLT | 3.70 | Table 1, theta18 |
| `lfrac_mm` | FRAC_MM | 2.72 | Table 1, theta19 |
| `lalph_mo` | alpha (1/month) | 0.167 1/month | Table 1, theta20 |
| `frac_t` | FRAC_T | 0.139 | Table 1, theta21 |
| `e_wt_cl` | WT on CL_SS | 0.73 | Table 2, theta22 |
| `e_wt_vc` | WT on Vc, Vp, Q (shared) | 0.50 | Table 2, theta23 |
| `e_sexf_vc` | sex on Vc (1/1.20) | 0.8333 | Table 2, theta24 (inverted) |
| `e_asian_vc` | Asian on Vc | 0.929 | Table 2, theta25 |
| `e_line1l_vc` | Treatment-naive on Vc | 1.20 | Table 2, theta26 |
| `e_sexf_cl` | sex on CL_SS (1/1.10) | 0.9091 | Table 2, theta27 (inverted) |
| `e_alb_cl` | ALB on CL_SS (power) | -0.247 | Table 2, theta28 |
| `e_combo_rg_cl` | R/G combo on CL_SS | 0.844 | Table 2, theta29 |
| `e_blbcell_cl` | B-cell on CL_SS (power) | 0.0212 | Table 2, theta30 |
| `e_tumsz_cl` | TUMSZ on CL_SS (linear) | 0.0521 | Table 2, theta31 |
| `e_line1l_kdes` | Treatment-naive on kdes | 3.38 | Table 2, theta32 |
| `e_combo_rg_kdes` | R/G combo on kdes | 0.932 | Table 2, theta33 |
| `e_line1l_cl_time` | Treatment-naive on CL_TIME | 3.53 | Table 2, theta34 |
| `tmbd50_cl_time` | TUMSZ50 on CL_TIME | 1150 mm^2 | Table 2, theta35 |
| `bcell_thr_cl_time` | B-cell threshold on CL_TIME | 121 cells/uL | Table 2, theta36 |
| `e_blbcell_cl_time` | B-cell power on CL_TIME | 0.578 | Table 2, theta37 |
| `e_wt_frac_mmae` | WT on FRAC_NS (power) | -0.467 | Table 2, theta38 |
| `e_sexf_frac_mmae` | sex on FRAC_NS (1/0.911) | 1.0977 | Table 2, theta39 (inverted) |
| `e_line1l_frac_mmae` | Treatment-naive on FRAC_NS | 0.756 | Table 2, theta40 |
| `e_combo_rg_frac_mmae` | R/G combo on FRAC_NS | 0.709 | Table 2, theta41 |
| `e_hepimp_frac_mmae` | NCI ODWG hep. impairment on FRAC_NS | 1.19 | Table 2, theta42 |
| `e_ecog_ge1_frac_mmae` | ECOG_GE1 on FRAC_NS (1/0.905) | 1.1050 | Table 2, theta43 (inverted) |
| `e_alb_frac_mmae` | ALB on FRAC_NS (power) | -0.613 | Table 2, theta44 |
| `propSd` | sqrt(sigma1^2) | 0.1594 | Table S3, Sigma11 = 0.0254 |
| `propSd_mmae` | sqrt(sigma2^2) | 0.2694 | Table S3, Sigma22 = 0.0726 |
| ODE system (4 states) | Lu 2019 Eq. 1 + supplement |  | Supplement \$DES + Equations |

Inter-individual variability (omega^2 stored as variances on the
log-normal scale; %CV = 100 \* sqrt(exp(omega^2) - 1)):

| Parameter | omega^2 (Table S3) |  %CV |
|-----------|-------------------:|-----:|
| CL_TIME   |               1.89 |  138 |
| CL_SS     |             0.0376 | 19.5 |
| Vc        |             0.0151 | 12.3 |
| Vp        |              0.107 | 32.7 |
| Q         |             0.0538 | 23.2 |
| Vmax (ac) |              0.462 | 67.9 |
| FRAC_NS   |             0.0972 | 31.2 |
| CL_MMAE   |              0.115 | 33.9 |
| V2_MMAE   |             0.0422 | 20.5 |

## Errata

The supplement’s printed *Equations* section contains a sign typo on the
acMMAE peripheral compartment ODE:

> Supplement Equations (as printed): `dA2/dt = -K12 A1 - K21 A2`

This is internally inconsistent with mass balance and with the
supplement’s own NONMEM control stream, which gives:

    DADT(2) = -K21*A(2) + K12*A(1)

i.e., `dA2/dt = +K12 A1 - K21 A2`. The packaged model adopts the
NONMEM-control-stream form (the +K12 sign), which conserves total acMMAE
mass under the linear distribution flux. The published Figure 2(a)
total-clearance curve and the cycle 1 / cycle 6 exposure values in Table
3 are derivable only with this sign convention.

## Virtual cohort

We simulate a 200-subject virtual cohort approximating the Lu 2019 study
population:

``` r

set.seed(20260426L)

n_subj   <- 200L
n_cycles <- 6L
cycle_hr <- 504    # 21-day cycle in hours

# Body weight: log-normal anchored at the published reference 75 kg, with
# spread chosen to reproduce the Figure 3 sensitivity-analysis 5th-95th
# percentile range of 48.7-118 kg.
wt    <- exp(rnorm(n_subj, mean = log(75), sd = 0.20))
wt    <- pmin(pmax(wt, 48.7), 118)

# Demographics distributions: sex 50/50; treatment-naive 10% (most studies
# enroll R/R patients with one Phase Ib/II first-line cohort, GO29044, n=45);
# Asian race 15% (representative of the multi-regional pool); albumin
# centered on 38 g/L with SD 4 (NHL clinical norm); ECOG_GE1 50/50 (a typical
# split in NHL studies); hepatic impairment 5% (rare in B-cell NHL trials);
# combo R/G 60% (most patients receive an anti-CD20 partner across the four
# studies); B-cell counts log-normally distributed centered at the threshold
# 121 cells/uL; tumor SPD log-normal centered on the 5,000 mm^2 reference.
sexf   <- rbinom(n_subj, 1, 0.50)
naive  <- rbinom(n_subj, 1, 0.10)
asian  <- rbinom(n_subj, 1, 0.15)
combo  <- rbinom(n_subj, 1, 0.60)
hepimp <- rbinom(n_subj, 1, 0.05)
ecog1  <- rbinom(n_subj, 1, 0.50)
alb    <- pmin(pmax(rnorm(n_subj, 38, 4), 25), 50)
bcell  <- pmin(pmax(exp(rnorm(n_subj, log(50), 1.2)), 0.5), 800)
tumsz  <- pmin(pmax(exp(rnorm(n_subj, log(5000), 0.7)), 200), 25000)

# MMAE-equivalent dose factor: pola_ug/kg * WT_kg * 3.65 * 718 / 145001
mmae_dose_factor <- 3.65 * 718 / 145001

make_subject <- function(id, wt_kg, sexf_ind, naive_ind, asian_ind,
                         combo_ind, hepimp_ind, ecog1_ind,
                         alb_val, bcell_val, tumsz_val) {
  pola_dose_ug <- 1.8 * 1000 * wt_kg                  # 1.8 mg/kg in ug
  mmae_eq_ug   <- pola_dose_ug * mmae_dose_factor

  ev <- rxode2::et() |>
    rxode2::et(amt = mmae_eq_ug, cmt = "central",
               ii = cycle_hr, addl = n_cycles - 1L, time = 0) |>
    rxode2::et(seq(0.5, n_cycles * cycle_hr, length.out = 100),
               cmt = "Cc") |>
    as.data.frame()
  ev$id         <- id
  ev$WT         <- wt_kg
  ev$SEXF       <- sexf_ind
  ev$LINE_1L    <- naive_ind
  ev$RACE_ASIAN <- asian_ind
  ev$COMBO_RG   <- combo_ind
  ev$HEPIMP    <- hepimp_ind
  ev$ECOG_GE1   <- ecog1_ind
  ev$ALB        <- alb_val
  ev$BLBCELL    <- bcell_val
  ev$TUMSZ      <- tumsz_val
  ev$pola_dose_mg <- pola_dose_ug / 1000
  ev$mmae_eq_ug   <- mmae_eq_ug
  ev
}

events <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
  make_subject(i, wt[i], sexf[i], naive[i], asian[i],
               combo[i], hepimp[i], ecog1[i],
               alb[i], bcell[i], tumsz[i])
}))

stopifnot(length(unique(events$id)) == n_subj)
```

## Typical-subject simulation (Figure 2 reproduction)

A single typical-subject curve replicates Lu 2019 Figure 2: a typical
patient receives 1.8 mg/kg q3w for 6 cycles (the labeled regimen). The
deterministic typical-value path uses
[`rxode2::zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html).

``` r

mod <- rxode2::rxode(readModelDb("Lu_2019_polatuzumab"))
#> ℹ parameter labels from comments will be replaced by 'label()'
typ <- mod |> rxode2::zeroRe()

typ_events <- rxode2::et() |>
  rxode2::et(amt = 1.8 * 1000 * 80 * mmae_dose_factor,
             cmt = "central",
             ii = cycle_hr, addl = n_cycles - 1L, time = 0) |>
  rxode2::et(seq(0.5, n_cycles * cycle_hr, length.out = 1500), cmt = "Cc")

# Lu 2019 Figure 2 caption: "typical patient ... bodyweight = 80 kg ... MMAE
# equivalent dose = 2,600 ug." We adopt the same WT = 80 kg reference subject.
typ_events$WT         <- 80
typ_events$SEXF       <- 0
typ_events$LINE_1L    <- 0
typ_events$RACE_ASIAN <- 0
typ_events$COMBO_RG   <- 0
typ_events$HEPIMP    <- 0
typ_events$ECOG_GE1   <- 1
typ_events$ALB        <- 35
typ_events$BLBCELL    <- 1
typ_events$TUMSZ      <- 5000

typ_sim <- rxode2::rxSolve(typ, typ_events, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalcl_time', 'etalcl', 'etalvc', 'etalvp', 'etalq', 'etalvmax', 'etalfrac_mmae', 'etalcl_mmae', 'etalvp_mmae'
```

``` r

typ_sim |>
  mutate(
    cl_mm = vmax / (km_ac + central / vc),  # CL_MM in L/hour
    cl_total = cl_ns + cl_t + cl_mm
  ) |>
  select(time, cl_ns, cl_t, cl_mm, cl_total) |>
  pivot_longer(c(cl_ns, cl_t, cl_mm, cl_total),
               names_to = "pathway", values_to = "cl_Lhr") |>
  mutate(pathway = factor(pathway,
                          levels = c("cl_ns", "cl_t", "cl_mm", "cl_total"),
                          labels = c("CL_NS", "CL_t", "CL_MM", "Total CL"))) |>
  ggplot(aes(time, cl_Lhr, colour = pathway)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = seq(0, by = cycle_hr, length.out = n_cycles),
             linetype = "dashed", alpha = 0.3) +
  labs(x = "Time (hour)", y = "Clearance of acMMAE (L/hour)",
       colour = NULL,
       caption = "Replicates Lu 2019 Figure 2(a).")
```

![Lu 2019 Figure 2(a) replication: time-dependent change of the three
acMMAE clearance pathways and their total under repeated 1.8 mg/kg q3w
dosing for 6 cycles in a typical 80 kg patient. CL_NS slowly declines
via the Hill function; CL_t rapidly decays toward zero; CL_MM
contributes negligibly at therapeutic doses (overlaps the x-axis). Black
arrow positions mark the dosing times at 0, 504, 1008, 1512, 2016, and
2520
hours.](Lu_2019_polatuzumab_files/figure-html/fig2a-clearances-1.png)

Lu 2019 Figure 2(a) replication: time-dependent change of the three
acMMAE clearance pathways and their total under repeated 1.8 mg/kg q3w
dosing for 6 cycles in a typical 80 kg patient. CL_NS slowly declines
via the Hill function; CL_t rapidly decays toward zero; CL_MM
contributes negligibly at therapeutic doses (overlaps the x-axis). Black
arrow positions mark the dosing times at 0, 504, 1008, 1512, 2016, and
2520 hours.

``` r

typ_sim |>
  mutate(
    cl_mm = vmax / (km_ac + central / vc),
    input_clns_ugh = frac_ns * cl_ns * central / vc,
    input_clt_ugh  = frac_ns * frac_clt * cl_t * central / vc,
    input_clmm_ugh = frac_ns * frac_mm * cl_mm * central / vc,
    input_total_ugh = input_clns_ugh + input_clt_ugh + input_clmm_ugh
  ) |>
  select(time, input_clns_ugh, input_clt_ugh, input_clmm_ugh, input_total_ugh) |>
  pivot_longer(c(input_clns_ugh, input_clt_ugh, input_clmm_ugh, input_total_ugh),
               names_to = "pathway", values_to = "rate_ugh") |>
  mutate(pathway = factor(pathway,
                          levels = c("input_clns_ugh", "input_clt_ugh",
                                     "input_clmm_ugh", "input_total_ugh"),
                          labels = c("CL_NS input", "CL_t input",
                                     "CL_MM input", "Total input"))) |>
  ggplot(aes(time, rate_ugh, colour = pathway)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = seq(0, by = cycle_hr, length.out = n_cycles),
             linetype = "dashed", alpha = 0.3) +
  labs(x = "Time (hour)", y = "Rate of input to MMAE central (ug/hour)",
       colour = NULL,
       caption = "Replicates Lu 2019 Figure 2(b).")
```

![Lu 2019 Figure 2(b) replication: time-dependent input rate (ug/hour)
to the unconjugated MMAE central compartment from each acMMAE
elimination pathway under repeated 1.8 mg/kg q3w dosing for 6 cycles.
The decline of the cumulative input rate per cycle is most pronounced
from cycle 1 to cycle 2, driven by the CL_t
pathway.](Lu_2019_polatuzumab_files/figure-html/fig2b-input-rates-1.png)

Lu 2019 Figure 2(b) replication: time-dependent input rate (ug/hour) to
the unconjugated MMAE central compartment from each acMMAE elimination
pathway under repeated 1.8 mg/kg q3w dosing for 6 cycles. The decline of
the cumulative input rate per cycle is most pronounced from cycle 1 to
cycle 2, driven by the CL_t pathway.

``` r

typ_sim |>
  select(time, Cc, Cc_mmae) |>
  pivot_longer(c(Cc, Cc_mmae), names_to = "analyte", values_to = "conc") |>
  mutate(analyte = recode(analyte,
                          Cc      = "acMMAE (ng/mL)",
                          Cc_mmae = "Unconjugated MMAE (ng/mL)")) |>
  ggplot(aes(time, conc)) +
  geom_line(colour = "steelblue", linewidth = 0.8) +
  facet_wrap(~ analyte, ncol = 1, scales = "free_y") +
  scale_y_log10() +
  labs(x = "Time (hour)", y = "Concentration (ng/mL)")
```

![Typical-subject acMMAE (top) and unconjugated MMAE (bottom) plasma
concentration-time profiles under 1.8 mg/kg q3w for 6 cycles in an 80 kg
patient. Cycle 6 acMMAE Cmax should approximate the published cycle 6
mean (721 ng/mL) and cycle 6 MMAE Cmax should approximate the published
mean (~2.27 ng/mL); see PKNCA cross-check
below.](Lu_2019_polatuzumab_files/figure-html/fig2-typical-conc-1.png)

Typical-subject acMMAE (top) and unconjugated MMAE (bottom) plasma
concentration-time profiles under 1.8 mg/kg q3w for 6 cycles in an 80 kg
patient. Cycle 6 acMMAE Cmax should approximate the published cycle 6
mean (721 ng/mL) and cycle 6 MMAE Cmax should approximate the published
mean (~2.27 ng/mL); see PKNCA cross-check below.

## Population (VPC-style) simulation

Full population simulation with all between-subject variability and
residual error preserved. Dropping
[`zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html)
enables both etas and proportional residual error from the model file.

``` r

pop_sim <- rxode2::rxSolve(
  mod, events,
  keep = c("WT", "SEXF", "LINE_1L", "RACE_ASIAN", "COMBO_RG",
           "HEPIMP", "ECOG_GE1", "ALB", "BLBCELL", "TUMSZ",
           "pola_dose_mg", "mmae_eq_ug"),
  returnType = "data.frame"
)
```

``` r

pop_sim |>
  filter(!is.na(Cc)) |>
  group_by(time) |>
  summarise(
    p05 = quantile(Cc, 0.05, na.rm = TRUE),
    p50 = quantile(Cc, 0.50, na.rm = TRUE),
    p95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time, p50)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "steelblue", alpha = 0.3) +
  geom_line(colour = "steelblue4", linewidth = 0.8) +
  scale_y_log10() +
  labs(x = "Time (hour)", y = "acMMAE concentration (ng/mL)")
```

![VPC-style 5th / 50th / 95th percentile envelope of simulated acMMAE
concentrations across 200 virtual subjects on 1.8 mg/kg q3w for 6
cycles. Replicates the structural shape of Lu 2019 Figure S1 (acMMAE
pcVPC).](Lu_2019_polatuzumab_files/figure-html/vpc-acmmae-1.png)

VPC-style 5th / 50th / 95th percentile envelope of simulated acMMAE
concentrations across 200 virtual subjects on 1.8 mg/kg q3w for 6
cycles. Replicates the structural shape of Lu 2019 Figure S1 (acMMAE
pcVPC).

``` r

pop_sim |>
  filter(!is.na(Cc_mmae)) |>
  group_by(time) |>
  summarise(
    p05 = quantile(Cc_mmae, 0.05, na.rm = TRUE),
    p50 = quantile(Cc_mmae, 0.50, na.rm = TRUE),
    p95 = quantile(Cc_mmae, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time, p50)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "tomato", alpha = 0.3) +
  geom_line(colour = "tomato4", linewidth = 0.8) +
  scale_y_log10() +
  labs(x = "Time (hour)", y = "Unconjugated MMAE concentration (ng/mL)")
```

![VPC-style 5th / 50th / 95th percentile envelope of simulated
unconjugated MMAE concentrations across 200 virtual subjects on 1.8
mg/kg q3w for 6 cycles. Replicates the structural shape of Lu 2019
Figure S1 (MMAE pcVPC), including the cycle-over-cycle decline in MMAE
Cmax driven by the time-decaying FRAC_NS
multiplier.](Lu_2019_polatuzumab_files/figure-html/vpc-mmae-1.png)

VPC-style 5th / 50th / 95th percentile envelope of simulated
unconjugated MMAE concentrations across 200 virtual subjects on 1.8
mg/kg q3w for 6 cycles. Replicates the structural shape of Lu 2019
Figure S1 (MMAE pcVPC), including the cycle-over-cycle decline in MMAE
Cmax driven by the time-decaying FRAC_NS multiplier.

## PKNCA validation by cycle vs. Lu 2019 Table 3

Lu 2019 Table 3 reports per-cycle mean +/- SD (CV%) for AUC0-tau, Cmax,
and Ctrough at cycles 1, 3, and 6 (and a cycle-30 hypothetical
steady-state) under 1.8 mg/kg q3w. We compute matching NCA on the
simulated population using PKNCA, with each observation assigned to its
cycle interval `[(c - 1) * 504, c * 504)` hours.

``` r

add_cycle <- function(df) {
  df |>
    mutate(cycle_n       = pmin(n_cycles,
                                pmax(1L,
                                     as.integer(floor((time - 1e-9) / cycle_hr) + 1L))),
           cycle         = factor(cycle_n, levels = 1:n_cycles),
           time_in_cycle = time - (cycle_n - 1) * cycle_hr)
}

acmmae_nca <- pop_sim |>
  filter(!is.na(Cc), time > 0) |>
  add_cycle() |>
  filter(cycle_n %in% c(1, 3, 6)) |>
  transmute(id, cycle, time = time_in_cycle, Cc)

mmae_nca <- pop_sim |>
  filter(!is.na(Cc_mmae), time > 0) |>
  add_cycle() |>
  filter(cycle_n %in% c(1, 3, 6)) |>
  transmute(id, cycle, time = time_in_cycle, Cc_mmae)

dose_df <- events |>
  filter(evid == 1) |>
  mutate(cycle_n = pmin(n_cycles,
                        pmax(1L, as.integer(round(time / cycle_hr) + 1L))),
         cycle   = factor(cycle_n, levels = 1:n_cycles),
         time_in_cycle = 0,
         amt = mmae_eq_ug) |>
  filter(cycle_n %in% c(1, 3, 6)) |>
  transmute(id, cycle, time = time_in_cycle, amt)
```

``` r

ac_conc <- PKNCA::PKNCAconc(acmmae_nca, Cc ~ time | cycle + id,
                            concu = "ng/mL", timeu = "hour")
ac_dose <- PKNCA::PKNCAdose(dose_df, amt ~ time | cycle + id,
                            doseu = "ug")

intervals_q3w <- data.frame(
  start    = 0,
  end      = cycle_hr,
  cmax     = TRUE,
  tmax     = TRUE,
  auclast  = TRUE,
  ctrough  = TRUE
)

ac_nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(ac_conc, ac_dose,
                                              intervals = intervals_q3w))
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#>  ■■■■■■■■■■■■■■■■■                 53% |  ETA:  2s
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Warning: Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
ac_nca_long <- as.data.frame(ac_nca_res$result)
```

``` r

mmae_conc <- PKNCA::PKNCAconc(mmae_nca, Cc_mmae ~ time | cycle + id,
                              concu = "ng/mL", timeu = "hour")
mmae_dose <- PKNCA::PKNCAdose(dose_df, amt ~ time | cycle + id,
                              doseu = "ug")
mmae_nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(mmae_conc, mmae_dose,
                                                intervals = intervals_q3w))
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.5) is not allowed
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.333333) is not allowed
#> Warning: Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      91% |  ETA:  0s
#> Warning: Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (15.3535) is not allowed
mmae_nca_long <- as.data.frame(mmae_nca_res$result)
```

### acMMAE: simulated vs. published

``` r

# Reduce per-subject NCA results to mean +/- SD by cycle, in the same units
# as Lu 2019 Table 3: AUC in ng*day/mL (= AUC_ng_h_mL / 24), Cmax in ng/mL,
# Ctrough in ng/mL.
summarize_cycle <- function(df, conv = c(auclast = 1 / 24,
                                          cmax    = 1,
                                          ctrough = 1)) {
  df |>
    filter(PPTESTCD %in% names(conv)) |>
    mutate(value = PPORRES * conv[PPTESTCD]) |>
    group_by(cycle, PPTESTCD) |>
    summarise(mean = mean(value, na.rm = TRUE),
              sd   = sd(value,   na.rm = TRUE),
              .groups = "drop") |>
    pivot_wider(names_from = PPTESTCD,
                values_from = c(mean, sd))
}

ac_summary <- summarize_cycle(ac_nca_long)
ac_summary_pretty <- ac_summary |>
  transmute(
    Cycle               = cycle,
    `AUC (ng*day/mL) sim`   = sprintf("%.0f +/- %.0f", mean_auclast, sd_auclast),
    `Cmax (ng/mL) sim`      = sprintf("%.0f +/- %.0f", mean_cmax,    sd_cmax),
    `Ctrough (ng/mL) sim`   = sprintf("%.1f +/- %.1f", mean_ctrough,    sd_ctrough)
  )

published_ac <- tibble::tibble(
  Cycle                       = factor(c(1, 3, 6), levels = 1:n_cycles),
  `AUC (ng*day/mL) published` = c("2,020 +/- 571",
                                   "2,640 +/- 522",
                                   "2,900 +/- 577"),
  `Cmax (ng/mL) published`    = c("690 +/- 116",
                                   "713 +/- 117",
                                   "721 +/- 118"),
  `Ctrough (ng/mL) published` = c("13.1 +/- 6.58",
                                   "23.7 +/- 10.1",
                                   "29.3 +/- 12.3")
)

knitr::kable(
  dplyr::left_join(ac_summary_pretty, published_ac, by = "Cycle"),
  caption = "acMMAE NCA -- simulated (200 virtual subjects, mixed covariates) vs. Lu 2019 Table 3 (n = 460, 1.8 mg/kg q3w EBE-based)."
)
```

| Cycle | AUC (ng\*day/mL) sim | Cmax (ng/mL) sim | Ctrough (ng/mL) sim | AUC (ng\*day/mL) published | Cmax (ng/mL) published | Ctrough (ng/mL) published |
|:---|:---|:---|:---|:---|:---|:---|
| 1 | NaN +/- NA | 846 +/- 160 | NaN +/- NA | 2,020 +/- 571 | 690 +/- 116 | 13.1 +/- 6.58 |
| 3 | NaN +/- NA | 869 +/- 162 | NaN +/- NA | 2,640 +/- 522 | 713 +/- 117 | 23.7 +/- 10.1 |
| 6 | NaN +/- NA | 670 +/- 107 | 25.8 +/- 13.8 | 2,900 +/- 577 | 721 +/- 118 | 29.3 +/- 12.3 |

acMMAE NCA – simulated (200 virtual subjects, mixed covariates) vs. Lu
2019 Table 3 (n = 460, 1.8 mg/kg q3w EBE-based). {.table}

### Unconjugated MMAE: simulated vs. published

``` r

mmae_summary <- summarize_cycle(mmae_nca_long)
mmae_summary_pretty <- mmae_summary |>
  transmute(
    Cycle                          = cycle,
    `AUC (ng*day/mL) sim`          = sprintf("%.1f +/- %.1f", mean_auclast, sd_auclast),
    `Cmax (ng/mL) sim`             = sprintf("%.2f +/- %.2f", mean_cmax,    sd_cmax),
    `Ctrough (ng/mL) sim`          = sprintf("%.3f +/- %.3f", mean_ctrough,    sd_ctrough)
  )

published_mmae <- tibble::tibble(
  Cycle                         = factor(c(1, 3, 6), levels = 1:n_cycles),
  `AUC (ng*day/mL) published`   = c("36.5 +/- 33.9",
                                     "26.4 +/- 24.5",
                                     "25.1 +/- 19.6"),
  `Cmax (ng/mL) published`      = c("4.08 +/- 3.99",
                                     "2.45 +/- 1.83",
                                     "2.27 +/- 1.57"),
  `Ctrough (ng/mL) published`   = c("0.269 +/- 0.820",
                                     "0.302 +/- 0.707",
                                     "0.308 +/- 0.527")
)

knitr::kable(
  dplyr::left_join(mmae_summary_pretty, published_mmae, by = "Cycle"),
  caption = "Unconjugated MMAE NCA -- simulated (200 virtual subjects, mixed covariates) vs. Lu 2019 Table 3 (n = 460, 1.8 mg/kg q3w EBE-based)."
)
```

| Cycle | AUC (ng\*day/mL) sim | Cmax (ng/mL) sim | Ctrough (ng/mL) sim | AUC (ng\*day/mL) published | Cmax (ng/mL) published | Ctrough (ng/mL) published |
|:---|:---|:---|:---|:---|:---|:---|
| 1 | NaN +/- NA | 6.14 +/- 4.12 | NaN +/- NA | 36.5 +/- 33.9 | 4.08 +/- 3.99 | 0.269 +/- 0.820 |
| 3 | NaN +/- NA | 3.62 +/- 1.86 | NaN +/- NA | 26.4 +/- 24.5 | 2.45 +/- 1.83 | 0.302 +/- 0.707 |
| 6 | NaN +/- NA | 3.40 +/- 1.76 | 0.409 +/- 0.395 | 25.1 +/- 19.6 | 2.27 +/- 1.57 | 0.308 +/- 0.527 |

Unconjugated MMAE NCA – simulated (200 virtual subjects, mixed
covariates) vs. Lu 2019 Table 3 (n = 460, 1.8 mg/kg q3w EBE-based).
{.table}

The simulated NCA values are not expected to match published values
exactly because:

1.  The published Table 3 mean +/- SD are over `n = 460` actual patients
    with empirical Bayes estimates of all etas, not over a 200-subject
    virtual cohort drawn from the prior.
2.  The virtual cohort here uses representative (not paper-replicated)
    distributions of weight, sex, treatment-naive status, R/G
    combination, ECOG, etc. Lu 2019’s analysis cohort weights individual
    patients by their actual covariates; covariates that systematically
    differ from the reference subject (especially WT, ALB, R/G combo,
    treatment-naive status) shift the means of these summaries.
3.  The model omits the eta-on-epsilon residual-error inflation present
    in the original NONMEM control stream (see Assumptions and
    deviations below); residual variability is therefore slightly
    underestimated for highly-variable parameters like CL_t and Cmax.

The cycle-over-cycle structure – acMMAE accumulating roughly 40% from
cycle 1 to cycle 6 while MMAE Cmax declines to roughly 56% of its cycle
1 value – is the published model’s primary qualitative claim and is
expected to be reproduced.

## Assumptions and deviations

- **Sex reference flip (SEXF vs. source SEX).** Lu 2019 NONMEM defines
  `SEX = SEXN - 1` so the source `SEX` indicator is 1 for males. The
  canonical SEXF reverses the value coding (`SEXF = 1 - source_SEX`).
  Multiplicative-ratio effects therefore invert: the paper’s
  `theta24 = 1.20` (V1 male:female) becomes
  `e_sexf_vc = 1/1.20 = 0.8333`; `theta27 = 1.10` becomes
  `e_sexf_cl = 1/1.10 = 0.9091`; `theta39 = 0.911` becomes
  `e_sexf_frac_mmae = 1/0.911 = 1.0977`. Effect magnitude and reference
  subject (male, SEXF = 0) are preserved.

- **ECOG_GE1 reference flip.** Lu 2019 NONMEM defines
  `ECOG0 = 1 if BECOG == 0`; the canonical ECOG_GE1 is the complement
  (`ECOG_GE1 = 1 - source_ECOG0`). The paper’s `theta43 = 0.905`
  (FRAC_NS_ECOG=0 / FRAC_NS_ECOG\>=1) becomes
  `e_ecog_ge1_frac_mmae = 1/0.905 = 1.1050`. The reference subject is
  ECOG_GE1 = 1 (i.e., ECOG \>= 1) – consistent with the package’s
  ECOG_GE1 canonical (reference 0 = “ECOG = 0”).

- **Hepatic-impairment indicator HEPIMP.** A new general-scope canonical
  was registered for the binary NCI ODWG hepatic-impairment indicator
  (`HEPIMP = 1` for groups \>= 2 = mild or worse; reference 0 = normal).
  Source column `BHPTGRPN` (categorical 1-4 with 9999 missing sentinel)
  decomposes to
  `HEPIMP = as.integer(BHPTGRPN > 1.5 & BHPTGRPN != 9999)`.

- **Anti-CD20 combination indicator COMBO_RG.** A new specific-scope
  canonical was registered for the binary “polatuzumab + rituximab OR
  obinutuzumab” indicator. The Lu 2019 NONMEM separately defines
  `RTX = as.integer(COMBO == 1)` and `GA101 = as.integer(COMBO == 2)`
  and applies effects as `theta^(RTX + GA101)`; because RTX and GA101
  are mutually exclusive, RTX+GA101 takes values 0 or 1 and the effect
  collapses to `theta^COMBO_RG`. The packaged model takes COMBO_RG as a
  single indicator; rituximab and obinutuzumab cannot be distinguished.

- **TUMSZ in mm^2 (SPD), not mm.** The canonical TUMSZ register pools
  RECIST sum-of-diameters (mm) and SPD (mm^2 for cHL/NHL) on a single
  column; for Lu 2019 the unit is mm^2 with reference 5,000 mm^2. The
  per-model `covariateData[[TUMSZ]]$units` and `notes` fields document
  the unit choice load-bearingly; users dosing with a TUMSZ in mm need
  to convert to mm^2 (TUMSZ_mm^2 = sum-of-products) before simulation.

- **B-cell threshold logic.** Lu 2019’s CL_SS and CL_TIME B-cell effects
  use floored-power forms `BLBCELL^0.0212` (with BLBCELL floored at 1)
  on CL_SS and `(BLBCELL/121)^0.578` (with the ratio floored at 1) on
  CL_TIME. Both flooring rules are implemented via `max(1, .)` inside
  [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html).
  For B-cell counts at or below the floor, both effects collapse to 1.

- **Eta-on-epsilon residual error not implemented.** The Lu 2019 NONMEM
  uses a per-individual scaling of residual error
  `Y = TY * (1 + EPS * exp(ETA10/11))` with a 2x2 omega block on
  ETA10/ETA11 (Table S3 omega(10,10), omega(10,11), omega(11,11) =
  0.0521, 0.038, 0.0427). nlmixr2lib does not natively support
  eta-on-epsilon; the packaged model collapses to a fixed proportional
  residual error per analyte using the Sigma point estimates (Sigma11 =
  0.0254 -\> propSd = 0.1594; Sigma22 = 0.0726 -\> propSd_mmae =
  0.2694). This slightly under-estimates residual variability for
  individuals with large positive ETA10/ETA11 and over-estimates for
  those with large negative values, but does not bias the typical-value
  parameters or the IIV terms.

- **Time variable.** `time` in rxode2 corresponds to the NONMEM `T` in
  the \$DES block (time since first dose / start of integration). For
  the CL_NS Hill function and FRAC_NS time-decay, repeated dosing of
  multiple patients is handled by setting up dose events at a common
  time origin; `time` is **not** reset at each dose. This matches the Lu
  2019 model’s intent: time-decaying clearance and conversion-fraction
  effects are cumulative across cycles.

- **Reference subject.** All covariate-effect ratios collapse to 1 when
  WT = 75 kg, ALB = 35 g/L, TUMSZ = 5,000 mm^2, BLBCELL = 1 cell/uL,
  SEXF = 0 (male), LINE_1L = 0 (R/R), RACE_ASIAN = 0, COMBO_RG = 0
  (single agent), HEPIMP = 0 (normal), ECOG_GE1 = 1 (ECOG \>= 1), so the
  typical-value parameters in Table 1 directly apply.

- **Race distribution.** Lu 2019 reports only the Asian indicator as
  retained in the final model (theta25). The 200-subject virtual cohort
  uses 15% Asian / 85% non-Asian, with race not finer-grained.

- **Treatment-naive distribution.** In the Lu 2019 cohort, the only
  exclusively first-line study was GO29044 (n = 45 of 460); the rest
  enrolled R/R patients. The virtual cohort uses 10% naive / 90% R/R,
  approximating the pooled mixture.

## Programmatic access

``` r

mod_meta <- rxode2::rxode(readModelDb("Lu_2019_polatuzumab"))
#> ℹ parameter labels from comments will be replaced by 'label()'
str(mod_meta$meta$population, max.level = 1)
#> List of 13
#>  $ n_subjects       : int 460
#>  $ n_studies        : int 4
#>  $ n_observations   : chr "4215 acMMAE + 4194 unconjugated MMAE concentration-time pairs (Lu 2019 Results section first paragraph)"
#>  $ age_range        : chr "Adults with NHL (Lu 2019 does not tabulate age in the main paper; the four constituent studies enrolled adult p"| __truncated__
#>  $ weight_range     : chr "5th-95th percentile 48.7-118 kg per Lu 2019 Figure 3 sensitivity-analysis annotation; reference 75 kg used for "| __truncated__
#>  $ sex_female_pct   : num NA
#>  $ race_ethnicity   : chr "Tested as Asian vs non-Asian indicator only in the final model (RACE_ASIAN); a multiplicative effect on Vc (e_a"| __truncated__
#>  $ disease_state    : chr "Relapsed/refractory or previously untreated B-cell non-Hodgkin lymphoma (NHL): diffuse large B-cell lymphoma (D"| __truncated__
#>  $ dose_range       : chr "Pola 0.1-2.4 mg/kg IV every 3 weeks (Q3W) as monotherapy or in combination with rituximab, obinutuzumab, bendam"| __truncated__
#>  $ regions          : chr "Multi-regional (the four constituent studies were global Phase I/Ib/II trials NCT01290549, NCT01691898, NCT0225"| __truncated__
#>  $ studies          : chr "DCS4968g (NCT01290549, Phase I/Ib single-agent and Pola+R run-in), GO27834 / ROMULUS (NCT01691898, Phase Ib/II "| __truncated__
#>  $ reference_subject: chr "75 kg, ALB 35 g/L, TUMSZ 5000 mm^2 SPD, B-cell 1 cell/uL (so max(1, BLBCELL) = 1), male, R/R, non-Asian, normal"| __truncated__
#>  $ notes            : chr "Population characteristics drawn from Lu 2019 main text Results, Figure 3 sensitivity-analysis annotations, and"| __truncated__
str(mod_meta$meta$covariateData, max.level = 1)
#> List of 10
#>  $ WT        :List of 6
#>  $ SEXF      :List of 6
#>  $ LINE_1L   :List of 6
#>  $ RACE_ASIAN:List of 6
#>  $ ALB       :List of 6
#>  $ BLBCELL   :List of 6
#>  $ TUMSZ     :List of 6
#>  $ ECOG_GE1  :List of 6
#>  $ HEPIMP    :List of 6
#>  $ COMBO_RG  :List of 6
```

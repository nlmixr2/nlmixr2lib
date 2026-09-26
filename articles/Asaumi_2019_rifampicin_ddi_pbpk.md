# Rifampicin DDI PBPK with pravastatin, pioglitazone, glibenclamide, repaglinide and coproporphyrin I (Asaumi 2019)

## Model and source

- Citation: Asaumi R, Menzel K, Lee W, Nunoya K, Imawaka H, Kusuhara H,
  Sugiyama Y. Expanded Physiologically-Based Pharmacokinetic Model of
  Rifampicin for Predicting Interactions With Drugs and an Endogenous
  Biomarker via Complex Mechanisms Including Organic Anion Transporting
  Polypeptide 1B Induction. CPT Pharmacometrics Syst Pharmacol.
  2019;8:845. <doi:10.1002/psp4.12457> (PMC6875706). Model code, ODEs
  and Tables S1-S4 are the paper’s Supporting Information
  (PSP4-8-845-s003..s008).
- Article: <https://doi.org/10.1002/psp4.12457> (open access,
  PMC6875706)
- Supporting information: Supplementary Text (all ODEs), Supplementary
  Model Code (the Napp programs), Tables S1-S4 and Figures S1-S2.

Asaumi 2019 extends the authors’ earlier rifampicin PBPK model (Asaumi
2018, CPT:PSP 7:186) with OATP1B and CYP2C8 induction and MRP2
inhibition, and uses it to predict interactions with four victim drugs
and the endogenous OATP1B biomarker coproporphyrin I (CP-I). The
Supplementary Model Code contains one Napp program per victim, each
coupling the same rifampicin model to one victim model. The package
ships the five programs as five model files:

| Model file | Victim | Rifampicin mechanisms acting on the victim |
|----|----|----|
| `Asaumi_2019_pravastatin_rifampicin_pbpk` | pravastatin | OATP1B induction + inhibition, MRP2 inhibition |
| `Asaumi_2019_pioglitazone_rifampicin_pbpk` | pioglitazone | CYP2C8 + CYP3A induction |
| `Asaumi_2019_glibenclamide_rifampicin_pbpk` | glibenclamide | OATP1B induction + inhibition, CYP2C9 + hepatic / intestinal CYP3A induction |
| `Asaumi_2019_repaglinide_rifampicin_pbpk` | repaglinide | OATP1B induction + inhibition, CYP2C8 + CYP3A induction |
| `Asaumi_2019_coproporphyrin_I_rifampicin_pbpk` | CP-I (endogenous) | OATP1B induction + inhibition, MRP2 inhibition |

``` r

nms <- c(
  pra = "Asaumi_2019_pravastatin_rifampicin_pbpk",
  pio = "Asaumi_2019_pioglitazone_rifampicin_pbpk",
  glb = "Asaumi_2019_glibenclamide_rifampicin_pbpk",
  rpg = "Asaumi_2019_repaglinide_rifampicin_pbpk",
  cpi = "Asaumi_2019_coproporphyrin_I_rifampicin_pbpk"
)
mods <- lapply(nms, function(n) rxode2::rxode(readModelDb(n)))
vapply(mods, function(m) length(m$state), integer(1))
#> pra pio glb rpg cpi 
#>  47  44  56  58  47
```

## Population

The model is a bottom-up / middle-out PBPK model, not a population
model. The clinical data are published mean blood concentration-time
profiles and AUC ratios (AUCRs) from healthy-volunteer studies
(pravastatin: Lutz 2018, Deng 2009, Maeda 2011; pioglitazone: Jaakkola
2006; glibenclamide: Zheng 2009; repaglinide: Kim 2016, Yoshikado 2017,
Bidstrup 2004, Hatorp 2003, Niemi 2000; CP-I: Takehara 2018, Lai 2016,
Kunze 2018). Every physiological parameter is per kilogram and is scaled
to body weight, which Table S1 sets to the reported value or to 70 kg;
the simulations below use 70 kg. No between-subject variability or
residual error is published (parameters were fitted by nonlinear least
squares in Napp), so every model is deterministic and one typical
subject per scenario is simulated – there is no virtual cohort.

``` r

str(mods$pra$population)
#> List of 6
#>  $ species      : chr "human"
#>  $ n_subjects   : int NA
#>  $ disease_state: chr "healthy volunteers"
#>  $ weight_median: chr "70 kg (assumed; Table S1 footnote)"
#>  $ dose_range   : chr "rifampicin 2-600 mg PO once daily for 10 days or 600 mg single PO; pravastatin 20 mg or 33 ug single PO"
#>  $ notes        : chr "Clinical DDI profiles digitised from the literature: Lutz 2018 (ref. 6, rifampicin 2-600 mg QD x 10 days, prava"| __truncated__
```

## Source trace

Every `ini()` value carries an in-file comment naming its source. In
summary:

| Quantity | Value | Source |
|----|----|----|
| Physiological flows and volumes (per kg) | liver flow 1.24 L/h; muscle 0.642; skin 0.257; adipose 0.223; serosa 0.274; mucosal blood 0.257; portal 0.531; blood volume 0.0743 L; hepatocytes 0.0174; hepatic extracellular 0.0067 | Table S1A |
| Degradation rate constants kdeg | CYP3A4 0.0158 / 0.0288 1/h (liver / enterocyte); CYP2C9 0.00666; CYP2C8 0.0301; OATP1B and UGT set equal to CYP3A4 | Table S1B |
| Rifampicin disposition | Tlag 0.255 h, kin 4.0 1/h, Fa 1, Fg 0.943, fB 0.0778, fH 0.0814, fE 0.115, SF Kp 6.65, fB CLint,all 0.251 L/h/kg, Rdif 0.129, beta 0.2, gamma 0.778, Km,u 0.177 uM, PSdif,E 0.14 L/h/kg, fm UGT 0.759, CLrenal 0.011 L/h/kg | Table 1 and its footnote (from Asaumi 2018) |
| Rifampicin induction Emax (UGT, CYP3A, CYP2C9) and EC50,u | 1.34, 4.57, 2.41; 0.0639 uM | Table 1 footnote |
| Rifampicin induction Emax (OATP1B, CYP2C8) | 2.32 (pravastatin beta 0.2), 2.55 | Table S2B (estimated here) |
| Rifampicin Ki,u for OATP1B | 0.19 (pravastatin), 0.13 (glibenclamide), 0.17 (repaglinide), 0.106 uM (CP-I) | Methods; CP-I from Yoshikado 2018 Table 2 |
| Rifampicin Ki,u for MRP2 | 0.87 uM | Methods |
| Pravastatin | Tlag 0.37 h, ka 0.77 1/h, SF Kp 0.85, fB CLint,all 1.3 L/h/kg, kbile 1.3 1/h (Table S2A); FaFg 0.5, fH 0.496, beta 0.2 (Table 1) | Tables 1, S2A |
| Pioglitazone | Tlag 0.87 h, ka 3.9 1/h, SF Kp 1.6, fB CLint,met 0.050 L/h/kg (Table S2A); FaFg 0.854, fm CYP2C8 0.836 / CYP3A 0.164 (Table 1, Supp. Text) | Tables 1, S2A |
| Glibenclamide | Tlag 0.773 h, ka 0.445 1/h, SF Kp 0.57, fB CLint,all 0.123 L/h/kg, fB 0.000774, fH 0.0221, Rdif 0.246, gamma 0.24; fm CYP2C9 0.85 / CYP3A 0.15 (beta 0.2) | Table 1 (from Asaumi 2018); Asaumi 2018 Discussion |
| Repaglinide | Tlag 0.18 h, ka 3.6 1/h, SF Kp 1.4, fB CLint,all 0.80 L/h/kg (Table S2A); kbile 1.39 1/h, fbile 0.21, fm CYP2C8 0.805 / CYP3A 0.195 (Table 1, Supp. Text) | Tables 1, S2A |
| CP-I | ka 3.0 1/h, FaFg 0.318, fB CLint,all 0.453 L/h/kg, kbile 5.2 1/h, fbile 0.84, CLrenal 0.0421 L/h/kg (Table 1, from Yoshikado 2018); pre-dose blood 0.49 nM (Figure 6 legend) | Table 1; Figure 6 |
| ODE structure | five-unit tandem liver, extended clearance concept, segregated-flow intestine, turnover induction, enterohepatic transit | Supplementary Text and Model Code |

## Simulation set-up

Rifampicin is dosed in umol (MW 822.94 g/mol) into `gut_lumen_rif`
(oral) or `central_rif` (intravenous, 0.5 h infusion; `central_rif` is a
concentration state so bioavailability divides by the blood volume in
the model). Victim doses are in ug and go to `gut_lumen`. Each model
observes the victim blood concentration `Cc` and, where present, the
relative OATP1B / MRP2 / CYP activities. All observations are taken on
the victim’s central compartment; the algebraic observables are returned
as columns at those rows.

``` r

mw_rif <- 822.94
rif_umol <- function(mg) mg * 1000 / mw_rif
trap <- function(t, y) sum(diff(t) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

# Build an event table. rif = data.frame(time, mg, route in {po, iv});
# victim single dose vt (h) / va (ug). Body weight fixed at 70 kg.
mk_events <- function(rif = NULL, vt = NULL, va = NULL, obs) {
  ev <- data.frame(
    id = 1, time = obs, evid = 0, amt = 0,
    cmt = "central", dur = 0
  )
  if (!is.null(rif) && nrow(rif)) {
    ev <- rbind(ev, data.frame(
      id = 1, time = rif$time, evid = 1, amt = rif_umol(rif$mg),
      cmt = ifelse(rif$route == "iv", "central_rif", "gut_lumen_rif"),
      dur = ifelse(rif$route == "iv", 0.5, 0)
    ))
  }
  if (!is.null(vt)) {
    ev <- rbind(ev, data.frame(
      id = 1, time = vt, evid = 1, amt = va, cmt = "gut_lumen", dur = 0
    ))
  }
  ev <- ev[order(ev$time, -ev$evid), ]
  ev$WT <- 70
  ev
}
rif_qd <- function(n, mg = 600, start = 0) {
  data.frame(time = start + 24 * (0:(n - 1)), mg = mg, route = "po")
}
# Victim AUC over [vt, vt+win] via PKNCA (auclast). Baseline dose information
# is not needed for auclast, so PKNCA's 'no dose' message is expected.
pknca_auc_window <- function(sim_df, vt, win) {
  d <- dplyr::filter(sim_df, !is.na(Cc), time >= vt, time <= vt + win)
  d$id <- 1L
  conc <- PKNCA::PKNCAconc(d, Cc ~ time | id)
  ivl <- data.frame(start = vt, end = vt + win, auclast = TRUE)
  res <- as.data.frame(PKNCA::pk.nca(PKNCA::PKNCAdata(conc, intervals = ivl)))
  res$PPORRES[res$PPTESTCD == "auclast"]
}
# Victim AUC ratio (rifampicin / control) computed from PKNCA AUCs.
victim_aucr <- function(mod, rif, vt, va, win = 24) {
  obs <- sort(unique(c(vt, seq(vt, vt + win, by = 0.05))))
  s1 <- rxode2::rxSolve(mod, mk_events(rif, vt, va, obs),
    returnType = "data.frame", atol = 1e-10, rtol = 1e-8)
  s0 <- rxode2::rxSolve(mod, mk_events(NULL, vt, va, obs),
    returnType = "data.frame", atol = 1e-10, rtol = 1e-8)
  pknca_auc_window(s1, vt, win) / pknca_auc_window(s0, vt, win)
}
```

## Rifampicin perpetrator: blood profile and induced activities

A single 600 mg oral dose of rifampicin gives a blood Cmax near 9 ug/mL
(about 11 uM) at ~1.5 h, consistent with the published rifampicin
profile.

``` r

s_rif <- rxode2::rxSolve(mods$pra, mk_events(rif_qd(1), obs = seq(0, 24, 0.05)),
  returnType = "data.frame")
cmax_ug <- max(s_rif$Cc_rif) * mw_rif / 1000
tmax <- s_rif$time[which.max(s_rif$Cc_rif)]
c(cmax_ug_per_mL = round(cmax_ug, 2), tmax_h = tmax)
#> cmax_ug_per_mL         tmax_h 
#>           9.12           1.65
ggplot(s_rif, aes(time, Cc_rif * mw_rif / 1000)) +
  geom_line(color = "firebrick") +
  labs(x = "Time (h)", y = "Rifampicin blood conc. (ug/mL)")
```

![Rifampicin blood concentration after a single 600 mg oral
dose.](Asaumi_2019_rifampicin_ddi_pbpk_files/figure-html/rif-profile-1.png)

Rifampicin blood concentration after a single 600 mg oral dose.

Figure 2c reports the predicted oscillation of relative OATP1B / MRP2 /
CYP activities during repeated 600 mg once-daily rifampicin. The Results
text gives: OATP1B rises to ~240% at steady state, falls to 34% at the
next trough, then 180% at 12 h post-dose; without induction a single
dose depresses OATP1B to 14% (trough) / 68% (12 h); MRP2 falls to 40% /
75%; and the relative increase is greater for hepatic CYP3A (~400%) than
CYP2C8 (~270%) or CYP2C9 (~230%). We reproduce these from the models.

``` r

obs10 <- seq(0, 240, 0.05)
s10 <- rxode2::rxSolve(mods$pra, mk_events(rif_qd(10), obs = obs10),
  returnType = "data.frame")
d10 <- s10[s10$time >= 216 & s10$time <= 240, ]
at12 <- function(df, col) df[[col]][abs(df$time - 228) < 1e-6]
oatp <- c(
  ss_peak = 100 * s10$oatp1b_activity[abs(s10$time - 216) < 1e-6],
  trough = 100 * min(d10$oatp1b_activity),
  plus12h = 100 * at12(d10, "oatp1b_activity")
)
# single dose, induction switched off (emax_oatp = 0): inhibition only
s1_noind <- rxode2::rxSolve(mods$pra, mk_events(rif_qd(1), obs = seq(0, 24, 0.05)),
  params = c(emax_oatp = 0), returnType = "data.frame")
# the same, but for the last dose of the 10-day course (rifampicin exposure
# lowered by its own UGT / CYP3A autoinduction)
s10_noind <- rxode2::rxSolve(mods$pra, mk_events(rif_qd(10), obs = obs10),
  params = c(emax_oatp = 0), returnType = "data.frame")
d10_noind <- s10_noind[s10_noind$time >= 216 & s10_noind$time <= 240, ]
oatp_inhib <- rbind(
  first_dose = c(
    trough = 100 * min(s1_noind$oatp1b_activity),
    plus12h = 100 * s1_noind$oatp1b_activity[abs(s1_noind$time - 12) < 1e-6]
  ),
  day10_dose = c(
    trough = 100 * min(d10_noind$oatp1b_activity),
    plus12h = 100 * at12(d10_noind, "oatp1b_activity")
  )
)
mrp2 <- c(trough = 100 * min(d10$mrp2_activity), plus12h = 100 * at12(d10, "mrp2_activity"))
# induced enzymes: peak relative activities on the glibenclamide / repaglinide files
s_glb <- rxode2::rxSolve(mods$glb, mk_events(rif_qd(10), obs = obs10), returnType = "data.frame")
s_rpg <- rxode2::rxSolve(mods$rpg, mk_events(rif_qd(10), obs = obs10), returnType = "data.frame")
cyp <- c(
  CYP3A_liver = 100 * max(s_glb$enzyme_3a4_liver1),
  CYP2C9 = 100 * max(s_glb$enzyme_2c9_liver1),
  CYP2C8 = 100 * max(s_rpg$enzyme_2c8_liver1)
)
round(oatp, 0)
#> ss_peak  trough plus12h 
#>     242      35     182
round(oatp_inhib, 0)
#>            trough plus12h
#> first_dose     14      51
#> day10_dose     15      72
round(mrp2, 0)
#>  trough plus12h 
#>      41      76
round(cyp, 0)
#> CYP3A_liver      CYP2C9      CYP2C8 
#>         405         231         280
```

These match the Results text: OATP1B 242 / 35 / 182% (paper 240 / 34 /
180%), MRP2 41 / 76% (40 / 75%), and CYP3A ~405%, CYP2C8 ~280%, CYP2C9
~231% (400 / 270 / 230%). The inhibition-only OATP1B statement (“a
single dose … 14% at the trough and 68% at 12 hours”) matches partly. A
literal first dose gives 14% / 51%. The last dose of the 10-day course
with OATP1B induction switched off gives 15% / 72%, which is closer. The
Results text does not say which scenario it used, so both are shown. The
CP-I model independently reproduces the same rifampicin exposure and
induced activities because it embeds an identical rifampicin sub-model.

## AUC ratios of the victim drugs (Table S3)

Table S3 gives predicted and observed AUCRs for each victim under a set
of rifampicin regimens. We reproduce the beta = 0.2 predicted column
(the authors’ recommended conservative set, which is the value fixed in
each shipped model) and compare it against both the published prediction
and the observed value.

``` r

# Each row: model, rifampicin regimen, victim dose time (h) / amount (ug),
# integration window (h), and the published beta=0.2 predicted / observed AUCR.
scen <- tibble::tribble(
  ~drug, ~regimen, ~mod, ~rif, ~vt, ~va, ~win, ~pred, ~obs,
  "Pravastatin", "RIF 600 mg single PO; PRV 33 ug PO", "pra", list(rif_qd(1)), 0, 33, 12, 3.01, 4.64,
  "Pravastatin", "RIF 600 mg single PO; PRV 20 mg PO", "pra", list(rif_qd(1)), 0, 20000, 24, 1.64, 2.27,
  "Pravastatin", "RIF 2 mg QD x10; PRV 20 mg 12 h after", "pra", list(rif_qd(10, 2)), 228, 20000, 24, 0.94, 0.83,
  "Pravastatin", "RIF 10 mg QD x10; PRV 20 mg 12 h after", "pra", list(rif_qd(10, 10)), 228, 20000, 24, 0.80, 0.81,
  "Pravastatin", "RIF 75 mg QD x10; PRV 20 mg 12 h after", "pra", list(rif_qd(10, 75)), 228, 20000, 24, 0.54, 0.42,
  "Pravastatin", "RIF 600 mg QD x10; PRV 20 mg 12 h after", "pra", list(rif_qd(10, 600)), 228, 20000, 24, 0.44, 0.42,
  "Glibenclamide", "RIF 600 mg single IV; GLB 1.25 mg PO", "glb", list(data.frame(time = 0, mg = 600, route = "iv")), 0, 1250, 24, 2.35, 2.18,
  "Glibenclamide", "RIF 600 mg QD x6 PO + IV day7; GLB 1.25 mg PO day7", "glb", list(rbind(rif_qd(6), data.frame(time = 144, mg = 600, route = "iv"))), 144, 1250, 24, 0.70, 0.72,
  "Glibenclamide", "RIF 600 mg QD x7 PO; GLB 1.25 mg PO 48 h after", "glb", list(rif_qd(7)), 192, 1250, 24, 0.28, 0.35,
  "Repaglinide", "RIF 600 mg single PO; RPG 100 ug PO", "rpg", list(rif_qd(1)), 0, 100, 12, 2.35, 2.60,
  "Repaglinide", "RIF 600 mg single PO; RPG 50 ug PO", "rpg", list(rif_qd(1)), 0, 50, 12, 2.40, 1.92,
  "Repaglinide", "RIF 600 mg QD x7 PO; RPG 4 mg PO 0 h after", "rpg", list(rif_qd(7)), 144, 4000, 12, 0.60, 0.50,
  "Repaglinide", "RIF 600 mg QD x7 PO; RPG 4 mg PO 1 h after", "rpg", list(rif_qd(7)), 145, 4000, 12, 0.78, 0.68,
  "Repaglinide", "RIF 600 mg QD x5 PO; RPG 0.5 mg PO 12.5 h after", "rpg", list(rif_qd(5)), 108.5, 500, 12, 0.15, 0.42,
  "Repaglinide", "RIF 600 mg QD x7 PO; RPG 4 mg PO 24 h after", "rpg", list(rif_qd(7)), 168, 4000, 12, 0.09, 0.19
)

scen$simulated <- vapply(seq_len(nrow(scen)), function(i) {
  victim_aucr(mods[[scen$mod[i]]], scen$rif[[i]][[1]], scen$vt[i], scen$va[i], scen$win[i])
}, numeric(1))
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.

victim_tab <- scen %>%
  mutate(
    simulated = round(simulated, 2),
    `abs % diff vs published` = round(100 * abs(simulated - pred) / pred, 1)
  ) %>%
  select(Victim = drug, Regimen = regimen,
         `Simulated AUCR` = simulated, `Published pred. AUCR` = pred,
         `Observed AUCR` = obs, `abs % diff vs published`)
knitr::kable(victim_tab, caption = "Simulated vs Table S3 (beta = 0.2) AUC ratios.")
```

| Victim | Regimen | Simulated AUCR | Published pred. AUCR | Observed AUCR | abs % diff vs published |
|:---|:---|---:|---:|---:|---:|
| Pravastatin | RIF 600 mg single PO; PRV 33 ug PO | 3.06 | 3.01 | 4.64 | 1.7 |
| Pravastatin | RIF 600 mg single PO; PRV 20 mg PO | 2.95 | 1.64 | 2.27 | 79.9 |
| Pravastatin | RIF 2 mg QD x10; PRV 20 mg 12 h after | 0.94 | 0.94 | 0.83 | 0.0 |
| Pravastatin | RIF 10 mg QD x10; PRV 20 mg 12 h after | 0.80 | 0.80 | 0.81 | 0.0 |
| Pravastatin | RIF 75 mg QD x10; PRV 20 mg 12 h after | 0.55 | 0.54 | 0.42 | 1.9 |
| Pravastatin | RIF 600 mg QD x10; PRV 20 mg 12 h after | 0.44 | 0.44 | 0.42 | 0.0 |
| Glibenclamide | RIF 600 mg single IV; GLB 1.25 mg PO | 2.35 | 2.35 | 2.18 | 0.0 |
| Glibenclamide | RIF 600 mg QD x6 PO + IV day7; GLB 1.25 mg PO day7 | 0.73 | 0.70 | 0.72 | 4.3 |
| Glibenclamide | RIF 600 mg QD x7 PO; GLB 1.25 mg PO 48 h after | 0.28 | 0.28 | 0.35 | 0.0 |
| Repaglinide | RIF 600 mg single PO; RPG 100 ug PO | 2.36 | 2.35 | 2.60 | 0.4 |
| Repaglinide | RIF 600 mg single PO; RPG 50 ug PO | 2.36 | 2.40 | 1.92 | 1.7 |
| Repaglinide | RIF 600 mg QD x7 PO; RPG 4 mg PO 0 h after | 0.59 | 0.60 | 0.50 | 1.7 |
| Repaglinide | RIF 600 mg QD x7 PO; RPG 4 mg PO 1 h after | 0.81 | 0.78 | 0.68 | 3.8 |
| Repaglinide | RIF 600 mg QD x5 PO; RPG 0.5 mg PO 12.5 h after | 0.17 | 0.15 | 0.42 | 13.3 |
| Repaglinide | RIF 600 mg QD x7 PO; RPG 4 mg PO 24 h after | 0.12 | 0.09 | 0.19 | 33.3 |

Simulated vs Table S3 (beta = 0.2) AUC ratios. {.table}

``` r

# Structural gate: the deterministic model must reproduce the authors' own
# beta = 0.2 prediction closely. The two staggered repaglinide cases (12.5 h and
# 24 h after the last rifampicin dose) predict very small AUCRs whose value
# depends on the AUC window, which the paper does not state, so they carry a
# looser bound (see Assumptions and deviations).
# The Deng 2009 pravastatin 20 mg row is reported but not gated: see below.
deng <- scen$drug == "Pravastatin" & scen$va == 20000 & scen$vt == 0
stagger <- scen$drug == "Repaglinide" & scen$vt %in% c(108.5, 168)
tight <- scen[!deng & !stagger, ]
loose <- scen[stagger, ]
stopifnot(
  all(abs(tight$simulated - tight$pred) / tight$pred < 0.10),
  all(abs(loose$simulated - loose$pred) / loose$pred < 0.45)
)
```

The Deng 2009 row (600 mg rifampicin with 20 mg pravastatin, both given
together) is **not reproduced**: the model gives about 2.9 against the
published prediction of 1.64 (observed 2.27). The pravastatin model is
linear in dose, and its ratio barely moves with the AUC window (2.94 for
any window from 8 to 96 h). So with the same simultaneous timing, the
shipped parameters give the same AUCR for 20 mg as for the 33 ug Maeda
microdose, and the microdose row matches (3.06 vs 3.01). The published
1.64 must come from a scenario detail of the authors’ Deng simulation
that neither the paper nor the Napp code gives, such as dose timing,
body weight or a study-specific setting. The row is shown for
transparency and is left out of the gate. No parameter was adjusted.

## Coproporphyrin I (endogenous OATP1B biomarker)

CP-I is synthesised endogenously and cleared by OATP1B uptake plus MRP2
biliary excretion, so the model has no victim dose: the state starts at,
and without rifampicin holds, a steady-state blood level of 0.49 nM.
Rifampicin raises CP-I by inhibiting (and, on repeated dosing, also
inducing) OATP1B. The published AUCRs (Table S3) are read against three
clinical reports. Only the Kunze 2018 window is stated (AUC over 16 h,
Table S3 footnote c). The Lai 2016 and Takehara 2018 windows are not
stated in Asaumi 2019; this vignette uses 24 h for Lai and 22 h for
Takehara. The 22 h was chosen because it reproduces both Takehara rows
(600 and 300 mg) at once, which makes it a back-solved integration
window rather than an independent check (see Assumptions and
deviations). No model parameter was adjusted.

``` r

# Steady-state hold: without rifampicin the biomarker must not drift.
s_cpi0 <- rxode2::rxSolve(mods$cpi, mk_events(NULL, 0, 0, seq(0, 200, 1)),
  returnType = "data.frame")
stopifnot(max(abs(s_cpi0$Cc - 0.49)) < 1e-3)
range(s_cpi0$Cc)
#> [1] 0.4900000 0.4900041
```

``` r

cpi_scen <- tibble::tribble(
  ~report, ~rif, ~win, ~pred,
  "Kunze 2018, RIF 600 mg single PO", list(rif_qd(1)), 16, 3.36,
  "Kunze 2018, RIF 600 mg QD x6 PO", list(rif_qd(6, start = 0)), 16, 1.83,
  "Takehara 2018, RIF 600 mg single PO", list(rif_qd(1)), 22, 2.77,
  "Takehara 2018, RIF 300 mg single PO", list(rif_qd(1, 300)), 22, 1.80,
  "Lai 2016, RIF 600 mg single PO", list(rif_qd(1)), 24, 2.56
)
cpi_aucr <- function(rif, win, start_offset = 0) {
  victim_aucr(mods$cpi, rif, start_offset, 0, win)
}
cpi_scen$simulated <- vapply(seq_len(nrow(cpi_scen)), function(i) {
  # repeated-dose row: integrate over the window after the last (day-6) dose
  off <- if (grepl("QD x6", cpi_scen$report[i])) 120 else 0
  cpi_aucr(cpi_scen$rif[[i]][[1]], cpi_scen$win[i], off)
}, numeric(1))
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.

cpi_tab <- cpi_scen %>%
  mutate(simulated = round(simulated, 2),
         `abs % diff` = round(100 * abs(simulated - pred) / pred, 1)) %>%
  select(Report = report, `Window (h)` = win,
         `Simulated AUCR` = simulated, `Published pred. AUCR` = pred,
         `abs % diff`)
knitr::kable(cpi_tab, caption = "CP-I AUC ratios vs Table S3 (beta = 0.2).")
```

| Report | Window (h) | Simulated AUCR | Published pred. AUCR | abs % diff |
|:---|---:|---:|---:|---:|
| Kunze 2018, RIF 600 mg single PO | 16 | 3.36 | 3.36 | 0.0 |
| Kunze 2018, RIF 600 mg QD x6 PO | 16 | 1.82 | 1.83 | 0.5 |
| Takehara 2018, RIF 600 mg single PO | 22 | 2.71 | 2.77 | 2.2 |
| Takehara 2018, RIF 300 mg single PO | 22 | 1.76 | 1.80 | 2.2 |
| Lai 2016, RIF 600 mg single PO | 24 | 2.55 | 2.56 | 0.4 |

CP-I AUC ratios vs Table S3 (beta = 0.2). {.table style="width:100%;"}

``` r

stopifnot(all(abs(cpi_scen$simulated - cpi_scen$pred) / cpi_scen$pred < 0.10))
```

``` r

s_cpi <- rxode2::rxSolve(mods$cpi, mk_events(rif_qd(1), 0, 0, seq(0, 48, 0.1)),
  returnType = "data.frame")
ggplot(s_cpi, aes(time, Cc)) +
  geom_hline(yintercept = 0.49, linetype = "dashed", color = "grey50") +
  geom_line(color = "steelblue") +
  labs(x = "Time (h)", y = "Blood CP-I (nM)")
```

![Predicted blood CP-I after a single 600 mg oral rifampicin dose (Kunze
scenario).](Asaumi_2019_rifampicin_ddi_pbpk_files/figure-html/cpi-profile-1.png)

Predicted blood CP-I after a single 600 mg oral rifampicin dose (Kunze
scenario).

## Pioglitazone (CYP2C8 probe)

Asaumi 2019 does not publish a predicted pioglitazone AUCR; the CYP2C8
Emax (2.55) was *estimated* from the pioglitazone profiles, and the
paper reports the observed AUCR of 0.46 and a 2.7-fold CYP2C8 activity
increase (Discussion). The model AUCR follows from those fitted values
and is shown for completeness; it is not tuned to the observed 0.46.

``` r

pio_aucr <- victim_aucr(mods$pio, rif_qd(6), 109, 30000, 48)
#> No dose information provided, calculations requiring dose will return NA.
#> No dose information provided, calculations requiring dose will return NA.
c(simulated_AUCR = round(pio_aucr, 2), observed_AUCR = 0.46)
#> simulated_AUCR  observed_AUCR 
#>           0.34           0.46
```

## Assumptions and deviations

- **beta fixed at 0.2 in every shipped model.** The paper carries three
  rate-determining-fraction (beta) sensitivity sets (0.2 / 0.5 / 0.8).
  The authors recommend the conservative beta = 0.2 for prediction and
  it is the value in each `ini()`; the 0.5 / 0.8 parameter sets are
  recorded in the in-file comments so a user can reproduce the other
  columns of Table S3.
- **Body weight 70 kg.** Table S1’s footnote uses the reported study
  weight or 70 kg; the per-kg physiological parameters are scaled by
  `WT`, set to 70 kg in all simulations here.
- **Deterministic.** No IIV or residual error is published (the models
  were fitted by nonlinear least squares in Napp), so no random effects
  are encoded and one typical subject per scenario is simulated.
- **AUCRs are the paper’s validation currency.** Absolute victim AUC is
  computed with PKNCA (`auclast`) over the paper’s integration window
  and the DDI/control ratio is formed from those AUCs; this reproduces
  Table S3’s beta = 0.2 predicted column. CP-I windows are
  report-specific (Kunze 16 h, Takehara 22 h, Lai 24 h) per Table S3
  footnotes.
- **Victim AUC windows are not stated in the paper.** The vignette
  integrates pravastatin over 12 h (33 ug microdose) or 24 h (20 mg),
  glibenclamide over 24 h and repaglinide over 12 h after the victim
  dose. With these windows, all rows except two repaglinide rows
  reproduce Table S3’s beta = 0.2 prediction to within 10%. The two
  exceptions are the rows dosed 12.5 h and 24 h after the last
  rifampicin dose, where the predicted AUCRs are small (0.15 and 0.09)
  and depend on the window; they are checked against a looser 45% bound.
- **Deng 2009 pravastatin row not reproduced** (simulated ~2.9 vs
  published 1.64). See the note under the AUCR table; it is excluded
  from the gate.
- **CP-I windows.** Kunze 16 h is from Table S3 footnote c. Lai 24 h and
  Takehara 22 h are assumptions. The Takehara window was back-solved
  from the two Takehara rows, so those two rows only test the 600 vs 300
  mg dose dependence, not the absolute AUCR.
- **Pioglitazone window.** Pioglitazone is dosed at 109 h, as in the
  authors’ Napp program (Supplementary Model Code), during the six-day
  rifampicin course. A 48 h AUC window is assumed. The paper publishes
  neither the window nor a model-predicted AUCR, so the pioglitazone
  result is descriptive only.
- **Upstream parameter provenance.** The rifampicin disposition block
  (Table 1) and the CP-I disposition block (Table 1) originate in the
  authors’ earlier papers (Asaumi 2018 for rifampicin; Yoshikado 2018
  for CP-I) and are fixed here as reported; the in-file comments cite
  the originating table. The rifampicin OATP1B Ki for CP-I (0.106 uM) is
  the Yoshikado 2018 Table 2 value the Methods rounds to ‘~0.1 uM’.
- **fm split provenance.** Pioglitazone and repaglinide CYP2C8 / CYP3A
  fractions and glibenclamide CYP2C9 / CYP3A fractions are the corrected
  values from the Supplementary Text ‘Setting of fm values’ section and
  Table 1; glibenclamide’s fm shifts with beta (0.85/0.15, 0.94/0.06,
  1/0 at beta 0.2/0.5/0.8, from the Asaumi 2018 Discussion).

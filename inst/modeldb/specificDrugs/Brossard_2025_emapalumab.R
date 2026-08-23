Brossard_2025_emapalumab <- function() {
  description <- "Two-compartment population PK model for emapalumab with allometric, age and total-bilirubin covariates and an interferon-gamma-dependent clearance component, linked to three parallel indirect-response (turnover) models for CXCL9, soluble IL-2 receptor alpha and ferritin in patients with macrophage activation syndrome associated with Still's disease"
  reference <- paste(
    "Brossard P. (2025) Emapalumab in Patients With Macrophage Activation",
    "Syndrome Associated With Still's Disease: A Population",
    "Pharmacokinetic/Pharmacodynamic Analysis.",
    "Clin Transl Sci 18(2):e70163. doi:10.1111/cts.70163.",
    "The structural PK model reproduced in Table S2 is the final model of the",
    "upstream analysis Brossard P, Laveille C. (2024) Population",
    "Pharmacokinetics of the Anti-Interferon-Gamma Monoclonal Antibody",
    "Emapalumab: An Updated Analysis. Rheumatol Ther 11:869-880,",
    "doi:10.1007/s40744-024-00669-y, whose Results text supplies the",
    "four-point clearance answer key used to recover the unprinted age and",
    "total-bilirubin reference values.",
    sep = " "
  )
  vignette <- "Brossard_2025_emapalumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # CXCL9, soluble IL-2 receptor alpha and ferritin are hyperinflammation
  # biomarkers specific to this paper's HLH / MAS setting; they are declared
  # here rather than registered as canonical compartments.
  paper_specific_compartments <- c("cxcl9", "sil2ra", "ferritin")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with reference weight 70 kg; exponents fixed at 0.75 (clearances, including Q) and 1 (volumes). Table S2 labels every structural typical value 'per 70 kg', which fixes the reference.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL (both the CLL and CLNL components) and Q with exponent +0.188, and on V1 with exponent -0.104; no age effect on V2 (upstream run log, run 013). Normalised to a 25-year reference age. The reference age is NOT printed in either paper; it was recovered by back-solving the upstream four-point clearance answer key and rounded to 25 y - see the vignette Errata and the derivation note in model().",
      source_name        = "AGE"
    ),
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect with exponent +0.162 on the linear clearance CLL and on both volumes V1 and V2, normalised to 12.8 umol/L. Does NOT act on CLNL or on Q (upstream run log, runs 016-019). The reference value is not printed; the back-solve returns a total-bilirubin factor of 1.0004 for the reference patient, identifying the reference as the pooled-dataset median 12.8 umol/L reported identically for primary HLH, MAS and All in Table S1.",
      source_name        = "TBIL"
    ),
    IFNG = list(
      description        = "Total (free plus emapalumab-bound) serum interferon-gamma concentration",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying driver of the target-mediated ('non-linear') clearance component CLNL, as a power function with exponent +0.542 normalised to 1e6 pg/mL. Above roughly 1e4 pg/mL it drives target-mediated disposition; below it emapalumab clearance is effectively linear. Inert for this model's own population, because DIS_MAS = 1 zeroes CLNL outright.",
      source_name        = "IFNG"
    ),
    DIS_MAS = list(
      description        = "Macrophage activation syndrome (secondary HLH in Still's disease) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (primary haemophagocytic lymphohistiocytosis)",
      notes              = "1 = MAS associated with Still's disease, 0 = primary HLH. The coefficient is fixed at -1, which sets the IFN-gamma-dependent clearance CLNL to exactly zero in MAS patients - the paper's central PK finding. This model's own population is DIS_MAS = 1 throughout; set DIS_MAS = 0 to recover the primary-HLH clearance behaviour.",
      source_name        = "MAS"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "emapalumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "emapalumab", units = "mg", specimen = "serum", verified = TRUE),
    cxcl9       = list(analyte = "CXCL9", units = "ng/L", specimen = "serum", verified = TRUE),
    sil2ra      = list(analyte = "soluble interleukin-2 receptor alpha", units = "ng/L", specimen = "serum", verified = TRUE),
    ferritin    = list(analyte = "ferritin", units = "ug/L", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 14,
    n_studies      = 1,
    age_range      = "2.1-25.4 years",
    age_median     = "11.5 years",
    weight_range   = "12.0-68.8 kg",
    weight_median  = "45.5 kg",
    sex_female_pct = 71.4,
    race_ethnicity = c(White = 78.6, Asian = 14.3, Unknown = 7.1),
    disease_state  = "macrophage activation syndrome associated with Still's disease (systemic juvenile idiopathic arthritis / adult-onset Still's disease)",
    dose_range     = "6 mg/kg IV loading dose, then 3 mg/kg every 3 days until Day 15 and twice weekly until Day 28",
    regions        = "North America and Europe",
    notes          = paste(
      "The PK/PD layer was estimated in the 14 MAS patients of NCT03311854.",
      "The structural PK model was developed on a pooled n = 58 dataset",
      "(44 primary HLH from NCT01818492 plus these 14 MAS patients, with",
      "long-term follow-up NCT02069899) and re-estimated on the MAS",
      "population; baseline characteristics are Table S1 of the source.",
      "Only one MAS patient received emapalumab beyond Day 28 (to Day 39)."
    )
  )

  ini({
    # --- Emapalumab PK: structural parameters, typical values per 70 kg -----
    # Table 1 (MAS re-estimate) and Table S2 (final pooled model) agree on
    # every structural value. Table S2 is in turn identical, value for value
    # and RSE for RSE, to the final-model table of the upstream Brossard &
    # Laveille 2024 analysis.
    lcl <- log(0.0143); label("Linear clearance CLL, independent of IFN-gamma (L/h per 70 kg)")                 # Table 1 CLL typical value 0.0143 (95% CI 0.0118-0.0169); Table S2 'CLL, L/h/70 kg'
    lclnl <- log(0.121); label("IFN-gamma-dependent clearance CLNL at IFN-gamma = 1e6 pg/mL (L/h per 70 kg)")   # Table 1 CLNL 'Intercept for IFNg = 1e6' 0.121 (95% CI 0.0926-0.149); Table S2 'CLNL, L/h/70 kg'
    lvc <- log(3.08); label("Central volume V1 (L per 70 kg)")                                                  # Table 1 V1 typical value 3.08 (95% CI 2.61-3.55); Table S2 'V1, L/70 kg'
    lq <- log(0.105); label("Inter-compartmental clearance Q (L/h per 70 kg)")                                  # Table 1 Q typical value 0.105 (95% CI 0.0815-0.128); Table S2 'Q, L/h/70 kg'
    lvp <- log(4.28); label("Peripheral volume V2 (L per 70 kg)")                                               # Table 1 V2 typical value 4.28 (95% CI 3.43-5.12); Table S2 'V2, L/70 kg'

    # --- Emapalumab PK: covariate effects -----------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")                                 # Table S2 CLs_BW, +0.75 (fixed)
    e_wt_vc <- fixed(1); label("Allometric exponent on V1 and V2 (unitless)")                                    # Table S2 Vs_BW, +1 (fixed)
    e_age_cl <- 0.188; label("Age power exponent on CL and Q (unitless)")                                       # Table 1 'CL, Q age effect, exponent' 0.188 (95% CI 0.124-0.251); Table S2 CLs_AGE reports 0.187 for the pooled fit
    e_age_vc <- -0.104; label("Age power exponent on V1 (unitless)")                                            # Table 1 V1 'Age effect, exponent' -0.104 (95% CI -0.156 to -0.0529); Table S2 V1_AGE
    e_tbili <- 0.162; label("Total-bilirubin power exponent on CLL, V1 and V2 (unitless)")                      # Table 1 'CL, V1, V2 bilirubin effect, exponent' 0.162 (95% CI 0.137-0.187); Table S2 CLs_V1_TBIL
    e_ifng_clnl <- 0.542; label("Total IFN-gamma power exponent on CLNL (unitless)")                            # Table 1 'IFNg effect, exponent' 0.542 (95% CI 0.503-0.580); Table S2 CLNL_IFNg
    e_mas_clnl <- fixed(-1); label("MAS effect on CLNL (-1 sets CLNL to zero in MAS patients)")                    # Table S2 CLNL_MAS, -1 (fixed)

    # --- Emapalumab PK: between-subject variability -------------------------
    # Table 1 reports these as omega SDs on the log scale, not variances:
    # CV = sqrt(exp(sd^2) - 1) maps 0.361 / 0.207 / 0.711 onto Table S2's
    # 37.3 / 20.9 / 81.2 CV% exactly, which fixes the scale. ini() takes
    # variances, so each entry below is the published SD squared.
    etalcl ~ 0.130321                                                                                          # Table 1 CL 'Inter-individual variability' 0.361 (95% CI 0.291-0.431) -> 0.361^2
    etalvc + etalvp ~ c(0.042849,
                        0.083744, 0.505521)                                                                    # Table 1 V1 IIV 0.207, V2 IIV 0.711, 'Correlation (V1, V2)' 0.569 -> 0.207^2, 0.569*0.207*0.711, 0.711^2

    # --- Emapalumab PK: residual error --------------------------------------
    propSd <- 0.301; label("Proportional residual error on emapalumab (fraction)")                              # Table 1 RUV 'Additive component, log' 0.301 (95% CI 0.293-0.309); additive-on-log in NONMEM is proportional in nlmixr2's linear space

    # --- PD: CXCL9 turnover model (Table 2, 'Original' estimates) -----------
    lrbase_cxcl9 <- log(8400); label("Baseline CXCL9 (ng/L)")                                                   # Table 2 CXCL9 Baseline typical value 8400 ng/L (bootstrap mean 8010, 95% CI 3420-12,600)
    lkdeg_cxcl9 <- log(0.414); label("CXCL9 degradation rate constant (1/day)")                                 # Table 2 CXCL9 Degradation rate typical value 0.414/day (bootstrap 0.401, 95% CI 0.242-0.559)
    logitimax_cxcl9 <- qlogis(0.983); label("Logit of maximum fractional inhibition of CXCL9 production (unitless)")       # Table 2 CXCL9 Imax typical value 98.3% (bootstrap 98.2%, 95% CI 97.2-99.2)
    propSd_cxcl9 <- 0.464; label("Proportional residual error on CXCL9 (fraction)")                             # Table 2 CXCL9 RUV 'Proportion, %' 46.4% (bootstrap 46.2%, 95% CI 38.7-53.6)

    # --- PD: soluble IL-2 receptor alpha turnover model (Table 2) -----------
    lrbase_sil2ra <- log(6630); label("Baseline soluble IL-2 receptor alpha (ng/L)")                            # Table 2 sIL-2Ra Baseline typical value 6630 ng/L; the abstract, Results text and Figure 1 all quote 6550 - see vignette Errata
    lkdeg_sil2ra <- log(0.112); label("sIL-2Ra degradation rate constant (1/day)")                              # Table 2 sIL-2Ra Degradation rate typical value 0.112/day (bootstrap 0.115, 95% CI 0.0682-0.162)
    logitimax_sil2ra <- qlogis(0.873); label("Logit of maximum fractional inhibition of sIL-2Ra production (unitless)")    # Table 2 sIL-2Ra Imax typical value 87.3% (bootstrap 86.9%, 95% CI 81.7-92.2); the abstract rounds this to 87%
    propSd_sil2ra <- 0.298; label("Proportional residual error on sIL-2Ra (fraction)")                          # Table 2 sIL-2Ra RUV 'Proportion, %' 29.8% (bootstrap 29.1%, 95% CI 21-37.2)

    # --- PD: ferritin turnover model (Table 2) ------------------------------
    lrbase_ferritin <- log(15300); label("Baseline ferritin (ug/L)")                                            # Table 2 Ferritin Baseline typical value 15,300 ug/L (bootstrap 15,500, 95% CI 6600-24,400)
    lkdeg_ferritin <- log(0.207); label("Ferritin degradation rate constant (1/day)")                           # Table 2 Ferritin Degradation rate typical value 0.207/day (bootstrap 0.208, 95% CI 0.165-0.252)
    logitimax_ferritin <- qlogis(0.996); label("Logit of maximum fractional inhibition of ferritin production (unitless)") # Table 2 Ferritin Imax typical value 99.6% (bootstrap 99.6%, 95% CI 99.2-99.9)
    propSd_ferritin <- 0.389; label("Proportional residual error on ferritin (fraction)")                       # Table 2 Ferritin RUV 'Proportion, %' 38.9% (bootstrap 39%, 95% CI 29.4-48.5)

    # --- PD: between-subject variability ------------------------------------
    # Table 2 labels the baseline and degradation-rate IIVs 'exponent'
    # (log-normal) and the Imax IIVs 'logit'. The values are taken as omega
    # SDs, the convention proven exactly against Table 1 / Table S2 for the PK
    # block above; ini() takes variances. Table 2 reports exactly three
    # correlations: CXCL9 baseline with baseline ferritin (0.675), sIL-2Ra
    # baseline with its own Imax (0.805), and ferritin baseline with its own
    # Imax (0.684). All three blocks are positive definite.
    etalrbase_cxcl9 + etalrbase_ferritin + etalogitimax_ferritin ~
      c(1.638400,
        1.296000, 2.250000,
        0.000000, 1.692900, 2.722500)                                                                          # Table 2: CXCL9 baseline IIV 1.28 -> 1.28^2; ferritin baseline IIV 1.5 -> 1.5^2; ferritin Imax IIV 1.65 -> 1.65^2; 0.675*1.28*1.5 = 1.296; 0.684*1.5*1.65 = 1.6929; no CXCL9-baseline/ferritin-Imax correlation is reported
    etalrbase_sil2ra + etalogitimax_sil2ra ~
      c(0.484416,
        0.499209, 0.793881)                                                                                    # Table 2: sIL-2Ra baseline IIV 0.696 -> 0.696^2; sIL-2Ra Imax IIV 0.891 -> 0.891^2; 0.805*0.696*0.891 = 0.499209
    etalogitimax_cxcl9 ~ 1.144900                                                                              # Table 2 CXCL9 Imax 'Interindividual variability, logit' 1.07 (bootstrap 1, 95% CI 0.606-1.4) -> 1.07^2; Table 2 reports no correlation between CXCL9 Imax and any other random effect, so it stays diagonal
    etalkdeg_cxcl9 ~ 0.481636                                                                                  # Table 2 CXCL9 Degradation rate IIV 0.694 -> 0.694^2
    etalkdeg_sil2ra ~ 0.556516                                                                                 # Table 2 sIL-2Ra Degradation rate IIV 0.746 -> 0.746^2
    etalkdeg_ferritin ~ 0.136900                                                                               # Table 2 Ferritin Degradation rate IIV 0.37 -> 0.37^2

    # --- PD: emapalumab potency ---------------------------------------------
    # Brossard 2025 reports Imax for each marker but NO potency (IC50 / EC50)
    # for the emapalumab-concentration-to-inhibition link, and prints no
    # equations at all. Methods 2.2 nevertheless states that "emapalumab
    # concentrations were assumed to affect (inhibit) PD marker production
    # rates", and Figure 1 draws the central compartment acting on all three
    # production arrows, so the link exists but its potency was never
    # published. The value below is imported from the upstream primary-HLH
    # analysis (Jacqmin et al. 2022, Br J Clin Pharmacol 88:2128-2139,
    # Supplementary Table 7: IC50 = 24.8 ng/mL = 0.0248 ug/mL) and is NOT
    # from this paper. Emapalumab troughs in this trial sit about three orders
    # of magnitude above it, so over the observed 0-56 day window this is
    # numerically indistinguishable from fully-saturated inhibition; the two
    # readings differ only on extrapolation to lower doses or to washout.
    # A single shared potency is used for all three markers, with the three
    # published Imax values carrying the marker-specific extent of inhibition.
    # It is held fixed() so it can never be read as a fit to this dataset.
    lec50 <- fixed(log(0.0248)); label("Emapalumab concentration giving half-maximal inhibition of PD marker production (ug/mL)")  # NOT FROM THIS PAPER - Jacqmin 2022 Supplementary Table 7; see vignette Errata
  })

  model({
    # 0. Constants and unit bookkeeping
    hoursPerDay <- 24    # PD degradation rates are published per day; the PK layer runs on hours

    # 1. Derived covariate terms
    #    References: 70 kg (Table S2 labels every typical value 'per 70 kg'),
    #    25 years and 12.8 umol/L total bilirubin. The age and bilirubin
    #    reference values are NOT printed in either publication. They were
    #    recovered from the upstream Brossard & Laveille 2024 Results text,
    #    which reports emapalumab CL = 0.00218, 0.00308, 0.00623 and 0.01718
    #    L/h at total IFN-gamma = 1e3, 1e4, 1e5 and 1e6 pg/mL "for a 1-year-old
    #    patient weighing 10 kg with primary HLH". Taking the 1e6 / 1e3 ratio
    #    cancels the age and weight factors and returns a total-bilirubin
    #    factor of 1.0004, i.e. the reference patient sits exactly at the
    #    bilirubin reference (the pooled median 12.8 umol/L, Table S1). The age
    #    factor then back-solves to (1/25.3)^0.188; rounding the reference to
    #    25 y reproduces all four published clearances to within 0.3%. The
    #    competing "centre at the dataset median age (1.7 y)" reading is off by
    #    83% and is decisively excluded. See the vignette Errata.
    wtCl <- (WT / 70)^e_wt_cl
    wtV <- (WT / 70)^e_wt_vc
    ageCl <- (AGE / 25)^e_age_cl
    ageVc <- (AGE / 25)^e_age_vc
    tbiliEff <- (TBILI / 12.8)^e_tbili

    # 2. Individual PK parameters
    #    Age acts on total CL (upstream run 014 shows age on CLL alone fits
    #    worse); total bilirubin acts on CLL, V1 and V2 only (upstream runs
    #    016-019); neither acts on CLNL beyond the shared age term.
    cll <- exp(lcl) * wtCl * ageCl * tbiliEff
    clnl <- exp(lclnl) * wtCl * ageCl * (IFNG / 1e6)^e_ifng_clnl * (1 + e_mas_clnl * DIS_MAS)
    cl <- (cll + clnl) * exp(etalcl)
    vc <- exp(lvc + etalvc) * wtV * ageVc * tbiliEff
    vp <- exp(lvp + etalvp) * wtV * tbiliEff
    q <- exp(lq) * wtCl * ageCl

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Emapalumab disposition (two-compartment, IV infusion into central)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central / vc

    # 5. PD: individual turnover parameters
    rbase_cxcl9 <- exp(lrbase_cxcl9 + etalrbase_cxcl9)
    kdeg_cxcl9 <- exp(lkdeg_cxcl9 + etalkdeg_cxcl9) / hoursPerDay
    imax_cxcl9 <- expit(logitimax_cxcl9 + etalogitimax_cxcl9)

    rbase_sil2ra <- exp(lrbase_sil2ra + etalrbase_sil2ra)
    kdeg_sil2ra <- exp(lkdeg_sil2ra + etalkdeg_sil2ra) / hoursPerDay
    imax_sil2ra <- expit(logitimax_sil2ra + etalogitimax_sil2ra)

    rbase_ferritin <- exp(lrbase_ferritin + etalrbase_ferritin)
    kdeg_ferritin <- exp(lkdeg_ferritin + etalkdeg_ferritin) / hoursPerDay
    imax_ferritin <- expit(logitimax_ferritin + etalogitimax_ferritin)

    # 6. PD: fractional target engagement driving the inhibition
    ec50 <- exp(lec50)
    inhib <- Cc / (Cc + ec50)

    # 7. PD: three parallel indirect-response (turnover) models
    #    Production is zero-order and inhibited by emapalumab; loss is first
    #    order. Each marker is assumed to be at its own steady-state baseline
    #    when treatment starts (Methods 2.2), so kin = rbase * kdeg. Figure 1
    #    shows the three turnover models acting in parallel off the central
    #    compartment, NOT sequentially (Results: the parallel model fitted
    #    better than the CXCL9 -> sIL-2Ra -> ferritin cascade).
    cxcl9(0) <- rbase_cxcl9
    sil2ra(0) <- rbase_sil2ra
    ferritin(0) <- rbase_ferritin

    d/dt(cxcl9) <- rbase_cxcl9 * kdeg_cxcl9 * (1 - imax_cxcl9 * inhib) - kdeg_cxcl9 * cxcl9
    d/dt(sil2ra) <- rbase_sil2ra * kdeg_sil2ra * (1 - imax_sil2ra * inhib) - kdeg_sil2ra * sil2ra
    d/dt(ferritin) <- rbase_ferritin * kdeg_ferritin * (1 - imax_ferritin * inhib) - kdeg_ferritin * ferritin

    # 8. Observations
    Cc ~ prop(propSd)
    cxcl9 ~ prop(propSd_cxcl9)
    sil2ra ~ prop(propSd_sil2ra)
    ferritin ~ prop(propSd_ferritin)
  })
}

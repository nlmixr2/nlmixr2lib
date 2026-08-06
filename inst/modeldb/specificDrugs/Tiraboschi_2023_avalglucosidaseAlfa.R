Tiraboschi_2023_avalglucosidaseAlfa <- function() {
  description <- "Cyclic (concatenated) three-compartment population PK model for intravenous avalglucosidase alfa in patients with late-onset and infantile-onset Pompe disease (Tiraboschi 2023), with one-way back-redistribution from the second peripheral compartment to central, parallel linear and Michaelis-Menten elimination from the central compartment, and time-varying body-weight allometric scaling on CL, Vc and Vmax."
  reference <- "Tiraboschi G, Marchionni D, Tuffal G, Fabre D, Martinez JM, An Haack K, Miossec P, Kittner B, Daba N, Hurbin F. Population pharmacokinetic modeling and dosing simulation of avalglucosidase alfa for selecting alternative dosing regimen in pediatric patients with late-onset pompe disease. J Pharmacokinet Pharmacodyn. 2023;50:461-474. doi:10.1007/s10928-023-09874-8"
  vignette <- "Tiraboschi_2023_avalglucosidaseAlfa"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    # Free (active) avalglucosidase alfa quantified by fluorometric enzyme
    # activity assay in plasma (Tiraboschi 2023, Bioanalysis). The paper does
    # not identify a biological matrix for either peripheral compartment --
    # they are mathematical distribution compartments of the concatenated
    # structure -- so those entries are unverified.
    central     = list(analyte = "avalglucosidase alfa", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "avalglucosidase alfa", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "avalglucosidase alfa", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying body weight, the only covariate retained in the final model.",
        "Power (allometric) effect on CL, V1 (= Vc) and Vmax with reference weight 70.5 kg,",
        "the median of the weight values in the analysis dataset (Tiraboschi 2023 Table 3 footnotes a-c:",
        "CL = TVCL * (WT/70.5)^theta10, V1 = TVV1 * (WT/70.5)^theta11, Vm = TVVm * (WT/70.5)^theta12).",
        "The paper used both baseline and time-varying weight during covariate screening and retained",
        "the time-varying form because paediatric patients grew during the studies (Results, Model development).",
        "No allometric factor was supported on the peripheral-compartment parameters (Q2, V2, Q3, V3),",
        "which were fixed; see supplementary Table 1 of the source for the model search.",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  # Covariates that Tiraboschi 2023 screened but did NOT retain in the final
  # model. Documented here (not in covariateData) because they are never
  # referenced in model(); see the vignette "Assumptions and deviations".
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous", reference_category = NULL,
      notes = "Screened by stepwise forward-inclusion / backward-elimination; not retained. Age and body weight were correlated, and the retained allometric weight model already captured the paediatric-versus-adult difference (Results, Model development)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male", units = "(binary)", type = "binary",
      reference_category = "0 (male)",
      notes = "Screened as 'gender'; not retained. 46.2% of the pooled cohort were female (Table 2)."
    ),
    CRCL = list(
      description = "Body-size-normalized creatinine clearance (renal function)", units = "mL/min/1.73 m^2",
      type = "continuous", reference_category = NULL,
      notes = paste(
        "The ONLY covariate that met the forward-inclusion criterion (BCLCRN on CL, dOFV = 15), but the authors",
        "deliberately did not retain it: across the observed range 50.2 to 528 mL/min/1.73 m^2 predicted CL moved",
        "only from 0.781 to 0.843 L/h versus 0.812 L/h at the median 161 mL/min/1.73 m^2, i.e. +/-4%, and CLCRN is",
        "not comparable between adults and children < 12 years (Results, Model development).",
        "Reported as CLCRN / BCLCRN in the source.", sep = " "
      )
    ),
    ALB = list(
      description = "Serum albumin", units = "g/L", type = "continuous", reference_category = NULL,
      notes = "Screened at baseline; not retained. Exposure was additionally stratified by albumin < / >= 45 g/L in Table 4 with no meaningful difference."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous", reference_category = NULL,
      notes = "Screened at baseline (reported in IU/L); not retained. Mean (SD) 83.8 (58.6) IU/L (Table 2)."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous", reference_category = NULL,
      notes = "Screened at baseline (reported in IU/L); not retained. Mean (SD) 98.8 (90.6) IU/L (Table 2)."
    ),
    ALP = list(
      description = "Alkaline phosphatase", units = "U/L", type = "continuous", reference_category = NULL,
      notes = "Screened at baseline (reported in IU/L); not retained."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous", reference_category = NULL,
      notes = "Screened at baseline; not retained. Exposure was additionally stratified by bilirubin < / >= 6.8 umol/L in Table 4 with no meaningful difference."
    ),
    CPK = list(
      description = "Creatine kinase", units = "U/L", type = "continuous", reference_category = NULL,
      notes = "Screened at baseline (reported in IU/L); not retained. Mean (SD) 769 (583) IU/L (Table 2); elevated as expected in Pompe disease."
    ),
    ADA_POS = list(
      description = "Anti-drug-antibody-positive status indicator", units = "(binary)", type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "ADA occurrence was investigated as a covariate using a three-tiered assay approach (Bioanalysis; Model development) and was not retained."
    ),
    DIS_IOPD = list(
      description = "Infantile-onset Pompe disease indicator (versus late-onset Pompe disease)",
      units = "(binary)", type = "binary", reference_category = "0 (late-onset Pompe disease, LOPD)",
      notes = paste(
        "Disease type was screened explicitly and did NOT meet the selection criteria -- one of the paper's key findings,",
        "since it means a single allometric model describes both IOPD (n = 16, 1-11 years) and LOPD (n = 75) patients",
        "(Results, Model development). Not a registered canonical in inst/references/covariate-columns.md:",
        "the final model does not use it, so no register entry is proposed by this extraction.", sep = " "
      )
    ),
    PRIOR_ALGLUCOSIDASE = list(
      description = "Prior treatment with alglucosidase alfa indicator", units = "(binary)", type = "binary",
      reference_category = "0 (treatment-naive)",
      notes = paste(
        "Previous alglucosidase alfa treatment status was screened and not retained; 67.0% of the pooled cohort were",
        "pre-treated (Table 2). Not a registered canonical in inst/references/covariate-columns.md: the final model",
        "does not use it, so no register entry is proposed by this extraction.", sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 91L,
    n_observations = 2242L,
    n_studies      = 4L,
    age_range      = "1-78.4 years",
    age_median     = "mean (SD) 39.2 (20.3) years overall; 46 (15.1) years in LOPD, 6.9 (3.2) years in IOPD",
    weight_range   = "9.9-129 kg",
    weight_median  = "70.5 kg (median of the weight values; mean (SD) 67.7 (26.3) kg)",
    sex_female_pct = 46.2,
    race_ethnicity = c(Caucasian = 83.5, Asian = 12.1, Black = 2.2, Other = 2.2),
    disease_state  = "Late-onset Pompe disease (LOPD, n = 75, including 1 adolescent) and infantile-onset Pompe disease (IOPD, n = 16, aged 1-11 years) showing clinical decline and suboptimal response to alglucosidase alfa.",
    dose_range     = "5, 10, 20 or 40 mg/kg intravenously every 2 weeks. Infusions were administered stepwise: 1 mg/kg/h for 30 min, then 3 mg/kg/h for 30 min, then 5 mg/kg/h for 30 min, then 7 mg/kg/h until the planned amount was delivered (total infusion 3.71 h for 20 mg/kg and 6.57 h for 40 mg/kg).",
    regions        = "International multi-regional (4 pooled clinical trials).",
    notes          = paste(
      "Pooled from NCT01898364 (phase 1, LOPD, N = 24), NCT02032524 (phase 1/2, LOPD, N = 19),",
      "NCT02782741 (phase 3 NEO/COMET, treatment-naive LOPD, N = 51) and NCT03019406 (phase 2 Mini-COMET, IOPD, N = 16);",
      "study designs in Tiraboschi 2023 Table 1, baseline demographics in Table 2.",
      "2498 plasma concentrations were collected; 241 below the LLOQ and 15 previously identified outliers were excluded,",
      "leaving 2242 (LOPD 2042, IOPD 200) concentration-time points.",
      "The fluorometric enzyme-activity assay for free avalglucosidase alfa was validated from 0.0125 (LLOQ) to 3.0 ug/mL",
      "for the diluted sample. NONMEM 7.4.1, FOCE with interaction.",
      "Pre-treatment with alglucosidase alfa: 67.0% of the pooled cohort.", sep = " "
    )
  )

  ini({
    # Structural parameters: Tiraboschi 2023 Table 3, "Estimate" column. Typical values
    # apply at the reference body weight of 70.5 kg (Table 3 footnotes a-c).
    # Cross-check of the reference weight and of the three allometric exponents against the
    # paper's own derived typical values at 68.1 kg (Results, Model development):
    #   0.808 * (68.1/70.5)^0.896 = 0.783 L/h   (paper: 0.783 L/h)
    #   3.37  * (68.1/70.5)^0.661 = 3.29  L     (paper: 3.29 L)
    #   12    * (68.1/70.5)^0.463 = 11.8  mg/h  (paper: 11.8 mg/h)
    lcl   <- log(0.808); label("Linear clearance CL (L/h)")                                 # Table 3, CL row (RSE 3.33%, 95% CI 0.755-0.862)
    lvc   <- log(3.37);  label("Central volume of distribution V1 (L)")                     # Table 3, V1 row (RSE 2.21%, 95% CI 3.22-3.52)
    lvmax <- log(12);    label("Maximum Michaelis-Menten elimination rate Vmax (mg/h)")     # Table 3, Vmax row (RSE 4.59%, 95% CI 10.9-13.1)
    lkm   <- log(0.541); label("Michaelis-Menten constant Km (ug/mL)")                      # Table 3, Km row (RSE 4.45%, 95% CI 0.493-0.589)

    # Concatenated (cyclic) distribution structure:
    #   central <-> peripheral1 -> peripheral2 -> central
    # Q2 is the only two-way inter-compartmental clearance; Q3 (peripheral1 -> peripheral2) and
    # Qpc (peripheral2 -> central, the "drug back redistribution") are one-way
    # (Results, Model development). Q2, V2, Q3 and V3 were fixed by the authors "to avoid
    # model overparameterization and non-identifiability"; only Qpc was estimated.
    lq       <- fixed(log(0.254)); label("Two-way inter-compartmental clearance Q2, central <-> peripheral1 (L/h)")   # Table 3, Q2 row (Fixed)
    lvp      <- fixed(log(296));   label("First peripheral volume of distribution V2 (L)")                            # Table 3, V2 row (Fixed)
    lq_p1_p2 <- fixed(log(1.87));  label("One-way inter-compartmental clearance Q3, peripheral1 -> peripheral2 (L/h)") # Table 3, Q3 row (Fixed)
    lvp2     <- fixed(log(1.31));  label("Second peripheral volume of distribution V3 (L)")                           # Table 3, V3 row (Fixed)
    lq_p2_c  <- log(0.0157);       label("One-way back-redistribution clearance Qpc, peripheral2 -> central (L/h)")   # Table 3, Qpc row (RSE 11.6%, 95% CI 0.0121-0.0194)

    # Allometric exponents on the time-varying body weight, reference 70.5 kg.
    # Table 3 rows "Effect of WT on CL / V1 / Vm" (theta10, theta11, theta12 in footnotes a-c).
    # All three were estimated, not fixed at the canonical 0.75 / 1: a sensitivity analysis
    # comparing optimized versus fixed exponents is reported in the Methods.
    e_wt_cl   <- 0.896; label("Allometric exponent of WT/70.5 on linear CL (unitless)")   # Table 3, Effect of WT on CL (RSE 7.91%, 95% CI 0.754-1.04)
    e_wt_vc   <- 0.661; label("Allometric exponent of WT/70.5 on Vc (unitless)")          # Table 3, Effect of WT on V1 (RSE 6.29%, 95% CI 0.578-0.744)
    e_wt_vmax <- 0.463; label("Allometric exponent of WT/70.5 on Vmax (unitless)")        # Table 3, Effect of WT on Vm (RSE 12%, 95% CI 0.352-0.574)

    # Inter-individual variability. The paper states "Log-normal inter-individual variability
    # of the parameters" (Methods, Model development), so the Table 3 "Inter-individual
    # variability omega^2" values ARE the log-scale variances and are used directly.
    # The parenthetical CV% in Table 3 is recovered as sqrt(exp(omega^2) - 1), which
    # confirms the log-normal reading of each value:
    #   CL   0.0907 -> sqrt(exp(0.0907) - 1) = 30.8%  (Table 3: 30.8%)
    #   V1   0.0184 -> 13.6%  (Table 3: 13.6%)
    #   Vmax 0.118  -> 35.4%  (Table 3: 35.4%)
    #   Km   0.243  -> 52.4%  (Table 3: 52.4%)
    #   Qpc  1.23   -> 156%   (Table 3: 156%)
    # Table 3 reports diagonal variances only; no off-diagonal covariances were published,
    # so the etas are encoded as independent.
    etalcl      ~ 0.0907  # Table 3, omega^2 CL   (RSE 18.5%, 95% CI 0.0578-0.124)
    etalvc      ~ 0.0184  # Table 3, omega^2 V1   (RSE 35.2%, 95% CI 0.00569-0.031)
    etalvmax    ~ 0.118   # Table 3, omega^2 Vm   (RSE 25.4%, 95% CI 0.0593-0.177)
    etalkm      ~ 0.243   # Table 3, omega^2 Km   (RSE 27.9%, 95% CI 0.11-0.376)
    etalq_p2_c  ~ 1.23    # Table 3, omega^2 Qpc  (RSE 26.2%, 95% CI 0.599-1.86); the allometric scaling reduced IIV for every parameter except Qpc

    # Residual error: proportional only ("The residual (intra-individual) variability, modeled
    # through a proportional error model was acceptable with a ~35% CV", Results). Table 3
    # reports sigma^2 = 0.12 with CV 34.6%, and sqrt(0.12) = 0.3464 = 34.6%, so the tabulated
    # value is a variance and the proportional SD is its square root.
    propSd <- 0.3464; label("Proportional residual error (fraction)")  # Table 3, residual variability sigma^2 = 0.12 (34.6%)
  })

  model({
    # Individual parameters. Allometric scaling on the time-varying body weight WT with
    # reference weight 70.5 kg (Tiraboschi 2023 Table 3 footnotes a-c). No allometric factor
    # is applied to the peripheral-compartment parameters -- the extensive model search
    # (supplementary Table 1) did not support one.
    cl   <- exp(lcl   + etalcl)   * (WT / 70.5)^e_wt_cl
    vc   <- exp(lvc   + etalvc)   * (WT / 70.5)^e_wt_vc
    vmax <- exp(lvmax + etalvmax) * (WT / 70.5)^e_wt_vmax
    km   <- exp(lkm   + etalkm)

    q       <- exp(lq)
    vp      <- exp(lvp)
    q_p1_p2 <- exp(lq_p1_p2)
    vp2     <- exp(lvp2)
    q_p2_c  <- exp(lq_p2_c + etalq_p2_c)

    # Concatenated three-compartment system with one-way back redistribution from the third
    # (second peripheral) compartment to central, and parallel linear + Michaelis-Menten
    # elimination from central (Methods, Model development; Results, Model development).
    # Amounts in mg and volumes in L, so central/vc is in mg/L = ug/mL, matching Km.
    # The Michaelis-Menten term uses the inlined state expression rather than a named
    # intermediate: in the nlmixr2 model-function form a named intermediate that references an
    # ODE state can silently evaluate to zero inside d/dt(), deleting the whole saturable arm.
    d/dt(central)     <- -(cl / vc) * central -
                          vmax * (central / vc) / (km + central / vc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1 +
                          (q_p2_c / vp2) * peripheral2
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1 -
                          (q_p1_p2 / vp) * peripheral1
    d/dt(peripheral2) <-  (q_p1_p2 / vp) * peripheral1 -
                          (q_p2_c / vp2) * peripheral2

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

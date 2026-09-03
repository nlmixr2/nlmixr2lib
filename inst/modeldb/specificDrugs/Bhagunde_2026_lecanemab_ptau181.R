Bhagunde_2026_lecanemab_ptau181 <- function() {
  description <- paste(
    "Indirect-response PK/PD model for absolute plasma tau phosphorylated",
    "at threonine 181 (p-tau181) in subjects with early Alzheimer's",
    "disease receiving the anti-protofibril monoclonal antibody",
    "lecanemab. Serum lecanemab concentration Cc inhibits p-tau181",
    "production through an Emax function:",
    "dR/dt = Kin * (1 - Emax * Cc / (EC50 + Cc)) - Kout * R, with",
    "Kin = baseline * Kout so the untreated pool sits at its baseline.",
    "Covariate effects act on the baseline only: power terms for body",
    "weight and baseline Mini-Mental State Examination score, and a ratio",
    "for APOE4-carrier status. Interindividual variability is exponential",
    "on all four PD parameters (baseline, Emax, EC50, Kout) with an",
    "estimated baseline ~ Emax correlation, and residual error is",
    "proportional. Fit by NONMEM FOCE-I to 7909 plasma p-tau181",
    "observations from 2179 subjects pooled across the lecanemab phase 2",
    "Study 201 (Core and open-label extension) and phase 3 Clarity AD",
    "(Study 301, Core and open-label extension) trials. The lecanemab",
    "serum-concentration driver is the two-compartment population PK model",
    "of Majid 2024 (zero-order IV input, linear elimination from the",
    "central compartment); its parameters are held fixed here because the",
    "present paper did not re-estimate them.",
    sep = " "
  )
  reference <- paste(
    "Bhagunde P, Penner N, Willis BA, Bell R, Sachdev P, Charil A,",
    "Irizarry MC, Hersch S, Reyderman L (2026).",
    "Pharmacokinetic/pharmacodynamic analyses of plasma pathophysiology",
    "biomarkers in subjects with early Alzheimer's disease following",
    "lecanemab treatment.",
    "Alzheimer's & Dementia: Translational Research & Clinical",
    "Interventions 12(2):e70246. doi:10.1002/trc2.70246.",
    "Lecanemab population PK driver (CL, V1, V2, Q and their covariate",
    "effects) taken from the cited upstream model:",
    "Majid O, Cao Y, Willis BA, et al. (2024) Population pharmacokinetics",
    "and exposure-response analyses of safety (ARIA-E and isolated ARIA-H)",
    "of lecanemab in subjects with early Alzheimer's disease.",
    "CPT Pharmacometrics Syst Pharmacol 13(12):2111-2123.",
    "doi:10.1002/psp4.13224.",
    sep = " "
  )
  vignette <- "Bhagunde_2026_lecanemab_biomarkers"

  paper_specific_compartments <- c("ptau181")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE
    ),
    ptau181 = list(
      analyte = "tau phosphorylated at threonine 181 (p-tau181)",
      units = "pg/mL", specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters twice. (1) On lecanemab CL and V1 as a power term",
        "normalised to 72 kg (Majid 2024 Table 1 equations). (2) On the",
        "baseline plasma p-tau181 as a power term (WT/72)^-0.235",
        "(Bhagunde 2026 Table 1). The 72 kg centring value is the",
        "reference weight printed in the Majid 2024 popPK equations and",
        "is independently recovered to 71.8-72.0 kg from the Bhagunde",
        "2026 Figure 4C body-weight forest panel; see the vignette",
        "Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "BW / Baseline Weight"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power term on lecanemab CL normalised to 43 g/L",
        "(Majid 2024 Table 1 equation).",
        sep = " "
      ),
      source_name        = "ALB"
    ),
    SEXF = list(
      description        = "Sex",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "1 = female, 0 = male. Acts on lecanemab CL (ratio 0.791) and V1",
        "(ratio 0.868) in the Majid 2024 popPK model. Sex was not",
        "retained as a covariate in the Bhagunde 2026 p-tau181 PD model",
        "(Table S4: females vs males on baseline, dOFV 0.165, not",
        "significant).",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = paste(
        "1 = ADA-positive. Ratio 1.13 on lecanemab CL (Majid 2024",
        "Table 1). Time-varying in the source popPK analysis. ADA status",
        "was tested but not retained in the Bhagunde 2026 p-tau181 PD",
        "model (Table S4).",
        sep = " "
      ),
      source_name        = "ADA"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = paste(
        "1 = Japanese. Acts on lecanemab V1 (ratio 0.920) and V2 (ratio",
        "0.671) in the Majid 2024 popPK model. Japanese race was tested",
        "but not retained in the Bhagunde 2026 p-tau181 PD model",
        "(Table S4).",
        sep = " "
      ),
      source_name        = "JPN"
    ),
    SCORE_MMSE = list(
      description        = "Baseline Mini-Mental State Examination total score",
      units              = "(SCORE_MMSE units, 0-30 score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power term on the baseline plasma p-tau181,",
        "(SCORE_MMSE/26)^-0.468 (Bhagunde 2026 Table 1). The paper does",
        "not print the normalisation value; 26 is recovered from the",
        "Bhagunde 2026 Figure 4C GFAP baseline-MMSE forest panel (22",
        "units -> 1.11 and 29 units -> 0.929 under the tabulated GFAP",
        "MMSE power of -0.63 imply 25.96 and 25.80 respectively). See the",
        "vignette Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "BMMSE"
    ),
    APOE4_CARRIER = list(
      description        = "APOE-epsilon4 carrier status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-carrier)",
      notes              = paste(
        "1 = carries at least one APOE-epsilon4 allele (heterozygous or",
        "homozygous), 0 = non-carrier. Applied as a ratio 1.07 on the",
        "baseline plasma p-tau181 (Bhagunde 2026 Table 1, explicitly",
        "labelled 'ratio'). The source paper tests carrier-vs-noncarrier",
        "only and does not separate heterozygotes from homozygotes, so",
        "the binary APOE4_CARRIER encoding is used rather than",
        "APOE4_COUNT.",
        sep = " "
      ),
      source_name        = "APOE4 carrier status"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2179,
    n_studies      = 2,
    n_observations = 7909,
    disease_state  = paste(
      "Early Alzheimer's disease (mild cognitive impairment due to",
      "Alzheimer's disease, or mild Alzheimer's disease dementia) with",
      "confirmed brain amyloid pathology.",
      sep = " "
    ),
    dose_range     = paste(
      "Lecanemab IV 2.5 mg/kg biweekly, 5 mg/kg monthly, 5 mg/kg",
      "biweekly, 10 mg/kg monthly, 10 mg/kg biweekly, or placebo",
      "(Bhagunde 2026 Table S1).",
      sep = " "
    ),
    regions        = "Global (North America, Europe, Asia including Japan)",
    notes          = paste(
      "Pooled Core and open-label-extension data from the lecanemab phase",
      "2 Study 201 and the phase 3 Clarity AD / Study 301 trial. Subject",
      "and observation counts per biomarker are Bhagunde 2026 Table S1.",
      "Plasma p-tau181 was measured with the Quanterix single-molecule",
      "array (Simoa) assay (Bhagunde 2026 Appendix A).",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------- #
    # Lecanemab population PK -- two-compartment, zero-order IV input,  #
    # linear elimination. All values FIXED: not re-estimated by         #
    # Bhagunde 2026, which used post hoc individual PK parameters.      #
    # Source: Majid 2024 Table 1 and the equations printed beneath it.  #
    # ---------------------------------------------------------------- #
    lcl <- fixed(log(0.0154)); label("Lecanemab clearance at the reference subject (CL, L/h)")        # Majid 2024 Table 1: CL = 0.0154 L/h (RSE 1.60%)
    lvc <- fixed(log(3.24));   label("Lecanemab central volume at the reference subject (V1, L)")     # Majid 2024 Table 1: V1 = 3.24 L (RSE 0.799%)
    lvp <- fixed(log(2.00));   label("Lecanemab peripheral volume at the reference subject (V2, L)")  # Majid 2024 Table 1: V2 = 2.00 L (RSE 4.09%)
    lq  <- fixed(log(0.00718)); label("Lecanemab intercompartmental clearance (Q, L/h)")              # Majid 2024 Table 1: Q = 0.00718 L/h (RSE 4.23%)

    e_wt_cl            <- fixed(0.353);  label("Power exponent for body weight on CL, (WT/72)^e (unitless)")            # Majid 2024 Table 1
    e_alb_cl           <- fixed(-0.374); label("Power exponent for serum albumin on CL, (ALB/43)^e (unitless)")         # Majid 2024 Table 1
    e_sexf_cl          <- fixed(0.791);  label("Ratio of CL in females to males, applied as e^SEXF (unitless)")         # Majid 2024 Table 1
    e_ada_pos_cl       <- fixed(1.13);   label("Ratio of CL in ADA-positive to ADA-negative subjects, applied as e^ADA_POS (unitless)")  # Majid 2024 Table 1
    e_wt_vc            <- fixed(0.513);  label("Power exponent for body weight on V1, (WT/72)^e (unitless)")            # Majid 2024 Table 1
    e_sexf_vc          <- fixed(0.868);  label("Ratio of V1 in females to males, applied as e^SEXF (unitless)")         # Majid 2024 Table 1
    e_race_japanese_vc <- fixed(0.920);  label("Ratio of V1 in Japanese to non-Japanese subjects, applied as e^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1
    e_race_japanese_vp <- fixed(0.671);  label("Ratio of V2 in Japanese to non-Japanese subjects, applied as e^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1

    # ---------------------------------------------------------------- #
    # Plasma p-tau181 indirect-response PD model                        #
    # (Bhagunde 2026 Section 2.2.2 equation and Table 1)                #
    # ---------------------------------------------------------------- #
    lrbase <- log(3.40);  label("Baseline plasma p-tau181 at the reference subject (pg/mL)")                    # Bhagunde 2026 Table 1: Baseline p-tau181 = 3.40 pg/mL (RSE 1.44%)
    lemax  <- log(0.480); label("Maximum fractional inhibition of p-tau181 production by lecanemab (Emax, unitless)")  # Bhagunde 2026 Table 1: Emax = 0.480 (RSE 6.67%)
    lec50  <- log(31.4);  label("Lecanemab serum concentration producing half of Emax (EC50, ug/mL)")           # Bhagunde 2026 Table 1: EC50 = 31.4 (RSE 59.2%). The table mg/mL unit label is a typographical error -- see the vignette Errata
    lkout  <- log(0.428); label("First-order degradation rate constant of the plasma p-tau181 pool (Kout, 1/year)")  # Bhagunde 2026 Table 1: Kout = 0.428 /year (RSE 21.9%), t1/2 = ln(2)/0.428 = 1.62 year, matching the paper approximately 1.6 years

    e_wt_rbase            <- -0.235; label("Power exponent for body weight on baseline p-tau181, (WT/72)^e (unitless)")           # Bhagunde 2026 Table 1: Body weight on baseline (exponent) = -0.235 (RSE 17.7%)
    e_mmse_rbase          <- -0.468; label("Power exponent for baseline MMSE on baseline p-tau181, (SCORE_MMSE/26)^e (unitless)") # Bhagunde 2026 Table 1: Baseline MMSE on baseline (exponent) = -0.468 (RSE 21.8%)
    e_apoe4_carrier_rbase <- 1.07;   label("Ratio of baseline p-tau181 in APOE4 carriers to non-carriers, applied as e^APOE4_CARRIER (unitless)")  # Bhagunde 2026 Table 1: APOE4 carrier status (vs noncarrier) on baseline (ratio) = 1.07 (RSE 1.81%)

    # Lecanemab PK interindividual variability, fixed from Majid 2024
    # (omega^2 = log(CV^2 + 1)).
    etalcl + etalvc ~ fixed(c(0.114936,
                              0.005934, 0.014774))   # Majid 2024 Table 1: IIV CL 34.9 %CV, IIV V1 12.2 %CV, correlation CL ~ V1 R = 0.144
    etalvp          ~ fixed(0.639175)                # Majid 2024 Table 1: IIV V2 94.6 %CV

    # Plasma p-tau181 interindividual variability (estimated here)
    etalrbase + etalemax ~ c(0.146357,
                             0.066004, 0.075478)     # Bhagunde 2026 Table 1: IIV baseline 39.7 %CV, IIV Emax 28.0 %CV, correlation baseline ~ Emax R = 0.628
    etalec50 ~ 1.994767                              # Bhagunde 2026 Table 1: IIV EC50 252 %CV
    etalkout ~ 0.336061                              # Bhagunde 2026 Table 1: IIV Kout 63.2 %CV

    propSd <- 0.192; label("Proportional residual error on plasma p-tau181 (fraction)")  # Bhagunde 2026 Table 1: Proportional residual variability 19.2 %CV (RSE 1.25%)
  })

  model({
    # Constants and reference (centring) values
    hrs_per_year <- 8766   # 365.25 days/year x 24 h/day; converts Kout, reported per year, to the model time unit (h)
    wt_ref       <- 72     # kg;  Majid 2024 popPK equations normalise WT to 72 kg
    alb_ref      <- 43     # g/L; Majid 2024 popPK equations normalise ALB to 43 g/L
    mmse_ref     <- 26     # MMSE units; normalisation recovered from Bhagunde 2026 Figure 4C (see covariateData notes)

    # Individual lecanemab PK parameters (Majid 2024 Table 1 equations)
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl * (ALB / alb_ref)^e_alb_cl *
      e_sexf_cl^SEXF * e_ada_pos_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc * e_sexf_vc^SEXF *
      e_race_japanese_vc^RACE_JAPANESE
    vp <- exp(lvp + etalvp) * e_race_japanese_vp^RACE_JAPANESE
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual plasma p-tau181 PD parameters
    rbase <- exp(lrbase + etalrbase) *
      (WT / wt_ref)^e_wt_rbase *
      (SCORE_MMSE / mmse_ref)^e_mmse_rbase *
      e_apoe4_carrier_rbase^APOE4_CARRIER
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    kout <- exp(lkout + etalkout) / hrs_per_year
    kin  <- rbase * kout   # Bhagunde 2026: Kin is set so the untreated pool sits at baseline

    # Lecanemab disposition
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Serum lecanemab concentration -- the PD driver. This is an individual
    # prediction (IPRED), matching the paper's use of post hoc individual PK
    # parameters; it therefore carries no residual error of its own.
    Cc <- central / vc

    # Plasma p-tau181 indirect response with Emax inhibition of production
    # dR/dt = Kin * (1 - Emax * Conc / (EC50 + Conc)) - R(t) * Kout
    # (Bhagunde 2026 Section 2.2.2)
    d/dt(ptau181) <- kin * (1 - emax * Cc / (ec50 + Cc)) - kout * ptau181
    ptau181(0)    <- rbase

    ptau181 ~ prop(propSd)
  })
}

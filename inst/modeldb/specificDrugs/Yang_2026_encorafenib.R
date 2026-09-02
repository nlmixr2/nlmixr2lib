Yang_2026_encorafenib <- function() {
  description <- "Semi-mechanistic global population PK model for encorafenib in healthy participants and patients with BRAF V600-mutant melanoma, metastatic colorectal cancer, non-small cell lung cancer, and other solid tumors: a two-compartment model with first-order oral absorption whose apparent clearance is multiplied by a unitary-baseline enzyme pool (CL/F * ENZ), the pool being driven by a concentration-dependent sigmoid Emax autoinduction of the enzyme production rate (Emax 1.861, EC50 9.097 ng/mL, Hill exponent fixed at 10, turnover half-life 64.3 h); CL/F therefore rises from 12.2 L/h on day 1 to 35 L/h at steady state (a 186% increase) within about 14 days of once-daily dosing. Age (power) and tumor type (metastatic CRC and a pooled other-tumor group versus a melanoma reference) act on CL/F day 1, and baseline body weight (power) acts on Vc/F."
  reference <- paste(
    "Yang D. Z., Hahn E. M., Piscitelli J., Pithavala Y. K., Hibma J. E. (2026).",
    "Population pharmacokinetic modeling of encorafenib in healthy participants",
    "and patients with BRAF V600-mutant solid tumors: a semi-mechanistic",
    "autoinduction model.",
    "Clinical Pharmacokinetics 65(5):693-707.",
    "doi:10.1007/s40262-025-01608-y.",
    "The enzyme-turnover autoinduction structure follows the same",
    "unitary-baseline enzyme-pool idiom used by",
    "modellib('Smythe_2012_rifampicin') and modellib('Clewe_2015_rifampicin'),",
    "with the addition of a Hill exponent on the Emax induction term.",
    sep = " "
  )
  vignette <- "Yang_2026_encorafenib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Baseline age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power model on CL/F day 1, centered on the analysis-population median of 58 years (Yang 2026 Eq. 6; Table 1 median age 58 years, range 19-94 years). Continuous covariates were entered as power models centered on the population median (Yang 2026 Methods 2.3.2).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power model on Vc/F, centered on the analysis-population median of 76 kg (Yang 2026 Eq. 7; Table 1 median body weight 76 kg, range 34-168 kg). Baseline (time-fixed), not time-varying; 3 of 1310 participants had body weight imputed to the population median (Yang 2026 Table 1 footnote a and Methods 2.3.5).",
      source_name        = "BWT"
    ),
    TUMTP_CRC = list(
      description        = "Metastatic colorectal cancer tumor-type indicator (1 = metastatic CRC, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = melanoma when TUMTP_OTHER is also 0 (melanoma is the reference tumor type in Yang 2026 Eq. 6)",
      notes              = "n = 240 / 1310 (18%) of the pooled cohort (Yang 2026 Table 1), drawn from studies C4221001 (CLGX818X2103) and BEACON (C4221009). Enters CL/F day 1 as the linear categorical multiplier (1 + e_tumtp_crc_cl * TUMTP_CRC) with e_tumtp_crc_cl = -0.175, i.e. mCRC lowers CL/F day 1 by 17.5% relative to the melanoma reference. A subject may have at most one of TUMTP_CRC / TUMTP_OTHER set to 1; both zero means melanoma.",
      source_name        = "mCRC"
    ),
    TUMTP_OTHER = list(
      description        = "Pooled 'other tumor type' indicator (1 = healthy participant, NSCLC, or other solid tumor; 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = melanoma when TUMTP_CRC is also 0 (melanoma is the reference tumor type in Yang 2026 Eq. 6)",
      notes              = "Per-paper composition: healthy participants n = 15 (1%), lung / NSCLC n = 98 (7%), and other solid tumors n = 65 (5%), pooled to n = 178 / 1310 (13%) (Yang 2026 Table 1 and Methods 2.3.5, 'Tumor type was evaluated as melanoma versus CRC versus all other tumor types, which included healthy (1%), lung (7%), and other (5%)'). Enters CL/F day 1 as the linear categorical multiplier (1 + e_tumtp_other_cl * TUMTP_OTHER) with e_tumtp_other_cl = -0.0938. Because the pool mixes healthy volunteers with two solid-tumor groups, this column is not interchangeable with another paper's TUMTP_OTHER.",
      source_name        = "other tumor"
    )
  )

  # Covariates the paper screened in the stepwise covariate model but did not
  # retain in the final model (Yang 2026 Methods 2.3.2, Results 3.2, and
  # Discussion paragraph 4). Documented here so the provenance of the covariate
  # screen survives without declaring unused entries in covariateData.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F, Vc/F, and Ka; not retained (Yang 2026 Results 3.2 and Discussion). Cohort was 594 (45%) female / 716 (55%) male (Table 1)."
    ),
    ALB = list(
      description = "Baseline serum albumin.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Tested on CL/F as a hepatic-function marker; not statistically significant (Yang 2026 Discussion paragraph 6). Cohort median 42 g/L (Table 1)."
    ),
    TPRO = list(
      description = "Baseline total serum protein.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Tested on CL/F; not retained (Yang 2026 Discussion paragraph 4). Cohort median 70 g/L (Table 1)."
    ),
    TBILI = list(
      description = "Baseline total bilirubin.",
      units       = "mg/dL as reported by Yang 2026 Table 1 (canonical register unit is umol/L; multiply by 17.104 to convert)",
      type        = "continuous",
      notes       = "Tested on CL/F as a hepatic-function marker; not statistically significant (Yang 2026 Discussion paragraph 6). Cohort median 0.47 mg/dL (Table 1)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL/F as a hepatic-function marker; not statistically significant (Yang 2026 Discussion paragraph 6). Cohort median 20 U/L (Table 1)."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "NOT tested, because of a correlation greater than 0.6 with AST (Yang 2026 Table 1 footnote d). Cohort median 17 U/L."
    ),
    LDH = list(
      description = "Baseline lactate dehydrogenase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL/F; not retained (Yang 2026 Discussion paragraph 4). Cohort median 207 U/L (Table 1)."
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate (MDRD formula).",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Tested on CL/F; not statistically significant, consistent with encorafenib's minimal renal clearance (Yang 2026 Discussion paragraph 5; Table 1 footnote c). Cohort median 91.62 mL/min/1.73 m^2."
    ),
    WHO_PS = list(
      description = "Baseline ECOG performance status.",
      units       = "(integer score 0-2)",
      type        = "categorical",
      notes       = "Tested on CL/F; not retained (Yang 2026 Results 3.2 and Discussion paragraph 4). Cohort: 0 in 797 (61%), 1 in 491 (37%), 2 in 7 (1%), missing in 15 (1%), with missing assigned to the most common value (Table 1, Methods 2.3.5)."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant moderate or strong CYP3A inhibitor use.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F as absence-or-weak versus moderate-or-strong; not retained. The paper attributes the null result to the small number of exposed participants (moderate n = 149, strong n = 62) rather than to an absent interaction, and notes that a dedicated study showed encorafenib AUCinf increased 180% with a strong CYP3A inhibitor (Yang 2026 Table 1, Methods 2.3.5, Discussion paragraph 4). Concomitant CYP3A INDUCER was not tested at all, because 99% of participants had absent or only weak inducer use."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "encorafenib",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "encorafenib",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "encorafenib",
      units    = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    enz_pool = list(
      analyte  = "encorafenib-metabolizing enzyme pool (predominantly CYP3A4)",
      units    = "fraction of baseline activity",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1310L,
    n_studies      = 9L,
    age_range      = "19-94 years",
    age_median     = "58 years",
    weight_range   = "34-168 kg",
    weight_median  = "76 kg",
    sex_female_pct = 45,
    race_ethnicity = c(White = 90, Black = 1, Asian = 6, Other = 3),
    disease_state  = "BRAF V600-mutant solid tumors plus healthy participants: melanoma n = 892 (68%), metastatic colorectal cancer n = 240 (18%), lung / NSCLC n = 98 (7%), other solid tumors n = 65 (5%), healthy n = 15 (1%) (Yang 2026 Table 1).",
    dose_range     = "Oral encorafenib 50-800 mg once daily (dose-escalation studies spanned 50-700 mg after a single dose), as monotherapy or in combination with binimetinib 45 mg twice daily or with cetuximab with or without mFOLFOX6. Approved regimens contributing most of the data are 450 mg once daily (melanoma / NSCLC with binimetinib) and 300 mg once daily (metastatic CRC with cetuximab).",
    regions        = "Global; nine phase 1, 2, and 3 studies (ARRAY-162-105 part 1, COLUMBUS, POLARIS, C4221010, LOGIC 2, C4221001, BEACON, PHAROS, C4221005)",
    notes          = "1310 participants enrolled; 1299 participants contributed 8651 evaluable observations after excluding 263 values below the limit of quantification or missing (from 11 participants) via the M1 method, and 3 further observations with CWRES > 6 were excluded from the final model (Yang 2026 Results 3.1 and 3.2). Combination therapy: none 372 (28%), binimetinib 689 (53%), other 249 (19%). Race was not evaluated as a covariate because 90% of participants were White."
  )

  ini({
    # =========================================================================
    # Structural PK. Two-compartment disposition with first-order oral
    # absorption, all parameters apparent (F-relative). Time unit is HOURS and
    # concentrations are ng/mL, matching Yang 2026 Table 2 (EC50 is reported in
    # ng/mL). Amounts are in mg, so the observation below carries the explicit
    # 1000 mg/L -> ng/mL conversion.
    # =========================================================================
    lka  <- log(0.954)   ; label("First-order absorption rate constant Ka (1/h)")                          # Yang 2026 Table 2: theta_Ka = 0.954 /h (RSE 3.3%; 95% CI 0.892-1.017; bootstrap 0.970, 0.877-1.083)
    lcl  <- log(12.238)  ; label("Apparent oral clearance CL/F on day 1, before autoinduction (L/h)")       # Yang 2026 Table 2: theta_CL/F day 1 = 12.238 L/h (RSE 3.0%; 95% CI 11.515-12.961; bootstrap 11.951, 11.066-12.653)
    lvc  <- log(61.73)   ; label("Apparent central volume of distribution Vc/F (L)")                        # Yang 2026 Table 2: theta_Vc/F = 61.73 L (RSE 2.6%; 95% CI 58.606-64.855; bootstrap 62.48, 58.78-66.99)
    lq   <- log(1.045)   ; label("Apparent intercompartmental clearance Q/F (L/h)")                         # Yang 2026 Table 2: theta_Q/F = 1.045 L/h (RSE 2.4%; 95% CI 0.996-1.094; bootstrap 1.095, 1.023-1.165)
    lvp  <- log(54.549)  ; label("Apparent peripheral volume of distribution Vp/F (L)")                     # Yang 2026 Table 2: theta_Vp/F = 54.549 L (RSE 2.4%; 95% CI 52.008-57.090; bootstrap 52.677, 49.125-90.174)

    # =========================================================================
    # Semi-mechanistic enzyme-turnover autoinduction (Yang 2026 Fig. 1 and
    # Eqs. 1-2):
    #   d(A_ENZ)/dT = kin * (1 + EFF) - kout * A_ENZ,  with kin = kout = kENZ
    #                 so that A_ENZ = 1 in the drug-free steady state,
    #   EFF         = Emax * Cp^gamma / (EC50^gamma + Cp^gamma),
    # and the apparent clearance consumed by the central ODE is CL/F * ENZ
    # (Yang 2026 Fig. 1 flux label).
    #
    # kENZ is entered here rather than the turnover half-life that Table 2
    # reports, following the enzyme-pool idiom already used by
    # Smythe_2012_rifampicin and Clewe_2015_rifampicin:
    #   kENZ = ln(2) / 64.279 h = 0.010784 /h,
    # which reproduces the K_ENZ = 0.0108 /h printed in Yang 2026 Results 3.3.
    # Table 2 places the IIV on the half-life; on the log scale an eta on kENZ
    # is exactly the negative of an eta on the half-life with the same
    # variance, so etalkenz below is distributionally identical to the
    # published IIV on theta_HL.
    # =========================================================================
    lkenz <- log(log(2) / 64.279) ; label("First-order turnover rate constant of the enzyme pool kENZ (1/h)")  # Yang 2026 Table 2: theta_Turnover HL = 64.279 h (RSE 5.9%; 95% CI 56.869-71.689; bootstrap 55.569, 52.816-60.591); Results 3.3 prints the equivalent K_ENZ = 0.0108 /h
    lemax <- log(1.861)           ; label("Maximal fractional increase in the enzyme production rate Emax (unitless)")  # Yang 2026 Table 2: theta_Emax = 1.861 (RSE 2.5%; 95% CI 1.768-1.953; bootstrap 1.931, 1.840-2.398); Results 3.3 reports the corresponding 186% rise in CL/F
    lec50 <- log(9.097)           ; label("Encorafenib plasma concentration giving half of Emax, EC50 (ng/mL)")         # Yang 2026 Table 2: theta_EC50 = 9.097 ng/mL (RSE 4.0%; 95% CI 8.381-9.813; bootstrap 8.386, 8.033-9.833); Results 3.3 rounds it to 9 ng/mL
    lhill <- fixed(log(10))       ; label("Hill (sigmoidicity) exponent of the autoinduction Emax function (unitless)")  # Yang 2026 Table 2: theta_Gamma = 10, reported with no RSE, no 95% CI and no bootstrap interval

    lfdepot <- fixed(log(1))      ; label("Oral bioavailability F (1 because CL/F, Vc/F, Q/F and Vp/F are apparent)")     # Structural anchor: Yang 2026 Fig. 1 doses the absorption compartment with F * Dose and Table 2 reports every disposition parameter as an apparent (/F) quantity

    # =========================================================================
    # Covariate effects (Yang 2026 Eqs. 6-7, coefficients in Table 2):
    #   TV(CL/F) = 12.2 * (AGE/58)^-0.326 * (1 - 0.175*mCRC) * (1 - 0.0938*other)
    #   TV(Vc/F) = 61.7 * (BWT/76)^0.588
    # The categorical terms are written here as (1 + e_* * IND) so that the
    # stored coefficient carries the published sign from Table 2 directly.
    # =========================================================================
    e_age_cl          <- -0.326  ; label("Power exponent on (AGE / 58 years) for CL/F day 1 (unitless)")          # Yang 2026 Table 2: theta_Age on CL/F day 1 = -0.326 (RSE 20.7%; 95% CI -0.459 to -0.194); Eq. 6 prints the same exponent with the 58-year centering value
    e_wt_vc           <-  0.588  ; label("Power exponent on (WT / 76 kg) for Vc/F (unitless)")                    # Yang 2026 Table 2: theta_BWT on Vc/F = 0.588 (RSE 14.5%; 95% CI 0.421-0.755); Eq. 7 prints the same exponent with the 76-kg centering value
    e_tumtp_crc_cl    <- -0.175  ; label("Fractional change in CL/F day 1 for metastatic colorectal cancer versus the melanoma reference")   # Yang 2026 Table 2: theta_mCRC tumor on CL/F day 1 = -0.175 (RSE 26.9%; 95% CI -0.267 to -0.083); Eq. 6 prints the matching (1 - 0.175 * mCRC) factor
    e_tumtp_other_cl  <- -0.0938 ; label("Fractional change in CL/F day 1 for the pooled other-tumor group versus the melanoma reference")   # Yang 2026 Eq. 6 prints (1 - 0.0938 * other tumor); Table 2 rounds the same theta to -0.094 (RSE 55.2%; 95% CI -0.195 to 0.008)

    # =========================================================================
    # Inter-individual variability. Yang 2026 Eq. 3 uses exponential random
    # effects and Eq. 4 defines the reported %CV as sqrt(omega^2) * 100, so the
    # Table 2 CV(%) column IS omega * 100 and omega^2 = (CV/100)^2. The
    # 15.81% entries reproduce exactly the 0.025 variance that Methods 2.2
    # states was assigned to random effects estimated near zero under
    # SAEM-IMP with mu-referencing; those, and the Ka IIV, are fixed.
    # Only the CL/F and Vc/F variances were estimated. The final omega matrix
    # is diagonal (Methods 2.3.1: all post-hoc eta correlations had R^2 < 0.6).
    # =========================================================================
    etalcl   ~ 0.319                # Yang 2026 Table 2: IIV on theta_CL/F = 56.48% CV (RSE 4.2%; shrinkage 13.0%; 95% CI 54.13-58.74) -> omega^2 = 0.5648^2 = 0.319
    etalvc   ~ 0.25                 # Yang 2026 Table 2: IIV on theta_Vc/F = 50% CV (RSE 5.7%; shrinkage 30.4%; 95% CI 47.12-52.73) -> omega^2 = 0.50^2 = 0.25
    etalka   ~ fixed(0.25)          # Yang 2026 Table 2: IIV on theta_Ka = 50% CV with no RSE or CI; Results 3.2 sets it to 50% because the absorption phase could not be characterized
    etalq    ~ fixed(0.025)         # Yang 2026 Table 2: IIV on theta_Q/F = 15.81% CV with no RSE or CI; Results 3.2 and Methods 2.2 give the variance 0.025
    etalvp   ~ fixed(0.025)         # Yang 2026 Table 2: IIV on theta_Vp/F = 15.81% CV with no RSE or CI; Results 3.2 and Methods 2.2 give the variance 0.025
    etalkenz ~ fixed(0.025)         # Yang 2026 Table 2: IIV on theta_HL = 15.81% CV with no RSE or CI; sign-flipped onto kENZ (see the kENZ note above), same variance 0.025
    etalemax ~ fixed(0.025)         # Yang 2026 Table 2: IIV on theta_Emax = 15.81% CV with no RSE or CI; Results 3.2 and Methods 2.2 give the variance 0.025
    etalec50 ~ fixed(0.025)         # Yang 2026 Table 2: IIV on theta_EC50 = 15.81% CV with no RSE or CI; Results 3.2 and Methods 2.2 give the variance 0.025
    etalhill ~ fixed(0.025)         # Yang 2026 Table 2: IIV on theta_Gamma = 15.81% CV with no RSE or CI; retained even though the typical Hill exponent itself was not estimated

    # =========================================================================
    # Residual error. Yang 2026 Eq. 5, as typeset on p. 696, is
    #   Y_ij = F_ij * (1 + W * eps_ij) + W * eps_ij
    # with sigma^2 fixed to 1. The SAME unsubscripted W appears in both terms
    # (verified against the typeset PDF equation, not only the text extract),
    # so the residual standard deviation is
    #   SD = W * (1 + F_ij) = W + W * F_ij,
    # i.e. a LINEAR sum of an additive and a proportional term sharing a single
    # epsilon. That is exactly nlmixr2's combined1() error model with
    # add == prop == W, and it is NOT the same as add() + prop() alone, which
    # would combine two independent epsilons in quadrature
    # (SD = W * sqrt(1 + F^2)).
    #
    # Table 2 publishes a single error theta (theta_Prop err = 0.589), which is
    # W. Both components are therefore set to that one published value; no
    # value is invented. The Methods prose calls this "a proportional model"
    # and says W "was estimated as two of the thetas", which the printed
    # equation does not resolve into two distinct symbols; the equation is
    # followed here (see the vignette Errata).
    # =========================================================================
    addSd  <- 0.589 ; label("Additive component of the combined1 residual error (ng/mL)")   # Yang 2026 Eq. 5: the additive term carries the same W as the proportional term; Table 2 theta_Prop err = 0.589 (RSE 1.0%; 95% CI 0.577-0.601; bootstrap 0.590, 0.577-0.603)
    propSd <- 0.589 ; label("Proportional component of the combined1 residual error (fraction)")  # Yang 2026 Table 2: theta_Prop err = 0.589 (RSE 1.0%; 95% CI 0.577-0.601; bootstrap 0.590, 0.577-0.603); Eq. 5 W
  })

  model({
    # --- 1. Individual disposition parameters. -------------------------------
    #        cl_base is the day-1 (pre-induction) apparent clearance carrying
    #        the age power term and the two tumor-type multipliers of Eq. 6;
    #        vc carries the body-weight power term of Eq. 7. Both centering
    #        values (58 years, 76 kg) are the analysis-population medians
    #        reported in Yang 2026 Table 1.
    ka      <- exp(lka + etalka)
    cl_base <- exp(lcl + etalcl) * (AGE / 58) ^ e_age_cl *
      (1 + e_tumtp_crc_cl * TUMTP_CRC) *
      (1 + e_tumtp_other_cl * TUMTP_OTHER)
    vc      <- exp(lvc + etalvc) * (WT / 76) ^ e_wt_vc
    q       <- exp(lq + etalq)
    vp      <- exp(lvp + etalvp)

    # --- 2. Autoinduction parameters. ----------------------------------------
    kenz <- exp(lkenz + etalkenz)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)

    # --- 3. Plasma concentration and the induction driver. -------------------
    #        Amounts are mg and vc is L, so central/vc is mg/L; the factor
    #        1000 converts to the ng/mL scale on which EC50 is reported.
    #        eff is Yang 2026 Eq. 2 with Cp taken as the central-compartment
    #        concentration (Fig. 1: "the EFF of encorafenib concentration in
    #        the central compartment increased the K_ENZ").
    Cc  <- 1000 * central / vc
    eff <- emax * Cc ^ hill / (ec50 ^ hill + Cc ^ hill)

    # --- 4. Time-varying apparent clearance (Yang 2026 Fig. 1: CL/F * ENZ). --
    #        enz_pool is unity at baseline, so cl equals cl_base before the
    #        first dose and approaches cl_base * (1 + Emax) once the pool has
    #        equilibrated at full induction.
    cl <- cl_base * enz_pool

    # --- 5. ODE system. ------------------------------------------------------
    #        Elimination and distribution are written as clearance times
    #        concentration rather than as micro-constants so that the
    #        autoinduction multiplier is read directly off cl; the two forms
    #        are algebraically identical (cl * central / vc == kel * central).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl * central / vc -
      q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <-  q * central / vc - q * peripheral1 / vp

    # --- 6. Enzyme pool (Yang 2026 Eq. 1 with kin = kout = kENZ). ------------
    d/dt(enz_pool)    <-  kenz * (1 + eff) - kenz * enz_pool

    # Normalized to unity at baseline (Yang 2026 Methods 2.3: "To normalize
    # the enzyme concentrations to unity at baseline, kin was set equal to
    # kout").
    enz_pool(0) <- 1.0

    # --- 7. Bioavailability anchor. ------------------------------------------
    f(depot) <- exp(lfdepot)

    # --- 8. Observation. -----------------------------------------------------
    #        combined1() gives SD = addSd + propSd * Cc, which with
    #        addSd == propSd == W reproduces Yang 2026 Eq. 5 exactly.
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}

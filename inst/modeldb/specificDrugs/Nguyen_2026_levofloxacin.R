Nguyen_2026_levofloxacin <- function() {
  description <- paste(
    "One-compartment oral population PK model for levofloxacin in Vietnamese adults treated for",
    "multidrug-resistant tuberculosis (Nguyen 2026), fitted jointly to paired plasma and saliva",
    "concentrations. Absorption is first order with a lag time, and both the absorption rate",
    "constant (4.18 1/h) and the lag time (0.95 h) were held constant. Saliva is carried as a",
    "kinetically distinct hypothetical effect compartment driven by the central compartment",
    "through a secretion rate constant (kin_saliva = 4.929 1/h) with irreversible salivary loss",
    "(kel_saliva = 5.084 1/h) and no reabsorption leg; the saliva state shares the central volume,",
    "so the steady-state saliva:plasma exposure ratio is the parameter-free constant",
    "kin_saliva/kel_saliva = 0.9695, matching the 0.928 scale factor the authors' competing",
    "scale-factor saliva model estimated. The authors selected this distinct-compartment structure",
    "over that scale-factor model, the same choice made for linezolid in Nguyen 2026 and the",
    "opposite of the choice made for busulfan in Xu 2023. No covariate was retained: neither total",
    "body weight nor fat-free mass improved the fit by allometric scaling, and age, sex, renal and",
    "hepatic function markers were all screened and rejected. Apparent volume of distribution",
    "(278.88 L) is about three times the commonly reported value, which the authors attribute to",
    "sparse sampling over 0-5 h post-dose. Interindividual variability is carried on apparent",
    "clearance and apparent volume. Combined additive plus proportional residual errors apply",
    "separately to plasma and saliva. The model underpins saliva-only limited sampling strategies",
    "for predicting plasma AUC(0-24).",
    sep = " "
  )
  reference <- paste(
    "Nguyen TA, Nguyen TP, Nguyen AT, Dinh LV, Nguyen HB, Vu HD, Nguyen TNB, Vu D, Fox GJ,",
    "Alffenaar JWC, Stocker SL. Single Saliva Sample Model-Informed Precision Dosing of",
    "Levofloxacin for Multidrug-Resistant Tuberculosis.",
    "Clin Pharmacokinet. 2026. doi:10.1007/s40262-026-01619-3",
    sep = " "
  )
  vignette <- "Nguyen_2026_levofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "levofloxacin", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "levofloxacin", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    # The saliva state is a hypothetical effect compartment (Nguyen 2026
    # Methods 2.3: 'a saliva bio-compartment (i.e., a hypothetical effect
    # compartment, which does not account for mass balance)'), so it holds a
    # notional amount that is rescaled by the central volume to give the
    # observed salivary concentration. It does not deplete the central
    # compartment -- see the model() note.
    saliva = list(
      analyte = "levofloxacin", units = "mg",
      specimen = "saliva", verified = TRUE
    )
  )

  # No covariate was retained in the final model, so covariateData is empty and
  # model() references no covariate column. Everything the authors screened is
  # recorded in covariatesDataExcluded below.
  covariateData <- list()

  # Covariates screened and NOT retained. Body size was tested explicitly as an
  # allometric descriptor (Methods 2.3, Results 3.2): total body weight with
  # fixed exponents gave dOFV = -1.93 and with estimated exponents dOFV = -0.64,
  # and fat-free mass made the fit worse (dOFV = +5.42 fixed, +4.9 estimated).
  # The remaining demographic, renal and hepatic markers were screened by
  # stepwise covariate modelling in Perl-speaks-NONMEM and none correlated with
  # clearance or volume (Results 3.2, Supplementary Figure S8).
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight at baseline.",
      units = "kg", type = "continuous",
      notes = paste(
        "Median 50 kg (IQR 45-56), Nguyen 2026 Table 1; a column in the Table S4 $INPUT list.",
        "Tested as an allometric size descriptor on CL/F and V/F with both fixed exponents",
        "(0.75 on clearance, 1 on volume; dOFV = -1.93 versus the base model) and estimated",
        "exponents (dOFV = -0.64). Neither reached the pre-specified 3.84-point forward-inclusion",
        "criterion, so no weight term is carried (Results 3.2)."
      )
    ),
    FFM = list(
      description = "Fat-free mass.",
      units = "kg", type = "continuous",
      notes = paste(
        "Median 41.4 kg (IQR 36.7-45.6), Nguyen 2026 Table 1. Derived from total body weight and",
        "height by the Janmahasatian / Anderson-Holford relation printed as Supplementary",
        "Information S2 Eq. 7, FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT), with WHSmax and",
        "WHS50 of 42.92 and 30.93 kg/m^2 for men and 37.99 and 35.98 kg/m^2 for women. Tested as",
        "an allometric size descriptor and made the fit WORSE (dOFV = +5.42 with fixed exponents,",
        "+4.9 with estimated exponents), so it is not carried (Results 3.2)."
      )
    ),
    AGE = list(
      description = "Age.", units = "years", type = "continuous",
      notes = paste(
        "Median 44 years (IQR 33-51), Nguyen 2026 Table 1. Screened against clearance and volume",
        "by stepwise covariate modelling and not retained (Methods 2.3, Results 3.2,",
        "Supplementary Figure S8)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)", type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "17 of 57 patients female (29.8%), Nguyen 2026 Table 1. Listed among the screened",
        "covariates (Methods 2.3, 'Covariates included age, sex, renal and hepatic function",
        "markers') and not retained (Results 3.2). Note that sex also enters the FFM relation",
        "above through its sex-specific WHSmax / WHS50 constants."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (renal-function marker).",
      units = "mL/min (raw, NOT BSA-normalized)", type = "continuous",
      notes = paste(
        "Median 71 mL/min (IQR 65.1-84.7), Nguyen 2026 Table 1. The estimating equation is not",
        "stated in the paper. All patients had renal function within normal limits (Results 3.1),",
        "so the cohort carries little information about renal impairment -- which matters because",
        "levofloxacin is predominantly renally cleared. Dropped in the Table S4 $INPUT list",
        "(CLCR=DROP) for the printed run. Screened and not retained (Results 3.2,",
        "Supplementary Figure S8)."
      )
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units = "umol/L", type = "continuous",
      notes = paste(
        "Median 77 umol/L (IQR 71-85), Nguyen 2026 Table 1, printed as 'SCR'. Carried into the",
        "Table S4 $INPUT list as CRE. Screened and not retained (Results 3.2,",
        "Supplementary Figure S8)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase (hepatic-function marker).",
      units = "U/L", type = "continuous",
      notes = paste(
        "Median 20 U/L (IQR 15-36), Nguyen 2026 Table 1. Within normal limits in all patients",
        "(Results 3.1). Screened and not retained (Results 3.2, Supplementary Figure S8)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase (hepatic-function marker).",
      units = "U/L", type = "continuous",
      notes = paste(
        "Median 24 U/L (IQR 19-32), Nguyen 2026 Table 1. Within normal limits in all patients",
        "(Results 3.1). Screened and not retained (Results 3.2, Supplementary Figure S8)."
      )
    )
    # Height (Table 1 median 162 cm), BMI (19.1 kg/m^2) and platelet count (a
    # dropped column in the Table S4 $INPUT list) are NOT listed here. Height
    # is reported as a baseline characteristic and enters the model only
    # indirectly, as an input to the FFM derivation above; BMI and platelets
    # are baseline characteristics that the paper never describes as screened
    # covariates. They are recorded in population instead of being claimed as
    # screened-and-rejected effects.
  )

  population <- list(
    species          = "human",
    n_subjects       = 57,
    n_studies        = 1,
    age_median       = "44 years (IQR 33-51)",
    age_range        = "adults aged 18 years and over; the full range is not reported",
    weight_median    = "50 kg (IQR 45-56)",
    weight_range     = "IQR 45-56 kg; the full range is not reported",
    height_median    = "162 cm (IQR 160-167)",
    ffm_median       = "41.4 kg (IQR 36.7-45.6)",
    bmi_median       = "19.1 kg/m^2 (IQR 17.5-20.8)",
    sex_female_pct   = 29.8,
    race_ethnicity   = c(Asian = 100),
    disease_state    = "multidrug-resistant pulmonary tuberculosis (MDR-TB)",
    renal_function   = "within normal limits in all patients; creatinine clearance median 71 mL/min (IQR 65.1-84.7), serum creatinine median 77 umol/L (IQR 71-85)",
    hepatic_function = "within normal limits in all patients; ALT median 20 U/L (IQR 15-36), AST median 24 U/L (IQR 19-32)",
    dose_range       = "oral levofloxacin at steady state, 750-1000 mg once daily (15-20 mg/kg/day, median 17.9 mg/kg/day)",
    regions          = "Vietnam (four provinces)",
    notes            = paste(
      "Pharmacokinetic sub-study of the V-SMART trial (ACTRN12620000681954), prospective and",
      "observational. Sixty patients met the inclusion criteria and 57 had evaluable",
      "pharmacokinetic data, contributing 342 paired plasma-saliva samples drawn at pre-dose, 2 h",
      "and 5 h post-dose after at least 7 days of treatment (i.e. at steady state). Three patients",
      "who had completed their levofloxacin course had below-limit-of-quantification samples and",
      "were excluded. Baseline demographics are Nguyen 2026 Table 1. The assay lower limit of",
      "quantification was 0.5 mg/L in both matrices; one plasma value (0.49 mg/L) and its paired",
      "saliva value (0.37 mg/L) fell below it and were retained 'as measured' rather than censored.",
      "The dose received by each individual is not published, only the 750-1000 mg range and the",
      "17.9 mg/kg/day median. All patients had renal and hepatic function markers within normal",
      "limits, so the model carries no information about organ impairment. The cohort is lean by",
      "international standards (median BMI 19.1 kg/m^2) and demographically uniform, which the",
      "authors flag as a limit on generalisability."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma disposition. Nguyen 2026 Table 2 'Final model' column
    # (objective function value = 695.794). Parameters are APPARENT
    # (CL/F, V/F) because absolute bioavailability was not identifiable
    # from the oral-only design; the Table S4 control stream carries no
    # bioavailability THETA at all, so F is structurally 1.
    # ------------------------------------------------------------------
    lcl <- log(10.311); label("Apparent clearance, CL/F (L/h)")                      # Table 2 final model theta_CL = 10.311 L/h (RSE 0.01%; SIR median 10.310, 95% CI 10.308-10.312); Table S4 $THETA (8, 10.3, 12)
    lvc <- log(278.88); label("Apparent central volume of distribution, V/F (L)")    # Table 2 final model theta_V = 278.88 L (RSE 0.008%; SIR median 278.88, 95% CI 278.82-278.92); Table S4 $THETA (250, 278, 300)

    # Absorption rate constant, held constant at a published value. Few
    # concentrations were observed in the absorption phase (sampling was
    # pre-dose, 2 h and 5 h only), so estimating Ka destabilised the fit --
    # in the no-lag base model it came out at 2.77 1/h with %RSE 105%
    # (Table S1). The authors adopted 4.18 1/h from published levofloxacin
    # analyses to stabilise convergence and gradients.
    lka <- fixed(log(4.18)); label("First-order absorption rate constant, Ka (1/h)") # Table 2 final model Ka = 4.18 (Fixed); Table S4 $THETA (4.18) FIX; Results 3.2 'fixed at 4.18 h-1 based on published estimates'

    # Absorption lag time, held constant. Estimated at 0.947 h in the base
    # plasma model (Table S1, RSE 3%) but %RSE rose to 141% when estimated
    # in the combined plasma-plus-saliva model, so it was held at the
    # rounded 0.95 h to maintain stability (Results 3.2).
    ltlag <- fixed(log(0.95)); label("Absorption lag time, Tlag (h)")                # Table 2 final model Tlag = 0.95 (Fixed); Table S4 $THETA (0.95) FIX; base-model estimate 0.947 h in Table S1

    # ------------------------------------------------------------------
    # Saliva effect-compartment rate constants. Nguyen 2026 Table 2 final
    # model; Figure 1 names them Kabs and Kel. Canonical names use the
    # registered kin_<compartment> tissue-exchange family plus
    # kel_<compartment> for the irreversible salivary loss -- the paper's
    # own symbol for that term is 'Kel', which collides with the canonical
    # central elimination rate constant, so the suffix is what keeps them
    # apart. Unlike the linezolid model from the same group
    # (modellib('Nguyen_2026_linezolid')), this model has NO saliva-to-
    # central reabsorption leg: Figure 1 draws a single dotted arrow into
    # the saliva compartment and a single dotted arrow out of it, and the
    # Table S4 $PK block defines only K23 and K30 (no K32).
    # ------------------------------------------------------------------
    lkin_saliva <- log(4.929); label("Central-to-saliva secretion rate constant, Kabs (1/h)")  # Table 2 final model K_abs = 4.929 (RSE 0.01%; SIR median 4.929, 95% CI 4.927-4.930); Table S4 K23, $THETA (1, 4.94, 6)
    lkel_saliva <- log(5.084); label("Irreversible salivary elimination rate constant, Kel (1/h)") # Table 2 final model K_el = 5.084 (RSE 0.02%; SIR median 5.084, 95% CI 5.082-5.086); Table S4 K30, $THETA (1, 5.08, 6)

    # ------------------------------------------------------------------
    # Interindividual variability. Exponential, P_i = P_TV * exp(eta_i)
    # (Supplementary Information S1.1 Eq. 1). Diagonal OMEGA -- Table S4
    # declares two separate $OMEGA elements with no BLOCK, so CL and V are
    # uncorrelated.
    #
    # VARIANCE CONVENTION. Table 2 labels the IIV rows 'omega (CV%)', which
    # is ambiguous: it can mean sqrt(variance)*100 or the log-normal
    # sqrt(exp(variance)-1)*100. Table S4 settles it, because its $THETA
    # initials are the final estimates rounded (10.3, 278, 4.94, 5.08
    # against 10.311, 278.88, 4.929, 5.084), so its $OMEGA and $SIGMA
    # initials are the final variances too. sqrt(0.166) = 40.74% against a
    # reported 41.2%, and sqrt(0.48) = 69.28% against a reported 69.7% --
    # the log-normal reading would give 42.50% and 78.49%, and 78.49% is
    # nowhere near 69.7%. The $SIGMA rows clinch it independently because
    # two of them are ADDITIVE errors in mg/L, where no CV convention
    # applies at all: sqrt(0.0195) = 0.1396 against a reported 0.1394 mg/L,
    # and sqrt(0.00364) = 0.0603 against a reported 0.06 mg/L. So the table
    # prints standard deviations and the control stream holds variances.
    # ------------------------------------------------------------------
    etalcl ~ 0.169744  # Table 2 final model omega_CL = 41.2 CV% -> 0.412^2 (RSE 62.9%, shrinkage 2%); Table S4 $OMEGA 0.166
    etalvc ~ 0.485809  # Table 2 final model omega_V = 69.7 CV% -> 0.697^2 (RSE 11.1%, shrinkage 12%); Table S4 $OMEGA 0.48

    # IIV on Ka, Tlag, Kabs and Kel is NOT carried. Table S4 fixes all four
    # to zero ($OMEGA '0 FIX' for IIV-KA, IIV-ALAG, IIV-K23, IIV-K30) and
    # Results 3.2 explains why for the saliva pair: 'The interindividual
    # variability on first-order saliva absorption rate (Kabs) and
    # elimination rate from the saliva compartment (Kel) could not be
    # estimated because including them rendered the model unstable'. They
    # are omitted rather than written as `~ fixed(0)` because a zero-variance
    # diagonal makes OMEGA singular and breaks the Cholesky sampler used by
    # rxSolve.

    # ------------------------------------------------------------------
    # Residual unexplained variability. Combined additive and proportional
    # in both matrices, modelled separately (Results 3.2; Supplementary
    # Information S1.2 Eq. 4, Y = IPRED*(1+EPS(1)) + EPS(2)). Values are the
    # Table 2 entries read as standard deviations per the convention proof
    # above.
    # ------------------------------------------------------------------
    propSd <- 0.224; label("Proportional residual SD for plasma Cc (fraction)")            # Table 2 final model sigma_Plasma_Prop = 22.4 CV% (RSE 0.5%, shrinkage 19%); Table S4 $SIGMA 0.0501 -> sqrt = 0.2238
    addSd <- 0.1394; label("Additive residual SD for plasma Cc (mg/L)")                    # Table 2 final model sigma_Plasma_Add = 0.1394 mg/L (RSE 0.02%, shrinkage 19%); Table S4 $SIGMA 0.0195 -> sqrt = 0.1396
    propSd_Csaliva <- 0.3453; label("Proportional residual SD for saliva Csaliva (fraction)") # Table 2 final model sigma_Saliva_Prop = 34.53 CV% (RSE 0.2%, shrinkage 8%); Table S4 $SIGMA 0.1240 -> sqrt = 0.3521
    addSd_Csaliva <- 0.06; label("Additive residual SD for saliva Csaliva (mg/L)")            # Table 2 final model sigma_Saliva_Add = 0.06 mg/L (RSE 0.03%, shrinkage 8%); Table S4 $SIGMA 0.00364 -> sqrt = 0.0603
  })

  model({
    # 1. Individual parameters. No covariate enters any of them -- the final
    #    model retained none (Results 3.2).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    kin_saliva <- exp(lkin_saliva)
    kel_saliva <- exp(lkel_saliva)

    # 2. Micro-constant. Table S4 $PK writes this as K20 = CL / V.
    kel <- cl / vc

    # 3. ODE system. NONMEM ADVAN13 TRANS1 with three compartments,
    #    COMP = (ABS) / (CENTRAL) / (SALIVA), Table S4.
    #
    #    The saliva state is a DRIVEN (non-depleting) hypothetical effect
    #    compartment. Methods 2.3 says so in as many words -- the plasma
    #    model was 'extended to include a saliva bio-compartment (i.e., a
    #    hypothetical effect compartment, which does not account for mass
    #    balance)' -- and Figure 1's caption repeats it. The central
    #    equation therefore carries no -kin_saliva*central loss term.
    #
    #    This is load-bearing, because Table S4 is headed 'Example of NONMEM
    #    model code' and omits the $DES block that its own ADVAN13 subroutine
    #    requires, while the K23 / K30 naming (borrowed from the general
    #    LINEAR ADVAN5/ADVAN7 subroutines, where those names do imply a
    #    mass-balance matrix) suggests the opposite reading. That reading is
    #    falsified by the paper's own numbers: at steady state the saliva
    #    state would sit at kin_saliva/kel_saliva = 0.9695 times the central
    #    amount, so salivary loss would be kel_saliva * 0.9695 = 4.929 per
    #    hour against a plasma kel of 10.311/278.88 = 0.036973 per hour --
    #    a 134-fold inflation of total clearance, collapsing the steady-state
    #    750 mg AUC(0-24) from 72.7 to about 0.5 mg*h/L and putting every
    #    simulated concentration two orders of magnitude below the 0.5-20 mg/L
    #    range plotted in Figure 3. See the vignette for the full arithmetic.
    #
    #    Same structure as this group's companion linezolid model
    #    (modellib('Nguyen_2026_linezolid')), minus the reabsorption leg.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(saliva) <- kin_saliva * central - kel_saliva * saliva

    # 4. Absorption lag. Table S4 $PK writes ALAG1 = TVALAG*EXP(ETA(4)) with
    #    the eta fixed to zero; compartment 1 is COMP = (ABS), the depot.
    alag(depot) <- tlag

    # 5. Observations. Both matrices are scaled by the central volume:
    #    Table S4 prints S2 = V and CP = A(2)/V for plasma, and CS = A(3)
    #    for saliva. The literal CS = A(3) (an implied 1 L saliva scale) is
    #    falsified by the paper: it puts the steady-state saliva:plasma
    #    exposure ratio at vc * kin_saliva/kel_saliva = 278.88 * 0.9695 =
    #    270, i.e. saliva running ~270x plasma, whereas Figure 3 plots both
    #    matrices on one shared 0-22 mg/L axis with near-identical medians
    #    (about 5.3 mg/L at 2 h post-dose in both panels). Dividing by vc
    #    gives the parameter-free steady-state ratio
    #    kin_saliva/kel_saliva = 4.929/5.084 = 0.9695, against the 0.928
    #    scale factor the authors' competing scale-factor saliva model
    #    estimated on the same data (Table S2) -- agreement to 4%. This is
    #    the same reading operator-ratified for the companion linezolid
    #    model (sidecar request-001 q2, 2026-09-02); see the vignette Errata.
    Cc <- central / vc
    Csaliva <- saliva / vc

    Cc ~ add(addSd) + prop(propSd)
    Csaliva ~ add(addSd_Csaliva) + prop(propSd_Csaliva)
  })
}

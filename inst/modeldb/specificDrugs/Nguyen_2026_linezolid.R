Nguyen_2026_linezolid <- function() {
  description <- paste(
    "One-compartment oral population PK model for linezolid in Vietnamese adults treated for",
    "multidrug-resistant tuberculosis (Nguyen 2026), fitted jointly to paired plasma and saliva",
    "concentrations. Saliva is carried as a kinetically distinct hypothetical effect compartment",
    "driven by the central compartment through a secretion rate constant (kin_saliva = 4.93 1/h),",
    "with reabsorption back towards plasma (kout_saliva = 1.84 1/h) and irreversible salivary loss",
    "(kel_saliva = 2.13 1/h); the saliva state shares the central volume, so the steady-state",
    "saliva:plasma exposure ratio is the parameter-free constant kin_saliva/(kout_saliva +",
    "kel_saliva) = 1.24. The authors selected this distinct-compartment structure over a",
    "scale-factor saliva model (d-2LL = -43.084, 2 df, p < 0.001), the opposite of the choice made",
    "for busulfan in Xu 2023. Apparent central volume increases with total body weight (power",
    "exponent 1.1 referenced to the population median 50 kg); absorption rate and bioavailability",
    "were both fixed. Interindividual variability was estimable only on apparent clearance (44.72%",
    "CV). Separate proportional residual errors apply to plasma (25.86%) and saliva (35.91%). The",
    "model underpins saliva-only limited sampling strategies for predicting plasma AUC(0-24).",
    sep = " "
  )
  reference <- paste(
    "Nguyen TA, Nguyen TP, Nguyen AT, Dinh LV, Nguyen HB, Vu HD, Nguyen TNB, Vu D, Fox GJ,",
    "Alffenaar JWC, Stocker SL. Single Saliva Sample Linezolid Dosing for Multidrug-Resistant",
    "Tuberculosis: A Population Pharmacokinetic Modelling of Plasma and Saliva.",
    "Clin Pharmacokinet. 2026. doi:10.1007/s40262-026-01626-4",
    sep = " "
  )
  vignette <- "Nguyen_2026_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "linezolid", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "linezolid", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    # The saliva state is a hypothetical effect compartment (Nguyen 2026
    # Figure 1 caption: "The saliva bio-compartment is a hypothetical effect
    # compartment"), so it holds a notional amount that is rescaled by the
    # central volume to give the observed salivary concentration. It does not
    # deplete the central compartment -- see the model() note.
    saliva = list(
      analyte = "linezolid", units = "mg",
      specimen = "saliva", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-scaled on the apparent central volume only, referenced to the population median",
        "50 kg (Nguyen 2026 Table 1; IQR 45-56 kg). The exponent 1.1 is estimated, not fixed at an",
        "allometric 1, and its 95% CI (0.62-1.46) spans 1; the paper calls it 'the modest",
        "association in our dataset' and notes it is consistent with a published linezolid popPK",
        "analysis in patients with TB (Zhou 2022). Weight was NOT retained on apparent clearance.",
        "The supplementary control stream (Table S7) writes this as VWT = ((WT/50)**THETA(8)).",
        "Note this cohort is lean by international standards (median BMI 19.1 kg/m^2), so the term",
        "was fitted over a narrow weight range and extrapolates poorly to heavier populations."
      ),
      source_name        = "WT"
    )
  )

  # Covariates screened in the stepwise covariate analysis but NOT retained in
  # the final model (Nguyen 2026 Results 3.2 and Supplementary Figure S8:
  # "other covariates, including age, renal function, and liver enzymes, were
  # not significant"). Documented here for provenance; they are deliberately
  # absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.", units = "years", type = "continuous",
      notes = paste(
        "Median 42 years (IQR 35-56), Nguyen 2026 Table 1. Screened against the structural",
        "parameters and not retained (Results 3.2, Supplementary Figure S8)."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (renal-function marker).",
      units = "mL/min (raw, NOT BSA-normalized)", type = "continuous",
      notes = paste(
        "Median 68.1 mL/min (IQR 62.8-79.3), Nguyen 2026 Table 1, printed as 'CLCR (mL/min)'.",
        "The estimating equation is not stated in the paper. All patients had normal renal",
        "function (Results 3.1), so the cohort carries little information about renal impairment.",
        "Screened and not retained (Results 3.2, Supplementary Figure S8)."
      )
    ),
    CREAT = list(
      description = "Serum creatinine.", units = "umol/L", type = "continuous",
      notes = paste(
        "Median 77.9 umol/L (IQR 75.9-84), Nguyen 2026 Table 1. Screened and not retained",
        "(Results 3.2, Supplementary Figure S8)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase (hepatic-function marker).",
      units = "U/L", type = "continuous",
      notes = paste(
        "Median 20 U/L (IQR 12.4-53), Nguyen 2026 Table 1. Screened and not retained",
        "(Results 3.2, Supplementary Figure S8)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase (hepatic-function marker).",
      units = "U/L", type = "continuous",
      notes = paste(
        "Median 24 U/L (IQR 16-46), Nguyen 2026 Table 1. Screened and not retained",
        "(Results 3.2, Supplementary Figure S8)."
      )
    )
    # Platelet count (Table 1 median 260 x10^9/L) and haemoglobin (a column in
    # the Table S7 $INPUT list) are NOT listed here. The paper describes the
    # covariate screen as covering "body weight, age, renal, and hepatic
    # function markers" (Methods 2.3) and reports that "age, renal function, and
    # liver enzymes" were not significant (Results 3.2) -- neither statement
    # covers the haematology markers, which are reported as baseline
    # characteristics and linezolid myelosuppression-safety markers rather than
    # as screened covariates. Recorded in population$notes instead of being
    # claimed as a screened-and-rejected effect.
  )

  population <- list(
    species        = "human",
    n_subjects     = 17,
    n_studies      = 1,
    age_median     = "42 years (IQR 35-56)",
    weight_median  = "50 kg (IQR 45-56)",
    weight_range   = "IQR 45-56 kg; the full range is not reported",
    height_median  = "162 cm (IQR 160-165)",
    bmi_median     = "19.1 kg/m^2 (IQR 16.5-20.5)",
    sex_female_pct = 17.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "multidrug-resistant pulmonary tuberculosis (MDR-TB)",
    renal_function = "normal in all patients; CLcr median 68.1 mL/min (IQR 62.8-79.3), serum creatinine median 77.9 umol/L (IQR 75.9-84)",
    hepatic_function = "normal in all patients; ALT median 20 U/L (IQR 12.4-53), AST median 24 U/L (IQR 16-46)",
    dose_range     = "oral linezolid at steady state, median daily dose 600 mg (range 450-600 mg)",
    regions        = "Vietnam (four provinces)",
    notes          = paste(
      "Pharmacokinetic sub-study of the V-SMART trial (ACTRN12620000681954), prospective and",
      "observational. Of 28 enrolled patients, 17 had evaluable PK data; 102 paired saliva-plasma",
      "samples were drawn at pre-dose, 2 h and 5 h post-dose after at least 7 days of treatment",
      "(i.e. at steady state). Baseline demographics are Nguyen 2026 Table 1. Assay LLOQ was",
      "0.5 mg/L in both matrices, with values below it treated as missing (Beal M1). Sex was NOT",
      "assessed as a covariate because of limited representation (14 male vs 3 female). Unbound",
      "concentrations were not measured; the paper assumes an unbound fraction of 87.2% when",
      "converting to the fAUC(0-24)/MIC >= 125 efficacy target. Total daily doses spanned",
      "450-600 mg, so the published mean NCA exposures in Table S6 correspond to a mixed-dose",
      "cohort rather than to a uniform 600 mg regimen."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma disposition. Nguyen 2026 Table 2 "Final model" column
    # (OFV = 275.423). Parameters are APPARENT (CL/F, V/F) because
    # absolute bioavailability was not identifiable from the oral-only
    # design; F was fixed to 1 as the anchor.
    # ------------------------------------------------------------------
    lcl <- log(4.35); label("Apparent clearance, CL/F (L/h)")                  # Table 2 final model theta_CL = 4.35 L/h (RSE 11%; SIR median 4.42, 95% CI 3.48-5.46)
    lvc <- log(46.8); label("Apparent central volume of distribution, V/F (L)") # Table 2 final model theta_V = 46.8 L (RSE 5%; SIR median 46.98, 95% CI 42.35-51.63)

    # Absorption rate constant, FIXED. Estimated at 1.51 1/h but with %RSE
    # 358% and high correlation with CL/F and V/F given no absorption-phase
    # data, so the authors fixed it to stabilise the model, matching the
    # bootstrap median of Alghamdi 2020 (1.53, 95% CI 0.98-3.08).
    # Supplementary Table S3 shows CL/F and V/F are insensitive to Ka over
    # 1.4-2.4 1/h.
    lka <- fixed(log(1.51)); label("First-order absorption rate constant, Ka (1/h)") # Table 2 final model Ka = 1.51 (Fixed); Table S7 $THETA (1.51) FIX

    # Oral bioavailability, FIXED to 1. This is the anchor that makes CL and
    # V apparent rather than absolute; it is not a claim that F is unity.
    lfdepot <- fixed(log(1)); label("Oral bioavailability, F (fraction)")       # Table 2 final model F = 1 (Fixed); Table S7 $THETA (1) FIX

    # ------------------------------------------------------------------
    # Saliva effect-compartment rate constants. Nguyen 2026 Table 2 final
    # model; Figure 1 names them Kabs / Kreabs / Kel. Canonical names use
    # the registered kin_<compartment> / kout_<compartment> tissue-exchange
    # family plus kel_<compartment> for the irreversible salivary loss --
    # the paper's own symbol for that term is "Kel", which collides with the
    # canonical central elimination rate constant, so the suffix is what
    # keeps them apart.
    # ------------------------------------------------------------------
    lkin_saliva <- log(4.93); label("Central-to-saliva secretion rate constant, Kabs (1/h)")   # Table 2 final model K_abs = 4.93 (RSE 2%; SIR median 4.92, 95% CI 4.80-4.99); Table S7 K23
    lkout_saliva <- log(1.84); label("Saliva-to-central reabsorption rate constant, Kreabs (1/h)") # Table 2 final model K_reabs = 1.84 (RSE 8%; SIR median 1.85, 95% CI 1.04-2.66); Table S7 K32
    lkel_saliva <- log(2.13); label("Irreversible salivary elimination rate constant, Kel (1/h)")  # Table 2 final model K_el = 2.13 (RSE 8%; SIR median 2.14, 95% CI 1.30-2.91); Table S7 K30

    # ------------------------------------------------------------------
    # Covariate effect. Table 2 prints the relation in the row header as
    # "V/F [L]  theta_V . (WT/50)^theta_COV".
    # ------------------------------------------------------------------
    e_wt_vc <- 1.1; label("Power exponent for body weight on V/F (unitless)")   # Table 2 final model theta_COV = 1.1 (RSE 4%; SIR median 1.09, 95% CI 0.62-1.46); d-2LL = -10.16, p < 0.01

    # ------------------------------------------------------------------
    # Interindividual variability. Exponential, P_i = P_TV * exp(eta_i)
    # (Supplementary Information S1.1 Eq. 1).
    #
    # VARIANCE CONVENTION. Table 2 labels the IIV row "omega_CL (CV%)" but
    # the tabulated number is sqrt(variance) * 100, NOT the log-normal
    # sqrt(exp(omega^2) - 1) * 100. Proof from Table S7, whose $THETA / $OMEGA
    # initials are the SIR medians of the final fit: $OMEGA 0.228 gives
    # sqrt(0.228) = 0.4775 -> 47.75%, against the reported SIR median of
    # 47.77%. The same identity holds for both $SIGMA rows (0.07 -> 26.46%
    # vs SIR 26.62%; 0.13 -> 36.06% vs SIR 37.01%). So the final-model
    # omega_CL variance is 0.4472^2 = 0.19999.
    # ------------------------------------------------------------------
    etalcl ~ 0.19999  # Table 2 final model omega_CL = 44.72 CV% -> 0.4472^2 (RSE 37%, shrinkage 0.1%)

    # IIV on V/F, Ka, F, Kabs, Kreabs and Kel is NOT carried. Table S7 fixes
    # all six to zero ($OMEGA "0 FIX" for BSVV, BSVKA, BSVF1, BSV K23,
    # BSV K32, BSV K30) and Results 3.2 explains why: "Interindividual
    # variability (IIV) could not be reliably estimated for saliva-specific
    # rate constants due to high shrinkage (> 55% for all of them)". They are
    # omitted rather than written as `~ fixed(0)` because a zero-variance
    # diagonal makes OMEGA singular and breaks the Cholesky sampler used by
    # rxSolve.

    # ------------------------------------------------------------------
    # Residual unexplained variability. Proportional in both matrices,
    # modelled separately (Methods 2.3: "Residual error in saliva
    # concentrations was modelled separately"). Values are the Table 2 CV%
    # read as sqrt(sigma) per the convention proof above.
    # ------------------------------------------------------------------
    propSd <- 0.2586; label("Proportional residual SD for plasma Cc (fraction)")              # Table 2 final model sigma_Plasma_Prop = 25.86 CV% (RSE 23%, shrinkage 9%)
    propSd_Csaliva <- 0.3591; label("Proportional residual SD for saliva Csaliva (fraction)") # Table 2 final model sigma_Saliva_Prop = 35.91 CV% (RSE 20%, shrinkage 4%)
  })

  model({
    # 1. Individual parameters. The weight term carries the population
    #    median 50 kg as its normalisation constant (Table 1; Table S7
    #    VWT = ((WT/50)**THETA(8))).
    wt_ref <- 50  # kg

    cl <- exp(lcl + etalcl)
    vc <- exp(lvc) * (WT / wt_ref)^e_wt_vc
    ka <- exp(lka)

    kin_saliva <- exp(lkin_saliva)
    kout_saliva <- exp(lkout_saliva)
    kel_saliva <- exp(lkel_saliva)

    # 2. Micro-constant
    kel <- cl / vc

    # 3. ODE system. NONMEM ADVAN13 TRANS1 with NCOMP = 3
    #    (DEPOT / CENT / SALIVA), Table S7.
    #
    #    The saliva state is a DRIVEN (non-depleting) hypothetical effect
    #    compartment, exactly as Figure 1's caption states: "The saliva
    #    bio-compartment is a hypothetical effect compartment". The central
    #    equation therefore carries no -kin_saliva*central loss and no
    #    +kout_saliva*saliva return.
    #
    #    This is load-bearing and was settled against the paper's own printed
    #    numbers, because Table S7 is headed "Example of NONMEM model code"
    #    and omits the $DES block that its own ADVAN13 subroutine requires.
    #    Had the three rate constants instead been read as a mass-balance
    #    exchange with central -- the reading the K23/K32/K30 naming suggests
    #    -- salivary loss would dominate total elimination: at steady state
    #    kel_saliva * A_saliva = 2.13 * 1.242 = 2.645 per hour against
    #    kel = 4.35/46.8 = 0.0929 per hour, giving a 28-fold larger total
    #    clearance. Solved in rxode2, the mass-balance reading returns a
    #    steady-state plasma AUC(0-24) of 4.7 mg*h/L and a Cmax of 2.8 mg/L,
    #    against the published 131.5 mg*h/L and 11.2 mg/L (Table S6). The
    #    driven reading below returns AUC(0-24) = Dose/(CL/F) = 137.9 mg*h/L
    #    and, on the paper's own hourly 0-24 h NCA grid, a cohort mean
    #    half-life of 8.03 h and Tmax of 1.98 h against the published 8.01 h
    #    and 1.98 h. See the vignette for the full arithmetic.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(saliva) <- kin_saliva * central - (kout_saliva + kel_saliva) * saliva

    # 4. Bioavailability (fixed to 1; retained so the apparent-parameter
    #    anchor is visible and so a user re-fitting the model has the hook).
    f(depot) <- exp(lfdepot)

    # 5. Observations. Both matrices are scaled by the central volume:
    #    Table S7 prints S2 = V and CP = A(2)/V for plasma, and CS = A(3)
    #    for saliva. The literal CS = A(3) (an implied 1 L saliva scale) is
    #    falsified by the paper: it puts the steady-state saliva:plasma
    #    exposure ratio at V * kin_saliva/(kout_saliva + kel_saliva) =
    #    46.8 * 1.242 = 58.1, i.e. saliva running ~58x plasma (~650 mg/L
    #    against a plasma Cmax of 11.2 mg/L), whereas Figure 3 plots both
    #    matrices on one shared 0-40 mg/L axis with the saliva median only
    #    slightly above plasma. Dividing by vc gives the parameter-free
    #    steady-state ratio kin_saliva/(kout_saliva + kel_saliva) =
    #    4.93/(1.84 + 2.13) = 1.2418, against the independently reported
    #    "mean saliva-to-plasma ratio ... 1.27 (95% CI 1.09-1.44)"
    #    (Discussion) -- agreement to 2%. Operator-ratified (sidecar
    #    request-001 q2, 2026-09-02); see the vignette Errata.
    Cc <- central / vc
    Csaliva <- saliva / vc

    Cc ~ prop(propSd)
    Csaliva ~ prop(propSd_Csaliva)
  })
}

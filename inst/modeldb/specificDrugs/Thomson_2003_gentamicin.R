Thomson_2003_gentamicin <- function() {
  description <- "One-compartment population PK model of intramuscular gentamicin in African infants with suspected severe sepsis (Thomson 2003). The 8 mg/kg i.m. dose is modelled as an IV bolus into the central compartment because first-order absorption could not be characterised from the sparse 1 h / next-morning sampling (the paper documents that ka estimates were poorly identified and absorption appeared complete by 1 h). Apparent clearance scales linearly with body weight and as a power function of (postnatal age + 1 day) normalised to the cohort median + 1 day; apparent volume of distribution scales linearly with body weight relative to the cohort median 3 kg. Reported CL and V are apparent values (CL/F, V/F) because all doses were administered by intramuscular injection and bioavailability could not be estimated."
  reference   <- "Thomson AH, Kokwaro GO, Muchohi SN, English M, Mohammed S, Edwards G. Population pharmacokinetics of intramuscular gentamicin administered to young infants with suspected severe sepsis in Kenya. Br J Clin Pharmacol. 2003;56(1):25-31. doi:10.1046/j.1365-2125.2003.01819.x"
  vignette    <- "Thomson_2003_gentamicin"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the time of the dose.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Thomson 2003 cohort median weight 3.05 kg (range 1.24-6.72 kg, Table 1). The volume-of-distribution effect uses a linear-deviation form anchored at 3 kg: V = 2.02 * (1 + 0.277 * (WT - 3)). The clearance effect is strictly linear in WT: CL = 0.0913 * WT * f_PNA. Time-varying within an episode is in principle allowed but the source data set carried a single baseline weight per infant.",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological since birth).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Thomson 2003 reports postnatal age in DAYS (cohort median 10 days, range 0-99 days, Table 1). The canonical PNA column is in MONTHS, so the model `model()` block reparameterises the source expression: the paper's `CL = theta1 * WT * ((PNA_days + 1) / 11)^theta2` is encoded as `CL = exp(lcl + etalcl) * WT * ((PNA_months + 1/30.4375) / (11/30.4375))^e_pna_cl`. Numerically identical because the 30.4375 days/month factor cancels in the ratio. The +1 day shift (in months: 1/30.4375 = 0.03285) was added by Thomson 2003 to keep the expression defined for PNA = 0 day cohort members. Users should supply PNA in months in the dataset.",
      source_name        = "PNA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 107L,
    n_studies      = 1L,
    n_observations = 203L,
    age_range      = "PNA 0-99 days",
    age_median     = "PNA 10 days",
    weight_range   = "1.24-6.72 kg",
    weight_median  = "3.05 kg",
    sex_female_pct = 36.4,
    race_ethnicity = c(Black_African = 100),
    disease_state  = "Young infants admitted to Kilifi District Hospital, Kenya, with suspected severe sepsis. Inclusion required body weight >= 1 kg, age < 3 months (changed to < 2 months from February 2001), no prior gentamicin exposure, no anuria for >= 24 h, and creatinine clearance on admission within an acceptable range. Children with primary diagnosis of tetanus or major life-threatening congenital malformations were excluded.",
    dose_range     = "Single intramuscular dose of 8 mg/kg gentamicin (the first dose of a once-daily regimen). The model was developed only on data following this first dose.",
    n_peak_samples       = 97L,
    n_nextday_samples    = 106L,
    peak_sample_time     = "median 1.05 h after dose (range 0.49-2.9 h)",
    nextday_sample_time  = "median 16.87 h after dose (range 8.35-32.85 h)",
    creatinine_range     = "17-173 micromol/L (median 52)",
    wbc_range            = "2.1-51.6 x 10^9/L (median 11.3)",
    haemoglobin_range    = "4.2-21.5 g/dL (median 14.5)",
    regions              = "Coastal Kenya (Kilifi District Hospital). Recruitment August 2000 - April 2001. Approval from the national Kenyan Ethical Committee.",
    notes                = "Baseline demographics from Thomson 2003 Table 1 (n = 107). Initial data set had 124 patients and 238 concentration measurements; 14 with missing creatinine, 2 with spurious results, and 1 below LLOQ on the only sample were excluded. The 8 mg/kg i.m. dose produced a median 1 h peak of 10.6 mg/L (range 3.0-19.8) and median next-day concentration of less than 2 mg/L (range 0.3-6.2). Gentamicin was assayed by FPIA (Abbott TDx FLx) with LLOQ 0.27 mg/L."
  )

  ini({
    # Final-model estimates from Thomson 2003 Table 4 ("Parameter estimates
    # from the best covariate model (covariance model results shown)").
    # The structural form is:
    #   CL/F (L/h) = theta1 * WT (kg) * ((PNA_days + 1) / 11)^theta2
    #   V/F  (L)   = theta3 * (1 + theta4 * (WT - 3 kg))
    # so the population estimate for a child of median weight (3 kg) and
    # median age (10 days) is CL/F = 0.274 L/h and V/F = 2.02 L (matches
    # the Abstract Results paragraph).
    #
    # Reported parameters are apparent (CL/F, V/F) because all doses were
    # i.m. and bioavailability could not be estimated. The vignette
    # documents the apparent-parameter caveat.
    #
    # Reference points for the canonical-PNA reparameterisation:
    #   pna_shift_months = 1 / 30.4375 = 0.03285 months  (1 day shift)
    #   pna_ref_months   = 11 / 30.4375 = 0.36139 months (median 10 d + 1 d shift)
    # The expression ((PNA_months + pna_shift_months) / pna_ref_months)
    # equals 1 at the cohort median (10 days = 0.32854 months) because
    # (0.32854 + 0.03285) / 0.36139 = 1.

    lcl      <- log(0.0913); label("Apparent clearance coefficient per kg body weight at reference postnatal age (cohort median 10 days + 1 day shift = 11 days) (L/h/kg)")  # Thomson 2003 Table 4: theta1 = 0.0913 L/h/kg (RSE 4.2%)
    e_pna_cl <- 0.130;       label("Power exponent of (PNA + 1 day) / 11 days on apparent clearance (unitless)")                                                            # Thomson 2003 Table 4: theta2 = 0.130 (RSE 24%)
    lvc      <- log(2.02);   label("Apparent volume of distribution intercept at WT = 3 kg (L)")                                                                            # Thomson 2003 Table 4: theta3 = 2.02 L (RSE 4.0%)
    e_wt_vc  <- 0.277;       label("Linear-deviation slope of body weight on apparent volume per kg above 3 kg (per kg)")                                                   # Thomson 2003 Table 4: theta4 = 0.277 (RSE 11%)

    # Inter-individual variability (Thomson 2003 Table 4 IIV column,
    # covariance-model results). The source reports IIV as a percentage
    # coefficient of variation under a log-normal model; the equivalent
    # log-scale variance is omega^2 = log(1 + CV^2).
    #   IIV CL : 42% CV -> log(1 + 0.42^2) = 0.16127
    #   IIV V  : 40% CV -> log(1 + 0.40^2) = 0.14842
    # The paper's text notes that the basic model "required a term for
    # covariance between CL and V" but the off-diagonal covariance value
    # is not reported numerically in Table 2 or Table 4; only the
    # diagonal IIVs are tabulated. The model file therefore encodes
    # independent etas and documents this gap in the vignette
    # Assumptions and deviations section.
    etalcl ~ 0.16127  # Thomson 2003 Table 4: IIV CL = 42% CV (RSE 12%)
    etalvc ~ 0.14842  # Thomson 2003 Table 4: IIV V  = 40% CV (RSE 15%)

    # Residual error. Thomson 2003 fixed the additive residual SD to a
    # negligible value (0.0001 mg/L) because the FOCE algorithm with at
    # most two samples per patient achieved perfect agreement between
    # measured and individual-predicted concentrations, which meant
    # residual error could not be characterised structurally. Fixing
    # residual error to a near-zero value preserved the OFV and parameter
    # estimates while allowing parameter standard errors to be returned.
    # The vignette flags the loss of any meaningful concentration-level
    # error term as an assumption.
    addSd <- fixed(0.0001); label("Additive residual error (mg/L) -- fixed to a negligible value in Thomson 2003 because residual error could not be characterised under FOCE with at most 2 samples per patient")  # Thomson 2003 Table 4: residual error 0.0001 mg/L FIXED
  })

  model({
    # ----- Derived covariate terms -----
    # Postnatal-age power factor on apparent clearance. Thomson 2003
    # expresses PNA in days with a 1-day shift and a reference of 11 days
    # (cohort median 10 d + 1 d shift). The canonical PNA column is in
    # months, so the same expression is encoded by converting the 1 d
    # shift and 11 d reference into months. Numerator and denominator
    # carry the same units factor, so the ratio is unchanged.
    pna_shift_months <- 1 / 30.4375
    pna_ref_months   <- 11 / 30.4375
    f_pna_cl <- ((PNA + pna_shift_months) / pna_ref_months) ^ e_pna_cl

    # ----- Individual PK parameters -----
    # Apparent clearance: linear in WT, power in (PNA + shift) / ref.
    cl <- exp(lcl + etalcl) * WT * f_pna_cl

    # Apparent volume of distribution: linear-deviation in WT anchored at
    # the cohort median 3 kg. At WT = 3 kg, V = exp(lvc + etalvc).
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 3))

    # Elimination microconstant.
    kel <- cl / vc

    # ----- One-compartment IV-bolus disposition -----
    # Thomson 2003 modelled the 8 mg/kg i.m. dose as an i.v. bolus input
    # to the central compartment because the first-order absorption rate
    # could not be characterised (absorption was complete by 1 h and the
    # sparse design provided too few early peak samples to identify ka).
    # Users dose into the central compartment in the event table.
    d/dt(central) <- -kel * central

    # ----- Observation and error -----
    # Plasma gentamicin concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}

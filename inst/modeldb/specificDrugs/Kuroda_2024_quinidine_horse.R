Kuroda_2024_quinidine_horse <- function() {
  description <- paste(
    "Veterinary (Thoroughbred horse). Three-compartment population PK model for the",
    "Vaughan-Williams class Ia antiarrhythmic quinidine (QND) in Thoroughbred",
    "horses, fitted by nonlinear mixed-effects modelling (Phoenix WinNonlin/NLME",
    "8.4, QRPEM engine) to pooled plasma data from 10 healthy horses (single",
    "5 mg/kg quinidine hydrochloride monohydrate as a 5-min IV infusion and single",
    "or twice-daily 20 mg/kg quinidine sulfate dihydrate by nasogastric tube, with",
    "a two-week washout) and 19 racehorses treated for naturally occurring atrial",
    "fibrillation (9.3-30.6 mg/kg quinidine sulfate dihydrate by nasogastric tube",
    "at veterinarian-chosen intervals). Absorption from the gastrointestinal depot",
    "is first order (Kabs) with a bioavailability factor F held on the logit scale;",
    "disposition is a central compartment (V1) connected to a rapidly equilibrating",
    "peripheral compartment (V2, distribution clearance CL2) and a slowly",
    "equilibrating peripheral compartment (V3, distribution clearance CL3), with",
    "elimination clearance CL from V1. Every structural parameter is reported per",
    "kilogram of body weight in Kuroda 2024 Table 1 (L/kg, L/kg/h), so the model",
    "scales each of them linearly with WT (exponent 1 fixed); doses are therefore",
    "given in mg and concentrations come out in ug/mL. All doses and parameters are",
    "expressed as QND BASE: the paper converts the administered salts with factors",
    "of 1.206 (sulfate dihydrate) and 1.168 (hydrochloride monohydrate). Between",
    "subject variability is exponential on the six disposition parameters and Kabs",
    "and logit-normal on F; Kuroda 2024 states a full OMEGA matrix was used but",
    "reports only the diagonal (the BSV% column), so only variances are packaged.",
    "The residual model is combined proportional plus additive with a separate pair",
    "estimated for the IV and the oral datasets; the ROUTE_IV indicator selects",
    "between them. No covariate (condition, age, body weight or sex) reached the",
    "paper's BIC threshold, so none is carried on any structural parameter."
  )
  reference <- paste(
    "Kuroda T, Minamijima Y, Kinman CK, Takahashi Y, Ebisuda Y, Inoue K,",
    "Ishikawa H, Mita H, Tamura N, Nukada T, Toutain P-L, Ohta M. (2024).",
    "Rational quinidine dosage regimen for atrial fibrillation in Thoroughbred",
    "racehorses based on population pharmacokinetics.",
    "Frontiers in Veterinary Science 11:1454342.",
    "doi:10.3389/fvets.2024.1454342.",
    sep = " "
  )
  vignette <- "Kuroda_2024_quinidine_horse"

  # Nonstandard residual-error SD names: Kuroda 2024 Table 1 estimates a
  # SEPARATE proportional and additive pair for the IV dataset (CMultStdev0,
  # Stdev0) and the oral dataset (CMultStdev1, Stdev1). Declared here so
  # checkModelConventions() accepts the route-suffixed names. Same pattern as
  # Ahmed_2015_topiramate.R.
  paper_specific_residual_sds <- c("propSdOral", "propSdIv", "addSdOral", "addSdIv")

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Kuroda 2024 section 2.5, which names
  # V1 the central compartment, V2 and V3 the peripheral compartments, and
  # adds Kabs and F "to the model for PO administration".
  compartmentData <- list(
    depot       = list(analyte = "quinidine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "quinidine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "quinidine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "quinidine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed adult-horse body weight. Kuroda 2024 Table 1 reports every structural volume and",
        "clearance per kilogram (L/kg, L/kg/h), so WT enters as a linear multiplier (exponent 1) on cl,",
        "vc, q, vp, q2 and vp2 rather than as an estimated allometric term; Kabs (1/h) and F (%) are not",
        "weight-normalized and are used as-is. This per-kilogram normalization IS the base model, not a",
        "covariate effect: the stepwise search described in section 2.5 tested body weight as an",
        "ADDITIONAL effect on top of it and rejected it (Discussion: 'None of the covariates explored",
        "(age, BW, sex, and presence/absence of AF) were significant with BIC <10.0'). Study horses",
        "weighed 473-563 kg (healthy) and 430-540 kg (atrial fibrillation).",
        sep = " "
      ),
      source_name        = "BW (body weight; Kuroda 2024 section 2.1)"
    ),
    ROUTE_IV = list(
      description        = "Indicator for intravenous administration of quinidine hydrochloride monohydrate",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral quinidine sulfate dihydrate by nasogastric tube)",
      notes              = paste(
        "Per-observation dosing-route indicator (1 = 5-min IV infusion, 0 = oral administration by",
        "nasogastric tube). Kuroda 2024 retained no route covariate effect on any structural parameter,",
        "but Table 1 estimates a separate proportional AND additive residual-error pair per dataset:",
        "CMultStdev0 = 0.0594 with Stdev0 = 0.0338 for IV, CMultStdev1 = 0.1571 with Stdev1 = 0.0479",
        "for oral. The model body selects between propSdIv/propSdOral and addSdIv/addSdOral via",
        "ROUTE_IV. Distinct from the rxode2 cmt event column (cmt = central for the IV infusion,",
        "cmt = depot for oral doses), which controls where the dose enters, not the residual magnitude.",
        sep = " "
      ),
      source_name        = "dataset index 0 (IV) / 1 (PO) in the Phoenix NLME residual-error parameter names Stdev0 / Stdev1"
    )
  )

  # Screened by the paper's stepwise BIC covariate search (section 2.5:
  # "covariates were tested for condition (healthy or AF), age, BW, and sex")
  # but NOT retained in the final model, and reported with no point estimate.
  # Body weight is deliberately absent from this list -- it IS referenced in
  # model() as the per-kilogram scaler, so it lives in covariateData above.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a candidate covariate on every structural parameter (section 2.5) but not retained; no effect size is reported. Horses were 2-7 years old (healthy) and 2-10 years old (atrial fibrillation)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened as a candidate covariate (section 2.5, 'sex') but not retained; no effect size is reported. The cohort was five stallions and five mares (healthy) and 12 stallions and seven mares (atrial fibrillation)."
    ),
    DIS_HEALTHY = list(
      description = "Healthy-horse indicator (complement is naturally occurring atrial fibrillation)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a candidate covariate (section 2.5, 'condition (healthy or AF)') but not retained; no effect size is reported. Because the disease state was rejected, the Discussion concludes 'the dosing regimen proposed in this study can probably be applied to the entire Thoroughbred population'."
    )
  )

  population <- list(
    species        = "horse (Thoroughbred)",
    n_subjects     = 29L,
    n_studies      = 2L,
    age_range      = "2-7 years (healthy); 2-10 years (atrial fibrillation)",
    weight_range   = "473-563 kg (healthy); 430-540 kg (atrial fibrillation)",
    sex_female_pct = 41.4,
    disease_state  = "healthy (n = 10) or naturally occurring atrial fibrillation under quinidine therapy (n = 19); 18 of the 19 converted to sinus rhythm",
    dose_range     = paste(
      "Healthy horses: single 5 mg/kg quinidine hydrochloride monohydrate (4.28 mg/kg QND base) as a",
      "5-min IV infusion into the right jugular vein in 500 mL sterile saline, and single 20 mg/kg",
      "quinidine sulfate dihydrate (16.58 mg/kg QND base) in 500 mL water by nasogastric tube with a",
      "500 mL water flush, in a two-way crossover with a two-week washout (n = 6); the same 20 mg/kg",
      "oral dose given twice q 6 h (n = 4). Atrial-fibrillation horses: 9.3-30.6 mg/kg quinidine",
      "sulfate dihydrate (7.7-25.4 mg/kg QND base) by nasogastric tube at veterinarian-chosen doses",
      "and intervals.",
      sep = " "
    ),
    sampling       = paste(
      "IV arm: pre-dose and 0, 5, 10, 20, 30, 45 min and 1, 2, 3, 4, 6, 8, 12 h. Single oral arm:",
      "pre-dose and 30 min, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 12, 24 h. Twice-oral arm: pre-dose and",
      "30 min, 1, 1.5, 2, 2.5, 3, 5.8, 6.5, 7, 7.5, 8, 8.5, 9, 12, 24 h after the first dose.",
      "Atrial-fibrillation horses: pre-dose and 0.5, 1, 2, 3, 4 h after administration, with exact",
      "dosing and sampling clock times recorded. Approximately 10 mL of blood from the right jugular",
      "vein into heparinized tubes; LC-MS/MS assay (Shimadzu LC / SCIEX MS, Acquity BEH column,",
      "quinidine-d3 internal standard) calibrated over 0.03-10.0 ug/mL with a lower limit of",
      "quantitation of 0.03 ug/mL. Values below the limit of quantitation, under 5% of the data,",
      "were excluded from the model.",
      sep = " "
    ),
    regions        = "Japan (Equine Research Institute, Japan Racing Association, Shimotsuke; Miho and Ritto Training Center racehorse hospitals)",
    notes          = paste(
      "Section 2.1 describes 10 healthy horses plus 19 horses with atrial fibrillation, and the",
      "Abstract states that 'the data from 29 horses were modeled', so n_subjects is 29. The Table 1",
      "caption instead says 'in 27 horses'; the paper never reconciles the two or names an exclusion,",
      "and 29 is supported by both the Abstract and the Animals section, so 29 is used here. See",
      "vignette Errata. Onset date of atrial fibrillation was unknown in 8 horses; the remaining 11",
      "started quinidine 4.6 +/- 1.5 days after onset. The therapeutic window established by the paper",
      "is 2.0-3.8 ug/mL: the median plasma quinidine concentration at conversion to sinus rhythm was",
      "2.0 ug/mL (range 0.5-2.7, n = 13 horses whose conversion fell inside the sampling window) and",
      "the median concentration at which adverse effects occurred was 3.8 ug/mL (range 1.6-5.1, two",
      "healthy horses after IV and four atrial-fibrillation horses after oral dosing). Ethics approvals",
      "22-8 and 23-6 (Institutional Animal Care and Use Committee, Equine Research Institute, Japan",
      "Racing Association); owner consent obtained. All eta shrinkage values were < 0.3.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Structural disposition parameters -- Kuroda 2024 Table 1
    # =====================================================================
    # Table 1 ("Bootstrap estimates of typical (median) population primary
    # and secondary parameters of quinidine") reports every primary
    # structural parameter normalized to body weight (L/kg for volumes,
    # L/kg/h for clearances). They are carried here on that per-kilogram
    # basis and multiplied by WT (exponent 1) in model(), which is exactly
    # what the per-kilogram normalization means.
    #
    # In Table 1 the CV% column is the BOOTSTRAP PRECISION of the typical
    # value (with the 2.5-97.5% percentile columns), NOT between-subject
    # variability; BSV is the separate final column. Both are quoted per
    # parameter below.
    #
    # Naming map (Kuroda 2024 Table 1 -> nlmixr2lib canonical):
    #   V1  -> vc            (central compartment)
    #   V2  -> vp            (rapidly equilibrating peripheral, via CL2)
    #   V3  -> vp2           (slowly equilibrating peripheral,   via CL3)
    #   CL  -> cl            (plasma clearance out of V1)
    #   CL2 -> q             (V1 <-> V2 distribution clearance)
    #   CL3 -> q2            (V1 <-> V3 distribution clearance)
    #   Kabs -> ka           (first-order oral absorption rate constant)
    #   F   -> fdepot        (oral bioavailability, on the logit scale)
    lvc <- log(0.63); label("Weight-normalized central volume of distribution V1 (L/kg)")                        # Table 1: V1 = 0.63 L/kg (bootstrap CV% 12.4; 2.5-97.5% 0.49-0.80)
    lvp <- log(0.59); label("Weight-normalized rapidly equilibrating peripheral volume V2 (L/kg)")               # Table 1: V2 = 0.59 L/kg (bootstrap CV% 12.1; 2.5-97.5% 0.49-0.72)
    lvp2 <- log(3.68); label("Weight-normalized slowly equilibrating peripheral volume V3 (L/kg)")               # Table 1: V3 = 3.68 L/kg (bootstrap CV% 7.0; 2.5-97.5% 3.09-4.13)
    lcl <- log(0.49); label("Weight-normalized plasma clearance CL (L/kg/h)")                                   # Table 1: CL = 0.49 L/kg/h (bootstrap CV% 6.4; 2.5-97.5% 0.43-0.56)
    lq <- log(2.87); label("Weight-normalized distribution clearance to peripheral1 CL2 (L/kg/h)")              # Table 1: CL2 = 2.87 L/kg/h (bootstrap CV% 15.2; 2.5-97.5% 2.08-3.68)
    lq2 <- log(2.44); label("Weight-normalized distribution clearance to peripheral2 CL3 (L/kg/h)")             # Table 1: CL3 = 2.44 L/kg/h (bootstrap CV% 13.9; 2.5-97.5% 1.91-3.22)

    # =====================================================================
    # Oral absorption -- Kuroda 2024 Table 1
    # =====================================================================
    # Section 2.5: "The absorption rate constant (Kabs) and the
    # bioavailability factor (F) were added to the model for PO
    # administration." Kabs cross-checks against the Table 1 secondary
    # parameter absorption half-life: log(2)/1.00 = 0.693 h vs 0.69 h.
    lka <- log(1.00); label("First-order oral absorption rate constant Kabs (1/h)")                             # Table 1: Kabs = 1.00 1/h (bootstrap CV% 24.5; 2.5-97.5% 0.67-1.70); secondary parameter absorption half-life = 0.69 h
    # Table 1 footnote: "For F, an ilogit transformation was used to prevent
    # estimates higher than 100%." F is therefore held on the logit scale
    # here and mapped back with expit() in model(), so that F stays inside
    # (0, 1) for every random draw; qlogis(0.364) = -0.55793. Encoding F as
    # exp(lfdepot + eta) instead would put roughly 0.1% of simulated horses
    # above 100% bioavailability, which is exactly what the paper's
    # transformation exists to prevent.
    logitfdepot <- qlogis(0.364); label("Logit of oral bioavailability F (unitless; F = 0.364)")                # Table 1: F = 36.4% (bootstrap CV% 2.9; 2.5-97.5% 34.1-38.1%)

    # =====================================================================
    # Between-subject variability -- Kuroda 2024 Table 1, BSV% column
    # =====================================================================
    # Section 2.5 Equation 1 gives an exponential BSV model,
    # theta_i = theta_tv * exp(eta_i), with eta ~ N(0, omega^2), and
    # Equation 2 defines the reported BSV% as
    #     CV(%) = 100 * sqrt(exp(omega^2) - 1)
    # so each variance below is recovered as
    #     omega^2 = log(1 + (BSV% / 100)^2).
    # Worked from the Table 1 BSV% column:
    #   V1   40.5% -> log(1 + 0.405^2) = 0.151884
    #   V2   37.8% -> log(1 + 0.378^2) = 0.133555
    #   V3   28.5% -> log(1 + 0.285^2) = 0.078095
    #   CL   25.6% -> log(1 + 0.256^2) = 0.063478
    #   CL2  74.9% -> log(1 + 0.749^2) = 0.445327
    #   CL3  47.3% -> log(1 + 0.473^2) = 0.201903
    #   Kabs 94.3% -> log(1 + 0.943^2) = 0.636179
    #   F    33.1% -> log(1 + 0.331^2) = 0.103964
    # F's variance is on the LOGIT scale, matching the ilogit
    # transformation of its typical value; Phoenix reports BSV% through
    # Equation 2 regardless of the structural transformation used.
    #
    # Section 2.5 also states "A full OMEGA matrix was used to determine the
    # random components of the model", but Table 1 reports ONLY the diagonal
    # (the BSV% column) -- no correlations, covariances or an OMEGA block
    # appear anywhere in the paper or its figures. Per the standing policy
    # for unreported variance components, only the reported variances are
    # packaged and the off-diagonals are left at zero rather than invented.
    # Simulated between-subject spread from this file is therefore slightly
    # narrower than the paper's Monte Carlo simulations, which used the full
    # matrix. See vignette Errata.
    etalvc ~ 0.151884         # Table 1 BSV% = 40.5 for V1
    etalvp ~ 0.133555         # Table 1 BSV% = 37.8 for V2
    etalvp2 ~ 0.078095        # Table 1 BSV% = 28.5 for V3
    etalcl ~ 0.063478         # Table 1 BSV% = 25.6 for CL
    etalq ~ 0.445327          # Table 1 BSV% = 74.9 for CL2
    etalq2 ~ 0.201903         # Table 1 BSV% = 47.3 for CL3
    etalka ~ 0.636179         # Table 1 BSV% = 94.3 for Kabs
    etalogitfdepot ~ 0.103964 # Table 1 BSV% = 33.1 for F (logit scale)

    # =====================================================================
    # Residual unexplained variability -- Kuroda 2024 Table 1
    # =====================================================================
    # Section 2.5 Equation 4 gives the Phoenix combined form
    #     C(t) = f(theta, Time) * (1 + eps1) + eps2
    # i.e. proportional plus additive. Table 1 estimates a separate pair per
    # dataset: suffix 0 is the IV dataset, suffix 1 the oral dataset. The
    # model body selects between them with the ROUTE_IV indicator.
    #
    # UNITS: Table 1 labels Stdev0/Stdev1 "ug/L", and section 2.5 says the
    # additive sigma was "reported as its standard deviation noted with the
    # same units as plasma concentration (ug/L)". That label is a typo -- the
    # assay's lower limit of quantitation is 0.03 ug/mL and every plotted and
    # tabulated concentration in the paper is in ug/mL, so an additive SD of
    # 0.0338 ug/L would sit three orders of magnitude below the
    # quantification limit and be physically meaningless. Read on the assay's
    # concentration scale, 0.0338 and 0.0479 ug/mL land essentially at the
    # limit of quantitation, which is what a well-behaved additive term does.
    # The values are therefore in ug/mL (= mg/L), which is also what
    # Phoenix's stdev parameter carries. The sibling model
    # Kuroda_2023_cephalothin.R carries the identical typo from the same
    # group and resolves it the same way. See vignette Errata.
    propSdIv <- 0.0594; label("Proportional residual error, intravenous dataset (fraction)")                    # Table 1: CMultStdev0 (residual, proportional for IV) = 0.0594 (bootstrap CV% 28.3; 2.5-97.5% 0.0058-0.0776)
    propSdOral <- 0.1571; label("Proportional residual error, oral dataset (fraction)")                         # Table 1: CMultStdev1 (residual, proportional for PO) = 0.1571 (bootstrap CV% 7.8; 2.5-97.5% 0.1380-0.1853)
    addSdIv <- 0.0338; label("Additive residual error, intravenous dataset (ug/mL)")                            # Table 1: Stdev0 (residual, additive for IV) = 0.0338 (bootstrap CV% 31.6; 2.5-97.5% 0.0137-0.0547); units read as ug/mL, not the printed ug/L
    addSdOral <- 0.0479; label("Additive residual error, oral dataset (ug/mL)")                                 # Table 1: Stdev1 (residual, additive for PO) = 0.0479 (bootstrap CV% 47.5; 2.5-97.5% 0.0008-0.0923); units read as ug/mL, not the printed ug/L
  })

  model({
    # 1. Individual parameters. Kuroda 2024 Table 1 is reported per kilogram
    #    of body weight, so every volume and clearance scales linearly with
    #    WT (exponent 1). Kabs (1/h) and F (%) are not weight-normalized in
    #    Table 1 and are used as-is.
    vc <- exp(lvc + etalvc) * WT
    vp <- exp(lvp + etalvp) * WT
    vp2 <- exp(lvp2 + etalvp2) * WT
    cl <- exp(lcl + etalcl) * WT
    q <- exp(lq + etalq) * WT
    q2 <- exp(lq2 + etalq2) * WT
    ka <- exp(lka + etalka)

    # Oral bioavailability: the random effect is added on the logit scale and
    # mapped back into (0, 1), matching the ilogit transformation recorded in
    # the Table 1 footnote.
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system (Kuroda 2024 section 2.5). An oral dose of quinidine
    #    sulfate dihydrate enters the absorption compartment `depot` and
    #    reaches the central compartment V1 with bioavailability F at
    #    first-order rate Kabs; the IV dose of quinidine hydrochloride
    #    monohydrate is infused directly into `central` over 5 min. Doses
    #    must be supplied as QND BASE in mg (see units and the vignette).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Bioavailability on the oral absorption compartment only; the IV
    #    infusion into `central` is unaffected.
    f(depot) <- fdepot

    # 5. Observation. Cc is total plasma quinidine, expressed as QND base in
    #    mg/L = ug/mL (dose in mg, volume in L).
    Cc <- central / vc

    # 6. Route-conditional residual error (Kuroda 2024 Table 1): each term
    #    collapses to the oral estimate when ROUTE_IV = 0 and to the
    #    intravenous estimate when ROUTE_IV = 1.
    propSd <- propSdOral + (propSdIv - propSdOral) * ROUTE_IV
    addSd <- addSdOral + (addSdIv - addSdOral) * ROUTE_IV

    Cc ~ prop(propSd) + add(addSd)
  })
}

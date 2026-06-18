# Semi-physiological population PK model for dextromethorphan and its three
# phase I metabolites (dextrorphan, 3-methoxymorphinan, 3-hydroxymorphinan)
# in adult breast-cancer patients on chronic tamoxifen, used as a dual
# CYP2D6 / CYP3A phenotypic probe. Extracted from ter Heine 2014 (BJCP)
# Table 2 plus Appendix 1.

TerHeine_2014_dextromethorphan <- function() {
  description <- paste(
    "Semi-physiological eight-compartment population PK model for",
    "dextromethorphan and its three phase I metabolites (dextrorphan,",
    "3-methoxymorphinan, 3-hydroxymorphinan) in adult breast-cancer",
    "patients receiving chronic oral tamoxifen, used as a dual",
    "CYP2D6 / CYP3A phenotypic probe (single 30 mg oral dose).",
    "Pre-systemic and systemic metabolism are integrated via a hypothetical",
    "hepatic metabolism compartment in rapid (quasi-steady-state) equilibrium",
    "with the dextromethorphan central compartment; the algebraic hepatic",
    "concentration drives parallel CYP2D6 (dextromethorphan -> dextrorphan)",
    "and CYP3A (dextromethorphan -> 3-methoxymorphinan) formation steps.",
    "Subsequent CYP3A-mediated conversion of dextrorphan and",
    "CYP2D6-mediated conversion of 3-methoxymorphinan both feed the",
    "terminal 3-hydroxymorphinan pool, which is eliminated by a single",
    "clearance to other species. All metabolite apparent volumes are",
    "fixed to 419 L (Abduljalil 2009 literature value) for",
    "identifiability. Individual post-hoc CYP2D6 (CL_CYP2D6,1) and",
    "CYP3A (CL_CYP3A,1) clearances serve as the phenotypic probe",
    "covariates used downstream by the companion tamoxifen / endoxifen",
    "model (TerHeine_2014_tamoxifen)."
  )
  reference <- paste(
    "ter Heine R, Binkhorst L, de Graan AJM, de Bruijn P, Beijnen JH,",
    "Mathijssen RHJ, Huitema ADR.",
    "Population pharmacokinetic modelling to assess the impact of CYP2D6",
    "and CYP3A metabolic phenotypes on the pharmacokinetics of tamoxifen",
    "and endoxifen.",
    "Br J Clin Pharmacol. 2014;78(3):572-586.",
    "doi:10.1111/bcp.12388.",
    sep = " "
  )
  vignette <- "TerHeine_2014_tamoxifen"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    # No estimated covariates in the dextromethorphan model itself; the
    # model is a typical-value semi-physiological structure used to
    # *derive* per-subject CYP2D6 and CYP3A activity scores for the
    # downstream tamoxifen model. Covariate column declaration is
    # therefore empty for this model.
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    age_range      = "22 - 71 years",
    age_median     = "53 years",
    weight_range   = "48.5 - 114 kg",
    weight_median  = "72.7 kg",
    height_range   = "1.56 - 1.79 m",
    height_median  = "1.69 m",
    sex_female_pct = 100,
    disease_state  = paste(
      "Hormone-receptor-positive breast cancer on chronic oral tamoxifen",
      "(20 mg QD in 70% of subjects, 40 mg QD in 30%). Patients on",
      "moderate or strong CYP3A inhibitors / inducers or ABCB1 / ABCG2",
      "modulators were excluded."
    ),
    dose_range     = paste(
      "Dextromethorphan: single 30 mg oral dose administered 2 h after",
      "the daily tamoxifen dose. Background chronic tamoxifen at 20 or",
      "40 mg QD."
    ),
    regions        = "The Netherlands (Erasmus MC Cancer Institute, Rotterdam)",
    notes          = paste(
      "40 patients enrolled, 1 dropped due to obstruction of a venous",
      "cannula; n = 39 evaluable. Dextromethorphan plasma samples at",
      "pre-dose and 0.5, 1, 1.5, 2, 4, 6, 10, 22 h post-dose. 293",
      "dextromethorphan, 296 dextrorphan, 139 3-methoxymorphinan and",
      "249 3-hydroxymorphinan concentrations available. Lower limit of",
      "quantitation 0.5 nM for dextromethorphan and metabolites; BLQ",
      "values with mass-spectrometer signal-to-noise > 5 were retained.",
      "Demographics from Table 1; cohort sex is 100% female because the",
      "indication is hormone-receptor-positive breast cancer.",
      "Dutch Trial Registry NTR1751."
    )
  )

  ini({
    # =================================================================
    # Structural parameters -- ter Heine 2014 Table 2 (dextromethorphan
    # and metabolites pharmacokinetic model). All clearance and flow
    # parameters are apparent (oral, /F).
    # =================================================================
    lka <- log(0.213)
    label("First-order absorption rate constant k12 (1/h)")
    # Table 2: k12 = 0.213 1/h, RSE 7.6%
    ltlag <- log(0.369)
    label("Absorption lag time (h)")
    # Table 2: tlag = 0.369 h, RSE 3.7%
    lvc <- log(188)
    label("Apparent central volume of distribution of dextromethorphan (L)")
    # Table 2: Central Vd dextromethorphan = 188 L, RSE 32.4%
    lvp <- log(1660)
    label("Apparent peripheral volume of distribution of dextromethorphan (L)")
    # Table 2: Peripheral Vd dextromethorphan = 1660 L, RSE 16.8%
    lvc_metab <- fixed(log(419))
    label("Apparent volume of distribution shared by all metabolites (L; FIXED)")
    # Methods, "Dextromethorphan and metabolites pharmacokinetic model
    # development": "we therefore assumed a volume of distribution of
    # 419 l for all dextromethorphan metabolites, as previously found
    # for dextrorphan [Abduljalil 2009]."
    lq <- log(415)
    label("Apparent hepatic-flow Q1 connecting dextromethorphan central and the algebraic hepatic compartment (L/h)")
    # Table 2: Q1 = 415 L/h, RSE 22.1%
    lq2 <- log(200)
    label("Apparent peripheral inter-compartmental clearance Q2 (L/h)")
    # Table 2: Q2 = 200 L/h, RSE 37.3%

    # ---- CYP2D6 / CYP3A mediated formation and onward conversion ----
    # The four "CL" terms are labelled in the paper by the catalysing
    # enzyme and the formation step (CL_CYP2D6,1 = dextromethorphan to
    # dextrorphan; CL_CYP2D6,2 = 3-methoxymorphinan to 3-hydroxymorphinan;
    # CL_CYP3A,1 = dextromethorphan to 3-methoxymorphinan; CL_CYP3A,2 =
    # dextrorphan to 3-hydroxymorphinan). The parameter names below
    # preserve the paper's enzyme-then-step numbering so the
    # source-trace is unambiguous.
    lcl_2d6_1 <- log(1560)
    label("CYP2D6-mediated formation clearance of dextrorphan from dextromethorphan (L/h)")
    # Table 2: CL_2D6,1 = 1560 L/h, RSE 27.8%
    lcl_2d6_2 <- log(362)
    label("CYP2D6-mediated formation clearance of 3-hydroxymorphinan from 3-methoxymorphinan (L/h)")
    # Table 2: CL_2D6,2 = 362 L/h, RSE 46.1%
    lcl_3a4_1 <- log(44.7)
    label("CYP3A-mediated formation clearance of 3-methoxymorphinan from dextromethorphan (L/h)")
    # Table 2: CL_3A4,1 = 44.7 L/h, RSE 26%
    lcl_3a4_2 <- log(1840)
    label("CYP3A-mediated formation clearance of 3-hydroxymorphinan from dextrorphan (L/h)")
    # Table 2: CL_3A4,2 = 1840 L/h, RSE 9.1%
    lcl_hm <- log(5730)
    label("Apparent elimination clearance of 3-hydroxymorphinan to other species (L/h)")
    # Table 2: CL_HM = 5730 L/h, RSE 11.4%

    # ---- Central-to-peripheral first-order rate constants for the
    # 3-methoxymorphinan and 3-hydroxymorphinan distribution chains.
    # The paper parameterises these peripheral compartments via rate
    # constants (1/h) rather than (Q, Vp) pairs; the peripheral volumes
    # are therefore not independently identifiable. The compartment
    # numbering 5-8 follows Figure 2 of the source paper (5 =
    # 3-methoxymorphinan central, 6 = 3-hydroxymorphinan central,
    # 7 = 3-methoxymorphinan peripheral, 8 = 3-hydroxymorphinan
    # peripheral).
    lk57 <- log(1.11)
    label("3-methoxymorphinan central -> peripheral rate constant k57 (1/h)")
    # Table 2: k57 = 1.11 1/h, RSE 19.5%
    lk75 <- log(0.0815)
    label("3-methoxymorphinan peripheral -> central rate constant k75 (1/h)")
    # Table 2: k75 = 0.0815 1/h, RSE 28.7%
    lk68 <- log(13.4)
    label("3-hydroxymorphinan central -> peripheral rate constant k68 (1/h)")
    # Table 2: k68 = 13.4 1/h, RSE 12.4%
    lk86 <- log(0.0637)
    label("3-hydroxymorphinan peripheral -> central rate constant k86 (1/h)")
    # Table 2: k86 = 0.0637 1/h, RSE 18.5%

    # =================================================================
    # Inter-individual variability -- Table 2. IIV is log-normal in the
    # paper's notation; CV% maps to log-variance via
    # omega^2 = log(1 + (CV/100)^2). The CYP2D6 clearances
    # CL_2D6,1 and CL_2D6,2 were >98% correlated and a SINGLE shared
    # variability term was estimated for both (paper Results,
    # "Dextromethorphan and metabolites pharmacokinetic model
    # development"). All other listed IIV terms are independent.
    # =================================================================
    etalka ~ log(1 + 0.31^2)
    # Table 2: IIV ka = 31% CV (shrinkage 4%), RSE 39.7%
    etalcl_2d6 ~ log(1 + 1.66^2)
    # Table 2: IIV CL_2D6 = 166% CV (shrinkage -0.9%), RSE 29.5%.
    # Applied to both lcl_2d6_1 and lcl_2d6_2 (shared eta).
    etalcl_3a4_1 ~ log(1 + 0.898^2)
    # Table 2: IIV CL_3A4,1 = 89.8% CV (shrinkage 19.8%), RSE 49.3%
    etalcl_3a4_2 ~ log(1 + 0.530^2)
    # Table 2: IIV CL_3A4,2 = 53% CV (shrinkage 2.4%), RSE 27.8%
    etalcl_hm ~ log(1 + 0.757^2)
    # Table 2: IIV CL_HM = 75.7% CV (shrinkage 6.2%), RSE 24.6%

    # =================================================================
    # Residual error -- proportional only (paper Methods, "Pharmacokinetic
    # analysis": "residual error was described with a proportional error
    # model"). One residual SD per measured species. The paper reports
    # the residuals as %CV, which equals the proportional-error SD on the
    # linear scale to first order.
    # =================================================================
    propSd <- 0.429
    label("Proportional residual SD for dextromethorphan plasma concentration (fraction)")
    # Table 2: Residual error dextromethorphan = 42.9% CV (shrinkage 5.2%), RSE 15%
    propSd_dxor <- 0.400
    label("Proportional residual SD for dextrorphan plasma concentration (fraction)")
    # Table 2: Residual error dextrorphan = 40% CV (shrinkage 7.2%), RSE 10.6%
    propSd_3mm <- 0.297
    label("Proportional residual SD for 3-methoxymorphinan plasma concentration (fraction)")
    # Table 2: Residual error 3-methoxymorphinan = 29.7% CV (shrinkage 9.2%), RSE 23.2%
    propSd_3hm <- 0.307
    label("Proportional residual SD for 3-hydroxymorphinan plasma concentration (fraction)")
    # Table 2: Residual error 3-hydroxymorphinan = 30.7% CV (shrinkage 7%), RSE 12%
  })

  model({
    # ---- Molecular masses (g/mol) used to convert the mg-scale input
    # dose to nmol at the depot. The paper states (Methods,
    # "Pharmacokinetic analysis") that all parent and metabolite
    # concentrations were converted to molar equivalents; tracking
    # internal compartment amounts in nmol therefore keeps the parent-
    # to-metabolite 1:1 molar conservation exact across the four
    # species (whose MWs differ by up to ~10%). Concentrations are
    # then directly in nmol/L = nM, matching the LOQ (0.5 nM) and the
    # paper's reporting convention.
    mdex <- 271.40  # dextromethorphan (C18H25NO)

    # Individual parameters. The single CL_2D6 IIV term is applied to
    # both lcl_2d6_1 and lcl_2d6_2 (paper Results).
    ka       <- exp(lka + etalka)
    vc       <- exp(lvc)
    vp       <- exp(lvp)
    vc_metab <- exp(lvc_metab)
    q        <- exp(lq)
    q2       <- exp(lq2)
    cl_2d6_1 <- exp(lcl_2d6_1 + etalcl_2d6)
    cl_2d6_2 <- exp(lcl_2d6_2 + etalcl_2d6)
    cl_3a4_1 <- exp(lcl_3a4_1 + etalcl_3a4_1)
    cl_3a4_2 <- exp(lcl_3a4_2 + etalcl_3a4_2)
    cl_hm    <- exp(lcl_hm + etalcl_hm)
    k57      <- exp(lk57)
    k75      <- exp(lk75)
    k68      <- exp(lk68)
    k86      <- exp(lk86)

    # Concentrations from compartment amounts. Depot mg-to-nmol scaling
    # is folded into the absorption rate so internal amounts are nmol.
    c_central <- central / vc
    c_perip   <- peripheral1 / vp
    c_dxor    <- central_dxor / vc_metab
    c_3mm     <- central_3mm / vc_metab
    c_3hm     <- central_3hm / vc_metab

    # ---- Hypothetical hepatic metabolism compartment in rapid
    # (quasi-steady-state) equilibrium with the dextromethorphan
    # central compartment (Appendix 1):
    #   C_L = (k12 * A1 + Q1 * C2) /
    #         (Q1 + CL_CYP2D6,1 + CL_CYP3A,1)
    # The absorption-rate factor ka * depot is in mg/h; the molar
    # unit conversion (1e6 / Mdex nmol per mg) brings it onto the
    # same nmol/h footing as the systemic concentration terms.
    ratein <- ka * depot * 1e6 / mdex
    c_liver <- (ratein + q * c_central) /
               (q + cl_2d6_1 + cl_3a4_1)

    # ---- ODE system (Appendix 1). Compartment numbering matches the
    # paper's Figure 2: 1 = depot, 2 = central dextromethorphan,
    # 3 = central dextrorphan, 4 = peripheral dextromethorphan,
    # 5 = central 3-methoxymorphinan, 6 = central 3-hydroxymorphinan,
    # 7 = peripheral 3-methoxymorphinan, 8 = peripheral
    # 3-hydroxymorphinan. Depot equation is in mg (-ka * depot, mg/h);
    # all other compartments accumulate in nmol via the molar-unit
    # conversion in `ratein`.
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  q * c_liver - q * c_central -
                               q2 * c_central + q2 * c_perip
    d/dt(peripheral1)      <-  q2 * c_central - q2 * c_perip
    d/dt(central_dxor)     <-  cl_2d6_1 * c_liver - cl_3a4_2 * c_dxor
    d/dt(central_3mm)      <-  cl_3a4_1 * c_liver + k75 * peripheral1_3mm -
                               cl_2d6_2 * c_3mm - k57 * central_3mm
    d/dt(central_3hm)      <-  cl_3a4_2 * c_dxor + cl_2d6_2 * c_3mm +
                               k86 * peripheral1_3hm -
                               cl_hm * c_3hm - k68 * central_3hm
    d/dt(peripheral1_3mm)   <-  k57 * central_3mm - k75 * peripheral1_3mm
    d/dt(peripheral1_3hm)   <-  k68 * central_3hm - k86 * peripheral1_3hm

    # ---- Absorption lag (standard NONMEM ALAG mechanism on depot).
    lag(depot) <- exp(ltlag)

    # ---- Observations. Internal compartment amounts are in nmol;
    # dividing by the apparent volume in L gives nmol/L = nM.
    Cc       <- c_central
    Cc_dxor  <- c_dxor
    Cc_3mm   <- c_3mm
    Cc_3hm   <- c_3hm

    Cc       ~ prop(propSd)
    Cc_dxor  ~ prop(propSd_dxor)
    Cc_3mm   ~ prop(propSd_3mm)
    Cc_3hm   ~ prop(propSd_3hm)
  })
}

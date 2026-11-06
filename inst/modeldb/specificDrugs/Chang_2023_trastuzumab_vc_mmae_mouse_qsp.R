Chang_2023_trastuzumab_vc_mmae_mouse_qsp <- function() {
  description <- "QSP (preclinical, mouse). Semi-mechanistic multi-analyte PK/PD model for the HER2-targeting antibody-drug conjugate trastuzumab-vc-MMAE (T-vc-MMAE, model payload MMAE) in NCI-N87 gastric-carcinoma xenograft-bearing nude mice, fit simultaneously to intravenous, subcutaneous and intratumoral administration. Carries four plasma analytes in a shared two-compartment disposition (intact conjugated antibody, naked antibody formed by deconjugation, conjugated MMAE, and unconjugated MMAE) plus a mechanistic Krogh-cylinder tumor-disposition module (vascular extravasation + surface diffusion, HER2 antigen binding/internalization, tumor-microenvironment-amplified deconjugation, tubulin binding of intracellular MMAE). Intracellular target engagement (%tubulin occupied) drives a nonlinear tumor-growth-inhibition model with an exponential-to-linear growth law and a four-compartment cell-death transit chain. Amounts in nmol, concentrations in nM, time in h."
  reference <- paste(
    "Chang HP, Le HK, Shah DK. Pharmacokinetics and Pharmacodynamics of",
    "Antibody-Drug Conjugates Administered via Subcutaneous and Intratumoral",
    "Routes. Pharmaceutics. 2023;15(4):1132. doi:10.3390/pharmaceutics15041132.",
    "Model structure: Figure 2 and Equations (1)-(24); parameters: Table 1.",
    sep = " "
  )
  vignette <- "Chang_2023_trastuzumab_vc_mmae"
  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # What each ODE state holds, in what amount units, in which matrix. Amounts in
  # nmol throughout; the tumor pools track the same analytes in the tumor
  # extracellular/intracellular space. verified = TRUE: matched to Figure 2 and
  # Equations (1)-(24) of the source.
  compartmentData <- list(
    depot_adc = list(
      analyte = "conjugated antibody (ADC)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_acpl = list(
      analyte = "conjugated MMAE payload",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    adc_central = list(analyte = "conjugated antibody (ADC)", units = "nmol", specimen = "plasma", verified = TRUE),
    adc_peripheral = list(analyte = "conjugated antibody (ADC)", units = "nmol", specimen = "tissue", verified = TRUE),
    mab_central = list(analyte = "naked antibody", units = "nmol", specimen = "plasma", verified = TRUE),
    mab_peripheral = list(analyte = "naked antibody", units = "nmol", specimen = "tissue", verified = TRUE),
    acpl_central = list(analyte = "conjugated MMAE payload", units = "nmol", specimen = "plasma", verified = TRUE),
    acpl_peripheral = list(analyte = "conjugated MMAE payload", units = "nmol", specimen = "tissue", verified = TRUE),
    mmae_central = list(analyte = "unconjugated MMAE", units = "nmol", specimen = "plasma", verified = TRUE),
    mmae_peripheral = list(analyte = "unconjugated MMAE", units = "nmol", specimen = "tissue", verified = TRUE),
    adc_free_tumor = list(
      analyte = "conjugated antibody (ADC), antigen-unbound",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_tumor = list(
      analyte = "conjugated antibody (ADC), HER2-bound",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    mab_free_tumor = list(
      analyte = "naked antibody, antigen-unbound",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    mab_bound_tumor = list(analyte = "naked antibody, HER2-bound", units = "nmol", specimen = "tumor", verified = TRUE),
    acpl_free_tumor = list(
      analyte = "conjugated MMAE payload, antigen-unbound",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    acpl_bound_tumor = list(
      analyte = "conjugated MMAE payload, HER2-bound",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    mmae_ext_tumor = list(analyte = "unconjugated MMAE", units = "nmol", specimen = "tumor", verified = TRUE),
    mmae_free_tumor = list(
      analyte = "unconjugated MMAE, intracellular free",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    mmae_tub_tumor = list(analyte = "MMAE bound to tubulin", units = "nmol", specimen = "tumor", verified = TRUE),
    tumor1 = list(analyte = "viable (proliferating) tumor volume", units = "mm^3", specimen = "tumor", verified = TRUE),
    tumor2 = list(analyte = "dying tumor volume (transit 1)", units = "mm^3", specimen = "tumor", verified = TRUE),
    tumor3 = list(analyte = "dying tumor volume (transit 2)", units = "mm^3", specimen = "tumor", verified = TRUE),
    tumor4 = list(analyte = "dying tumor volume (transit 3)", units = "mm^3", specimen = "tumor", verified = TRUE)
  )

  # Covariate: per-dose-record intratumoral-route indicator. 1 = intratumoral
  # (IT) dose, 0 = intravenous or subcutaneous. Gates the FIT retained-fraction
  # bioavailability on IT doses placed in the central compartment (the escaped
  # ~25% of an IT dose that immediately reaches plasma). See
  # inst/references/covariate-columns.md (ROUTE_IT).
  covariateData <- list(
    ROUTE_IT = list(
      description = "Intratumoral-vs-IV administration route indicator (1 = IT, 0 = IV/SC)",
      units = "(binary)",
      type = "binary",
      source_name = "administration route (Section 2.3 / 2.5)"
    ),
    TUM_VOL = list(
      description = "Measured initial (baseline) tumor volume, sets the viable tumor state tumor1(0)",
      units = "mm^3",
      type = "continuous",
      source_name = "initial tumor volume TV0 (Table 1, Experimental; Table S1 per-arm means 243-296 mm3)"
    )
  )

  covariatesDataExcluded <- list()

  population <- list(
    species = "mouse (male athymic nude, NCI-N87 gastric-carcinoma subcutaneous xenograft)",
    n_subjects = 49L,
    n_studies = 1L,
    disease_state = "HER2-overexpressing NCI-N87 gastric-carcinoma flank xenograft",
    age_range = "7 weeks at cell implantation (5 weeks on arrival + 2 weeks acclimation)",
    weight_range = "24-36 g body weight",
    sex_female_pct = 0,
    dose_range = "10 mg/kg single dose for PK (IV n=15, SC n=15, IT n=19); efficacy at high/mid/low: IV & IT 10/3/1 mg/kg, SC 20/6/2 mg/kg single dose",
    regions = "Preclinical (University at Buffalo)",
    notes = paste(
      "Average drug-to-antibody ratio (DAR) ~3.5 (HIC, Figure S1).",
      "Plasma mAb PK parameters (CLmAb, CLDmAb, V1mAb, V2mAb) were estimated in a",
      "prior step by fitting trastuzumab plasma data in nude mice (ref [37]) and",
      "then fixed; unconjugated-MMAE plasma PK (CLPL, CLDPL, V1PL, V2PL) were",
      "estimated from free-MMAE mouse data (ref [47]) and fixed. Tumor-disposition",
      "constants are literature-derived (Table 1 refs [40-49]). Trastuzumab does",
      "not bind mouse HER2, so there is no systemic antigen sink. Tumor volumes",
      "<50 mm3 were treated as interval-censored during fitting. PK fit in ADAPT 5",
      "(additive + proportional error); PD fit in Monolix (SAEM). Initial tumor",
      "volume ~150-200 mm3 (Table S1 per-arm means 243-296 mm3)."
    )
  )

  # Every d/dt() state below is a mechanistic, model-specific amount pool (nmol),
  # not a canonical library compartment; whitelist them so checkModelConventions()
  # does not flag them. Order fixes the compartment slots used by the event table.
  paper_specific_compartments <- c(
    "depot_adc",
    "depot_acpl",
    "adc_central",
    "adc_peripheral",
    "mab_central",
    "mab_peripheral",
    "acpl_central",
    "acpl_peripheral",
    "mmae_central",
    "mmae_peripheral",
    "adc_free_tumor",
    "adc_bound_tumor",
    "mab_free_tumor",
    "mab_bound_tumor",
    "acpl_free_tumor",
    "acpl_bound_tumor",
    "mmae_ext_tumor",
    "mmae_free_tumor",
    "mmae_tub_tumor",
    "tumor1",
    "tumor2",
    "tumor3",
    "tumor4"
  )

  ini({
    # --------------------------------------------------------------------- #
    # Plasma / systemic PK (Table 1).  Volumes are tabulated in mL and      #
    # clearances in mL/h; they are converted to L / L/h in model() so that  #
    # amount(nmol)/volume(L) yields nM to match the nM tumor module.        #
    # --------------------------------------------------------------------- #
    # Conjugated antibody (ADC) / naked mAb disposition. Estimated by fitting
    # trastuzumab mouse plasma data (ref [37]) then fixed for the ADC fit.
    lclmab <- fixed(log(0.006)); label("Plasma clearance of ADC/mAb CLmAb (mL/h)") # Table 1: 0.006 (10.8% RSE), fixed from ref [37] fit
    lcldmab <- fixed(log(0.040)); label("Distribution clearance of ADC/mAb CLDmAb (mL/h)") # Table 1: 0.040 (28.9%)
    lv1mab <- fixed(log(1.41)); label("Central volume of ADC/mAb V1mAb (mL)") # Table 1: 1.41 (6.48%)
    lv2mab <- fixed(log(0.861)); label("Peripheral volume of ADC/mAb V2mAb (mL)") # Table 1: 0.861 (13.3%)

    # Unconjugated MMAE disposition. Estimated from free-MMAE mouse data
    # (ref [47]) then fixed; CLPL via NCA.
    lclpl <- fixed(log(15.9)); label("Plasma clearance of unconjugated MMAE CLPL (mL/h)") # Table 1: 15.9, fixed from ref [47]
    lcldpl <- fixed(log(0.811)); label("Distribution clearance of unconjugated MMAE CLDPL (mL/h)") # Table 1: 0.811 (20.0%)
    lv1pl <- fixed(log(2.25)); label("Central volume of unconjugated MMAE V1PL (mL)") # Table 1: 2.25 (9.68%)
    lv2pl <- fixed(log(5.60)); label("Peripheral volume of unconjugated MMAE V2PL (mL)") # Table 1: 5.60 (25.7%)

    # Deconjugation rate constants (estimated in the ADC fit).
    lkdecadc <- log(0.00344); label("Deconjugation rate forming naked mAb kdec,ADC (1/h)") # Table 1: 0.00344 (9.65%)
    lkdecpl <- log(0.00905); label("Deconjugation rate releasing MMAE kdec,PL (1/h)") # Table 1: 0.00905 (6.48%)

    # --------------------------------------------------------------------- #
    # Tumor-disposition module (literature-derived; Table 1 refs [40-49]).  #
    # Permeabilities in mm/h, diffusivities in mm^2/h, radii in mm.         #
    # --------------------------------------------------------------------- #
    lp_ab <- fixed(log(0.01)); label("Tumor vascular permeability of ADC PADC (mm/h)") # Table 1: 0.01 [40]
    lp_mmae <- fixed(log(0.0875)); label("Tumor vascular permeability of MMAE PMMAE (mm/h)") # Table 1: 0.0875 [41]
    ld_ab <- fixed(log(0.00054)); label("Diffusion rate of ADC DADC (mm^2/h)") # Table 1: 0.00054 [42]
    ld_mmae <- fixed(log(1.04)); label("Diffusion rate of MMAE DMMAE (mm^2/h)") # Table 1: 1.04 [41]
    lrcap <- fixed(log(0.008)); label("Tumor blood-capillary radius Rcap (mm)") # Table 1: 0.008 [43]
    lrkrogh <- fixed(log(0.075)); label("Average inter-capillary (Krogh) distance Rkrough (mm)") # Table 1: 0.075 [43]
    leps_ab <- fixed(log(0.24)); label("Tumor void fraction for ADC epsilon_ADC (unitless)") # Table 1: 0.24 [41]
    leps_mmae <- fixed(log(0.44)); label("Tumor void fraction for MMAE epsilon_MMAE (unitless)") # Table 1: 0.44 [41]

    lkon_ag <- fixed(log(1.25)); label("ADC-HER2 association rate kon,Antigen (1/nM/h)") # Table 1: 1.25 [44]
    lkoff_ag <- fixed(log(2.26)); label("ADC-HER2 dissociation rate koff,Antigen (1/h)") # Table 1: 2.26 [44]
    lkint_cplx <- log(0.112); label("Internalization of ADC-antigen complex kint,ADC (1/h)") # Table 1: 0.112 (19.6%), estimated
    lkint_mab <- fixed(log(0.027)); label("Internalization of mAb-antigen complex kint,mAb (1/h)") # Table 1: 0.027 [45]
    lantigen <- fixed(log(1799)); label("Total HER2 antigen concentration Antigen (nM)") # Table 1: 1799 [44]

    lkon_tub <- fixed(log(0.00187)); label("MMAE-tubulin association rate kon,Tubulin (1/nM/h)") # Table 1: 0.00187 [46,47]
    lkoff_tub <- fixed(log(0.545)); label("MMAE-tubulin dissociation rate koff,Tubulin (1/h)") # Table 1: 0.545 [38,46]
    ltubulin <- fixed(log(500)); label("Total tubulin concentration Tubulin (nM)") # Table 1: 500 [48]
    lkin <- fixed(log(0.075)); label("MMAE nonspecific cellular uptake kin (1/h)") # Table 1: 0.075 [49]
    lkout <- log(0.0116); label("MMAE cellular efflux kout (1/h)") # Table 1: 0.0116 (30.9%), estimated
    lfdec <- log(30.8); label("Fold-increase of deconjugation in tumor fdec (unitless)") # Table 1: 30.8 (20.7%), estimated

    # --------------------------------------------------------------------- #
    # Route-specific parameters (estimated).                                #
    # --------------------------------------------------------------------- #
    lfit <- log(0.75); label("Fraction of IT dose retained in tumor FIT (fraction)") # Table 1: 75.0% (4.61%)
    lfsc <- log(0.476); label("SC bioavailability F (fraction)") # Table 1: 47.6% (0.836%)
    lka <- log(0.0498); label("SC first-order absorption rate ka (1/h)") # Table 1: 0.0498 (1.96%)
    ldar0 <- fixed(log(3.5)); label("Characterized drug-to-antibody ratio DAR0 (unitless)") # Results 3.1 / Figure S1: average DAR ~3.5

    # --------------------------------------------------------------------- #
    # Tumor growth / kill (PD; Table 1).  Growth estimated from untreated    #
    # controls (step 4); kill from all treated arms (step 5).                #
    # --------------------------------------------------------------------- #
    lkgl <- log(3.08); label("Zero-order (linear-phase) tumor growth kgl (mm^3/h)") # Table 1: 3.08 (24.7%)
    lkge <- log(0.0018); label("First-order (exponential-phase) tumor growth kge (1/h)") # Table 1: 0.0018 (3.21%)
    lkmax <- log(0.00673); label("Maximum cell-killing rate Kmax (1/h)") # Table 1: 0.00673 (6.18%)
    lec50 <- log(17.9); label("Percent tubulin engagement giving 50% kill EC50 (%)") # Table 1: 17.9 (1.17%)
    ltau <- log(19.3); label("Mean transit time of the cell-death chain tau (h)") # Table 1: 19.3 (21.6%)

    # Between-mouse variability (Monolix SAEM population variability, Table 1,
    # reported as CV%); omega^2 = log(1 + CV^2) for the log-normal random effect.
    etalkgl ~ log(1 + 0.538^2); label("IIV on kgl (Table 1 Var_kgl 53.8%)") # 0.254
    etalkmax ~ log(1 + 0.432^2); label("IIV on Kmax (Table 1 Var_Kmax 43.2%)") # 0.171
    etaltau ~ log(1 + 1.31^2); label("IIV on tau (Table 1 Var_tau 131%)") # 0.999

    # Residual error. Per-analyte residual magnitudes are NOT reported in the
    # paper (PK: ADAPT additive + proportional variance; PD: Monolix); a nominal
    # proportional error is supplied for simulation only (see vignette Errata).
    propSd <- 0.15; label("Proportional residual SD on plasma total antibody (simulation default; not from paper)")
  })

  model({
    # ---- back-transform parameters ----
    # Volumes mL -> L and clearances mL/h -> L/h so amount(nmol)/volume(L) = nM.
    vc_mab <- exp(lv1mab) / 1000
    vp_mab <- exp(lv2mab) / 1000
    cl_mab <- exp(lclmab) / 1000
    cld_mab <- exp(lcldmab) / 1000
    vc_pl <- exp(lv1pl) / 1000
    vp_pl <- exp(lv2pl) / 1000
    cl_pl <- exp(lclpl) / 1000
    cld_pl <- exp(lcldpl) / 1000
    kdec_adc <- exp(lkdecadc)
    kdec_pl <- exp(lkdecpl)

    p_adc <- exp(lp_ab)
    p_mmae <- exp(lp_mmae)
    d_adc <- exp(ld_ab)
    d_mmae <- exp(ld_mmae)
    rcap <- exp(lrcap)
    rkrogh <- exp(lrkrogh)
    eps_adc <- exp(leps_ab)
    eps_mmae <- exp(leps_mmae)
    kon_ag <- exp(lkon_ag)
    koff_ag <- exp(lkoff_ag)
    kint_cplx <- exp(lkint_cplx)
    kint_mab <- exp(lkint_mab)
    antigen <- exp(lantigen)
    kon_tub <- exp(lkon_tub)
    koff_tub <- exp(lkoff_tub)
    tubulin <- exp(ltubulin)
    kin <- exp(lkin)
    kout <- exp(lkout)
    fdec <- exp(lfdec)

    fit <- exp(lfit)
    fsc <- exp(lfsc)
    ka <- exp(lka)
    dar0 <- exp(ldar0)

    kgl <- exp(lkgl + etalkgl)
    kge <- exp(lkge)
    kmax <- exp(lkmax + etalkmax)
    ec50 <- exp(lec50)
    tau <- exp(ltau + etaltau)

    # ---- geometry: total tumor volume and derived radius ----
    # Vtumor tracks the total (viable + dying) tumor volume; mm^3 -> L for conc.
    tv_mm3 <- tumor1 + tumor2 + tumor3 + tumor4
    vtumor <- tv_mm3 * 1e-6
    rtumor <- (3 * tv_mm3 / (4 * 3.141592653589793))^(1 / 3) # mm, from spherical volume

    # ---- vascular (2*P*Rcap/Rkrogh^2) + surface (6*D/Rtumor^2) exchange, 1/h ----
    kx_adc <- 2 * p_adc * rcap / rkrogh^2 + 6 * d_adc / rtumor^2
    kx_mmae <- 2 * p_mmae * rcap / rkrogh^2 + 6 * d_mmae / rtumor^2

    # ---- plasma concentrations (nM) ----
    c_adc <- adc_central / vc_mab
    c_mab <- mab_central / vc_mab
    c_acpl <- acpl_central / vc_mab
    c_mmae <- mmae_central / vc_pl

    # ---- time-varying DAR in plasma (Eq 9) and in tumor (Eq 19) ----
    dar_plasma <- acpl_central / (adc_central + mab_central + 1e-12)
    dar_tumor <- (acpl_free_tumor + acpl_bound_tumor) /
      (adc_free_tumor + adc_bound_tumor + mab_free_tumor + mab_bound_tumor + 1e-12)

    # ---- Krogh-cylinder extravasation + surface-diffusion fluxes (nmol/h) ----
    # flux = kx * (C_plasma - C_tumor_extracellular / epsilon) * Vtumor
    ex_adc <- kx_adc * (c_adc - (adc_free_tumor / vtumor) / eps_adc) * vtumor
    ex_mab <- kx_adc * (c_mab - (mab_free_tumor / vtumor) / eps_adc) * vtumor
    ex_acpl <- kx_adc * (c_acpl - (acpl_free_tumor / vtumor) / eps_adc) * vtumor
    ex_mmae <- kx_mmae * (c_mmae - (mmae_ext_tumor / vtumor) / eps_mmae) * vtumor

    # ---- HER2 antigen binding in tumor extracellular space (nmol/h) ----
    # Each species competes for antigen reduced only by its own bound form, as
    # written in Eqs (10)/(14); the printed Eqs (12)/(13) antigen terms use the
    # unbound/cross-species forms, read here as transcription slips (see Errata).
    bind_adc <- kon_ag * (adc_free_tumor / vtumor / eps_adc) *
      (antigen - adc_bound_tumor / vtumor) * vtumor - koff_ag * adc_bound_tumor
    bind_mab <- kon_ag * (mab_free_tumor / vtumor / eps_adc) *
      (antigen - mab_bound_tumor / vtumor) * vtumor - koff_ag * mab_bound_tumor
    bind_acpl <- kon_ag * (acpl_free_tumor / vtumor / eps_adc) *
      (antigen - acpl_bound_tumor / vtumor) * vtumor - koff_ag * acpl_bound_tumor

    # ---- intracellular MMAE-tubulin binding (nmol/h) ----
    bind_tub <- kon_tub * (mmae_free_tumor / vtumor) *
      (tubulin - mmae_tub_tumor / vtumor) * vtumor - koff_tub * mmae_tub_tumor

    # ---- SC absorption depots (Eqs 1/5, +ka*Xa term) ----
    d/dt(depot_adc) <- -ka * depot_adc
    d/dt(depot_acpl) <- -ka * depot_acpl

    # ---- plasma / peripheral disposition (Eqs 1-8) ----
    d/dt(adc_central) <- -cld_mab * c_adc + cld_mab * (adc_peripheral / vp_mab) -
      cl_mab * c_adc - kdec_adc * adc_central - ex_adc + ka * depot_adc
    d/dt(adc_peripheral) <- cld_mab * c_adc - cld_mab * (adc_peripheral / vp_mab) -
      kdec_adc * adc_peripheral

    d/dt(mab_central) <- -cld_mab * c_mab + cld_mab * (mab_peripheral / vp_mab) -
      cl_mab * c_mab + kdec_adc * adc_central - ex_mab
    d/dt(mab_peripheral) <- cld_mab * c_mab - cld_mab * (mab_peripheral / vp_mab) +
      kdec_adc * adc_peripheral

    # Conjugated-payload tracker: loses payload with the ADC it rides on
    # (CLmAb*c_adc*DAR), deconjugates at kdec,PL, extravasates, absorbs from SC.
    d/dt(acpl_central) <- -cld_mab * c_acpl + cld_mab * (acpl_peripheral / vp_mab) -
      cl_mab * c_adc * dar_plasma - kdec_pl * acpl_central - ex_acpl + ka * depot_acpl
    d/dt(acpl_peripheral) <- cld_mab * c_acpl - cld_mab * (acpl_peripheral / vp_mab) -
      kdec_pl * acpl_peripheral

    # Unconjugated MMAE plasma: gains payload from ADC catabolism (CLmAb*DAR) and
    # from conjugated-payload deconjugation (kdec,PL).
    d/dt(mmae_central) <- -cld_pl * c_mmae + cld_pl * (mmae_peripheral / vp_pl) -
      cl_pl * c_mmae + cl_mab * c_adc * dar_plasma + kdec_pl * acpl_central - ex_mmae
    d/dt(mmae_peripheral) <- cld_pl * c_mmae - cld_pl * (mmae_peripheral / vp_pl) +
      kdec_pl * acpl_peripheral

    # ---- tumor extracellular ADC / naked mAb free and antigen-bound (Eqs 10-15) ----
    d/dt(adc_free_tumor) <- ex_adc - bind_adc - fdec * kdec_adc * adc_free_tumor
    d/dt(adc_bound_tumor) <- bind_adc - kint_cplx * adc_bound_tumor -
      fdec * kdec_adc * adc_bound_tumor
    d/dt(mab_free_tumor) <- ex_mab - bind_mab + fdec * kdec_adc * adc_free_tumor
    d/dt(mab_bound_tumor) <- bind_mab - kint_mab * mab_bound_tumor +
      fdec * kdec_adc * adc_bound_tumor
    d/dt(acpl_free_tumor) <- ex_acpl - bind_acpl - fdec * kdec_pl * acpl_free_tumor
    d/dt(acpl_bound_tumor) <- bind_acpl - kint_cplx * acpl_bound_tumor -
      fdec * kdec_pl * acpl_bound_tumor

    # ---- tumor MMAE: extracellular, intracellular free, tubulin-bound (Eqs 16-18) ----
    # Extracellular free MMAE gains the payload released by tumor deconjugation of
    # conjugated payload (both free and antigen-bound); sign per mass balance.
    d/dt(mmae_ext_tumor) <- ex_mmae - kin * mmae_ext_tumor + kout * mmae_free_tumor +
      fdec * kdec_pl * (acpl_free_tumor + acpl_bound_tumor)
    d/dt(mmae_free_tumor) <- kin * mmae_ext_tumor - kout * mmae_free_tumor -
      bind_tub + kint_cplx * acpl_bound_tumor * dar_tumor
    d/dt(mmae_tub_tumor) <- bind_tub

    # ---- intracellular target engagement (Eq 20): % tubulin occupied ----
    pct_te <- (mmae_tub_tumor / vtumor) / tubulin * 100

    # ---- tumor growth-inhibition model (Eqs 21-24) ----
    # Exponential-to-linear Koch growth on the viable pool (tumor1), scaled by the
    # viable fraction tumor1/Vtumor; the printed Eq (21) denominator kge is read as
    # kgl so the linear-phase limit is kgl mm^3/h (see Errata). Kill routes viable
    # cells through a three-stage death transit chain (tumor2-4) at rate 1/tau.
    kill <- kmax * pct_te / (ec50 + pct_te)
    growth <- 2 * kge * kgl * tumor1 * (tumor1 / tv_mm3) / (kgl + 2 * kge * tumor1)
    d/dt(tumor1) <- growth - kill * tumor1
    d/dt(tumor2) <- kill * tumor1 - tumor2 / tau
    d/dt(tumor3) <- (tumor2 - tumor3) / tau
    d/dt(tumor4) <- (tumor3 - tumor4) / tau

    # ---- route-specific bioavailability ----
    # SC: first-order absorption from the depots with bioavailability F.
    # IV (ROUTE_IT = 0): full dose into the central compartments (F = 1).
    # IT (ROUTE_IT = 1): FIT retained in tumor free pools; (1 - FIT) escapes to
    #   plasma. The conjugated-payload amount is DAR0x the antibody amount and is
    #   supplied as such in the event table, so f() carries no DAR0 factor here.
    f(depot_adc) <- fsc
    f(depot_acpl) <- fsc
    f(adc_central) <- 1 - ROUTE_IT * fit
    f(acpl_central) <- 1 - ROUTE_IT * fit
    f(adc_free_tumor) <- fit
    f(acpl_free_tumor) <- fit

    # ---- initial tumor volume: viable pool starts at the measured baseline ----
    tumor1(0) <- TUM_VOL

    # ---- observations (nM unless noted); plasma total antibody is the primary ----
    Cc <- (adc_central + mab_central) / vc_mab # plasma total antibody (nM)
    plasma_total_mmae <- acpl_central / vc_mab + mmae_central / vc_pl
    plasma_free_mmae <- mmae_central / vc_pl
    tumor_total_antibody <- (adc_free_tumor + adc_bound_tumor +
      mab_free_tumor + mab_bound_tumor) / vtumor
    tumor_free_mmae <- (mmae_ext_tumor + mmae_free_tumor + mmae_tub_tumor) / vtumor
    tumor_total_mmae <- (acpl_free_tumor + acpl_bound_tumor) / vtumor + tumor_free_mmae
    tumor_volume <- tv_mm3 # total tumor volume (mm^3)
    target_engagement <- pct_te # % tubulin occupied by MMAE

    Cc ~ prop(propSd)
  })
}

Mallick_2020_pyrethroids_pbpk <- function() {
  description <- "PBPK (whole-body, life-stage, plasma plus six tissue compartments; recoded from acslX to R by the authors, deposited as Supplementary 4). Generic life-stage human PBPK model for the pyrethroid insecticide class of Mallick et al. 2020, parameterised here for the reference compound deltamethrin (DLM) in a 25-year-old adult male. A single generic structure is applied to eight pyrethroids (deltamethrin, cis- and trans-permethrin, esfenvalerate, cyphenothrin, cyhalothrin, cyfluthrin, bifenthrin); only the hepatic intrinsic clearance clint_liver and molecular weight are compound-specific, so a different pyrethroid or age is simulated by overriding clint_liver and the age-specific physiology rather than by a new structure. Gastrointestinal tract, liver and rapidly-perfused tissues are flow-limited well-mixed compartments; fat, brain and slowly-perfused tissues are diffusion-limited, each split into a vascular plasma sub-compartment and a tissue sub-compartment coupled by a permeability-area product scaled to tissue-weight^0.75. Oral dosing enters a gut-lumen depot absorbed first order; a lymphatic fraction bypasses hepatic first pass and enters plasma directly, the remainder enters the GI portal path to the liver. Hepatic elimination is first-order restrictive clearance: the metabolic rate is clint_liver times liver weight times the free liver concentration divided by an empirical free-concentration adjustment factor kmf. Deterministic typical-value model (no IIV, no residual error): the authors built it by IVIVE from expressed-enzyme in vitro clearances and enzyme ontogeny and evaluated internal target-tissue (brain) exposure across ages rather than fitting individual data. The published inhalation, dermal and drinking-water routes are omitted here because every published simulation is a single daily oral dose. Concentrations are in molar units (umol/L); the paper's molecular weights, needed only for the ng/mL display conversion, are not on disk (see vignette Errata)."
  reference <- paste(
    "Mallick P, Moreau M, Song G, Efremenko AY, Pendse SN, Creek MR, Osimitz TG,",
    "Hines RN, Hinderliter P, Clewell HJ, Lake BG, Yoon M. (2020).",
    "Development and Application of a Life-Stage Physiologically Based",
    "Pharmacokinetic (PBPK) Model to the Assessment of Internal Dose of",
    "Pyrethroids in Humans.",
    "Toxicological Sciences 173(1):86-99.",
    "doi:10.1093/toxsci/kfz211. PMCID PMC6944222.",
    "Structure and parameters: Supplementary 4 (R model code) and Supplementary 2,",
    "Table 1S (PBPK parameters); age-specific hepatic clearances from Supplementary 3.",
    sep = " "
  )
  vignette <- "Mallick_2020_pyrethroids"

  units <- list(
    time = "h",
    dosing = "umol",
    concentration = "umol/L",
    amount = "umol",
    weight = "kg"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "All tissue volumes are computed as a fraction of body weight (fv_* * WT).",
        "The reference value is the 25-year-old adult male, 81.73748 kg, from",
        "Supplementary 2, Table 1S. To simulate another life stage, set WT to that",
        "age's body weight AND override the age-specific physiology (cardiac_output,",
        "hct, the fv_* volume fractions and the fq_* flow fractions) and the compound",
        "clearance clint_liver, all of which are tabulated per age in Table 1S and",
        "Supplementary 3. The six tabulated ages are 0.5, 2, 5, 12, 19 and 25 years.",
        sep = " "
      ),
      source_name = "BODYWT"
    )
  )

  compartmentData <- list(
    gut_lumen = list(analyte = "Deltamethrin", units = "umol", specimen = "not applicable", verified = TRUE),
    a_gut = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_metabolized = list(analyte = "Deltamethrin", units = "umol", specimen = "not applicable", verified = TRUE),
    a_fat_plasma = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_slowly_perfused_plasma = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_slowly_perfused = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_brain_plasma = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    a_brain = list(analyte = "Deltamethrin", units = "umol", specimen = "tissue", verified = TRUE),
    plasma = list(analyte = "Deltamethrin", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "6 months to 25 years (life-stage model; tabulated at 0.5, 2, 5, 12, 19, 25 years)",
    weight_range = "7.79 kg (6 months) to 81.74 kg (25 years), male (Supplementary 2, Table 1S)",
    sex_female_pct = 0,
    race_ethnicity = NA_character_,
    disease_state = "healthy (virtual life-stage population)",
    dose_range = paste(
      "Published simulations use a single daily oral dose of 1 mg/kg for 120 days to",
      "steady state, in males of 1, 5, 19 and 25 years (Monte-Carlo, 1000 subjects per",
      "age group); a 14-day single-daily-oral profile is also shown (Figure S5).",
      sep = " "
    ),
    regions = "United States (physiology adapted to NHANES 2005-2006 males)",
    notes = paste(
      "This is a forward IVIVE-PBPK simulation, not a fit, so n_subjects = 0. The",
      "generic structure and all chemical-specific parameters are taken from the",
      "authors' rat pyrethroid model (Song et al. 2019); age-specific physiology is",
      "adapted from published life-stage models (Clewell et al. 2004, Ruark et al.",
      "2017, Song et al. 2016, Wu et al. 2015) originally for women and rescaled to",
      "men. Age-specific hepatic clearance was built by IVIVE from expressed-enzyme",
      "in vitro intrinsic clearances (CYP1A2, 2B6, 2C8, 2C9, 2C19, 2D6, 3A4, 3A5, CES1,",
      "CES2) scaled by enzyme abundance and non-linear ontogeny curves; the resulting",
      "age-specific total hepatic clearances are tabulated in Supplementary 3 and are",
      "carried here as the compound parameter clint_liver rather than re-derived. The",
      "central finding is that brain (target-tissue) Cmax in children is comparable to",
      "or lower than in adults for the same exposure, because efficient CES-mediated",
      "hydrolysis matures rapidly after birth and clearance becomes liver-blood-flow",
      "limited across ages. The Monte-Carlo interindividual variability of the paper",
      "(Table 1: CVs on body weight, hematocrit, cardiac output, unbound fraction,",
      "brain flow and partition, liver volume and flow, metabolic constant, fat volume,",
      "lymphatic fraction) is documentation only here and is not encoded as IIV.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Every value in this model is a fixed literature/IVIVE constant; nothing
    # was estimated from individual data, so every parameter is fixed().
    # Age-specific physiology is given here at the 25-year-old adult male
    # reference (Supplementary 2, Table 1S, last column); younger ages are
    # simulated by overriding these through rxSolve(params = ...) with the
    # tabulated column for that age (see the vignette).
    # =====================================================================

    # --- Age-specific whole-body physiology (adult 25Y reference) ---------
    cardiac_output <- fixed(236.63) ; label("Total plasma flow / cardiac output, adult (L/h)") # Table 1S, CARDOUTPC at 25Y
    hct <- fixed(0.441) ; label("Hematocrit, adult (unitless)") # Table 1S, HCT at 25Y

    # Tissue volumes as a fraction of body weight (adult 25Y).
    fv_brain <- fixed(0.017) ; label("Brain volume, fraction of body weight (adult)") # Table 1S, VOLBRAINC at 25Y
    fv_fat <- fixed(0.239) ; label("Fat volume, fraction of body weight (adult)") # Table 1S, VOLFATC at 25Y
    fv_gut <- fixed(0.016) ; label("GI-tract volume, fraction of body weight (adult)") # Table 1S, VOLGIC at 25Y
    fv_liver <- fixed(0.0197) ; label("Liver volume, fraction of body weight (adult)") # Table 1S, VOLLIVERC at 25Y
    fv_rapidly_perfused <- fixed(0.0411) ; label("Rapidly-perfused volume, fraction of body weight (adult)") # Table 1S, VOLRPC at 25Y
    fv_slowly_perfused <- fixed(0.4521) ; label("Slowly-perfused volume, fraction of body weight (adult)") # Table 1S, VOLSPC at 25Y
    fv_blood <- fixed(0.0553) ; label("Blood volume, fraction of body weight (adult)") # Table 1S, VOLBLOODC at 25Y

    # Tissue plasma flows as a fraction of cardiac output (adult 25Y). The
    # five perfused-tissue fractions sum to 1.0 at every tabulated age.
    fq_brain <- fixed(0.1155) ; label("Brain plasma flow, fraction of cardiac output (adult)") # Table 1S, FRBRNC at 25Y
    fq_fat <- fixed(0.049) ; label("Fat plasma flow, fraction of cardiac output (adult)") # Table 1S, FRFATC at 25Y
    fq_liver <- fixed(0.215) ; label("Total liver plasma flow, fraction of cardiac output (adult)") # Table 1S, FRLIVC at 25Y
    fq_rapidly_perfused <- fixed(0.215) ; label("Rapidly-perfused plasma flow, fraction of cardiac output (adult)") # Table 1S, FRRPC at 25Y
    fq_slowly_perfused <- fixed(0.4055) ; label("Slowly-perfused plasma flow, fraction of cardiac output (adult)") # Table 1S, FRSPC at 25Y
    fq_liver_arterial <- fixed(0.05) ; label("Hepatic-arterial plasma flow, fraction of cardiac output") # Table 1S, FRLIVH (all ages)

    f_vascular_tissue <- fixed(0.05) ; label("Fraction of a diffusion-limited tissue that is vascular plasma") # Table 1S, VTBC (all ages)

    # --- Compound-independent partition and permeability parameters -------
    # (rat pyrethroid model, adapted to human; single values for all 8 pyrethroids)
    pc_fat <- fixed(68.7) ; label("Fat:plasma partition coefficient (unitless)") # Table 1S, PFAT
    pc_brain <- fixed(0.44) ; label("Brain:plasma partition coefficient (unitless)") # Table 1S, PBRN
    pc_slowly_perfused <- fixed(3.94) ; label("Slowly-perfused:plasma partition coefficient (unitless)") # Table 1S, PSP
    pa_fat_coef <- fixed(1.5) ; label("Fat permeability-area coefficient (L/h/kg^0.75)") # Table 1S, PAFC
    pa_brain_coef <- fixed(0.095) ; label("Brain permeability-area coefficient (L/h/kg^0.75)") # Table 1S, PABC
    pa_slowly_perfused_coef <- fixed(0.05) ; label("Slowly-perfused permeability-area coefficient (L/h/kg^0.75)") # Table 1S, PASPC

    # --- Absorption, binding and restrictive-clearance parameters ---------
    k_uptake <- fixed(5) ; label("Oral uptake rate constant (1/h)") # Table 1S, KA
    f_lymphatic <- fixed(0.086) ; label("Fraction of the oral dose absorbed via lymph, bypassing hepatic first pass (unitless)") # Table 1S, LYMPHSWTCH
    fu_plasma <- fixed(0.1) ; label("Fraction unbound in plasma (unitless)") # Table 1S, FuPLS
    kmf <- fixed(5) ; label("Empirical free-concentration adjustment factor for restrictive clearance (unitless)") # Table 1S, KMF (fitted)

    # --- Compound-specific hepatic clearance (deltamethrin, adult 25Y) ----
    clint_liver <- fixed(7418.58) ; label("Total hepatic intrinsic clearance, deltamethrin adult (L/h/kg liver)") # Supplementary 3, DLM total Clint at 25Y (= Table 1S VKM1C DLM)
  })

  model({
    # -------------------------------------------------------------------
    # 1. Calculated volumes (L) and plasma flows (L/h). Supplementary 4
    #    scenario file: volumes are fractions of body weight, flows are
    #    fractions of cardiac output (total plasma flow). Liver density is
    #    taken as 1 kg/L so liver weight (kg) equals liver volume (L).
    # -------------------------------------------------------------------
    BW <- WT
    VOLBRAIN <- fv_brain * BW
    VOLFAT <- fv_fat * BW
    VOLGI <- fv_gut * BW
    VOLLIVER <- fv_liver * BW
    VOLRP <- fv_rapidly_perfused * BW
    VOLSP <- fv_slowly_perfused * BW
    VOLBLOOD <- fv_blood * BW
    VOLPLS <- VOLBLOOD * (1 - hct) # plasma volume = blood volume x (1 - hematocrit)

    QLIV <- fq_liver * cardiac_output
    QLIVH <- fq_liver_arterial * cardiac_output # hepatic arterial inflow
    QLIVGI <- QLIV - QLIVH # portal inflow (drains the GI compartment)
    QFAT <- fq_fat * cardiac_output
    QRP <- fq_rapidly_perfused * cardiac_output
    QSP <- fq_slowly_perfused * cardiac_output
    QBRN <- fq_brain * cardiac_output

    # Liver, GI and rapidly-perfused tissues share the liver:plasma partition.
    # Table 1S footnote: scaled PLIV = 0.7524 * exp(log(4.451/0.7524) * exp(-0.1555 * kmf))
    PLIV <- 0.7524 * exp(log(4.451 / 0.7524) * exp(-0.1555 * kmf))
    PGI <- PLIV
    PRP <- PLIV

    # Permeability-area products for the diffusion-limited tissues, scaled to
    # tissue weight^0.75 (Chemical-Specific Parameters: scaled to volume^0.75).
    PAF <- pa_fat_coef * VOLFAT^0.75
    PAB <- pa_brain_coef * VOLBRAIN^0.75
    PASP <- pa_slowly_perfused_coef * VOLSP^0.75

    # Total hepatic clearance = intrinsic clearance per kg liver x liver weight.
    VKM1 <- clint_liver * VOLLIVER

    # -------------------------------------------------------------------
    # 2. Concentrations (umol/L). Free/bound split by fu_plasma; venous
    #    effluent of each compartment carries free tissue plus bound plasma.
    # -------------------------------------------------------------------
    CPLS <- plasma / VOLPLS
    CPLSf <- CPLS * fu_plasma # unbound arterial plasma
    CPLSb <- CPLS * (1 - fu_plasma) # bound plasma (rides through unchanged)

    CGIf <- (a_gut / VOLGI) * fu_plasma / PGI
    CLf <- (a_liver / VOLLIVER) * fu_plasma / PLIV # free liver concentration
    CVL <- CLf + CPLSb

    CVFf <- a_fat_plasma / (f_vascular_tissue * VOLFAT)
    CFf <- (a_fat / ((1 - f_vascular_tissue) * VOLFAT)) * fu_plasma / pc_fat
    CVF <- CPLSb + CVFf

    CRPf <- (a_rapidly_perfused / VOLRP) * fu_plasma / PRP
    CVRP <- CPLSb + CRPf

    CVSPf <- a_slowly_perfused_plasma / (f_vascular_tissue * VOLSP)
    CSPf <- (a_slowly_perfused / ((1 - f_vascular_tissue) * VOLSP)) * fu_plasma / pc_slowly_perfused
    CVSP <- CPLSb + CVSPf

    CVBf <- a_brain_plasma / (f_vascular_tissue * VOLBRAIN)
    CBf <- (a_brain / ((1 - f_vascular_tissue) * VOLBRAIN)) * fu_plasma / pc_brain
    CVB <- CPLSb + CVBf

    # Mixed venous return: five perfused tissues (GI drains into liver, not
    # directly to mixed venous, so it is not summed here).
    CV <- (QLIV * CVL + QFAT * CVF + QRP * CVRP + QSP * CVSP + QBRN * CVB) / cardiac_output

    # -------------------------------------------------------------------
    # 3. Absorption and metabolism.
    # -------------------------------------------------------------------
    avlbl_dose <- k_uptake * gut_lumen # rate of oral absorption (umol/h)
    LRAM <- VKM1 * (CLf / kmf) # first-order restrictive hepatic metabolism (umol/h)

    # -------------------------------------------------------------------
    # 4. Mass-balance ODE system (Supplementary 4, Model Equations).
    # -------------------------------------------------------------------
    d/dt(gut_lumen) <- -k_uptake * gut_lumen
    d/dt(a_gut) <- QLIVGI * (CPLSf - CGIf) + (1 - f_lymphatic) * avlbl_dose
    d/dt(a_liver) <- QLIVH * CPLSf + QLIVGI * CGIf - QLIV * CLf - LRAM
    d/dt(a_metabolized) <- LRAM
    d/dt(a_fat_plasma) <- QFAT * (CPLSf - CVFf) + PAF * (CFf - CVFf)
    d/dt(a_fat) <- PAF * (CVFf - CFf)
    d/dt(a_rapidly_perfused) <- QRP * (CPLSf - CRPf)
    d/dt(a_slowly_perfused_plasma) <- QSP * (CPLSf - CVSPf) + PASP * (CSPf - CVSPf)
    d/dt(a_slowly_perfused) <- PASP * (CVSPf - CSPf)
    d/dt(a_brain_plasma) <- QBRN * (CPLSf - CVBf) + PAB * (CBf - CVBf)
    d/dt(a_brain) <- PAB * (CVBf - CBf)
    d/dt(plasma) <- cardiac_output * (CV - CPLS) + f_lymphatic * avlbl_dose

    # -------------------------------------------------------------------
    # 5. Outputs (umol/L). Cc is plasma; Cbrain is the target-tissue
    #    (brain) total concentration the paper reports as its dose metric.
    # -------------------------------------------------------------------
    Cc <- CPLS
    Cbrain <- (a_brain_plasma + a_brain) / VOLBRAIN
    Cliver <- a_liver / VOLLIVER
    # Mass balance: total drug in the body plus metabolised should equal the
    # cumulative absorbed dose; useful as a structural check in simulation.
    Amass <- gut_lumen + a_gut + a_liver + a_metabolized +
      a_fat_plasma + a_fat + a_rapidly_perfused +
      a_slowly_perfused_plasma + a_slowly_perfused +
      a_brain_plasma + a_brain + plasma
  })
}

VanWart_2025_tobramycin <- function() {
  description <- paste(
    "Two-compartment population PK model for tobramycin in serum with zero-order",
    "intravenous input, first-order elimination, and a linked epithelial-lining-fluid",
    "(ELF) effect compartment, developed on time-matched steady-state serum and",
    "urea-corrected ELF (bronchoalveolar-lavage) concentrations from 16 adult patients",
    "with pneumonia originally published by Carcas 1999 (Van Wart 2025). The ELF",
    "compartment is driven by the serum concentration and does not remove drug from the",
    "central compartment, so it is a partitioned effect compartment rather than a",
    "distribution compartment: the pseudo-partition coefficient ppc = k13 / k30 = 0.49",
    "sets the steady-state ELF:serum ratio and the equilibration rate ke0 = k30 gives an",
    "equilibration half-life of about 12 minutes. Body-surface-area-normalised creatinine",
    "clearance, computed inside the model from weight, height, age, sex, and serum",
    "creatinine, acts on clearance as a power function referenced to 90 mL/min/1.73 m2.",
    "Central and peripheral volumes were constrained equal during fitting (Vss 20.8 L).",
    sep = " "
  )
  reference <- paste(
    "Van Wart SA, Trang M, Safir MC, Santulli AR, Rubino CM, Bhavnani SM.",
    "Population pharmacokinetic analysis of tobramycin in serum and ELF using data",
    "from patients with pneumonia.",
    "Antimicrob Agents Chemother. 2025;69(5):e0090824. doi:10.1128/aac.00908-24.",
    "Parameter estimates are from Table 1; the structural equations, the creatinine",
    "clearance derivation, and the residual-error structure are from the NONMEM control",
    "stream printed on page 5 of the supplemental material (aac.00908-24-s0001.pdf).",
    "The underlying serum and ELF concentration data were published by",
    "Carcas AJ, Garcia-Satue JL, Zapater P, Frias-Iniesta J. Tobramycin penetration into",
    "epithelial lining fluid of patients with pneumonia.",
    "Clin Pharmacol Ther. 1999;65(3):245-250. doi:10.1016/S0009-9236(99)70103-7.",
    sep = " "
  )
  vignette <- "VanWart_2025_tobramycin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at study entry.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the model only through the Cockcroft-Gault creatinine clearance and the",
        "Gehan and George body-surface-area formula computed in the supplemental control",
        "stream $PK block; there is no separate allometric term on any structural",
        "parameter. Weight ranged 50-84 kg (median 68 kg) across the 16 analysis subjects",
        "listed in the supplemental NONMEM data set (column WTKG). Dosing was 1 mg/kg, so",
        "weight also sets the dose amount in the analysis data set.",
        sep = " "
      ),
      source_name        = "WTKG"
    ),
    HT = list(
      description        = "Body height at study entry.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to compute body surface area by the Gehan and George formula",
        "BSA = 0.0235 * WT^0.51456 * HT^0.42246, which in turn normalises the male",
        "creatinine clearance to 1.73 m2 (supplemental control stream $PK). Height ranged",
        "147-190 cm across the 16 analysis subjects (column HTCM).",
        sep = " "
      ),
      source_name        = "HTCM"
    ),
    AGE = list(
      description        = "Age at study entry.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters only through the (140 - AGE) numerator of the Cockcroft-Gault creatinine",
        "clearance in the supplemental control stream $PK block. Age ranged 23-72 years",
        "(median 43 years) across the 16 analysis subjects (column AGE).",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Selects between the two creatinine clearance branches of the supplemental",
        "control stream $PK block. The source column is already named SEXF with the",
        "canonical 1 = female orientation, so no recoding is needed. The two branches are",
        "NOT symmetric as published: males get the Cockcroft-Gault estimate normalised to",
        "1.73 m2 by (1.73 / BSA) and no sex factor, while females get the 0.85 sex factor",
        "and no BSA normalisation, even though both are commented",
        "'mL/min/1.73 m^2' in the control stream. The asymmetry is reproduced verbatim",
        "here because the reported fixed effects are conditional on it; see the vignette",
        "Errata. 5 of the 16 analysis subjects (31.2%) were female.",
        sep = " "
      ),
      source_name        = "SEXF"
    ),
    CREAT = list(
      description        = "Serum creatinine at study entry.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Denominator of the Cockcroft-Gault creatinine clearance in the supplemental",
        "control stream $PK block, which uses the 72 constant appropriate to mg/dL.",
        "Serum creatinine ranged 0.6-1.2 mg/dL across the 16 analysis subjects",
        "(column SCR). Supplying CREAT in umol/L instead of mg/dL would inflate the",
        "denominator about 88-fold and is the single easiest way to misuse this model.",
        sep = " "
      ),
      source_name        = "SCR"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description        = "Body-surface-area-normalised creatinine clearance.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "CRCL is the covariate that actually acts on clearance, but it is not a data",
        "column for this model: the supplemental control stream derives it inside $PK",
        "from WT, HT, AGE, SEXF, and CREAT, and that derivation is reproduced in",
        "model() as the local variable crcl. It is documented here so the renal-function",
        "dependency is discoverable, and because a user who has CRCL but not its five",
        "inputs can bypass the derivation by supplying WT/HT/AGE/SEXF/CREAT values that",
        "reproduce their CRCL. Reference value 90 mL/min/1.73 m^2 (Table 1 reports CL",
        "'for a typical 90 mL/min/1.73m2 patient').",
        sep = " "
      ),
      source_name        = "CLCR"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "tobramycin", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "tobramycin", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    elf = list(
      analyte = "tobramycin", units = "mg/L",
      specimen = "epithelial lining fluid", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    age_range      = "23-72 years (median 43 years)",
    weight_range   = "50-84 kg (median 68 kg)",
    height_range   = "147-190 cm",
    sex_female_pct = 100 * 5 / 16,
    disease_state  = "Adults with pneumonia undergoing bronchoscopy with bronchoalveolar lavage.",
    renal_function = paste(
      "Serum creatinine 0.6-1.2 mg/dL. Derived creatinine clearance spans roughly",
      "60-160 mL/min/1.73 m^2 using the control stream's own sex-specific formulas; no",
      "subject had severe renal impairment.",
      sep = " "
    ),
    dose_range     = paste(
      "1 mg/kg tobramycin as a 30-minute IV infusion every 8 h, with dose adjustment as",
      "needed to reach peak and trough concentrations of approximately 8 and < 2 mg/L;",
      "bronchoscopy was performed at least 2 days after any adjustment so serum",
      "concentrations were at steady state.",
      sep = " "
    ),
    regions        = "Spain (the source clinical study of Carcas 1999).",
    n_observations = paste(
      "63 serum concentrations and 16 urea-corrected ELF concentrations (1 per patient,",
      "4 patients at each of 0.5, 2, 4, and 8 h), counted from the 95-row analysis data",
      "set printed on pages 7-9 of the supplemental material (16 dose records + 79",
      "observations). Serum was sampled at 0.5, 2, 4, and 8 h after the previous dose in",
      "every patient except subject 6, who has no 0.5 h serum record, so the count is 63",
      "and not the 64 that the Results text implies.",
      sep = " "
    ),
    notes          = paste(
      "Van Wart 2025 is a re-analysis: it fits a population PK model to individual",
      "dosing, demographic, renal-function, and time-matched serum/ELF concentration data",
      "published by Carcas 1999, which was chosen over the gentamicin and netilmicin",
      "candidates because it staggered the bronchoalveolar-lavage sampling times across",
      "the dosing interval and so can resolve the serum-to-ELF hysteresis. The complete",
      "16-subject analysis data set is printed on pages 7-9 of the supplemental material;",
      "the demographic ranges above are computed from it. Race and ethnicity are not",
      "reported.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural serum PK. Table 1 reports untransformed final estimates;
    # the supplemental control stream $PK block parameterises CL and Vc as
    # exp(log(theta) + eta), i.e. log-normal about the printed typical value.
    # ---------------------------------------------------------------------
    lcl <- log(3.26)
    label("Clearance at CRCL = 90 mL/min/1.73 m^2 (L/h)")                        # Table 1: CL (L/h) for a typical 90 mL/min/1.73m2 patient = 3.26 (%RSE 6.41); control stream THETA(1)
    lvc <- log(10.4)
    label("Central volume of distribution (L)")                                  # Table 1: Vc (L) = 10.4 (%RSE 7.03); control stream THETA(2)
    lq <- log(0.518)
    label("Intercompartmental clearance between serum and peripheral1 (L/h)")    # Table 1: CLd (L/h) = 0.518 (%RSE 25.8); control stream THETA(3)

    # Vp is deliberately absent from ini(). Table 1 footnote b: "Due to model
    # parameter identifiability concerns, Vc and Vp were constrained to have
    # the same value during model fitting", and the control stream implements
    # that as VP = TVVC -- the TYPICAL central volume, without eta. Encoding a
    # separate fixed(log(10.4)) would break the constraint the moment lvc were
    # re-estimated, so model() derives vp from lvc directly. The resulting
    # Vss is 20.8 L (0.297 L/kg at 70 kg), inside the 0.28-0.45 L/kg range of
    # the literature cited in that footnote.

    # ---------------------------------------------------------------------
    # ELF effect compartment. The control stream $DES writes
    #   DADT(3) = K13*A(1) - K30*A(3),  with S3 = S1 = VC,
    # which in concentration terms is d Celf / dt = k13 * Cc - k30 * Celf.
    # That is the canonical partitioned effect compartment
    #   d Celf / dt = ke0 * (ppc * Cc - Celf)
    # with ke0 = k30 (equilibration rate) and ppc = k13 / k30 (steady-state
    # ELF:serum ratio). Reported k13 / k30 = 1.81 / 3.69 = 0.4905, which
    # reproduces the paper's median simulated penetration ratio of 0.49
    # (Table 2), and log(2) / 3.69 = 0.188 h = 11.3 min reproduces the
    # "approximately 12 minutes" equilibrium half-life stated in the text.
    # ---------------------------------------------------------------------
    lke0 <- log(3.69)
    label("ELF equilibration rate constant k30 (1/h)")                           # Table 1: k30 (h-1) = 3.69 (%RSE 38.9); control stream THETA(5)
    lppc <- log(1.81 / 3.69)
    label("ELF:serum steady-state pseudo-partition coefficient k13/k30 (fraction)") # Table 1: k13 (h-1) = 1.81 (%RSE 23.1) and k30 (h-1) = 3.69 (%RSE 38.9); control stream THETA(4) / THETA(5)

    # ---------------------------------------------------------------------
    # Covariate effect.
    # ---------------------------------------------------------------------
    e_crcl_cl <- 0.685
    label("Power exponent of CRCL on clearance (unitless)")                      # Table 1: CL-CLcr exponent = 0.685 (%RSE 31.0); control stream THETA(6) in MU_1 = LOG(THETA(1)) + THETA(6)*LOG(CLCR/90.0)

    # ---------------------------------------------------------------------
    # Interindividual variability. Table 1 footnote c: "For model reduction
    # purposes, the magnitude of the interindividual variability for CL, Vc,
    # and k30 was constrained to have the same estimated value during model
    # fitting", implemented in the control stream as three separate
    # $OMEGA BLOCK(1) with SAME -- i.e. one shared variance and no
    # covariances between the three etas. nlmixr2 has no equality constraint
    # across omegas, so the same value is written three times; re-estimating
    # this model would free them.
    # ---------------------------------------------------------------------
    etalcl ~ 0.0488                                                              # Table 1: omega^2 CL = 0.0488 (22.4% CV, %RSE 31.6); control stream $OMEGA BLOCK(1) 0.1 (initial)
    etalvc ~ 0.0488                                                              # Table 1: omega^2 Vc = 0.0488 (22.4% CV, %RSE 31.6); control stream $OMEGA BLOCK(1) SAME
    etalke0 ~ 0.0488                                                             # Table 1: omega^2 k30 = 0.0488 (22.4% CV, %RSE 31.6); control stream $OMEGA BLOCK(1) SAME

    # ---------------------------------------------------------------------
    # Residual error. Control stream $ERROR:
    #   Y = F + F*EPS(1)*(1-FLGELF) + EPS(2)*FLGELF
    # so serum records get a purely proportional error and ELF records a
    # purely additive one; neither output carries both.
    # ---------------------------------------------------------------------
    propSd <- sqrt(0.0299)
    label("Proportional residual error, serum (fraction)")                       # Table 1: sigma^2 CCV for serum data = 0.0299 (17.3% CV, %RSE 36.3); sqrt(0.0299) = 0.173
    addSd_Celf <- sqrt(0.264)
    label("Additive residual error, ELF (mg/L)")                                 # Table 1: sigma^2 Additive for ELF data = 0.264 (0.514 SD, %RSE 61.6); sqrt(0.264) = 0.514
  })

  model({
    # 1. Derived covariate terms. Reproduced verbatim from the supplemental
    #    control stream $PK block:
    #
    #      BSA = 0.0235*(WTKG**0.51456)*(HTCM**0.42246)          ; Gehan and George
    #      IF(SEXF.EQ.0) THEN
    #       CLCR = (140.0-AGE)*WTKG/(72.0*SCR)*(1.73/BSA)        ; Males
    #      ELSE
    #       CLCR = ((140.0-AGE)*WTKG/(72.0*SCR))*0.85            ; Females
    #      ENDIF
    #
    #    The male branch normalises to 1.73 m^2 and applies no sex factor;
    #    the female branch applies the 0.85 sex factor and no BSA
    #    normalisation. That asymmetry is as published and is preserved here
    #    because the fixed effects were estimated under it -- see the
    #    vignette Errata. SEXF is 0/1, so the two branches are composed
    #    arithmetically rather than with an if/else.
    bsa <- 0.0235 * WT^0.51456 * HT^0.42246
    crcl <- (140 - AGE) * WT / (72 * CREAT) *
      ((1 - SEXF) * (1.73 / bsa) + SEXF * 0.85)

    # 2. Individual parameters.
    cl <- exp(lcl + etalcl) * (crcl / 90)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    # VP = TVVC in the control stream: the peripheral volume equals the
    # TYPICAL central volume and therefore carries no eta of its own.
    vp <- exp(lvc)

    ke0 <- exp(lke0 + etalke0)
    # Only k30 carries an eta in the source $PK block; k13 is a pure fixed
    # effect. Since k13 = ke0 * ppc, the eta that inflates ke0 must deflate
    # ppc by exactly the same factor, which makes ke0 * ppc =
    # exp(lke0 + lppc) = 1.81 1/h for every subject while the individual
    # penetration ratio ppc_i = k13 / k30_i varies as exp(-etalke0). This
    # is what generates the between-patient spread in the ELF:serum ratio
    # reported in Table 2 (SD 0.12 about a mean of 0.51).
    ppc <- exp(lppc - etalke0)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. $DES of the supplemental control stream:
    #      DADT(1) = K21*A(2) - K12*A(1) - K10*A(1)
    #      DADT(2) = K12*A(1) - K21*A(2)
    #      DADT(3) = K13*A(1) - K30*A(3)
    #    The ELF compartment is driven by the central compartment but does
    #    not return drug to it and does not appear in DADT(1): the paper
    #    states the structural model "used serum drug concentrations to drive
    #    appearance in and subsequent removal of drug from the ELF without
    #    impacting the serum PK data fitting". The ELF state holds a
    #    CONCENTRATION, not an amount, which is the algebraic equivalent of
    #    the control stream's S3 = VC scaling of an amount state.
    Cc <- central / vc

    d/dt(central) <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(elf) <- ke0 * (ppc * Cc - elf)

    # 5. Observations and residual error. Doses are administered to the
    #    central compartment as 30-minute zero-order infusions (rate set in
    #    the event table), matching RATE = 2 * AMT in the analysis data set.
    Celf <- elf

    Cc ~ prop(propSd)
    Celf ~ add(addSd_Celf)
  })
}

vanHest_2005_mycophenolic_acid <- function() {
  description <- "Two-compartment population PK model with time-lagged first-order absorption for mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) in adult renal transplant recipients co-treated with ciclosporin (van Hest 2005). Apparent clearance CL/F carries four covariates (creatinine clearance, serum albumin, ciclosporin daily dose, and sex) and apparent central volume V1/F carries two (creatinine clearance and serum albumin), all as centered power terms except sex, which is a multiplicative power-of-indicator factor of 1.11 for women. Log-normal inter-individual variability on Ka, CL/F, V1/F and V2/F, plus inter-occasion variability on Ka, CL/F and V1/F. Transcribed from the executable mrgsolve source distributed as S1 File of Maizaud 2025, which re-implements the van Hest 2005 model; the van Hest 2005 primary publication was not available at extraction time (see vignette Errata)."
  reference <- paste(
    "van Hest RM, van Gelder T, Vulto AG, Mathot RAA. Population pharmacokinetics",
    "of mycophenolic acid in renal transplant recipients.",
    "Clin Pharmacokinet. 2005;44(10):1083-1096.",
    "doi:10.2165/00003088-200544100-00006.",
    "Parameter values transcribed from the mrgsolve implementation in S1 File of",
    "Maizaud F, Arraki-Zava S, Sayadi H, Fromage Y, Marquet P, Woillard J-B,",
    "Monchaud C. Population pharmacokinetic modeling of missed mycophenolate",
    "mofetil doses: impact on exposure and exploration of mitigation strategies.",
    "PLoS One. 2025;20(8):e0330854. doi:10.1371/journal.pone.0330854.",
    sep = " "
  )
  vignette <- "Maizaud_2025_missed_mycophenolate_doses"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "mycophenolic acid (MPA)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "mycophenolic acid (MPA)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance in raw mL/min (NOT BSA-normalized). Covariate on both",
        "MPA apparent clearance CL/F and apparent central volume V1/F."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centered at 48 mL/min in both power terms, per the Maizaud 2025 S1 File",
        "'[MAIN]' block: pow(CLCR/48, thetaCLcR) on CL and pow(CLCR/48, thetaCLCR) on",
        "VC. Exponents are -0.12 on CL/F and -0.62 on V1/F, so lower renal function",
        "raises both apparent clearance and apparent central volume -- the expected",
        "direction for MPA, whose accumulating MPAG metabolite displaces MPA from",
        "albumin in renal impairment and so increases the free fraction. The mrgsolve",
        "'@covariates' simulation default is 60 mL/min. The estimating equation used",
        "in van Hest 2005 is not stated in any on-disk source; the 48 mL/min centering",
        "constant is presumably the cohort median."
      ),
      source_name = "CLCR"
    ),
    ALB = list(
      description = paste(
        "Serum albumin concentration. Covariate on both MPA apparent clearance CL/F",
        "and apparent central volume V1/F."
      ),
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centered at 30 g/L in both power terms, per the Maizaud 2025 S1 File",
        "'[MAIN]' block: pow(ALB/30, thetaALBCL) on CL and pow(ALB/30, thetaALB) on",
        "VC. The source values are already in SI g/L (the canonical unit as of the",
        "2026-06-19 register standardization), so no inline conversion is applied.",
        "Exponents are -1.07 on CL/F and -1.13 on V1/F: MPA is about 97% albumin-bound,",
        "so hypoalbuminaemia raises the free fraction and therefore both the apparent",
        "clearance and the apparent volume of the total-concentration model. The",
        "mrgsolve '@covariates' simulation default is 40 g/L."
      ),
      source_name = "ALB"
    ),
    CONMED_CSA_DOSE = list(
      description = paste(
        "Total daily dose of concomitant ciclosporin (cyclosporine A). Covariate on",
        "MPA apparent clearance CL/F."
      ),
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centered at 450 mg/day, per the Maizaud 2025 S1 File '[MAIN]' block:",
        "pow(CIC/450, thetaCIC) with exponent 0.31. Ciclosporin inhibits MRP2-mediated",
        "biliary efflux of MPAG and so suppresses enterohepatic recirculation of MPA;",
        "a higher ciclosporin dose therefore raises apparent MPA clearance. All",
        "subjects in the van Hest 2005 cohort received ciclosporin, so the companion",
        "binary CONMED_CSA is 1 throughout and is not carried as a separate covariate",
        "here. The mrgsolve '@covariates' simulation default is 300 mg/day. Setting",
        "CONMED_CSA_DOSE = 0 is NOT meaningful in this model: the power form is",
        "singular at zero and the cohort contains no ciclosporin-free subjects."
      ),
      source_name = "CIC"
    ),
    SEXF = list(
      description = "1 = female, 0 = male. Covariate on MPA apparent clearance CL/F.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "The Maizaud 2025 S1 File codes the source column as 'Gender : 0 : Gender",
        "(0 = male, 1 = female)' and applies pow(thetaGender, Gender) with",
        "thetaGender = 1.11, which matches the canonical SEXF orientation directly --",
        "no value inversion and no sign flip. Women have 11% higher apparent MPA",
        "clearance than men in this model. Encoded here as e_sexf_cl^SEXF. The S1",
        "simulation draws sex from rbinom(1000, 1, 0.50), i.e. a 50/50 split."
      ),
      source_name = "Gender"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for inter-occasion variability on Ka, CL/F and V1/F.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "van Hest 2005 reports inter-occasion variability (the paper's kappa terms,",
        "printed as the 'x' rows of Maizaud 2025 Table 1) on Ka, V1/F and CL/F. The",
        "NUMBER of occasions in the original analysis is not stated in any on-disk",
        "source, so this model declares the minimal two-occasion expansion following",
        "the Chen_2023_nemonoxacin.R / Jonsson_2011_ethambutol.R precedent, with the",
        "occasion-2 variances fix()ed equal to occasion 1 per the NONMEM",
        "'$OMEGA BLOCK(1) SAME' convention every shipped IOV model uses. Decomposed",
        "inside model() into binary indicators oc1 / oc2. For single-occasion",
        "simulation pass OCC = 1, which reproduces the Maizaud 2025 S1 File behaviour",
        "exactly (S1 adds one extra subject-level draw per IOV parameter inside the",
        "same exp(), i.e. it does not vary them by occasion). Users with genuine",
        "multi-occasion data should extend the pattern to their occasion count. See",
        "vignette Errata."
      ),
      source_name = "(not present in the Maizaud 2025 S1 File; see notes)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 140L,
    n_studies = 1L,
    sex_female_pct = 50,
    disease_state = "adult renal transplant recipients",
    co_medication = "ciclosporin (cyclosporine A) as the concomitant calcineurin inhibitor",
    renal_function = "creatinine clearance centered at 48 mL/min (presumed cohort median)",
    dose_range = "750-1250 mg MMF BID in the Maizaud 2025 simulations",
    notes = paste(
      "Population described in Maizaud 2025 Methods 2.1: 'The second model, by van",
      "Hest et al. [28], was based on data from 140 renal transplant recipients",
      "receiving ciclosporin.' The 50% female figure is the split used by the",
      "Maizaud 2025 S1 File virtual cohort (rbinom(1000, 1, 0.50)), not a reported",
      "cohort characteristic. Detailed baseline demographics (age, weight, race,",
      "time post-transplant) are in the van Hest 2005 primary publication, which was",
      "not available at extraction time."
    )
  )

  ini({
    # Structural parameters. Maizaud 2025 S1 File, van Hest '[PARAM] @annotated'
    # block; cross-checked against Maizaud 2025 Table 1 column 'Van Hest et al. [28]'.
    ltlag <- log(0.21); label("Absorption lag time (h)")                          # S1 TVTlag = 0.21; Table 1 'Tlag=0.21 h'
    lka   <- log(4.1);  label("Absorption rate constant (1/h)")                   # S1 TVKA = 4.1; Table 1 'Ka=4.1 h-1'
    lcl   <- log(32.5); label("Apparent clearance (CL/F, L/h)")                   # S1 TVCL = 32.5; Table 1 rounds this to 'CL/F=33 L/h'
    lvc   <- log(91);   label("Apparent central volume (V1/F, L)")                # S1 TVVC = 91; Table 1 'V1/F=91 L'
    lq    <- log(35);   label("Apparent intercompartmental clearance (Q/F, L/h)") # S1 TVQ = 35; Table 1 'Q/F=35 L/h'
    lvp   <- log(237);  label("Apparent peripheral volume (V2/F, L)")             # S1 TVVP = 237; Table 1 'V2/F=237 L'

    # Covariate effects. All continuous covariates enter as CENTERED power terms
    # (CRCL/48, ALB/30, CONMED_CSA_DOSE/450); sex enters as a power of the binary
    # indicator, e_sexf_cl^SEXF.
    e_crcl_cl            <- -0.12; label("Power exponent of creatinine clearance on CL/F (unitless)")   # S1 'thetaCLcR : -0.12 : Creatinine clearance effect on CL'
    e_alb_cl             <- -1.07; label("Power exponent of serum albumin on CL/F (unitless)")          # S1 'thetaALBCL : -1.07 : Albumin effect on CL'
    e_conmed_csa_dose_cl <-  0.31; label("Power exponent of ciclosporin daily dose on CL/F (unitless)") # S1 'thetaCIC : 0.31 : Ciclosporine dose effect on CL'
    e_sexf_cl            <-  1.11; label("Multiplicative factor on CL/F for women (unitless)")          # S1 'thetaGender : 1.11 : Gender covariate', applied as pow(thetaGender, Gender)
    e_crcl_vc            <- -0.62; label("Power exponent of creatinine clearance on V1/F (unitless)")   # S1 'thetaCLCR : -0.62 : Creatinine clearance effect on VC'
    e_alb_vc             <- -1.13; label("Power exponent of serum albumin on V1/F (unitless)")          # S1 'thetaALB : -1.13 : Albumin effect on VC'

    # IIV. S1 '[OMEGA] @annotated' entries are VARIANCES; Maizaud 2025 Table 1
    # reports the same quantities as SDs (omega). All four cross-check to two
    # decimals: 0.89^2 = 0.792 vs 0.802, 0.30^2 = 0.090 vs 0.091,
    # 0.77^2 = 0.593 vs 0.602, 0.84^2 = 0.706 vs 0.712. The S1 variances are used
    # here because they carry the extra significant figure.
    etalka ~ 0.802 # S1 'ETAKA : 0.802 : IIV on Ka'; Table 1 'omega Ka=0.89'
    etalcl ~ 0.091 # S1 'ETACL : 0.091 : IIV on CL'; Table 1 'omega CL/F=0.30'
    etalvc ~ 0.602 # S1 'ETAVC : 0.602 : IIV on VC'; Table 1 'omega V1/F=0.77'
    etalvp ~ 0.712 # S1 'ETAVP : 0.712 : IIV on VP'; Table 1 'omega V2/F=0.84'

    # IOV (the paper's kappa terms, printed as the 'x' rows of Table 1). Table 1
    # cross-check: 0.93^2 = 0.865 vs 0.86, 0.62^2 = 0.384 vs 0.384,
    # 0.29^2 = 0.084 vs 0.084. Occasion 2 is fix()ed equal to occasion 1 per the
    # NONMEM '$OMEGA BLOCK(1) SAME' convention; the number of occasions in the
    # original analysis is not on disk (see the OCC covariateData notes).
    etaiov_ka_1 ~ 0.86      # S1 'ETA_KA : 0.86 : IPV on Ka'; Table 1 'x Ka=0.93'
    etaiov_ka_2 ~ fix(0.86) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_1 ~ 0.384      # S1 'ETA_VC : 0.384 : IPV on VC'; Table 1 'x V1/F=0.62'
    etaiov_vc_2 ~ fix(0.384) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_1 ~ 0.084      # S1 'ETA_CL : 0.084 : IPV on CL'; Table 1 'x CL/F=0.29'
    etaiov_cl_2 ~ fix(0.084) # SAME-equivalent: equal to the occasion-1 IOV variance

    # Residual error. This is the ORIGINAL van Hest 2005 value as reported in
    # Maizaud 2025 Table 1, not the deliberately-decreased sigma = 0.001 that
    # Maizaud 2025 substituted for its own simulations (Methods 2.1).
    addSd <- 0.45; label("Additive residual error (mg/L)") # Table 1 'Original additive error sigma =0.45 mg/L'
  })

  model({
    # Decompose the integer occasion column into binary indicators multiplexing the
    # IOV etas. Pass OCC = 1 for single-occasion simulation.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # Individual parameters, exactly as written in the Maizaud 2025 S1 File '[MAIN]'
    # block:
    #   double CL = ((TVCL * pow(CLCR/48, thetaCLcR)) * pow(ALB/30, thetaALBCL)
    #                * pow(CIC/450, thetaCIC) * pow(thetaGender, Gender))
    #               * exp(ETACL + ETA_CL);
    #   double VC = TVVC * pow(CLCR/48, thetaCLCR) * pow(ALB/30, thetaALB)
    #               * exp(ETAVC + ETA_VC);
    # Q and Tlag carry neither IIV nor IOV in this model.
    cl <- exp(lcl + etalcl + iov_cl) *
      (CRCL / 48)^e_crcl_cl *
      (ALB / 30)^e_alb_cl *
      (CONMED_CSA_DOSE / 450)^e_conmed_csa_dose_cl *
      e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc + iov_vc) *
      (CRCL / 48)^e_crcl_vc *
      (ALB / 30)^e_alb_vc
    ka   <- exp(lka + etalka + iov_ka)
    q    <- exp(lq)
    vp   <- exp(lvp + etalvp)
    tlag <- exp(ltlag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODE system. MMF is dosed into the depot; MPA appears in the central
    # compartment. The MMF-to-MPA molar conversion is absorbed into the apparent
    # parameters CL/F and V/F, so dose amounts stay in mg of MMF.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}

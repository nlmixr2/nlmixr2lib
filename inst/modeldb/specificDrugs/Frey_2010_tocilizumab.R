Frey_2010_tocilizumab <- function() {
  description <- "Two-compartment population PK model for tocilizumab in adults with moderate-to-severe rheumatoid arthritis (Frey 2010), with parallel first-order linear and Michaelis-Menten elimination from the central compartment."
  reference <- "Frey N, Grange S, Woodworth T. Population pharmacokinetic analysis of tocilizumab in patients with rheumatoid arthritis. J Clin Pharmacol. 2010;50(7):754-766. doi:10.1177/0091270009350623"
  vignette <- "Frey_2010_tocilizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear CL with reference 1.8 m^2 per Frey 2010 Table III equation CL = 0.3 * (BSA/1.8)^0.7. Frey 2010 does not state which BSA formula was used; assume DuBois unless the source data dictionary states otherwise. BSA, BMI, and body weight were highly correlated in the dataset (Discussion p763); BSA was the body-size descriptor retained.",
      source_name        = "BSA"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Frey 2010 Methods (p757) parameterizes sex with male as reference (TVP = theta_P for males; TVP = theta_P * (1 + theta_SEX) for females) and reports theta_SEX = -0.16, i.e. CL is 16% lower in women than in men (Table III, p760). Maps directly to the canonical SEXF column.",
      source_name        = "SEX"
    ),
    HDLC = list(
      description        = "Baseline serum high-density lipoprotein cholesterol",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear CL with reference 54 mg/dL per Frey 2010 Table III equation CL = 0.3 * (HDL-C/54)^-0.2. Time-fixed at baseline. The paper interprets the HDLC effect on CL as a body-size surrogate (Discussion p763) rather than mechanistic.",
      source_name        = "HDL-C"
    ),
    RHEUMATOID_FACTOR = list(
      description        = "Baseline serum rheumatoid factor (autoantibody against the Fc portion of IgG)",
      units              = "U/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Frey 2010 fits the covariate on the natural-log scale: CL = 0.3 * (LRF/4.7)^0.1, where LRF = log(RHEUMATOID_FACTOR) and the reference LRF = 4.7 corresponds to RF = exp(4.7) ~= 110 U/mL (Table III, p760). The canonical column carries the raw RF value in U/mL; the log transform is applied inside model() as (log(RHEUMATOID_FACTOR) / log(110))^e_lrf_cl, which is algebraically identical to the paper's (LRF/4.7)^0.1 form.",
      source_name        = "RF"
    ),
    TPRO = list(
      description        = "Baseline total serum protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on V1 with reference 74 g/L per Frey 2010 Table III equation V1 = 3.5 * (PROT/74)^-1.1. Frey 2010 retains both TPRO (negative exponent) and ALB (positive exponent) on V1 and notes there is no clear mechanistic explanation; the joint effect may reflect serum-volume modifications (Discussion p763).",
      source_name        = "PROT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects on both V1 (positive exponent) and Vm (negative exponent) with reference 38 g/L per Frey 2010 Table III equations V1 = 3.5 * (ALBU/38)^0.7 and VM = 7.5 * (ALBU/38)^-0.4.",
      source_name        = "ALBU"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (measured Cockcroft-Gault method per the source paper's clinical-chemistry panel)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Vm with reference 106 mL/min per Frey 2010 Table III equation VM = 7.5 * (CRCL/106)^0.2. Frey 2010 does not BSA-normalize the creatinine clearance; the column carries the raw CrCl in mL/min, NOT the canonical BSA-normalized mL/min/1.73 m^2 form. Documented here so downstream simulation does not double-correct for body size (BSA already enters CL separately).",
      source_name        = "CRCL"
    ),
    SMOKE = list(
      description        = "Current-smoker indicator at baseline, 1 = smoker, 0 = non-smoker",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-smoker)",
      notes              = "Multiplicative fractional effect on Vm: VM_smoker = VM * (1 + 0.11) per Frey 2010 Table III row 'Smoking on VM' (effect +11%, p760). About 18% of the pooled cohort were smokers (Frey 2010 Results p759).",
      source_name        = "Smoking"
    )
  )

  population <- list(
    n_subjects     = 1793L,
    n_observations = 7415L,
    n_studies      = 4L,
    age_range      = "18-89 years",
    age_median     = "52 years",
    weight_range   = "38-150 kg",
    weight_median  = "70 kg",
    sex_female_pct = 82,
    race_ethnicity = c(White = 77, Asian = 8, "American Indian or Alaskan native" = 7, Black = 4, Other = 4),
    disease_state  = "Moderate-to-severe rheumatoid arthritis (adults). Patients with inadequate response to methotrexate (OPTION), inadequate response to traditional DMARDs (TOWARD), inadequate response to anti-TNF therapy (RADIATE), or as monotherapy (AMBITION).",
    dose_range     = "4 or 8 mg/kg by 1-hour intravenous infusion every 4 weeks for 24 weeks (per body weight).",
    regions        = "International multi-regional (4 phase III studies pooled).",
    notes          = "Baseline demographics from Frey 2010 Table I (4 phase III studies: OPTION N=396, TOWARD N=718, RADIATE N=341, AMBITION N=338; total 1793 subjects, 7415 PK samples). Smoking status: ~82% non-smokers / ~18% smokers across the pooled cohort. Concomitant medications: methotrexate (1227 patients), folic acid (1512), corticosteroids (879), NSAIDs (1010), other DMARDs / immunosuppressants in smaller subgroups (Frey 2010 Results p757-758); none of these concomitant drugs were retained as PK covariates."
  )

  ini({
    # Structural PK parameters - Frey 2010 Table II final-model estimates (p759). Reference
    # covariate values for the typical subject: male, BSA = 1.8 m^2, HDL-C = 54 mg/dL,
    # log(RF) = 4.7 (RF ~= 110 U/mL), total protein = 74 g/L, albumin = 38 g/L,
    # creatinine clearance = 106 mL/min, non-smoker.
    lcl <- log(0.3); label("Linear clearance CL (L/day)")                                       # Frey 2010 Table II, CL row
    lvc <- log(3.5); label("Central volume of distribution V1 (L)")                             # Frey 2010 Table II, V1 row
    lq  <- log(0.2); label("Inter-compartmental clearance Q (L/day)")                           # Frey 2010 Table II, Q row
    lvp <- log(2.9); label("Peripheral volume of distribution V2 (L)")                          # Frey 2010 Table II, V2 row

    # Parallel Michaelis-Menten elimination from the central compartment.
    lvmax <- log(7.5); label("Maximum Michaelis-Menten elimination rate Vmax (mg/day)")           # Frey 2010 Table II, VM row
    lkm   <- log(2.7); label("Michaelis-Menten constant Km (ug/mL)")                              # Frey 2010 Table II, KM row

    # Covariate exponents and multiplicative effects - Frey 2010 Table II ("Covariate effects")
    # and Table III equations (p760).
    e_bsa_cl  <-  0.7;  label("Power exponent of BSA/1.8 on linear CL (unitless)")               # Frey 2010 Table II / Table III: BSA on CL
    e_sexf_cl <- -0.16; label("Fractional change in linear CL for female sex (unitless)")        # Frey 2010 Table II: Sex on CL = -16%
    e_hdlc_cl <- -0.2;  label("Power exponent of HDLC/54 on linear CL (unitless)")               # Frey 2010 Table II / Table III: HDL-C on CL
    e_lrf_cl  <-  0.1;  label("Power exponent of log(RF)/log(110) on linear CL (unitless)")      # Frey 2010 Table II / Table III: log(RF) on CL
    e_tpro_vc   <- -1.1;  label("Power exponent of TPRO/74 on Vc (unitless)")                      # Frey 2010 Table II / Table III: total protein on V1
    e_alb_vc    <-  0.7;  label("Power exponent of ALB/38 on Vc (unitless)")                       # Frey 2010 Table II / Table III: albumin on V1
    e_alb_vmax  <- -0.4;  label("Power exponent of ALB/38 on Vmax (unitless)")                     # Frey 2010 Table II / Table III: albumin on VM
    e_crcl_vmax <-  0.2;  label("Power exponent of CRCL/106 on Vmax (unitless)")                   # Frey 2010 Table II / Table III: creatinine CL on VM
    e_smk_vmax  <-  0.11; label("Fractional change in Vmax for current smoker (unitless)")         # Frey 2010 Table II: Smoking on VM = +11%

    # Inter-individual variability: Frey 2010 Table II reports CV% (linear-scale) and the
    # off-diagonal correlation coefficients r between the four etas (CL, V1, V2, Vm).
    # Convert each CV% to NONMEM-style log-normal variance via omega^2 = log(1 + CV^2):
    #   CL  CV 39% -> omega^2 = log(1 + 0.39^2) = 0.1416, omega = 0.3762
    #   V1  CV 37% -> omega^2 = log(1 + 0.37^2) = 0.1284, omega = 0.3583
    #   V2  CV 66% -> omega^2 = log(1 + 0.66^2) = 0.3614, omega = 0.6012
    #   Vm  CV 54% -> omega^2 = log(1 + 0.54^2) = 0.2562, omega = 0.5062
    # Off-diagonal covariances cov_ij = r_ij * omega_i * omega_j with the six r values
    # reported in Table II ("Interindividual variability" section):
    #   r(CL,V1)= 0.6 -> cov =  0.6 * 0.3762 * 0.3583 =  0.0809
    #   r(CL,V2)=-0.1 -> cov = -0.1 * 0.3762 * 0.6012 = -0.0226
    #   r(CL,Vm)=-0.5 -> cov = -0.5 * 0.3762 * 0.5062 = -0.0952
    #   r(V1,V2)= 0.5 -> cov =  0.5 * 0.3583 * 0.6012 =  0.1077
    #   r(V1,Vm)= 0.2 -> cov =  0.2 * 0.3583 * 0.5062 =  0.0363
    #   r(V2,Vm)= 0.2 -> cov =  0.2 * 0.6012 * 0.5062 =  0.0609
    etalcl + etalvc + etalvp + etalvmax ~ c(0.1416,
                                            0.0809,  0.1284,
                                           -0.0226,  0.1077, 0.3614,
                                           -0.0952,  0.0363, 0.0609, 0.2562)

    # Residual error - Frey 2010 Table II ("Residual error" section, p759). Combined
    # additive + proportional model: Cobs = Cpred * (1 + eps_prop) + eps_add.
    # NOTE: the published Table II header swaps the Greek symbols
    #   "Additive (sigma_prop), ug/mL = 2.4" and
    #   "Proportional (sigma_add), %    = 22"
    # the row labels and units are correct (Additive 2.4 ug/mL, Proportional 22%) but the
    # parenthetical sigma_prop / sigma_add subscripts are swapped relative to the equation
    # in the Methods (p756), where eps1 carries sigma^2_prop and eps2 carries sigma^2_add.
    # The values 2.4 ug/mL and 22% are unambiguous; only the subscripts are mislabelled.
    # See vignette "Errata" section.
    propSd <- 0.22; label("Proportional residual error (fraction)")                              # Frey 2010 Table II "Proportional" row, 22%
    addSd  <- 2.4;  label("Additive residual error (ug/mL)")                                    # Frey 2010 Table II "Additive" row, 2.4 ug/mL
  })
  model({
    # Individual PK parameters with the Frey 2010 final-model covariate equations.
    cl <- exp(lcl + etalcl) *
          (BSA / 1.8)^e_bsa_cl *
          (1 + e_sexf_cl * SEXF) *
          (HDLC / 54)^e_hdlc_cl *
          (log(RHEUMATOID_FACTOR) / log(110))^e_lrf_cl
    vc <- exp(lvc + etalvc) *
          (TPRO / 74)^e_tpro_vc *
          (ALB / 38)^e_alb_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
    vmax <- exp(lvmax + etalvmax) *
          (ALB / 38)^e_alb_vmax *
          (CRCL / 106)^e_crcl_vmax *
          (1 + e_smk_vmax * SMOKE)
    km <- exp(lkm)

    # Two-compartment IV PK with parallel linear and Michaelis-Menten elimination
    # from the central compartment. Concentration in ug/mL; central in mg
    # (1 mg/L = 1 ug/mL).
    Cc <- central / vc

    d/dt(central)     <- -(cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    Cc ~ add(addSd) + prop(propSd)
  })
}

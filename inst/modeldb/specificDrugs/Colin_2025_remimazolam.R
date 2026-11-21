Colin_2025_remimazolam <- function() {
  description <- paste(
    "Joint population pharmacokinetic model for intravenous remimazolam and",
    "its inactive carboxylic-acid metabolite CNS 7054 in 31 children and",
    "adolescents aged 6 to 18 years undergoing procedural sedation (Colin",
    "2025, NCT04851717). Remimazolam follows a three-compartment model and",
    "is assumed to be eliminated entirely by conversion to CNS 7054, which",
    "reaches its one-compartment central volume through a chain of five",
    "transit compartments and is eliminated linearly. Every clearance and",
    "volume is scaled a priori by body weight (exponents 0.75 and 1,",
    "reference 70 kg); age, sex and fentanyl coadministration were tested on",
    "remimazolam clearance and none was retained.",
    sep = " "
  )
  reference <- paste(
    "Colin PJ, Bichajian LH, Curt VR, Koomen JV, Stohr T, Struys MMRF,",
    "Mason KP. Pharmacokinetics and Pharmacodynamics of Remimazolam for",
    "Procedural Sedation in Children and Adolescents. Anesthesiology.",
    "2025;143(2):368-382. doi:10.1097/ALN.0000000000005560. PMC12227205.",
    sep = " "
  )
  vignette <- "Colin_2025_remimazolam"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds, verified against the Supplemental Digital
  # Content 1 'final PopPK model' control stream ($MODEL COMP(CENTRAL)
  # COMP(PERIPH) COMP(PERIPH2) COMP(MET) COMP(TRANS1)..COMP(TRANS5);
  # $ERROR PARENT = (A(1)/V1)*1000 'Conversion from ug.mL-1 to ng.mL-1').
  # Doses are in mg and volumes in L. The transit and CNS 7054 states carry
  # remimazolam-mass equivalents: the converted remimazolam amount moves
  # through the chain unchanged and the molecular-weight factor is applied
  # only when the CNS 7054 concentration is formed.
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_cns7054 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "A priori allometric scaling of every clearance (exponent 0.75) and volume (exponent 1) to a 70 kg reference (Methods, Population Pharmacokinetic Modeling; Supplemental Digital Content 1 'FSIZE = WEIGHT/70'). The mean transit time is not scaled. Cohort mean 50.6 kg (SD 21.2, range 20.5-106) (Table 3).",
      source_name = "WEIGHT"
    )
  )

  # Tested on remimazolam clearance by likelihood-ratio test and not
  # retained (Results, Population Pharmacokinetic Model): age (dOFV -2.21,
  # P = 0.13), sex (dOFV -2.38, P = 0.12) and fentanyl coadministration
  # (dOFV -0.40, P = 0.52). No point estimates are reported.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Tested on remimazolam clearance; not retained (dOFV -2.21, P = 0.13). Cohort mean 12.3 years (SD 3.6, range 6-18)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = "Tested on remimazolam clearance; not retained (dOFV -2.38, P = 0.12). The paper does not state its coding. 12 of 31 patients were girls."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "6-18 years",
    age_mean = "12.3 years (SD 3.6)",
    weight_range = "20.5-106 kg",
    weight_mean = "50.6 kg (SD 21.2)",
    sex_female_pct = 100 * 12 / 31,
    disease_state = "Children and adolescents (ASA physical status I-III) undergoing elective diagnostic or therapeutic procedures under procedural sedation of less than 2 h with spontaneous breathing: lumbar puncture, endoscopy, bone marrow aspiration and other procedures.",
    dose_range = "Intravenous remimazolam titrated to a University of Michigan Sedation Scale score of 2-3: initial bolus 125-150 ug/kg (max 5-7 mg) or initial infusion 5.0-7.5 ug/kg/min, with top-up boluses at least 3 min apart (Table 2, revised protocol). 21 of 31 patients also received intravenous fentanyl 2-5 min before the first remimazolam bolus.",
    regions = "United States",
    n_observations = "170 remimazolam and 160 CNS 7054 venous plasma concentrations above the lower limit of quantification (2 and 20 ng/mL) from 171 samples.",
    notes = "First age cohort (6-18 years) of the phase 2/3 open-label trial CNS7056-026 (EudraCT 2020-004118-37, NCT04851717), four treatment arms: bolus (n = 5), bolus + fentanyl (n = 11), infusion (n = 5), infusion + fentanyl (n = 10) (Table 3). Fitted in NONMEM 7.4 with FOCE-I. The pharmacodynamic model for the sedation scale is modellib('Colin_2025_remimazolam_umss')."
  )

  ini({
    # ---------------------------------------------------------------------
    # Remimazolam disposition, typical values for a 70 kg patient. The
    # NONMEM model estimates every parameter as THETA on the log scale
    # (Supplemental Digital Content 1: 'V1 = EXP(THETA(1)+ETA(1))*FSIZE').
    # The values below are the final log-scale estimates that the
    # supplement's UMSS control stream carries forward as '$THETA ... FIX';
    # each back-transforms to the Table 4 estimate quoted alongside.
    # ---------------------------------------------------------------------
    lvc <- 0.223
    label("Remimazolam central volume of distribution V1 for 70 kg (L)") # SDC1 UMSS stream THETA(1) = 0.223; Table 4 V1 = 1.25 L/70 kg (95% CI 0.77-2.10)
    lvp <- 2.38
    label("Remimazolam fast peripheral volume V2 for 70 kg (L)") # SDC1 UMSS stream THETA(2) = 2.38; Table 4 V2 = 10.8 L/70 kg (7.45-18.1)
    lvp2 <- 2.51
    label("Remimazolam slow peripheral volume V3 for 70 kg (L)") # SDC1 UMSS stream THETA(3) = 2.51; Table 4 V3 = 12.3 L/70 kg (8.10-19.3)
    lcl <- -0.359
    label("Remimazolam elimination clearance CL for 70 kg (L/min)") # SDC1 UMSS stream THETA(4) = -0.359; Table 4 CL = 0.70 L/min/70 kg (0.52-0.94)
    lq <- 0.519
    label("Remimazolam fast distributional clearance Q2 for 70 kg (L/min)") # SDC1 UMSS stream THETA(5) = 0.519; Table 4 Q2 = 1.68 L/min/70 kg (1.00-2.74)
    lq2 <- -1.47
    label("Remimazolam slow distributional clearance Q3 for 70 kg (L/min)") # SDC1 UMSS stream THETA(6) = -1.47; Table 4 Q3 = 0.23 L/min/70 kg (0.15-0.33)

    # ---------------------------------------------------------------------
    # CNS 7054. The fraction of remimazolam metabolised to CNS 7054 is fixed
    # to 1 (Methods: 'we fixed the fraction metabolized to 1.00'), so the
    # whole remimazolam elimination flux k10 * A(1) enters the transit chain
    # (SDC1 $DES 'DADT(5) = k10*A(1) - ktr*A(5)').
    # ---------------------------------------------------------------------
    lvc_cns7054 <- 2.5
    label("CNS 7054 volume of distribution Vm for 70 kg (L)") # SDC1 UMSS stream THETA(7) = 2.5; Table 4 Vm = 12.2 L/70 kg (10.2-15.7)
    lcl_cns7054 <- -2.69
    label("CNS 7054 elimination clearance CLm for 70 kg (L/min)") # SDC1 UMSS stream THETA(8) = -2.69 (L/min); Table 4 CLm = 4.07 L/h/70 kg (2.52-5.88), estimated in L/min (Table 4 footnote); exp(-2.69) * 60 = 4.07
    lmtt <- 1.38
    label("Mean transit time of CNS 7054 formation (min)") # SDC1 UMSS stream THETA(9) = 1.38; Table 4 MTT = 3.97 min (3.85-5.09)

    # ---------------------------------------------------------------------
    # Between-subject variability. Table 4 footnote: 'IIV ... expressed as
    # a coefficient of variation (CV%) calculated as the square root of the
    # variance estimated by NONMEM', so omega^2 = (CV/100)^2. V1, V3, Q2, Q3
    # and MTT carry an IIV that was not estimated but 'fixed to 10%'
    # (Results), coded as '$OMEGA 1E-2 FIX'.
    # ---------------------------------------------------------------------
    etalvc ~ fixed(0.01) # SDC1 $OMEGA 1E-2 FIX 'IIV_V1'; Results: held at 10%
    etalvp ~ 0.853776 # Table 4 IIV V2 = 92.4% (64.3-137); 0.924^2
    etalvp2 ~ fixed(0.01) # SDC1 $OMEGA 1E-2 FIX 'IIV_V3'
    etalcl ~ 0.398161 # Table 4 IIV CL = 63.1% (45.8-91.0); 0.631^2
    etalq ~ fixed(0.01) # SDC1 $OMEGA 1E-2 FIX 'IIV_Q2'
    etalq2 ~ fixed(0.01) # SDC1 $OMEGA 1E-2 FIX 'IIV_Q3'
    etalvc_cns7054 ~ 0.245025 # Table 4 IIV Vm = 49.5% (37.5-72.0); 0.495^2
    etalcl_cns7054 ~ 0.417316 # Table 4 IIV CLm = 64.6% (43.6-101); 0.646^2
    etalmtt ~ fixed(0.01) # SDC1 $OMEGA 1E-2 FIX 'BSV_MTT'

    # ---------------------------------------------------------------------
    # Residual error. SDC1 $ERROR: SD = SQRT(IPRED^2 * PROP^2 + ADD^2) with
    # $SIGMA 1 FIX, i.e. a combined error whose proportional and additive
    # SDs add in variance. PROP = EXP(THETA(10)), ADD = EXP(THETA(11)) and
    # so on; Table 4 reports them as SDs. The remimazolam additive SD is
    # EXP(0 FIX) = 1 ng/mL (';; SD = 1 ng.mL-1').
    # ---------------------------------------------------------------------
    propSd <- 0.411
    label("Proportional residual error on remimazolam (fraction)") # Table 4 PropRMZ = 41.1% (35.2-49.0); SDC1 UMSS stream THETA(10) = -0.888, exp(-0.888) = 0.411
    addSd <- fixed(1)
    label("Additive residual error on remimazolam (ng/mL)") # SDC1 THETA(11) = 0 FIX ';; SD = 1 ng.mL-1', exp(0) = 1
    propSd_cns7054 <- 0.175
    label("Proportional residual error on CNS 7054 (fraction)") # Table 4 PropCNS7054 = 17.5% (14.8-21.1); SDC1 UMSS stream THETA(12) = -1.74, exp(-1.74) = 0.176
    addSd_cns7054 <- 22.9
    label("Additive residual error on CNS 7054 (ng/mL)") # Table 4 AddCNS7054 = 22.9 ng/mL (14.4-37.2); SDC1 UMSS stream THETA(13) = 3.13, exp(3.13) = 22.9
  })

  model({
    # A priori allometric scaling to 70 kg (SDC1 'FSIZE = WEIGHT/70').
    fsize <- WT / 70

    vc <- exp(lvc + etalvc) * fsize
    vp <- exp(lvp + etalvp) * fsize
    vp2 <- exp(lvp2 + etalvp2) * fsize
    cl <- exp(lcl + etalcl) * fsize^0.75
    q <- exp(lq + etalq) * fsize^0.75
    q2 <- exp(lq2 + etalq2) * fsize^0.75
    vc_cns7054 <- exp(lvc_cns7054 + etalvc_cns7054) * fsize
    cl_cns7054 <- exp(lcl_cns7054 + etalcl_cns7054) * fsize^0.75
    mtt <- exp(lmtt + etalmtt)

    # Five transit compartments with ktr = (N + 1) / MTT, N = 5 (SDC1
    # 'N = 5', 'ktr = (N+1)/MTT'), so MTT spans the five transit states and
    # the CNS 7054 central compartment.
    ktr <- (5 + 1) / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_cns7054 <- cl_cns7054 / vc_cns7054

    d/dt(central) <- -k12 * central + k21 * peripheral1 - k13 * central +
      k31 * peripheral2 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(transit1) <- kel * central - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(central_cns7054) <- ktr * transit5 - kel_cns7054 * central_cns7054

    # Concentrations in ng/mL (mg/L * 1000). The CNS 7054 concentration is
    # multiplied by 439.3 / 425.3 (remimazolam / CNS 7054 molecular weights)
    # exactly as the SDC1 $ERROR block codes it ('MET =
    # (A(4)/Vm)*1000*439.3/425.3 ; to account for difference in MW'); Vm
    # was estimated with this factor in place.
    Cc <- central / vc * 1000
    Cc_cns7054 <- central_cns7054 / vc_cns7054 * 1000 * 439.3 / 425.3

    Cc ~ prop(propSd) + add(addSd)
    Cc_cns7054 ~ prop(propSd_cns7054) + add(addSd_cns7054)
  })
}

Salem_2023_paclitaxel <- function() {
  description <- paste(
    "Semi-physiological (Friberg) PK/PD model of paclitaxel- and",
    "carboplatin-induced myelosuppression fitted to first-cycle absolute",
    "neutrophil counts from electronic health records of adults treated",
    "with weekly paclitaxel, weekly paclitaxel/carboplatin or every-3-weeks",
    "paclitaxel/carboplatin in routine care (Salem 2023, n = 405). The PK",
    "layers are fixed literature models: a three-compartment paclitaxel",
    "model with saturable elimination and saturable distribution to the",
    "first peripheral compartment (covariates body surface area, sex, total",
    "bilirubin and age), and a two-compartment linear carboplatin model.",
    "Both drug concentrations inhibit proliferation through multiplicative",
    "linear effects (1 - SLOPE * C); the PD parameters (baseline ANC, mean",
    "transit time, feedback exponent and both slopes) were estimated with",
    "log-normal IIV on baseline ANC, MTT and both slopes and a proportional",
    "residual error on ANC."
  )
  reference <- paste(
    "Salem AM, Dvergsten E, Karovic S, Maitland ML, Gopalakrishnan M.",
    "Model-based approach to identify predictors of paclitaxel-induced",
    "myelosuppression in 'real-world' administration. CPT Pharmacometrics",
    "Syst Pharmacol. 2023;12(7):929-940. doi:10.1002/psp4.12963.",
    "Model code and final estimates in Appendix S1 (Pumas code and fit",
    "output); equations, Table S1 and Figure S2 in Appendix S2.",
    "The fixed PK layers reproduce Joerger et al. Clin Cancer Res.",
    "2006;12(7 Pt 1):2150-2157 (paclitaxel) and Joerger et al. Clin Cancer",
    "Res. 2007;13(21):6410-6418 (carboplatin), as transcribed in Appendix S1.",
    sep = " "
  )
  vignette <- "Salem_2023_paclitaxel"
  units <- list(
    time = "h",
    dosing = "umol",
    dosing_carboplatin = "mg",
    concentration = "umol/L",
    concentration_carboplatin = "mg/L",
    anc = "10^9 cells/L"
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value. Power scaling (BSA / 1.8)^e on the paclitaxel second-peripheral volume, VMEL, VMTR and Q (Appendix S1 @pre block). Cohort mean (SD) 1.8 (0.2) m^2 in patients < 70 years and 1.8 (0.3) m^2 in patients >= 70 years (Table 1).",
      source_name = "BSA"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "female (SEXF = 1)",
      notes = "The source codes SEX as 0 = female, 1 = male (Appendix S1 @covariates comment), so SEX = 1 - SEXF. The male multipliers 1.2 (VMEL, VMTR), 2.11 (KMTR) and 0.893 (K21) are kept as printed and applied as factor^(1 - SEXF).",
      source_name = "SEX"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value. Power scaling (TBILI / 7)^-0.167 on the paclitaxel VMEL (Appendix S1 @pre block; the covariate comment gives the unit as 'micromole', i.e. umol/L). Missing values were imputed with the dataset mean (Results). Cohort mean (SD) 8.4 (4.2) umol/L in patients < 70 years and 9.2 (5.4) in patients >= 70 years (Table 1).",
      source_name = "BILIRUBIN"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling (AGE / 56)^-0.352 on the paclitaxel VMEL (Appendix S1 @pre block). 319 patients were < 70 years and 86 were >= 70 years (Table 1).",
      source_name = "AGE"
    )
  )

  compartmentData <- list(
    central = list(analyte = "paclitaxel", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paclitaxel", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "paclitaxel", units = "umol", specimen = "plasma", verified = TRUE),
    central_carboplatin = list(analyte = "carboplatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_carboplatin = list(analyte = "carboplatin", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(analyte = "neutrophils", units = "10^9 cells/L", specimen = "not applicable", verified = TRUE),
    transit1 = list(analyte = "neutrophils", units = "10^9 cells/L", specimen = "not applicable", verified = TRUE),
    transit2 = list(analyte = "neutrophils", units = "10^9 cells/L", specimen = "not applicable", verified = TRUE),
    transit3 = list(analyte = "neutrophils", units = "10^9 cells/L", specimen = "not applicable", verified = TRUE),
    circ = list(analyte = "neutrophils", units = "10^9 cells/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 405L,
    n_studies = 1L,
    n_observations = "2274 ANC measurements (median 5 per patient, range 3-39)",
    age_range = "18 years and older; 319 patients < 70 years, 86 (21.2%) >= 70 years and 20 (5%) >= 80 years",
    sex_female_pct = 86.2,
    race_ethnicity = "Self-reported: 221 White (European ancestry), 55 Asian, 45 African American; 184 non-White in total (Table 1)",
    disease_state = "Adults with cancer treated with paclitaxel-containing chemotherapy in routine care: breast 58.7%, ovarian 8.1%, non-small cell lung cancer 7.4%, other sites the remainder",
    dose_range = "Weekly paclitaxel (median 80 mg/m^2, 1-h infusion; n = 167), weekly paclitaxel/carboplatin (paclitaxel median 80 mg/m^2; n = 179) or every-3-weeks paclitaxel/carboplatin (paclitaxel median 175 mg/m^2, 3-h infusion; n = 59); carboplatin as a 30-min infusion",
    regions = "United States (Inova Schar Cancer Institute, Fairfax, Virginia)",
    notes = "Retrospective electronic health record data, 2015-2019; first treatment cycle only. Baseline ANC mean (SD) 7.4 (5.4) x 10^9/L in patients < 70 years and 6.9 (4.4) in patients >= 70 years (Table 1). Infusion durations were assigned from institutional standards and doses were assumed to start 30 min after a same-day ANC sample (Methods)."
  )

  ini({
    # ---- Paclitaxel PK (fixed; Appendix S1 @pre block, Figure S2A) ----
    # Typical values for a woman with BSA 1.8 m^2, total bilirubin 7 umol/L
    # and age 56 years. Amounts in umol, concentrations in umol/L.
    lvc <- fixed(log(12.8)); label("Paclitaxel central volume Vc (L)") # Appendix S1: Vc = 12.8
    lvp2 <- fixed(log(252)); label("Paclitaxel second peripheral volume Vp at BSA 1.8 m^2 (L)") # Appendix S1: Vp = 252 * (BSA/1.8)^1.17
    lq2 <- fixed(log(20.1)); label("Paclitaxel intercompartmental clearance Q, central <-> second peripheral, at BSA 1.8 m^2 (L/h)") # Appendix S1: Q = 20.1 * (BSA/1.8)^0.724
    lvmax <- fixed(log(37.4)); label("Paclitaxel maximum elimination rate VMEL for the reference woman (umol/h)") # Appendix S1: VMEL = 37.4 * ...
    lkm <- fixed(log(0.53)); label("Paclitaxel Michaelis-Menten constant of elimination KMEL (umol/L)") # Appendix S1: KMEL = 0.53
    lvmtr <- fixed(log(169)); label("Paclitaxel maximum distribution rate VMTR, central -> first peripheral, for the reference woman (umol/h)") # Appendix S1: VMTR = 169 * ...
    lkmtr <- fixed(log(0.83)); label("Paclitaxel Michaelis-Menten constant of distribution KMTR in women (umol/L)") # Appendix S1: KMTR = 0.83 * 2.11^SEX
    lk21 <- fixed(log(1.15)); label("Paclitaxel first-order return rate K21, first peripheral -> central, in women (1/h)") # Appendix S1: K21 = 1.15 * 0.893^SEX

    e_bsa_vp2 <- fixed(1.17); label("Power exponent of (BSA / 1.8) on Vp (unitless)") # Appendix S1: (BSA/1.8)^1.17
    e_bsa_q2 <- fixed(0.724); label("Power exponent of (BSA / 1.8) on Q (unitless)") # Appendix S1: (BSA/1.8)^0.724
    e_sexm_vmax <- fixed(1.2); label("Multiplicative male effect on VMEL, applied as e_sexm_vmax^(1 - SEXF) (unitless)") # Appendix S1: 1.2^SEX
    e_bsa_vmax <- fixed(0.842); label("Power exponent of (BSA / 1.8) on VMEL (unitless)") # Appendix S1: (BSA/1.8)^0.842
    e_tbili_vmax <- fixed(-0.167); label("Power exponent of (TBILI / 7 umol/L) on VMEL (unitless)") # Appendix S1: (BILIRUBIN/7)^(-0.167)
    e_age_vmax <- fixed(-0.352); label("Power exponent of (AGE / 56 years) on VMEL (unitless)") # Appendix S1: (AGE/56)^(-0.352)
    e_sexm_vmtr <- fixed(1.2); label("Multiplicative male effect on VMTR, applied as e_sexm_vmtr^(1 - SEXF) (unitless)") # Appendix S1: 1.2^SEX
    e_bsa_vmtr <- fixed(0.911); label("Power exponent of (BSA / 1.8) on VMTR (unitless)") # Appendix S1: (BSA/1.8)^0.911
    e_sexm_kmtr <- fixed(2.11); label("Multiplicative male effect on KMTR, applied as e_sexm_kmtr^(1 - SEXF) (unitless)") # Appendix S1: 2.11^SEX
    e_sexm_k21 <- fixed(0.893); label("Multiplicative male effect on K21, applied as e_sexm_k21^(1 - SEXF) (unitless)") # Appendix S1: 0.893^SEX

    # ---- Carboplatin PK (fixed; Appendix S1 @pre block, Figure S2A) ----
    lcl_carboplatin <- fixed(log(7.38)); label("Carboplatin clearance CL (L/h)") # Appendix S1: cl = 7.38
    lvc_carboplatin <- fixed(log(11.9)); label("Carboplatin central volume Vc (L)") # Appendix S1: vc = 11.9
    lq_carboplatin <- fixed(log(90.5)); label("Carboplatin intercompartmental clearance Q (L/h)") # Appendix S1: q = 90.5
    lvp_carboplatin <- fixed(log(8.23)); label("Carboplatin peripheral volume Vp (L)") # Appendix S1: vp = 8.23

    # ---- Myelosuppression PD (estimated; Appendix S1 fit output) ----
    lcirc0 <- log(6.2787); label("Baseline absolute neutrophil count CIRC0 (10^9 cells/L)") # Appendix S1 estimate tvPOP_CIRC0 = 6.2787; Table S1 6.3 (RSE 3%)
    lmtt <- log(143.71); label("Mean transit time MTT through the maturation chain (h)") # Appendix S1 estimate tvPOP_MTT = 143.71; Table S1 143.7 (RSE 3%)
    gamma <- 0.19585; label("Feedback exponent gamma on (CIRC0 / circ) (unitless)") # Appendix S1 estimate tvPOP_GAMMA = 0.19585; Table S1 0.2 (RSE 6%)
    lslope <- log(5.8194); label("Linear paclitaxel effect SLOPE on proliferation (L/umol)") # Appendix S1 estimate tvPOP_SLOP = 5.8194; Table S1 5.8 (RSE 5.1%)
    lslope_carboplatin <- log(0.0058817); label("Linear carboplatin effect SLOPEC on proliferation (L/mg)") # Appendix S1 estimate tvPOP_SLOPC = 0.0058817; Table S1 0.006 (RSE 34.2%)

    # IIV: Appendix S1 Omega_1 = diag(CIRC0, MTT, SLOPE) and Omega_2 = SLOPEC,
    # printed as variances. Table S1 reports sqrt(omega^2) x 100 as '% CV'
    # (48.3, 17.7, 33.3, 182.2).
    etalcirc0 ~ 0.23317 # Appendix S1 estimate Omega_1(1,1); Table S1 BSV CIRC0 48.3%
    etalmtt ~ 0.031576 # Appendix S1 estimate Omega_1(2,2); Table S1 BSV MTT 17.7%
    etalslope ~ 0.11131 # Appendix S1 estimate Omega_1(3,3); Table S1 BSV SLOPE 33.3%
    etalslope_carboplatin ~ 3.3201 # Appendix S1 estimate Omega_2(1,1); Table S1 BSV SLOPEC 182.2%

    # Appendix S1 @derived: ANC ~ Normal(E, E * sigma_prop), so sigma_prop is
    # the proportional SD.
    propSd <- 0.36922; label("Proportional residual error on ANC (fraction)") # Appendix S1 estimate sigma_prop = 0.36922; Table S1 36.9% (RSE 3.6%)
  })

  model({
    # ---- Paclitaxel PK (Appendix S1 @pre; reference BSA 1.8 m^2,
    # bilirubin 7 umol/L, age 56 years; SEX = 1 - SEXF) ----
    sexm <- 1 - SEXF
    vc <- exp(lvc)
    vp2 <- exp(lvp2) * (BSA / 1.8)^e_bsa_vp2
    q2 <- exp(lq2) * (BSA / 1.8)^e_bsa_q2
    vmax <- exp(lvmax) * e_sexm_vmax^sexm * (BSA / 1.8)^e_bsa_vmax *
      (TBILI / 7)^e_tbili_vmax * (AGE / 56)^e_age_vmax
    km <- exp(lkm)
    vmtr <- exp(lvmtr) * e_sexm_vmtr^sexm * (BSA / 1.8)^e_bsa_vmtr
    kmtr <- exp(lkmtr) * e_sexm_kmtr^sexm
    k21 <- exp(lk21) * e_sexm_k21^sexm

    # ---- Carboplatin PK ----
    cl_carboplatin <- exp(lcl_carboplatin)
    vc_carboplatin <- exp(lvc_carboplatin)
    q_carboplatin <- exp(lq_carboplatin)
    vp_carboplatin <- exp(lvp_carboplatin)

    Cc <- central / vc
    Cc_carboplatin <- central_carboplatin / vc_carboplatin

    # Saturable elimination and saturable central -> first-peripheral
    # transfer; the return is first order on the peripheral amount. The
    # transfer term enters d/dt(central) with a minus sign as in the
    # Appendix S1 code; Appendix S2 Eq. 1 prints it with a plus sign, which
    # would not conserve mass.
    elim <- vmax * Cc / (km + Cc)
    transfer <- vmtr * Cc / (kmtr + Cc)

    d/dt(central) <- -elim + k21 * peripheral1 - transfer + q2 / vp2 * peripheral2 - q2 / vc * central
    d/dt(peripheral1) <- -k21 * peripheral1 + transfer
    d/dt(peripheral2) <- -q2 / vp2 * peripheral2 + q2 / vc * central

    d/dt(central_carboplatin) <- -cl_carboplatin / vc_carboplatin * central_carboplatin +
      q_carboplatin / vp_carboplatin * peripheral1_carboplatin -
      q_carboplatin / vc_carboplatin * central_carboplatin
    d/dt(peripheral1_carboplatin) <- q_carboplatin / vc_carboplatin * central_carboplatin -
      q_carboplatin / vp_carboplatin * peripheral1_carboplatin

    # ---- Myelosuppression PD (Appendix S1 @vars / @dynamics; Appendix S2
    # Eqs. 4-11). Three transit compartments, so ktr = (3 + 1) / MTT. ----
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt <- exp(lmtt + etalmtt)
    slope <- exp(lslope + etalslope)
    slope_carboplatin <- exp(lslope_carboplatin + etalslope_carboplatin)
    ktr <- 4 / mtt

    edrug <- 1 - slope * Cc
    edrug_carboplatin <- 1 - slope_carboplatin * Cc_carboplatin
    feed <- abs(circ0 / circ)^gamma

    d/dt(prol) <- ktr * prol * edrug * feed * edrug_carboplatin - ktr * prol
    d/dt(transit1) <- ktr * prol - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(circ) <- ktr * transit3 - ktr * circ

    prol(0) <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    circ(0) <- circ0

    ANC <- circ
    ANC ~ prop(propSd)
  })
}

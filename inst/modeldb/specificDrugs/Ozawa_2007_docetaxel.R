Ozawa_2007_docetaxel <- function() {
  description <- "Three-compartment IV PK coupled with a modified Friberg-style semimechanistic-physiological PK/PD model for docetaxel-induced neutropenia in Japanese cancer patients (Ozawa 2007). The PD layer extends Friberg 2002 with an additional zero-order input compartment that captures the transient ANC increase attributable to dexamethasone premedication; alpha-1 acid glycoprotein modulates the linear drug-effect slope on the proliferating compartment via a power-law form. Per-subject baseline ANC is supplied as a covariate and is used to initialise the proliferation, transit, and circulating compartments."
  reference <- paste(
    "Ozawa K, Minami H, Sato H. (2007).",
    "Population pharmacokinetic and pharmacodynamic analysis for time courses",
    "of docetaxel-induced neutropenia in Japanese cancer patients.",
    "Cancer Sci 98(12):1985-1992. doi:10.1111/j.1349-7006.2007.00615.x.",
    "PD structure extends Friberg LE et al. (2002) J Clin Oncol 20(24):4713-4721",
    "(see modellib('Friberg_2002_paclitaxel') for the leukocyte arm of the original).",
    sep = " "
  )
  vignette <- "Ozawa_2007_docetaxel"
  paper_specific_compartments <- c("input")

  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L",
    anc           = "10^9 cells/L",
    aag           = "g/L"
  )

  covariateData <- list(
    AAG = list(
      description        = "Serum alpha-1 acid glycoprotein concentration, time-fixed at baseline.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Multiplicative power-form effect on the drug-effect slope: SLOPE = theta_SLOPE * (AAG / 0.94)^e_aag_slope (Ozawa 2007 Equation 11). Reference value 0.94 g/L (= 94 mg/dL) is the population-median constant used inside the published NONMEM control stream (Appendix I, `AGPm = 94`); Table 1 reports the cohort median as 90 mg/dL but the model fit was carried out with AGPm = 94. Source NM-TRAN column `AGP1` (Appendix I $INPUT) is reported in mg/dL; conversion to canonical g/L is `AAG_g_per_L = AGP_mg_per_dL / 100`. Cohort range 0.51-2.41 g/L (Table 1 `51-241 mg/dL`). Time-fixed per subject.",
      source_name        = "AGP1"
    ),
    NEUT = list(
      description        = "Per-subject baseline absolute neutrophil count (observed value before drug administration), used as the fixed initial condition for the proliferation, transit, and circulating compartments.",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source NM-TRAN column `BASE` (Appendix I $INPUT). Ozawa 2007 fixes Circ(t=0) at the per-subject observed baseline ANC (Methods page 1986: 'Circ (t = 0) was fixed at its observed value') rather than estimating a population baseline. The proliferating and three transit compartments are also initialised at this value (Methods: 'Prol (t = 0) = Transit1 (t = 0) = Transit2 (t = 0) = Transit3 (t = 0) = Circ (t = 0)'); the same value enters the feedback term as Circ0 in the (Circ0/Circ)^gamma1 ratio. The canonical NEUT register units are cells/mm^3; this model documents the per-paper unit override of 10^9 cells/L. To convert, NEUT_per_mm3 = NEUT_10e9_per_L * 1000.",
      source_name        = "BASE"
    )
  )

  population <- list(
    n_subjects     = 62L,
    n_studies      = 1L,
    age_range      = "21-77 years",
    age_median     = "55.5 years",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 82.3,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Japanese adult cancer patients (44 breast, 10 non-small cell lung, 3 head and neck, 5 other) including subjects with liver dysfunction or poor performance status who would normally be excluded from drug-development trials. ECOG performance status distribution 0/1/2/3 = 14/36/7/5.",
    dose_range     = "Docetaxel 30-60 mg/m^2 IV over 1 h every 3 weeks (the approved Japanese dose at the time was 60 mg/m^2; attending physicians could reduce based on liver function, performance status, or the extent of prior chemotherapy). Per-PS-stratum medians (Table 1): PS0 60 mg/m^2 (range 39-60), PS1 60 (40-63), PS2 60 (58-62), PS3 30 (26-60).",
    regions        = "Japan (National Cancer Center Hospital East, Kashiwa, Chiba)",
    co_medication  = "Dexamethasone administered before each docetaxel infusion to prevent emesis. The model's `input` compartment captures the resulting transient ANC increase but the dexamethasone dose itself is not encoded as an explicit covariate. Subjects who received granulocyte-colony stimulating factors after docetaxel were excluded from the analysis.",
    notes          = "Baseline characteristics from Ozawa 2007 Table 1. Median (range) baseline laboratory values: albumin 3.7 (2.6-4.5) g/dL, total bilirubin 0.6 (0.2-1.2) mg/dL, AST 24 (11-310) IU/L, ALT 20 (6-140) IU/L, alkaline phosphatase 249 (93-1382) IU/L, creatinine 0.6 (0.4-1.3) mg/dL, AGP 90 (51-241) mg/dL. Forty-six subjects had < 3 prior chemotherapy regimens, 16 had >= 3. Body weight range and median are not reported in Table 1 (dose is given as mg/m^2 only). 395 ANC observations across 62 subjects fed the PD analysis. The PK and PD parameters were estimated in two sequential steps (popPK first using FOCEI, then PD using FO with individual posthoc PK parameters as input columns) but the `ini()` block below encodes the joint typical-value model so simulation reproduces both layers in one rxode2 solve."
  )

  ini({
    # ---- Three-compartment IV PK (Ozawa 2007 Table 2; ADVAN11 TRANS4) -------
    # Parameter values are population means (theta) reported in Table 2.
    # The paper's PK structure mirrors the upstream popPK model of Bruno et al. 1996
    # (J Pharmacokin Biopharm 24:153) but the values were re-estimated on the
    # 62-subject Japanese cohort.
    lcl  <- log(35.7); label("CL: docetaxel clearance (L/h)")                                            # Table 2: theta_CL  = 35.7  (SE 1.30)
    lvc  <- log(6.94); label("V1: central volume of distribution (L)")                                   # Table 2: theta_V1  = 6.94  (SE 0.303)
    lq   <- log(5.58); label("Q2: intercompartmental clearance central<->peripheral1 (L/h)")             # Table 2: theta_Q2  = 5.58  (SE 0.356)
    lvp  <- log(7.39); label("V2: peripheral1 volume of distribution (L)")                               # Table 2: theta_V2  = 7.39  (SE 1.08)
    lq2  <- log(12.5); label("Q3: intercompartmental clearance central<->peripheral2 (L/h)")             # Table 2: theta_Q3  = 12.5  (SE 1.22)
    lvp2 <- log(225);  label("V3: peripheral2 volume of distribution (L)")                               # Table 2: theta_V3  = 225   (SE 46.2)

    # PK IIV (Table 2). Reported as CV%; converted to log-eta variance via log(1 + CV^2).
    etalcl ~ 0.1478   # log(1 + 0.399^2) = 0.1478; Table 2 omega_CL = 39.9% CV
    etalvc ~ 0.0673   # log(1 + 0.264^2) = 0.0673; Table 2 omega_V1 = 26.4% CV

    # PK proportional residual error (Table 2 sigma = 26.8% CV; Methods text: "proportional model").
    propSd <- 0.268; label("Proportional residual SD on docetaxel concentration (fraction)")          # Table 2: sigma     = 26.8% CV

    # ---- Friberg-extension myelosuppression PD (Ozawa 2007 Table 3) ----------
    # Six PD compartments: Prol (proliferating), Transit1..3 (maturation chain), Input (dexamethasone-
    # induced ANC bump), Circ (observed circulating neutrophils). See Equations 1-7 + Equation 11.
    # NONMEM control stream in Appendix I confirms the structure (ADVAN6, $MODEL block COMP=POOL,
    # DELAY1..3, INPUT2, OBS) and the constraint kprol = ktr = kcirc = 4/MTT (Methods).
    lmtt        <- log(113);   label("MTT: mean transit time through proliferation -> circulation chain (h)")    # Table 3: theta_MTT    = 113   (SE 4.62)
    lslope      <- log(17.9);  label("SLOPE: linear drug-effect slope on proliferation (1/(mg/L))")              # Table 3: theta_SLOPE  = 17.9  (SE 1.75)
    lgamma1     <- log(0.196); label("gamma1: feedback exponent on (Circ0/Circ) (unitless)")                     # Table 3: theta_gamma1 = 0.196 (SE 0.013)
    lmit        <- log(35.5);  label("MIT: mean input time of the dexamethasone-induced ANC bump (h)")           # Table 3: theta_MIT    = 35.5  (SE 5.08)
    lip0        <- log(5.19);  label("IP0: initial input-compartment amount, dexamethasone-driven (10^9 cells/L)") # Table 3: theta_IP0  = 5.19  (SE 1.51)
    e_aag_slope <- -1.38;      label("gamma2: AGP power-law exponent on SLOPE (unitless, negative)")             # Table 3: theta_gamma2 = -1.38 (SE 0.287)

    # PD IIV (Table 3). Reported as CV%; converted to log-eta variance via log(1 + CV^2).
    # Appendix I $OMEGA fixes ETA(3) on gamma1 and ETA(4) on MIT to 0, so no IIV on those.
    etalmtt   ~ 0.01225  # log(1 + 0.111^2) = 0.01225; Table 3 omega_MTT   = 11.1% CV
    etalslope ~ 0.3818   # log(1 + 0.682^2) = 0.3818;  Table 3 omega_SLOPE = 68.2% CV
    etalip0   ~ 0.7930   # log(1 + 1.10^2)  = 0.7930;  Table 3 omega_IP0   = 110% CV

    # PD residual error: Appendix I $ERROR Y = F * EXP(ERR(1)) -- log-normal residual.
    # Table 3 reports sigma as CV%; mapped to lnorm() addSd on the log scale via the
    # small-CV approximation addSd ~= CV/100 (paper's convention, since the source labels
    # the value 'CV%').
    addSd_ANC <- 0.291; label("Lognormal residual SD on log(ANC) (unitless on log scale)")               # Table 3: sigma = 29.1% CV
  })

  model({
    # ---- PK individual parameters ----
    # CL and V1 carry IIV (Table 2); Q2, V2, Q3, V3 are typical-only (no eta in Table 2).
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    # PK micro-constants (Appendix I $PK)
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ---- PD individual parameters ----
    # MTT, SLOPE, IP0 carry IIV (Table 3); gamma1 and MIT typical-only (Appendix I $OMEGA fixes
    # their etas to 0).
    mtt    <- exp(lmtt + etalmtt)
    gamma1 <- exp(lgamma1)
    mit    <- exp(lmit)
    ip0    <- exp(lip0 + etalip0)

    # AGP covariate effect on SLOPE (Equation 11). The published model normalises by AGPm = 94 mg/dL
    # (= 0.94 g/L) per Appendix I; Table 1 reports the cohort median as 90 mg/dL, but the fit was
    # carried out with AGPm = 94 -- value 0.94 g/L is therefore the source-of-truth here.
    slope <- exp(lslope + etalslope) * (AAG / 0.94)^e_aag_slope

    # KTR = (n + 1) / MTT with n = 3 transits (Appendix I `KTR = 4/MTT`; Methods text)
    ktr <- 4 / mtt
    # KIN = 1 / MIT (Methods: MIT = 1/kin)
    kin <- 1 / mit

    # Concentration in central compartment (mg/L, since dose is in mg and Vc in L)
    Cc <- central / vc

    # Linear drug effect (Equation 7)
    edrug <- slope * Cc

    # ---- ODE system ----
    # Three-compartment IV PK (Appendix I $DES DADT(1..3); compartment order matches NONMEM
    # CENTRAL=1, PERIPH1=2, PERIPH2=3 so dose records may use cmt = "central" or cmt = 1).
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Dexamethasone-induced ANC bump (Equation 5).
    # `input` is a paper-specific (non-canonical) compartment name retained from Ozawa 2007's
    # nomenclature; it represents the zero-order ANC increase attributed to dexamethasone
    # premedication that feeds the circulating compartment with rate constant kin.
    d/dt(input) <- -kin * input

    # Friberg myelosuppression chain with feedback and dex-induced input (Equations 1-4, 6).
    # Compartment naming follows the Friberg-family convention used by Friberg_2002_paclitaxel
    # and Netterberg_2017_docetaxel: precursor1 = Prol (proliferating pool), precursor2..4 =
    # Transit1..3 (maturation chain), circ = Circ (circulating). The feedback uses
    # (Circ0 / Circ)^gamma1 with Circ0 = NEUT (per-subject baseline ANC, supplied as the
    # `NEUT` covariate). Self-renewing proliferation rate kprol = ktr per the steady-state
    # assumption (Methods); clearance rate from the circulating pool kcirc = ktr per the
    # parameter-reduction assumption (Methods).
    d/dt(precursor1) <- ktr * precursor1 * (1 - edrug) * (NEUT / circ)^gamma1 - ktr * precursor1
    d/dt(precursor2) <- ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <- ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr * precursor3 - ktr * precursor4
    d/dt(circ)       <- ktr * precursor4 + kin * input - ktr * circ

    # Initial conditions:
    # - precursor1 = precursor2..4 = circ initialised at the per-subject baseline ANC (NEUT
    #   covariate), per the Methods assumption that the system is at steady state at t = 0
    #   with no drug.
    # - Input(0) = IP0 (Methods text 'Input (t = 0) = IP0'; Appendix I $PK 'F4 = IP0' combined
    #   with TIME=0 AMT=1 record on CMT=4 is the NONMEM-style way of producing the same
    #   effect; in rxode2 we set the state directly).
    precursor1(0) <- NEUT
    precursor2(0) <- NEUT
    precursor3(0) <- NEUT
    precursor4(0) <- NEUT
    circ(0)       <- NEUT
    input(0)      <- ip0

    # ---- Observations ----
    Cc  ~ prop(propSd)
    ANC <- circ
    ANC ~ lnorm(addSd_ANC)
  })
}

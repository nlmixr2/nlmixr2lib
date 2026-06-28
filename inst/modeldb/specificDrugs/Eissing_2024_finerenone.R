Eissing_2024_finerenone <- function() {
  description <- "Two-compartment population PK model with a 4-transit-compartment delayed first-order absorption for finerenone in adults with chronic kidney disease and type 2 diabetes (FIGARO-DKD final PK model; Eissing 2024)"
  reference <- paste(
    "Eissing T, Goulooze SC, van den Berg P, et al. Pharmacokinetics and",
    "pharmacodynamics of finerenone in patients with chronic kidney",
    "disease and type 2 diabetes: Insights based on FIGARO-DKD and",
    "FIDELIO-DKD. Diabetes Obes Metab. 2024;26(3):924-936.",
    "doi:10.1111/dom.15387",
    sep = " "
  )
  vignette <- "Eissing_2024_finerenone"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on Vc/F with reference 87 kg (median in FIGARO-DKD PK dataset, per supplement NONMEM code header).",
      source_name        = "BW0"
    ),
    CRCL = list(
      description        = "CKD-EPI estimated glomerular filtration rate (eGFR)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F and inversely on F1 with reference 60.6 mL/min/1.73 m^2 (FIGARO-DKD median). Derived per CKD Epidemiology Collaboration formula (eGFR-EPI). Source column EGFREP in the NONMEM supplement.",
      source_name        = "EGFREP"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F only (does NOT enter F1) with reference 73 U/L (FIGARO-DKD median). Replaces serum creatinine and gamma-glutamyl transferase as the retained hepatic-secretion / liver-function descriptor in the FIGARO-DKD final PK model.",
      source_name        = "ALP"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F and inversely on F1 with reference 1.97 m^2 (FIGARO-DKD median). Selected in FIGARO-DKD in preference to body height (which was the analogous descriptor retained in FIDELIO-DKD).",
      source_name        = "BSA"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F and inversely on F1 with reference 71 g/L (= 7.1 g/dL FIGARO-DKD median in the supplement code). The supplement reports the reference and coefficient in g/dL; the model converts the canonical TPRO g/L input to g/dL inline via tpro_gdL <- TPRO * 0.1 so the structural coefficient remains aligned with its original calibration (per inst/references/covariate-columns.md TPRO Section g/L canonical convention).",
      source_name        = "PROT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3102,
    n_studies      = 1,
    age_range      = "median (5th-95th percentile) not separately tabulated for the PK subpopulation; FIGARO-DKD overall median age 64 years",
    age_median     = "65 years (FIGARO-DKD population)",
    weight_range   = "5th-95th percentile 58.0-123.6 kg",
    weight_median  = "87 kg (FIGARO-DKD PK dataset; reference in NONMEM code header)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Adults with chronic kidney disease (CKD) and type 2 diabetes (T2D); eGFR 25-90 mL/min/1.73 m^2, urine albumin-creatinine ratio 30-5000 mg/g, and serum potassium <= 4.8 mmol/L at screening. Median (5th-95th percentile) eGFR-EPI was 67.6 (33.3-103) mL/min/1.73 m^2 and median UACR was 308 (36.8-2115) mg/g.",
    dose_range     = "10 or 20 mg oral once daily, with potassium- and eGFR-guided dose titration during the trial (starting 10 mg if baseline eGFR < 60 mL/min/1.73 m^2 and 20 mg if baseline eGFR >= 60); average dose level over time was 17.5 mg in FIGARO-DKD.",
    regions        = "Global multinational phase 3 trial (NCT02540993; FIGARO-DKD).",
    notes          = "PK dataset: 8142 valid finerenone plasma concentrations from 3102 patients (Results, Pharmacokinetics, paragraph 1). Median exposure parameters were 10% (AUC) and 7% (Cmax) lower than the more advanced-disease FIDELIO-DKD cohort. The full FIDELITY pool (FIGARO-DKD + FIDELIO-DKD) randomized 13026 participants. This model file extracts the FIGARO-DKD final population PK model only; the paper's downstream serum-potassium, UACR, eGFR, and parametric time-to-event PK/PD models use individual posthoc CL/F and F1 from this PK fit as fixed inputs and are out of scope for the compact nlmixr2lib fittable popPK/PD format."
  )

  ini({
    # Reference participant (FIGARO-DKD PK dataset typical):
    #   BW = 87 kg, eGFR-EPI = 60.6 mL/min/1.73 m^2, ALP = 73 U/L,
    #   BSA = 1.97 m^2, total serum protein = 7.1 g/dL.
    # Parameter values: supplement Item S2, FIGARO-DKD population PK
    # model NONMEM code, $THETA block (final estimates).

    lka   <- log(20.0);  label("Transit-compartment absorption / chain rate constant (1/h)")   # supplement $THETA TH2 (4 sequential transit steps DEPOT->BUFFER->BUFFER2->BUFFER3->CENTRAL with shared rate Ka)
    lcl   <- log(35.1);  label("Apparent clearance, reference participant (CL/F, L/h)")        # supplement $THETA TH3
    lvc   <- log(126);   label("Apparent central volume of distribution, reference participant (Vc/F, L)")  # supplement $THETA TH4
    lq    <- log(0.357); label("Apparent inter-compartmental clearance (Q/F, L/h)")             # supplement $THETA TH5

    # Structural fixed parameters
    ltlag    <- fixed(log(0.215)); label("Absorption lag time (h)")           # supplement $THETA TH7 (FIX)
    lfdepot  <- fixed(log(1));     label("Relative bioavailability anchor (unitless)")  # supplement $THETA TH8 (F1 = 1 FIX); covariate-driven correction applied in model()
    # Vp/F : Vc/F ratio is fixed to 1 in the FIGARO-DKD final model
    # (supplement $THETA TH6 = 1 FIX); encoded as vp <- vc inline.

    # Covariate effects (power-form exponents)
    e_crcl_cl <-  0.124;   label("Power exponent of eGFR-EPI on CL/F and inverse on F1 (unitless)")  # supplement $THETA TH10
    e_alp_cl  <- -0.0954;  label("Power exponent of alkaline phosphatase on CL/F only (unitless)")    # supplement $THETA TH11
    e_bsa_cl  <-  0.492;   label("Power exponent of BSA on CL/F and inverse on F1 (unitless)")        # supplement $THETA TH12
    e_tpro_cl <- -0.246;   label("Power exponent of total serum protein on CL/F and inverse on F1 (unitless)")  # supplement $THETA TH13
    e_wt_vc   <-  0.642;   label("Power exponent of body weight on Vc/F (unitless)")                  # supplement $THETA TH9

    # IIV: 2x2 omega block on lcl and lvc (variance on log scale).
    #   var(etalcl) = 0.0985 ; var(etalvc) = 0.0925 ; cov = 0.0330
    # corresponds to CV% ~ 32.0% (CL/F), 31.0% (Vc/F), correlation 0.35.
    etalcl + etalvc ~ c(
      0.0985,
      0.0330,  0.0925
    )                                                                                              # supplement $OMEGA BLOCK(2)

    # Residual error: NONMEM proportional with W = SQRT(THETA(1)) * IPRED
    # and $SIGMA 1 FIX gives propSd = sqrt(0.05) ~ 0.2236.
    propSd <- 0.2236; label("Proportional residual error (fraction)")                              # supplement $THETA TH1 (= 0.05, on log-CV scale)
  })
  model({
    # Reference values (FIGARO-DKD PK dataset typical, supplement NONMEM code header)
    ref_wt   <- 87
    ref_crcl <- 60.6
    ref_alp  <- 73
    ref_bsa  <- 1.97
    ref_tpro <- 7.1   # g/dL (paper-calibrated reference)

    # Convert canonical g/L TPRO input to the g/dL scale at which the
    # supplement coefficient was calibrated (1 g/dL = 10 g/L).
    tpro_gdL <- TPRO * 0.1

    # Covariate factor shared between CL/F and the inverse F1 correction
    # (eGFR, BSA, TPRO): in the NONMEM supplement these enter CL/F as
    # CV2 * CV4 * CV5 and F1 as 1 / (CV2 * CV4 * CV5).
    covfac_cl_f <- (CRCL     / ref_crcl)^e_crcl_cl *
                   (BSA      / ref_bsa )^e_bsa_cl  *
                   (tpro_gdL / ref_tpro)^e_tpro_cl

    # Covariate factor that enters CL/F only (ALP); supplement CV3.
    covfac_cl_only <- (ALP / ref_alp)^e_alp_cl

    cl <- exp(lcl + etalcl) * covfac_cl_f * covfac_cl_only
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    vp <- vc                   # supplement: V3 = THETA(6) * V2 with THETA(6) = 1 FIX
    q  <- exp(lq)
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Four sequential transit compartments share the same rate Ka
    # (supplement: K14 = K45 = K56 = K62 = THETA(2)). The chain is
    #   depot -> transit1 -> transit2 -> transit3 -> central.
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot     - ka  * transit1
    d/dt(transit2)    <-  ka * transit1  - ka  * transit2
    d/dt(transit3)    <-  ka * transit2  - ka  * transit3
    d/dt(central)     <-  ka * transit3  - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Bioavailability and lag time on the depot compartment
    # (supplement: ALAG1 and F1 on compartment 1 = DEPOT).
    f(depot)    <- exp(lfdepot) / covfac_cl_f
    alag(depot) <- exp(ltlag)

    # Concentration: dose in mg, Vc in L -> mg/L; multiply by 1000 to
    # express finerenone plasma concentration in ug/L (= ng/mL), the
    # unit used in the source paper's PK figures and labelling.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

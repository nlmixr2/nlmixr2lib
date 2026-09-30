Trang_2021_vaborbactam <- function() {
  description <- "Two-compartment population PK model for intravenous vaborbactam (given with meropenem) in noninfected adults and adults with complicated urinary tract, bloodstream or other serious infections, with a sigmoidal Hill relationship between renal clearance and MDRD eGFR, power effects of height on clearance and body surface area on central volume, proportional study-phase shifts on clearance and both volumes, and a cumulative urine compartment for the urinary concentrations"
  reference <- paste(
    "Trang M, Griffith DC, Bhavnani SM, Loutit JS, Dudley MN, Ambrose PG,",
    "Rubino CM. 2021. Population pharmacokinetics of meropenem and",
    "vaborbactam based on data from noninfected subjects and infected",
    "patients. Antimicrob Agents Chemother 65:e02606-20.",
    "doi:10.1128/AAC.02606-20. Covariate-equation forms and reference",
    "values from the FDA Clinical Pharmacology review of NDA 209776",
    "(Vabomere, 2017), section 4.2 (vaborbactam covariate equations).",
    sep = " "
  )
  vignette <- "Trang_2021_meropenem_vaborbactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vaborbactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vaborbactam", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "vaborbactam", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste0(
        "Estimated glomerular filtration rate by the Modification of Diet in ",
        "Renal Disease equation, BSA-normalised (mL/min/1.73 m^2); time-varying"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Materials and Methods, Demographics: MDRD eGFR recalculated at every ",
        "serum creatinine measurement and treated as time-varying, with serum ",
        "creatinine capped at a lower bound of 0.5 mg/dL. Drives the renal ",
        "clearance arm through the sigmoidal Hill term (no normalising value). ",
        "Pooled median 90.1, range 4.50-338 mL/min/1.73 m^2 (Table 1)."
      ),
      source_name = "eGFR"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters total CL as (HT / 168)^2.24 (Table 3). The 168 cm reference is ",
        "printed in the FDA review vaborbactam CL equation and equals the ",
        "Table 1 pooled median. Range 145-193 cm (Table 1)."
      ),
      source_name = "HTCM"
    ),
    BSA = list(
      description = "Body surface area by the DuBois and DuBois method",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters Vc as (BSA / 1.88)^1.50 (Table 3). The 1.88 m^2 reference is ",
        "printed in the FDA review vaborbactam Vc equation. Pooled median 1.84, ",
        "range 1.27-2.83 m^2 (Table 1)."
      ),
      source_name = "BSA"
    ),
    STUDY_PHASE3 = list(
      description = "Phase 3 study (infected patient) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (phase 3 infected patient, studies 505 and 506)",
      notes = paste0(
        "The source indicator 'Phase' is 1 for phase 1 noninfected subjects ",
        "(studies 501 and 504, including the renally impaired subjects of 504) ",
        "and 0 for phase 3 infected patients (FDA review: 'Phase is an ",
        "indicator variable with values equal to 0 for Phase III patients and ",
        "1 for Phase I subjects'), so Phase = 1 - STUDY_PHASE3. The phase 3 ",
        "patients are the typical-value reference; the proportional shifts ",
        "(1 + theta * Phase) raise CL and Vp and lower Vc in phase 1 subjects."
      ),
      source_name = "Phase"
    ),
    URINE_VOL_INTERVAL = list(
      description = "Urine volume collected in the urine collection interval containing the current urinary observation",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Observation-equation divisor only (Curine = urine amount / collected ",
        "volume). Urine was collected in the phase 1 studies over 0-4, 4-8, ",
        "8-12, 12-24, 24-48 and 48-72 h intervals (supplement Table S3), so ",
        "the urine state must be reset to zero at each interval boundary ",
        "(evid = 5, amt = 0, cmt = 'urine'). Not needed for plasma-only ",
        "simulation; any positive value may then be supplied."
      ),
      source_name = "urine volume"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 414L,
    n_studies = 4L,
    age_range = "18-92 years",
    age_median = "53 years",
    weight_range = "40.0-177 kg",
    weight_median = "75.0 kg",
    sex_female_pct = 55.5,
    race_ethnicity = "Not tabulated; race was screened and not retained (Figure S7).",
    disease_state = paste(
      "Pooled phase 1 noninfected adults (study 501, healthy volunteers; study",
      "504, normal renal function to end-stage renal disease) and phase 3",
      "infected patients (study 505 TANGO I, complicated urinary tract",
      "infection or acute pyelonephritis; study 506 TANGO II, infections due",
      "to confirmed or suspected carbapenem-resistant Enterobacterales,",
      "including bacteremia, HABP/VABP, cIAI and cUTI)."
    ),
    dose_range = paste(
      "Vaborbactam 0.25-2 g (0.5 g in severe renal impairment in study 506) as",
      "3-h intravenous infusions (1 h in study 501 group 6), single dose or",
      "every 8 h, co-administered with meropenem; renal dose adjustments per",
      "supplement Table S3."
    ),
    regions = "Multinational (phase 3 TANGO I and TANGO II).",
    renal_function = "eGFR (MDRD) 4.50-338 mL/min/1.73 m^2; median 90.1.",
    notes = paste(
      "Final data set: 4,082 plasma concentrations from 93 noninfected subjects",
      "and 321 infected patients plus 746 urine concentrations from 75",
      "noninfected subjects (Results). Table 1 demographics are for all 431",
      "enrolled subjects/patients (239 female, 55%). Assay LC-MS/MS, 0.02-100",
      "mg/L; BLQ handled with Beal M3. NONMEM 7.2, FOCE-I."
    )
  )

  ini({
    # Structural parameters -- Trang 2021 Table 3 (final vaborbactam model),
    # for a phase 3 infected patient of height 168 cm and BSA 1.88 m^2.
    # CL = (CLNR + CLR,max * eGFR^h / (eGFR50^h + eGFR^h))
    #      * (1 + theta_CL * Phase) * (HTCM/168)^theta_HT
    # (functional form from the FDA review; see covariateData).
    lcl_nonren <- log(0.157)
    label("Nonrenal clearance, phase 3 patient, 168 cm (L/h)") # Table 3, 'CL NR' = 0.157 L/h (%SEM 12.7)
    lcl_renal_max <- log(8.86)
    label("Maximal renal clearance, phase 3 patient, 168 cm (L/h)") # Table 3, 'CL R,max' = 8.86 L/h (%SEM 3.50)
    lcrcl50 <- log(49.7)
    label("eGFR at half-maximal renal clearance (mL/min/1.73 m^2)") # Table 3, 'eGFR50' = 49.7 (%SEM 3.10)
    lhill <- log(2.25)
    label("Hill coefficient of the renal clearance-eGFR relationship (unitless)") # Table 3, 'Hill coefficient' = 2.25 (%SEM 3.40)
    lvc <- log(17.1)
    label("Central volume of distribution, phase 3 patient, BSA 1.88 m^2 (L)") # Table 3, 'V c' = 17.1 L (%SEM 2.90)
    lq <- log(2.75)
    label("Distributional clearance (L/h)") # Table 3, 'CL d' = 2.75 L/h (%SEM 10.5)
    lvp <- log(1.77)
    label("Peripheral volume of distribution, phase 3 patient (L)") # Table 3, 'V p' = 1.77 L (%SEM 11.0)

    # Covariate effects.
    e_ht_cl <- 2.24
    label("Power exponent of height on CL (unitless)") # Table 3, 'Power coefficient of HTCM on CL' = 2.24 (%SEM 22.3)
    e_phase1_cl <- 0.517
    label("Proportional shift in CL for phase 1 noninfected subjects (fraction)") # Table 3, 'Proportional shift with Phase on CL' = 0.517 (%SEM 30.8)
    e_bsa_vc <- 1.50
    label("Power exponent of body surface area on Vc (unitless)") # Table 3, 'Power coefficient of BSA on V c' = 1.50 (%SEM 13.0)
    e_phase1_vc <- -0.215
    label("Proportional shift in Vc for phase 1 noninfected subjects (fraction)") # Table 3, 'Proportional shift with Phase on V c' = -0.215 (%SEM 42.0)
    e_phase1_vp <- 1.28
    label("Proportional shift in Vp for phase 1 noninfected subjects (fraction)") # Table 3, 'Proportional shift with Phase on V p' = 1.28 (%SEM 22.3)

    # IIV. Table 3 reports %CV; omega^2 = log(CV^2 + 1). The final model
    # fitted a full covariance matrix whose off-diagonal elements are not
    # published, so the etas are entered as independent (vignette errata).
    etalcl ~ 0.18891 # Table 3, IIV 45.6 %CV (printed on the 'CL R,max' row; applied to total CL, see vignette)
    etalvc ~ 0.14430 # Table 3, 'V c' IIV 39.4 %CV
    etalq ~ 0.11246 # Table 3, 'CL d' IIV 34.5 %CV
    etalvp ~ 0.05155 # Table 3, 'V p' IIV 23.0 %CV

    # Residual error. Table 3 prints unlabelled sigma values; they are
    # variances (sigma^2), as for the meropenem model: read as SDs the plasma
    # proportional error would be 3.7%, inconsistent with the +/-2 IWRES
    # spread and observed-versus-IPRED scatter of supplement Figure S2.
    propSd <- sqrt(0.0372)
    label("Plasma proportional residual SD (fraction)") # Table 3, 'Plasma proportional error' = 0.0372 (variance); SD 0.1929
    addSd <- sqrt(0.0287)
    label("Plasma additive residual SD (mg/L)") # Table 3, 'Plasma additive error' = 0.0287 (variance); SD 0.1694 mg/L
    propSd_Curine <- sqrt(0.115)
    label("Urine proportional residual SD (fraction)") # Table 3, 'Urine proportional error' = 0.115 (variance); SD 0.3391
    addSd_Curine <- sqrt(5.46)
    label("Urine additive residual SD (mg/L)") # Table 3, 'Urine additive error' = 5.46 (variance); SD 2.337 mg/L
  })

  model({
    # 1. Derived covariate terms. Source 'Phase' indicator = 1 for phase 1
    #    noninfected subjects.
    phase1 <- 1 - STUDY_PHASE3
    hill <- exp(lhill)
    crcl50 <- exp(lcrcl50)

    # Nonrenal and renal clearance arms.
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal_max) * CRCL^hill / (crcl50^hill + CRCL^hill)
    frac_renal <- cl_renal / (cl_nonren + cl_renal)

    # 2. Individual PK parameters. Phase and height scale total CL; the
    #    renal fraction is carried unchanged into the urine arm.
    cl <- (cl_nonren + cl_renal) * (1 + e_phase1_cl * phase1) *
      (HT / 168)^e_ht_cl * exp(etalcl)
    vc <- exp(lvc + etalvc) * (1 + e_phase1_vc * phase1) * (BSA / 1.88)^e_bsa_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (1 + e_phase1_vp * phase1)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Zero-order intravenous infusion into central via the
    #    event-table rate; renally cleared drug accumulates in urine.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(urine) <- frac_renal * kel * central

    # 5. Observations. Urine concentration = amount excreted in the
    #    collection interval / collected volume (mL converted to L; guarded
    #    against a zero volume).
    Cc <- central / vc
    urine_volume <- max(URINE_VOL_INTERVAL / 1000, 0.001)
    Curine <- urine / urine_volume

    Cc ~ add(addSd) + prop(propSd)
    Curine ~ add(addSd_Curine) + prop(propSd_Curine)
  })
}

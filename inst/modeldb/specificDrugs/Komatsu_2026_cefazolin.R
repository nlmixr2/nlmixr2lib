Komatsu_2026_cefazolin <- function() {
  description <- "Two-compartment population PK model for intravenous cefazolin surgical prophylaxis in men undergoing robotic-assisted radical prostatectomy (Komatsu 2026, n = 67). The model is fitted to UNBOUND serum concentrations, so Cc = central / vc is the unbound concentration and CL, Vc, Q and Vp are unbound-referenced. Periprostatic adipose tissue concentration is linked to the central compartment by a rapid-equilibrium partition coefficient (Cadipose = kp_adipose * Cc), with no transfer rate, because each patient contributed a single tissue sample. Unbound clearance scales with Cockcroft-Gault creatinine clearance as a power function centred at the cohort median of 69 mL/min. Log-normal IIV on CL, Vc, Vp and the partition coefficient (IIV on Q fixed to zero); proportional residual error. Estimated in NONMEM 7.6 by FOCE-I."
  reference <- "Komatsu T, Takahashi Y, Takayama Y, Akamada Y, Ikeda M, Tsumura H, Ishii D, Matsumoto K, Iwamura M, Okamoto H, Hanaki H, Otori K. Unbound and Periprostatic Adipose Tissue Cefazolin Pharmacokinetics in Robotic-Assisted Radical Prostatectomy. Antibiotics (Basel). 2026;15(2):181. doi:10.3390/antibiotics15020181. PMCID PMC12937387. All parameter estimates from Table 2; model equation from Results Section 2; tissue-link and covariate forms from Methods Section 4.3."
  vignette <- "Komatsu_2026_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "cefazolin",
      units = "mg",
      specimen = "serum",
      verified = TRUE,
      notes = "Holds the administered cefazolin amount. Dividing by the unbound-referenced central volume Vc gives the UNBOUND serum concentration (ultrafiltrate through a 10 kDa Nanosep device, Methods Section 4.2). No total concentration was modelled. The unbound-referenced Vc of 44.3 L and CL of 25.8 L/h are several-fold larger than total-referenced cefazolin values for the same reason (cefazolin is about 80-90% protein bound, Introduction)."
    ),
    peripheral1 = list(
      analyte = "cefazolin",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "Unbound-referenced peripheral distribution compartment (Vp = 52.2 L). It is NOT the periprostatic adipose tissue: the tissue concentration is an algebraic function of the central concentration."
    )
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation, raw mL/min and NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 69 mL/min, the cohort median (Table 1, IQR 60-88), printed as the denominator of the Results Section 2 equation CL = 25.8 * (CLcr/69)^0.735. Cockcroft-Gault per Methods Section 4.3; the paper's Monte Carlo simulations (Figure 4) extrapolate the term across 5-120 mL/min.",
      source_name = "CLcr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods Section 4.3) and not retained. Cohort median 69 years (IQR 63-73)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods Section 4.3) and not retained; weight still set the dose (1 g below 80 kg, 2 g at or above 80 kg). Cohort median 67.2 kg (IQR 60.4-74.2)."
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened (Methods Section 4.3) and not retained. Cohort median 17.1 umol/L (IQR 13.7-23.1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Methods Section 4.3) and not retained, even though cefazolin binding is albumin-dependent; the model is fitted to unbound concentrations only. Cohort median 32 g/L (IQR 30-34)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods Section 4.3) and not retained separately from CRCL. Cohort median 0.92 mg/dL (IQR 0.83-1.06)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 67,
    n_studies = 1,
    age_range = "median 69 years (IQR 63-73)",
    weight_range = "median 67.2 kg (IQR 60.4-74.2)",
    sex_female_pct = 0,
    race_ethnicity = "Japanese (single-centre study at Kitasato University Hospital, Kanagawa, Japan)",
    disease_state = "Men with localized prostate cancer undergoing robotic-assisted radical prostatectomy (RARP), receiving cefazolin surgical prophylaxis. Excluded: known cephalosporin hypersensitivity.",
    dose_range = "1 g (weight < 80 kg, n = 59) or 2 g (weight >= 80 kg, n = 8) cefazolin as a 10-minute IV infusion within 60 min before incision, redosed every 3 h intraoperatively at the same dose.",
    renal_function = "Cockcroft-Gault creatinine clearance median 69 mL/min (IQR 60-88); serum creatinine median 0.92 mg/dL (IQR 0.83-1.06).",
    regions = "Japan",
    notes = "Prospective cohort, December 2021 to February 2025. 295 unbound serum samples (at incision, during prostate removal, at redosing, before skin closure) and 67 periprostatic adipose tissue samples (one per patient, at prostate removal; tissue density assumed 1 kg/L). Albumin median 32 g/L (IQR 30-34); total bilirubin median 17.1 umol/L (IQR 13.7-23.1). Table 1."
  )

  ini({
    lcl <- log(25.8); label("Unbound clearance at CRCL = 69 mL/min (L/h)") # Table 2 theta1 = 25.8 (RSE 5.47%)
    lvc <- log(44.3); label("Unbound central volume of distribution (L)") # Table 2 theta3 = 44.3 (RSE 10.18%)
    lq <- log(44.0); label("Unbound intercompartmental clearance (L/h)") # Table 2 theta4 = 44.0 (RSE 11.84%)
    lvp <- log(52.2); label("Unbound peripheral volume of distribution (L)") # Table 2 theta5 = 52.2 (RSE 8.85%)
    lkp_adipose <- log(0.638); label("Periprostatic adipose tissue to unbound serum partition coefficient FPA (unitless)") # Table 2 theta6 = 0.638 (RSE 6.68%)
    e_crcl_cl <- 0.735; label("Power exponent on (CRCL/69) for unbound clearance (unitless)") # Table 2 theta2 = 0.735 (RSE 19.46%)

    # Table 2 reports omega as the square root of the inter-individual
    # variance in percent (table footnote), so the variance is its square.
    # IIV on Q was fixed to zero in the final model (Results Section 2).
    etalcl ~ 0.283^2 # Table 2 omega CL = 28.3%
    etalvc ~ 0.529^2 # Table 2 omega Vc = 52.9%
    etalvp ~ 0.463^2 # Table 2 omega Vp = 46.3%
    etalkp_adipose ~ 0.45^2 # Table 2 omega FPA = 45.0%

    # Table 2 lists a single proportional residual error for the joint
    # serum + tissue fit; it is carried for both outputs.
    propSd <- 0.216; label("Proportional residual SD, unbound serum Cc (fraction)") # Table 2 epsilon pro = 21.6% (RSE 8.79%)
    propSd_Cadipose <- 0.216; label("Proportional residual SD, periprostatic adipose Cadipose (fraction)") # Table 2 epsilon pro = 21.6% (RSE 8.79%), the single reported residual term
  })

  model({
    cl <- exp(lcl + etalcl) * (CRCL / 69)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    kp_adipose <- exp(lkp_adipose + etalkp_adipose)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Unbound serum concentration; tissue in rapid equilibrium with it
    # (Methods Section 4.3), mg/kg read as mg/L at a density of 1 kg/L.
    Cc <- central / vc
    Cadipose <- kp_adipose * Cc

    Cc ~ prop(propSd)
    Cadipose ~ prop(propSd_Cadipose)
  })
}

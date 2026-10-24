Rose_2016b_follitropin_delta <- function() {
  description <- paste(
    "One-compartment population PK model for FE 999049 (recombinant human",
    "FSH; INN follitropin delta) after repeated daily subcutaneous dosing in",
    "24 gonadotropin-downregulated healthy women, with first-order absorption",
    "through a single transit compartment and first-order elimination. The",
    "central compartment holds TOTAL serum FSH: exogenous FE 999049 plus",
    "endogenous FSH from a zero-order turnover input that is suppressed over",
    "time by observed serum inhibin B through an Imax function. The",
    "endogenous FSH baseline decreases with baseline progesterone (power",
    "function), and body weight scales CL/F and V/F allometrically",
    "(exponents fixed at 0.75 and 1, reference 65 kg)."
  )
  reference <- paste(
    "Rose TH, Roshammar D, Erichsen L, Grundemar L, Ottesen JT.",
    "Characterisation of Population Pharmacokinetics and Endogenous",
    "Follicle-Stimulating Hormone (FSH) Levels After Multiple Dosing of a",
    "Recombinant Human FSH (FE 999049) in Healthy Women.",
    "Drugs R D. 2016;16(2):165-172. doi:10.1007/s40268-016-0126-z.",
    sep = " "
  )
  vignette <- "Rose_2016b_follitropin_delta"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "follitropin delta", units = "ug", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "follitropin delta", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "FSH (endogenous FSH + follitropin delta)",
      units = "ug",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1),",
        "both exponents fixed, reference weight 65 kg (Section 3 and Table 2",
        "footnote). Cohort mean 71.5 kg, range 46.1-86.6 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    PROG_BASE = list(
      description = "Observed baseline (pre-dose) serum progesterone",
      units = "ug/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the endogenous FSH baseline as the un-normalised power",
        "function rbase = FSHbl * PROG_BASE^Progblef, i.e. the reference is",
        "1 ug/L: Figure 2's typical-prediction curve passes through",
        "FSHbl = 0.162 ug/L at progesterone 1.0 ug/L. Cohort mean 0.652",
        "ug/L, range 0.08-1.49 ug/L (Table 1). Subjects with a missing",
        "baseline were given the median population baseline (Section 2.1)."
      ),
      source_name = "Prog"
    ),
    INHIBINB = list(
      description = "Observed serum inhibin B, time-varying",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Suppresses the zero-order endogenous FSH input through the",
        "fractional Imax term 1 - INHIBINB / (INHIBINB + ec50_inhibinb),",
        "complete inhibition possible (Eq. 8). Measured once daily; the",
        "paper carried the last measured value forward between samples",
        "(Section 2.1), so solve with covsInterpolation = 'locf'.",
        "Baseline cohort mean 22.0 pg/mL, range 4.0-89.0 pg/mL (Table 1);",
        "rose to the 1000-6000 pg/mL range during dosing (Figure 5)."
      ),
      source_name = "InhB(t)"
    ),
    INHIBINB_BASE = list(
      description = "Observed baseline (pre-dose) serum inhibin B, time-fixed",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The t = 0 value of the time-varying INHIBINB column. Needed because",
        "the endogenous production rate kin is back-calculated per subject",
        "from the pre-dose steady state (Section 2.2), which requires the",
        "inhibin B value at t = 0 at every later time. Cohort mean 22.0",
        "pg/mL, range 4.0-89.0 pg/mL (Table 1)."
      ),
      source_name = "InhB(0)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "21.5-38.7 years (mean 31.1; Table 1)",
    weight_range = "46.1-86.6 kg (mean 71.5; Table 1)",
    sex_female_pct = 100,
    disease_state = paste(
      "Healthy women downregulated with two doses of a 1-month GnRH-agonist",
      "depot (LUPRON DEPOT) on days -28 and -14. Downregulation was",
      "incomplete: baseline serum FSH 0.211 ug/L (range 0.089-0.376),",
      "baseline progesterone 0.652 ug/L, baseline inhibin B 22.0 pg/mL",
      "(Table 1)."
    ),
    dose_range = "225 IU FE 999049 subcutaneously once daily for 7 days",
    regions = "Not reported",
    notes = paste(
      "FE 999049 arm (24 of 49 women) of a randomised, double-blind,",
      "active-control multiple-dose trial versus GONAL-F (Olsson 2014).",
      "672 serum FSH samples (MSD electrochemiluminescence immunoassay,",
      "LLOQ 0.075 ug/L); 3 (0.4%) were below the LLOQ and excluded. FSH was",
      "sampled pre-dose, daily for 15 days, and every 4 h from the last dose",
      "on day 6 to day 8. NONMEM 7.2.0."
    )
  )

  ini({
    # Structural parameters (Table 2). CL/F and V/F are typical values for a
    # 65 kg woman; F is not identified from the SC-only data.
    lcl <- log(0.423); label("Apparent clearance CL/F at WT = 65 kg (L/h)") # Table 2: CL/F 0.423 L/h (RSE 3.9%)
    lvc <- log(24.3); label("Apparent volume of distribution V/F at WT = 65 kg (L)") # Table 2: V/F 24.30 L (RSE 4.6%)
    lktr <- log(0.329); label("Absorption rate constant from the dosing site to the transit compartment (1/h)") # Table 2: ktr 0.329 1/h (RSE 17.0%)
    lka <- log(0.148); label("Absorption rate constant from the transit compartment to central (1/h)") # Table 2: ka 0.148 1/h (RSE 13.2%)

    # Endogenous FSH (Table 2). FSHbl is the endogenous FSH baseline
    # concentration at progesterone = 1 ug/L; Progblef is the power exponent
    # of baseline progesterone on FSHbl; InhBef is the inhibin B level giving
    # 50% inhibition of endogenous FSH production.
    lrbase <- log(0.162); label("Endogenous FSH baseline concentration FSHbl at PROG_BASE = 1 ug/L (ug/L)") # Table 2: FSHbl 0.162 ug/L (RSE 9.1%)
    e_prog_base_rbase <- -0.246; label("Power exponent of baseline progesterone on FSHbl (unitless)") # Table 2: Progblef -0.246 (RSE 34.7%)
    lec50_inhibinb <- log(100); label("Inhibin B level giving half-maximal inhibition of endogenous FSH production (pg/mL)") # Table 2: InhBef 100 (RSE 37.2%)

    # Allometric exponents, held at the canonical allometric values (Section 3).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL/F (unitless)") # Section 3: power exponent fixed to 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of WT on V/F (unitless)") # Section 3: power exponent fixed to 1

    # IIV, Table 2 'IIV CV%', converted with omega^2 = log(CV^2 + 1). The
    # CL/F - V/F covariance is reported as positive and significant (Table 3,
    # dOFV 11.77) but its magnitude is not printed, so the block is diagonal.
    etalcl ~ 0.024045 # Table 2: CL/F IIV 15.6% CV
    etalvc ~ 0.033296 # Table 2: V/F IIV 18.4% CV
    etalktr ~ 0.52801 # Table 2: ktr IIV 83.4% CV
    etalrbase ~ 0.074443 # Table 2: FSHbl IIV 27.8% CV

    # Combined additive + proportional residual error (Section 3); neither
    # magnitude is reported in Table 2 or the text.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)") # Section 3: combined error model, no estimate published
    addSd <- fixed(0); label("Additive residual SD (ug/L; 0 -- not reported in the source)") # Section 3: combined error model, no estimate published
  })
  model({
    # Individual parameters; allometric WT scaling with reference 65 kg.
    cl <- exp(lcl + etalcl) * (WT / 65)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 65)^e_wt_vc
    ktr <- exp(lktr + etalktr)
    ka <- exp(lka)
    kel <- cl / vc

    # Endogenous FSH baseline: un-normalised power function of baseline
    # progesterone in ug/L (Figure 2 curve passes FSHbl at Prog = 1 ug/L).
    rbase <- exp(lrbase + etalrbase) * PROG_BASE^e_prog_base_rbase
    ec50_inhibinb <- exp(lec50_inhibinb)

    # Zero-order endogenous FSH input back-calculated from the pre-dose
    # steady state of Eq. 8 (Section 2.2): at t = 0 the central amount is
    # rbase * vc and d/dt(central) = 0, so
    #   kin * ec50 / (InhB(0) + ec50) = kel * rbase * vc.
    kin <- kel * rbase * vc * (INHIBINB_BASE + ec50_inhibinb) / ec50_inhibinb
    inh_inhibinb <- INHIBINB / (INHIBINB + ec50_inhibinb)

    # Eqs. 6-8: rhFSH_DS -> depot, rhFSH_TR -> transit1, total FSH -> central.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- kin * (1 - inh_inhibinb) + ka * transit1 - kel * central
    central(0) <- rbase * vc

    # Total (endogenous + exogenous) serum FSH: ug / L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

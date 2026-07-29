Willmann_2018_rivaroxaban <- function() {
  description <- "Paediatric population PK model for oral rivaroxaban in children aged 0.5-18 years (Willmann 2018, EINSTEIN-Jr phase I). Linear two-compartment model with first-order absorption from a depot and first-order elimination from the central compartment; CL and central V allometrically scaled to body weight (CL exponent 0.323 estimated; V exponent 1 fixed); ka shifted between the undiluted-oral-suspension formulation and the tablet / diluted-oral-suspension reference; relative bioavailability F1 reduced for the 20 mg-equivalent body-weight-adjusted dose relative to the 10 mg-equivalent reference."
  reference   <- "Willmann S, Thelen K, Kubitza D, Lensing AWA, Frede M, Coboeken K, Stampfuss J, Burghaus R, Mueck W, Lippert J. Pharmacokinetics of rivaroxaban in children using physiologically based and population pharmacokinetic modelling: an EINSTEIN-Jr phase I study. Thromb J. 2018;16:32. doi:10.1186/s12959-018-0185-1"
  vignette    <- "Willmann_2018_rivaroxaban"
  units       <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; the EINSTEIN-Jr phase I single-dose study reports a single per-subject weight).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL (estimated exponent 0.323) and on central V (fixed exponent 1) with reference weight 70 kg per Willmann 2018 Table 1; Q and Vp are not weight-scaled per Willmann 2018 Results paragraph 1.",
      source_name        = "WT"
    ),
    FORM_UNDIL_SUSP = list(
      description        = "Formulation indicator at the dose record: 1 = undiluted oral suspension, 0 = tablet or diluted oral suspension (the model lumps these because Willmann 2018 found their ka statistically indistinguishable).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet or diluted oral suspension; ka 0.717 1/h, the typical-value reference).",
      notes              = "Per-dose-record covariate. Shifts ka from 0.717 1/h (reference) to 0.208 1/h (undiluted suspension); the in-vitro dissolution result reported in Willmann 2018 Discussion paragraph 5 attributes the slower undiluted-suspension ka to excipient-driven dissolution suppression at low pH.",
      source_name        = "Formulation (Tablet / Diluted suspension / Undiluted suspension)"
    ),
    DOSE_HIGH_RIV = list(
      description        = "Dose-level indicator at the dose record: 1 = 20 mg-equivalent body-weight-adjusted rivaroxaban dose, 0 = 10 mg-equivalent body-weight-adjusted dose.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (10 mg-equivalent body-weight-adjusted dose; relative bioavailability fixed to 1 by definition).",
      notes              = "Per-dose-record covariate. The 20 mg-equivalent dose has F1 = 0.648 relative to the 10 mg-equivalent reference (Willmann 2018 Table 1); the saturable-bioavailability behaviour is consistent with the adult patient popPK literature [reference 21 of Willmann 2018, Mueck 2014].",
      source_name        = "Dose level (10 mg-equivalent / 20 mg-equivalent)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59,
    n_studies      = 1,
    age_range      = "0.5-18 years (four cohorts: 0.5-<2, 2-<6, 6-<12, 12-<18 years)",
    weight_range   = "approximately 5-80 kg across the 0.5-18-year cohorts (the paper reports detailed weight ranges in the accompanying study report, reference 18)",
    sex_female_pct = NA_real_,
    disease_state  = "Children who had completed treatment for venous thromboembolism (VTE); EINSTEIN-Jr phase I single-dose pharmacokinetic study.",
    dose_range     = "Body-weight-adjusted single oral doses targeting adult exposure of either rivaroxaban 10 mg or 20 mg, delivered as tablet (children >=6 years at investigator discretion; required for children >=12 years), undiluted oral suspension, or diluted oral suspension (required for children <6 years).",
    regions        = "18 sites in 7 countries: Australia, Austria, Canada, France, Israel, Italy, United States.",
    n_observations = "206 plasma concentrations from 59 children; 7 below LLOQ (0.500 ug/L) excluded, leaving 199 for the popPK fit.",
    trial_id       = "NCT01145859 (EINSTEIN-Jr phase I).",
    notes          = "Demographic and dosing summary from Willmann 2018 Methods (Participants and study design) and Results (paragraph 1); detailed baseline demographics are reported in the accompanying clinical study publication (reference 18 of Willmann 2018, Monagle et al.). Sampling windows were age-dependent: 90 min-5 h and 20-24 h in children aged 0.5-2 years, with an additional 8-12 h sample in 2-6 year-olds; children >=6 years contributed up to four samples within 24 h post-dose."
  )

  ini({
    # Structural parameters - reference body weight 70 kg.
    # Final estimates from Willmann 2018 Table 1 ("Population estimates for the
    # first pediatric PopPK model of rivaroxaban").
    lka <- log(0.717);  label("Absorption rate constant for tablet or diluted oral suspension (1/h)")  # Willmann 2018 Table 1: ka(tablet/diluted suspension) = 0.717 1/h, RSE 21.3%
    lcl <- log(7.26);   label("Apparent clearance CL/F at 70 kg (L/h)")                                 # Willmann 2018 Table 1: CL = 7.26 L/h, RSE 9.38%
    lvc <- log(50.9);   label("Apparent central volume of distribution V/F at 70 kg (L)")              # Willmann 2018 Table 1: V = 50.9 L, RSE 12%
    lq  <- log(0.928);  label("Apparent inter-compartmental clearance Q/F (L/h)")                       # Willmann 2018 Table 1: Q = 0.928 L/h, RSE 17.5%
    lvp <- log(13.5);   label("Apparent peripheral volume of distribution Vp/F (L)")                    # Willmann 2018 Table 1: Vp = 13.5 L, RSE 51.5%

    # Formulation effect on ka: log-shift for undiluted suspension relative to
    # the tablet / diluted-suspension reference.
    e_undilsusp_ka <- log(0.208 / 0.717);  label("Log-ratio shift on ka for undiluted oral suspension (unitless)")  # Willmann 2018 Table 1: ka(undiluted suspension) = 0.208 1/h, RSE 15.4%; log-shift = log(0.208/0.717) = -1.238

    # Relative bioavailability shift for the 20 mg-equivalent dose vs the
    # 10 mg-equivalent reference (which has F1 = 1 by definition in the paper).
    lfdepot <- log(0.648);  label("Log relative bioavailability F1 for the 20 mg-equivalent dose vs the 10 mg-equivalent reference (unitless)")  # Willmann 2018 Table 1: F1 = 0.648, RSE 9.03%

    # Allometric exponents - relative to a body weight of 70 kg.
    # CL exponent is estimated; V exponent is fixed to 1 per Willmann 2018
    # Results paragraph ("the scaling exponent of V with body weight was
    # estimated to be not significantly different from 1; therefore, it was
    # fixed to 1, consistent with the allometric theory").
    allo_cl <- 0.323;         label("Allometric exponent on CL with body weight (unitless)")  # Willmann 2018 Table 1: CL exponent = 0.323, RSE 27.1%
    allo_vc <- fixed(1.0);    label("Allometric exponent on central V with body weight (unitless, fixed)")  # Willmann 2018 Table 1: V exponent fixed at 1, consistent with allometric theory

    # IIV - Willmann 2018 reports CV% (Table 1 footer "b"; exponential IIV
    # model per Methods "exponential inter-individual variability (IIV)").
    # Convert CV% to log-normal variance via omega^2 = log(1 + CV^2):
    #   ka: log(1 + 0.397^2) = log(1.15761) = 0.14640
    #   CL: log(1 + 0.262^2) = log(1.06864) = 0.06637
    etalka ~ 0.14640                                                        # Willmann 2018 Table 1: CV 39.7% on ka, RSE 63.9%
    etalcl ~ 0.06637                                                        # Willmann 2018 Table 1: CV 26.2% on CL, RSE 39.2%

    # Residual error - proportional (Willmann 2018 Results paragraph: "The
    # residual error was described by a proportional error model").
    propSd <- 0.466;  label("Proportional residual error (fraction)")  # Willmann 2018 Table 1: residual error = 46.6%, RSE 14.1%
  })

  model({
    # Reference body weight for allometric scaling (Willmann 2018 Table 1
    # description column: "Clearance for a subject with a body weight of
    # 70 kg" and "Volume of distribution ... for a subject with a body
    # weight of 70 kg").
    ref_wt <- 70

    # Absorption rate constant: tablet / diluted-suspension typical value,
    # log-shifted for the undiluted-suspension arm.
    ka <- exp(lka + etalka + e_undilsusp_ka * FORM_UNDIL_SUSP)

    # Individual PK parameters with allometric scaling on CL and central V
    # (Willmann 2018 Methods "Development of a first PopPK model in
    # children" paragraph 3: "CL and V were allometrically scaled with body
    # weight relative to a body weight of 70 kg"; "Scaling of Vp and Q with
    # body weight did not improve the fit and led to implausibly large
    # values of Vp").
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^allo_cl
    vc <- exp(lvc) * (WT / ref_wt)^allo_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # Relative bioavailability on the depot:
    #   DOSE_HIGH_RIV = 0 (10 mg-equivalent): exp(0) = 1 (reference)
    #   DOSE_HIGH_RIV = 1 (20 mg-equivalent): exp(lfdepot) = 0.648
    fdepot <- exp(lfdepot * DOSE_HIGH_RIV)

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                  k12 * central - k21 * peripheral1

    # Bioavailability on the depot.
    f(depot) <- fdepot

    # Observation: central amount is in mg (dose units), vc in L, so
    # central/vc is mg/L. The Willmann 2018 bioanalytical method reports
    # rivaroxaban plasma concentration in ug/L (Methods "Rivaroxaban plasma
    # concentrations were determined using ... calibration range from 0.500
    # ug/L (LLOQ) to 500 ug/L"); convert mg/L -> ug/L via x 1000.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}

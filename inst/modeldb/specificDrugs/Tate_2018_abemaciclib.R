Tate_2018_abemaciclib <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "linear elimination for oral abemaciclib in 224 adults with advanced",
    "cancer from the phase I first-in-human study NCT01394016 (Tate 2018).",
    "Apparent clearance and volume are fixed per patient; all exposure",
    "variability is carried by a relative bioavailability (Frel) that is",
    "both dose- and time-dependent. Frel starts at Frel,i = 1 (fixed) at",
    "the first dose and decays exponentially (half-life 14.7 days) towards",
    "a dose-dependent steady-state value Frel,ss = D50 / (D50 + Dose) with",
    "D50 = 101 mg, so repeated dosing at higher doses absorbs proportionally",
    "less drug. Serum albumin (power) and alkaline phosphatase (linear",
    "deviation) act on the global Frel. Log-normal IIV on Frel,i, Frel,ss,",
    "ka and CL/F, inter-occasion variability on CL/F, and combined",
    "additive + proportional residual error."
  )
  reference <- paste(
    "Tate SC, Sykes AK, Kulanthaivel P, Chan EM, Turner PK, Cronier DM.",
    "A Population Pharmacokinetic and Pharmacodynamic Analysis of",
    "Abemaciclib in a Phase I Clinical Trial in Cancer Patients.",
    "Clin Pharmacokinet. 2018;57(3):335-344.",
    "doi:10.1007/s40262-017-0559-8.",
    sep = " "
  )
  vignette <- "Tate_2018_abemaciclib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "abemaciclib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "abemaciclib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the global relative bioavailability,",
        "(ALB / 39)^-1.32, normalised to the 39 g/L cohort median (Tate",
        "2018 Table 1; range 24-48 g/L). Tate 2018 Sect. 2.5 treats albumin",
        "as time-varying, updated at each visit with last observation",
        "carried forward. Lower albumin gives higher Frel and exposure."
      ),
      source_name = "serum albumin"
    ),
    ALP = list(
      description = "Serum alkaline phosphatase activity",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear-deviation effect on the global relative bioavailability,",
        "(1 + 0.00197 * (ALP - 81)), centred on the 81 IU/L cohort median",
        "(Tate 2018 Table 1; range 36-1175 IU/L; IU/L and U/L are the same",
        "unit). Time-varying with last observation carried forward (Sect.",
        "2.5). Higher ALP gives higher Frel and exposure."
      ),
      source_name = "alkaline phosphatase"
    ),
    T_FIRSTDOSE = list(
      description = "Time elapsed since the first abemaciclib dose of the treatment course",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The t of Tate 2018 Eq. 5, which drives the decay of Frel from its",
        "initial value towards its steady-state value. Counted from the",
        "very first dose -- in NCT01394016 the single lead-in dose on day",
        "-3, three days before continuous dosing began -- and NOT reset at",
        "later doses. The paper does not define the origin of t in words;",
        "the day -3 origin is the one that reproduces all three end-of-",
        "cycle-1 numbers in Sect. 3.2 (Frel ~0.55 at 150 mg, ~0.50 at 200",
        "mg, a 21% rise in absorbed amount for a 33% dose increase); see",
        "the vignette. For a simulation whose first dose is at time zero,",
        "T_FIRSTDOSE equals the model time."
      ),
      source_name = "t"
    ),
    OCC = list(
      description = "Integer-valued pharmacokinetic sampling-occasion indicator for inter-occasion variability on CL/F",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Tate 2018 Table 3 reports a single IOV magnitude on CL/F (52.1%)",
        "but does not define an occasion or state their number. This file",
        "encodes four occasions, one per cycle-1 PK sampling visit listed",
        "in Sect. 2.3 (day -3, day 15, day 22, day 28), as an",
        "extraction-side construct. The index is decomposed inside model()",
        "into mutually exclusive indicators oc1..oc4 multiplexing four",
        "per-occasion CL etas; occasions 2-4 are fixed equal to occasion 1",
        "because a single shared IOV variance was estimated. Records with",
        "OCC outside 1..4 carry no IOV. For a single-occasion simulation",
        "pass OCC = 1."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5, ESM Fig. 3) but not retained; no estimate reported."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Tate 2018 Sect. 2.5) but not retained; no estimate reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 224L,
    n_studies = 1L,
    n_observations = "4012 abemaciclib plasma concentrations (Tate 2018 Table 2)",
    age_range = "24-85 years (median 61)",
    weight_range = "43.6-175 kg (median 70.4)",
    sex_female_pct = 67,
    race_ethnicity = c(White = 94, Black = 2, Asian = 4),
    disease_state = paste(
      "Adults with advanced cancer: mixed tumours (dose escalation), NSCLC,",
      "glioblastoma, breast cancer (including HR+ with fulvestrant),",
      "melanoma and colorectal cancer"
    ),
    dose_range = "Oral abemaciclib capsules 50-225 mg q24h or 75-275 mg q12h",
    regions = "Not reported (multicentre)",
    notes = paste(
      "Phase I NCT01394016 (Tate 2018 Sect. 2.1, Tables 1-2). PK model built",
      "on cycle-1 data: single dose on day -3, continuous dosing days 1-28.",
      "54 dose reductions in 48 patients during cycle 1, mostly 200 to 150",
      "mg q12h."
    )
  )

  ini({
    # ----- Disposition and absorption (Tate 2018 Table 3, final model) -----
    lka <- log(0.197); label("First-order absorption rate constant ka (1/h)")  # Table 3 final model: ka = 0.197 /h (6.6 %SEE)
    lcl <- log(35.9); label("Apparent clearance CL/F (L/h)")  # Table 3 final model: CL/F = 35.9 L/h (5.8 %SEE)
    lvc <- log(1050); label("Apparent volume of distribution Vd/F (L)")  # Table 3 final model: Vd/F = 1050 L (5.6 %SEE)

    # ----- Relative bioavailability (Tate 2018 Eqs. 4-5, Table 3) -----
    # Frel(t) = Frel,ss - (Frel,ss - Frel,i) * exp(-lambda * t), with
    # Frel,ss = 1 - Dose / (D50 + Dose). Frel,ss has no typical value of its
    # own; it is fully determined by D50 and the dose.
    lfdepot <- fixed(log(1)); label("Initial relative bioavailability Frel,i (unitless)")  # Table 3 final model: Frel,i = 1 (fixed)
    led50 <- log(101); label("Dose at which Frel,ss equals 0.5, D50 (mg)")  # Table 3 final model: D50 = 101 mg (20 %SEE)
    lfdepot_exp_kdes <- log(0.00197); label("Rate constant lambda of the transition from Frel,i to Frel,ss (1/h)")  # Table 3 final model: lambda = 0.00197 /h (9.5 %SEE)

    # ----- Covariate effects on the global Frel (Table 3 footnote b) -----
    # Frel = TV_Frel * (ALB / 39)^theta1 * (1 + theta2 * (ALP - 81))
    e_alb_fdepot <- -1.32; label("Power exponent of serum albumin on Frel (unitless)")  # Table 3 final model: theta1 serum albumin on Frel = -1.32 (20 %SEE)
    e_alp_fdepot <- 0.00197; label("Linear coefficient of alkaline phosphatase on Frel (per U/L)")  # Table 3 final model: theta2 alkaline phosphatase on Frel = 0.00197 (29 %SEE)

    # ----- Inter-individual variability (Table 3 final model, IIV column) -----
    # The tabulated percentages are read as omega x 100, so the variance is
    # (P/100)^2. The paper does not state the convention; the vignette shows
    # that this reading reproduces the width of the Fig. 1 VPC band and the
    # log(1 + CV^2) reading does not.
    etalfdepot ~ 0.538756  # Table 3 final model: IIV Frel,i 73.4 -> 0.734^2
    etalfdepot_exp_inf ~ 1.5129  # Table 3 final model: IIV Frel,ss 123 -> 1.23^2
    etalka ~ 0.602176  # Table 3 final model: IIV ka 77.6 -> 0.776^2
    etalcl ~ 0.087616  # Table 3 final model: IIV CL/F 29.6 -> 0.296^2

    # ----- Inter-occasion variability on CL/F (Table 3 final model, IOV column) -----
    # One shared variance (52.1 -> 0.521^2); occasions 2-4 repeat it as fixed.
    etaiov_cl_1 ~ 0.271441  # Table 3 final model: IOV CL/F 52.1 -> 0.521^2
    etaiov_cl_2 ~ fixed(0.271441)
    etaiov_cl_3 ~ fixed(0.271441)
    etaiov_cl_4 ~ fixed(0.271441)

    # ----- Residual error (Table 3 final model) -----
    addSd <- 14.8; label("Additive residual error (ng/mL)")  # Table 3 final model: additive 14.8 ng/mL (27 %SEE)
    propSd <- 0.171; label("Proportional residual error (fraction)")  # Table 3 final model: proportional 17.1% (8.9 %SEE)
  })

  model({
    # ----- 1. Occasion indicators for IOV on CL/F -----
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # ----- 2. Individual disposition and absorption -----
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc)
    kel <- cl / vc

    # ----- 3. Dose- and time-dependent relative bioavailability -----
    # Eq. 4: Frel,ss = 1 - Dose / (D50 + Dose) = D50 / (D50 + Dose), where
    # Dose is the amount of the dose being absorbed, read off the event
    # record with podo(depot) so a dose reduction changes Frel,ss at once.
    ed50 <- exp(led50)
    fdepot_i <- exp(lfdepot + etalfdepot)
    fdepot_exp_inf <- ed50 / (ed50 + podo(depot)) * exp(etalfdepot_exp_inf)
    # Eq. 5 rewritten in the exp-decay form P(t) = P_inf + P_component * exp(-kdes * t):
    # Frel = Frel,ss - (Frel,ss - Frel,i) * exp(-lambda * t).
    fdepot_exp_component <- fdepot_i - fdepot_exp_inf
    fdepot_exp_kdes <- exp(lfdepot_exp_kdes)
    frel_typ <- fdepot_exp_inf + fdepot_exp_component * exp(-fdepot_exp_kdes * T_FIRSTDOSE)
    # Table 3 footnote b: covariates act on the global Frel.
    fdepot <- frel_typ * (ALB / 39)^e_alb_fdepot * (1 + e_alp_fdepot * (ALP - 81))

    # ----- 4. ODE system -----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- fdepot

    # ----- 5. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

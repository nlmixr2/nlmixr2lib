Brill_2015_midazolam <- function() {
  description <- paste(
    "Three-compartment population PK model for midazolam with a five-transit-",
    "compartment oral absorption chain (Ka = Ktr), supporting oral and",
    "intravenous dosing, in 20 morbidly obese patients studied before and 18",
    "of them again one year after bariatric surgery (Roux-en-Y gastric bypass",
    "or sleeve gastrectomy). Post-surgery status multiplies clearance by 1.68",
    "and inter-compartmental clearance Q by 3.22 and has its own transit rate",
    "constant; before surgery, total body weight scales the central volume",
    "linearly and both peripheral volumes by a power function (reference",
    "127 kg); the second peripheral volume is a fixed fraction of the first."
  )
  reference <- paste(
    "Brill MJ, van Rongen A, van Dongen EP, van Ramshorst B, Hazebroek EJ,",
    "Darwich AS, Rostami-Hodjegan A, Knibbe CA (2015). The Pharmacokinetics",
    "of the CYP3A Substrate Midazolam in Morbidly Obese Patients Before and One",
    "Year After Bariatric Surgery. Pharm Res 32(12):3927.",
    "doi:10.1007/s11095-015-1752-9."
  )
  vignette <- "Brill_2015_midazolam"
  dosing <- c("depot", "central")
  units <- list(time = "min", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Brill 2015 Table II: 'BW body weight (median = 127 kg for all data)'.",
        "Enters only while PRIOR_BARIATRIC = 0 (before surgery): linearly on",
        "the central volume, Vc = 37.3 * (1 + 0.0435 * (WT - 127)), and as a",
        "power on the first (and therefore the second) peripheral volume,",
        "Vp = 106 * (WT / 127)^3.93. After surgery the volumes are the",
        "WT-independent typical values (Results: 'without a significant",
        "influence of (lean) bodyweight within this group'). The linear Vc",
        "term turns negative below WT = 104 kg, so the pre-surgery branch is",
        "only valid within the observed pre-surgery range of 112-186 kg.",
        "Record the weight measured at each occasion."
      ),
      source_name = "BW"
    ),
    PRIOR_BARIATRIC = list(
      description = "Bariatric surgery performed before this record (1 = post-surgery, 0 = before surgery)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (morbidly obese patient before bariatric surgery)",
      notes = paste(
        "Brill 2015 covariate 'bariatric surgery' (Methods Eq. 3, Z^COVARIATE).",
        "0 = occasion 1 (day of surgery, sampled before and during the",
        "operation); 1 = occasion 2, 52 +/- 3 weeks after laparoscopic",
        "Roux-en-Y gastric bypass (16 patients) or sleeve gastrectomy (2).",
        "Time-varying within subject: the same patient carries 0 at the first",
        "occasion and 1 at the second, with one set of random effects shared",
        "across both (no inter-occasion variability was reported)."
      ),
      source_name = "bariatric surgery"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "midazolam", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "Before surgery 26-57 years (mean 43.6, SD 7.6); after surgery 27-58 years (mean 45.5, SD 7.4)",
    weight_range = "Before surgery 112-186 kg (mean 144.4, SD 21.7); after surgery 62-138 kg (mean 98.3, SD 18.0)",
    weight_median = "127 kg (pooled over both occasions; covariate reference)",
    sex_female_pct = 60,
    disease_state = paste(
      "Morbid obesity (BMI before surgery 40-68 kg/m^2, mean 47.1) in adults",
      "undergoing laparoscopic Roux-en-Y gastric bypass or sleeve gastrectomy;",
      "restudied one year after surgery (BMI 24-50 kg/m^2, mean 31.9; mean",
      "weight loss 44.5 kg). No CYP3A inducers or inhibitors, no renal",
      "insufficiency (MDRD4 eGFR >= 60 mL/min)."
    ),
    dose_range = paste(
      "At each occasion 7.5 mg oral midazolam tablet followed by 5 mg IV",
      "160 +/- 48 min later (at induction of anaesthesia on occasion 1, at",
      "150 min on occasion 2)."
    ),
    regions = "Netherlands",
    notes = paste(
      "Brill 2015 Table I: 20 patients (12 F / 8 M) before surgery, 18 of them",
      "(11 F / 7 M) at 51.8 +/- 2.5 weeks after surgery; 22 +/- 3 and 21 +/- 1",
      "samples per patient per occasion. NTC01519726 / EudraCT",
      "2011-003293-93. Occasion 1 is the same cohort published as Brill 2014",
      "(doi:10.1007/s40262-014-0166-x), refitted here jointly with occasion 2."
    )
  )

  ini({
    # Final model of the simultaneous analysis, Brill 2015 Table II
    # ('Final model of simultaneous analysis', Value (RSE)).

    # Transit / absorption rate (Ka = Ktr), estimated separately per occasion.
    lka_presurg <- log(0.117); label("Ka = Ktr before bariatric surgery (1/min)") # Table II 'Ka = Ktr Morbidly obese' 0.117 (15%)
    lka_postsurg <- log(0.267); label("Ka = Ktr after bariatric surgery (1/min)") # Table II 'Ka = Ktr Bariatric patients' 0.267 (19%)

    lfdepot <- log(0.537); label("Oral bioavailability (fraction)") # Table II F 0.537 (9%)

    lcl <- log(0.385); label("Clearance before bariatric surgery (L/min)") # Table II 'CL Morbidly obese' 0.385 (11%)
    e_prior_bariatric_cl <- 1.68; label("Multiplicative factor on CL after bariatric surgery (unitless)") # Table II 'fCL Bariatric patients' 1.68 (7%)

    lq <- log(0.669); label("Inter-compartmental clearance to peripheral1 before surgery (L/min)") # Table II 'Q Morbidly obese' 0.669 (24%)
    e_prior_bariatric_q <- 3.22; label("Multiplicative factor on Q after bariatric surgery (unitless)") # Table II 'fQ Bariatric patients' 3.22 (32%)

    lq2 <- log(0.551); label("Inter-compartmental clearance to peripheral2 (L/min)") # Table II Q2 0.551 (23%)

    lvc <- log(37.3); label("Central volume at WT = 127 kg before surgery, and after surgery (L)") # Table II 'Vmedian BW' 37.3 (18%) = 'Vcentral Bariatric patients' 37.3
    e_wt_vc <- 0.0435; label("Linear WT slope on Vc before surgery (1/kg)") # Table II X 0.0435 (92%)

    lvp <- log(106); label("Peripheral1 volume at WT = 127 kg before surgery, and after surgery (L)") # Table II 'Vmedian BW' 106 (17%) = 'V1st peripheral Bariatric patients' 106
    e_wt_vp <- 3.93; label("Power exponent of WT/127 on peripheral volumes before surgery (unitless)") # Table II Y 3.93 (20%)

    # Table II estimates the second peripheral volume as a fraction fV of the
    # first (fV = 0.359, '*V 1st peripheral = 38 L' i.e. 0.359 * 106 L).
    # With no IIV on the peripheral volumes this equals Vp2 = 38.05 L scaled
    # by the same WT power as Vp, which is how model() applies it.
    lvp2 <- log(0.359 * 106); label("Peripheral2 volume at WT = 127 kg before surgery, and after surgery (L)") # Table II 'fV2nd peripheral' 0.359 (27%) x Vmedian BW 106 L

    # IIV, Table II CV%; omega^2 = log(1 + CV^2). No IIV on Q, Q2 or the
    # peripheral volumes ('V1st peripheral 0 FIX').
    etalka ~ log(1 + 0.424^2) # Table II IIV 'Ktr = Ka' 42.4% (15%), shared by both occasions
    etalcl ~ log(1 + 0.197^2) # Table II IIV CL 19.7% (38%)
    etalfdepot ~ log(1 + 0.334^2) # Table II IIV F 33.4% (18%)
    etalvc ~ log(1 + 0.537^2) # Table II IIV Vcentral 53.7% (39%)

    propSd <- 0.421; label("Proportional residual error (fraction)") # Table II proportional residual error 42.1% (5%)
  })
  model({
    # PRIOR_BARIATRIC = 1 after surgery, 0 before (Methods Eq. 3-8).
    pre <- 1 - PRIOR_BARIATRIC

    ka <- exp(lka_presurg * pre + lka_postsurg * PRIOR_BARIATRIC + etalka)
    ktr <- ka # Methods: transit rates equalized to Ka; MTT = (N + 1) / Ktr with N = 5
    fdepot <- exp(lfdepot + etalfdepot)
    cl <- exp(lcl + etalcl) * e_prior_bariatric_cl^PRIOR_BARIATRIC
    q <- exp(lq) * e_prior_bariatric_q^PRIOR_BARIATRIC
    q2 <- exp(lq2)

    # Body weight acts only before surgery (Table II equations with
    # median BW = 127 kg); after surgery the volumes are the reference values.
    vc <- exp(lvc + etalvc) * (PRIOR_BARIATRIC + pre * (1 + e_wt_vc * (WT - 127)))
    wt_vp <- PRIOR_BARIATRIC + pre * (WT / 127)^e_wt_vp
    vp <- exp(lvp) * wt_vp
    vp2 <- exp(lvp2) * wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(central) <- ktr * transit5 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot) <- fdepot

    # ug / L = ng/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

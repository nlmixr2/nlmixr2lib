vanRongen_2017_midazolam <- function() {
  description <- paste(
    "Two-compartment population PK model for midazolam in 19 obese",
    "adolescents (12-18.9 years, total body weight 62-149.8 kg, BMI",
    "24.8-55 kg/m^2) and 20 morbidly obese adults (26-57 years, total",
    "body weight 112.3-186.3 kg, BMI 39.9-67.6 kg/m^2), with a",
    "five-transit-compartment first-order oral absorption chain (Ka = Ktr)",
    "supporting oral and intravenous dosing (van Rongen 2017 Final model,",
    "Table 2). Study population (adolescent vs morbidly obese adult)",
    "separates clearance into two cohort-specific values (CL_104.7 kg in",
    "adolescents with an estimated TBW power on top, CL fixed-across-WT in",
    "morbidly obese adults). The same V_141.8 kg central value of the",
    "peripheral compartment is shared between cohorts but TBW power",
    "scaling applies only to morbidly obese adults. Central volume,",
    "inter-compartmental clearance Q, transit-absorption rate Ka = Ktr,",
    "and oral bioavailability F are shared across both cohorts. Oral",
    "data were collected only in morbidly obese adults; adolescents",
    "received only IV bolus doses."
  )
  reference <- paste(
    "van Rongen A, Brill MJE, Vaughns JD, Valitalo PAJ, van Dongen EPA,",
    "van Ramshorst B, Barrett JS, van den Anker JN, Knibbe CAJ (2018).",
    "Higher Midazolam Clearance in Obese Adolescents Compared with",
    "Morbidly Obese Adults.",
    "Clin Pharmacokinet 57(5):601-611.",
    "doi:10.1007/s40262-017-0579-4.",
    "PK structure builds on the morbidly-obese-adult midazolam model of",
    "Brill et al. (2014); see modellib('Brill_2014_midazolam').",
    "A sibling extraction of the same paper sourced from the DDMORE",
    "Foundation Repository bundle (DDMODEL00000250) and re-fit on the",
    "bundle's simulated 9-subject dataset is available as",
    "modellib('vanRongen_2018_midazolam'); this entry uses the published",
    "Table 2 Final-model point estimates from the full 39-subject pooled",
    "dataset instead.",
    sep = " "
  )
  vignette <- "vanRongen_2017_midazolam"
  units    <- list(time = "minute", dosing = "microgram", concentration = "microgram/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'TBW' (total body weight, kg) maps to the canonical WT.",
        "Time-fixed at baseline. Enters the Final model in two cohort-specific",
        "ways (van Rongen 2017 Table 2 footnotes): power-law on CL in obese",
        "adolescents with reference 104.7 kg (the adolescent-cohort median TBW)",
        "and power-law on peripheral volume in morbidly obese adults with",
        "reference 141.8 kg (the adult-cohort median TBW). Adolescent peripheral",
        "volume and adult clearance are flat across WT in the Final model. The",
        "pooled-cohort observed TBW range is 62-186.3 kg."
      ),
      source_name        = "TBW"
    ),
    ADOLESCENT = list(
      description        = "Obese-adolescent cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = morbidly obese adult (BMI > 40 kg/m^2; 26-57 years; 112.3-186.3 kg)",
      notes              = paste(
        "1 = obese adolescent (overweight BMI-for-age 85th-95th percentile or",
        "obese BMI-for-age >= 95th percentile; 12-18.9 years; 62-149.8 kg).",
        "0 = morbidly obese adult (BMI > 40 kg/m^2; 26-57 years; 112.3-186.3 kg).",
        "Per van Rongen 2017 Methods 2.5, the variable 'study population'",
        "(adolescents [27] vs adults [14]) was retained as a binary covariate",
        "on clearance after stepwise forward inclusion and backward exclusion.",
        "Source paper uses no specific column name; the indicator encodes the",
        "study-cohort partition described in Tables 1 and 2. This canonical",
        "encoding pairs the indicator with WT to recover the cohort-specific",
        "Final-model parametrisation. The model is intended for subjects in",
        "one of these two cohorts; using it for non-obese subjects or",
        "non-adolescent / non-adult age groups is outside the calibration",
        "envelope."
      ),
      source_name        = NA_character_
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 2L,
    age_range      = "Obese adolescents 12.5-18.9 years (mean 15.9, SD 1.6); morbidly obese adults 26-57 years (mean 43.6, SD 7.6).",
    weight_range   = "Obese adolescents 62-149.8 kg (mean 102.7, SD 24.9); morbidly obese adults 112.3-186.3 kg (mean 144.4, SD 21.7).",
    weight_median  = "104.7 kg (adolescent reference for CL power), 141.8 kg (adult reference for Vp power).",
    sex_female_pct = 64,
    race_ethnicity = paste(
      "Obese adolescents (n = 19): 5 Caucasian, 9 African American, 5 Hispanic",
      "(per Discussion). Morbidly obese adults (n = 20): 19 Caucasian,",
      "1 African American."
    ),
    disease_state  = paste(
      "Obese adolescents (n = 19): three overweight (BMI-for-age 85th-95th",
      "percentile) and 16 obese (BMI-for-age >= 95th percentile) adolescents",
      "scheduled for general surgery (orthopaedic surgery, tonsillectomy,",
      "bariatric surgery). Morbidly obese adults (n = 20): BMI > 40 kg/m^2",
      "scheduled for bariatric surgery."
    ),
    dose_range     = paste(
      "Obese adolescents: single IV bolus 2 or 3 mg midazolam given a few",
      "minutes before transfer to the operating room. Morbidly obese adults:",
      "7.5 mg oral midazolam followed by a 5 mg IV bolus at induction of",
      "anaesthesia (159 +/- 67 min after the oral dose)."
    ),
    regions        = paste(
      "Obese adolescents: Children's National Health System, Washington DC,",
      "USA (IRB Protocol No. 4718). Morbidly obese adults: St. Antonius",
      "Hospital, Nieuwegein, The Netherlands (VCMO NL35861.100.11, EudraCT",
      "2011-003293-93)."
    ),
    notes          = paste(
      "Pooled analysis of 530 midazolam plasma concentration observations",
      "from 19 obese adolescents (129 samples, all above the 0.5 ng/mL LLOQ)",
      "and 20 morbidly obese adults (401 samples retained; 33 of 434 were",
      "below 0.3 ng/mL LOD and excluded, while 9 samples between 0.3-0.8",
      "ng/mL were kept). Demographics per Table 1. Sex split: 13 F / 6 M",
      "adolescents, 12 F / 8 M morbidly obese adults (25/39 = 64% female",
      "overall). NONMEM 7.2 first-order conditional estimation with",
      "interaction; internal validation by 1000-replicate stratified",
      "bootstrap (93% successful) and 1000-simulation normalised prediction",
      "distribution error analysis (Methods 2.4)."
    )
  )

  ini({
    # Final model parameter values from van Rongen 2017 Table 2 ("Final model"
    # column with RSE%). The Final model retains study population as a binary
    # covariate on CL (delta-OFV -8.0, p < 0.01), TBW as a power covariate on
    # CL within the adolescent cohort (delta-OFV -10.6, p < 0.01), and TBW as
    # a power covariate on V_peripheral within the morbidly obese adult cohort
    # (delta-OFV -10.9, p < 0.001) per Results 3.2. The covariate-screening
    # round eliminated age, BMI, lean body weight (Janmahasatian formula),
    # race, and sex (all p > 0.05). Bioavailability is only meaningful for the
    # oral arm of the morbidly obese adult cohort; adolescents received only
    # IV doses.

    # Clearance for obese-adolescent reference subject (TBW = 104.7 kg).
    lcl_adol     <- log(0.71)
    label("Clearance in obese adolescents at TBW = 104.7 kg (L/min)")        # Table 2 final CL_104.7 kg = 0.71 L/min (RSE 7%)

    # TBW power exponent on CL in obese adolescents.
    e_wt_cl_adol <- 1.2
    label("Power exponent of (WT / 104.7) on CL in obese adolescents")       # Table 2 final W = 1.2 (RSE 31%)

    # Clearance for morbidly obese adult cohort (no TBW dependency in Final
    # model; "No significant trend was found for TBW and clearance in the
    # morbidly obese adults" per Results 3.2).
    lcl_adult    <- log(0.44)
    label("Clearance in morbidly obese adults (L/min)")                      # Table 2 final CL_morbidly obese adults = 0.44 L/min (RSE 11%)

    # Oral bioavailability (shared across cohorts; only relevant for adult
    # oral doses since adolescents received no oral midazolam).
    lfdepot      <- log(0.562)
    label("Oral bioavailability F (unitless)")                                # Table 2 final F = 0.562 (RSE 12%)

    # Absorption rate constant Ka = Ktr (transit-compartment rate constant
    # equalised to absorption rate constant per Methods 2.4).
    lka          <- log(0.115)
    label("First-order absorption / transit rate constant Ka = Ktr (1/min)") # Table 2 final K_a = K_tr = 0.115 1/min (RSE 11%)

    # Central volume of distribution (shared across cohorts; no WT covariate
    # retained in the Final model).
    lvc          <- log(55.2)
    label("Central volume of distribution Vc (L)")                            # Table 2 final V_central = 55.2 L (RSE 11%)

    # Peripheral volume reference. The Final model shares the same
    # V_141.8 kg = 172 L between the two cohorts (Table 2 reports the same
    # 172 (13%) entry for both rows). TBW power scaling applies only to the
    # morbidly obese adult cohort.
    lvp          <- log(172)
    label("Peripheral volume of distribution Vp at TBW = 141.8 kg (L)")       # Table 2 final V_peripheral obese adolescents / V_141.8 kg = 172 L (RSE 13%)

    # TBW power exponent on Vp in morbidly obese adults only (not applied
    # to adolescents; "Nor for TBW and peripheral volume of distribution of
    # midazolam in the obese adolescents (p > 0.05)" per Results 3.2).
    e_wt_vp_adult <- 3.3
    label("Power exponent of (WT / 141.8) on Vp in morbidly obese adults")   # Table 2 final X = 3.3 (RSE 33%)

    # Inter-compartmental clearance (shared across cohorts).
    lq           <- log(1.14)
    label("Inter-compartmental clearance Q (L/min)")                          # Table 2 final Q = 1.14 L/min (RSE 12%)

    # Inter-individual variability. Van Rongen 2017 Table 2 reports IIVs as
    # CV%; internal variance on the log-eta scale is omega^2 = log(CV^2 + 1).
    #   CL   :  21.0% -> log(0.210^2 + 1) = 0.04316
    #   F    :  39.2% -> log(0.392^2 + 1) = 0.14296
    #   Ka   :  49.5% -> log(0.495^2 + 1) = 0.21915
    #   Vc   :  58.5% -> log(0.585^2 + 1) = 0.29440
    #   Vp   :  42.2% -> log(0.422^2 + 1) = 0.16395
    #   Q    :  42.4% -> log(0.424^2 + 1) = 0.16539
    # No off-diagonal correlations are reported in Table 2; etas are
    # treated as independent log-normal variates as in the source NONMEM
    # implementation (the paper's text references only a diagonal OMEGA).
    etalcl_adol     ~ 0.04316    # Table 2 final IIV CL 21% (RSE 26%); applied to both adolescent and adult typical CL via a shared eta
    etalfdepot      ~ 0.14296    # Table 2 final IIV F 39.2% (RSE 21%)
    etalka          ~ 0.21915    # Table 2 final IIV K_a = K_tr 49.5% (RSE 17%)
    etalvc          ~ 0.29440    # Table 2 final IIV V_central 58.5% (RSE 11%)
    etalvp          ~ 0.16395    # Table 2 final IIV V_peripheral 42.2% (RSE 22%)
    etalq           ~ 0.16539    # Table 2 final IIV Q 42.4% (RSE 25%)

    # Proportional residual error (Final model; van Rongen 2017 Methods 2.4
    # "Residual variability was tested using proportional, additive, or
    # combined proportional and additive error models" with the Final model
    # retaining proportional error only per Table 2).
    propSd <- 0.297
    label("Proportional residual error (fraction)")                           # Table 2 final proportional error 29.7% (RSE 9%)
  })

  model({
    # Cohort-specific typical CL. ADOLESCENT == 1 selects the obese-
    # adolescent parametrisation with reference 104.7 kg and estimated TBW
    # power exponent W = e_wt_cl_adol; ADOLESCENT == 0 selects the
    # morbidly-obese-adult constant CL with no WT dependence (van Rongen
    # 2017 Table 2 Final model). The same eta is applied to both branches
    # to retain a single IIV omega for CL across the pooled cohort
    # (matching the paper's single CL row in the IIV section of Table 2).
    cl_typ_adol  <- exp(lcl_adol)  * (WT / 104.7)^e_wt_cl_adol
    cl_typ_adult <- exp(lcl_adult)
    cl_typ       <- ADOLESCENT * cl_typ_adol + (1 - ADOLESCENT) * cl_typ_adult
    cl           <- cl_typ * exp(etalcl_adol)

    # Cohort-specific typical Vp. Adolescents have no TBW dependence; adults
    # carry the (WT / 141.8)^X = (WT / 141.8)^e_wt_vp_adult power. Both
    # branches share the same reference value V_141.8 kg = 172 L (Table 2).
    vp_typ <- exp(lvp) * (WT / 141.8)^(e_wt_vp_adult * (1 - ADOLESCENT))
    vp     <- vp_typ * exp(etalvp)

    # Shared parameters across cohorts.
    fdepot <- exp(lfdepot + etalfdepot)
    ka     <- exp(lka     + etalka)
    vc     <- exp(lvc     + etalvc)
    q      <- exp(lq      + etalq)
    ktr    <- ka                                                              # van Rongen 2017 Methods 2.4: "K_tr was equalised to K_a"

    # Micro-constants for the two-compartment central-disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Compartment ordering: depot (oral input) -> 5 transit compartments
    # (van Rongen 2017 Methods 2.4: "a five transit absorption compartment
    # model") -> central. IV doses (bolus) load central directly via
    # cmt = central in the event table; oral doses load depot.
    d/dt(depot)       <- -ka  * depot
    d/dt(transit1)    <-  ka  * depot     - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2  - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3  - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4  - ktr * transit5
    d/dt(central)     <-  ktr * transit5  - kel * central -
                          k12 * central   + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central   - k21 * peripheral1

    # Oral bioavailability on the depot (IV doses bypass).
    f(depot) <- fdepot

    # Observed plasma midazolam concentration. Dose in microgram, volume in
    # L gives microgram/L which is numerically equal to ng/mL (the units
    # used in van Rongen 2017 for the LLOQ and concentration reporting).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

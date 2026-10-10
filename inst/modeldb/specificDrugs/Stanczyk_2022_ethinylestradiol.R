Stanczyk_2022_ethinylestradiol <- function() {
  description <- "One-compartment population PK model for ethinyl estradiol (EE) delivered by the once-weekly levonorgestrel/ethinyl estradiol transdermal delivery system (TDS) in healthy women (Stanczyk 2022). Each 7-day patch splits its bioavailable dose between a zero-order input directly into the central compartment over the 168-h wear period (fraction Frel) and a bolus into a first-order absorption depot (fraction 1 - Frel). Overall bioavailability F = 30% is fixed. Dose each patch as TWO records of the patch EE content (2300 ug), one to `central` with `rate = -2` and one to `depot`; f() splits the content between the two routes."
  reference <- "Stanczyk FZ, Archer DF, Lohmer LRL, Pirone J, Previtera M, Korner P. Extended regimen of a levonorgestrel/ethinyl estradiol transdermal delivery system: Predicted serum hormone levels using a population pharmacokinetic model. PLoS ONE. 2022;17(12):e0279640. doi:10.1371/journal.pone.0279640"
  vignette <- "Stanczyk_2022_levonorgestrel_ethinylestradiol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "ethinyl estradiol", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethinyl estradiol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "18-45 years (enrolment criterion)",
    sex_female_pct = 100,
    bmi_range = "18-32 kg/m^2 (enrolment criterion)",
    disease_state = "Healthy, normotensive, nonsmoking women with 24- to 35-day menstrual intervals",
    dose_range = "LNG/EE TDS (2.6 mg levonorgestrel, 2.3 mg ethinyl estradiol per patch; nominal delivery 120 ug/day LNG and 30 ug/day EE), one patch per week for 3 weeks followed by a 1-week patch-free period",
    regions = "United States (single centre)",
    notes = "Phase 1 open-label randomized study ATI-CL14 (NCT01243580), conducted 2009. The popPK model was fitted only to cycle-2 data from the 18 women who wore the TDS in both cycles 1 and 2 (weeks 1 and 3 of cycle 2 sampled; 23 samples per subject). No covariate analysis was performed. Stanczyk 2022 Methods 'Study design' and 'PopPK models'."
  )

  ini({
    lka <- log(0.0204); label("First-order absorption rate constant from the depot Ka (1/h)") # Stanczyk 2022 S1 Table: Ka = 0.0204 1/h (CV 14.1%)
    lcl <- log(107); label("Apparent clearance CL (L/h)") # Stanczyk 2022 S1 Table: Clearance = 107 L/h (CV 7.37%)
    lvc <- log(2490); label("Apparent central volume of distribution V (L)") # Stanczyk 2022 S1 Table: Volume = 2490 L (CV 8.16%)
    logitfrel <- logit(0.732); label("Logit of the fraction of the bioavailable dose delivered by the zero-order input into central (Frel)") # Stanczyk 2022 S1 Table: relative bioavailability for zero-order absorption Frel = 73.2 (CV 8.24%); no IIV
    lfdepot <- fixed(log(0.30)); label("Overall bioavailability F of the patch EE content (fraction)") # Stanczyk 2022 S1 Table: Bioavailability F = 30.0, fixed; Methods 'EE popPK models'
    ld1 <- fixed(log(168)); label("Duration of the zero-order input into central D1 (h)") # Stanczyk 2022 Methods 'EE popPK models': duration of the zero-order absorption fixed to 168 h

    # S1 Table reports IIV as CV%; omega^2 = log(1 + CV^2) for Phoenix
    # exponential (log-normal) random effects. No correlations reported.
    etalvc ~ log(1 + 0.434^2) # Stanczyk 2022 S1 Table: IIV Volume 43.4% CV -> log(1 + 0.434^2)
    etalcl ~ log(1 + 0.346^2) # Stanczyk 2022 S1 Table: IIV Clearance 34.6% CV -> log(1 + 0.346^2)
    etalka ~ log(1 + 0.531^2) # Stanczyk 2022 S1 Table: IIV Ka 53.1% CV -> log(1 + 0.531^2)

    propSd <- 0.269; label("Proportional residual error (fraction)") # Stanczyk 2022 S1 Table: multiplicative residual error 26.9 (CV 9.17%); Cobs = C*(1 + CEps)
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    frel <- expit(logitfrel)
    fbio <- exp(lfdepot)
    d1 <- exp(ld1)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # The patch content is dosed to both states: Frel * F enters central as
    # a zero-order input over the wear period, (1 - Frel) * F enters the
    # first-order depot as a bolus.
    f(central) <- fbio * frel
    dur(central) <- d1
    f(depot) <- fbio * (1 - frel)

    # ug / L = ng/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

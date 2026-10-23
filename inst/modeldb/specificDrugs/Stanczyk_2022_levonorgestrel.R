Stanczyk_2022_levonorgestrel <- function() {
  description <- "One-compartment population PK model for levonorgestrel (LNG) delivered by the once-weekly levonorgestrel/ethinyl estradiol transdermal delivery system (TDS) in healthy women (Stanczyk 2022). Each 7-day patch splits its bioavailable dose between a short zero-order infusion into a first-order absorption depot (fraction Frel) and a zero-order input directly into the central compartment over the 168-h wear period (fraction 1 - Frel). Overall bioavailability F = 27% is fixed. Clearance and volume are lower from the third consecutive weekly patch onwards (multipliers 0.87 and 0.85). Dose each patch as TWO records of the patch LNG content (2600 ug): one to `depot` with `rate = -2` (modelled duration D1) and one to `central` with an explicit `dur = 168`, which equals the model's fixed D2 and avoids an rxode2 5.1.8 mis-solve of back-to-back modelled-duration infusions when CYCLE changes on the dose row (see vignette); f() splits the content between the two routes."
  reference <- "Stanczyk FZ, Archer DF, Lohmer LRL, Pirone J, Previtera M, Korner P. Extended regimen of a levonorgestrel/ethinyl estradiol transdermal delivery system: Predicted serum hormone levels using a population pharmacokinetic model. PLoS ONE. 2022;17(12):e0279640. doi:10.1371/journal.pone.0279640"
  vignette <- "Stanczyk_2022_levonorgestrel_ethinylestradiol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    CYCLE = list(
      description = "Consecutive weekly patch number since the last patch-free week (1 = first patch, 2 = second, 3 = third, ...)",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = "Time-varying; set it to the number of the patch currently being worn and supply it on every row. The paper's 'week 3 effect' was estimated from week-1 and week-3 data of a 3-patch cycle; it is applied here to CYCLE >= 3, which reproduces the paper's own week-3 and week-12 simulations of the 12-week extended regimen (S4 Table, Table 3), including the week-3 pre-dose concentration that requires weeks 1 and 2 to use the reference parameters.",
      source_name = "week (Week 3 effect)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "levonorgestrel", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levonorgestrel", units = "ug", specimen = "plasma", verified = TRUE)
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
    lka <- log(0.0281); label("First-order absorption rate constant from the depot Ka (1/h)") # Stanczyk 2022 S2 Table: Ka = 0.0281 1/h (CV 31.2%)
    lcl <- log(2.46); label("Apparent clearance CL for patches 1 and 2 (L/h)") # Stanczyk 2022 S2 Table: Clearance = 2.46 L/h (CV 10.9%)
    lvc <- log(163); label("Apparent central volume of distribution V for patches 1 and 2 (L)") # Stanczyk 2022 S2 Table: Volume = 163 L (CV 15.0%)
    e_cycle3_cl <- 0.87; label("Multiplier on CL from the third consecutive patch onwards (unitless)") # Stanczyk 2022 S2 Table: Week 3 effect on clearance = 87.0 (CV 46.7%); read as CL_week3 = 0.870 * CL, see vignette
    e_cycle3_vc <- 0.85; label("Multiplier on V from the third consecutive patch onwards (unitless)") # Stanczyk 2022 S2 Table: Week 3 effect on volume = 85.0 (CV 57.0%); read as V_week3 = 0.850 * V, see vignette
    lfrel <- log(0.378); label("Fraction of the bioavailable dose routed to the first-order depot (Frel)") # Stanczyk 2022 S2 Table: relative bioavailability for first-order absorption Frel = 37.8 (CV 23.2%)
    lfdepot <- fixed(log(0.27)); label("Overall bioavailability F of the patch LNG content (fraction)") # Stanczyk 2022 S2 Table: Overall bioavailability F = 27.0, fixed; Methods 'LNG popPK model'
    ld1 <- log(0.5); label("Duration of the zero-order infusion into the depot D1 (h)") # NOT PRINTED: Methods 'LNG popPK model' says this duration was estimated but S2 Table omits it; back-solved by the maintainers from the S4 Table week-1 6-h and 12-h means and medians (least squares over 0.1-2 h, optimum 0.5 h); see vignette
    ld2 <- fixed(log(168)); label("Duration of the zero-order input into central D2 (h)") # Stanczyk 2022 Methods 'LNG popPK model': duration of the zero-order absorption into central fixed to 168 h

    # S2 Table reports IIV as CV%; omega^2 = log(1 + CV^2) for Phoenix
    # exponential (log-normal) random effects. No correlations reported.
    etalvc ~ log(1 + 0.389^2) # Stanczyk 2022 S2 Table: IIV Volume 38.9% CV
    etalcl ~ log(1 + 0.444^2) # Stanczyk 2022 S2 Table: IIV Clearance 44.4% CV
    etalfrel ~ log(1 + 0.235^2) # Stanczyk 2022 S2 Table: IIV Frel 23.5% CV

    # Phoenix mixed error Cobs = C + Ceps*sqrt(1 + C^2*(CMultStdev/sigma)^2),
    # i.e. Var = addSd^2 + (propSd*C)^2. The printed multiplicative SD is
    # -0.206; its sign is immaterial because only its square enters.
    propSd <- 0.206; label("Proportional residual error (fraction)") # Stanczyk 2022 S2 Table: multiplicative residual error SD = -0.206
    addSd <- 0.086; label("Additive residual error (ng/mL)") # Stanczyk 2022 S2 Table: additive residual error 0.0860 ug/L (CV 16.5%)
  })
  model({
    cycle3 <- (CYCLE >= 3)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * e_cycle3_cl^cycle3
    vc <- exp(lvc + etalvc) * e_cycle3_vc^cycle3
    frel <- exp(lfrel + etalfrel)
    fbio <- exp(lfdepot)
    d1 <- exp(ld1)
    d2 <- exp(ld2)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # The patch content is dosed to both states: Frel * F is infused into the
    # depot over D1, (1 - Frel) * F enters central as a zero-order input over
    # the wear period.
    f(depot) <- fbio * frel
    dur(depot) <- d1
    f(central) <- fbio * (1 - frel)
    dur(central) <- d2

    # ug / L = ng/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

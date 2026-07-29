Rose_2016_follitropin_delta <- function() {
  description <- paste(
    "One-compartment population PK model for FE 999049 (recombinant human",
    "FSH; INN follitropin delta) with first-order subcutaneous absorption",
    "through a single transit compartment and first-order elimination, in",
    "27 healthy pituitary-suppressed female subjects after a single",
    "subcutaneous dose of 37.5-450 IU (2.2-26.3 ug). Body weight enters as",
    "an allometric covariate on apparent clearance (exponent 0.75) and",
    "apparent volume of distribution (exponent 1) with reference weight",
    "65 kg."
  )
  reference <- paste(
    "Rose TH, Roshammar D, Erichsen L, Grundemar L, Ottesen JT.",
    "Population Pharmacokinetic Modelling of FE 999049, a Recombinant",
    "Human Follicle-Stimulating Hormone, in Healthy Women After Single",
    "Ascending Doses. Drugs R D. 2016 Jun;16(2):173-180.",
    "doi:10.1007/s40268-016-0129-9.",
    sep = " "
  )
  vignette <- "Rose_2016_follitropin_delta"
  units    <- list(time = "hour", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL/F (exponent 0.75) and V/F (exponent 1)",
        "with reference weight 65 kg (Rose 2016 Section 3, typical-value",
        "definition). Body-weight range in the trial was 51.6-90.0 kg",
        "(Table 1)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 27L,
    n_studies      = 1L,
    age_range      = "21-35 years (Table 1, range across all five dose groups)",
    weight_range   = "51.6-90.0 kg (Table 1, range across all five dose groups)",
    weight_median  = "65 kg (reference weight used by the paper for the typical-value parameters)",
    sex_female_pct = 100,
    disease_state  = paste(
      "Healthy pituitary-suppressed female volunteers. All subjects switched",
      "to a single high-dose combined oral contraceptive (OGESTREL 0.5/50,",
      "ethinyl estradiol 50 ug + norgestrel 0.5 mg) 14 days prior to dosing",
      "to suppress endogenous FSH. BMI 18-29 kg/m^2; menstrual cycle",
      "24-35 days; no history of ovarian dysfunction."
    ),
    dose_range     = paste(
      "Single subcutaneous abdominal injection at 37.5, 75, 150, 225, or",
      "450 IU (converted to 2.2, 4.4, 8.8, 13.1, or 26.3 ug using the",
      "drug's specific activity; Section 2.2)."
    ),
    regions        = "Caucasian study (Olsson et al. 2014 reference 8)",
    notes          = paste(
      "First-in-human single ascending dose trial of FE 999049 (INN",
      "follitropin delta, commercial: Rekovelle). The analysis pooled the",
      "5 active dose groups; 27 of the original 40 enrolled subjects were",
      "included after excluding 8 placebo subjects, 2 subjects in the 4.4",
      "ug group with a late endogenous-FSH peak, and 1 subject in the",
      "8.8 ug group with the same pattern. 594 serum FSH samples were",
      "analysed; 258 (43%) were below the LLOQ of 0.075 ug/L and were",
      "handled with the M3 method during estimation (Section 2.3)."
    )
  )

  ini({
    # Structural parameters - typical values for a 65 kg woman (Rose 2016
    # Table 2). The paper reports CL/F and V/F (apparent clearance and
    # apparent volume); F is not separately identified from the SC-only
    # data and is therefore not modelled as a distinct parameter.
    lcl  <- log(0.430); label("Apparent clearance CL/F at WT = 65 kg (L/h)")        # Rose 2016 Table 2 (CL/F = 0.430 L/h, final)
    lvc  <- log(28.0);  label("Apparent volume of distribution V/F at WT = 65 kg (L)") # Rose 2016 Table 2 (V/F = 28.0 L, final)
    lka  <- log(0.160); label("First-order absorption rate from transit to central (1/h)") # Rose 2016 Table 2 (ka = 0.160 1/h, final)
    lktr <- log(0.517); label("First-order transit rate from dosing site to transit compartment (1/h)") # Rose 2016 Table 2 (ktr = 0.517 1/h, final)

    # Allometric exponents on body weight. The paper held these fixed at
    # canonical allometric values (Section 3: "with the power exponent
    # fixed to allometric values").
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL/F (unitless, fixed)") # Rose 2016 Section 3 (fixed allometric exponent)
    e_wt_vc <- fixed(1.00); label("Allometric exponent of WT on V/F (unitless, fixed)")  # Rose 2016 Section 3 (fixed allometric exponent)

    # IIV. Rose 2016 Table 2 reports IIV as CV%; conversion to the
    # log-eta variance scale uses omega^2 = log(CV^2 + 1).
    #   CL/F: 28.2% CV -> log(0.282^2 + 1) = 0.07652
    #   V/F : 44.3% CV -> log(0.443^2 + 1) = 0.17917
    #   ka  : 23.3% CV -> log(0.233^2 + 1) = 0.05286
    # The paper notes (Section 3) a positive correlation between CL/F and
    # V/F was identified in the final model but the magnitude of the
    # covariance is not reported. To avoid fabricating a value, the
    # correlation is omitted here and IIVs are encoded as diagonal; see
    # the vignette's Assumptions and deviations section for context.
    etalcl ~ 0.07652  # Rose 2016 Table 2 (CL/F IIV 28.2% CV)
    etalvc ~ 0.17917  # Rose 2016 Table 2 (V/F IIV 44.3% CV)
    etalka ~ 0.05286  # Rose 2016 Table 2 (ka IIV 23.3% CV)

    # Combined additive-plus-proportional residual error. The proportional
    # term is unitless even though Table 2 prints the column header as
    # "ug/L" -- the additive term is what carries the ug/L units (Section
    # 3 and Table 2 footnote pattern).
    addSd  <- 0.038; label("Additive residual error (ug/L)")                          # Rose 2016 Table 2 (additive 0.038 ug/L, final)
    propSd <- 0.033; label("Proportional residual error (fraction)")                  # Rose 2016 Table 2 (proportional 0.033, final)
  })
  model({
    # Individual PK parameters. Allometric body-weight scaling applies to
    # CL/F and V/F per Rose 2016 Section 3. The reference weight is 65 kg,
    # the typical-value weight reported beneath Table 2.
    cl  <- exp(lcl + etalcl) * (WT / 65)^e_wt_cl
    vc  <- exp(lvc + etalvc) * (WT / 65)^e_wt_vc
    ka  <- exp(lka + etalka)
    ktr <- exp(lktr)
    kel <- cl / vc

    # ODE system reproducing Rose 2016 Equations (2)-(4): dosing-site
    # (depot), transit, and central compartments. The dose enters the
    # depot at rate ktr -> transit -> ka -> central -> elimination at kel.
    # A1 in the paper maps to `depot`, A2 to `transit1`, A3 to `central`.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot - ka * transit1
    d/dt(central)  <-  ka * transit1 - kel * central

    # Observed serum FE 999049 concentration: dose in ug, volume in L ->
    # ug/L which is numerically equal to ng/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

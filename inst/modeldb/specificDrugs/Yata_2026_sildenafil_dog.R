Yata_2026_sildenafil_dog <- function() {
  description <- paste(
    "Veterinary (dog). One-compartment population PK model with first-order absorption",
    "for orally administered sildenafil in 20 client-owned dogs with naturally occurring",
    "pulmonary hypertension (Yata 2026). Dogs had been maintained on 0.8-3.6 mg/kg by",
    "mouth every 8 h for at least 3 days and were sampled sparsely (three plasma samples",
    "per dog, drawn from a randomised allocation of nine nominal times) across a single",
    "8 h dosing interval at steady state; the fit used the steady-state input option in",
    "Phoenix NLME. The disposition is parameterised as the authors fitted it -- an",
    "absorption rate constant Ka, an apparent volume of distribution V/F and an",
    "elimination rate constant Ke, each carrying its own independent exponential",
    "inter-individual random effect -- rather than being reparameterised to the",
    "algebraically equivalent CL/F plus V/F form, which would require a correlated eta",
    "block the authors did not fit. V/F is reported in L/kg, so the model scales it",
    "linearly by body weight; because clinical doses are also prescribed per kg, that",
    "weight dependence cancels out of the predicted concentrations. There is no",
    "intravenous arm, so bioavailability is not identifiable and V/F and the derived CL/F",
    "are apparent values that absorb F. Body weight, age, mg/kg dose, serum creatinine,",
    "alkaline phosphatase, alanine aminotransferase, concurrent diuretic administration",
    "and concurrent pimobendan administration were all screened as covariates and none",
    "met the authors' significance criterion, so the final model carries no covariate",
    "effects. The very large and unexplained inter-individual variability in absorption",
    "(Ka CV 94%) is the paper's central finding. The residual-error magnitude is not",
    "reported anywhere in the paper and is encoded as a zero placeholder."
  )
  reference <- paste(
    "Yata M, DeFrancesco TC, Bonagura JD, Papich MG.",
    "Population Pharmacokinetics of Sildenafil in Dogs With Naturally Occurring",
    "Pulmonary Hypertension.",
    "Journal of Veterinary Pharmacology and Therapeutics. 2026;49:267-273.",
    "doi:10.1111/jvp.70057"
  )
  vignette <- "Yata_2026_sildenafil_dog"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Subject-level, time-fixed. Structural, not a fitted covariate effect: Yata 2026",
        "Table 1 reports V/F in L/kg, so the apparent central volume is the per-kg value",
        "multiplied by body weight. This is the linear (exponent 1) scaling implied by the",
        "reported unit, NOT an estimated allometric exponent -- the paper fits no exponent.",
        "Ke is a rate constant (1/h) and is not weight-scaled, so the apparent clearance",
        "CL/F = Ke * V/F inherits the same linear weight scaling (0.97 L/kg/h).",
        "Body weight was ALSO screened as an additional covariate on Ka, Ke and V/F in the",
        "exploratory step (Yata 2026 Methods 2.6) and was not retained (Results).",
        "Because clinical sildenafil doses are prescribed in mg/kg, dose and volume scale",
        "with WT together and the predicted concentration-time profile is independent of",
        "body weight; WT matters only when an absolute (mg) dose is simulated.",
        "Individual body weights are not tabulated in the paper; dogs weighing < 4 kg were",
        "excluded (Methods 2.2)."
      ),
      source_name = "WT"
    )
  )

  # Covariates Yata 2026 screened in the exploratory covariate step (Methods 2.6)
  # but did NOT retain in the final model: "Exploration of covariates (body weight,
  # dose, creatinine, ALKP, ALT, age, diuretic administration, pimobendan
  # administration) did not identify any significant factors that could explain this
  # variability using our criteria" (Results). No point estimate is published for any
  # of them, so none can be encoded; they are recorded here for provenance only.
  # WT is deliberately absent from this list because it IS referenced in model(),
  # structurally, via the L/kg unit of V/F.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened on Ka, Ke and V/F; not retained. Median 11 years, range 1-16",
        "(Yata 2026 Results)."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = paste(
        "Screened on Ka, Ke and V/F; not retained. Measured on the study day as part of a",
        "complete biochemical panel (Yata 2026 Methods 2.3). The paper reports neither the",
        "unit nor the distribution; mg/dL is the convention of the reporting laboratory",
        "(North Carolina State University clinical pathology) and of US veterinary practice,",
        "and is recorded here so a user supplying this column knows which scale the",
        "screening was done on. No effect is encoded, so the unit is documentation only."
      )
    ),
    ALP = list(
      description = "Serum alkaline phosphatase activity",
      units = "U/L",
      type = "continuous",
      notes = paste(
        "Screened on Ka, Ke and V/F; not retained. Written 'ALKP' in Yata 2026 (Methods 2.6",
        "and Results), the abbreviation used by US veterinary biochemistry panels for",
        "alkaline phosphatase. Included as a marker of hepatic and biliary function because",
        "sildenafil undergoes extensive hepatic first-pass metabolism and the authors",
        "hypothesised that hepatic congestion secondary to pulmonary hypertension might",
        "alter it (Discussion). No distribution is reported."
      )
    ),
    ALT = list(
      description = "Serum alanine aminotransferase activity",
      units = "U/L",
      type = "continuous",
      notes = paste(
        "Screened on Ka, Ke and V/F; not retained. Hepatocellular-injury marker, screened",
        "for the same reason as ALP. No distribution is reported."
      )
    ),
    CONMED_DIURETIC = list(
      description = "Concurrent diuretic administration indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on Ka, Ke and V/F; not retained. Class membership for this paper is",
        "furosemide (n = 2), torsemide (n = 1) and spironolactone (n = 2) among the 20 dogs",
        "(Yata 2026 Results, concurrent-medication list), i.e. loop plus potassium-sparing",
        "agents. Per the CONMED_DIURETIC register entry the pooled class composition is",
        "paper-specific and is enumerated here."
      )
    ),
    CONMED_PIMOBENDAN = list(
      description = "Concurrent pimobendan administration indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on Ka, Ke and V/F; not retained. Pimobendan is an inodilator",
        "(calcium sensitiser plus PDE3 inhibitor) used in canine cardiac disease; it was the",
        "most common concurrent medication, given to 7 of the 20 dogs (Yata 2026 Results).",
        "Screened separately from CONMED_DIURETIC because the authors listed it as its own",
        "binary categorical covariate (Methods 2.6)."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sildenafil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sildenafil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "dog (client-owned, mixed and pure breeds)",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = paste(
      "Three plasma sildenafil samples per dog (60 total) drawn from randomised",
      "allocations of nine nominal times (0, 20 min, 40 min, 1, 1.5, 2, 3, 4, 6, 8 h)",
      "so that 4-7 dogs contributed at each time point; one dog (Dog 5) was sampled at",
      "30 min instead of the scheduled 20 min. All measured concentrations were above",
      "the 5 ng/mL limit of quantitation (Yata 2026 Results)."
    ),
    age_range = "1-16 years",
    age_median = "11 years",
    weight_range = "not reported; dogs weighing < 4 kg were excluded (Methods 2.2)",
    weight_median = "not reported",
    sex_female_pct = 55,
    disease_state = paste(
      "Naturally occurring moderate-to-severe pulmonary hypertension of mixed aetiology,",
      "classified per the 2020 ACVIM consensus statement: undetermined (Group 6, n = 6),",
      "heartworm disease (Group 5, n = 4), lung or airway disease (Group 3, n = 4),",
      "congenital heart disease (Groups 1 and 2, n = 2), mixed heartworm plus lung or",
      "airway disease (Groups 3 and 5, n = 2), chronic lung disease with confirmed",
      "pulmonary thromboembolism (Groups 3 and 4, n = 1), and ACVIM stage D myxomatous",
      "mitral valve disease (Group 2, n = 1). Median tricuspid regurgitation pressure",
      "gradient 76.5 mmHg (range 49-120) at diagnosis in the 18 dogs with a measurable",
      "gradient, and 49.4 mmHg (range 19-117) on sildenafil at the time of the study in",
      "16 dogs."
    ),
    dose_range = paste(
      "0.8-3.6 mg/kg orally every 8 h (Sildenafil 20 mg tablets, Amneal Pharmaceuticals)",
      "for at least 3 days before the study day; median dose on the study day 1.8 mg/kg",
      "with a CV of 32%, median total daily dose 5.2 mg/kg (range 2.4-10.7)."
    ),
    regions = "United States (North Carolina State University Veterinary Hospital, January 2021 to January 2022)",
    sex_breakdown = "9 spayed female, 2 intact female, 7 castrated male, 2 intact male",
    breeds = paste(
      "Mixed breed (n = 7), West Highland White Terrier (n = 2), Pekingese (n = 2),",
      "Shih Tzu (n = 2), and one each of Pembroke Corgi, American Staffordshire Terrier,",
      "Miniature Dachshund, Boxer, Terrier, Miniature Pinscher and Chihuahua."
    ),
    co_medication = paste(
      "Sildenafil was monotherapy in 8 dogs and part of multi-drug therapy in 12.",
      "Concurrent medications were pimobendan (n = 7), enalapril (n = 4),",
      "spironolactone (n = 2), furosemide (n = 2), torsemide (n = 1), prednisolone",
      "(n = 3), doxycycline (n = 2), hydrocodone (n = 2), and one each of fluoxetine,",
      "taurine, theophylline, levothyroxine, diphenhydramine, cetirizine, tacrolimus eye",
      "ointment, grapiprant, atenolol, telmisartan, clopidogrel and carprofen. Drugs",
      "affecting gastrointestinal absorption (H2-receptor blockers, sucralfate, proton",
      "pump inhibitors) or hepatic metabolism (phenobarbital, rifampin, ketoconazole)",
      "were exclusion criteria."
    ),
    notes = paste(
      "Prospective, open-label, steady-state population PK study (IACUC # 20-430); owner",
      "consent obtained for every dog. Dogs were fasted overnight and usual morning",
      "medications were withheld; food and other morning medications were given 3 h after",
      "the sildenafil dose, so the sampled interval is a fasted one. Trazodone",
      "(1.6-5.1 mg/kg) was permitted at home before presentation to reduce stress.",
      "Plasma assayed by LC-MS (m/z 475.3), calibration range 5-1000 ng/mL, LOQ 5 ng/mL.",
      "Median duration of sildenafil therapy before the study was 57 days (range 14-829)."
    )
  )

  ini({
    # ========================================================================
    # Structural model -- Yata 2026 Table 1, 'Estimate' column. The three
    # primary parameters are the typical values (thetas) of a one-compartment
    # model with first-order input, Equation (1):
    #
    #   C(T) = D * Ka / (V/F * (Ka - Ke)) * [exp(-Ke*T) - exp(-Ka*T)]
    #
    # No RSE or confidence interval is published for any of them; Table 1
    # reports only the estimate, the eta shrinkage and the between-subject
    # CV%. None is flagged as fixed, so all three are estimated.
    #
    # The parameterisation is Ka / V/F / Ke exactly as fitted. It is NOT
    # reparameterised to lcl + lvc: Yata 2026 places an independent
    # exponential random effect on each of Ka, V/F and Ke, and the equivalent
    # clearance form would need etalcl = etalkel + etalvc, a correlated eta
    # block the authors did not fit. See the `lkel` register entry, case (b).
    # ========================================================================
    lka  <- log(0.62); label("Absorption rate constant Ka (1/h)")                                                 # Yata 2026 Table 1: theta Ka = 0.62 /h (shrinkage 0.12); Results gives the matching absorption half-life 1.12 h = log(2)/0.62
    lvc  <- log(4.03); label("Apparent central volume of distribution per body weight, V/F = exp(lvc) * WT (L/kg)")  # Yata 2026 Table 1: theta V/F = 4.03 L/kg (shrinkage 0.14). Apparent: no intravenous arm, so V/F absorbs bioavailability (Methods 2.5)
    lkel <- log(0.24); label("Elimination rate constant Ke (1/h)")                                                # Yata 2026 Table 1: theta Ke = 0.24 /h (shrinkage 0.11); Results gives the matching elimination half-life 2.89 h = log(2)/0.24, and CL/F = Ke * V/F = 0.97 L/kg/h reproduces Table 1

    # ========================================================================
    # Inter-individual variability -- Yata 2026 Table 1, 'CV%' column.
    #
    # Equation (2) is the exponential (log-normal) form P_i = theta_P *
    # exp(eta_i,P), with eta ~ N(0, omega^2), so the reported CV% converts to
    # the internal variance scale as omega^2 = log(1 + CV^2).
    #
    # The CV% column is the between-subject variability, not the precision of
    # the estimate: the Abstract ('but variable (CV 94%)'), the Results
    # ('highly variable with CV for Ka of 94%') and the Discussion ('the
    # highest degree of variability was for the parameter of Ka, the
    # absorption rate (CV 94%)') all read it as variability, and it sits
    # beside an eta-shrinkage column, which only a random effect has.
    #
    # All three etas are independent: Equation (2) and the accompanying text
    # state the eta values 'were assumed to be independent'. No covariance is
    # reported, so no block is fitted here.
    # ========================================================================
    etalka  ~ log(1 + 0.9368^2)  # Yata 2026 Table 1: Ka CV% = 93.68 -> omega^2 = 0.6300 (eta shrinkage 0.12)
    etalvc  ~ log(1 + 0.1548^2)  # Yata 2026 Table 1: V/F CV% = 15.48 -> omega^2 = 0.0237 (eta shrinkage 0.14)
    etalkel ~ log(1 + 0.2067^2)  # Yata 2026 Table 1: Ke CV% = 20.67 -> omega^2 = 0.0418 (eta shrinkage 0.11)

    # ========================================================================
    # Residual unexplained variability.
    #
    # The FORM is proportional. Yata 2026 Equation (3) prints
    #   Cobs = Cpred * (1 + epsilon)
    # which is Phoenix NLME's multiplicative residual model. The surrounding
    # prose says 'An additive model described the residual random variability'
    # and 'Cpred is the model predicted concentration plus the error value',
    # which contradicts the displayed equation. The equation is taken as
    # authoritative; it is also the only reading consistent with a data set
    # spanning 5 to 1000 ng/mL on a single proportional error term.
    #
    # The MAGNITUDE is not reported. Table 1 lists no sigma, the paper has no
    # supplementary parameter table (Table S1 is the blood-sampling allocation
    # schedule only), and no residual-error value appears anywhere in the text.
    # Encoded as a zero placeholder rather than invented; simulations with this
    # model therefore return IPRED-level concentrations with between-subject
    # variability but no residual scatter. See the vignette Errata.
    # ========================================================================
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- magnitude not reported in the source)")  # Yata 2026 Equation (3) gives the proportional form Cobs = Cpred * (1 + epsilon) but no value for sigma; Table 1 reports no residual-error row
  })

  model({
    # Individual parameters. Equation (2): exponential IIV on each of the
    # three primary parameters, independently.
    ka  <- exp(lka + etalka)
    # V/F is reported in L/kg (Table 1), so the apparent central volume is the
    # per-kg value multiplied by body weight. Linear scaling implied by the
    # unit, not a fitted allometric exponent.
    vc  <- exp(lvc + etalvc) * WT
    kel <- exp(lkel + etalkel)

    # One-compartment model with first-order input, the ODE form of Equation (1).
    # Bioavailability is not identifiable (oral data only), so no f(depot) term
    # is applied; F is absorbed into the apparent V/F and CL/F instead.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # central is in mg and vc in L, so central/vc is mg/L; the factor of 1000
    # converts to the ng/mL of Table 1 and of the LC-MS assay.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

Faisal_2019_enalapril_odmt <- function() {
  description <- paste(
    "Orodispersible-minitablet-only fit of the Faisal 2019 simultaneous",
    "semi-mechanistic population PK model for enalapril and its active",
    "metabolite enalaprilat in serum and urine (Faisal 2019 Supplementary",
    "Table Supp-II). Same structure as Faisal_2019_enalapril.R -- an Erlang",
    "transit absorption chain, a one-compartment enalapril disposition with",
    "parallel urinary (kurine) and metabolic (kmet) exits, a 2-compartment",
    "Erlang delay on enalaprilat formation, and a two-compartment",
    "enalaprilat disposition cleared only into urine -- but estimated on the",
    "24 subjects' ODMT period alone, with no formulation covariate and with",
    "SIX rather than eight absorption transit compartments, which is the one",
    "structural difference the ODMT data supported. This is the test arm of",
    "the paper's second analysis strategy: the two formulation-specific fits",
    "were run so the individual parameter estimates could be compared by a",
    "paired Wilcoxon signed-rank test, which found a difference only in the",
    "enalapril mean transit time (0.484 h here versus 0.570 h for the",
    "reference tablet, p = 0.03, about 5 min earlier appearance). DOSE",
    "BASIS: supply the dose as the LABELLED mass of enalapril MALEATE (the",
    "trial's 10 mg = 10000 ug), not as enalapril free base."
  )
  reference <- paste(
    "Faisal M, Cawello W, Burckhardt BB, de Hoon J, Laer S; LENA Consortium.",
    "Simultaneous Semi-Mechanistic Population Pharmacokinetic Modeling",
    "Analysis of Enalapril and Enalaprilat Serum and Urine Concentrations",
    "From Child Appropriate Orodispersible Minitablets.",
    "Front Pediatr. 2019;7:281. doi:10.3389/fped.2019.00281",
    "(Supplementary Material, Table Supp-II)"
  )
  vignette <- "Faisal_2019_enalapril"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")
  dosing <- c("transit1")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Figure 1 (model schematic); the
  # transit chain is six compartments long in this formulation-specific fit.
  compartmentData <- list(
    transit1 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "enalapril", units = "ug", specimen = "serum", verified = TRUE),
    urine = list(analyte = "enalapril", units = "ug", specimen = "urine", verified = TRUE),
    transit1_enaat = list(analyte = "enalaprilat", units = "ug", specimen = "administration site", verified = TRUE),
    transit2_enaat = list(analyte = "enalaprilat", units = "ug", specimen = "administration site", verified = TRUE),
    central_enaat = list(analyte = "enalaprilat", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1_enaat = list(analyte = "enalaprilat", units = "ug", specimen = "serum", verified = TRUE),
    urine_enaat = list(analyte = "enalaprilat", units = "ug", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Enters the enalapril central volume of",
        "distribution through Equation 6, TV = theta_TV * (WT_ind /",
        "WT_ref)^theta, with WT_ref = 69.76 kg, the mean body weight of the",
        "analysed cohort (Table 2), and the exponent fixed at 1 per the",
        "Methods statement that theta 'was tested with a fixed value of 1",
        "for the volume of distribution'. Supplementary Table Supp-II is",
        "headed 'FINAL (FULL) MODEL', the same label Table 1 uses for the",
        "covariate-carrying pooled model, and a fixed exponent produces no",
        "table row, so the weight effect is carried here as well. Because",
        "the exponent is fixed at 1 this term is the identity at the",
        "reference weight. Observed weight range 51.8-95.6 kg."
      ),
      source_name = "Weight (kg) (Table 2; Equation 6)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    age_range = "22.08-47.16 years",
    age_mean = "28.00 years",
    age_median = "24.40 years",
    weight_range = "51.8-95.6 kg",
    weight_mean = "69.76 kg",
    weight_median = "67.60 kg",
    height_range = "153.0-189.0 cm",
    total_body_water_range = "32.86-53.70 L",
    race_ethnicity = NULL,
    disease_state = "Healthy adult volunteers",
    dose_range = paste(
      "Single 10 mg oral dose of enalapril maleate as ten 1 mg",
      "child-appropriate orodispersible minitablets (ODMT) with 240 mL of",
      "water -- the test period of the two-treatment, two-period crossover."
    ),
    regions = "Belgium (University Hospitals Leuven / KU Leuven)",
    notes = paste(
      "Demographics from Table 2; parameter estimates from Supplementary",
      "Table Supp-II. Same 24 subjects as the pooled model, restricted to",
      "the ODMT study period. Results: '6 transits were added for the ODMT",
      "formulation (transits = 6)', against eight for the reference and",
      "pooled analyses. NONMEM 7.4.0, ADVAN6, FOCE with interaction;",
      "validated by SIR and a non-parametric bootstrap (Supplementary",
      "Figure Supp-2 carries the goodness-of-fit plots)."
    )
  )

  ini({
    # ================================================================
    # Enalapril (parent) -- Supplementary Table Supp-II,
    # 'FINAL (FULL) MODEL' column, ODMT formulation only.
    # ================================================================
    lka <- log(7.71)
    label("Enalapril absorption rate constant from the last transit (1/h)")
    # Table Supp-II: KA = 7.71 1/h (26.0% RSE)

    lvc <- log(50.70)
    label("Enalapril central volume of distribution at 69.76 kg (L)")
    # Table Supp-II: VC = 50.70 L (5.00% RSE)

    lfdepot <- log(0.589)
    label("Enalapril absolute bioavailability F1 (fraction)")
    # Table Supp-II: F1 = 0.589 (4.00% RSE)

    lmtt <- log(0.484)
    label("Enalapril mean absorption transit time (h)")
    # Table Supp-II: MTT1 = 0.484 h (7.00% RSE). The reference-tablet fit
    # gives 0.570 h; the 0.086 h difference is the paper's headline
    # '5 min early appearance of enalapril from ODMT'.

    lkurine <- log(0.298)
    label("Enalapril urinary excretion rate constant (1/h)")
    # Table Supp-II: KREN = 0.298 1/h (5.00% RSE)

    lkmet <- log(0.693)
    label("Enalaprilat formation rate constant from enalapril (1/h)")
    # Table Supp-II: KM = 0.693 1/h (6.00% RSE)

    # ================================================================
    # Enalaprilat (active metabolite) -- Supplementary Table Supp-II.
    # ================================================================
    lvc_enaat <- log(47.60)
    label("Enalaprilat central volume of distribution (L)")
    # Table Supp-II: VM = 47.60 L (7.00% RSE)

    lk12_enaat <- log(0.060)
    label("Enalaprilat central-to-peripheral rate constant (1/h)")
    # Table Supp-II: KQ1 = 0.060 1/h (5.00% RSE)

    lk21_enaat <- log(0.051)
    label("Enalaprilat peripheral-to-central rate constant (1/h)")
    # Table Supp-II: KQ2 = 0.051 1/h (14.0% RSE)

    lkurine_enaat <- log(0.175)
    label("Enalaprilat urinary excretion rate constant (1/h)")
    # Table Supp-II: KME = 0.175 1/h (6.00% RSE). This is the value the
    # main-text Discussion quotes as 'the estimated value of the rate
    # constant of enalaprilat elimination was 0.175 1/h', although the
    # pooled Table 1 estimate is 0.184 1/h.

    lmtt_enaat <- log(0.873)
    label("Enalaprilat mean formation transit time (h)")
    # Table Supp-II: MTT2 = 0.873 h (11.0% RSE)

    # ================================================================
    # Covariate effect.
    # ================================================================
    e_wt_vc <- fixed(1)
    label("Body-weight exponent on the enalapril central volume (unitless)")
    # Methods, 'Covariate Modeling Analysis', Equation 6: theta 'was tested
    # with a fixed value of 1 for the volume of distribution'.

    # ================================================================
    # Between-subject variability -- Supplementary Table Supp-II,
    # 'Interindividual variability (IIV)' block. Exponential model
    # (Equation 3), so the tabulated values are log-scale variances.
    # ================================================================
    etalka ~ 0.779 # Table Supp-II row 'IIV_KA' = 0.779 (50.0% RSE)
    etalvc ~ 0.047 # Table Supp-II row 'IIV_VC' = 0.047 (34.0% RSE)
    etalfdepot ~ 0.025 # Table Supp-II row 'IIV_F1' = 0.025 (31.0% RSE)
    etalmtt ~ 0.269 # Table Supp-II row 'IIV_MTT1' = 0.269 (30.0% RSE)
    etalkurine ~ 0.056 # Table Supp-II row 'IIV_KREN' = 0.056 (31.0% RSE)
    etalkmet ~ 0.067 # Table Supp-II row 'IIV_KM' = 0.067 (31.0% RSE)
    etalvc_enaat ~ 0.087 # Table Supp-II row 'IIV_VM' = 0.087 (32.0% RSE)
    etalkurine_enaat ~ 0.071 # Table Supp-II row 'IIV_KME' = 0.071 (31.0% RSE)
    etalmtt_enaat ~ 0.094 # Table Supp-II row 'IIV_ MTT2' = 0.094 (31.0% RSE)

    # ================================================================
    # Residual unexplained variability -- Supplementary Table Supp-II.
    # Proportional rows are labelled 'sigma 2' (variances, so the
    # nlmixr2 SD is their square root); additive rows are in ug/L and
    # are already on the SD scale.
    # ================================================================
    propSd <- sqrt(0.010)
    label("Proportional residual SD, enalapril serum (fraction)")
    # Table Supp-II, Serum Enalapril: proportional error variance 0.010
    # (12.0% RSE) -> SD 0.1

    addSd <- 0.186
    label("Additive residual SD, enalapril serum (ug/L)")
    # Table Supp-II, Serum Enalapril: additive error 0.186 ug/L (23.0% RSE)

    propSd_enaat <- sqrt(0.016)
    label("Proportional residual SD, enalaprilat serum (fraction)")
    # Table Supp-II, Serum Enalaprilat: proportional error variance 0.016
    # (13.0% RSE) -> SD 0.1265

    addSd_enaat <- 0.220
    label("Additive residual SD, enalaprilat serum (ug/L)")
    # Table Supp-II, Serum Enalaprilat: additive error 0.220 ug/L (17.0% RSE)

    propSd_urineEna <- sqrt(0.026)
    label("Proportional residual SD, cumulative enalapril in urine (fraction)")
    # Table Supp-II, Urine Enalapril: proportional error variance 0.026
    # (13.0% RSE) -> SD 0.1612

    propSd_urineEnaat <- sqrt(0.005)
    label("Proportional residual SD, cumulative enalaprilat in urine (fraction)")
    # Table Supp-II, Urine Enalaprilat: proportional error variance 0.005
    # (15.0% RSE) -> SD 0.0707
  })

  model({
    # Reference body weight: mean weight of the analysed cohort (Table 2).
    wt_ref <- 69.76

    # Results: '6 transits were added for the ODMT formulation
    # (transits = 6)'; two transits for enalaprilat formation.
    n_transit <- 6
    n_transit_enaat <- 2

    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc
    fdepot <- exp(lfdepot + etalfdepot)
    mtt <- exp(lmtt + etalmtt)
    kurine <- exp(lkurine + etalkurine)
    kmet <- exp(lkmet + etalkmet)
    vc_enaat <- exp(lvc_enaat + etalvc_enaat)
    k12_enaat <- exp(lk12_enaat)
    k21_enaat <- exp(lk21_enaat)
    kurine_enaat <- exp(lkurine_enaat + etalkurine_enaat)
    mtt_enaat <- exp(lmtt_enaat + etalmtt_enaat)

    # MTT = (N + 1) / KTR as printed in the Methods text and Figure 1.
    ktr <- (n_transit + 1) / mtt
    ktr_enaat <- (n_transit_enaat + 1) / mtt_enaat

    # ODE system (Figure 1) with a six-compartment absorption chain.
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ka * transit6

    d/dt(central) <- ka * transit6 - kurine * central - kmet * central
    d/dt(urine) <- kurine * central

    d/dt(transit1_enaat) <- kmet * central - ktr_enaat * transit1_enaat
    d/dt(transit2_enaat) <- ktr_enaat * transit1_enaat - ktr_enaat * transit2_enaat
    d/dt(central_enaat) <- ktr_enaat * transit2_enaat -
      kurine_enaat * central_enaat -
      k12_enaat * central_enaat + k21_enaat * peripheral1_enaat
    d/dt(peripheral1_enaat) <- k12_enaat * central_enaat - k21_enaat * peripheral1_enaat
    d/dt(urine_enaat) <- kurine_enaat * central_enaat

    f(transit1) <- fdepot

    Cc <- central / vc
    Cc_enaat <- central_enaat / vc_enaat
    urineEna <- urine
    urineEnaat <- urine_enaat

    Cc ~ add(addSd) + prop(propSd)
    Cc_enaat ~ add(addSd_enaat) + prop(propSd_enaat)
    urineEna ~ prop(propSd_urineEna)
    urineEnaat ~ prop(propSd_urineEnaat)
  })
}

Faisal_2019_enalapril_reference <- function() {
  description <- paste(
    "Reference-tablet-only fit of the Faisal 2019 simultaneous",
    "semi-mechanistic population PK model for enalapril and its active",
    "metabolite enalaprilat in serum and urine (Faisal 2019 Supplementary",
    "Table Supp-I). Identical structure to Faisal_2019_enalapril.R -- an",
    "8-compartment Erlang transit absorption chain, a one-compartment",
    "enalapril disposition with parallel urinary (kurine) and metabolic",
    "(kmet) exits, a 2-compartment Erlang delay on enalaprilat formation,",
    "and a two-compartment enalaprilat disposition cleared only into urine",
    "-- but estimated on the 24 subjects' reference-formulation period",
    "alone, with no formulation covariate. This is the comparator arm of the",
    "paper's second analysis strategy: the two formulation-specific fits",
    "were run so the individual parameter estimates could be compared by a",
    "paired Wilcoxon signed-rank test, which found a difference only in the",
    "enalapril mean transit time (p = 0.03). DOSE BASIS: supply the dose as",
    "the LABELLED mass of enalapril MALEATE (the trial's 10 mg = 10000 ug),",
    "not as enalapril free base."
  )
  reference <- paste(
    "Faisal M, Cawello W, Burckhardt BB, de Hoon J, Laer S; LENA Consortium.",
    "Simultaneous Semi-Mechanistic Population Pharmacokinetic Modeling",
    "Analysis of Enalapril and Enalaprilat Serum and Urine Concentrations",
    "From Child Appropriate Orodispersible Minitablets.",
    "Front Pediatr. 2019;7:281. doi:10.3389/fped.2019.00281",
    "(Supplementary Material, Table Supp-I)"
  )
  vignette <- "Faisal_2019_enalapril"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")
  dosing <- c("transit1")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Figure 1 (model schematic).
  compartmentData <- list(
    transit1 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "enalapril", units = "ug", specimen = "administration site", verified = TRUE),
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
        "for the volume of distribution'. Supplementary Table Supp-I is",
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
      "Single 10 mg oral dose of enalapril maleate as two 5 mg",
      "market-authorised conventional tablets (Renitec) with 240 mL of",
      "water -- the reference period of the two-treatment, two-period",
      "crossover."
    ),
    regions = "Belgium (University Hospitals Leuven / KU Leuven)",
    notes = paste(
      "Demographics from Table 2; parameter estimates from Supplementary",
      "Table Supp-I. Same 24 subjects as the pooled model, restricted to",
      "the reference-formulation study period. Eight absorption transit",
      "compartments, as in the pooled analysis. NONMEM 7.4.0, ADVAN6,",
      "FOCE with interaction; validated by SIR and a non-parametric",
      "bootstrap (Supplementary Figure Supp-1 carries the goodness-of-fit",
      "plots)."
    )
  )

  ini({
    # ================================================================
    # Enalapril (parent) -- Supplementary Table Supp-I,
    # 'FINAL (FULL) MODEL' column, reference formulation only.
    # ================================================================
    lka <- log(7.780)
    label("Enalapril absorption rate constant from the last transit (1/h)")
    # Table Supp-I: KA = 7.780 1/h (37.0% RSE)

    lvc <- log(51.60)
    label("Enalapril central volume of distribution at 69.76 kg (L)")
    # Table Supp-I: VC = 51.60 L (6.00% RSE)

    lfdepot <- log(0.630)
    label("Enalapril absolute bioavailability F1 (fraction)")
    # Table Supp-I: F1 = 0.630 (5.00% RSE)

    lmtt <- log(0.570)
    label("Enalapril mean absorption transit time (h)")
    # Table Supp-I: MTT1 = 0.570 h (10.00% RSE). Compare the ODMT fit's
    # 0.484 h: the 0.086 h gap is the paper's '5 min early appearance'.

    lkurine <- log(0.312)
    label("Enalapril urinary excretion rate constant (1/h)")
    # Table Supp-I: KREN = 0.312 1/h (5.00% RSE)

    lkmet <- log(0.683)
    label("Enalaprilat formation rate constant from enalapril (1/h)")
    # Table Supp-I: KM = 0.683 1/h (7.00% RSE)

    # ================================================================
    # Enalaprilat (active metabolite) -- Supplementary Table Supp-I.
    # ================================================================
    lvc_enaat <- log(44.70)
    label("Enalaprilat central volume of distribution (L)")
    # Table Supp-I: VM = 44.70 L (5.00% RSE)

    lk12_enaat <- log(0.060)
    label("Enalaprilat central-to-peripheral rate constant (1/h)")
    # Table Supp-I: KQ1 = 0.060 1/h (6.0% RSE)

    lk21_enaat <- log(0.057)
    label("Enalaprilat peripheral-to-central rate constant (1/h)")
    # Table Supp-I: KQ2 = 0.057 1/h (15.0% RSE)

    lkurine_enaat <- log(0.192)
    label("Enalaprilat urinary excretion rate constant (1/h)")
    # Table Supp-I: KME = 0.192 1/h (6.00% RSE)

    lmtt_enaat <- log(0.942)
    label("Enalaprilat mean formation transit time (h)")
    # Table Supp-I: MTT2 = 0.942 h (13.0% RSE)

    # ================================================================
    # Covariate effect.
    # ================================================================
    e_wt_vc <- fixed(1)
    label("Body-weight exponent on the enalapril central volume (unitless)")
    # Methods, 'Covariate Modeling Analysis', Equation 6: theta 'was tested
    # with a fixed value of 1 for the volume of distribution'.

    # ================================================================
    # Between-subject variability -- Supplementary Table Supp-I,
    # 'Interindividual variability (IIV)' block. Exponential model
    # (Equation 3), so the tabulated values are log-scale variances.
    # ================================================================
    etalka ~ 1.310 # Table Supp-I row 'IIV_KA' = 1.310 (53.0% RSE)
    etalvc ~ 0.069 # Table Supp-I row 'IIV_VC' = 0.069 (34.0% RSE)
    etalfdepot ~ 0.057 # Table Supp-I row 'IIV_F1' = 0.057 (31.0% RSE)
    etalmtt ~ 0.203 # Table Supp-I row 'IIV_MTT1' = 0.203 (31.0% RSE)
    etalkurine ~ 0.056 # Table Supp-I row 'IIV_KREN' = 0.056 (33.0% RSE)
    etalkmet ~ 0.088 # Table Supp-I row 'IIV_KM' = 0.088 (32.0% RSE)
    etalvc_enaat ~ 0.048 # Table Supp-I row 'IIV_VM' = 0.048 (35.0% RSE)
    etalkurine_enaat ~ 0.053 # Table Supp-I row 'IIV_KME' = 0.053 (33.0% RSE)
    etalmtt_enaat ~ 0.330 # Table Supp-I row 'IIV_ MTT2' = 0.330 (32.0% RSE)

    # ================================================================
    # Residual unexplained variability -- Supplementary Table Supp-I.
    # Proportional rows are labelled 'sigma 2' (variances, so the
    # nlmixr2 SD is their square root); additive rows are in ug/L and
    # are already on the SD scale.
    # ================================================================
    propSd <- sqrt(0.010)
    label("Proportional residual SD, enalapril serum (fraction)")
    # Table Supp-I, Serum Enalapril: proportional error variance 0.010
    # (12.0% RSE) -> SD 0.1

    addSd <- 0.189
    label("Additive residual SD, enalapril serum (ug/L)")
    # Table Supp-I, Serum Enalapril: additive error 0.189 ug/L (21.0% RSE)

    propSd_enaat <- sqrt(0.021)
    label("Proportional residual SD, enalaprilat serum (fraction)")
    # Table Supp-I, Serum Enalaprilat: proportional error variance 0.021
    # (14.0% RSE) -> SD 0.1449

    addSd_enaat <- 0.220
    label("Additive residual SD, enalaprilat serum (ug/L)")
    # Table Supp-I, Serum Enalaprilat: additive error 0.220 ug/L (22.0% RSE)

    propSd_urineEna <- sqrt(0.011)
    label("Proportional residual SD, cumulative enalapril in urine (fraction)")
    # Table Supp-I, Urine Enalapril: proportional error variance 0.011
    # (14.0% RSE) -> SD 0.1049

    propSd_urineEnaat <- sqrt(0.005)
    label("Proportional residual SD, cumulative enalaprilat in urine (fraction)")
    # Table Supp-I, Urine Enalaprilat: proportional error variance 0.005
    # (16.0% RSE) -> SD 0.0707
  })

  model({
    # Reference body weight: mean weight of the analysed cohort (Table 2).
    wt_ref <- 69.76

    # Results: 'Eight transits (transits = 8) were added for reference and
    # pooled data analysis'; two transits for enalaprilat formation.
    n_transit <- 8
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

    # ODE system (Figure 1). Dose enters transit 1 with bioavailability F1.
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - ktr * transit7
    d/dt(transit8) <- ktr * transit7 - ka * transit8

    d/dt(central) <- ka * transit8 - kurine * central - kmet * central
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

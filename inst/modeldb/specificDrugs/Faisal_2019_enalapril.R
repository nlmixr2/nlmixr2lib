Faisal_2019_enalapril <- function() {
  description <- paste(
    "Simultaneous semi-mechanistic population PK model for the prodrug",
    "enalapril and its active diacid metabolite enalaprilat, fitted jointly",
    "to serum concentrations AND cumulative urinary excretion amounts of both",
    "analytes in 24 healthy adults given a single 10 mg oral dose of enalapril",
    "maleate (Faisal 2019, LENA consortium). Enalapril absorption is an",
    "8-compartment Erlang transit chain (rate ktr) feeding the central",
    "compartment at rate ka; enalapril leaves the central compartment by two",
    "parallel first-order routes, direct urinary excretion (kurine) and",
    "metabolic conversion to enalaprilat (kmet). The enalaprilat formation",
    "phase is delayed by a 2-compartment Erlang transit chain (rate",
    "ktr_enaat) before entering a two-compartment enalaprilat disposition",
    "model that is eliminated only by urinary excretion (kurine_enaat).",
    "Cumulative urinary amounts of both analytes are carried as explicit",
    "excretion compartments, which is what makes the absolute bioavailability",
    "F1 identifiable (estimated 0.606). Body weight scales the enalapril",
    "central volume with the exponent FIXED at 1 and referenced to the study",
    "mean 69.76 kg; the orodispersible-minitablet formulation multiplies the",
    "enalapril mean transit time by 0.730 relative to the reference tablet.",
    "DOSE BASIS: supply the dose as the LABELLED mass of enalapril MALEATE",
    "(the trial's 10 mg = 10000 ug), not as enalapril free base. This differs",
    "from Steichert_2025_enalapril_enalaprilat_pediatric.R, which requires",
    "free base; to move a dose between the two models multiply or divide by",
    "376.45/492.52 = 0.76433."
  )
  reference <- paste(
    "Faisal M, Cawello W, Burckhardt BB, de Hoon J, Laer S; LENA Consortium.",
    "Simultaneous Semi-Mechanistic Population Pharmacokinetic Modeling",
    "Analysis of Enalapril and Enalaprilat Serum and Urine Concentrations",
    "From Child Appropriate Orodispersible Minitablets.",
    "Front Pediatr. 2019;7:281. doi:10.3389/fped.2019.00281"
  )
  vignette <- "Faisal_2019_enalapril"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")
  dosing <- c("transit1")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Figure 1 (model schematic) and the
  # Methods 'Population Pharmacokinetic Model Structure' section.
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
        "analysed cohort (Table 2). The paper states the exponent theta 'was",
        "tested with a fixed value of 1 for the volume of distribution and",
        "0.75 for clearance'; weight was retained only on a volume, so the",
        "exponent is encoded as fixed(1). Adding this covariate dropped the",
        "objective function by 18.2 and backward elimination increased it",
        "significantly, so it was retained in the final model (Results,",
        "'Model Evaluation Results'). Observed weight range 51.8-95.6 kg."
      ),
      source_name = "Weight (kg) (Table 2; Equation 6)"
    ),
    FORM_ODMT = list(
      description = paste(
        "Orodispersible minitablet (ODMT) formulation indicator:",
        "1 = the dose was given as child-appropriate orodispersible",
        "minitablets, 0 = the market-authorised reference tablet."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference tablet, Renitec, 2 x 5 mg enalapril maleate)",
      notes = paste(
        "Per-dose-record indicator. In this trial the ODMT arm was 10",
        "minitablets of 1 mg strength and the reference arm was two 5 mg",
        "Renitec tablets; both were swallowed with 240 mL of water, so the",
        "contrast is the multiple-minitablet presentation rather than",
        "in-mouth disintegration. Enters the enalapril mean transit time",
        "multiplicatively through Equation 7, TV = theta_X * theta_12^FORM *",
        "exp(eta), with theta_12 = 0.730, i.e. the ODMT shortens MTT1 by",
        "27% relative to the reference tablet. Adding the effect dropped the",
        "objective function by 6.51 and it was retained in the final model.",
        "No other model parameter carried a formulation effect."
      ),
      source_name = "FORM (Equation 7; Results, 'Pharmacokinetics Comparison of ODMT and Reference Formulation')"
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
      "Single 10 mg oral dose of enalapril maleate in each of two periods",
      "(two-treatment, two-period crossover). Reference treatment: two 5 mg",
      "market-authorised conventional tablets (Renitec) with 240 mL water.",
      "Test treatment: ten 1 mg child-appropriate orodispersible minitablets",
      "(ODMT) with 240 mL water."
    ),
    regions = "Belgium (University Hospitals Leuven / KU Leuven)",
    notes = paste(
      "Demographics from Table 2. Phase I relative-bioavailability trial",
      "NCT02252692 within the EU FP7 LENA project (Labeling of Enalapril",
      "from Neonate to Adolescence, grant 602295). 2,208 serum and urine",
      "concentrations. Serum sampling at 0.17, 0.33, 0.5, 0.75, 1, 1.25,",
      "1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 12, 24 and 48 h; urine",
      "collected over 0-2, 2-4, 4-8, 8-12, 12-24, 24-36 and 36-48 h and",
      "converted to cumulative excreted amounts. NONMEM 7.4.0, ADVAN6,",
      "FOCE with interaction. Serum LLOQ 0.195 ng/mL (enalapril) and",
      "0.180 ng/mL (enalaprilat); samples below LLOQ were excluded.",
      "Validated by VPC, a 200-sample non-parametric bootstrap and a",
      "sampling-importance-resampling (SIR) uncertainty analysis; all eta",
      "shrinkages were below 25%."
    )
  )

  ini({
    # ================================================================
    # Enalapril (parent) disposition and absorption -- Table 1,
    # 'FINAL (FULL) MODEL' column, pooled two-formulation analysis.
    # ================================================================
    lka <- log(6.010)
    label("Enalapril absorption rate constant from the last transit (1/h)")
    # Table 1: KA = 6.010 1/h (15.0% RSE); Figure 1 shows KA as the transfer
    # from transit 8 into the enalapril central compartment.

    lvc <- log(51.10)
    label("Enalapril central volume of distribution at 69.76 kg (L)")
    # Table 1: VC = 51.10 L (4.00% RSE).

    lfdepot <- log(0.606)
    label("Enalapril absolute bioavailability F1 (fraction)")
    # Table 1: F1 = 0.606 (3.00% RSE). Identifiable because the cumulative
    # urinary amounts of both analytes are fitted; the paper notes the
    # estimate matches the published 60% absorbed fraction.

    lmtt <- log(0.558)
    label("Enalapril mean absorption transit time, reference tablet (h)")
    # Table 1: MTT1 = 0.558 h (9.00% RSE). This is the FORM = 0 (reference
    # tablet) value; the ODMT value is obtained via e_form_odmt_mtt.

    lkurine <- log(0.305)
    label("Enalapril urinary excretion rate constant (1/h)")
    # Table 1: KREN = 0.305 1/h (4.00% RSE).

    lkmet <- log(0.688)
    label("Enalaprilat formation rate constant from enalapril (1/h)")
    # Table 1: KM = 0.688 1/h (4.00% RSE).

    # ================================================================
    # Enalaprilat (active metabolite) -- Table 1.
    # ================================================================
    lvc_enaat <- log(46.10)
    label("Enalaprilat central volume of distribution (L)")
    # Table 1: VM = 46.10 L (4.00% RSE).

    lk12_enaat <- log(0.060)
    label("Enalaprilat central-to-peripheral rate constant (1/h)")
    # Table 1: KQ1 = 0.060 1/h (4.00% RSE). Table 1 footnote defines KQ1 as
    # 'Rate constant of enalaprilat distribution from central to peripheral
    # compartment'.

    lk21_enaat <- log(0.054)
    label("Enalaprilat peripheral-to-central rate constant (1/h)")
    # Table 1: KQ2 = 0.054 1/h (10.0% RSE). Table 1 footnote defines KQ2 as
    # 'Rate constant of enalaprilat distribution from peripheral to central
    # compartment'.

    lkurine_enaat <- log(0.184)
    label("Enalaprilat urinary excretion rate constant (1/h)")
    # Table 1: KME = 0.184 1/h (4.00% RSE). Figure 1 labels the same arrow
    # KMEL. The Discussion quotes 0.175 1/h, which is the ODMT-only value
    # from Supplementary Table Supp-II, not this pooled estimate.

    lmtt_enaat <- log(0.910)
    label("Enalaprilat mean formation transit time (h)")
    # Table 1: MTT2 = 0.910 h (8.00% RSE).

    # ================================================================
    # Covariate effects.
    # ================================================================
    e_wt_vc <- fixed(1)
    label("Body-weight exponent on the enalapril central volume (unitless)")
    # Methods, 'Covariate Modeling Analysis', Equation 6: 'The parameter
    # theta was tested with a fixed value of 1 for the volume of
    # distribution and 0.75 for clearance.' Weight was retained only on VC,
    # a volume, so the exponent is fixed at 1 and carries no RSE in Table 1.

    e_form_odmt_mtt <- 0.730
    label("Multiplicative effect of the ODMT formulation on MTT1 (unitless)")
    # Table 1: 'THETA (X)' = 0.730 (12.0% RSE), read as theta_12 of
    # Equation 7 (TV = theta_X * theta_12^FORM * exp(eta)). See the vignette
    # Errata for the three lines of evidence supporting this reading.

    # ================================================================
    # Between-subject variability -- Table 1, 'INTERINDIVIDUAL
    # VARIABILITY (IIV)' block. Exponential model (Equation 3,
    # Pi = TVp * exp(ETAi)), so the tabulated values are variances of
    # the eta on the log scale. No IIV was estimated on KQ1 or KQ2.
    # ================================================================
    etalka ~ 0.688 # Table 1 row 'IIV_KA' = 0.688 (31.0% RSE)
    etalvc ~ 0.058 # Table 1 row 'IIV_VC' = 0.058 (24.0% RSE)
    etalfdepot ~ 0.041 # Table 1 row 'IIV_F1' = 0.041 (22.0% RSE)
    etalmtt ~ 0.151 # Table 1 row 'IIV_MTT1' = 0.151 (22.0% RSE)
    etalkurine ~ 0.058 # Table 1 row 'IIV_KREN' = 0.058 (24.0% RSE)
    etalkmet ~ 0.078 # Table 1 row 'IIV_KM' = 0.078 (22.0% RSE)
    etalvc_enaat ~ 0.069 # Table 1 row 'IIV_VM' = 0.069 (23.0% RSE)
    etalkurine_enaat ~ 0.063 # Table 1 row 'IIV_KME' = 0.063 (23.0% RSE)
    etalmtt_enaat ~ 0.296 # Table 1 row 'IIV_ MTT2' = 0.296 (22.0% RSE)

    # ================================================================
    # Residual unexplained variability -- Table 1, 'RESIDUAL
    # UNEXPLAINED VARIABILITY (RUV)' block. The proportional rows are
    # labelled 'Proportional error (sigma 2)', i.e. variances, so the
    # nlmixr2 SD is their square root. The additive rows are labelled
    # 'Additive error (ug/l)', i.e. already on the SD scale (and both
    # land at the respective serum LLOQ, 0.195 and 0.180 ng/mL).
    # ================================================================
    propSd <- sqrt(0.010)
    label("Proportional residual SD, enalapril serum (fraction)")
    # Table 1, Serum Enalapril: proportional error variance 0.010 (8.00%
    # RSE) -> SD 0.1, a 10% CV.

    addSd <- 0.188
    label("Additive residual SD, enalapril serum (ug/L)")
    # Table 1, Serum Enalapril: additive error 0.188 ug/L (15.0% RSE).

    propSd_enaat <- sqrt(0.018)
    label("Proportional residual SD, enalaprilat serum (fraction)")
    # Table 1, Serum Enalaprilat: proportional error variance 0.018 (9.00%
    # RSE) -> SD 0.1342.

    addSd_enaat <- 0.220
    label("Additive residual SD, enalaprilat serum (ug/L)")
    # Table 1, Serum Enalaprilat: additive error 0.220 ug/L (13.0% RSE).

    propSd_urineEna <- sqrt(0.019)
    label("Proportional residual SD, cumulative enalapril in urine (fraction)")
    # Table 1, Urine Enalapril: proportional error variance 0.019 (9.00%
    # RSE) -> SD 0.1378.

    propSd_urineEnaat <- sqrt(0.005)
    label("Proportional residual SD, cumulative enalaprilat in urine (fraction)")
    # Table 1, Urine Enalaprilat: proportional error variance 0.005 (11.0%
    # RSE) -> SD 0.0707.
  })

  model({
    # Reference body weight for the allometric term: the mean weight of the
    # analysed cohort (Table 2), which is what Equation 6 normalises by.
    wt_ref <- 69.76

    # Erlang transit-chain lengths. Results: 'Eight transits (transits = 8)
    # were added for reference and pooled data analysis'; two transits were
    # added for the enalaprilat formation phase (Figure 1, N = 8 and N = 2).
    n_transit <- 8
    n_transit_enaat <- 2

    # Individual parameters. Equation 3 gives the exponential BSV model;
    # Equation 6 the weight effect on VC; Equation 7 the formulation effect
    # on MTT1, written exactly as the paper does (a power of the binary
    # indicator, so FORM_ODMT = 0 recovers the reference-tablet value).
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc
    fdepot <- exp(lfdepot + etalfdepot)
    mtt <- exp(lmtt + etalmtt) * e_form_odmt_mtt^FORM_ODMT
    kurine <- exp(lkurine + etalkurine)
    kmet <- exp(lkmet + etalkmet)
    vc_enaat <- exp(lvc_enaat + etalvc_enaat)
    k12_enaat <- exp(lk12_enaat)
    k21_enaat <- exp(lk21_enaat)
    kurine_enaat <- exp(lkurine_enaat + etalkurine_enaat)
    mtt_enaat <- exp(lmtt_enaat + etalmtt_enaat)

    # Transit rate constants. The paper prints MTT1 = (N + 1) / KTR in both
    # the Methods text and Figure 1, and estimates KA separately as the
    # transfer out of the last transit; the (N + 1) convention is therefore
    # used as written rather than the bare N / MTT form.
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

    # Enalapril central compartment: two parallel first-order exits, direct
    # urinary excretion and conversion to enalaprilat. The paper assumes no
    # other elimination route for either analyte.
    d/dt(central) <- ka * transit8 - kurine * central - kmet * central
    d/dt(urine) <- kurine * central

    # Enalaprilat formation delay (two Erlang transits) then two-compartment
    # disposition with urinary excretion only.
    d/dt(transit1_enaat) <- kmet * central - ktr_enaat * transit1_enaat
    d/dt(transit2_enaat) <- ktr_enaat * transit1_enaat - ktr_enaat * transit2_enaat
    d/dt(central_enaat) <- ktr_enaat * transit2_enaat -
      kurine_enaat * central_enaat -
      k12_enaat * central_enaat + k21_enaat * peripheral1_enaat
    d/dt(peripheral1_enaat) <- k12_enaat * central_enaat - k21_enaat * peripheral1_enaat
    d/dt(urine_enaat) <- kurine_enaat * central_enaat

    f(transit1) <- fdepot

    # Observations. The two serum outputs are concentrations; the two urine
    # outputs are the cumulative amounts excreted, which is the quantity the
    # paper fitted after converting measured urine concentrations.
    Cc <- central / vc
    Cc_enaat <- central_enaat / vc_enaat
    urineEna <- urine
    urineEnaat <- urine_enaat

    # Equations 4 and 5: combined additive plus proportional error on the
    # serum outputs, proportional only on the urine outputs.
    Cc ~ add(addSd) + prop(propSd)
    Cc_enaat ~ add(addSd_enaat) + prop(propSd_enaat)
    urineEna ~ prop(propSd_urineEna)
    urineEnaat ~ prop(propSd_urineEnaat)
  })
}

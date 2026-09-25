Li_2018_guanfacine <- function() {
  description <- paste(
    "One-compartment oral population pharmacokinetic reduction of the",
    "Simcyp (version 14) physiologically based pharmacokinetic model",
    "for guanfacine extended-release (GXR) in healthy adults, with the",
    "CYP3A4 drug-drug-interaction layer carried as coadministered-",
    "perpetrator covariates (Li 2018). The platform model itself cannot",
    "be encoded here - the virtual North European Caucasian population",
    "file, the per-organ physiology and the perpetrator compound files",
    "are Simcyp database content and are not published. What IS fully",
    "reported is the guanfacine compound layer (Table 1 plus the",
    "derivation chain in Section 2.2.1), and it is sufficient to",
    "reconstruct the disposition as an ordinary compartmental model",
    "with NO fitted parameters: first-order absorption into a depot,",
    "a single distribution volume from the reported Vss of 8.0 L/kg,",
    "and additive renal plus CYP3A4-mediated non-renal elimination",
    "whose reported 12.6 and 12.2 L/h sum to the observed intravenous",
    "clearance of 24.8 L/h. Oral bioavailability is not a free",
    "parameter either: the gut and hepatic first-pass availabilities",
    "are pinned by the paper's own well-stirred arithmetic and their",
    "product reproduces the observed CL/F of 37.8 L/h exactly.",
    "The reduction reproduces the paper's own predicted Cmax to within",
    "8.6% at both simulated doses, with Tmax 5.6 h against an observed",
    "6 h and a terminal half-life of 15.7 h against an observed 17 h.",
    "Coadministered CYP3A4 inhibitors and inducers act through a single",
    "relative CYP3A4 activity term that scales gut first-pass",
    "extraction, hepatic first-pass extraction and non-renal clearance",
    "together. Each perpetrator's activity was back-solved from its",
    "published AUC ratio ALONE; the published Cmax ratio was held out",
    "and is then reproduced within 10% for all six arms (within 5% for",
    "five of them), which is out-of-sample evidence that the reduction",
    "has the right shape and not merely the right exposure totals.",
    "This is a typical-value simulation model: the source is a PBPK",
    "analysis and reports no inter-individual variance components and",
    "no residual-error model, so there are no etas and propSd is fixed",
    "at zero. It is an ADULT model - the source states explicitly that",
    "no pediatric or adolescent pharmacokinetic data were used in model",
    "development - so no weight scaling is carried; see the vignette.",
    sep = " "
  )
  reference <- paste(
    "Li A, Yeo K, Welty D, Rong H. (2018). Development of Guanfacine",
    "Extended-Release Dosing Strategies in Children and Adolescents",
    "with ADHD Using a Physiologically Based Pharmacokinetic Model to",
    "Predict Drug-Drug Interactions with Moderate CYP3A4 Inhibitors or",
    "Inducers. Paediatr Drugs 20(1):19-28.",
    "doi:10.1007/s40272-017-0270-0.",
    sep = " "
  )
  vignette <- "Li_2018_guanfacine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in
  # what biological matrix. Verified against Li 2018 Figure 1b (the
  # Simcyp PBPK schematic: an absorption site feeding the gut wall and
  # portal vein, then the liver, then the systemic circulation) and
  # Section 2.2.1, which describes first-order absorption with a single
  # steady-state distribution volume.
  compartmentData <- list(
    depot = list(
      analyte = "guanfacine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "guanfacine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  # The perpetrators Li 2018 simulates. All flags zero (and
  # DOSE_EFAVIRENZ_MG zero) is GXR monotherapy. Each coefficient in
  # ini() is the log relative CYP3A4 activity back-solved from that
  # perpetrator's published AUC geometric mean ratio.
  covariateData <- list(
    CONMED_KETOCONAZOLE = list(
      description = "Concomitant ketoconazole 400 mg once daily (strong CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 perpetrator coadministered (GXR monotherapy)",
      notes = paste(
        "Li 2018 Table 2. GXR 4 mg was given on day 3 of 6 days of",
        "ketoconazole dosing. This arm has observed clinical data",
        "(Shire study SPD503-106, n = 20) as well as a prediction; the",
        "coefficient here is back-solved from the PREDICTED AUC ratio",
        "of 2.56, because it is the model that is being reproduced.",
        "Perpetrator inhibition constants were the Simcyp version 14",
        "compound-file defaults and are not published, so the",
        "coefficient is an empirical relative activity rather than a",
        "transcribed Ki or kinact."
      ),
      source_name = "ketoconazole"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampicin 600 mg once daily (strong CYP3A4 inducer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 perpetrator coadministered (GXR monotherapy)",
      notes = paste(
        "Li 2018 Table 2. GXR 4 mg was given on day 8 of 11 days of",
        "rifampicin dosing (Shire study SPD503-108, n = 20). The source",
        "reports two candidate rifampicin models, differing only in the",
        "maximum fold induction Indmax. The coefficient here is",
        "back-solved from the Indmax = 8 prediction (AUC ratio 0.35),",
        "which is the model the authors adopted: Section 3.2 and",
        "Figure 5 report that Indmax = 16 overpredicted induction",
        "(AUC ratio 0.16 against an observed 0.31). The rejected",
        "Indmax = 16 variant is reproduced in the vignette but is not",
        "encoded here."
      ),
      source_name = "rifampicin"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole 200 mg once daily (moderate CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 perpetrator coadministered (GXR monotherapy)",
      notes = paste(
        "Li 2018 Table 3. GXR 4 mg was given on day 3 of 6 days of",
        "fluconazole dosing, after a 400 mg loading dose on day 1.",
        "Prospective prediction only - no clinical drug-drug-interaction",
        "study has been conducted for this pair."
      ),
      source_name = "fluconazole"
    ),
    CONMED_ERYTHROMYCIN = list(
      description = "Concomitant erythromycin 500 mg three times daily (moderate CYP3A4 inhibitor)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 perpetrator coadministered (GXR monotherapy)",
      notes = paste(
        "Li 2018 Table 3. GXR 4 mg was given on day 3 of 6 days of",
        "erythromycin dosing. Prospective prediction only."
      ),
      source_name = "erythromycin"
    ),
    CONMED_EFV = list(
      description = "Concomitant efavirenz (moderate CYP3A4 inducer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no CYP3A4 perpetrator coadministered (GXR monotherapy)",
      notes = paste(
        "Li 2018 Table 3. GXR 4 mg was given on day 10 of 14 days of",
        "efavirenz dosing. Here the reference category is simply 'no",
        "efavirenz' rather than the alternative-antiretroviral-regimen",
        "reference used by the antiretroviral population-PK models that",
        "share this column. Set DOSE_EFAVIRENZ_MG alongside this flag",
        "to select which of the two simulated efavirenz dose levels",
        "applies; this flag alone carries the 400 mg effect.",
        "Prospective prediction only."
      ),
      source_name = "efavirenz"
    ),
    DOSE_EFAVIRENZ_MG = list(
      description = "Daily dose of coadministered efavirenz (mg)",
      units = "mg",
      type = "continuous",
      reference_category = "0 = no efavirenz coadministered",
      notes = paste(
        "Li 2018 simulated two efavirenz dose levels, 400 mg and 600 mg",
        "once daily (Section 2.5, Table 3), and reports a materially",
        "different guanfacine AUC ratio for each (0.58 and 0.33). The",
        "source gives no dose-response function linking them, so this",
        "column is used only to select between the two published",
        "levels: a value of 500 mg or more selects the 600 mg",
        "coefficient and any lower nonzero value selects the 400 mg",
        "coefficient. Do not interpolate. The 600 mg model is the one",
        "published by Ke et al.; the 400 mg model was built from it by",
        "adjusting efavirenz oral clearance to 17.7 L/h."
      ),
      source_name = "efavirenz dose"
    )
  )

  # Reported by the source but deliberately not carried.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Li 2018 Table 1 gives Vss in L/kg, so the Simcyp",
        "model does scale distribution volume with weight, and this",
        "reduction fixes it at the 70 kg reference instead. Weight is",
        "NOT carried as a covariate because the source publishes no",
        "weight relationship for clearance - Simcyp derives hepatic",
        "clearance from liver weight and microsomal protein content",
        "that are population-file content - so scaling volume alone",
        "would shorten the half-life in a lighter subject with no",
        "support from the paper. Section 4 states explicitly that no",
        "pediatric pharmacokinetic data were used in model development",
        "and that the weight-scaling claim rests on a separate",
        "population-PK analysis (reference 18) that is not reproduced",
        "here. Treat this as an adult model."
      ),
      units = "kg",
      type = "continuous",
      notes = "Reference weight 70 kg assumed; the Simcyp population mean is not published."
    ),
    SEXF = list(
      description = paste(
        "Sex. Each simulated trial was sex-matched to its comparator",
        "clinical study (46.2%, 65%, 40% and 50% women across the",
        "arms), but Li 2018 reports no sex effect on any guanfacine",
        "parameter; sex enters only through the Simcyp population file's",
        "demographic sampling."
      ),
      units = "(binary)",
      type = "binary",
      notes = "Trial-matching input to the virtual population, not a model covariate."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 3L,
    age_range = "18-54 years (model-development study); 18-53 years across the validation studies",
    sex_female_pct = 46.2,
    disease_state = paste(
      "Healthy adults. The paper's dosing recommendations are for",
      "children and adolescents aged 6-17 years with",
      "attention-deficit/hyperactivity disorder, but no pediatric data",
      "were used in model development or validation (Section 4,",
      "limitation 1)."
    ),
    dose_range = "Single oral doses of guanfacine extended-release 2 mg and 4 mg",
    regions = "North European Caucasian virtual population (Simcyp version 14 default, Howgate 2006), demographics matched per trial",
    studies = paste(
      "Swearingen et al. (reference 19), a randomized crossover trial",
      "in 52 healthy adults of whom 49 completed, each receiving a",
      "single oral dose of GXR 2 mg or 4 mg once a week for 4 weeks",
      "with sampling to 96 h - this study supplied the observed Cmax,",
      "AUC and CL/F the compound file was optimized against. Shire",
      "study SPD503-106 (n = 20, 19-50 years, 65% women) supplied the",
      "ketoconazole drug-drug-interaction data and SPD503-108 (n = 20,",
      "18-53 years, 40% women) the rifampicin data; both are data on",
      "file with their summary statistics printed in Table 2. Renal",
      "and total intravenous clearance came from a separate",
      "intravenous study (reference 20)."
    ),
    notes = paste(
      "n_subjects records the 49 completers of the Swearingen study,",
      "the data the disposition parameters were optimized against;",
      "n_studies counts the three clinical studies above. This is a",
      "PBPK analysis rather than a population-PK fit, so there is no",
      "pooled analysis dataset and no estimated variance components.",
      "The between-trial spread visible in Figure 2 is Simcyp virtual-",
      "population output driven by demographic and CYP3A4-abundance",
      "distributions that are not published, so it is not encoded as",
      "an omega here."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter is fixed: this is a typical-value simulation
    # model. Nothing below was fitted to anything. Each value is either
    # a verbatim Li 2018 Table 1 input, an arithmetic consequence of the
    # derivation chain the paper prints in Section 2.2.1, or (for the
    # perpetrator terms) back-solved from a published AUC ratio.
    #
    # The Section 2.2.1 chain, reproduced here so the arithmetic below
    # can be checked line by line:
    #   CL_IV  = 24.8 L/h   total clearance after intravenous dosing
    #   CL_R   = 12.6 L/h   renal clearance, 'estimated to be 50% of
    #                       total body clearance'
    #   B:P    = 1.45       blood-to-plasma ratio
    #   Q_H    = 90 L/h     hepatic blood flow
    #   CL/F   = 37.8 L/h   observed apparent oral clearance
    #   CL_H,B = (24.8 - 12.6) / 1.45 = 8.4138 L/h   (paper prints 8.41)
    #   F_H    = 1 - 8.4138 / 90     = 0.90651       (paper prints 0.91)
    # Both printed roundings match, which is the confirmation that the
    # transcription of CL_IV, CL_R, B:P and Q_H is correct.
    # ------------------------------------------------------------------

    # -- Absorption ----------------------------------------------------
    lka <- fixed(log(0.465))
    label("First-order absorption rate constant ka (1/h)")
    # Li 2018 Table 1, row 'First-order absorption rate constant (ka)'
    # = 0.465 /h, footnoted as derived from in vivo data. Section 2.2.1
    # states 0.459 /h for the same quantity; the discrepancy is real in
    # the published PDF and is not a conversion artifact. Table 1 is the
    # compound-file input table, so 0.465 is the value the simulations
    # actually used and is the one taken here. The choice is immaterial
    # - it moves Tmax from 5.59 h to 5.64 h - but it is recorded rather
    # than silently resolved. The authors obtained it by sensitivity
    # analysis, as the value recovering the observed Tmax of about 6 h.

    fa <- fixed(1)
    label("Fraction of the dose absorbed")
    # Li 2018 Table 1, row 'Fraction absorbed (fa)' = 1, from Caco-2
    # permeability of 26e-6 cm/s (data on file, V00977-SPD503).

    # -- Distribution --------------------------------------------------
    lvc <- fixed(log(560))
    label("Central volume of distribution Vc (L)")
    # Li 2018 Table 1, row 'Volume of distribution at steady state
    # (Vss)' = 8.0 L/kg, times the 70 kg reference weight assumed here
    # (the Simcyp population mean weight is not published). The model
    # has a single distribution volume, so Vc = Vss. Table 1 footnote a
    # records that 8.0 L/kg was chosen to recover the observed Cmax
    # values of 1.57 and 3.58 ng/mL, and that it is consistent with the
    # Rodgers and Rowland prediction of 9.08 L/kg and with in vivo Vss
    # values of 6.3-8.6 L/kg.

    # -- Elimination ---------------------------------------------------
    lcl_renal <- fixed(log(12.6))
    label("Renal clearance CL_R (L/h)")
    # Li 2018 Table 1, row 'Mean renal clearance (CL_R)' = 12.6 L/h,
    # from observed intravenous data (reference 20). Not affected by
    # CYP3A4 perpetrators.

    lcl_nonren <- fixed(log(12.2))
    label("Non-renal (CYP3A4-mediated) clearance CL_met (L/h)")
    # Section 2.2.1: CL_IV = 24.8 L/h with 'equal contributions from
    # metabolic and renal clearances', so CL_met = 24.8 - 12.6 = 12.2
    # L/h. Section 4 confirms the metabolic share is attributed to
    # CYP3A4 ('a CYP3A4 fraction of metabolism of 50%'). This is the
    # only clearance term the perpetrator covariates act on. Carrying
    # the split explicitly rather than a single 24.8 L/h matters: a
    # complete CYP3A4 inhibitor can at most halve total clearance,
    # which is what bounds the attainable AUC ratio.

    # -- First-pass extraction -----------------------------------------
    # The two extraction ratios below are fully determined by the
    # paper's own numbers; neither is fitted.
    #   F   = CL_IV / (CL/F) = 24.8 / 37.8 = 0.656085
    #   F_H = 0.906513 (chain above), so E_H = 0.093487
    #   F_G = F / (fa * F_H) = 0.656085 / 0.906513 = 0.723745,
    #         so E_gut = 0.276255
    # F_G is the quantity Section 2.2.1 says Q_G was refined to 1 L/h to
    # produce, via Equation 1 and F = fa * F_G * F_H. Equation 1 does
    # not survive the trim as text, but it is the standard well-stirred
    # gut form F_G = Q_G / (Q_G + fu_G * CLu_G,int) with fu_G = 1, and
    # the extraction-ratio parameterisation used here is algebraically
    # identical to it while avoiding the unpublished CLu_G,int.
    eh <- fixed(0.0934865900)
    label("Baseline hepatic first-pass extraction ratio")

    egut <- fixed(0.2762548807)
    label("Baseline gut-wall first-pass extraction ratio")

    # -- Coadministered CYP3A4 perpetrator effects ---------------------
    # Each coefficient is log(relative CYP3A4 activity) for that
    # perpetrator, back-solved from its published AUC geometric mean
    # ratio ALONE. Activity below 1 is inhibition, above 1 induction.
    # The published Cmax ratio was deliberately held out; the model then
    # reproduces it as an out-of-sample check (vignette Table 5):
    #   ketoconazole  AUCR 2.56 (Table 2)  Cmax 1.74 -> predicted 1.57
    #   rifampicin    AUCR 0.35 (Table 2)  Cmax 0.54 -> predicted 0.53
    #   fluconazole   AUCR 1.98 (Table 3)  Cmax 1.45 -> predicted 1.41
    #   erythromycin  AUCR 2.31 (Table 3)  Cmax 1.58 -> predicted 1.50
    #   efavirenz 400 AUCR 0.58 (Table 3)  Cmax 0.72 -> predicted 0.73
    #   efavirenz 600 AUCR 0.33 (Table 3)  Cmax 0.50 -> predicted 0.51
    # The activities ladder monotonically with the perpetrators' FDA
    # potency classification - ketoconazole is the strongest inhibitor
    # and rifampicin the strongest inducer - which is itself a
    # consistency check on the back-solving.
    e_conmed_ketoconazole_cyp3a4 <- fixed(-2.2315747467)
    label("log relative CYP3A4 activity with ketoconazole 400 mg once daily")

    e_conmed_rifampicin_cyp3a4 <- fixed(1.0043012336)
    label("log relative CYP3A4 activity with rifampicin 600 mg once daily")

    e_conmed_fluconazole_cyp3a4 <- fixed(-1.1875914226)
    label("log relative CYP3A4 activity with fluconazole 200 mg once daily")

    e_conmed_erythromycin_cyp3a4 <- fixed(-1.6997279711)
    label("log relative CYP3A4 activity with erythromycin 500 mg three times daily")

    e_conmed_efv_cyp3a4 <- fixed(0.5751309406)
    label("log relative CYP3A4 activity with efavirenz 400 mg once daily")

    e_dose_efavirenz_mg_cyp3a4 <- fixed(0.4750705574)
    label("Additional log relative CYP3A4 activity when the efavirenz dose is 600 mg once daily")
    # Increment on top of e_conmed_efv_cyp3a4, so the 600 mg activity is
    # exp(0.5751309 + 0.4750706) = 2.8582, back-solved from the
    # published AUC ratio of 0.33.

    propSd <- fixed(0)
    label("Proportional residual error (none reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Typical-value parameters. No random effects.
    # ------------------------------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc)
    cl_renal <- exp(lcl_renal)
    cl_nonren <- exp(lcl_nonren)

    # ------------------------------------------------------------------
    # 2. Relative CYP3A4 activity contributed by a coadministered
    # perpetrator. All flags zero (GXR monotherapy) gives a3a4 = 1. The
    # terms are additive on the log scale, so the flags are
    # multiplicative on activity; the source only ever simulates one
    # perpetrator at a time.
    #
    # efv_hi selects the 600 mg efavirenz coefficient. It is written
    # with a single comparison so that any nonzero efavirenz dose below
    # 500 mg falls through to the 400 mg coefficient rather than
    # silently contributing nothing.
    # ------------------------------------------------------------------
    efv_hi <- CONMED_EFV * (DOSE_EFAVIRENZ_MG >= 500)

    a3a4 <- exp(
      e_conmed_ketoconazole_cyp3a4 * CONMED_KETOCONAZOLE +
        e_conmed_rifampicin_cyp3a4 * CONMED_RIFAMPICIN +
        e_conmed_fluconazole_cyp3a4 * CONMED_FLUCONAZOLE +
        e_conmed_erythromycin_cyp3a4 * CONMED_ERYTHROMYCIN +
        e_conmed_efv_cyp3a4 * CONMED_EFV +
        e_dose_efavirenz_mg_cyp3a4 * efv_hi
    )

    # Well-stirred scaling of each extraction ratio. Multiplying an
    # intrinsic clearance by a factor A takes an extraction ratio E to
    # E*A / (1 - E + E*A), so the corresponding availability 1 - E is
    # divided by (1 - E + E*A). These two denominators are all that is
    # needed. Both the gut wall and the non-renal hepatic pathway are
    # treated as entirely CYP3A4, which is the paper's own assumption
    # (Section 2.2.1 assigns the whole gut intrinsic clearance to CYP3A
    # and Section 4 assigns the whole metabolic fraction to CYP3A4).
    gfac <- 1 - egut + egut * a3a4
    hfac <- 1 - eh + eh * a3a4

    # ------------------------------------------------------------------
    # 3. Clearance. Renal clearance is untouched by a CYP3A4
    # perpetrator; the non-renal arm follows the hepatic extraction
    # ratio. With a3a4 = 1 this is 12.6 + 12.2 = 24.8 L/h, the observed
    # intravenous clearance.
    # ------------------------------------------------------------------
    cl <- cl_renal + cl_nonren * a3a4 / hfac
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 4. ODE system. One distribution compartment, matching the single
    # Vss the source reports.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ------------------------------------------------------------------
    # 5. Bioavailability, F = fa * F_G * F_H, with each first-pass
    # availability raised by a CYP3A4 inhibitor and lowered by an
    # inducer. With a3a4 = 1 this is 1 * 0.723745 * 0.906513 = 0.656085,
    # so CL/F = 24.8 / 0.656085 = 37.8 L/h, the observed value the
    # source calibrated Q_G against.
    # ------------------------------------------------------------------
    fdepot <- fa * ((1 - egut) / gfac) * ((1 - eh) / hfac)
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 6. Observation. Doses are in mg and vc is in L, so central / vc is
    # in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    # throughout Li 2018 Tables 2 and 3.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

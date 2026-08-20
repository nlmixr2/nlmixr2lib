Yang_2024_meropenem_pbpk <- function() {
  description <- paste(
    "PBPK-derived reduced one-compartment model for meropenem in critically ill",
    "(severe pneumonia with or without sepsis) adult Asian ICU patients. The source",
    "paper built an 18-tissue whole-body PBPK model in PK-Sim v11.2 for healthy",
    "adults and extrapolated it to critically ill patients by scaling albumin,",
    "alpha-1-acid glycoprotein, hematocrit and GFR; that whole-body structure is a",
    "PK-Sim platform model whose ODEs, organ volumes, blood flows and partition",
    "coefficients are NOT written out in the publication and are therefore not",
    "reproduced here. What IS fully reported, and is what this file encodes, is the",
    "reduced disposition model the authors carried out of the PBPK and into their",
    "Monte Carlo target-attainment analysis: Vd 23.21 L, CL 12.07 L/h and unbound",
    "fraction 0.98 for Asian patients with severe infection (Results, Monte Carlo",
    "simulations). Intravenous infusion dosing; linear elimination; no covariates.",
    "The PK/PD index is the fraction of the dosing interval during which the free",
    "concentration exceeds the MIC (f%T>MIC), with 40%fT>MIC as the primary target",
    "and 100%fT>MIC as a stricter alternative; f%T>MIC is computed from the",
    "simulated profile rather than integrated as a model state, following the",
    "Minichmayr_2024_ceftaroline precedent. Interindividual variability was assumed",
    "lognormal on Vd and CL in the Monte Carlo but its magnitude is not reported",
    "anywhere in the paper or its supplement, so both etas are encoded as fixed(0);",
    "no residual error model was reported.",
    sep = " "
  )
  reference <- paste(
    "Yang Y, Wang Y, Zeng W, Zhou J, Xu M, Lan Y, Liu L, Shen J, Zhang C, He Q.",
    "Physiologically-based pharmacokinetic/pharmacodynamic modeling of meropenem",
    "in critically ill patients. Sci Rep. 2024;14:19249.",
    "doi:10.1038/s41598-024-64223-0. PMCID PMC11335869.",
    "The reduced disposition parameters encoded here (Vd 23.21 L, CL 12.07 L/h,",
    "f 0.98 in Asians with severe infection) are stated in the Results section",
    "'Monte carlo simulations'; Vd 23.21 L also appears in Table 1 as the",
    "parameter-identification result. The PK/PD target 40%fT>MIC, the MIC panel",
    "(1, 2, 4, 8, 16 mg/L) and the dosing regimens (1 g q12h and 1 g q8h at 0.5,",
    "4 and 6 h infusion) are Methods 'Monte carlo simulations'; the resulting PTA",
    "curves are Figure 8B (Asians) and Figure S1 (100%fT>MIC, Asians). Supplement",
    "Table S2 gives the critically-ill physiological scaling factors applied inside",
    "PK-Sim; Table S1 lists the 40 literature PK studies used to build and verify",
    "the PBPK model; Table S5 lists the 31 prospective ICU patients.",
    sep = " "
  )
  vignette <- "Yang_2024_meropenem_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 31,
    n_studies      = 40,
    age_range      = "38-95 years (prospective ICU cohort, Table S5)",
    age_median     = "66 years (mean; Results 'Prospective analysis of meropenem')",
    weight_range   = "42-85 kg (prospective ICU cohort, Table S5)",
    weight_median  = "57 kg (mean; Results 'Prospective analysis of meropenem')",
    sex_female_pct = 45.2,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Critically ill adults with severe pneumonia with or without sepsis,",
      "requiring ICU treatment; patients with cirrhosis or any liver damage,",
      "and patients in renal failure, were excluded (Methods 'Clinical data",
      "collection'). Six of the 31 prospective patients (19.35%) were on CRRT.",
      sep = " "
    ),
    renal_function = paste(
      "Mean serum creatinine 87.06 umol/L in the prospective cohort (range",
      "25.8-184 umol/L, Table S5). Inside PK-Sim the critically-ill",
      "extrapolation applied a GFR scaling factor of 0.50 relative to the",
      "PK-Sim healthy control value (Table S2).",
      sep = " "
    ),
    dose_range = paste(
      "Prospective cohort: 0.5-1 g every 8 or 12 h as 0.5, 4 or 6 h intravenous",
      "infusions (Table S5). Monte Carlo target-attainment analysis: 1 g q12h and",
      "1 g q8h, each at 0.5, 4 and 6 h infusion (Methods 'Monte carlo",
      "simulations'). The literature studies used to build and verify the PBPK",
      "model span 0.25-3 g by bolus, 0.5-3 h infusion and 24 h continuous",
      "infusion (Table S1).",
      sep = " "
    ),
    regions = "China (ICU of the Third People's Hospital of Chengdu); literature cohorts from Europe, North America and Asia",
    notes = paste(
      "Three distinct populations appear in the paper and should not be",
      "conflated. (1) The PBPK model was developed for a typical 30-year-old",
      "European male using PK-Sim's virtual population, then verified against 12",
      "healthy-adult literature PK studies. (2) It was extrapolated to critically",
      "ill patients and verified against 10 literature studies in severe",
      "infection; the extrapolation scaled four physiological parameters relative",
      "to PK-Sim healthy control values (Table S2, with the sensitivity index the",
      "paper computed for each): albumin 0.54 (sensitivity 0.78), alpha-1-acid",
      "glycoprotein 1.62 (0.73), hematocrit 0.77 (0.54) and GFR 0.50 (0.87). The",
      "absolute control values those factors multiply are PK-Sim internals and are",
      "not printed in the paper, which is one reason the whole-body layer is not",
      "reproducible here. (3) A prospective TDM cohort of 92 plasma samples from",
      "31 Chengdu ICU patients was used as an external check; 71 of 92 samples",
      "(77.17%) fell inside the predicted 5th-95th percentile band (Figure 7).",
      "The n_subjects field records the prospective cohort (31); n_studies records",
      "the 40 literature PK study arms in Table S1.",
      "The Monte Carlo target-attainment analysis simulated 10,000 virtual",
      "patients with Pseudomonas aeruginosa infection.",
      "Sex split: 14 of 31 prospective patients were female (Table S5).",
      sep = " "
    )
  )

  ini({
    # Reduced disposition parameters carried out of the PK-Sim whole-body PBPK
    # model and into the Monte Carlo target-attainment analysis. All three are
    # wrapped in fixed() because none was estimated by fitting this reduced
    # model: they are PBPK outputs held constant as Monte Carlo inputs
    # ("The Vd, CL, and f were defined as 23.21 L, 12.07 L/h, and 0.98 in
    # Asians with severe infection from PBPK model").
    #
    # Internal consistency check on CL: 12.07 L/h is 2.87 mL/min/kg at 70 kg,
    # which is exactly the predicted critically-ill CL that Table 2 reports for
    # the 1 g intravenous regimens (Cohen 2002 and Frippiat 2015 rows).
    lcl <- fixed(log(12.07)); label("Clearance (L/h)")                    # Yang 2024 Results 'Monte carlo simulations' (CL = 12.07 L/h, Asians with severe infection)
    lvc <- fixed(log(23.21)); label("Volume of distribution (L)")         # Yang 2024 Results 'Monte carlo simulations' and Table 1 (Vd = 23.21 L, parameter identification)

    # Unbound fraction in plasma. Table 1 lists fu_p = 0.98 as the PBPK input
    # (from Martins 2020); the Monte Carlo re-used the same value, sampling it
    # from a uniform distribution whose bounds the paper does not state.
    fu <- fixed(0.98); label("Fraction unbound in plasma (fraction)")     # Yang 2024 Table 1 (fu p = 0.98) and Results 'Monte carlo simulations' (f = 0.98)

    # Interindividual variability. Methods 'Monte carlo simulations' states
    # that Vd and CL were assumed to follow a log-gaussian distribution, but
    # the magnitude of that variability is reported NOWHERE in the paper or its
    # supplement. These are encoded as fixed(0) to record that the authors
    # declared lognormal IIV on exactly these two parameters while making clear
    # that no variance was published; fixed(0) here means "magnitude
    # unreported", NOT "the authors estimated zero variability". Do not treat
    # simulations from this model as reproducing the published PTA spread
    # without supplying an omega. See the vignette Errata.
    etalcl ~ fixed(0)   # Yang 2024 Methods 'Monte carlo simulations' (lognormal IIV on CL declared; magnitude not reported)
    etalvc ~ fixed(0)   # Yang 2024 Methods 'Monte carlo simulations' (lognormal IIV on Vd declared; magnitude not reported)

    # No residual error model was reported: the PBPK model is a forward
    # simulation evaluated by MFE / GMFE against literature means, not a
    # likelihood-based popPK fit, so no sigma exists to transcribe.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported by the source)")   # Yang 2024 - no residual error model reported
  })

  model({
    # Individual disposition parameters. The IIV terms are declared by the
    # paper but carry zero variance here (see ini()).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment disposition with linear elimination. Meropenem is given
    # as an intravenous infusion, so doses go directly to central; set rate =
    # or dur = on the dose records to reproduce the 0.5, 4 and 6 h infusions.
    d/dt(central) <- -kel * central

    # Total plasma meropenem concentration (mg/L).
    Cc <- central / vc

    # Free (unbound) plasma concentration (mg/L). This is the quantity the
    # PK/PD index is defined on: the target is the fraction of the dosing
    # interval with Cfree > MIC (40%fT>MIC primary, 100%fT>MIC alternative).
    # f%T>MIC is computed from the simulated profile in the vignette rather
    # than integrated as a model state, matching how the authors computed it
    # (Oracle Crystal Ball Monte Carlo over sampled PK parameters) and
    # following Minichmayr_2024_ceftaroline.
    Cfree <- fu * Cc

    Cc ~ prop(propSd)
  })
}

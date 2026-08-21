Izat_2025_ziprasidone <- function() {
  description <- paste(
    "One-compartment intravenous pharmacokinetic reduction of the Simcyp",
    "full-PBPK model for the antipsychotic ziprasidone in healthy adults",
    "(Izat 2025). The source model was built in the Simcyp simulator",
    "(version 21) and its whole-body mass-balance equations are not",
    "published, so the platform model itself cannot be encoded here. What",
    "IS fully reported is the ziprasidone compound layer: Table 1 gives",
    "the steady-state volume of distribution recovered by the fitted",
    "global tissue-partition scalar, and the observed intravenous plasma",
    "clearance. No parameter is fitted here and none is imported from a",
    "Simcyp population file: every value is either a Table 1 entry or an",
    "arithmetic consequence of one at the assumed 70 kg reference body",
    "weight. Only the intravenous arm is reproduced; the oral arm needs a",
    "hepatic blood flow that the paper never prints. Because the source",
    "used a full-PBPK distribution model, only a single lumped",
    "steady-state volume is reported and this reduction is",
    "mono-exponential, which underpredicts the observed peak",
    "concentration by about 15 percent; see the vignette for the",
    "quantitative comparison. This is a typical-value simulation model:",
    "the source reports no inter-individual variance components and no",
    "residual-error model, so there are no etas and propSd is fixed at",
    "zero. The fraction metabolized by aldehyde oxidase and by CYP3A4",
    "that is the paper's main contribution is a property of the",
    "unpublished platform layer and is NOT reproducible from this model;",
    "see the vignette for the full list of deviations.",
    sep = " "
  )
  reference <- paste(
    "Izat N, Bolleddula J, Carione P, Huertas Valentin L, Jones RS,",
    "Kulkarni P, Moss D, Peterkin VC, Tian DD, Treyer A, Venkatakrishnan K,",
    "Zientek MA, Barber J, Houston JB, Galetin A, Scotcher D. (2025).",
    "Establishing a physiologically based pharmacokinetic framework for",
    "aldehyde oxidase and dual aldehyde oxidase-CYP substrates.",
    "CPT Pharmacometrics Syst Pharmacol 14(1):164-178.",
    "doi:10.1002/psp4.13255.",
    sep = " "
  )
  vignette <- "Izat_2025_aldehyde_oxidase_substrates"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what the single ODE state holds, in what amount units, in
  # what biological matrix. Verified against Izat 2025 Table 1 (ziprasidone
  # column, "Distribution model: Full PBPK") and Appendix S1 Supplement 2,
  # which states that the global Kp scalar was fitted "to recover the
  # observed volume of distribution at steady state (Vss)" against observed
  # plasma concentrations.
  compartmentData <- list(
    central = list(
      analyte = "ziprasidone", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  # Body weight is implicit in the volume input but is never printed.
  # Recorded as screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Izat 2025 Table 1 expresses Vss in L/kg, so the",
        "Simcyp model does scale distribution volume with weight. This",
        "reduction fixes a 70 kg reference weight instead of carrying a",
        "weight term, because the corresponding weight scaling of",
        "clearance lives in the unpublished Simcyp Healthy Volunteers",
        "population file. Scaling volume with weight while holding",
        "clearance fixed would produce an internally inconsistent model,",
        "so neither is scaled."
      ),
      units = "kg",
      type  = "continuous",
      notes = paste(
        "Implicit in the L/kg volume input; not carried. No body weight",
        "is printed anywhere in the paper or Appendix S1 -- Table S3",
        "reports only age range and female fraction. The 70 kg value is",
        "the standing rounded-standard assumption; see the vignette",
        "Errata."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "19-37 years",
    weight_median  = "70 kg (assumed reference weight for the L/kg volume input; not reported)",
    sex_female_pct = 0,
    disease_state  = "Healthy adult volunteers.",
    dose_range     = "Single 5 mg intravenous infusion.",
    route          = "intravenous",
    studies        = paste(
      "Miceli JJ, Wilner KD, Swan SK, Tensfeldt TG. Pharmacokinetics,",
      "safety, and tolerability of intramuscular ziprasidone in healthy",
      "volunteers. J Clin Pharmacol 2005;45(6):620-630, cited as",
      "reference 23 of Appendix S1 Table S3. Twelve healthy men aged",
      "19-37 years received 5 mg ziprasidone intravenously (a separate",
      "cohort of thirteen received 20 mg orally). The intravenous arm is",
      "the model-development dataset for the distribution parameter",
      "reproduced here."
    ),
    notes          = paste(
      "This is a PBPK analysis rather than a population-PK fit, so there",
      "is no pooled analysis dataset and no estimated variance",
      "components. The percent-CV figures in Table 1 and Appendix S1",
      "Table S10 are the spread of a Simcyp virtual population driven by",
      "unpublished population files, not estimated omegas.",
      "Observed fractions metabolized from the human mass-balance study",
      "(Appendix S1 Table S2) are fmAO 0.67 and fmCYP3A4 0.33, with renal",
      "clearance assumed negligible; those fractions are recorded here",
      "for provenance but are not encoded as parameters because they do",
      "not affect the plasma concentration-time profile and their pathway",
      "split is a property of the unpublished platform layer. Table 1",
      "footnote f records that no reaction-phenotyping data with an",
      "aldehyde oxidase inhibitor were available for ziprasidone, so its",
      "hepatocyte AO intrinsic clearance is structurally unavailable",
      "rather than merely unreported."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are either verbatim Izat 2025 Table 1
    # entries or arithmetic consequences of them. The one assumption is
    # the 70 kg reference body weight needed to turn the L/kg volume
    # input into litres; no body weight is printed anywhere in the paper
    # (see covariatesDataExcluded$WT and the vignette Errata).
    #
    # WHY THIS IS A ONE-COMPARTMENT MODEL. Izat 2025 Table 1 classifies
    # ziprasidone under "Distribution model: Full PBPK" with a fitted
    # global Kp scalar of 0.56. In a full-PBPK layout the multi-tissue
    # structure IS the distribution model, and the only distribution
    # quantity Table 1 reports is the lumped steady-state volume that the
    # Kp scalar was tuned to recover; no individual tissue volumes,
    # blood flows or partition coefficients are printed, and Table 1
    # gives neither a Vsac nor a Q for this compound. The reduction is
    # therefore mono-exponential. Unlike zaleplon -- whose reported
    # inter-compartmental clearance of zero makes its reduction exact --
    # this IS an approximation of a genuinely multi-exponential profile,
    # and it costs about 15 percent on the peak concentration. That
    # deviation is reported in the vignette rather than tuned away.
    # ------------------------------------------------------------------

    # Central volume. Table 1 gives the steady-state volume of
    # distribution Vss = 1.01 L/kg. At the 70 kg reference weight:
    #   1.01 L/kg * 70 kg = 70.7 L
    # There is no Vsac to subtract for this compound.
    lvc <- fixed(log(70.7))
    label("Central compartment volume vc (L) at the 70 kg reference weight")
    # Derived from Izat 2025 Table 1 (Vss 1.01 L/kg, recovered by the
    # optimized global Kp scalar of 0.56).

    # Systemic plasma clearance. Table 1: CLiv = 23.04 L/h [14% CV]
    # (i.v.). Renal clearance is reported as negligible, so this is
    # entirely non-renal and is carried as a single lcl rather than split
    # into lcl_renal / lcl_nonren.
    lcl <- fixed(log(23.04))
    label("Systemic plasma clearance CL (L/h)")
    # Izat 2025 Table 1, ziprasidone column, "CLiv/oral (L/h) [CV%]":
    # 23.04 [14] (i.v.). Cross-checked against Appendix S1 Table S10,
    # which reports observed CLiv 23.04 (14) and observed AUCinf 217
    # ng.h/mL after 5 mg, consistent with 5 mg / 23.04 L/h = 217 ng.h/mL.

    # Izat 2025 is a PBPK simulation analysis, not a population-PK fit.
    # It reports no residual-error model and no inter-individual variance
    # components. Rather than invent a variance, the residual error is
    # fixed at zero, which makes this a deterministic typical-value
    # simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    # ------------------------------------------------------------------
    vc <- exp(lvc)
    cl <- exp(lcl)

    # ------------------------------------------------------------------
    # 2. Micro-constant for elimination.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system. A single compartment, as justified in the ini()
    #    header. Amounts are in mg. Intravenous administration only --
    #    there is no depot, because the oral bioavailability of this
    #    compound cannot be reconstructed from the paper (the hepatic
    #    first-pass term needs a hepatic blood flow the paper never
    #    prints).
    # ------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ------------------------------------------------------------------
    # 4. Observation. Doses are in mg and vc is in L, so central / vc is
    #    in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units used
    #    throughout Izat 2025 Appendix S1 Table S10.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

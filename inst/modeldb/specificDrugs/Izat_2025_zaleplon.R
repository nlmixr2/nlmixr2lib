Izat_2025_zaleplon <- function() {
  description <- paste(
    "One-compartment intravenous pharmacokinetic reduction of the Simcyp",
    "minimal-PBPK model for the sedative-hypnotic zaleplon in healthy",
    "adults (Izat 2025). The source model was built in the Simcyp",
    "simulator (version 21) and its whole-body mass-balance equations are",
    "not published, so the platform model itself cannot be encoded here.",
    "What IS fully reported is the zaleplon compound layer: Table 1 gives",
    "the steady-state volume of distribution, the volume of the single",
    "adjusting compartment (SAC), the inter-compartmental clearance, and",
    "the observed intravenous plasma clearance separately in women and in",
    "men. Because the reported inter-compartmental clearance is exactly",
    "zero and the SAC volume is 1e-05 L/kg, the SAC is kinetically",
    "disconnected and the source model reduces without approximation to a",
    "single compartment. Clearance is sex-dependent, which is the only",
    "covariate effect anywhere in the paper. No parameter is fitted here",
    "and none is imported from a Simcyp population file: every value is",
    "either a Table 1 entry or an arithmetic consequence of one at the",
    "assumed 70 kg reference body weight. Only the intravenous arm is",
    "reproduced; the oral arm needs a hepatic blood flow that the paper",
    "never prints. This is a typical-value simulation model: the source",
    "reports no inter-individual variance components and no residual-error",
    "model, so there are no etas and propSd is fixed at zero. The fraction",
    "metabolized by aldehyde oxidase and by CYP3A4 that is the paper's main",
    "contribution is a property of the unpublished platform layer and is",
    "NOT reproducible from this model; see the vignette for the full list",
    "of deviations.",
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
  # what biological matrix. Verified against Izat 2025 Table 1 (zaleplon
  # column, "Distribution model: Minimal PBPK") and Appendix S1 Supplement 2,
  # which describes the central compartment of the minimal-PBPK layout as a
  # non-physiological systemic compartment whose concentration is compared
  # against observed plasma data.
  compartmentData <- list(
    central = list(
      analyte = "zaleplon", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    SEXF = list(
      description = paste(
        "Biological sex, 1 = female and 0 = male. Izat 2025 Table 1",
        "reports two separate observed intravenous plasma clearances for",
        "zaleplon, 52.51 L/h in women and 71.58 L/h in men, and the",
        "Simcyp simulations were run separately for a female and a male",
        "cohort (Appendix S1 Table S10 tabulates every zaleplon result",
        "under an F and an M row). This is the only covariate effect",
        "anywhere in the paper. Men are the reference category, so the",
        "effect parameter e_sexf_cl carries the female-to-male clearance",
        "ratio."
      ),
      units = "(binary)",
      type  = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Not an estimated covariate coefficient. The two clearances are",
        "the observed values from the source pharmacokinetic study",
        "(Rosen 1999, Appendix S1 Table S3 reference 22), entered into",
        "Simcyp as sex-specific inputs. The steady-state volume of",
        "distribution is NOT sex-split in Table 1, consistent with the",
        "near-identical observed peak concentrations in the two cohorts",
        "(62.7 versus 61.4 ng/mL, Appendix S1 Table S10)."
      ),
      source_name = "female / male cohort labels F and M"
    )
  )

  # Body weight is implicit in every volume input but is never printed.
  # Recorded as screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Izat 2025 Table 1 expresses Vss and Vsac in L/kg,",
        "so the Simcyp model does scale distribution volume with weight.",
        "This reduction fixes a 70 kg reference weight instead of",
        "carrying a weight term, because the corresponding weight",
        "scaling of clearance lives in the unpublished Simcyp Healthy",
        "Volunteers population file. Scaling volume with weight while",
        "holding clearance fixed would produce an internally",
        "inconsistent model, so neither is scaled."
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
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "19-32 years",
    weight_median  = "70 kg (assumed reference weight for the L/kg volume input; not reported)",
    sex_female_pct = 50,
    disease_state  = "Healthy adult volunteers.",
    dose_range     = "Single 5 mg intravenous infusion.",
    route          = "intravenous",
    studies        = paste(
      "Rosen AS, Fournie P, Darwish M, Danjou P, Troy SM. Zaleplon",
      "pharmacokinetics and absolute bioavailability. Biopharm Drug",
      "Dispos 1999;20(3):171-175, cited as reference 22 of Appendix S1",
      "Table S3. Ten healthy volunteers aged 19-32 years received 5 mg",
      "zaleplon intravenously and orally; Appendix S1 Table S3 records",
      "the female fraction as 0 or 1, i.e. separate all-female and",
      "all-male cohorts, which is why Table 1 carries two clearances.",
      "The intravenous arm is the model-development dataset for the",
      "distribution parameters reproduced here."
    ),
    notes          = paste(
      "This is a PBPK analysis rather than a population-PK fit, so there",
      "is no pooled analysis dataset and no estimated variance",
      "components. The percent-CV figures in Table 1 and Appendix S1",
      "Table S10 are the spread of a Simcyp virtual population driven by",
      "unpublished population files, not estimated omegas.",
      "Observed fractions metabolized from the human mass-balance study",
      "(Appendix S1 Table S2) are fmAO 0.57-0.74 and fmCYP3A4 0.26-0.43,",
      "with renal clearance assumed negligible; those fractions are",
      "recorded here for provenance but are not encoded as parameters",
      "because they do not affect the plasma concentration-time profile",
      "and their pathway split is a property of the unpublished platform",
      "layer."
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
    # zaleplon under "Distribution model: Minimal PBPK", which Appendix
    # S1 Supplement 2 describes as "a non-physiological central
    # compartment excluding the liver and portal vein coupled with a
    # non-physiological non-central compartment ... (SAC)". For zaleplon
    # the reported SAC parameters are Vsac = 0.00001 L/kg and Q = 0 L/h.
    # An inter-compartmental clearance of exactly zero means no drug ever
    # enters the SAC, so the two-compartment layout collapses to a single
    # compartment with no approximation at all. The reduction below is
    # therefore exact with respect to the reported compound layer, not a
    # simplification of it.
    # ------------------------------------------------------------------

    # Systemic volume. Table 1 gives the whole-body steady-state volume
    # Vss = 1.16 L/kg and the SAC volume Vsac = 0.00001 L/kg. Because the
    # SAC never fills (Q = 0), the volume the drug actually distributes
    # into is the remainder, at the 70 kg reference weight:
    #   (1.16 - 0.00001) L/kg * 70 kg = 81.1993 L
    # The liver and portal-vein volumes that the minimal-PBPK layout
    # excludes from the central compartment are not reported and are not
    # subtracted; because they re-enter the observable steady-state
    # volume they do not affect a one-compartment reduction. See the
    # vignette Errata.
    lvc <- fixed(log(81.1993))
    label("Central compartment volume vc (L) at the 70 kg reference weight")
    # Derived from Izat 2025 Table 1 (Vss 1.16 L/kg; Vsac 0.00001 L/kg).

    # Systemic plasma clearance in MEN, the reference category of SEXF.
    # Table 1: CLiv = 71.58 L/h [20.2% CV] (i.v., male). Renal clearance
    # is reported as negligible, so this is entirely non-renal and is
    # carried as a single lcl rather than split into lcl_renal /
    # lcl_nonren.
    lcl <- fixed(log(71.58))
    label("Systemic plasma clearance CL in men (L/h)")
    # Izat 2025 Table 1, zaleplon column, "CLiv/oral (L/h) [CV%]":
    # 71.58 [20.2] (i.v., male). Cross-checked against Appendix S1
    # Table S10, which reports observed CLiv 71.6 (20.2) for the male arm.

    # Female-to-male clearance ratio, applied multiplicatively:
    #   52.51 / 71.58 = 0.733585  ->  log(0.733585) = -0.309812
    # Table 1 gives CLiv = 52.51 L/h [22.5% CV] (i.v., female).
    e_sexf_cl <- fixed(log(0.733585))
    label("Multiplicative effect of female sex on CL (log ratio; women/men)")
    # Derived from Izat 2025 Table 1, zaleplon column: 52.51 L/h female
    # divided by 71.58 L/h male. Cross-checked against Appendix S1
    # Table S10 observed CLiv 52.5 (22.5) female and 71.6 (20.2) male.

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
    # 1. Individual parameters. One binary covariate, no random effects.
    #    SEXF = 1 for women, 0 for men; men are the reference category.
    # ------------------------------------------------------------------
    vc <- exp(lvc)
    cl <- exp(lcl + e_sexf_cl * SEXF)

    # ------------------------------------------------------------------
    # 2. Micro-constant for elimination.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system. A single compartment, as justified in the ini()
    #    header: the reported SAC has Q = 0 L/h and so is unreachable.
    #    Amounts are in mg. Intravenous administration only -- there is
    #    no depot, because the oral bioavailability of this compound
    #    cannot be reconstructed from the paper (the hepatic first-pass
    #    term needs a hepatic blood flow the paper never prints).
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

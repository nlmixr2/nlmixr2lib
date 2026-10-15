Wu_2022_treprostinil <- function() {
  description <- paste(
    "One-compartment intravenous pharmacokinetic reduction of the Simcyp",
    "full-PBPK model for the prostacyclin analogue treprostinil in healthy",
    "adults (Wu 2022). The source model was built in the Simcyp simulator",
    "(version 17 release 1) and its whole-body mass-balance equations are",
    "not published, so the platform model itself cannot be encoded here.",
    "What IS reported is the treprostinil intravenous disposition layer:",
    "the total in vivo intravenous clearance of 43 L/h that the retrograde",
    "clearance model was built to reproduce, its printed renal component",
    "of 0.9 L/h (Table 4), and the final model's predicted steady-state",
    "volume of distribution of 0.42 L/kg (Table 3). The elimination is",
    "carried as an explicit renal plus non-renal sum. No parameter is",
    "fitted here and none is imported from a Simcyp population file: every",
    "value is either a printed entry or an arithmetic consequence of one",
    "at the assumed 70 kg reference body weight. Because the source used a",
    "full-PBPK distribution model, only a single lumped steady-state volume",
    "is reported and this reduction is mono-exponential. It reproduces the",
    "paper's own predicted Cmax and AUC for all three intravenous infusion",
    "studies (Supplemental Table 1) to within 14 percent, with no fitted",
    "parameter. Only the intravenous arm is reproduced: the extended-release",
    "oral tablet arm, and with it the hepatic-impairment and patient",
    "extrapolations, rests on an absorption model and an oral-model",
    "clearance the paper does not print; see the vignette for the",
    "quantitative reason. This is a typical-value simulation model: the",
    "source reports no inter-individual variance components and no",
    "residual-error model, so there are no etas and propSd is fixed at",
    "zero.",
    sep = " "
  )
  reference <- paste(
    "Wu X, Zhang X, Xu R, Shaik IH, Venkataramanan R. (2022).",
    "Physiologically based pharmacokinetic modelling of treprostinil after",
    "intravenous injection and extended-release oral tablet administration",
    "in healthy volunteers: An extrapolation to other patient populations",
    "including patients with hepatic impairment.",
    "Br J Clin Pharmacol 88(2):587-599.",
    "doi:10.1111/bcp.14966.",
    sep = " "
  )
  vignette <- "Wu_2022_treprostinil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what the single ODE state holds, in what amount units, in
  # what biological matrix. Verified against Wu 2022 section 2.2 ("Full
  # PBPK model with Rodgers and Rowland method was used to predict the
  # tissue:plasma partition coefficient") and Table 3, whose only lumped
  # distribution quantity is the predicted Vss.
  compartmentData <- list(
    central = list(
      analyte = "treprostinil",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  # Body weight is implicit in the volume input but is never printed.
  # Recorded as screened-but-not-carried rather than silently dropped.
  covariatesDataExcluded <- list(
    WT = list(
      description = paste(
        "Body weight. Wu 2022 Table 3 expresses Vss in L/kg, and two of",
        "the three intravenous studies dosed per kilogram (Table 1), so",
        "the Simcyp model does scale distribution volume with weight. This",
        "reduction fixes a 70 kg reference weight instead of carrying a",
        "weight term, because the corresponding weight scaling of",
        "clearance lives in the unpublished Simcyp Healthy Volunteers",
        "population file; the printed clearance of 43 L/h is an absolute",
        "value. Scaling volume with weight while holding clearance fixed",
        "would produce an internally inconsistent model, so neither is",
        "scaled."
      ),
      units = "kg",
      type = "continuous",
      notes = paste(
        "Implicit in the L/kg volume input; not carried. No body weight is",
        "printed anywhere in the paper or its supplement -- Table 1",
        "reports only race, sex split and age range. The 70 kg value is",
        "the standing rounded-standard assumption; see the vignette",
        "Assumptions and deviations."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 90L,
    n_studies = 3L,
    age_range = "18-63 years across the three intravenous studies (Table 1; the third study's range is unreported and was defaulted to 18-65 years in the simulations)",
    weight_median = "70 kg (assumed reference weight for the L/kg volume input; not reported)",
    sex_female_pct = "40-50 across the three intravenous studies (Table 1)",
    race_ethnicity = "Predominantly white (Table 1: 8 white, 3 black, 4 Hispanic; 26 white, 8 black, 17 other; third study unreported); simulated as the Simcyp healthy Caucasian population",
    disease_state = "Healthy adult volunteers.",
    dose_range = paste(
      "Single intravenous infusions of 0.00225 mg/kg over 2.5 h,",
      "0.0432 mg/kg over 72 h and 0.2 mg over 4 h (Table 1)."
    ),
    route = "intravenous",
    regions = "Simcyp virtual healthy (Caucasian) population (section 2.2).",
    studies = paste(
      "Three published intravenous studies in healthy volunteers",
      "(Table 1). Wade et al. (reference 10; n = 15, 0.00225 mg/kg over",
      "2.5 h) was used to build the model; Laliberte et al. (reference",
      "11; n = 51, 0.0432 mg/kg over 72 h) and the treprostinil",
      "extended-release tablet package insert (reference 4; n = 24, 0.2",
      "mg over 4 h) were the model-naive verification studies. The 43",
      "L/h clearance is the average intravenous clearance of references",
      "4, 10 and 11."
    ),
    notes = paste(
      "This is a PBPK analysis rather than a population-PK fit, so there",
      "is no pooled analysis dataset and no estimated variance",
      "components; each simulation was 10 virtual trials of 10 subjects",
      "matched to the cohort demographics, and the 5th-95th percentile",
      "bands in Figures 1 and 3-5 are the spread of that virtual",
      "population. Elimination as built in the source (section 2.2,",
      "Table 4): renal clearance 0.9 L/h (2% of 43 L/h), additional",
      "systemic clearance 2.6 L/h (6%), biliary clearance 1% included",
      "in the hepatic component, and CYP2C8 / CYP2C9 contributing 90% /",
      "10% of metabolism through retrograde intrinsic clearances of 20.5",
      "and 0.75 uL/min/pmol. Those fractions are recorded for provenance;",
      "the renal component is encoded as lcl_renal, the rest is the",
      "lumped lcl_nonren. Reported but NOT reproducible from this model:",
      "the extended-release oral tablet absorption (ADAM model with a",
      "digitised dissolution profile, jejunal Peff 0.46e-4 cm/s, colon",
      "absorption scalar 0.05, plasma fraction unbound adjusted to 0.04 in",
      "the oral model; predicted fa 0.48, Fg 0.87, Fh 0.60 and ka 0.20",
      "1/h), the Child-Pugh A/B/C hepatic-impairment extrapolation",
      "(predicted AUC 2.4-, 3.8- and 4.6-fold healthy, from Simcyp",
      "cirrhosis population files), and the predicted lung:plasma",
      "partition coefficient of 0.17."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter below is fixed: nothing was estimated in building
    # this reduction. Values are either verbatim Wu 2022 entries or
    # arithmetic consequences of them. The one assumption is the 70 kg
    # reference body weight needed to turn the L/kg volume into litres;
    # no body weight is printed anywhere in the paper or its supplement
    # (see covariatesDataExcluded$WT and the vignette).
    #
    # WHY THIS IS A ONE-COMPARTMENT MODEL. Section 2.2 states that a
    # full PBPK model with Rodgers and Rowland partition coefficients was
    # used, with the Kp scalar and the adipose Kp optimized. Table 3
    # prints the twelve tissue Kp values but no organ volumes or blood
    # flows, so the multi-tissue structure cannot be rebuilt; the only
    # lumped distribution quantity is the predicted Vss. The reduction
    # is therefore mono-exponential. It reproduces the paper's own
    # predicted Cmax and AUC for all three intravenous studies to within
    # 14 percent (Supplemental Table 1); that comparison is in the
    # vignette and nothing was tuned to reach it.
    # ------------------------------------------------------------------

    # Central volume. Table 3: 'Vss (L/kg) predicted' = 0.42. At the
    # 70 kg reference weight: 0.42 L/kg * 70 kg = 29.4 L.
    lvc <- fixed(log(29.4))
    label("Central compartment volume vc (L) at the 70 kg reference weight")

    # Renal arm of clearance. Table 4: 'CL R (L/h)' = 0.9; section 2.2:
    # renal clearance assumed to be 2% of the total in vivo clearance
    # (0.02 * 43 = 0.86, printed as 0.9 L/h).
    lcl_renal <- fixed(log(0.9))
    label("Renal plasma clearance CL_renal (L/h)")

    # Non-renal arm of clearance, by difference from the total
    # intravenous clearance of section 2.2 ('The clearance of 43 L/h was
    # the average clearance value after intravenous administration
    # obtained from 3 reports'), which the retrograde model was built to
    # reproduce: 43 - 0.9 = 42.1 L/h. It lumps the CYP2C8 and CYP2C9
    # hepatic arms, biliary clearance and the 2.6 L/h additional systemic
    # clearance of Table 4.
    lcl_nonren <- fixed(log(42.1))
    label("Non-renal plasma clearance CL_nonren (L/h)")

    # Wu 2022 is a PBPK simulation analysis, not a population-PK fit. It
    # reports no residual-error model and no inter-individual variance
    # components. Rather than invent a variance, the residual error is
    # fixed at zero, which makes this a deterministic typical-value
    # simulation model.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariates and no random effects.
    #    Total plasma clearance is the sum of the two arms and reproduces
    #    the 43 L/h of section 2.2 by construction.
    # ------------------------------------------------------------------
    vc         <- exp(lvc)
    cl_renal   <- exp(lcl_renal)
    cl_nonren  <- exp(lcl_nonren)
    cl         <- cl_renal + cl_nonren

    # ------------------------------------------------------------------
    # 2. Micro-constant for elimination.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system. A single compartment, as justified in the ini()
    #    header. Amounts are in mg. Intravenous administration only
    #    (infusions are given through the event table's rate or dur).
    # ------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ------------------------------------------------------------------
    # 4. Observation. Doses are in mg and vc is in L, so central / vc is
    #    in mg/L = ug/mL; multiply by 1000 to report ng/mL, the units of
    #    Wu 2022 Figure 1 and Supplemental Table 1.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

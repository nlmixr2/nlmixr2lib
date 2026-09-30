Tse_2020_tofacitinib <- function() {
  description <- paste(
    "One-compartment oral and intravenous pharmacokinetic reduction of the",
    "Simcyp (version 15, release 1) minimal-PBPK model for the Janus kinase",
    "inhibitor tofacitinib in healthy adults (Tse 2020). The source model",
    "uses first-order absorption, the Simcyp minimal-PBPK distribution",
    "option with no single adjusting compartment and a user-entered",
    "steady-state volume of 1.24 L/kg, and a clearance built from the",
    "observed intravenous clearance (24.7 L/h) split into a renal arm",
    "(7.62 L/h) and a hepatic CYP3A4 + CYP2C19 arm. With no adjusting",
    "compartment the plasma profile is one-compartmental, so the model is",
    "encoded as depot + central with the reported ka, Vss, renal and",
    "non-renal clearance, and the reported oral bioavailability of 74%.",
    "No parameter is fitted. The reduction reproduces the paper's own",
    "predicted intravenous Cmax to 0.4%, the predicted oral Cmax and",
    "AUC0-inf at all six single doses from 1 to 100 mg to within 7%, and",
    "the predicted 14% AUC increase when active renal secretion is",
    "abolished (see the validation vignette). This is a typical-value",
    "simulation model: the source reports no inter-individual variance",
    "and no residual-error model. The drug-interaction (fluconazole,",
    "ketoconazole, rifampicin) and renal / hepatic impairment predictions",
    "depend on proprietary Simcyp compound and population files and are",
    "not reproducible from this model.",
    sep = " "
  )
  reference <- paste(
    "Tse S, Dowty ME, Menon S, Gupta P, Krishnaswami S. (2020).",
    "Application of Physiologically Based Pharmacokinetic Modeling to",
    "Predict Drug Exposure and Support Dosing Recommendations for",
    "Potential Drug-Drug Interactions or in Special Populations: An",
    "Example Using Tofacitinib.",
    "J Clin Pharmacol 60(12):1617-1628. doi:10.1002/jcph.1679.",
    sep = " "
  )
  vignette <- "Tse_2020_tofacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales the distribution volume linearly, because the source enters",
        "Vss per kilogram (Table 1, 'Vss (L/kg) mode: User' = 1.24). The",
        "Simcyp virtual population's weights are not reported, so 70 kg is",
        "used as the reference weight. Clearance is NOT weight-scaled: the",
        "source enters it as an absolute systemic clearance in L/h."
      ),
      source_name = "body weight (Simcyp virtual population)"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance / glomerular filtration rate",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Renal impairment is simulated in the source by swapping in the",
        "Simcyp GFR 30-60 and GFR < 30 population files (Table 4,",
        "Supplemental Table S3). The resulting change in renal clearance is",
        "internal to those population files and is not reported, so no",
        "renal-function covariate relationship is reproducible here."
      )
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (Child-Pugh A)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Simulated in the source with the Simcyp Liver Cirrhosis CP-A",
        "population file (Table 4). The physiological changes that file",
        "applies are not reported."
      )
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (Child-Pugh B)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Simulated in the source with the Simcyp Liver Cirrhosis CP-B",
        "population file (Table 4). The physiological changes that file",
        "applies are not reported."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 9L,
    age_range = "19-54 years across the healthy-volunteer simulations (matched to each clinical study)",
    weight_median = "70 kg (reference weight for the L/kg volume input; the Simcyp population weights are not reported)",
    sex_female_pct = 0,
    disease_state = "Healthy adult volunteers.",
    dose_range = paste(
      "Single 10 mg intravenous infusion over 0.5 h; single oral doses",
      "of 1 to 100 mg; 15 mg orally twice daily for 14 days."
    ),
    regions = "Simcyp Healthy Volunteers population file (version 15, release 1).",
    studies = paste(
      "Model development: absolute bioavailability study (10 mg IV and",
      "PO, n = 12; Gupta 2011) and the 14C mass-balance study (10 mg PO,",
      "n = 6; Dowty 2014). Verification: single ascending-dose studies",
      "(1 mg, n = 6, Suzuki 2017; 3 to 100 mg, n = 7-9 per dose,",
      "Krishnaswami 2015) and a 15 mg twice-daily multiple-dose study",
      "(n = 23; Lawendy 2009). Supplemental Table S1."
    ),
    notes = paste(
      "This is a PBPK simulation analysis, not a population-PK fit.",
      "n_subjects records the virtual cohort size of every source",
      "simulation (10 trials x 10 subjects, 100% male for the",
      "healthy-volunteer scenarios); n_studies counts the clinical",
      "studies in Supplemental Table S1 with data used for development or",
      "verification (ascending-dose cohorts counted separately). The DDI",
      "(Supplemental Table S2) and organ-impairment (Supplemental",
      "Table S3) simulations are outside the scope of this model."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Every parameter is fixed: nothing was estimated in building this
    # reduction. Values are Tse 2020 Table 1 / Supplemental Figure S1
    # inputs, or arithmetic consequences of them.
    # ------------------------------------------------------------------

    lka <- fixed(log(5.7))
    label("First-order absorption rate constant ka (1/h)")
    # Tse 2020 Table 1, 'ka (per h)' = 5.7, estimated with the Simcyp
    # Parameter Estimation and Automated Sensitivity Analysis modules.

    # Distribution volume. Table 1 enters Vss = 1.24 L/kg (user mode,
    # from the absolute-bioavailability study, Gupta 2011) with the
    # minimal-PBPK distribution model and NO single adjusting
    # compartment, so Vss is the one-compartment plasma volume:
    #   1.24 L/kg * 70 kg = 86.8 L
    lvc <- fixed(log(86.8))
    label("Central volume vc at the 70 kg reference weight (L)")
    # Tse 2020 Table 1, 'Vss (L/kg) b mode: User' = 1.24, times 70 kg.

    lcl_renal <- fixed(log(7.62))
    label("Renal clearance (L/h)")
    # Tse 2020 Table 1, 'Cl R (L/h)' = 7.62 (mean observed renal
    # clearance across clinical studies). The Methods text also writes
    # it as 'ClIV x 0.29', which evaluates to 7.16 L/h; the Table 1
    # value is used because it is the model input and because the
    # renal-secretion scenario (ClIV 24.7 -> 21.68 L/h when ClR goes
    # 7.62 -> 4.6 L/h) closes only with 7.62.

    # Non-renal (hepatic CYP3A4 + CYP2C19) clearance. Table 1 enters the
    # observed intravenous clearance ClIV = 24.7 L/h (footnote c) and
    # retrograde-calculates the CYP intrinsic clearances from the
    # remainder after renal clearance:
    #   24.7 - 7.62 = 17.08 L/h
    # Supplemental Figure S1 splits total clearance as CYP3A4 54%,
    # CYP2C19 17% and renal 29% (passive 18%, active 11%). The CYP split
    # matters only for the DDI simulations, whose perpetrator models are
    # proprietary Simcyp compound files, so it is not carried here.
    lcl_nonren <- fixed(log(17.08))
    label("Non-renal (hepatic CYP3A4 + CYP2C19) clearance (L/h)")
    # Tse 2020 Table 1 footnote c (ClIV 24.7 L/h) minus Table 1 ClR.

    # Oral bioavailability. The source model obtains F as
    # fa * fg * fh with fa = 0.93 (Table 1) and fg, fh computed inside
    # Simcyp; the baseline fg and fh values are not printed. The
    # disposition scheme the model was built to reproduce
    # (Supplemental Figure S1: 93% absorbed, 74% reaching the systemic
    # circulation) gives F = 0.74, the same figure the Discussion cites.
    lfdepot <- fixed(log(0.74))
    label("Oral bioavailability F (fraction)")
    # Tse 2020 Supplemental Figure S1 ('74%', footnote b Gupta 2011)
    # and Discussion ('high oral bioavailability of 74%').

    # Tse 2020 is a PBPK simulation analysis, not a population-PK fit;
    # the SDs in Tables 2-4 describe a Simcyp virtual population driven
    # by unpublished population files, not estimated omegas. The
    # residual error is fixed at zero rather than invented.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model reported by the source)")
  })

  model({
    # 1. Individual parameters. Body weight scales the volume only.
    ka <- exp(lka)
    vc <- exp(lvc) * WT / 70
    cl_renal <- exp(lcl_renal)
    cl_nonren <- exp(lcl_nonren)
    fdepot <- exp(lfdepot)

    # 2. Total systemic clearance = renal + non-renal (Table 1).
    cl <- cl_renal + cl_nonren
    kel <- cl / vc

    # 3. ODEs. The source's portal-vein and liver compartments carry
    # first-pass extraction, which is represented by the bioavailability
    # on the depot; their volumes are not reported and are lumped into
    # the systemic compartment. Amounts are in mg.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Bioavailability on the oral depot. Intravenous doses go to
    # central and are unaffected.
    f(depot) <- fdepot

    # 5. Observation: mg / L * 1000 = ng/mL (the units of Tables 2-4).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

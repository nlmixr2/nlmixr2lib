Gidal_2018_eslicarbazepine_ssf <- function() {
  description <- paste0(
    "Emax exposure-efficacy model for the total STANDARDIZED SEIZURE ",
    "FREQUENCY (SSF, seizures per 28 days) during the maintenance period ",
    "in adults with focal-onset seizures taking adjunctive ",
    "eslicarbazepine acetate (ESL) (Gidal 2018, phase 3 trials 2093-301, ",
    "2093-302 and 2093-304). The model works on the log scale with a ",
    "0.33 offset that keeps zero-seizure patients finite: ",
    "ln(SSF + 0.33) = BSLN - 0.276*placebo + (1 - placebo) * ",
    "[Emax * Cav-ss / (3530 + Cav-ss)], where ",
    "BSLN = 2.19 + 0.228*WesternEurope + 0.310*LatinAmerica ",
    "+ 0.460*NorthAmerica - 0.00922*(age - 37) and ",
    "Emax = -0.822 + 0.150*baselineCBZ + 0.242*WesternEurope ",
    "(Gidal 2018 Appendix S1 Eqs. E-7, E-8 and E-9, Table S-7). Note the ",
    "switch structure: the constant placebo effect applies ONLY to the ",
    "placebo arm and the Emax drug effect ONLY to the active arms, so ",
    "the two are mutually exclusive rather than additive. Emax is ",
    "negative because the effect is a reduction in seizure frequency, ",
    "and both retained Emax covariates SHRINK that reduction -- patients ",
    "taking carbamazepine at baseline and patients from Western Europe ",
    "are predicted to benefit less. Rest of World is the region ",
    "reference. The EC50 of 3,530 ng/mL is distinct from the 9,450 ng/mL ",
    "EC50 of the companion weekly-seizure-count model; the 9.5 ug/mL ",
    "value quoted in the paper's main text is the latter. Linear, ",
    "log-linear and saturable exposure forms were screened and the ",
    "saturable one was selected. There is no PK layer and no ODE: ",
    "exposure enters as the static per-patient column CAV, an ",
    "empirical-Bayes prediction from ",
    "modellib('Gidal_2018_eslicarbazepine'). Minimum objective function ",
    "918.822."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Parameter table and equations are in Appendix S1 (supporting",
    "information), Table S-7 and Equations E-7, E-8 and E-9.",
    "Exposure metric produced by modellib('Gidal_2018_eslicarbazepine').",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"
  units <- list(
    time = "n/a (landmark maintenance-period seizure count standardized to 28 days; the model carries no time term)",
    dosing = "n/a (no dose events; exposure enters as the covariate CAV)",
    concentration = "lnssf (ln of standardized seizure frequency plus 0.33; the natural-scale ssf, seizures per 28 days, is derived)"
  )

  covariateData <- list(
    CAV = list(
      description = paste(
        "Individual predicted average steady-state eslicarbazepine plasma",
        "concentration over the once-daily dosing interval."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical-Bayes prediction from the population PK model of the",
        "same paper; compute it as dose / (24 * CL/F) with",
        "modellib('Gidal_2018_eslicarbazepine'), equivalently AUC0-24 / 24.",
        "The ng/mL scale is fixed by the EC50 of 3,530 ng/mL that it is",
        "compared against inside the Emax term and by the companion",
        "probability-of-response model, which centres the same metric at a",
        "median of 10,205 ng/mL. Set to 0 for placebo patients; the",
        "(1 - PLACEBO) switch already removes the drug term for them, so",
        "the two encodings agree. NOT centred and NOT scaled -- it enters",
        "raw into the Emax denominator."
      ),
      source_name = "C_av-ss_i (Eq. E-7)"
    ),
    PLACEBO = list(
      description = "Randomised placebo-arm membership; 1 = placebo, 0 = active eslicarbazepine acetate.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (active treatment arm)",
      notes = paste(
        "Acts as an exclusive SWITCH rather than as an ordinary additive",
        "covariate. Eq. E-7 multiplies the constant placebo effect by",
        "PLACEBO and the Emax drug term by (1 - PLACEBO), so a placebo",
        "patient receives the -0.276 shift and NO drug term, while an",
        "active patient receives the drug term and NO placebo shift. That",
        "is the structure the paper printed and it is reproduced",
        "literally here. A consequence worth noting when simulating: the",
        "model does not predict a continuous approach to the placebo",
        "response as Cav-ss goes to zero -- an active patient with zero",
        "exposure sits 0.276 log units above a placebo patient."
      ),
      source_name = "plac_i (Eq. E-7)"
    ),
    REGION_WESTERNEUROPE = list(
      description = "Western European study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with REGION_LATINAMERICA and REGION_NORTHAMERICA also 0 this selects the REST OF WORLD reference group",
      notes = paste(
        "The only covariate in this model that acts on BOTH the baseline",
        "and the Emax: +0.228 on BSLN (Eq. E-8) and +0.242 on Emax",
        "(Eq. E-9). Because Emax is negative, a positive shift SHRINKS the",
        "seizure reduction, so Western European patients are predicted to",
        "have both a higher baseline seizure frequency and a smaller",
        "treatment benefit. Gidal 2018 Discussion attributes this to",
        "demographic and clinical differences between regions and notes",
        "the study was not designed to test it. The three region",
        "indicators are mutually exclusive; all three 0 selects Rest of",
        "World. NOTE this four-way split (Western Europe / Latin America /",
        "North America against a Rest-of-World reference) differs from the",
        "split used by the paper's TEAE models, where Europe is the",
        "reference -- the two are not interchangeable."
      ),
      source_name = "WEU_i (Eqs. E-8 and E-9)"
    ),
    REGION_LATINAMERICA = list(
      description = "Latin American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; see REGION_WESTERNEUROPE for the shared Rest-of-World reference",
      notes = paste(
        "Acts on the BASELINE only (+0.310 on BSLN, Eq. E-8); it was not",
        "retained on Emax. Mutually exclusive with the other two region",
        "indicators."
      ),
      source_name = "LA_i (Eq. E-8)"
    ),
    REGION_NORTHAMERICA = list(
      description = "North American study site; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; see REGION_WESTERNEUROPE for the shared Rest-of-World reference",
      notes = paste(
        "Acts on the BASELINE only (+0.460 on BSLN, Eq. E-8), the largest",
        "of the three regional baseline shifts; it was not retained on",
        "Emax. Study 2093-304 was the North American trial. Mutually",
        "exclusive with the other two region indicators."
      ),
      source_name = "NA_i (Eq. E-8)"
    ),
    AGE = list(
      description = "Patient age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the BASELINE linearly, centred at the population median of",
        "37 years; the negative slope means older patients are predicted",
        "to have a lower baseline seizure frequency. The figures of Gidal",
        "2018 are all drawn at this median age. Not retained on Emax."
      ),
      source_name = "Age_i (Eq. E-8)"
    ),
    CONMED_CBZ = list(
      description = "Carbamazepine use during the baseline period; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (took other antiepileptic drugs at baseline)",
      notes = paste(
        "Acts on the EMAX only (+0.150, Eq. E-9); it was not retained on",
        "the baseline. Because Emax is negative, the positive shift shrinks",
        "the predicted seizure reduction by about 18%, so baseline",
        "carbamazepine users benefit less. About 49% of the analysis",
        "population was on carbamazepine. Gidal 2018 Discussion presents",
        "this as confirmation of earlier post hoc analyses of the same",
        "trials, and combines it with the PK finding that carbamazepine",
        "raises eslicarbazepine CL/F to recommend an ESL dose increase",
        "when the two are taken together."
      ),
      source_name = "BCAR_i (Eq. E-9)"
    )
  )

  population <- list(
    species = "human",
    n_studies = 3L,
    n_subjects = 1152L,
    age_median = "37 years (the centring value used by Eq. E-8)",
    disease_state = "adults with focal-onset seizures, at least four in the 4 weeks before screening despite 1-3 concomitant antiepileptic drugs",
    dose_range = "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily, or placebo",
    regions = "Rest of World (reference), Western Europe, Latin America and North America",
    baseline_seizure_frequency = "2-412 seizures per 28 days during the 8-week baseline period (Gidal 2018 Results)",
    co_medication = "about 49% were receiving concomitant carbamazepine",
    notes = paste(
      "The paper reports model-predicted standardized seizure frequencies",
      "of 6.5 (placebo), 5.4 (400 mg), 4.6 (800 mg) and 4.3 (1,200 mg)",
      "seizures per 28 days. The placebo value reproduces exactly from",
      "Eqs. E-7 and E-8 at the reference patient: exp(2.19 - 0.276) - 0.33",
      "= 6.45. The active-arm values reproduce to within rounding when",
      "Cav-ss is taken as the median for each dose implied by the",
      "companion probability-of-response model's 10,205 ng/mL centring at",
      "800 mg; see the vignette source trace. The number of subjects is",
      "the safety analysis set size, which Appendix S1 does not restate",
      "separately for the seizure-frequency analysis."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-7 and Equations E-7, E-8, E-9:
    #
    #   ln(SSF + 0.33)_i = BSLN_i - 0.276*plac_i
    #                      + (1 - plac_i) * [Emax_i * Cav / (3530 + Cav)]
    #   BSLN_i = 2.19 + 0.228*WEU + 0.310*LA + 0.460*NA
    #                 - 0.00922*(Age - 37)
    #   Emax_i = -0.822 + 0.150*BCAR + 0.242*WEU
    #
    # The 0.33 offset inside the log is part of the transform the authors
    # fitted, not a nuisance constant: it keeps ln finite for patients
    # with zero seizures during maintenance.
    # ==================================================================

    # ----- Baseline (Eq. E-8), on the ln(SSF + 0.33) scale -----
    lrbase <- 2.19 ; label("Baseline ln(standardized seizure frequency + 0.33) for a 37-year-old Rest-of-World patient (unitless, log scale)")  # Table S-7, 'Baseline SSF' 2.19, 1.3% SEM; Eq. E-8. exp(2.19) - 0.33 = 8.61 seizures per 28 days
    e_region_westerneurope_rbase <- 0.228 ; label("Shift in baseline ln(SSF + 0.33) for a Western European versus a Rest-of-World study site (unitless, log scale)")  # Table S-7, 'Effect of Western European region on baseline SSF' 0.228, 30.6% SEM; Eq. E-8
    e_region_latinamerica_rbase <- 0.310 ; label("Shift in baseline ln(SSF + 0.33) for a Latin American versus a Rest-of-World study site (unitless, log scale)")  # Table S-7, 'Effect of Latin American region on baseline SSF' 0.310, 18.4% SEM; Eq. E-8
    e_region_northamerica_rbase <- 0.460 ; label("Shift in baseline ln(SSF + 0.33) for a North American versus a Rest-of-World study site (unitless, log scale)")  # Table S-7, 'Effect of North American region on baseline SSF' 0.460, 15.3% SEM; Eq. E-8
    e_age_rbase <- -0.00922 ; label("Change in baseline ln(SSF + 0.33) per year of age above 37 years (1/year, log scale)")  # Table S-7, 'Slope of age effect on baseline SSF' -0.00922, 18.3% SEM; Eq. E-8

    # ----- Placebo effect (Eq. E-7), applied to the placebo arm only ----
    e_placebo_lnssf <- -0.276 ; label("Constant shift in ln(SSF + 0.33) during maintenance for patients randomised to placebo (unitless, log scale)")  # Table S-7, 'Constant placebo effect' -0.276, 13.1% SEM; Eq. E-7. Check: exp(2.19 - 0.276) - 0.33 = 6.45, the 6.5 seizures per 28 days the paper predicts for placebo

    # ----- Emax drug effect (Eqs. E-7 and E-9) -----
    # Bare emax, not lemax: the value is NEGATIVE (a reduction in seizure
    # frequency on the log scale) and cannot be log-transformed.
    emax <- -0.822 ; label("Maximum reduction in ln(SSF + 0.33) attributable to eslicarbazepine, for a patient not taking baseline carbamazepine outside Western Europe (unitless, log scale)")  # Table S-7, 'Emax at the baseline SSF of 2.4' -0.822, 13.9% SEM; Eq. E-9. The table's row label quotes the baseline SSF at which the parameter was estimated; Eq. E-9 carries NO baseline-dependent scaling and is what is encoded here
    e_conmed_cbz_emax <- 0.150 ; label("Shift in the eslicarbazepine Emax for baseline carbamazepine use (unitless, log scale)")  # Table S-7, 'Effect of baseline carbamazepine use on Emax' 0.150, 37.5% SEM; Eq. E-9. Positive, so it shrinks the negative Emax and reduces the benefit
    e_region_westerneurope_emax <- 0.242 ; label("Shift in the eslicarbazepine Emax for a Western European versus a Rest-of-World study site (unitless, log scale)")  # Table S-7, 'Effect of Western European region on Emax' 0.242, 27.9% SEM; Eq. E-9. Positive, so it shrinks the benefit
    lec50 <- log(3530) ; label("Average steady-state eslicarbazepine concentration giving half the maximum seizure-frequency reduction (ng/mL)")  # Table S-7, 'EC50 (ng/mL)' 3,530, 51.0% SEM; Eq. E-7. Distinct from the 9,450 ng/mL EC50 of the companion weekly-seizure model, which is the value quoted as 9.5 ug/mL in the main text

    # ----- Interindividual variability -----
    # Table S-7 reports these three as VARIANCES on the additive log
    # scale (footnotes a, b and c), not as %CV. The footnote SDs confirm
    # each: sqrt(0.544) = 0.74, sqrt(0.503) = 0.71, sqrt(1.80) = 1.34
    # (printed as '134 %CV').
    etalrbase ~ 0.544          # Table S-7, baseline SSF IIV 0.544, 5.4% SEM; footnote a: 'The estimate presented is a variance. The IIV in the baseline SSF is described by an SD of 0.74'
    etae_placebo_lnssf ~ 0.503 # Table S-7, placebo effect IIV 0.503, 14.1% SEM; footnote b: 'The estimate presented is a variance. The IIV in the placebo effect is described by an SD of 0.71'
    etaemax ~ 1.80             # Table S-7, Emax IIV 1.80, 11.1% SEM; footnote c: 'The estimate presented is a variance. The IIV in the Emax is described by a 134 %CV' -- sqrt(1.80) = 1.342

    # ----- Residual error -----
    addSd_lnssf <- 0.1597 ; label("Additive residual SD on ln(standardized seizure frequency + 0.33) (unitless, log scale)")  # Table S-7, 'Additive RV' 0.0255 as a variance, 58.0% SEM; footnote d: 'The estimate presented is a variance. The RV is described by an SD of 0.16'. sqrt(0.0255) = 0.15969
  })

  model({
    # ----- Baseline on the ln(SSF + 0.33) scale (Eq. E-8) -----
    # IIV is additive on this scale because the scale is already
    # logarithmic; the printed variances are for exactly this form.
    bsln <- lrbase +
      e_region_westerneurope_rbase * REGION_WESTERNEUROPE +
      e_region_latinamerica_rbase * REGION_LATINAMERICA +
      e_region_northamerica_rbase * REGION_NORTHAMERICA +
      e_age_rbase * (AGE - 37) +
      etalrbase

    # ----- Placebo effect, placebo arm only (Eq. E-7) -----
    plac_eff <- e_placebo_lnssf + etae_placebo_lnssf

    # ----- Emax drug effect, active arms only (Eqs. E-7 and E-9) -----
    emax_i <- (emax +
      e_conmed_cbz_emax * CONMED_CBZ +
      e_region_westerneurope_emax * REGION_WESTERNEUROPE) *
      exp(etaemax)
    ec50 <- exp(lec50)
    drug_eff <- emax_i * CAV / (ec50 + CAV)

    # ----- Prediction (Eq. E-7) -----
    # The PLACEBO switch is exclusive: the placebo shift and the drug
    # term never both contribute to the same subject.
    lnssf <- bsln + PLACEBO * plac_eff + (1 - PLACEBO) * drug_eff

    # ----- Natural-scale standardized seizure frequency -----
    # Inverts the ln(SSF + 0.33) transform the model was fitted on.
    ssf <- exp(lnssf) - 0.33

    # ----- Observation -----
    # The residual is additive on the transformed scale, matching the
    # paper's goodness-of-fit plots of measured versus predicted ln(SSF).
    lnssf ~ add(addSd_lnssf)
  })
}

Park_2023_everolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for everolimus in patients with refractory seizures associated with focal cortical dysplasia (FCD) type II. Apparent oral clearance rises LINEARLY (additively) with body surface area, CL/F = 12.5 + 9.71 * (BSA / 1.5) L/h, where 1.5 m2 is the cohort median BSA; the apparent central volume (293 L) and the absorption rate constant (0.585 /h) carry no covariate. Inter-individual variability on the absorption rate constant was fixed to zero by the authors. Everolimus is assayed in whole blood and concentrations are reported in ng/mL."
  reference <- paste(
    "Park J, Kim SH, Hahn J, Kang H-C, Lee S-G, Kim HD, Chang MJ.",
    "Population pharmacokinetics of everolimus in patients with seizures",
    "associated with focal cortical dysplasia.",
    "Front Pharmacol. 2023;14:1197549. doi:10.3389/fphar.2023.1197549.",
    sep = " "
  )
  vignette <- "Park_2023_everolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Everolimus partitions extensively into erythrocytes and
  # the Park 2023 assay was run on whole blood (Methods 2.3; the authors cite
  # the whole-blood assay as their reason for rejecting RBC count as a
  # covariate in Results 3.2), so `central` is a whole-blood compartment.
  compartmentData <- list(
    depot   = list(analyte = "everolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "everolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model (Results 3.2:",
        "forward selection dOFV = -6.356; backward elimination dOFV = +6.92,",
        "above the 6.63 retention threshold). Enters apparent oral clearance",
        "as a LINEAR / additive term rather than the more usual multiplicative",
        "power form: TVCL = 12.5 + 9.71 * (BSA / 1.5) L/h (Table 2 header and",
        "Abstract). 1.5 m2 is the cohort median BSA (Table 1), so the",
        "normalisation is median-centred and CL/F at the median BSA is",
        "12.5 + 9.71 = 22.21 L/h. Observed range 0.6-2 m2 (Table 1); the",
        "authors simulated BSA of 0.5, 0.7, 1, 1.5, 1.7 and 2 m2",
        "(Supplementary Table S2). Baseline value, not time-varying.",
        "The authors note the linear form agrees with the FDA clinical",
        "pharmacology review (CDER 2009), in which everolimus clearance also",
        "increases linearly with BSA."
      ),
      source_name        = "BSA"
    )
  )

  # Screened in the stepwise covariate search (Methods 2.4 lists 14 candidate
  # covariates) but NOT retained in the final model, so none is referenced in
  # model(). Recorded here to preserve the provenance of the covariate screen.
  # Note: the three CYP3A4 co-medication class flags below are documentation-
  # only labels for this paper's screened predictors; they are not entries in
  # inst/references/covariate-columns.md and mint no canonical name.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.", units = "years", type = "continuous",
      notes = "Screened (Methods 2.4); not retained. Median 13.5 y, range 4-32 y (Table 1).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator.", units = "(binary)", type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened (Methods 2.4); not retained. The paper codes sex as",
        "0 = male, 1 = female, which already matches the canonical SEXF",
        "polarity, so no value inversion would be required. 9 of 22 (40.9%)",
        "were male (Table 1)."
      ),
      source_name = "sex"
    ),
    WT = list(
      description = "Body weight.", units = "kg", type = "continuous",
      notes = "Screened (Methods 2.4); not retained. Median 50 kg, range 13-86 kg (Table 1).",
      source_name = "weight"
    ),
    ALB = list(
      description = "Serum albumin.", units = "g/dL", type = "continuous",
      notes = paste(
        "Screened (Methods 2.4). Selected in the SECOND forward-selection step",
        "as a proportional effect on clearance (dOFV = -4.941) and carried in",
        "the full model, but dropped at backward elimination because removing",
        "it raised OFV by only 4.941, below the 6.63 retention threshold",
        "(Results 3.2). The authors rationalised the candidate effect by",
        "everolimus's 75% protein binding. No point estimate for the albumin",
        "coefficient is reported anywhere in the paper or supplement, so the",
        "full model cannot be reconstructed. Median 4.6 g/dL, range 4.1-5.2",
        "g/dL (Table 1)."
      ),
      source_name = "albumin"
    ),
    CREAT = list(
      description = "Serum creatinine.", units = "mg/dL", type = "continuous",
      notes = "Screened (Methods 2.4); not retained. Median 0.67 mg/dL, range 0.33-0.92 mg/dL (Table 1).",
      source_name = "serum creatinine"
    ),
    ALT = list(
      description = "Alanine aminotransferase.", units = "IU/L", type = "continuous",
      notes = "Screened (Methods 2.4); not retained. Median 11 IU/L, range 5-67 IU/L (Table 1).",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase.", units = "IU/L", type = "continuous",
      notes = "Screened (Methods 2.4); not retained. Median 16 IU/L, range 10-44 IU/L (Table 1).",
      source_name = "AST"
    ),
    HGB = list(
      description = "Blood hemoglobin concentration.", units = "g/dL", type = "continuous",
      notes = "Screened (Methods 2.4); not retained. Median 13.15 g/dL, range 10.3-16.5 g/dL (Table 1).",
      source_name = "hemoglobin"
    ),
    HCT = list(
      description = "Hematocrit.", units = "%", type = "continuous",
      notes = "Screened (Methods 2.4); not retained. Median 38.95%, range 35-47.8% (Table 1).",
      source_name = "hematocrit"
    ),
    RBC = list(
      description = "Red blood cell count.", units = "10^6 cells/uL", type = "continuous",
      notes = paste(
        "Screened (Methods 2.4). A proportional effect of RBC on clearance",
        "(dOFV = -5.6) and a power effect of RBC on volume (dOFV = -4.095)",
        "both reduced the OFV, but the authors deliberately EXCLUDED RBC on",
        "mechanistic grounds: the everolimus bioassay is run on whole blood,",
        "so erythrocyte count does not influence the measured whole-blood",
        "concentration (Results 3.2). Median 4.48, range 3.89-5.19",
        "10^6 cells/uL (Table 1)."
      ),
      source_name = "RBC"
    ),
    CYP3A4_STRONG_INDUCER = list(
      description = "Presence of at least one concomitant strong CYP3A4 inducer (carbamazepine, phenytoin).",
      units = "(binary)", type = "binary", reference_category = "0 (absent)",
      notes = paste(
        "Screened (Methods 2.4); not retained. Only 1 of 22 patients (4.5%)",
        "took a strong inducer, and the one patient on a moderate inducer",
        "dropped out without valid concentrations (Results 3.1), so the",
        "authors state the inducer/inhibitor effect could not be estimated",
        "(Discussion; Limitations). Inducer list in Supplementary Table S1."
      ),
      source_name = "CYP3A4 strong inducer"
    ),
    CYP3A4_WEAK_INDUCER = list(
      description = "Presence of at least one concomitant weak CYP3A4 inducer (rufinamide, topiramate, clobazam, perampanel, oxcarbazepine).",
      units = "(binary)", type = "binary", reference_category = "0 (absent)",
      notes = "Screened (Methods 2.4); not retained. 14 of 22 (63.6%) exposed (Table 1); list in Supplementary Table S1.",
      source_name = "CYP3A4 weak inducer"
    ),
    CYP3A4_WEAK_INHIBITOR = list(
      description = "Presence of at least one concomitant weak CYP3A4 inhibitor (valproic acid, perampanel, ranitidine).",
      units = "(binary)", type = "binary", reference_category = "0 (absent)",
      notes = "Screened (Methods 2.4); not retained. 14 of 22 (63.6%) exposed (Table 1); list in Supplementary Table S1.",
      source_name = "CYP3A4 weak inhibitor"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    n_observations = "152 everolimus whole-blood concentrations (Results 3.2).",
    age_range      = "4-32 years",
    age_median     = "13.5 years (IQR 12-17.75; Table 1)",
    weight_range   = "13-86 kg",
    weight_median  = "50 kg (IQR 35.5-60; Table 1)",
    bsa_range      = "0.6-2 m^2",
    bsa_median     = "1.5 m^2 (IQR 1.23-1.7; Table 1)",
    sex_female_pct = 59.1,
    race_ethnicity = "Not reported. Single-centre Korean cohort (Severance Hospital, Seoul).",
    disease_state  = paste(
      "Focal cortical dysplasia (FCD) type II with refractory seizures",
      "despite more than two antiepileptic drugs; at least three seizures per",
      "month over two months, and no response to vagus nerve stimulation or",
      "dietary treatment (Supplementary Material inclusion criteria).",
      "All patients continued at least one concomitant antiepileptic drug,",
      "unchanged through the baseline and core phases."
    ),
    dose_range     = paste(
      "Everolimus (Afinitor disperz, tablet for oral suspension) once daily.",
      "Initial dose 4.5 mg/m2/day, then TDM-guided adjustment in 2 mg steps",
      "to a target trough of 5-15 ng/mL (Methods 2.2)."
    ),
    regions        = "Republic of Korea (Severance Hospital, Seoul).",
    notes          = paste(
      "Data from a double-blinded, placebo-controlled, crossover randomised",
      "clinical trial (IRB 4-2017-0299; additional post-dose sampling under",
      "IRB 4-2019-1232), September 2017 to May 2020. Phases: 4-week baseline,",
      "12-week core I, 12-week core II (crossover), 29-week extension.",
      "Sampling was sparse and predominantly trough-only: pre-dose at weeks",
      "2, 3, 4 and 8 (core I), 14, 15, 16 and 20 (core II) and 25, 28 and 40",
      "(extension), plus additional 1-4 h post-dose samples at weeks 24, 28",
      "or 40 (Methods 2.3). Assay HPLC-MS/MS on whole blood, LLOQ 1.1 ng/mL,",
      "validated 1.1-41.6 ng/mL, inter- and intra-assay CV below 7.3%.",
      "The sparse sampling is the authors' stated reason for preferring a",
      "one-compartment structure over the two-compartment models reported for",
      "everolimus elsewhere (Discussion)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Final-model parameter estimates, Park 2023 Table 2 ("Final model"
    # column) and the Table 2 header equation
    #   TVCL = CL + theta_BSA_on_CL * (BSA / 1.5);  TVV = V;  TVKA = KA
    # The same three values are restated in the Abstract Results sentence
    # ("TVCL = 12.5 + 9.71 x (BSA/1.5), TVV = 293, and TVKA = 0.585").
    # The "Structural model" column of Table 2 (CL 21.9, V 302, KA 0.573)
    # is the covariate-free base model, an intermediate development step;
    # this file packages the final model, per the base-vs-final policy.
    # ------------------------------------------------------------------
    lka      <- log(0.585); label("Absorption rate constant Ka (1/h)")                                        # Table 2 final KA = 0.585 /h (RSE 30%)
    lcl      <- log(12.5);  label("BSA-independent component of apparent oral clearance CL/F (L/h)")          # Table 2 final CL = 12.5 L/h (RSE 26%)
    e_bsa_cl <- 9.71;       label("Linear slope of BSA on apparent oral clearance CL/F, per unit BSA/1.5 (L/h)") # Table 2 final theta_BSA on CL = 9.71 (RSE 33%)
    lvc      <- log(293);   label("Apparent central volume of distribution V/F (L)")                          # Table 2 final V = 293 L (RSE 30%)

    # ------------------------------------------------------------------
    # Inter-individual variability, exponential (Methods 2.4:
    # theta_i = theta_POP * exp(eta_i)).
    #
    # Table 2 tabulates omega on the STANDARD-DEVIATION scale, not as a
    # variance. Proof from the paper's own numbers: Results 3.2 gives the
    # base-model variances as omega_CL^2 = 0.0409 and omega_V^2 = 0.0602,
    # and the base-model Table 2 cells are 0.2022 = sqrt(0.0409) and
    # 0.2454 = sqrt(0.0602). The final-model cells (0.1652, 0.2083) are
    # therefore also SDs, and are squared here to the variance scale that
    # nlmixr2 expects.
    # ------------------------------------------------------------------
    etalcl ~ 0.02729  # Table 2 final omega_CL = 0.1652 (log-scale SD) -> variance 0.1652^2 = 0.02729
    etalvc ~ 0.04339  # Table 2 final omega_V  = 0.2083 (log-scale SD) -> variance 0.2083^2 = 0.04339

    # The authors did not estimate IIV on ka (Results 3.2, omega_KA^2), because
    # "there was not much information about the absorption phase from our
    # data ... with the assumption that absorption rates are the same among
    # the patients" (Discussion).
    etalka ~ fixed(0) # Results 3.2 / Discussion: no absorption-phase information in the sparse data

    # ------------------------------------------------------------------
    # Residual error: proportional only (Results 3.2, "The residual
    # variability was best explained by the proportional error model").
    # Table 2 reports sigma on the same SD scale as the omega rows above,
    # so 0.2943 is a 29.4% proportional error. Confirmed against
    # Supplementary Table S3: the simulated trough CV there is ~50%, which
    # is what IIV (omega as SD) plus a 29.4% proportional residual gives;
    # reading these cells as variances would put the CV far above 50%.
    # ------------------------------------------------------------------
    propSd <- 0.2943; label("Proportional residual error (fraction)")  # Table 2 final sigma_proportional = 0.2943 (RSE 14%)
  })

  model({
    # Absorption. etalka is fixed at zero (see ini()), so ka is the same
    # 0.585 /h in every subject; the eta is carried explicitly to record
    # that the authors modelled and then fixed it rather than omitting it.
    ka <- exp(lka + etalka)

    # Apparent oral clearance. Park 2023 parameterises the BSA effect as a
    # LINEAR / ADDITIVE term, TVCL = CL + theta_BSA_on_CL * (BSA / 1.5),
    # not as the more common multiplicative power form -- so e_bsa_cl
    # carries units of L/h rather than being a unitless exponent. The
    # normalising constant 1.5 m2 is the cohort median BSA (Table 1); at
    # that BSA, CL/F = 12.5 + 9.71 = 22.21 L/h. The exponential IIV
    # multiplies the whole typical value (Methods 2.4).
    cl <- (exp(lcl) + e_bsa_cl * (BSA / 1.5)) * exp(etalcl)

    vc  <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Doses are in mg and vc is in L, so central / vc is mg/L; everolimus
    # whole-blood concentrations are reported in ng/mL (1 mg/L = 1000 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

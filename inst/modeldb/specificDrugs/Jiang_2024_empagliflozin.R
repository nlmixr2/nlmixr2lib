Jiang_2024_empagliflozin <- function() {
  description <- paste(
    "Two-compartment population PK model with Savic transit-compartment",
    "absorption for oral empagliflozin in healthy Korean adults (Jiang",
    "2024). Pools two randomized, open-label, two-period crossover phase 1",
    "studies that compared a novel empagliflozin L-proline cocrystal",
    "formulation (CKD-370) against the conventional empagliflozin",
    "formulation at 25 mg, and as a 5 mg/1000 mg empagliflozin/metformin",
    "fixed-dose combination. An analytical transit chain (Ktr = 9.19 1/h,",
    "Mtt = 0.63 h, derived N = 4.79) feeds first-order absorption (Ka) into",
    "a two-compartment disposition model. Log-transformed body weight",
    "scales apparent clearance (exponent 0.64) and the apparent peripheral",
    "volume (exponent 0.57). Interindividual variability on Ka, CL, and the",
    "peripheral volume; interoccasion variability across the two crossover",
    "periods on Ktr, Mtt, Ka, and CL. Formulation (cocrystal vs",
    "conventional) was screened as a covariate and was NOT retained: the",
    "paper's central finding is that the L-proline cocrystal does not alter",
    "the absorption phase or any PK parameter, so this single model",
    "describes both formulations.",
    sep = " "
  )
  reference <- paste(
    "Jiang X, Yu KS, Nam DH, Oh J. A Population Pharmacokinetic Study to",
    "Compare a Novel Empagliflozin L-Proline Formulation with Its",
    "Conventional Formulation in Healthy Subjects. Pharmaceuticals (Basel).",
    "2024 Apr 18;17(4):522. doi:10.3390/ph17040522.",
    sep = " "
  )
  vignette <- "Jiang_2024_empagliflozin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the Jiang 2024 final model. Entered",
        "as log-transformed body weight on apparent clearance and on the",
        "apparent peripheral volume of distribution (Results 2.3: 'the",
        "inclusion of log-transformed body weight on clearance (CL) and the",
        "volume of distribution in the peripheral compartment (V2) ...",
        "significantly improved the model fit'). Table 2 reports the",
        "coefficients as beta_CL_logWT = 0.64 (RSE 24.5%) and",
        "beta_V2_logWT = 0.57 (RSE 30.3%). A coefficient on",
        "log(WT / WT_ref) added to a log-scale parameter is algebraically",
        "identical to the power form (WT / WT_ref)^beta used in model(),",
        "so the model is encoded in the power form that is standard across",
        "nlmixr2lib. IMPORTANT: Jiang 2024 does NOT report the centering",
        "(reference) weight used for the log transformation, and Monolix",
        "does not fix it by convention. A reference of 70 kg is used here",
        "per the nlmixr2lib default for an unreported centering value; it",
        "also coincides with the pooled cohort mean weight of 69.69 kg",
        "(Table 1, n = 54). The pooled cohort MEDIAN of 72.1 kg is the other",
        "plausible reading; because the covariate enters as a power, the",
        "choice only rescales the typical values, and over the observed",
        "55.6-83.4 kg range the two readings differ by at most",
        "(72.1/70)^0.64 = 1.9% on CL and 1.7% on the peripheral volume.",
        "See the vignette 'Assumptions and deviations' section.",
        sep = " "
      ),
      source_name        = "Weight"
    ),
    OCC = list(
      description        = "Crossover-period (occasion) indicator: 1 = period 1, 2 = period 2",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Both source studies were randomized, open-label, two-period,",
        "two-sequence crossover studies, so every subject contributes two",
        "single-dose occasions (one per formulation). Jiang 2024 Discussion:",
        "'Since PK modeling was conducted in this study by using data from",
        "two crossover studies, the effects of period 1 and period 2 were",
        "reflected by the IOV.' Decomposed inside model() into binary",
        "indicators oc1 and oc2 that multiplex the per-occasion IOV etas on",
        "log-Ktr, log-Mtt, log-Ka, and log-CL. Table 2 reports a single IOV",
        "standard deviation per parameter (not one per occasion), so the",
        "occasion-2 variances are fix()'d equal to the occasion-1 variances,",
        "mirroring the NONMEM '$OMEGA BLOCK(1) SAME' idiom already used by",
        "Smythe_2013_gatifloxacin.R and Goggin_2004_emfilermin.R. Note that",
        "occasion is NOT the same thing as formulation: because formulation",
        "was screened and rejected as a covariate, period 1 and period 2",
        "differ only through the IOV random effects.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  # Covariates screened by Jiang 2024 but NOT retained in the final model.
  # Methods 4.4 lists all 13 tested covariates: "treatment, age, height,
  # weight, BMI, sex, glucose, protein, aspartate aminotransferase (AST),
  # alkaline phosphatase (ALP), alanine transaminase (ALT), lactate
  # dehydrogenase (LDH), and estimated glomerular filtration rate (eGFR)".
  # Only weight survived, so the other 12 are recorded here as documentation
  # (they are deliberately never referenced in model()). Ranges are from
  # Table 1, pooled column (n = 54).
  covariatesDataExcluded <- list(
    FORM_EMPA_LPROLINE = list(
      description        = "Formulation indicator: 1 = empagliflozin L-proline cocrystal (CKD-370), 0 = conventional empagliflozin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (conventional empagliflozin formulation)",
      notes              = paste(
        "Screened as the 'treatment' covariate (Methods 4.4) and NOT",
        "retained -- this null result is the paper's primary conclusion.",
        "Results 2.3 / Discussion: 'the results of the population PK",
        "analysis showed that formulation differences did not affect PK",
        "parameters as covariates, which suggested that the difference in",
        "the absorption phase is not significantly meaningful and can be",
        "ignored.' Consequently the packaged model is formulation-agnostic",
        "and predicts identical profiles for both products; the vignette",
        "asserts that invariance numerically. Retained here purely to",
        "document that the comparison was made. Study A compared 25 mg",
        "CKD-370 vs 25 mg conventional empagliflozin (NCT03849495); Study B",
        "compared 5 mg/1000 mg empagliflozin L-proline/metformin vs",
        "5 mg/1000 mg empagliflozin/metformin fixed-dose combinations",
        "(NCT03848637).",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (Methods 4.4), not retained. Table 1 pooled median 27.5 years [20-50], mean 30.26 +/- 7.48. Study eligibility required 19-50 years."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened (Methods 4.4), not retained. Table 1 pooled median 171.5 cm [155.2-186.4], mean 171.21 +/- 7.28."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened (Methods 4.4), not retained. Table 1 pooled median 24.1 kg/m^2 [19.5-26.9], mean 23.75 +/- 2.03. Study eligibility required 18-27 kg/m^2."
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the cohort majority, 42/54)",
      notes              = paste(
        "Screened (Methods 4.4). Results 2.3 reports that sex on the",
        "peripheral volume of distribution DID significantly improve the",
        "fit during covariate screening, but it 'was also excluded from the",
        "final model due to poor estimation precision (relative standard",
        "error greater than 70%)'. No point estimate is published for the",
        "sex effect, so it cannot be encoded even optionally. Cohort: 22",
        "male / 5 female in Study A, 20 male / 7 female in Study B.",
        sep = " "
      )
    ),
    FPG = list(
      description = "Fasting plasma glucose (screening clinical chemistry)",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened (Methods 4.4), not retained. Table 1 pooled median 87 mg/dL [78-101], mean 87.70 +/- 5.97. Reported in mg/dL by the source; the canonical FPG unit is mmol/L (divide by 18.02)."
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened as 'protein' (Methods 4.4), not retained. Table 1 pooled median 6.7 g/dL [6.2-7.4], mean 6.69 +/- 0.25. Reported in g/dL by the source; the canonical TPRO unit is g/L (multiply by 10)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.4), not retained. Table 1 pooled median 17 IU/L [12-35], mean 17.61 +/- 4.47."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "IU/L",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods 4.4). Results 2.3 reports that ALP on",
        "intercompartmental clearance (Q) DID significantly improve the fit",
        "during covariate screening, but 'ALP was excluded because it was",
        "not clinically relevant to Q'. No point estimate is published, so",
        "it cannot be encoded. Table 1 pooled median 55.5 IU/L [30-83],",
        "mean 56.54 +/- 13.18.",
        sep = " "
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.4), not retained. Table 1 pooled median 16.5 IU/L [7-64], mean 19.00 +/- 10.60."
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened (Methods 4.4), not retained. Table 1 pooled median 147 IU/L [113-219], mean 150.06 +/- 21.88."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (MDRD equation)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened as 'eGFR' (Methods 4.4), not retained. Table 1 pooled median 112.35 mL/min/1.73 m^2 [74.2-143.9], mean 111.17 +/- 15.62. Stored under the canonical CRCL, whose register entry lists MDRD eGFR as an accepted assay form."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species         = "human",
    n_subjects      = 54L,
    n_studies       = 2L,
    n_observations  = 1832L,
    age_range       = "20-50 years (Table 1 pooled: median 27.5 years, mean 30.26 +/- 7.48); study eligibility 19-50 years",
    weight_range    = "55.6-83.4 kg (Table 1 pooled: median 72.1 kg, mean 69.69 +/- 7.82)",
    bmi_range       = "19.5-26.9 kg/m^2 (Table 1 pooled: median 24.1); study eligibility 18-27 kg/m^2",
    renal_function  = "Normal; MDRD eGFR 74.2-143.9 mL/min/1.73 m^2 (Table 1 pooled median 112.35)",
    sex_female_pct  = 22.2,
    race_ethnicity  = c(Asian = 100),
    disease_state   = "Healthy volunteers (no clinically significant medical history, physical examination, vital signs, 12-lead ECG, or clinical laboratory findings)",
    dose_range      = paste(
      "Single oral doses. Study A (NCT03849495): 25 mg empagliflozin",
      "L-proline vs 25 mg conventional empagliflozin. Study B",
      "(NCT03848637): 5 mg/1000 mg empagliflozin L-proline/metformin vs",
      "5 mg/1000 mg empagliflozin/metformin fixed-dose combinations, i.e.",
      "5 mg of the empagliflozin moiety.",
      sep = " "
    ),
    regions         = "Republic of Korea (single centre, Seoul National University Hospital)",
    co_medication   = "Study B only: metformin 1000 mg co-administered as a fixed-dose combination. Jiang 2024 Discussion argues a PK interaction is unlikely because 'metformin and empagliflozin have no common pathways for metabolism or any common transporters', and reports that PK parameters were similar with or without the Study B data.",
    notes           = paste(
      "27 subjects per study (54 total), all healthy Korean adults.",
      "1832 pooled plasma empagliflozin concentrations: 864 from Study A",
      "(432 per formulation) and 968 from Study B (485 conventional, 483",
      "L-proline). Sampling in Study A at 0, 0.33, 0.67, 1, 1.5, 2, 2.5, 3,",
      "4, 6, 8, 10, 12, 24, 34, and 48 h (16 samples); Study B added 3.5 and",
      "5 h (18 samples). LC-MS/MS calibration range 2-1500 ng/mL (Study A)",
      "and 0.5-150 ng/mL (Study B). Empagliflozin L-proline is a cocrystal",
      "that breaks down to empagliflozin in the digestive system and is",
      "absorbed as empagliflozin, so only empagliflozin was assayed. Model",
      "built in Monolix Suite 2021R1 with SAEM; evaluated by GOF plots,",
      "NPDE, VPC (n = 500), and bootstrap (n = 1000). Only an internal",
      "evaluation was performed -- the paper notes the absence of an",
      "independent external dataset as a limitation, along with the",
      "healthy-subject-only and Korean-only composition of the cohort.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. All values are FINAL ESTIMATES from Jiang 2024
    # Table 2 ("Estimates of population pharmacokinetics parameters"),
    # "Fixed effects" block. Every value is corroborated by the bootstrap
    # median in the same row (1000 resamples). Because the source studies
    # are oral-only with no intravenous reference arm, CL, V1, V2, and Q are
    # all APPARENT (/F) quantities; Jiang 2024 does not estimate or report a
    # bioavailability term, so F is implicitly 1 throughout.
    # ---------------------------------------------------------------------

    # Absorption. Jiang 2024 parameterizes the transit chain by BOTH the
    # transit rate constant Ktr and the mean transit time Mtt, and derives
    # the number of transit compartments N from Formula (1), Mtt = (N+1)/Ktr
    # -> N = Ktr*Mtt - 1 = 9.19*0.63 - 1 = 4.79, exactly the value quoted in
    # Results 2.3 ("The number of transit compartments was confirmed to be
    # 4.79 using Formula (1)"). Ktr and Mtt are therefore the two free
    # absorption-shape parameters and N is a derived quantity, computed in
    # model() rather than stored here. (That Ktr and Mtt are genuinely free
    # -- rather than one being a deterministic function of the other -- is
    # established by Table 2 reporting two separately estimated IOV
    # standard deviations for them, 1.17 with RSE 7.43% and 0.42 with RSE
    # 7.15%; two deterministically linked parameters could not support two
    # identifiable IOV variances.)
    lktr <- log(9.19);  label("Transit rate constant Ktr of the absorption chain (1/h)")               # Jiang 2024 Table 2: Ktr = 9.19 1/h (RSE 12.4%; bootstrap 8.94 [6.766-11.728])
    lmtt <- log(0.63);  label("Mean transit time Mtt for absorption (h)")                              # Jiang 2024 Table 2: Mtt = 0.63 h (RSE 4.26%; bootstrap 0.64 [0.526-0.725]). Discussion quotes this as "an Mtt of 37.8 min" (0.63 h * 60 = 37.8 min)
    lka  <- log(0.26);  label("First-order absorption rate constant Ka into the central compartment (1/h)") # Jiang 2024 Table 2: Ka = 0.26 1/h (RSE 3.42%; bootstrap 0.25 [0.23-0.275])

    # Disposition. Jiang 2024 uses the Monolix {Cl, V1, Q, V2}
    # parameterization, in which V1 is the CENTRAL and V2 the PERIPHERAL
    # volume. The paper is explicit and consistent about this: the abstract
    # and Results 2.3 both describe V2 as "the volume of distribution in the
    # peripheral compartment (V2)". V2 therefore maps to the canonical lvp
    # (and carries the weight covariate and the IIV), while V1 maps to lvc.
    lcl  <- log(8.21); label("Apparent clearance CL/F at the reference weight (L/h)")                  # Jiang 2024 Table 2: CL = 8.21 L/h (RSE 1.78%; bootstrap 8.21 [7.945-8.465])
    lvc  <- log(0.6);  label("Apparent central volume of distribution V1/F (L)")                       # Jiang 2024 Table 2: V1 = 0.6 L (RSE 6.67%; bootstrap 0.64 [0-1.889] -- see vignette Errata, this split is near-degenerate)
    lvp  <- log(44.6); label("Apparent peripheral volume of distribution V2/F at the reference weight (L)") # Jiang 2024 Table 2: V2 = 44.6 L (RSE 2.51%; bootstrap 43.93 [41.255-47.171])
    lq   <- log(4.92); label("Apparent intercompartmental clearance Q/F (L/h)")                        # Jiang 2024 Table 2: Q = 4.92 L/h (RSE 5.09%; bootstrap 4.66 [3.979-5.456])

    # Covariate effects. Log-transformed body weight on CL and on the
    # peripheral volume; these are exponents of a power model on WT/70 (see
    # covariateData$WT for the algebraic equivalence and for the
    # centering-weight assumption).
    e_wt_cl <- 0.64; label("Exponent of body weight on apparent clearance (unitless)")                 # Jiang 2024 Table 2: beta_CL_logWT = 0.64 (RSE 24.5%; bootstrap 0.64 [0.366-0.926])
    e_wt_vp <- 0.57; label("Exponent of body weight on apparent peripheral volume (unitless)")         # Jiang 2024 Table 2: beta_V2_logWT = 0.57 (RSE 30.3%; bootstrap 0.54 [0.221-0.851])

    # ---------------------------------------------------------------------
    # Interindividual variability. Table 2 footnote b: "The IIV and IOV are
    # presented as SD (CV)" -- i.e. the first number in each cell is the
    # standard deviation of the log-normal random effect and the
    # parenthesised percentage is the corresponding coefficient of
    # variation. That reading is confirmed exactly, for all seven random
    # effects, by CV = sqrt(exp(SD^2) - 1): e.g. 0.17 -> 17.12%,
    # 0.12 -> 12.04%, 1.17 -> 171.20%. The nlmixr2 omega entries below are
    # therefore VARIANCES, i.e. the published SD squared -- NOT
    # log(CV^2 + 1) applied to the parenthesised percentage (which would
    # double-count the transformation). Jiang 2024 Methods 4.3: "all the
    # individual parameters considered to be log-normally distributed",
    # exponential random effects per Formulas (2) and (3).
    # Note there is no IIV on V1 or Q in the final model.
    # ---------------------------------------------------------------------
    etalka ~ 0.0289;  # Jiang 2024 Table 2 IIV: omega_Ka = 0.17 SD (17.12% CV, RSE 10.8%) -> variance 0.17^2 = 0.0289
    etalcl ~ 0.0144;  # Jiang 2024 Table 2 IIV: omega_CL = 0.12 SD (12.04% CV, RSE 10.9%) -> variance 0.12^2 = 0.0144
    etalvp ~ 0.0169;  # Jiang 2024 Table 2 IIV: omega_V2 = 0.13 SD (13.06% CV, RSE 12.0%) -> variance 0.13^2 = 0.0169

    # ---------------------------------------------------------------------
    # Interoccasion variability across the two crossover periods, on Ktr,
    # Mtt, Ka, and CL (Results 2.3). Table 2 reports ONE IOV standard
    # deviation per parameter, shared across occasions, so occasion 2 is
    # fix()'d equal to occasion 1 (the NONMEM "$OMEGA BLOCK(1) SAME" idiom;
    # nlmixr2 has no SAME shortcut). Variances are the published SD squared,
    # as for the IIV block above.
    # ---------------------------------------------------------------------
    etaiov_ktr_1 ~ 1.3689;        # Jiang 2024 Table 2 IOV: gamma_Ktr = 1.17 SD (171.20% CV, RSE 7.43%) -> variance 1.17^2 = 1.3689 (occasion 1)
    etaiov_ktr_2 ~ fix(1.3689);   # IOV on log-Ktr, occasion 2; variance held equal to occasion 1 per the source's single-SD IOV reporting
    etaiov_mtt_1 ~ 0.1764;        # Jiang 2024 Table 2 IOV: gamma_Mtt = 0.42 SD (43.92% CV, RSE 7.15%) -> variance 0.42^2 = 0.1764 (occasion 1)
    etaiov_mtt_2 ~ fix(0.1764);   # IOV on log-Mtt, occasion 2; variance held equal to occasion 1
    etaiov_ka_1  ~ 0.000841;      # Jiang 2024 Table 2 IOV: gamma_Ka = 0.029 SD (2.90% CV, RSE 34.3%) -> variance 0.029^2 = 0.000841 (occasion 1)
    etaiov_ka_2  ~ fix(0.000841); # IOV on log-Ka, occasion 2; variance held equal to occasion 1
    etaiov_cl_1  ~ 0.001849;      # Jiang 2024 Table 2 IOV: gamma_CL = 0.043 SD (4.30% CV, RSE 16.7%) -> variance 0.043^2 = 0.001849 (occasion 1)
    etaiov_cl_2  ~ fix(0.001849); # IOV on log-CL, occasion 2; variance held equal to occasion 1

    # Residual variability. Jiang 2024 Results 2.2: "A proportional error
    # model was found to best describe the residuals." Table 2 reports the
    # single residual parameter as b = 0.16 -- Monolix's proportional error
    # model is y = f + b*f*eps with eps ~ N(0,1), so b is directly the
    # proportional standard deviation and maps onto nlmixr2's prop().
    propSd <- 0.16; label("Proportional residual error (fraction)")                                    # Jiang 2024 Table 2 Residual variability: b = 0.16 (RSE 2.01%; bootstrap 0.16 [0.141-0.175])
  })

  model({
    # 1. Occasion indicators (binary decomposition of the OCC column).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    # 2. Per-occasion IOV etas. Each affected parameter picks up one of two
    # etas depending on the current occasion; the two variances are equal
    # (see ini()), so this reproduces a single shared IOV magnitude.
    iov_ktr <- oc1 * etaiov_ktr_1 + oc2 * etaiov_ktr_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_ka  <- oc1 * etaiov_ka_1  + oc2 * etaiov_ka_2
    iov_cl  <- oc1 * etaiov_cl_1  + oc2 * etaiov_cl_2

    # 3. Individual structural parameters. The weight terms are the power
    # form of Jiang 2024's log-transformed-weight covariate model:
    # exp(beta * log(WT/70)) == (WT/70)^beta. Reference weight 70 kg (see
    # covariateData$WT -- the source does not report the centering value).
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc)
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q  <- exp(lq)
    ka <- exp(lka + etalka + iov_ka)

    # Absorption-chain shape. N is derived from Ktr and Mtt via Jiang 2024
    # Formula (1), Mtt = (N+1)/Ktr, so N = Ktr*Mtt - 1; at the typical
    # values this is 9.19*0.63 - 1 = 4.79, the value reported in Results
    # 2.3. rxode2's transit() re-derives the internal rate constant as
    # (nn+1)/mtt, which returns exactly ktr, so the parameterization is
    # self-consistent.
    ktr <- exp(lktr + iov_ktr)
    mtt <- exp(lmtt + iov_mtt)
    nn  <- ktr * mtt - 1

    # 4. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 5. ODE system. The analytical Savic transit-compartment chain feeds the
    # depot, which then drains into central by first-order ka. The dose
    # record targets depot and the bolus is suppressed with f(depot) <- 0, so
    # the entire dose enters through the gamma-density input rate below.
    #
    # The gamma density is written out explicitly rather than via rxode2's
    # transit() macro. Both forms are algebraically identical -- the line
    # below is exactly what transit(nn, mtt, 1) expands to, and the vignette
    # asserts the two agree to <1e-8 on Cmax / Tmax / AUC -- but the
    # transit() macro silently returns ZERO input for a model in nlmixr2 UI
    # (ini()/model()) form under rxode2 5.1.3 when combined with the
    # conventional f(depot) <- 0 bolus suppression, which would make the
    # whole dose vanish without any error. podo(depot) is not
    # bioavailability-adjusted, so it still returns the full dose amount
    # under f(depot) <- 0. Jiang 2024 estimates no bioavailability term (CL
    # and the volumes are apparent /F quantities), so no F multiplies podo.
    d/dt(depot)       <- exp(log(podo(depot)) + log(ktr) + nn * log(ktr * tad(depot)) -
                              ktr * tad(depot) - lgamma(nn + 1)) - ka * depot
    d/dt(central)     <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- 0

    # 6. Observation and error. Dose in mg and volumes in L give Cc in mg/L;
    # Jiang 2024 reports plasma concentrations in ng/mL, which is 1000x this
    # value (1 mg/L = 1000 ng/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

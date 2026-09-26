Mehta_2018_fluticasoneFuroate <- function() {
  description <- "Two-compartment population PK model with first-order absorption for inhaled fluticasone furoate in patients with COPD receiving single-inhaler fluticasone furoate/umeclidinium/vilanterol triple therapy, refit on the FULFIL study pooled with the historical fluticasone furoate/vilanterol program"
  reference <- "Mehta R, Pefani E, Beerahee M, Brealey N, Barnacle H, Birk R, Zhu CQ, Lipson DA. Population Pharmacokinetic Analysis of Fluticasone Furoate/Umeclidinium/Vilanterol via a Single Inhaler in Patients with COPD. J Clin Pharmacol. 2018;58(11):1461-1467. doi:10.1002/jcph.1253"
  vignette <- "Mehta_2018_fluticasoneFuroate_umeclidinium_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses are entered in ug and volumes are in L, so `Cc` is in
  # ug/L == ng/mL. Mehta 2018 reports fluticasone furoate concentrations and
  # exposures in pg/mL and pg*h/mL (assay LLQ 10 pg/mL, Sect. 'Pharmacokinetic
  # Assessments'); multiply `Cc` by 1000 to compare against the published
  # values. The source control stream (Supplementary Table 1) carries the same
  # factor as `S2 = V2/1000`; it is NOT reproduced inside the model here so
  # that dose / volume / clearance stay mutually consistent. This matches the
  # sibling extraction `Siederer_2016_fluticasoneFuroate.R`.

  covariateData <- list(
    RACE_ASIAN_EAST_SE = list(
      description = "Composite East Asian / Japanese / South East Asian heritage indicator (source RACE1 = 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White/Caucasian (RACE1 = 1), the model reference category",
      notes = "1 if the subject's race is East Asian, Japanese, or South East Asian heritage; 0 otherwise. Mutually exclusive with RACE_BLACK and RACE_ASIAN_CENTRAL_ARABIC_AMIND_OTH; all three = 0 selects the White/Caucasian reference. Mehta 2018 states that 'a covariate analysis was not planned' and that 'the same covariate relationship was assumed' as in the historical fluticasone furoate model, so the RACE1 structure is carried over from Siederer 2016 but NO re-estimated coefficient is reported for the combined dataset. The coefficient is therefore held at a structural zero here; see the ini() comment and the vignette 'Assumptions and deviations' section. The entire FULFIL PK population was White (Table 3, 74/74 = 100%), i.e. the reference category.",
      source_name = "RACE1 = 2"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator (source RACE1 = 3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White/Caucasian (RACE1 = 1), the model reference category",
      notes = "Carried over from the historical fluticasone furoate model as part of the RACE1 grouping; no coefficient is re-estimated or reported for the combined dataset, so it is held at a structural zero. No FULFIL PK subject was Black (Table 3).",
      source_name = "RACE1 = 3"
    ),
    RACE_ASIAN_CENTRAL_ARABIC_AMIND_OTH = list(
      description = "Composite Central/South Asian, White-Arabic/North African, American Indian/Alaska Native and 'other' heritage indicator (source RACE1 = 4)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White/Caucasian (RACE1 = 1), the model reference category",
      notes = "1 if the subject's race is Central/South Asian, White-Arabic/North African, American Indian/Native Alaskan, or 'other'; 0 otherwise. Mutually exclusive with RACE_ASIAN_EAST_SE and RACE_BLACK. Carried over from the historical model; no coefficient is re-estimated or reported for the combined dataset, so it is held at a structural zero.",
      source_name = "RACE1 = 4"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "fluticasone furoate", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fluticasone furoate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fluticasone furoate", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 74,
    n_studies = 1,
    age_mean = "64 years",
    weight_mean = "81 kg",
    sex_female_pct = 26,
    race_ethnicity = c(`White` = 100),
    disease_state = "symptomatic chronic obstructive pulmonary disease (mean 45% predicted FEV1)",
    dose_range = "fluticasone furoate 100 ug once daily by oral inhalation, as the fluticasone furoate/umeclidinium/vilanterol 100/62.5/25 ug single-inhaler triple combination (Ellipta)",
    regions = "global (162 centers in 15 countries: Russian Federation, Ukraine, Mexico, Germany, Greece, Czech Republic, Romania, Bulgaria, China, Estonia, Hungary, Italy, Poland, Republic of Korea, Slovakia)",
    notes = "Demographics are the FULFIL (CTT116853, NCT02345161) PK population of 74 patients randomized to fluticasone furoate/umeclidinium/vilanterol who provided serial (n = 10) or sparse (n = 64) samples at weeks 12 and 24 (Table 3); mean BMI 28 kg/m2, mean height 171 cm. The PARAMETER ESTIMATES in this file were obtained on a COMBINED dataset that pools these FULFIL data with the historical fluticasone furoate program data used to build the Siederer 2016 model; the size of the historical half is not restated in Mehta 2018, so n_subjects records the FULFIL contribution only. Data below the 10 pg/mL quantification limit were treated as censored using the NONMEM M3 full-likelihood approach (Ahn 2008)."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters. Mehta 2018 Table 1 reports every fluticasone
    # furoate THETA on two scales: the estimated log scale ('Combined Model Ln
    # Estimates') and the back-transformed value ('Model Parameter Estimates
    # With Combined Dataset'). The LOG-scale column is used here because it is
    # the scale NONMEM estimated on (Supplementary Table 1: MU_1 = THETA(1),
    # CL = EXP(MU_1 + ETA(1)), ...) and because it carries three significant
    # figures where the back-transformed column sometimes carries only two
    # (KA 0.053 vs Ln -2.94). The back-transformed value is quoted in each
    # trailing comment for cross-checking.
    #
    # These are the COMBINED-dataset estimates (FULFIL + historical data), not
    # the 'Historical Model' column of the same table -- Table 1's historical
    # column is what the deposited control stream carries as its $THETA
    # initial estimates.
    # -----------------------------------------------------------------------
    lka <- -2.94; label("Apparent first-order absorption rate constant after oral inhalation (log 1/h)") # Table 1 'KA (h-1)' Combined Model Ln Estimate -2.94 (95% CI -3 to -2.88); back-transformed 0.053 1/h (0.049 to 0.056)
    lcl <- 5.43; label("Apparent inhaled clearance CL/F for a White/Caucasian subject (log L/h)") # Table 1 'CL/F (L/h)' Combined Model Ln Estimate 5.43 (5.38 to 5.48); back-transformed 228 L/h (217 to 240)
    lvc <- fixed(0.31); label("Apparent central volume of distribution V2/F (log L)") # Table 1 'V2/F (L)' Combined Model Ln Estimate 0.31 '(Fixed)'; back-transformed 1.36 L '(Fixed)'. Supplementary Table 1 '$THETA 0.31 FIX ; V2'
    lq <- 5.74; label("Apparent intercompartmental clearance Q/F (log L/h)") # Table 1 'Q/F (L/h)' Combined Model Ln Estimate 5.74 (5.54 to 5.94); back-transformed 311 L/h (255 to 380)
    lvp <- 4.66; label("Apparent peripheral volume of distribution V3/F (log L)") # Table 1 'V3 /F (L)' Combined Model Ln Estimate 4.66 (4.45 to 4.87); back-transformed 106 L (86 to 130)

    # -----------------------------------------------------------------------
    # Race effect on CL/F. Supplementary Table 1 retains the four-level RACE1
    # grouping in $PK (the CLRACE1-DEFINITION block, reference RACE1 = 1) and
    # three race THETAs, so the covariate is part of the estimated model.
    #
    # Its MAGNITUDE, however, is not reported for the combined dataset:
    #   * Table 1 lists only the five structural rows -- there is no RACE1 row.
    #   * Methods states 'a covariate analysis was not planned. The same
    #     covariate relationship was assumed'.
    #   * The $THETA values in the deposited stream (-0.211, -0.0602, -0.265)
    #     are INITIAL estimates for this run, seeded from the historical model:
    #     the same block's structural $THETA (5.44, 0.31 FIX, 5.59, 4.71,
    #     -2.95) reproduces Table 1's 'Historical Model Ln Estimates' column
    #     exactly, not its combined column.
    #   * That block is additionally self-inconsistent as printed -- it defines
    #     CLRACE1 = 1 for the reference and then adds it to the log-scale
    #     MU_1, which would give a typical CL/F of exp(5.43 + 1) = 620 L/h
    #     against the 228 L/h Table 1 reports -- and its RACE1 = 3 sign is
    #     opposite to the +0.0602 published in Siederer 2016 Table 1.
    # The coefficients are therefore held at a structural zero rather than
    # transcribed from an initial estimate. For the published historical
    # coefficients use `Siederer_2016_fluticasoneFuroate`, which is fit to the
    # data those coefficients came from. See the vignette 'Assumptions and
    # deviations' section.
    # -----------------------------------------------------------------------
    e_race_asian_east_se_cl <- fixed(0); label("Log-scale effect of East Asian/Japanese/South East Asian heritage on CL/F (unitless)") # Supplementary Table 1 '$PK IF(RACE1.EQ.2) CLRACE1 = ( 1 + THETA(7))'; no combined-dataset estimate reported in Table 1
    e_race_black_cl <- fixed(0); label("Log-scale effect of Black/African American race on CL/F (unitless)") # Supplementary Table 1 '$PK IF(RACE1.EQ.3) CLRACE1 = ( 1 + THETA(8))'; no combined-dataset estimate reported in Table 1
    e_race_asian_central_arabic_amind_oth_cl <- fixed(0); label("Log-scale effect of Central Asian/White-Arabic/American Indian/other heritage on CL/F (unitless)") # Supplementary Table 1 '$PK IF(RACE1.EQ.4) CLRACE1 = ( 1 + THETA(9))'; no combined-dataset estimate reported in Table 1

    # -----------------------------------------------------------------------
    # Inter-individual variability. Supplementary Table 1 places an exponential
    # ETA on all five structural parameters (CL = EXP(MU_1 + ETA(1)), ... ,
    # KA = EXP(MU_5 + ETA(5))), and Results confirms individual MAP Bayes
    # parameter estimates were obtained for all 74 patients -- so the etas are
    # real. Their MAGNITUDES for the combined dataset are not reported: Table 1
    # has no OMEGA rows, and the stream's $OMEGA block (0.367, 3.2, 0.598,
    # 0.46, 0.204) is the initial-estimate half of the same $THETA block shown
    # above to hold the HISTORICAL, not combined, values. They are declared at
    # fixed(0) so the structure is preserved without adopting an initial
    # estimate as a result.
    # -----------------------------------------------------------------------
    etalka ~ fixed(0) # Supplementary Table 1 'KA=EXP(MU_5+ETA(5))'; combined-dataset variance not reported
    etalcl ~ fixed(0) # Supplementary Table 1 'CL=EXP(MU_1+ETA(1))'; combined-dataset variance not reported
    etalvc ~ fixed(0) # Supplementary Table 1 'V2=EXP(MU_2+ETA(2))'; combined-dataset variance not reported
    etalq ~ fixed(0) # Supplementary Table 1 'Q=EXP(MU_3+ETA(3))'; combined-dataset variance not reported
    etalvp ~ fixed(0) # Supplementary Table 1 'V3=EXP(MU_4+ETA(4))'; combined-dataset variance not reported

    # -----------------------------------------------------------------------
    # Residual error. Supplementary Table 1 is explicit that fluticasone
    # furoate uses a pure ADDITIVE (constant) error on the untransformed
    # concentration scale: '$PROB TWO COMPARTMENT FIRST ORDER CONST ERROR',
    # 'SIG=THETA(6) ; RESIDUAL ERROR SD', 'Y=IPRED+SIG*ERR(1)' with
    # '$SIGMA 1 FIX'. This differs from the umeclidinium and vilanterol
    # streams, which use a combined additive-plus-proportional SD.
    #
    # The magnitude is again not reported for the combined dataset (the
    # stream's initial 8.72 pg/mL == 0.00872 ng/mL sits in the historical
    # $THETA block), so it is declared at fixed(0).
    # -----------------------------------------------------------------------
    addSd <- fixed(0); label("Additive residual error (ng/mL)") # Supplementary Table 1 '$ERROR ... Y=IPRED+SIG*ERR(1)' with 'SIG=THETA(6)'; combined-dataset magnitude not reported
  })

  model({
    # Race enters log-additively on CL/F, the form the reference extraction
    # Siederer_2016_fluticasoneFuroate uses for the same RACE1 grouping. All
    # three indicators are 0 for the White/Caucasian reference, which is the
    # whole FULFIL PK population.
    cl <- exp(
      lcl + etalcl +
        e_race_asian_east_se_cl * RACE_ASIAN_EAST_SE +
        e_race_black_cl * RACE_BLACK +
        e_race_asian_central_arabic_amind_oth_cl * RACE_ASIAN_CENTRAL_ARABIC_AMIND_OTH
    )
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}

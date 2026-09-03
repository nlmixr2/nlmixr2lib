Maleki_2024_brepocitinib <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption for oral brepocitinib (PF-06700841, a dual TYK2/JAK1 inhibitor) in healthy participants and patients with alopecia areata, hidradenitis suppurativa, psoriasis, psoriatic arthritis, ulcerative colitis, or vitiligo. Apparent (oral) disposition parameters carry fixed allometric body-weight scaling; clearance additionally depends on baseline aspartate aminotransferase and Asian race, apparent central volume on sex, and the absorption rate constant on ulcerative-colitis status. An absorption lag time applies only to the tablet formulation, and relative bioavailability is 30 percent higher above 100 mg/day. Interindividual variability on CL/F and Vc/F is Box-Cox transformed (Petersson 2009 form) and its magnitude, together with the combined additive-plus-proportional residual error, differs between healthy participants and patients."
  reference <- paste(
    "Maleki F, Clark E, Banfield C, Byon W, Nicholas T.",
    "Population pharmacokinetic modeling of oral brepocitinib in healthy",
    "volunteers and patients with immuno-inflammatory diseases.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(4):551-562.",
    "doi:10.1002/psp4.13100",
    sep = " "
  )
  vignette <- "Maleki_2024_brepocitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline body weight. Allometric scaling referenced to 70 kg,",
        "with exponents fixed at 0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F",
        "(Maleki 2024 Table 3 rows 'Effect of BWT on CL/F and Q/F' and 'Effect of",
        "BWT on Vc/F and Vp/F', both marked Fixed). The reference weight of 70 kg",
        "is a deliberate rounding: the observed median was 81 kg (Table 1), and",
        "the Discussion states 'the typical values of 70 kg and 22 U/L were used",
        "for BWT and BAST'. Observed range 42-204 kg."
      ),
      source_name        = "BWT"
    ),
    AST = list(
      description        = "Baseline aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value entered as a power function referenced to",
        "22 U/L (Maleki 2024 Table 3 footnote e). As with body weight, 22 U/L is",
        "a rounded reference rather than the observed median of 19 U/L (Table 1);",
        "the Discussion states the rounded values were used deliberately.",
        "Observed range 8-133 U/L; the paper notes the normal range is 8-33 U/L."
      ),
      source_name        = "BAST"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; the pooled White and Other groups)",
      notes              = paste(
        "Maleki 2024 Table 3 footnote e: RACE is 1 for the Asian population and 0",
        "for others. Asians were 7% of the analysis population (52 of 775);",
        "White 85% and Other 8% form the reference group. The effect is on CL/F only."
      ),
      source_name        = "RACE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) for this model",
      notes              = paste(
        "The source column is a male indicator (Maleki 2024 Table 3 footnote f:",
        "'SEX is 1 for male population'), so SEXF = 1 - SEX. Females are the",
        "paper's reference category -- the typical Vc/F of 88.52 L is a female",
        "value -- so the effect is applied as exp(e_sexf_vc * (1 - SEXF)) to",
        "preserve the published female reference and the published positive",
        "coefficient. Reported as 13% higher Vc/F in males. Note that the paper's",
        "reference category was 'selected alphabetically' (Discussion), not by",
        "group size: males were 59% of the cohort."
      ),
      source_name        = "SEX"
    ),
    DIS_UC = list(
      description        = "Ulcerative colitis disease-state indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participants and the AA / HS / PsA / PsO / vitiligo cohorts pooled)",
      notes              = paste(
        "Maleki 2024 Table 3 footnote d: TA is 1 for patients with UC and 0 for",
        "everyone else. Patients with UC were 19% of the analysis population",
        "(145 of 775). The effect is additive on Ka (Ka is a normally distributed",
        "parameter in Table 2), not multiplicative, and is the only structural",
        "effect of disease state on the typical parameters."
      ),
      source_name        = "TA"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient; the pooled AA / HS / PsA / PsO / UC / vitiligo cohorts)",
      notes              = paste(
        "Healthy participants were 12% of the analysis population (92 of 775).",
        "This covariate does not change any typical structural parameter. It",
        "gates (i) the magnitude of the interindividual variability on CL/F and",
        "Vc/F -- Maleki 2024 Table 3 reports the healthy-participant omegas with a",
        "separate multiplicative 'Effect of patients' row for each -- and (ii)",
        "which pair of residual-error magnitudes applies, since Table 3 reports",
        "the additive and proportional residual SDs separately for healthy",
        "participants and for patients. The source encoded the covariate as a",
        "patient indicator ('Health status (healthy vs patient)', Table 2), so",
        "both effects are applied on (1 - DIS_HEALTHY)."
      ),
      source_name        = "Health status"
    ),
    FORM_TABLET = list(
      description        = "Tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral suspension)",
      notes              = paste(
        "Maleki 2024 Table 3 footnote b: 'Tlag = 0 for suspension formulation'.",
        "The lag time of 0.26 h therefore applies only when FORM_TABLET = 1.",
        "Tablet was 92% of the analysis population; the non-tablet comparator is",
        "the oral suspension used in the phase I relative-bioavailability arm",
        "(NCT02310750)."
      ),
      source_name        = "Drug formulation"
    ),
    DOSE_HIGH = list(
      description        = "Daily-dose-above-100-mg indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (administered daily dose of 100 mg or less)",
      notes              = paste(
        "Maleki 2024 Table 3 footnote c: 'D100 is 1 when the administered dose",
        ">100 and 0 for other dose levels', with dose in mg/day. In this analysis",
        "the >100 mg/day cohorts were the 175 and 200 mg phase I arms, 5% of the",
        "population (41 of 775); the 100 mg cohort (4%) sits in the reference",
        "group. The effect is on relative bioavailability, giving 30% higher F",
        "and hence 1.3-fold higher dose-normalized AUCtau and Cmax above 100 mg",
        "(Discussion). Time-fixed per subject in this design."
      ),
      source_name        = "D100"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Reported as a baseline demographic in Maleki 2024 Table 1 (median 43,",
        "range 18-75 years) and covered by the univariate covariate screen, but",
        "not retained in the final covariate model (Table 2 lists only the",
        "covariates carried forward to the multivariate analysis, and age is not",
        "among them). Documentation only; not referenced in model()."
      )
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Reported as a baseline demographic in Maleki 2024 Table 1 (median 4.5,",
        "range 3.2-5.6 g/dL) but not retained in the final covariate model.",
        "Documentation only; not referenced in model()."
      )
    ),
    CRCL = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Reported as a baseline demographic in Maleki 2024 Table 1 (median 119,",
        "range 41.4-357 mL/min; Table 1 footnote a cites Cockcroft and Gault) but",
        "not retained in the final covariate model. Not BSA-normalized in the",
        "source. Documentation only; not referenced in model()."
      )
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "brepocitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "brepocitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "brepocitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 775,
    n_studies      = 9,
    n_observations = 8552,
    age_range      = "18-75 years",
    age_median     = "43 years",
    weight_range   = "42-204 kg",
    weight_median  = "81 kg",
    sex_female_pct = 41,
    race_ethnicity = c(White = 85, Asian = 7, Other = 8),
    disease_state  = paste(
      "Healthy participants (12%) pooled with patients with alopecia areata (8%),",
      "hidradenitis suppurativa (5%), psoriatic arthritis (24%), plaque psoriasis",
      "(26%), ulcerative colitis (19%), and non-segmental vitiligo (7%)"
    ),
    dose_range     = "1-200 mg oral, single dose or once daily (one 50 mg twice-daily phase I arm); tablet (92%) or suspension",
    regions        = "Not reported by region; includes a dedicated phase I study in healthy Japanese participants (NCT03236493)",
    notes          = paste(
      "Baseline demographics from Maleki 2024 Table 1; study-level detail from",
      "Table S1 of the supplement. Three phase I studies (NCT02310750,",
      "NCT03236493, NCT03656952) and six phase II studies (NCT02969018,",
      "NCT02974868, NCT03963401, NCT02958865, NCT03715829, NCT04092452).",
      "8552 plasma concentrations from 775 individuals were used for the",
      "analysis reported in the Abstract; the Results section reports 9532 data",
      "records from the same 775 individuals, the difference being dosing and",
      "other non-observation records. LLOQ was 0.2 ng/mL and fewer than 10% of",
      "samples were below it, so no BLQ imputation was applied. Placebo records",
      "and screening observations were excluded. Estimation was by SAEM in",
      "Monolix 2021R1, not by NONMEM."
    )
  )

  ini({
    # ==================================================================
    # Structural parameters (apparent / oral) for the paper's reference
    # individual: 70 kg, non-Asian, female, non-UC, baseline AST 22 U/L.
    #
    # SOURCE-CONFLICT NOTE (see the vignette's Assumptions and deviations
    # section for the full reconciliation): Maleki 2024 Table 3 reports an
    # "Estimate" column whose values disagree with the covariate equations
    # printed in the same table's own footnotes c, d, e, and f. Where a
    # footnote equation exists it is used here, because the footnote values
    # -- and not the Estimate column -- are the ones reproduced in the
    # Abstract and Discussion (CL/F 17.5 L/h, Vc/F 88.5 L, 30% higher F
    # above 100 mg, ~10% lower CL/F in Asians, 13% higher Vc/F in males).
    # Footnote d is the single exception: its coefficient of -0.61 on Ka
    # falls outside the bootstrap 95% CI (-1.58, -0.68) printed on the same
    # row and implies 23% slower absorption against the paper's stated 46%,
    # so the Table 3 Estimate is used for that row instead. Parameters with
    # no footnote equation take the Table 3 Estimate column, their only source.
    # ==================================================================
    lka     <- log(2.66)     ; label("Absorption rate constant for non-UC subjects (Ka, 1/h)")             # Table 3 row 'Ka (1/h)', Estimate column
    lcl     <- log(17.49)    ; label("Apparent oral clearance (CL/F, L/h)")                                # Table 3 footnote e (CL/F = 17.49 L/h); Abstract and Discussion give 17.5
    lvc     <- log(88.52)    ; label("Apparent central volume of distribution (Vc/F, L)")                  # Table 3 footnote f (Vc/F = 88.52 L); Abstract gives 88.5
    lq      <- log(0.55)     ; label("Apparent inter-compartmental clearance (Q/F, L/h)")                  # Table 3 row 'Q/F (L/h)', Estimate column
    lvp     <- log(232.5)    ; label("Apparent peripheral volume of distribution (Vp/F, L)")               # Table 3 row 'Vp/F (L)', Estimate column; poorly identified, bootstrap 95% CI (12, 455)
    ltlag   <- log(0.26)     ; label("Absorption lag time of the tablet formulation (h)")                  # Table 3 row 'Tlag of tablet formulation (h)', Estimate column; footnote b: Tlag = 0 for suspension
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability at daily doses of 100 mg or less (unitless)") # Table 3 row 'Frel' = 1 with no bootstrap entry, i.e. a structural anchor

    # ==================================================================
    # Covariate effects
    # ==================================================================
    e_wt_cl_q          <- fixed(0.75)  ; label("Allometric body-weight exponent shared by CL/F and Q/F (unitless)")   # Table 3 row 'Effect of BWT on CL/F and Q/F' = 0.75, marked Fixed
    e_wt_vc_vp         <- fixed(1)     ; label("Allometric body-weight exponent shared by Vc/F and Vp/F (unitless)")  # Table 3 row 'Effect of BWT on Vc/F and Vp/F' = 1, marked Fixed
    e_ast_cl           <- -0.17        ; label("Power exponent of baseline AST on CL/F, referenced to 22 U/L (unitless)") # Table 3 footnote e: (BAST / 22)^-0.17
    e_race_asian_cl    <- -0.11        ; label("Log-scale effect of Asian race on CL/F (unitless)")                   # Table 3 footnote e: exp(-0.11 * RACE); Discussion reports ~10% slower clearance, exp(-0.11) = 0.896
    e_sexf_vc          <- 0.12         ; label("Log-scale effect of male sex on Vc/F, applied on (1 - SEXF) (unitless)") # Table 3 footnote f: exp(0.12 * SEX) with SEX = 1 for males; Discussion reports 13% higher Vc/F, exp(0.12) = 1.128
    e_dis_uc_ka        <- -1.16        ; label("Additive shift of ulcerative colitis on Ka (1/h)")                    # Table 3 row 'Effect of UC patients on Ka' = 1.16 with bootstrap 95% CI (-1.58, -0.68); sign taken from the CI. Footnote d's -0.61 is outside that CI and is not used
    e_dose_high_fdepot <- 0.3          ; label("Fractional increase in relative bioavailability above 100 mg/day (unitless)") # Table 3 footnote c: F = 1 * (1 + 0.3 * D100); Discussion reports 30% higher bioavailability

    # ==================================================================
    # Interindividual variability
    #
    # Maleki 2024 Methods, "Stochastic model development": IIV enters as
    #   P_i = theta_P * exp(eta_i_transformed),
    #   eta_i_transformed = ((exp(eta_i))^theta_s - 1) / theta_s
    # i.e. the Box-Cox transformation of the random effect described by
    # Petersson et al. (2009) Pharm Res 26:2174-2185. theta_s is the shape
    # parameter, tabulated below as boxcox_cl and boxcox_vc; theta_s -> 0
    # recovers the usual log-normal IIV. The transformation is applied
    # explicitly inside model() because rxode2 has no native Box-Cox eta
    # distribution.
    #
    # The omegas below are the healthy-participant values; patients carry a
    # multiplicative uplift on the eta magnitude (see e_dis_healthy_eta*).
    # ==================================================================
    etalcl + etalvc ~ c(
      # var(etalcl)                              # Table 3 row 'omega2 CL/F (SD)' = 0.42
      0.42 * 0.42,
      # cov(etalcl, etalvc); rho = 0.71          # Table 3 row 'rho CL/F, Vc/F' = 0.71
      0.71 * 0.42 * 0.19,
      # var(etalvc)                              # Table 3 row 'omega2 Vc/F (SD)' = 0.19
      0.19 * 0.19
    )

    boxcox_cl             <- 0.26 ; label("Box-Cox shape parameter of the CL/F random effect (unitless)")   # Table 3 row 'Box-Cox transformation shape parameter on omega2 CL/F' = 0.26
    boxcox_vc             <- 5.86 ; label("Box-Cox shape parameter of the Vc/F random effect (unitless)")   # Table 3 row 'Box-Cox transformation shape parameter on omega2 Vc/F' = 5.86
    e_dis_healthy_etalcl  <- 0.44 ; label("Fractional uplift of the CL/F IIV magnitude in patients (unitless)") # Table 3 row 'Effect of patients on omega2 CL/F' = 0.44; applied as (1 + 0.44) on (1 - DIS_HEALTHY)
    e_dis_healthy_etalvc  <- 0.29 ; label("Fractional uplift of the Vc/F IIV magnitude in patients (unitless)") # Table 3 row 'Effect of patients on omega2 Vc/F' = 0.29; applied as (1 + 0.29) on (1 - DIS_HEALTHY)

    # ==================================================================
    # Residual unexplained variability
    #
    # Maleki 2024 Methods prints the combined error as
    #   DV_ij = IPRED_ij + sqrt(sigma_add^2 + (sigma_pro * IPRED_ij)^2) * eps_ij
    # which is exactly nlmixr2's default combined2 form for
    # `add(addSd) + prop(propSd)`. The additive SDs below are therefore in
    # ng/mL, the units of the analysis dataset (LLOQ 0.2 ng/mL). Table 3
    # reports a separate pair for healthy participants and for patients.
    # ==================================================================
    addSdHv   <- 0.89 ; label("Additive residual SD in healthy participants (ng/mL)")       # Table 3 row 'Additive residual error of HV (SD)' = 0.89
    propSdHv  <- 0.1  ; label("Proportional residual SD in healthy participants (fraction)") # Table 3 row 'Proportional residual error of HV (SD)' = 0.1
    addSdPat  <- 0.45 ; label("Additive residual SD in patients (ng/mL)")                     # Table 3 row 'Additive residual error of patients (SD)' = 0.45
    propSdPat <- 0.05 ; label("Proportional residual SD in patients (fraction)")              # Table 3 row 'Proportional residual error of patients (SD)' = 0.05
  })

  model({
    # ---- 1. Covariate-derived multipliers ----
    # Allometry referenced to 70 kg, exponents fixed by the source.
    allomCl <- (WT / 70.0)^e_wt_cl_q
    allomVc <- (WT / 70.0)^e_wt_vc_vp

    # ---- 2. Box-Cox transformed random effects ----
    # The eta magnitude itself carries the health-status covariate: Ka and
    # the omegas are normally distributed parameters in Maleki 2024 Table 2,
    # for which the paper's categorical-covariate form is
    # P = theta_P * (1 + theta_COV) when the covariate indicator is 1. The
    # tabulated omegas are the healthy-participant values, so the uplift is
    # applied on the patient indicator (1 - DIS_HEALTHY). Scaling a
    # mean-zero normal deviate multiplies its SD by the same factor.
    etalclScaled <- etalcl * (1.0 + e_dis_healthy_etalcl * (1.0 - DIS_HEALTHY))
    etalvcScaled <- etalvc * (1.0 + e_dis_healthy_etalvc * (1.0 - DIS_HEALTHY))

    # Petersson 2009 Box-Cox transformation of the random effect.
    etalclBc <- (exp(boxcox_cl * etalclScaled) - 1.0) / boxcox_cl
    etalvcBc <- (exp(boxcox_vc * etalvcScaled) - 1.0) / boxcox_vc

    # ---- 3. Individual parameters ----
    cl <- exp(lcl) * allomCl *
      (AST / 22.0)^e_ast_cl *
      exp(e_race_asian_cl * RACE_ASIAN) *
      exp(etalclBc)
    vc <- exp(lvc) * allomVc *
      exp(e_sexf_vc * (1.0 - SEXF)) *
      exp(etalvcBc)
    q  <- exp(lq)  * allomCl
    vp <- exp(lvp) * allomVc

    # Ka carries an ADDITIVE ulcerative-colitis shift, matching the paper's
    # covariate form for normally distributed parameters and the sign of the
    # bootstrap CI on that row.
    ka <- exp(lka) + e_dis_uc_ka * DIS_UC

    # Relative bioavailability and the tablet-only absorption lag.
    fdepot <- exp(lfdepot) * (1.0 + e_dose_high_fdepot * DOSE_HIGH)
    tlag   <- exp(ltlag) * FORM_TABLET

    # ---- 4. Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 5. ODE system ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- 6. Dose input ----
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ---- 7. Observation and residual error ----
    # Amounts are mg and volumes L, so central / vc is mg/L; the factor of
    # 1000 converts to ng/mL, the concentration unit of the analysis dataset
    # (LLOQ 0.2 ng/mL) in which the residual-error SDs are reported.
    Cc <- 1000.0 * central / vc

    addSdCohort  <- addSdHv  * DIS_HEALTHY + addSdPat  * (1.0 - DIS_HEALTHY)
    propSdCohort <- propSdHv * DIS_HEALTHY + propSdPat * (1.0 - DIS_HEALTHY)

    Cc ~ prop(propSdCohort) + add(addSdCohort)
  })
}

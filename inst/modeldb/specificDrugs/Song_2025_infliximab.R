Song_2025_infliximab <- function() {
  description <- "Two-compartment population PK model of intravenous and subcutaneous infliximab in Korean adults with inflammatory bowel disease on maintenance therapy, with body-mass-index, serum-albumin, C-reactive-protein and quantitative anti-drug-antibody-concentration effects on clearance and on subcutaneous bioavailability (Song 2025)"
  reference <- "Song JH, Hong SN, Kim MG, Kim M, Kim SK, Kim ER, Chang DK, Kim YH. Population pharmacokinetic model for the use of intravenous or subcutaneous infliximab in patients with inflammatory bowel disease: real-world data from a prospective cohort study. Gut Liver. 2025;19(3):376-387. doi:10.5009/gnl240503"
  vignette <- "Song_2025_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Song 2025 Methods section 6 (two
  # compartments with first-order absorption from the subcutaneous site and
  # first-order elimination) and section 5 (serum infliximab assay).
  compartmentData <- list(
    depot       = list(analyte = "infliximab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    BMI = list(
      description        = "Body mass index (time-varying; recorded at every outpatient visit)",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on BOTH clearance and subcutaneous bioavailability, normalised to 22.5 kg/m^2: CL includes (BMI/22.5)^0.360 and F includes (BMI/22.5)^-0.832 (Song 2025 p. 380 final-model equations; Table 2 theta BMI-CL and theta BMI-CL). Song 2025 evaluated body weight and BMI head-to-head on CL and reports BMI as the better descriptor (Discussion: 'BMI provides a superior explanation of infliximab PK'), so WT is not in the final model. The 22.5 normalising value is the cohort mean per Methods section 6 ('Continuous covariates were normalized to the mean value of the data') and coincides with the BMI implied by the Table 1 medians (65.0 kg / 1.700 m^2 = 22.5 kg/m^2); BMI itself is not tabulated in Table 1.",
      source_name        = "BMI"
    ),
    ALB = list(
      description        = "Serum albumin (time-varying; measured at every outpatient visit)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance only: (ALB_gdL/4.4)^-0.372 (Song 2025 p. 380 final-model equation; Table 2 theta ALB-CL = -0.372). Song 2025 reports albumin in US-convention g/dL (Table 1 median 4.4 g/dL, IQR 4.2-4.6) and calibrated the exponent on that scale, so model() applies the inline conversion alb_gdL <- ALB * 0.1 required by the ALB register entry's Units note; the canonical column stays SI g/L, i.e. the reference 4.4 g/dL is 44 g/L. The 4.4 g/dL normalising value is the cohort mean per Methods section 6 and happens to equal the Table 1 median.",
      source_name        = "ALB"
    ),
    CRP = list(
      description        = "C-reactive protein, standard (non-high-sensitivity) assay (time-varying; measured at every outpatient visit)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance only: (CRP_mgdL/0.18)^0.022 (Song 2025 p. 380 final-model equation; Table 2 theta CRP-CL = 0.022). Song 2025 reports CRP in US-convention mg/dL (Table 1 median 0.06 mg/dL, IQR 0.06-0.15) and calibrated the exponent on that scale, so model() applies the inline conversion crp_mgdL <- CRP * 0.1; the canonical column stays SI mg/L, i.e. the reference 0.18 mg/dL is 1.8 mg/L. The 0.18 mg/dL normalising value is the cohort MEAN per Methods section 6 and is materially above the Table 1 median of 0.06 mg/dL, as expected for a right-skewed acute-phase reactant -- do not substitute the median. Standard assay: Song 2025 Methods section 4 lists CRP among routine laboratory findings with no high-sensitivity qualifier.",
      source_name        = "CRP"
    ),
    CONC_ADA_NGML = list(
      description        = "Total (free plus drug-bound) anti-infliximab antibody concentration measured by a drug-tolerant precipitation / acid-dissociation immunoassay (time-varying; drawn with each pre-dose PK sample)",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on BOTH clearance and subcutaneous bioavailability, normalised to 10 ng/mL: CL includes (CONC_ADA_NGML/10)^0.022 and F includes (CONC_ADA_NGML/10)^-0.213 (Song 2025 p. 380 final-model equations; Table 2 theta ADA-CL = 0.022, theta ADA-F = -0.213). This column is a calibrated MASS CONCENTRATION, not a dilution titer: every subject carries a strictly positive measured value, ADA-negatives included (Song 2025 Results: median 8.0 ng/mL, IQR 7.1-9.0 for ADA-negative vs median 12.6 ng/mL, IQR 11.0-19.0 for ADA-positive; overall Table 1 median 9.3 ng/mL, IQR 7.6-11.6). There is therefore NO zero- or one-encoding for ADA-negative subjects, and the power form (CONC_ADA_NGML/10)^theta is undefined at 0 -- which is why this is not ADA_TITER. Positivity, where the paper reports it, is a derived dichotomy at the assay cutoff of 10 ng/mL (Methods section 5), the same value the authors used as the normalising constant. Song 2025 treats the continuous encoding as the paper's central methodological contribution (Discussion: 'the distinctive feature of our model was that ADA levels were chosen as covariates rather than the presence of ADA ... while binary classification is intuitive, it can lead to information loss'), so do NOT dichotomise onto ADA_POS. Assay is drug-tolerant, so the measurement is not confounded by circulating infliximab.",
      source_name        = "ADA"
    )
  )

  # Covariates that Song 2025 screened in the stepwise covariate search but did
  # NOT retain in the final model (Methods section 6: covariates evaluated for
  # CL were age, sex, BW, height, BMI, diagnosis, ADA, ALB, CRP and WBC;
  # covariates evaluated for F were BW, BMI and ADA). Only BMI, ALB, CRP and
  # ADA survived forward selection / backward elimination, and the paper
  # reports no point estimate for any rejected covariate (the per-step
  # objective-function changes are in Supplementary Table 1, which is not on
  # disk). Recorded here so the provenance of the covariate screen is not lost;
  # documentation only, never referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 36 years (IQR 30-47) per Song 2025 Table 1."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL; not retained. Song 2025 Table 1 reports 130/181 (71.8%) male, i.e. 28.2% female."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on CL and on F; not retained in either. Song 2025 evaluated body weight and BMI head-to-head and reports BMI as the better descriptor of infliximab PK (Discussion), so the retained covariate is BMI. Note weight still enters the DOSING regimen (5 or 10 mg/kg IV) even though it is absent from the parameter model. Cohort median 65.0 kg (IQR 55.0-73.6) per Table 1."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on CL; not retained (its information enters through BMI). Cohort median 170.0 cm (IQR 162.8-174.2) per Song 2025 Table 1."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 5.7 x 10^9/L (IQR 4.7-6.8) per Song 2025 Table 1."
    ),
    DIS_CD = list(
      description = "Crohn's disease indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Diagnosis (CD vs UC) was screened on CL and not retained; Song 2025 Discussion reads this as 'no significant PK differences between the diseases'. 149/181 (82.3%) CD per Table 1."
    ),
    DIS_UC = list(
      description = "Ulcerative colitis indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Complementary level of the screened diagnosis covariate; not retained. 32/181 (17.7%) UC per Song 2025 Table 1."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 181L,
    n_studies        = 1L,
    age_range        = "Age at diagnosis <16 years in 7.7%, 17-40 years in 79.0%, >40 years in 13.3% (adults at the time of sampling)",
    age_median       = "36 years (IQR 30-47)",
    weight_range     = "IQR 55.0-73.6 kg (full range not reported)",
    weight_median    = "65.0 kg",
    height_median    = "170.0 cm (IQR 162.8-174.2)",
    bmi_reference    = "22.5 kg/m^2 (cohort mean used as the covariate normalising value; BMI is not tabulated in Table 1)",
    sex_female_pct   = 28.2,
    race_ethnicity   = "All Korean (single-centre Korean cohort); the authors list the resulting limited generalisability to other ethnicities as a study limitation.",
    disease_state    = "Inflammatory bowel disease on maintenance infliximab: Crohn's disease 149/181 (82.3%), ulcerative colitis 32/181 (17.7%). Remnant disease burden at enrolment: transmural healing (global sMaRIA = 0) in 16/106 CD, complete endoscopic remission (SES-CD = 0) in 8/43 CD, endoscopic remission (Mayo endoscopic score = 0) in 16/32 UC.",
    dose_range       = "Maintenance only. IV: 5 mg/kg every 8 weeks (139/181 patients at enrolment) or 10 mg/kg every 8 weeks after dose escalation (42/181 at enrolment). SC: Remsima SC 120 mg fixed dose by prefilled pen every 2 weeks. Induction (5 mg/kg IV at weeks 0, 2 and 6) preceded enrolment and is NOT covered by this model -- the authors list inapplicability to induction as a study limitation.",
    regions          = "Single centre, Samsung Medical Center, Seoul, Korea; enrolment February 2020 to December 2022.",
    concomitant_immunomodulator = "95/181 (52.5%) at enrolment.",
    infliximab_product = "Remicade 112/181 (61.9%), Remsima 58/181 (32.0%), Remaloce 11/181 (6.1%) at enrolment; SC product is Remsima SC.",
    prior_surgery    = "53/181 (29.3%) with previous intestinal surgery.",
    smoking          = "25/181 (13.8%) current smokers.",
    ada_incidence    = "ADA positive (>10 ng/mL) at least once in 171 patients, covering 861 of the 2,132 measurements.",
    n_observations   = 2132L,
    sampling         = "Sparse, pre-dose trough samples only, drawn under a proactive therapeutic-drug-monitoring protocol: IV troughs at 8 +/- 2 weeks after the last dose, SC troughs at 14 +/- 4 days after the last dose. Concentrations below the 3 ng/mL limit of quantitation were substituted with 1.5 ng/mL.",
    notes            = "Baseline demographics from Song 2025 Table 1 (n = 181, 2,132 samples). 212 patients were screened; 31 excluded (1 withdrew consent, 5 drug holiday >4 months, 25 switched to SC maintenance immediately after IV induction). Reference covariate values for the typical patient (the cohort MEANS per Methods section 6, which differ from the Table 1 medians for CRP): ALB = 4.4 g/dL (44 g/L), CRP = 0.18 mg/dL (1.8 mg/L), ADA = 10 ng/mL, BMI = 22.5 kg/m^2. Because only sparse trough data were available, Vp and Q were fixed to literature values taken from Hanzel 2021 (see modellib('Hanzel_2021_infliximab'))."
  )

  ini({
    # Structural parameters. These are the typical values AT the reference
    # covariate vector (ALB 4.4 g/dL, CRP 0.18 mg/dL, ADA 10 ng/mL,
    # BMI 22.5 kg/m^2), because every covariate term in the published
    # equations is normalised to those values and so evaluates to 1 there.
    lcl     <- log(0.248); label("Clearance at the reference covariate vector (CL, L/day)")                      # Song 2025 Table 2: CL = 0.248 L/day (RSE 20%); bootstrap 0.249 (0.173-0.344)
    lvc     <- log(1.87);  label("Central volume of distribution (Vc, L)")                                       # Song 2025 Table 2: Vc = 1.87 L (RSE 61%); bootstrap 1.88 (0.23-4.36); also stated in the Results text
    lka     <- log(0.083); label("First-order subcutaneous absorption rate constant (ka, 1/day)")                # Song 2025 Table 2: ka = 0.083 /day (RSE 26%); bootstrap 0.083 (0.051-0.151); also stated in the Results text
    lfdepot <- log(0.667); label("Subcutaneous bioavailability at the reference covariate vector (F, fraction)")  # Song 2025 Table 2: F = 0.667 (RSE 19%); bootstrap 0.677 (0.477-0.921)

    # Vp and Q were NOT estimated. Song 2025 Methods section 6: "only sparsely
    # measured concentrations were utilized, necessitating the fixation of Vp
    # and Q to literature-based values (1.96 L and 0.599 L/day, respectively)",
    # citing Hanzel 2021 (reference 29 of Song 2025), which is itself in this
    # library as Hanzel_2021_infliximab (Vp = 1.93 L, Q = 0.598 L/day at its
    # 70 kg reference patient -- the same values to rounding). Table 2 prints
    # both as "fixed", so they are wrapped in fixed() here; the log() goes
    # INSIDE fixed().
    lvp     <- fixed(log(1.96));  label("Peripheral volume of distribution (Vp, L), from Hanzel 2021")     # Song 2025 Table 2: Vp = 1.96 fixed; Methods section 6
    lq      <- fixed(log(0.599)); label("Inter-compartmental clearance (Q, L/day), from Hanzel 2021")      # Song 2025 Table 2: Q = 0.599 fixed; Methods section 6

    # Covariate effects. Song 2025 Methods section 6: "Continuous covariates
    # were normalized to the mean value of the data, while both continuous and
    # discrete covariates were incorporated into the base model using power
    # functions." The two published final-model equations (p. 380, column 1;
    # rendered as images in the PDF and recovered with `pdftotext -layout`) are
    #   CL = 0.248 * (ALB/4.4)^-0.372 * (CRP/0.18)^0.022
    #               * (ADA/10)^0.022  * (BMI/22.5)^0.360
    #   F  = 0.667 * (ADA/10)^-0.213  * (BMI/22.5)^-0.832
    # ADA enters as the quantitative ng/mL concentration, NOT as a positivity
    # flag -- see covariateData[[CONC_ADA_NGML]]$notes.
    e_alb_cl      <- -0.372; label("Power exponent of serum albumin on CL ((ALB_gdL/4.4)^e_alb_cl)")                        # Song 2025 Table 2: theta ALB-CL = -0.372 (RSE 30%); bootstrap -0.383 (-0.592 to -0.154)
    e_crp_cl      <-  0.022; label("Power exponent of C-reactive protein on CL ((CRP_mgdL/0.18)^e_crp_cl)")                 # Song 2025 Table 2: theta CRP-CL = 0.022 (RSE 33%); bootstrap 0.023 (0.010-0.036)
    e_conc_ada_cl <-  0.022; label("Power exponent of ADA concentration on CL ((CONC_ADA_NGML/10)^e_conc_ada_cl)")          # Song 2025 Table 2: theta ADA-CL = 0.022 (RSE 42%); bootstrap 0.021 (0.001-0.042)
    e_bmi_cl      <-  0.360; label("Power exponent of body mass index on CL ((BMI/22.5)^e_bmi_cl)")                         # Song 2025 Table 2: theta BMI-CL = 0.360 (RSE 31%); bootstrap 0.354 (0.150-0.627)
    e_conc_ada_f  <- -0.213; label("Power exponent of ADA concentration on SC bioavailability ((CONC_ADA_NGML/10)^e_conc_ada_f)") # Song 2025 Table 2: theta ADA-F = -0.213 (RSE 19%); bootstrap -0.209 (-0.298 to -0.057)
    e_bmi_f       <- -0.832; label("Power exponent of body mass index on SC bioavailability ((BMI/22.5)^e_bmi_f)")          # Song 2025 Table 2: theta BMI-F = -0.832 (RSE 33%); bootstrap -0.838 (-1.450 to -0.300)

    # Inter-individual variability. Song 2025 Methods section 6: "The
    # inter-individual variability was modeled using exponential error model.
    # Inter-individual variability was evaluated for all PK parameters, but it
    # was only included for CL and F." -- so there is deliberately NO IIV on
    # Vc or ka, and none is invented here.
    #
    # Table 2 reports the two diagonals as %CV and the off-diagonal as a
    # covariance on the same (log-scale variance) footing:
    #   omega CL = 20.3 %CV  ->  var = 0.203^2 = 0.041209
    #   omega F  = 21.6 %CV  ->  var = 0.216^2 = 0.046656
    #   Cov CL-F = 0.017     ->  correlation = 0.017 / sqrt(0.041209 * 0.046656) = 0.388
    # The standard NONMEM reporting convention CV% = 100 * sqrt(omega^2) is
    # used. The log-normal-exact alternative omega^2 = log(1 + CV^2) would give
    # 0.040404 / 0.045648 and a correlation of 0.396 -- under a 2% difference
    # in the variances, so the two readings are not discriminable here; see the
    # vignette Errata.
    etalcl + etalfdepot ~ c(0.041209,
                            0.017,    0.046656)  # Song 2025 Table 2: omega CL 20.3 %CV, omega F 21.6 %CV, Cov CL-F 0.017

    # Residual error. Song 2025 Methods section 6: "The residual variability
    # was modeled using combined proportional and additive error model."
    # Table 2 labels the proportional term with the sigma symbol and gives the
    # additive term in mg/L, i.e. both are on the SD scale (a variance reading
    # of 0.333 would imply a 58% proportional SD, implausible for these data
    # and inconsistent with the explicit sigma label).
    addSd  <- 0.561; label("Additive residual error (mg/L)")            # Song 2025 Table 2: additive error = 0.561 mg/L (RSE 12%); bootstrap 0.571 (0.438-0.732)
    propSd <- 0.333; label("Proportional residual error (fraction)")    # Song 2025 Table 2: proportional error sigma = 0.333 (RSE 6%); bootstrap 0.328 (0.290-0.368)
  })
  model({
    # Unit conversions. The canonical ALB and CRP data columns are SI (g/L and
    # mg/L); Song 2025 reports and calibrated against US-convention g/dL and
    # mg/dL, so convert here to keep the published exponents aligned with the
    # scale they were fitted on. Required by the ALB register entry's Units
    # note; applied to CRP for the same reason.
    alb_gdL  <- ALB * 0.1   # SI g/L  -> US-convention g/dL  (reference 4.4 g/dL  = 44 g/L)
    crp_mgdL <- CRP * 0.1   # SI mg/L -> US-convention mg/dL (reference 0.18 mg/dL = 1.8 mg/L)

    # Individual PK parameters, exactly as published (Song 2025 p. 380).
    # Every covariate term is normalised to the cohort mean, so at the
    # reference vector cl reduces to exp(lcl) and fdepot to exp(lfdepot).
    cl <- exp(lcl + etalcl) *
      (alb_gdL / 4.4)^e_alb_cl *
      (crp_mgdL / 0.18)^e_crp_cl *
      (CONC_ADA_NGML / 10)^e_conc_ada_cl *
      (BMI / 22.5)^e_bmi_cl
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)

    # Bioavailability carries its own covariate model and its own eta,
    # correlated with the clearance eta (Song 2025 Table 2 Cov CL-F).
    fdepot <- exp(lfdepot + etalfdepot) *
      (CONC_ADA_NGML / 10)^e_conc_ada_f *
      (BMI / 22.5)^e_bmi_f

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # Bioavailability applies to the subcutaneous depot only. IV doses go
    # straight to central with an implied F of 1, which is the identifiability
    # anchor that lets Song 2025 estimate the SC F from a joint IV + SC fit
    # (Methods section 6: "Both IV and SC data were used simultaneously in
    # estimating CL, Vc, Vp, Q, F and ka").
    f(depot) <- fdepot

    # Dose in mg, volume in L -> mg/L, which is numerically identical to the
    # ug/mL used throughout Song 2025.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

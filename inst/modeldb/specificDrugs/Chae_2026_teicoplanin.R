Chae_2026_teicoplanin <- function() {
  description <- paste(
    "Two-compartment IV population PK model for teicoplanin in 405 Korean",
    "neutropenic adults after haematopoietic stem cell transplantation (HSCT),",
    "built from routine therapeutic-drug-monitoring trough data (Chae 2026).",
    "CKD-EPI eGFR and serum albumin both enter clearance as POWER terms",
    "normalised to the cohort means:",
    "CL = 1.26 * (eGFR / 113.31)^0.574 * (albumin / 3.10)^-0.793 L/h,",
    "with albumin in g/dL. Between-subject variability is exponential on CL",
    "and on the PERIPHERAL volume V2 only; the central volume V1 and the",
    "inter-compartmental clearance Q carry no random effect. Residual error is",
    "combined additive plus proportional. Because only trough samples were",
    "available, the authors supplemented the 568 observed troughs with 360",
    "early post-dose concentrations simulated from two previously published",
    "Korean teicoplanin models, then obtained the final estimates by stochastic",
    "simulation and estimation (SSE) - re-drawing the simulated backbone and",
    "re-fitting 1000 times. The clearance of 1.26 L/h is roughly double the",
    "0.63-0.69 L/h reported in non-HSCT Korean cohorts, which the authors",
    "attribute to the supranormal mean eGFR of this population",
    "(113 mL/min/1.73 m^2), itself partly an artefact of the reduced muscle",
    "mass and low serum creatinine typical of HSCT recipients."
  )
  reference <- paste(
    "Chae H, Cha HJ, Kang M, Han S, Lee DG. Pharmacokinetic model based on",
    "stochastic simulation and estimation for therapeutic drug monitoring of",
    "teicoplanin in Korean neutropenic hematopoietic stem cell transplant",
    "recipients. Drug Des Devel Ther. 2026. doi:10.2147/DDDT.S550736"
  )
  vignette <- "Chae_2026_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate calculated with the CKD-EPI formula; retained on clearance as a power term",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eGFR. BSA-NORMALISED renal function in mL/min/1.73 m^2",
        "derived with the CKD-EPI equation (Chae 2026 Methods, 'Patient",
        "Population and Observed Data': 'Estimated glomerular filtration rate",
        "(eGFR) was calculated using the CKD-EPI (Chronic Kidney Disease",
        "Epidemiology Collaboration) formula'), stored under the canonical CRCL",
        "column per inst/references/covariate-columns.md, which accepts a",
        "BSA-normalised creatinine-based estimate alongside measured GFR. Supply",
        "the value on that normalised scale; a raw Cockcroft-Gault mL/min value",
        "silently rescales the renal term.",
        "Enters the final model as the power term printed above Chae 2026",
        "Table 3: CL = theta1 * (eGFR / 113.31)^theta5 * (Albumin / 3.10)^theta6.",
        "The normalising constant 113.31 mL/min/1.73 m^2 is the cohort mean; the",
        "Table 1 demographics row prints 113.38 +/- 20.80, a 0.06% difference",
        "attributable to the 12 subjects excluded for missing height. The",
        "equation's 113.31 is used here because it is the constant against which",
        "theta5 was estimated.",
        "This reference value is well above the register's previously observed",
        "range (76-141 mL/min/1.73 m^2) at its upper end, and deliberately so:",
        "the Discussion explains that HSCT recipients have reduced muscle mass",
        "and correspondingly low serum creatinine (cohort mean 0.63 mg/dL), so a",
        "creatinine-based eGFR runs supranormal in this population and should not",
        "be read as genuine hyperfiltration.",
        "Treated as time-fixed at the subject level. The paper records eGFR among",
        "the demographic variables captured per patient and gives no",
        "within-subject time course for it."
      ),
      source_name        = "eGFR"
    ),
    ALB = list(
      description        = "Serum albumin concentration; retained on clearance as a power term with a negative exponent, consistent with teicoplanin's high plasma protein binding",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'Serum albumin level', reported by Chae 2026 in g/dL",
        "(Table 1: 3.04 +/- 0.41 g/dL). The canonical ALB column is g/L (SI) per",
        "inst/references/covariate-columns.md, so model() applies the inline",
        "conversion alb_gdL <- ALB * 0.1 before the power term, exactly as",
        "Roepcke_2023_rezafungin.R does. Populate ALB in g/L: the reference",
        "subject is 31.0 g/L, i.e. the paper's 3.10 g/dL.",
        "Enters the final model as the power term printed above Chae 2026",
        "Table 3: CL = theta1 * (eGFR / 113.31)^theta5 * (Albumin / 3.10)^theta6,",
        "with theta6 = -0.793 (SSE final, Table 4). The NEGATIVE exponent is the",
        "expected direction for a drug that is roughly 90% albumin-bound: higher",
        "albumin means a smaller unbound fraction and therefore a lower apparent",
        "total-drug clearance. The Abstract states the covariate pair reflects",
        "'the effects of renal function and protein binding', and the sign is",
        "consistent in both the base model (theta6 = -0.801, Table 3) and the SSE",
        "final model (-0.793, Table 4), whose bootstrap 95% CI (-1.11 to -0.488)",
        "excludes zero.",
        "The normalising constant 3.10 g/dL differs slightly from the Table 1",
        "cohort mean of 3.04 g/dL (2% higher), the same small offset seen for",
        "eGFR; the equation's constant is used because theta6 was estimated",
        "against it.",
        "Treated as time-fixed at the subject level; no within-subject albumin",
        "time course is reported."
      ),
      source_name        = "Serum albumin level"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Cohort mean 164.82 +/- 9.23 cm (Chae 2026 Table 1; 12 subjects with missing height were excluded). Listed among the covariates tested by GAM screening and PsN stepwise covariate modelling (Methods, 'Base Model Building': 'Tested covariates included height, weight, sex, albumin, serum creatinine level, and estimated creatinine clearance') and not retained in the final model."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Cohort mean 64.06 +/- 12.22 kg (Chae 2026 Table 1). Screened by GAM plus PsN stepwise covariate modelling and not retained. No allometric scaling appears anywhere in the final model, so the published volumes and clearances are absolute values for a roughly 64 kg adult rather than per-70 kg values."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "195 of 417 initially screened subjects were female (222 male / 195 female, Chae 2026 Table 1). Screened by GAM plus PsN stepwise covariate modelling and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Cohort mean 0.63 +/- 0.30 mg/dL (Chae 2026 Table 1). Screened by GAM plus PsN stepwise covariate modelling and not retained; the CKD-EPI eGFR derived from it was retained instead, so creatinine enters the model only through CRCL. The Discussion flags the low cohort creatinine as a consequence of the reduced muscle mass of HSCT recipients and therefore as the reason the derived eGFR runs high."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 405L,
    n_studies        = 1L,
    n_concentrations = 568L,
    age_mean         = "49.29 +/- 14.86 years (Chae 2026 Table 1)",
    age_range        = "Adults aged 18 years and older (Methods inclusion criterion); no upper bound reported",
    weight_mean      = "64.06 +/- 12.22 kg (Chae 2026 Table 1)",
    height_mean      = "164.82 +/- 9.23 cm (Chae 2026 Table 1; 12 subjects with missing height excluded)",
    sex_female_pct   = 46.8,
    race_ethnicity   = "Korean (single-centre Korean cohort; the paper frames the model throughout as population-specific to Koreans and cites inter-ethnic variability in glycopeptide PK as the reason non-Korean models are not transportable here)",
    disease_state    = paste(
      "Adults with haematologic malignancies who had undergone haematopoietic",
      "stem cell transplantation, were neutropenic, were hospitalised for",
      "infection and received intravenous teicoplanin. Teicoplanin was given",
      "for (i) positive Gram-positive cultures, (ii) severe sepsis or shock with",
      "blood cultures pending, (iii) a history of MRSA infection or",
      "colonisation, (iv) skin and soft-tissue infection, or (v) suspected",
      "catheter-related infection. Patients whose teicoplanin PK was judged",
      "unrepresentative of the group - specifically those on massive fluid",
      "therapy or renal replacement therapy - were EXCLUDED, so the model",
      "carries no dialysis-clearance term and should not be applied to patients",
      "on RRT."
    ),
    dose_range       = paste(
      "Teicoplanin intravenously. The standard institutional regimen was three",
      "400 mg loading doses at 12 h intervals followed by 400 mg once daily,",
      "but the actual dose was at the physician's discretion and the standard",
      "regimen was not always applied; the complete real dosing history and",
      "sampling times were reconstructed and used in the analysis. The paper",
      "does not state the infusion duration."
    ),
    regions          = "Korea (Bone Marrow Transplantation Center, now Catholic Hematology Hospital, Seoul St. Mary's Hospital, The Catholic University of Korea; TDM data collected 2015-2017)",
    renal_function   = paste(
      "Supranormal on a creatinine basis. CKD-EPI eGFR mean 113.38 +/- 20.80",
      "mL/min/1.73 m^2 and serum creatinine mean 0.63 +/- 0.30 mg/dL (Chae 2026",
      "Table 1). The Discussion compares this with eGFR means of 64 and 103",
      "mL/min/1.73 m^2 in the comparator literature and attributes the",
      "difference partly to reduced muscle mass and low serum creatinine in HSCT",
      "recipients rather than to true hyperfiltration. Patients receiving renal",
      "replacement therapy were excluded from the analysis."
    ),
    nutritional_state = "Hypoalbuminaemic on average: serum albumin mean 3.04 +/- 0.41 g/dL (30.4 g/L), below the usual 3.5-5.0 g/dL reference interval, which is the clinical context for the retained albumin effect on clearance.",
    screened_covariates = paste(
      "Tested by generalized additive model screening plus PsN stepwise",
      "covariate modelling with forward selection (p < 0.05, dOFV > 3.84) and",
      "backward elimination (p < 0.01, dOFV > 6.63): height, weight, sex,",
      "albumin, serum creatinine and estimated creatinine clearance (Chae 2026",
      "Methods, 'Base Model Building'). Only eGFR and albumin, both on CL, were",
      "retained. Chae 2026 Table 2 gives the OFV ladder against the",
      "no-IIV two-compartment base model (OFV 4369.751): IIV on CL -486.746,",
      "IIV on CL and V2 -565.546, plus eGFR on CL -600.598, plus albumin on CL",
      "-588.965, plus both -634.415."
    ),
    notes            = paste(
      "SINGLE-CENTRE RETROSPECTIVE TDM STUDY WITH A SIMULATED BACKBONE. 417",
      "patients were initially included; after excluding 12 with missing",
      "covariates and 1 with an erroneously recorded concentration, 405 patients",
      "contributed 568 valid observations. Sampling times ranged from 19 to 2879",
      "h after the first dose; TDM samples were drawn about 1 h before a",
      "scheduled dose (roughly 23 h after the previous maintenance dose) after",
      "at least 48-96 h on a stable regimen, though 7.6% of records came from",
      "patients maintained under 48 h. Assay: UHPLC-MS/MS, calibration range",
      "3.9-52.9 mg/L with a separately validated LLOQ of 0.72 mg/L.",
      "BECAUSE ONLY TROUGHS WERE AVAILABLE, the observed data could not identify",
      "the distribution phase. The authors therefore simulated 360 backbone",
      "concentrations (6 timepoints - 0.3, 1, 1.5, 2.5, 4 and 6 h post-dose - in",
      "each of 60 virtual subjects, 30 drawn from each of the two previously",
      "published Korean teicoplanin models: a two-compartment model in 15",
      "critically ill elderly patients and a three-compartment model in 12",
      "healthy volunteers) and merged them with the observed troughs. Covariates",
      "for the virtual subjects were drawn from multivariate normal",
      "distributions matched to those published models using R 4.4.0 and",
      "tmvtnorm 1.6.0.",
      "The FINAL parameter estimates come from the SSE procedure: the",
      "simulate-merge-refit cycle was repeated 1000 times in NONMEM 7.5 with",
      "FOCE-I; 88.1% of runs both minimised successfully and returned reliable",
      "estimates, and the mean and 90% CI across those runs are Chae 2026",
      "Table 4. RSEs ranged 0.522-8.88%, with 2.09% on CL and 6.81% on V2.",
      "A three-compartment model fitted the pooled data better (lower OFV,",
      "better residual normality and homoscedasticity) but destabilised the SSE",
      "loop because many subjects were adequately described by two compartments;",
      "the authors selected the two-compartment structure for estimation",
      "stability and explicitly acknowledge the resulting structural",
      "misspecification. Qualification: 1000-sample bootstrap and VPC on the base",
      "model (Table 3, Figures 1-2) and a 32-bin prediction-corrected VPC on the",
      "SSE final model (Figure 3), which the authors describe as performing best",
      "within 720 h of the first dose - the window in which teicoplanin TDM is",
      "actually done - with prediction intervals inflating beyond it where data",
      "are sparse."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects. Chae 2026 Table 4 ('Parameter Estimates of the
    # Final Pharmacokinetic Model Using SSE (N = 1000)'), Median column. The
    # Median column is used throughout rather than the Mean column because it
    # is the one the paper's own narrative quotes: the Discussion states
    # 'V1 was estimated at 4.71 L' and 'clearance was estimated at 1.26 L/h',
    # matching Table 4 medians 4.71 and 1.26 (the means are 4.71 and 1.25).
    # Table 3 gives the corresponding BASE-model estimates, which the SSE
    # procedure then refined; the SSE values are the paper's final model.
    #
    # The clearance equation printed above Chae 2026 Table 3 is
    #   CL = theta1 * (eGFR / 113.31)^theta5 * (Albumin / 3.10)^theta6
    # so lcl is the clearance of a REFERENCE subject at the cohort-mean
    # eGFR of 113.31 mL/min/1.73 m^2 and albumin of 3.10 g/dL (31.0 g/L).
    # ------------------------------------------------------------------
    lcl <- log(1.26); label("Clearance at eGFR 113.31 mL/min/1.73 m^2 and albumin 3.10 g/dL (L/h)")  # Chae 2026 Table 4 theta1 median 1.26 L/h (mean 1.25 +/- 0.01, RSE 0.522%); base-model Table 3 gave 1.25 (RSE 1.98%), bootstrap median 1.24, 95% CI 1.2-1.3
    lvc <- log(4.71); label("Central volume V1 (L)")                                                 # Chae 2026 Table 4 theta2 median 4.71 L (mean 4.71 +/- 0.114, RSE 5.8%); base-model Table 3 gave 4.58 (RSE 5.61%), bootstrap median 4.57, 95% CI 4.06-5.09
    lq  <- log(5.50); label("Inter-compartmental clearance Q (L/h)")                                 # Chae 2026 Table 4 theta4 median 5.50 L/h (mean 5.50 +/- 0.051, RSE 2.59%); base-model Table 3 gave 5.53 (RSE 2.62%), bootstrap median 5.53, 95% CI 5.23-5.84
    lvp <- log(46.4); label("Peripheral volume V2 (L)")                                              # Chae 2026 Table 4 theta3 median 46.4 L (mean 46.3 +/- 0.058, RSE 2.97%); base-model Table 3 gave 47.1 (RSE 10.64%), bootstrap median 46.6, 95% CI 33.2-60.9

    # ------------------------------------------------------------------
    # Covariate effects on CL, both POWER exponents on ratios to the cohort
    # means, per the equation printed above Chae 2026 Table 3. Applied in
    # model() as (CRCL / 113.31)^e_crcl_cl * (alb_gdL / 3.10)^e_alb_cl, which
    # reproduces that equation verbatim.
    #
    # e_alb_cl is NEGATIVE as printed: teicoplanin is about 90% albumin-bound,
    # so a higher albumin lowers the unbound fraction and hence the apparent
    # total clearance. The sign is stable across the base model (-0.801) and
    # the SSE final model (-0.793), and the base-model bootstrap 95% CI
    # (-1.11 to -0.488) excludes zero.
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.574;  label("Power exponent of CKD-EPI eGFR on CL (unitless)")   # Chae 2026 Table 4 theta5 median 0.574 (mean 0.582 +/- 0.05, RSE 2.55%); base-model Table 3 gave 0.596 (RSE 20.3%), bootstrap median 0.609, 95% CI 0.339-0.853
    e_alb_cl  <- -0.793; label("Power exponent of serum albumin on CL (unitless)")  # Chae 2026 Table 4 theta6 median -0.793 (mean -0.795 +/- 0.027, RSE 1.39%); base-model Table 3 gave -0.801 (RSE 19.9%), bootstrap median -0.804, 95% CI -1.11 to -0.488

    # ------------------------------------------------------------------
    # Between-subject variability, log-normal (Methods, 'Base Model Building':
    # 'The individual variability of each PK parameter was described using a
    # lognormal distribution ... eta_i is the between-subject variability
    # (BSV), which follows a normal distribution with a mean of zero and a
    # variance of omega^2').
    #
    # Chae 2026 Table 4 reports these as NONMEM omega^2 VARIANCES, not as CV%.
    # Table 3 reports the same two random effects under a header that says
    # explicitly '(as CV%)' - 35.7% for CL and 78.9% for V2 - while Table 4's
    # header carries no such qualifier and prints 0.118 and 0.458. Converting
    # the Table 4 variances the log-normal way confirms the identification:
    #   sqrt(exp(0.118) - 1) = 35.4% vs Table 3's 35.7%
    #   sqrt(exp(0.458) - 1) = 76.2% vs Table 3's 78.9%
    # (the small residual gap is the base-model-versus-SSE difference). Reading
    # them as SDs instead would give 34.4% and 67.7%, and the V2 figure is then
    # 11 points off its Table 3 counterpart.
    #
    # IIV is on CL and on the PERIPHERAL volume only. Chae 2026 Table 2 names
    # the selected structure '2-compartment with IIV in Cl, V2', and Table 3
    # labels theta3 (the parameter V2 carries an eta on) 'Peripheral volume
    # (L)' while theta2 under the header 'V1' is 'Central volume (L)'. So the
    # eta sits on lvp, not on lvc. No off-diagonal covariances were published,
    # although the Methods state an OMEGA BLOCK was evaluated.
    # ------------------------------------------------------------------
    etalcl ~ 0.118  # Chae 2026 Table 4 omega_CL median 0.118 (mean 0.118 +/- 0.041, RSE 2.09%); equals 35.4% CV, matching Table 3's printed 35.7% CV for the base model
    etalvp ~ 0.458  # Chae 2026 Table 4 omega_V2 median 0.458 (mean 0.452 +/- 0.133, RSE 6.81%); equals 76.2% CV, matching Table 3's printed 78.9% CV for the base model

    # ------------------------------------------------------------------
    # Residual variability: combined additive plus proportional (Methods:
    # 'Proportional, additive, and combined error models (epsilon) were
    # assessed for residual variability, which follows a normal distribution
    # with a mean of zero and a variance of sigma^2'; both an additive and a
    # proportional row appear in the final Table 4).
    #
    # SCALE: Chae 2026 prints sigma_add = 0.672 and sigma_prop = 0.245, which
    # are NONMEM $SIGMA VARIANCES, so the SDs encoded here are their square
    # roots. Three things establish this.
    #  (i) Table 4 demonstrably reports raw NONMEM output for its random
    #      effects: its omega rows are variances (see the BSV note above), and
    #      the residual rows sit in the same table under the same unqualified
    #      header.
    # (ii) Table 3 flags the ONE place the authors transformed anything -- its
    #      BSV subsection header reads '(as CV%)' -- and leaves the residual
    #      subsection header bare.
    #(iii) The reported RSE is on the same scale as the printed number, i.e.
    #      NONMEM's own variance-scale RSE. For sigma_prop in Table 3, the
    #      estimate 0.236 with RSE 6.18% implies an asymptotic 95% interval of
    #      0.236 +/- 1.96 * 0.0146 = 0.207-0.265, which reproduces the printed
    #      bootstrap 95% CI of 0.204-0.268 almost exactly. Had the authors
    #      converted the estimate to an SD they would have had to delta-method
    #      the RSE too, and the two intervals would not line up.
    # The alternative reading -- that 0.245 is already an SD, i.e. a 24.5%
    # proportional residual -- is the more typical published magnitude and
    # matches the 24.3% of the sibling Korean teicoplanin model
    # Wi_2017_teicoplanin.R, so it is recorded here as the competing
    # interpretation and in the vignette's Assumptions section. It is not
    # adopted: a roughly 50% proportional residual is consistent with what this
    # dataset actually is -- retrospectively reconstructed dose histories and
    # real-world sampling times, fitted with a two-compartment structure the
    # authors themselves describe as misspecified relative to the
    # three-compartment fit.
    # ------------------------------------------------------------------
    addSd  <- 0.8198; label("Additive residual error (mg/L)")          # Chae 2026 Table 4 sigma_add median 0.672 (mean 0.671 +/- 0.174, RSE 8.88%) read as a variance -> SD = sqrt(0.672) = 0.8198 mg/L. Base-model Table 3 gave 0.741 (RSE 25.9%), bootstrap median 0.72, 95% CI 0.247-1.23. For scale, the assay LLOQ is 0.72 mg/L.
    propSd <- 0.4950; label("Proportional residual error (fraction)")  # Chae 2026 Table 4 sigma_prop median 0.245 (mean 0.245 +/- 0.054, RSE 2.74%) read as a variance -> SD = sqrt(0.245) = 0.4950. Base-model Table 3 gave 0.236 (RSE 6.18%), bootstrap median 0.234, 95% CI 0.204-0.268.
  })

  model({
    # Individual PK parameters. CKD-EPI eGFR (canonical column CRCL, in
    # mL/min/1.73 m^2) and serum albumin (canonical column ALB, in SI g/L)
    # both act on CL as power terms normalised to the cohort means printed in
    # the Chae 2026 clearance equation. Albumin is converted from the canonical
    # SI g/L to the g/dL the paper calibrated against, as required by
    # inst/references/covariate-columns.md; the reference subject is
    # ALB = 31.0 g/L = 3.10 g/dL.
    alb_gdL <- ALB * 0.1  # canonical SI g/L -> the paper's US-convention g/dL

    cl <- exp(lcl + etalcl) * (CRCL / 113.31)^e_crcl_cl * (alb_gdL / 3.10)^e_alb_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Intravenous dosing into the central compartment, first-order
    # distribution to peripheral1 and first-order elimination from central.
    # Chae 2026 does not report the infusion duration, so the route is left to
    # the event table: dose central as a bolus (amt only) or as an infusion
    # (amt plus rate or dur). Dose in mg and volumes in L give central / vc in
    # mg/L, matching the mg/L concentrations the paper reports.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

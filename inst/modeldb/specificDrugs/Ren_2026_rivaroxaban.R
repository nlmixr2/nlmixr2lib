Ren_2026_rivaroxaban <- function() {
  description <- paste(
    "Joint population pharmacokinetic-pharmacodynamic model for rivaroxaban",
    "in real-world Chinese adults treated for pulmonary embolism (Ren 2026).",
    "The PK layer is a one-compartment oral model with first-order absorption",
    "and elimination; creatinine clearance enters CL/F and age enters V/F, both",
    "as median-normalized power terms. The PD layer is a direct linear",
    "concentration-prothrombin time relationship, PT = baseline + slope * Cc,",
    "with alanine aminotransferase entering the slope as a median-normalized",
    "power term. The two endpoints (plasma rivaroxaban concentration and",
    "prothrombin time) are declared jointly; the authors fitted the PD layer",
    "sequentially on the PK posterior and then re-estimated the combined",
    "PK/PD model. Interindividual variability was not supported on ka",
    "(99.8% shrinkage) and was dropped by the authors.",
    sep = " "
  )
  reference <- paste(
    "Ren J, Li Y, Zheng X, Wu Z, Shi J, Han X.",
    "Population pharmacokinetic and pharmacodynamic analysis of rivaroxaban in",
    "real-world patients with pulmonary embolism in China.",
    "Res Pract Thromb Haemost. 2026;10:106618.",
    "doi:10.1016/j.rpth.2026.106618.",
    sep = " "
  )
  vignette <- "Ren_2026_rivaroxaban"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Rivaroxaban was given orally with a meal, 5-20 mg once
  # daily (Methods 2.1), and measured in citrated plasma by UPLC-MS/MS over
  # 1-500 ng/mL (Methods 2.2).
  compartmentData <- list(
    depot   = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance, Cockcroft-Gault, raw (NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Computed with the Cockcroft-Gault equation (Methods 2.1); the paper",
        "reports it in raw mL/min and never normalizes to 1.73 m^2, so a user",
        "supplying BSA-normalized eGFR would mis-scale CL/F. Enters CL/F as",
        "(CRCL / 88.3)^0.270, where 88.3 mL/min is the cohort median (Results",
        "3.2, text following Equation 14). Renal strata in the cohort: >= 80",
        "mL/min 58.3%, 50-79 mL/min 34.2%, 30-49 mL/min 7.5% (Table 1);",
        "patients with CrCL < 15 mL/min were excluded (Methods 2.1), so the",
        "model carries no information about severe renal impairment or",
        "dialysis. The paper's own dose-adjustment simulations evaluate the",
        "term at 30, 50 and 80 mL/min.",
        sep = " "
      ),
      source_name        = "CrCL"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters V/F as (AGE / 63)^-0.231, where 63 years is the cohort median",
        "(Results 3.2, text following Equation 13). Observed range 21-92 years,",
        "mean 58.8 +/- 15.8 (Table 1). Age was also picked up on CL/F and on ka",
        "during forward inclusion but both were removed in backward elimination",
        "(OFV drops of only 5.935 and 1.362); only the V/F effect is in the",
        "final model.",
        sep = " "
      ),
      source_name        = "Age"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "PD covariate only: enters the concentration-PT slope as",
        "(ALT / 19)^-0.201, where 19 IU/L is the cohort median (Results 3.3,",
        "text following Equation 15). The source reports IU/L; the canonical",
        "register uses U/L, which is numerically interchangeable. Observed",
        "range 10-157 IU/L, mean 25.3 +/- 19.1 (Table 1); patients with",
        "Child-Pugh B/C hepatic impairment were excluded (Methods 2.1). The",
        "authors judge the effect not clinically significant -- raising ALT",
        "from 10 to 157 IU/L lowers average steady-state PT by only ~7%",
        "(Results 3.5, Figure 6) -- but it is retained in the final model, so",
        "it is retained here.",
        sep = " "
      ),
      source_name        = "ALT"
    )
  )

  # Screened in the stepwise covariate search but NOT retained in the final
  # model (Methods 2.3 lists the screened set; Results 3.2 and 3.3 report the
  # retained ones). ALT is absent from this list because it WAS retained, on
  # the PD slope.
  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.3). Median 164 cm, range 153-192 (Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.3). Median 70 kg, range 49-124 (Table 1). The Discussion notes body weight as a likely driver of the V/F difference against Caucasian cohorts, but it did not survive the paper's own screen."
    ),
    SEXF = list(
      description = "Biological sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (paper covariate 'gender') and not retained (Methods 2.3). Cohort 121/187 female (64.7%)."
    ),
    HGB = list(
      description = "Hemoglobin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.3). Median 132 g/L, range 84-177 (Table 1)."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.3). Median 39.5%, range 24.1-52 (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.3). Median 43 g/L, range 32-49 (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.3). Median 10.7 umol/L, range 4.1-50.4 (Table 1)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and not retained (Methods 2.3). Median 3.3 umol/L, range 0.9-16.6 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened and not retained as a covariate in its own right (Methods 2.3); renal function entered the final model through the Cockcroft-Gault CRCL derived from it. Median 68 umol/L, range 36-155 (Table 1)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 187L,
    n_observations   = "517 rivaroxaban plasma concentrations and 376 prothrombin time assays (Results 3.2). The two counts differ because residual plasma volume was sometimes insufficient for the coagulation assay (Results 3.1).",
    n_studies        = 1L,
    age_range        = "21-92 years (median 63, mean 58.8 +/- 15.8); 48.1% aged 65 years or older, 12.8% aged 75 years or older",
    weight_range     = "49-124 kg (median 70, mean 70.8 +/- 14.7)",
    sex_female_pct   = 64.7,
    race_ethnicity   = c(Asian = 100),
    disease_state    = "Acute pulmonary embolism (a manifestation of venous thromboembolism), treated in routine clinical care. Patients with Child-Pugh B/C hepatic impairment, creatinine clearance < 15 mL/min, pregnancy or breastfeeding, spontaneous bleeding tendency or platelets < 20 x 10^9/L, contraindications to other factor Xa inhibitors, or hereditary thrombophilia / antiphospholipid syndrome were excluded (Methods 2.1).",
    dose_range       = "Oral rivaroxaban 5 mg (3.7%), 10 mg (34.2%), 15 mg (8.6%) or 20 mg (53.5%) once daily, taken with a meal for at least 5 days before sampling; dose was chosen by the treating physician (Methods 2.1, Table 1).",
    regions          = "China (single center: Peking Union Medical College Hospital, Beijing), prospective cohort enrolled April 2021 to August 2024.",
    renal_function   = "Cockcroft-Gault creatinine clearance: >= 80 mL/min in 58.3%, 50-79 mL/min in 34.2%, 30-49 mL/min in 7.5%; cohort median 88.3 mL/min. 41.7% had mild-to-moderate renal impairment.",
    hepatic_function = "ALT median 19 IU/L (10-157); albumin median 43 g/L (32-49); total bilirubin median 10.7 umol/L (4.1-50.4). Moderate/severe hepatic impairment was an exclusion criterion.",
    co_medication    = "Statins 26.2%, metformin 5.3%, metoprolol 5.3%, aspirin 3.7%, strong P-glycoprotein inhibitors 3.2%, strong BCRP inhibitors 1.6% (Table 1). Too few patients received strong P-gp or BCRP inhibitors for the interaction to be modelled (Discussion).",
    notes            = paste(
      "Baseline prothrombin time in the cohort: median 13.8 s, range 9.6-31.9",
      "(Table 1); PT was assayed with Thromborel S reagent on a Sysmex CS5100",
      "with within-run and inter-run precision < 3.5%. Sparse sampling: 2-4",
      "samples per patient at predose (within 0.5 h) and 2, 4, 6, 12 h",
      "(+/- 0.5 h) after a dose. Comorbidities: cancer 22.9%, hypertension",
      "16.0%, hyperlipidemia 11.2%, diabetes 10.2%, coronary artery disease",
      "4.3%. The model was fitted in Phoenix WinNonlin 8.3 by first-order",
      "conditional estimation extended least squares. PT sensitivity to",
      "rivaroxaban is reagent-dependent; this model's slope is specific to the",
      "Thromborel S / Sysmex CS5100 system used here (Discussion).",
      sep = " "
    )
  )

  ini({
    # =================================================================
    # PK structural parameters (Ren 2026 Table 2). Time in h, dose in
    # mg, CL/F in L/h, V/F in L. ka was ESTIMATED here (CV 13.5%,
    # bootstrap 0.922-1.58), unlike the two earlier Chinese rivaroxaban
    # models that fixed it at the Japanese value of 0.617 1/h -- the
    # Discussion argues explicitly against reusing that fixed value.
    # =================================================================
    lka <- log(1.25);  label("Absorption rate constant (ka, 1/h)")         # Table 2: Ka = 1.25 (CV 13.5%, bootstrap 0.922-1.58)
    lvc <- log(53.0);  label("Apparent volume of distribution (V/F, L)")   # Table 2: V/F = 53.0 (CV 3.91%, bootstrap 49.0-57.1)
    lcl <- log(6.13);  label("Apparent clearance (CL/F, L/h)")             # Table 2: CL/F = 6.13 (CV 3.18%, bootstrap 5.75-6.52)

    # =================================================================
    # PK covariate effects. Methods Equation 4 gives the continuous-
    # covariate form P_i = P_pop * (COV / COV_median)^theta, instantiated
    # in Results 3.2 as
    #   Equation 13:  V/F  (L)   = 53.0 * (Age  / 63  )^-0.231
    #   Equation 14:  CL/F (L/h) = 6.13 * (CrCL / 88.3)^ 0.270
    # NOTE both signs: the age exponent is NEGATIVE. The Elsevier PDF's
    # symbol font drops the minus sign under some text extractors (the
    # preprocessed _trimmed.md renders Table 2 as a bare "0.231"), but
    # the typeset Equation 13, the Results prose ("the coefficient (theta)
    # of age on V/F value was -0.231, indicating that the V/F value
    # decreased with age") and the bootstrap CI (-0.391 to -0.070) all
    # agree that it is negative.
    # =================================================================
    e_age_vc  <- -0.231; label("Age power exponent on V/F (unitless)")                    # Table 2 / Equation 13: -0.231 (CV 34.7%, bootstrap -0.391 to -0.070)
    e_crcl_cl <-  0.270; label("Creatinine clearance power exponent on CL/F (unitless)")  # Table 2 / Equation 14: 0.270 (CV 23.9%, bootstrap 0.142-0.391)

    # =================================================================
    # PD structural parameters (Ren 2026 Table 3). Methods Equation 5
    # (the retained linear form, chosen over the near-linear Equation 6
    # and the Emax / sigmoid-Emax Equations 8-9):
    #   Equation 15: PT (s) = 11.3 + 0.0184 * (ALT / 19)^-0.201 * C
    # with C in ug/L (= ng/mL), which is why Cc is scaled to ng/mL in
    # model() below.
    #
    # SLOPE VALUE -- the paper is self-inconsistent. Results 3.3 prose
    # says "a slope of 0.00184 s per 1 ug/L", but Table 3, Equation 15
    # and the cross-study comparison in Table 5 all print 0.0184. The
    # prose figure is a typo: at 0.00184 the model predicts an average
    # steady-state PT of ~11.6 s at every renal stratum, which sits below
    # the 1st percentile of every box in the paper's own Figure 6,
    # whereas 0.0184 reproduces the Figure 6 panel centres (14.6 s at
    # CrCL 30, 13.8 s at CrCL 80, both at ALT 20 IU/L). 0.0184 is used.
    # =================================================================
    lrbase <- log(11.3);   label("Baseline prothrombin time with no rivaroxaban present (s)")     # Table 3: Baseline = 11.3 s (CV 0.711%, bootstrap 11.1-11.5)
    lslope <- log(0.0184); label("Linear slope of prothrombin time on concentration (s per ug/L)") # Table 3 / Equation 15: Slope = 0.0184 (CV 3.20%, bootstrap 0.0174-0.0191)

    # ALT exponent is negative -- same dropped-minus-sign caveat as the
    # age exponent above; Equation 15, the Results prose ("the
    # coefficient (theta) of ALT on slope was -0.201") and the bootstrap
    # CI (-0.305 to -0.100) all agree.
    e_alt_slope <- -0.201; label("Alanine aminotransferase power exponent on the PT slope (unitless)")  # Table 3 / Equation 15: -0.201 (CV 26.1%, bootstrap -0.305 to -0.100)

    # =================================================================
    # Interindividual variability. Methods Equation 1 is exponential
    # (P_i = P_pop * exp(eta_i)); Tables 2 and 3 head the block
    # "Interindividual variability (%)" and report CV%, mapped here via
    # omega^2 = log(CV^2 + 1):
    #   V/F      CV 12.5%   -> 0.0155042
    #   CL/F     CV 14.6%   -> 0.0210920
    #   Baseline CV  0.242% -> 0.00000585638
    #   Slope    CV 10.4%   -> 0.0107579
    # No eta on ka: "Due to the high shrinkage (99.8%) in the
    # interindividual variation of Ka values, the random effect of Ka
    # values among individuals was removed from the model" (Results 3.2),
    # and Table 2 prints "-" on that row.
    #
    # The baseline-PT eta is essentially zero (CV 0.242%, i.e. +/- 0.03 s
    # on an 11.3 s baseline) and is kept only for fidelity to Table 3.
    # See the vignette's "Assumptions and deviations" section for why the
    # 0.242 is read as 0.242% rather than as a fraction: at 24.2% the
    # simulated steady-state PT distribution would span roughly 9.0-18.7 s
    # (5th-95th), about three times the spread the authors actually show
    # in Figure 6, whereas 0.242% reproduces it.
    # =================================================================
    etalvc    ~ 0.0155042        # Table 2: IIV on V/F, CV 12.5% (CV of estimate 16.8%, shrinkage 15.9%)
    etalcl    ~ 0.0210920        # Table 2: IIV on CL/F, CV 14.6% (CV of estimate 13.0%, shrinkage 6.6%)
    etalrbase ~ 0.00000585638    # Table 3: IIV on Baseline, CV 0.242% (CV of estimate 35.4%, shrinkage 46.7%)
    etalslope ~ 0.0107579        # Table 3: IIV on Slope, CV 10.4% (CV of estimate 30.5%, shrinkage 31.7%)

    # =================================================================
    # Residual error. Methods Equation 2 is proportional on both
    # endpoints (Y = IPRED * (1 + eps)). Tables 2 and 3 label the row
    # "Proportional error (%)" but print the standard deviation as a
    # FRACTION, not a percentage: read literally, 0.221% would mean the
    # sparse real-world concentration data were fitted to within a
    # quarter of a percent, which is finer than the assay's own
    # inter/intraday precision of < 8.5% (Methods 2.2). The same applies
    # to the PD row against the < 3.5% PT-assay precision. Encoded as
    # 22.1% and 6.69% CV respectively.
    # =================================================================
    propSd    <- 0.221;  label("Proportional residual SD on rivaroxaban concentration (fraction)")  # Table 2: proportional error 0.221 (CV of estimate 8.50%, shrinkage 28.1%)
    propSd_PT <- 0.0669; label("Proportional residual SD on prothrombin time (fraction)")           # Table 3: proportional error 0.0669 (CV of estimate 7.61%, shrinkage 21.0%)
  })

  model({
    # =================================================================
    # PK: one-compartment oral disposition with first-order absorption
    # and elimination (Results 3.2). Covariates enter as the median-
    # normalized power terms of Equations 13 and 14.
    # =================================================================
    ka <- exp(lka)
    vc <- exp(lvc + etalvc) * (AGE  / 63)^e_age_vc      # Equation 13
    cl <- exp(lcl + etalcl) * (CRCL / 88.3)^e_crcl_cl   # Equation 14

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x1000 converts to ng/mL, which is
    # the unit the assay was validated in (1-500 ng/mL, Methods 2.2) and
    # is numerically identical to the ug/L in which the PD slope is
    # reported.
    Cc <- central / vc * 1000

    # =================================================================
    # PD: direct linear concentration-prothrombin time relationship with
    # no effect compartment and no time delay (Methods Equation 5,
    # instantiated as Equation 15). Cc here is the individual prediction
    # -- residual error is applied to each endpoint separately below, so
    # the PT prediction is driven by the error-free concentration, which
    # matches the authors' sequential PK-then-PD fitting strategy
    # (Methods 2.4).
    # =================================================================
    rbase <- exp(lrbase + etalrbase)
    slope <- exp(lslope + etalslope) * (ALT / 19)^e_alt_slope   # Equation 15

    PT <- rbase + slope * Cc

    Cc ~ prop(propSd)
    PT ~ prop(propSd_PT)
  })
}

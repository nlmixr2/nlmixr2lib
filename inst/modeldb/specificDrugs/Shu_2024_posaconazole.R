Shu_2024_posaconazole <- function() {
  description <- "Population PK model for posaconazole oral suspension (Noxafil) in Chinese hematopoietic stem cell transplantation (HSCT) recipients (Shu 2024). One-compartment disposition with first-order absorption and first-order elimination, parameterised on the apparent (oral) scale as CL/F and V/F because the therapeutic-drug-monitoring data were oral-only and bioavailability was not identifiable. The absorption rate constant Ka was fixed to 0.4 1/h from earlier posaconazole-suspension popPK studies, because almost every sample was a pre-dose trough and the absorption phase could not support an estimate. Creatinine clearance enters apparent clearance as a power function centred on the cohort median 103.81 mL/min; body weight enters apparent volume as a power function centred on the cohort median 45.85 kg; and concomitant proton-pump-inhibitor use multiplies apparent volume by 3.83. Inter-individual variability is exponential on CL/F and V/F and is very large (omega 1.118 and 0.826 on the log scale). Residual variability is proportional. The model was used for Monte Carlo dose optimisation of weight-banded BID and TID regimens against steady-state trough targets of 0.7 ug/mL for prophylaxis and 1.0 ug/mL for treatment of invasive fungal disease."
  reference   <- "Shu YS, Dong ZH, Yang YL, Li SW, Yi QY, Wang P, Shi YP, Zhang YY, Shi HY. Individualized regimen of Posaconazole oral suspension in Chinese HSCT patients based on population pharmacokinetic model. Sci Rep. 2024;14:20288. doi:10.1038/s41598-024-70955-w."
  vignette    <- "Shu_2024_posaconazole"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "posaconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "posaconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance calculated with the Cockcroft-Gault equation. NOT body-surface-area normalised: the source reports raw Cockcroft-Gault mL/min (Shu 2024 Materials and methods, 'Study design and patients': 'Ccr was calculated using the Cockcroft-Gault equation'). The paper separately tabulates a CKD-EPI eGFR column, which was screened but not retained in the final model.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject in the source analysis. Retained on apparent clearance in the final model",
        "as a power function centred on the cohort median: CL/F = 11.5 * (CCR / 103.81)^0.68",
        "(Shu 2024 Table 2 footnote). Cohort median 103.81 mL/min, range 14.44-240.64 mL/min",
        "(Shu 2024 Table 1); do not extrapolate outside that range.",
        "The raw (non-BSA-normalised) Cockcroft-Gault form matches the CRCL register's Delattre 2010 amikacin",
        "precedent; the register requires the assay form to be documented per model, which this note does.",
        "The authors flag this covariate as a departure from earlier posaconazole popPK models",
        "(Shu 2024 Discussion): only about 13 percent of an oral suspension dose is excreted renally,",
        "so the CCR effect is more plausibly a marker of general physiological reserve in HSCT recipients",
        "than a mechanistic renal-elimination pathway.",
        sep = " "
      ),
      source_name        = "CCR"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject in the source analysis. Retained on apparent volume of distribution",
        "in the final model as a power function centred on the cohort median:",
        "V/F proportional to (WT / 45.85)^1.78 (Shu 2024 Table 2 footnote).",
        "Cohort median 45.85 kg, range 13.65-80 kg (Shu 2024 Table 1) -- a wide range because the cohort",
        "spans ages 3-65 years and therefore includes children. The estimated exponent 1.78 is far above",
        "the allometric 1.0 expected for a volume term and reflects posaconazole's high lipophilicity and",
        "adipose-tissue distribution (Shu 2024 Discussion), plus the confounding of V/F with bioavailability",
        "in oral-only data. Do not extrapolate outside the observed weight range.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CONMED_PPI = list(
      description        = "Concomitant proton-pump-inhibitor use (1 = patient co-administered a PPI, 0 = no PPI)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant PPI; 51 of 62 patients, Shu 2024 Table 1 'Proton pump inhibitor (No. of patients) 11')",
      notes              = paste(
        "Subject-level ever-versus-never indicator as reported in Shu 2024 Table 1; the paper does not",
        "define a per-record time-varying flag or a minimum-duration threshold, so the operational",
        "definition is looser than the Goel 2016 sonidegib precedent (>= 80 percent of the PK assessment phase).",
        "Retained on apparent volume of distribution: Shu 2024 Table 2 reports Vppi = 3.34 L for PPI = 0",
        "and 12.8 L for PPI = 1, entering the printed V/F equation as a multiplicative factor.",
        "Mechanistically the effect is a bioavailability effect -- PPIs raise gastric pH, reducing the",
        "solubility and hence the absorption of the posaconazole suspension (Shu 2024 Discussion) --",
        "which in an oral-only dataset appears on the apparent volume V/F.",
        "See the ini() block and the vignette Assumptions and deviations section for the reference-level",
        "arithmetic used to encode the two-level Vppi.",
        sep = " "
      ),
      source_name        = "PPI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search but not retained in the final model (Shu 2024 Results, 'Model building'). Cohort median 20.5 years, range 3-65 (Shu 2024 Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). 30 of 62 patients were male (Shu 2024 Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 19.17 kg/m^2, range 10.74-28.34 (Shu 2024 Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 27.15 U/L, range 2.7-188.3 (Shu 2024 Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 19.2 U/L, range 2.5-234.1 (Shu 2024 Table 1)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 107 U/L, range 45-243 (Shu 2024 Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 40.1 g/L, range 5.1-49.3 (Shu 2024 Table 1)."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 32.2 U/L, range 9-327 (Shu 2024 Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 11 umol/L, range 4.00-134.7 (Shu 2024 Table 1)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). Cohort median 4.9 umol/L, range 1.37-114.9 (Shu 2024 Table 1)."
    ),
    CONMED_PHENYTOIN = list(
      description = "Concomitant phenytoin sodium use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). 12 of 62 patients (Shu 2024 Table 1)."
    ),
    CONMED_METOCLOPRAMIDE = list(
      description = "Concomitant metoclopramide use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Shu 2024 Results, 'Model building'). 7 of 62 patients (Shu 2024 Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 62L,
    n_studies      = 1L,
    age_range      = "3-65 years (median 20.5; Shu 2024 Table 1)",
    age_median     = "20.5 years",
    weight_range   = "13.65-80 kg (median 45.85; Shu 2024 Table 1)",
    weight_median  = "45.85 kg",
    sex_female_pct = 51.6,
    race_ethnicity = "Chinese (single-centre cohort at the First Affiliated Hospital of Shandong First Medical University, Jinan)",
    disease_state  = "Patients undergoing hematopoietic stem cell transplantation (HSCT) receiving posaconazole oral suspension for prophylaxis or treatment of invasive fungal disease, for at least 7 days before sampling.",
    dose_range     = "150-600 mg/day posaconazole oral suspension (Noxafil), given twice daily (BID) or three times daily (TID), usually with meals at approximately 08:00, 11:30 and 18:30.",
    regions        = "China (Jinan, Shandong), October 2021 to April 2023",
    n_observations = "103 posaconazole plasma concentrations from 62 patients (1-3 samples per patient, at least one trough each). Samples were drawn 0.5 h before the first dose of the day (defined as Cmin) and 2 or 4 h after that dose; 98 of the 103 points were predominantly Cmin. Observed concentrations 0.1-4.03 ug/mL, median 0.706 ug/mL (Shu 2024 Results, 'Demographic characteristics').",
    renal_function = "Cockcroft-Gault creatinine clearance median 103.81 mL/min (range 14.44-240.64); CKD-EPI eGFR median 129.98 mL/min (range 33.38-282.08) (Shu 2024 Table 1).",
    hepatic_function = "ALT median 27.15 U/L (2.7-188.3), AST median 19.2 U/L (2.5-234.1), ALP median 107 U/L (45-243), GGT median 32.2 U/L (9-327), total bilirubin median 11 umol/L (4.00-134.7), direct bilirubin median 4.9 umol/L (1.37-114.9), albumin median 40.1 g/L (5.1-49.3) (Shu 2024 Table 1).",
    co_medication  = "Proton pump inhibitor 11 of 62 patients; phenytoin sodium 12 of 62; metoclopramide 7 of 62 (Shu 2024 Table 1). Only PPI was retained as a covariate.",
    notes          = paste(
      "Retrospective, single-centre, therapeutic-drug-monitoring study (Ethics approval R202312200209).",
      "Posaconazole was quantified by a validated HPLC-UV assay (261 nm) with a calibration range of",
      "0.1-10 ug/mL and a lower limit of quantification of 0.1 ug/mL (Shu 2024 Materials and methods).",
      "Model estimated in NONMEM 7.3.0; internal validation by 1000-sample nonparametric bootstrap",
      "(83.7 percent success), NPDE, and prediction-corrected VPC. No external validation was performed,",
      "which the authors name as the main limitation of the study.",
      "Baseline demographics are in Shu 2024 Table 1.",
      sep = " "
    )
  )

  ini({
    # ====================================================================
    # Structural PK parameters (Shu 2024 Table 2).
    #
    # One-compartment model on the APPARENT (oral) scale. Posaconazole was
    # given only as an oral suspension, so bioavailability F is not
    # identifiable and the reported disposition parameters are CL/F and
    # V/F (Shu 2024 Results, 'Model building'). No separate bioavailability
    # term is encoded: the nominal dose enters the depot compartment and F
    # is absorbed into lcl and lvc.
    # ====================================================================

    # Ka was NOT estimated. Shu 2024 Materials and methods, 'Population
    # pharmacokinetic modeling': 'The absorption rate constant (Ka) was
    # fixed to 0.4 /h'. Table 2 accordingly prints '0.4 FIX' with no RSE
    # and no bootstrap summary. The Discussion explains why: almost every
    # sample was a pre-dose trough, 'Due to the small sample size of the
    # absorption phase collected in this study, it was not possible to
    # robustly estimate the value of the absorption rate constant Ka', and
    # the value 0.4 /h was taken from earlier posaconazole-suspension
    # popPK literature.
    lka <- fixed(log(0.4)); label("First-order absorption rate constant Ka (1/h; taken from earlier posaconazole-suspension popPK studies)")  # Shu 2024 Table 2: Ka = 0.4 1/h FIX

    lcl <- log(11.5);       label("Apparent clearance CL/F at the cohort-median CCR of 103.81 mL/min (L/h)")  # Shu 2024 Table 2: CLpop = 11.5 L/h, RSE 22.3 percent, bootstrap median 10.63 (95% CI 4.55-17.2)

    # --------------------------------------------------------------------
    # lvc: the reference-level apparent volume, 2768.86 L.
    #
    # ARITHMETIC AND WHY IT IS NOT THE PRINTED 829. The Shu 2024 Table 2
    # footnote prints the final volume model as
    #
    #     V/F (L) = 829 * e^0.826 * (WT/45.85)^1.78 * Vppi
    #
    # where Table 2 gives Vppi = 3.34 for PPI = 0 and Vppi = 12.8 for
    # PPI = 1. (The 'e^0.826' factor is the inter-individual random effect
    # written with omega substituted for eta -- 0.826 is exactly the
    # IIVV row, 82.6 percent; the same substitution appears as 'e^1.12' in
    # the CL/F equation, matching IIVCL = 111.8 percent. Both are eta terms
    # with mean 0, i.e. a factor of 1 at the typical value.)
    #
    # Taken literally, the typical V/F of a median-weight patient NOT on a
    # PPI is therefore 829 * 3.34 = 2768.86 L, and Vpop = 829 alone is not
    # attained by anyone. Because PPI = 0 is the reference category, the
    # reference-level factor 3.34 is folded into lvc here and only the
    # PPI = 1 / PPI = 0 RATIO is carried as the covariate effect. This is
    # exact algebra, not a reparameterisation of the fit: 829 and Vppi are
    # confounded in the printed model -- they only ever appear as a product.
    #
    # CONFLICT WITHIN TABLE 2 ITSELF. The Table 2 legend defines Vpop as
    # 'volume of distribution of the typical subject' -- the same wording it
    # uses for CLpop, where the label IS consistent with the equation because
    # the CCR factor equals 1 at the cohort median. For V/F it is not, because
    # Vppi is 3.34 at its reference level and never 1. The Discussion repeats
    # the legend's number ('The typical value of the V/F population of HSCT
    # patients evaluated was 829 L', compared there with 548 L for the tablet
    # formulation). So the legend and the equation of the same table disagree
    # and no reading satisfies both. The printed EQUATION is implemented here,
    # per the standing convention that where a source's equation and its
    # description conflict, the equation is the implementable statement.
    #
    # An additive reading of Vppi is eliminated outright: its '(L)' unit would
    # give 829 + 3.34 = 832 L versus 829 + 12.8 = 842 L, a 1.2 percent effect
    # that could not survive a backward elimination requiring dOFV > 10.83.
    #
    # The strongest corroboration is internal to Table 2: Vpop is
    # reported with an RSE of 5.2 percent, which implies an interval of
    # roughly 745-915 L, yet its own bootstrap 95% CI is 123.55-1520 L --
    # more than an order of magnitude, and irreconcilable with a 5.2 percent
    # RSE unless Vpop is confounded with another estimated parameter. Vppi
    # is the only candidate, and its CI (0.99-16.53 for PPI = 0) is
    # correspondingly wide, while the two bootstrap MEDIANS multiply to
    # 763 * 4.05 = 3090 L -- close to 2768.86 L and nowhere near 829 L.
    # That is the signature of a fit in which Vpop and Vppi were both
    # estimated as free multiplicative thetas, i.e. of the literal equation.
    #
    # WHAT DOES *NOT* SETTLE IT, to save a future reader the detour:
    #   (1) Cross-model V/F comparison is uninformative. Kohl 2010, the only
    #       other posaconazole oral-suspension model here, has V/F = 2250 L,
    #       superficially close to 2768.86 L -- but its paired CL/F is
    #       67 L/h, so its half-life is 23 h against 167 h here. An apparent
    #       volume is only interpretable together with the clearance it is
    #       paired with, and on the half-life axis the comparison actually
    #       favours the 829 L reading (50 h). Neither reading reproduces
    #       posaconazole's usual 20-35 h, because trough-only data cannot
    #       inform V/F at all -- the same unidentifiability that produced
    #       the confounding above and the implausible WT exponent of 1.78.
    #   (2) Reproducing the paper's Monte Carlo regimens does not
    #       discriminate either, because the target quantity is a
    #       steady-state TROUGH and the trough of a one-compartment model is
    #       governed by dose/(CL*tau) once the half-life is long relative to
    #       the interval. At 3 mg/kg TID and the reference covariates the two
    #       readings differ by under 2 percent (about 1.48 vs 1.46 ug/mL).
    #       The vignette runs the whole PTA grid under BOTH readings to show
    #       this explicitly rather than claiming the replication as support.
    #       That run recovers IDENTICAL recommendations in all eight strata,
    #       so the ambiguity, while real and unresolvable from the source,
    #       changes nothing the paper concludes. A user who prefers the other
    #       reading gets it in one line: rxode2::ini(mod, lvc = log(829)).
    # --------------------------------------------------------------------
    lvc <- log(2768.86);    label("Apparent volume of distribution V/F at the cohort-median weight of 45.85 kg without concomitant PPI (L)")  # Shu 2024 Table 2: Vpop = 829 L (RSE 5.2%) x Vppi(PPI=0) = 3.34 (RSE 8.7%) = 2768.86 L

    # ====================================================================
    # Covariate effects (Shu 2024 Table 2 and its footnote equations).
    # ====================================================================

    # CL/F = 11.5 * (CCR / 103.81)^0.68, centred on the Table 1 cohort
    # median CCR. Table 2 labels this row 'CLccr, magnitude of the effect
    # of creatinine clearance on CL'.
    e_crcl_cl <- 0.68;  label("Power exponent of creatinine clearance on apparent clearance (unitless; CL/F ~ (CRCL/103.81)^e_crcl_cl)")  # Shu 2024 Table 2: CLccr = 0.68, RSE 32.1 percent, bootstrap median 0.88 (95% CI 0.4-3.21)

    # V/F proportional to (WT / 45.85)^1.78, centred on the Table 1 cohort
    # median weight. Table 2 labels this row 'Vwt, magnitude of the effect
    # of bodyweight on V'.
    e_wt_vc   <- 1.78;  label("Power exponent of body weight on apparent volume of distribution (unitless; V/F ~ (WT/45.85)^e_wt_vc)")  # Shu 2024 Table 2: Vwt = 1.78, RSE 16.2 percent, bootstrap median 1.58 (95% CI 0.82-2.39)

    # PPI effect on V/F, expressed as the ratio of the two printed Vppi
    # level values because PPI = 0 is the reference category and its value
    # (3.34) has been folded into lvc above:
    #     12.8 / 3.34 = 3.832335...
    # Applied as e_conmed_ppi_vc^CONMED_PPI so that PPI = 0 contributes a
    # factor of 1 and PPI = 1 multiplies V/F by 3.83.
    e_conmed_ppi_vc <- 3.832335; label("Fold-change in apparent volume of distribution with concomitant PPI versus no PPI (unitless)")  # Shu 2024 Table 2: Vppi = 12.8 (PPI=1, RSE 20.0%) / 3.34 (PPI=0, RSE 8.7%)

    # ====================================================================
    # Inter-individual variability (Shu 2024 Table 2, rows IIVCL and IIVV).
    #
    # SCALE. The paper describes IIV with 'a power exponential error model'
    # (Materials and methods, Eq. 1) -- i.e. theta_i = theta_pop *
    # exp(eta_i) with eta ~ N(0, omega^2) -- and reports IIVCL and IIVV as
    # percentages. The Table 2 footnote equations print those same numbers
    # inside the exponential as 'e^1.12' (CL/F) and 'e^0.826' (V/F), which
    # identifies the reported percentages as omega itself (the standard
    # deviation of the log-scale random effect) expressed as a percentage,
    # NOT as a CV% requiring the omega^2 = log(CV^2 + 1) conversion. The
    # nlmixr2 variances below are therefore simply (IIV% / 100)^2.
    #
    # The magnitudes are very large -- omega = 1.118 on CL/F corresponds to
    # a 95 percent range spanning roughly a 80-fold spread of clearance --
    # and are consistent with the paper's own VPC commentary ('the large
    # shaded area in the graph indicates a high degree of variability') and
    # with the well-documented erratic absorption of the suspension
    # formulation. Reported eta-shrinkage was 25.2 percent (CL) and
    # 32.1 percent (V), so these are only moderately informed by the data.
    # ====================================================================
    etalcl ~ 1.118^2  # Shu 2024 Table 2: IIVCL = 111.8 percent, RSE 17.3 percent, shrinkage 25.2 percent; variance = 1.249924
    etalvc ~ 0.826^2  # Shu 2024 Table 2: IIVV  =  82.6 percent, RSE 12.7 percent, shrinkage 32.1 percent; variance = 0.682276

    # ====================================================================
    # Residual error. Shu 2024 Results, 'Model building': 'The model
    # follows the Residual unknown variability (RUV) of the proportional
    # error model'. Four candidate residual models (additive,
    # proportional, power exponential, mixed) were compared by OFV and
    # diagnostic plots (Materials and methods, Eqs. 2-5).
    # ====================================================================
    propSd <- 0.137; label("Proportional residual error (fraction)")  # Shu 2024 Table 2: Proportional error = 13.7 percent, RSE 13.6 percent, shrinkage 35.2 percent
  })

  model({
    # ====================================================================
    # 1. Individual PK parameters.
    #
    # Both covariate relationships are power functions centred on the
    # Table 1 cohort medians, exactly as printed in the Shu 2024 Table 2
    # footnote:
    #     CL/F = 11.5 * (CCR/103.81)^0.68            * exp(eta_CL)
    #     V/F  = 829  * (WT /45.85 )^1.78 * Vppi     * exp(eta_V)
    # with the PPI = 0 level of Vppi folded into lvc (see ini()).
    # ====================================================================
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 103.81)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 45.85)^e_wt_vc * e_conmed_ppi_vc^CONMED_PPI

    # ====================================================================
    # 2. Micro-constants
    # ====================================================================
    kel <- cl / vc

    # ====================================================================
    # 3. ODE system -- one compartment with first-order absorption.
    # Oral suspension dose records target cmt = depot with the nominal
    # dose in mg; bioavailability is absorbed into the apparent parameters
    # CL/F and V/F, so no f(depot) term is applied.
    # ====================================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ====================================================================
    # 4. Observation and residual error.
    # Cc (ug/mL) = central (mg) / vc (L), since mg/L == ug/mL.
    # ====================================================================
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

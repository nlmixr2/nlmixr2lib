Lai_2026_rivaroxaban <- function() {
  description <- "One-compartment population PK model for rivaroxaban in Asian (Taiwanese) adults with atrial fibrillation sampled under real-world therapeutic drug monitoring, with a creatinine-clearance power effect and a CYP3A4/P-gp-inhibitor comedication effect on CL/F and a lean-body-weight effect on V/F (Lai 2026)"
  reference <- "Lai NS, Lin CJ, Kuo CH, Peng YF, Tang SC, Huang CF, Lin SY, Lin SW. Development and Clinical Application of a Real-World Population Pharmacokinetic Model of Rivaroxaban in Asian Patients with Atrial Fibrillation. Clin Pharmacokinet. 2026. doi:10.1007/s40262-026-01650-4"
  vignette <- "Lai_2026_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Rivaroxaban is given orally (10 or 15 mg once daily,
  # Table 1 "Rivaroxaban regimen and concentration") and measured in plasma by
  # validated UHPLC-MS/MS (Methods 2.2 "Study Procedure").
  compartmentData <- list(
    depot   = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the Cockcroft-Gault equation",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT BSA-normalized to",
        "mL/min/1.73 m^2. Methods 2.2: 'To estimate creatinine clearance, the",
        "Cockcroft-Gault equation was utilized'. Enters CL/F as a power of the ratio to",
        "the cohort median 54 mL/min; the median is named in the Discussion ('An",
        "increase in CrCL from the median (54 mL/min) to the maximum (133 mL/min)",
        "resulted in a 75% increase in CL/F, while a decrease to the minimum",
        "(19.8 mL/min) led to a 46% reduction'). Cohort mean 53.7 +/- 18.9 mL/min",
        "(Table 1), observed range 19.8-133 mL/min (Discussion); 44.7% of patients had",
        "CrCL < 50 mL/min (Results 3.1). MDRD-4 eGFR and serum creatinine were screened",
        "as alternative renal descriptors and rejected: 'CrCL was selected as a better",
        "predictor of the impact on renal function than CRE or eGFR' (Results 3.3.2).",
        "Patients with end-stage renal disease requiring dialysis were excluded",
        "(Methods 2.1), so the model carries no information about dialysis."
      ),
      source_name        = "CrCL"
    ),
    FFM = list(
      description        = "Lean body weight (fat-free mass) from the Janmahasatian equation",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The source calls this column 'LBW' and derives it with the Janmahasatian",
        "formula, which is the equation this register files under FFM: Methods 2.2",
        "'To estimate fat-free body mass, we evaluated two descriptors: the",
        "Janmahasatian formula (lean body weight, LBW) and the James formula (lean body",
        "mass, LBM) ... The former is considered the contemporary standard in PPK",
        "studies'. The James-formula LBM is the distinct canonical LBM and was screened",
        "as an alternative and rejected (Results 3.3.1). Requires body weight, height",
        "and sex to compute; the source's Table S2 holds the printed equations and is",
        "not on disk, but the Janmahasatian form is fixed by the cited reference",
        "(Janmahasatian et al., Clin Pharmacokinet 2005;44:1051-1065). Enters V/F as an",
        "exponential (Monolix log-linear 'proportional shift') deviation from the cohort",
        "median 48 kg, which the Discussion names ('Increasing LBW from the median",
        "(48 kg) to the maximum (81.6 kg) ... reducing LBW to the minimum (24.9 kg)').",
        "Observed range 24.9-81.6 kg (Discussion). Total body weight, BSA, IBW, the",
        "James LBM and sex were all screened on V/F; LBW gave the largest OFV drop",
        "(Results 3.3.1) and total body weight was explicitly not significant",
        "(Discussion: 'actual bodyweight did not significantly affect V/F')."
      ),
      source_name        = "LBW"
    ),
    CONMED_CYP3A4_PGP_INH = list(
      description        = "Concomitant CYP3A4 or P-glycoprotein inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant CYP3A4 or P-gp inhibitor)",
      notes              = paste(
        "1 = the patient was receiving a CYP3A4 OR a P-glycoprotein inhibitor at the",
        "time of rivaroxaban concentration measurement, 0 = neither. The source pools",
        "the two inhibition mechanisms into a single indicator (Methods 2.2:",
        "'Concurrent medications were classified into four categories: cytochrome P450",
        "3A4 (CYP3A4) or P-glycoprotein (P-gp) inhibitors, ...'), so the column is the",
        "UNION of the two mechanisms and is not interchangeable with either",
        "CONMED_CYP3A4_INH or CONMED_PGP_INH alone. 86 of 226 patients (38.1%) were",
        "positive (Table 1). The agents were amiodarone (n = 56), dronedarone (n = 15),",
        "diltiazem (n = 12) and verapamil (n = 6) (Results 3.1); all four inhibit P-gp,",
        "and dronedarone, diltiazem and verapamil are additionally moderate CYP3A4",
        "inhibitors while amiodarone is a weak one, so in THIS cohort the two mechanisms",
        "coincide for every positive patient. The full classification list is the",
        "source's Table S1, which is not on disk. The Discussion attributes most of the",
        "effect to the two antiarrhythmics: 'these effects were primarily driven by",
        "amiodarone and dronedarone, which accounted for 82.6% of co-medication cases'.",
        "Concomitant BCRP inhibitors (benzbromarone, febuxostat, sulfasalazine; 22",
        "patients, 9.7%) were screened as a SEPARATE covariate on CL/F and eliminated",
        "in the backward step (Results 3.3.2), and are therefore NOT part of this",
        "indicator. The CYP3A4/P-gp inducer arm (levetiracetam, 4 patients, 1.8%) is",
        "likewise not part of this indicator and was not retained."
      ),
      source_name        = "CYP3A4 or P-gp inhibitors"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final model.
  # Results 3.3.1 (on V/F) and 3.3.2 (on CL/F) name each one; the paper reports
  # only OFV drops for these (Tables S4 and S5, which are not on disk) and no
  # point estimate, so none can be encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F and entered the full model, but removed in backward",
        "elimination (Results 3.3.2: 'the full model incorporated age, IBW, CrCL,",
        "CYP3A4/P-gp inhibitor use, BCRP inhibitor use, and sex ... the removal of CrCL",
        "or the concomitant use of CYP3A4/P-gp inhibitors led to significant increases",
        "in the OFV'). Mean 73.2 +/- 8.2 years, 38.1% older than 75 years (Table 1,",
        "Results 3.1)."
      )
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened on BOTH V/F and CL/F and significant on each in forward inclusion,",
        "but dropped from the final model for collinearity with the body-size",
        "descriptors (Discussion: 'Although sex was also identified as a significant",
        "factor influencing V/F, it was excluded from the final model because of its",
        "high collinearity with LBW'). 124 of 226 patients male (54.9%), i.e. 45.1%",
        "female (Table 1). Sex is still needed indirectly, as an input to the",
        "Janmahasatian FFM equation."
      )
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened and NOT retained on V/F (Discussion: 'Interestingly, actual",
        "bodyweight did not significantly affect V/F, possibly because of the",
        "moderate-to-low tissue affinity of rivaroxaban, making LBW a more accurate",
        "predictor of distribution'). Mean 65.0 +/- 13.1 kg (Table 1). Still needed",
        "indirectly, as an input to the Janmahasatian FFM equation."
      )
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Not itself a screened covariate, but a required input to the Janmahasatian FFM",
        "equation and to the Cockcroft-Gault / MDRD-4 derivations. Recorded here because",
        "a user supplying real data must have it: one validation-cohort patient was",
        "excluded precisely for 'missing essential covariates (height)' (Results 3.4).",
        "The paper does not tabulate a height distribution."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Screened on V/F and significant in forward inclusion, but rejected in favour of",
        "LBW (Results 3.3.1: 'As BSA, IBW, LBW, and LBM were highly correlated, LBW",
        "resulted in the greatest reduction in OFV and was therefore selected for further",
        "modeling'). No distribution tabulated."
      )
    ),
    IBW = list(
      description = "Ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened on BOTH V/F and CL/F. On V/F it lost to LBW (Results 3.3.1). On CL/F",
        "it was actually the winning body-size descriptor in forward inclusion ('Among",
        "IBW, LBM, and LBW, IBW was selected as the best body size descriptor on the",
        "basis of the greatest reduction in OFV') and entered the full model, but it was",
        "then removed in backward elimination -- only CrCL and CYP3A4/P-gp-inhibitor use",
        "survived on CL/F (Results 3.3.2). The final model therefore carries no body-size",
        "term on clearance at all. No distribution tabulated, and the IBW formula variant",
        "is in the source's Table S2, which is not on disk."
      )
    ),
    LBM = list(
      description = "Lean body mass from the James equation",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "The James-formula lean body mass, screened on BOTH V/F and CL/F as the",
        "'traditional body-size descriptor' alternative to the Janmahasatian LBW",
        "(Methods 2.2) and rejected on both (Results 3.3.1, 3.3.2). Distinct from the",
        "retained FFM column, which is the Janmahasatian quantity. No distribution",
        "tabulated."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Reported as 'CRE' in Methods 2.2. Screened on CL/F as an alternative renal",
        "descriptor and rejected in favour of CrCL (Results 3.3.2). Mean 1.1 +/- 0.3",
        "mg/dL with 0% missing (Table 1). Still needed indirectly, as the input to the",
        "Cockcroft-Gault CrCL used by the retained CRCL column. MDRD-4 eGFR was the",
        "third renal descriptor screened and is not listed separately here because it is",
        "an alias of the retained CRCL canonical."
      )
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Not screened in the primary covariate search because 30.1% of values were",
        "missing (Table 1 footnote a; Methods 2.8). The sensitivity analysis found it to",
        "be the ONLY missing-data laboratory parameter significantly correlated with",
        "CL/F, and the authors deliberately declined to add it because the retained CRCL",
        "already carries renal function (Results 3.6: 'since CrCL, a primary indicator of",
        "renal function, was already successfully incorporated into the final PPK model,",
        "BUN was not further evaluated as an additional covariate'). Mean 21.5 +/- 8.7",
        "mg/dL. Recorded here because it is the one screened-out covariate the paper's",
        "own diagnostics flag as carrying residual signal."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 226,
    n_studies      = 1,
    n_observations = 452,
    age_range      = "over 20 years (eligibility); mean 73.2 +/- 8.2 years",
    age_median     = "not reported; mean 73.2 years",
    weight_range   = "not reported; mean 65.0 +/- 13.1 kg",
    weight_median  = "not reported; mean 65.0 kg",
    sex_female_pct = 45.1,
    race_ethnicity = c(Asian = 100),
    disease_state  = "atrial fibrillation, on rivaroxaban for at least 3 days",
    dose_range     = "10 mg (37.6%) or 15 mg (62.4%) orally once daily",
    regions        = "Taiwan (single center: National Taiwan University Hospital, Taipei)",
    renal_function = "Cockcroft-Gault CrCL mean 53.7 +/- 18.9 mL/min, median 54, range 19.8-133 mL/min; 44.7% below 50 mL/min. Patients with end-stage renal disease requiring dialysis were excluded.",
    co_medication  = "CYP3A4 or P-gp inhibitors in 86 patients (38.1%): amiodarone 56, dronedarone 15, diltiazem 12, verapamil 6. CYP3A4/P-gp inducer (levetiracetam) in 4 (1.8%). BCRP inhibitors in 22 (9.7%): benzbromarone 13, febuxostat 8, sulfasalazine 1.",
    notes          = paste(
      "Prospective observational DOAC-T (direct oral anticoagulant-Taiwan) cohort",
      "enrolled January 2016 to June 2023; ClinicalTrials.gov NCT05333666. Table 1",
      "baseline demographics. 258 patients enrolled, 226 analysed after exclusions",
      "(Fig. 1). Sampling is deliberately sparse and paired: exactly two steady-state",
      "samples per patient, one peak drawn 2-4 h after a pharmacist-supervised dose",
      "(observed mean sampling time 2.11 +/- 0.40 h) and one trough drawn immediately",
      "before the next scheduled dose (observed mean 24.6 +/- 3.14 h since the last",
      "dose). Observed trough concentration median 27.1 ng/mL (IQR 15.6-50.5); observed",
      "peak median 232.8 ng/mL (IQR 162.7-343.2). 12 of the 464 collected samples were",
      "excluded before model fitting -- 11 below the limit of quantification (10 trough,",
      "1 peak) and 1 peak outlier of 1225.1 ng/mL -- leaving 452 (Results 3.2).",
      "Risk scores: CHA2DS2-VASc median 3 (IQR 2-5), HAS-BLED median 2 (IQR 1-2).",
      "A separate 67-patient validation cohort (65 analysable) was drawn from the same",
      "database and is NOT part of the 226; individual-level MAPE was 10.24% and",
      "population-level MAPE 91.1% (Results 3.4). Only one post-dose sample per patient",
      "lies in the absorption phase, which is why ka could not be estimated.",
      "Approximately 30% of the cohort received off-label underdose regimens",
      "(Discussion)."
    )
  )

  ini({
    # Structural parameters. Time in h, dose in mg, CL/F in L/h, V/F in L.
    # ka was NOT estimated: "Given that only one post-dose peak sample was
    # collected and absorption kinetics may not be adequately characterized, the
    # absorption rate constant (Ka) was fixed at 0.821 h-1 according to
    # literature findings [8]" (Methods 2.3). Table 2 prints the rounded 0.82 and
    # marks it "(fixed)" in both the base and final columns; Methods 2.3 carries
    # the un-rounded 0.821, which is the value used here.
    lka <- fixed(log(0.821)); label("Absorption rate constant (ka, 1/h), from the literature (Lai 2026 reference 8)")  # Methods 2.3; Table 2 row "Ka (h-1) 0.82 (fixed)"

    # Typical values AT the covariate reference point, i.e. CrCL = 54 mL/min,
    # no CYP3A4/P-gp inhibitor, LBW = 48 kg. The Abstract states the same two
    # numbers ("The estimated apparent clearance (CL/F) and the volume of
    # distribution (V/F) were 6.13 L/h and 45.57 L, respectively").
    lcl <- log(6.13);  label("Apparent clearance (CL/F, L/h) at CrCL 54 mL/min with no CYP3A4/P-gp inhibitor")  # Table 2 final model CL/F 6.13 (RSE 6.21%); Abstract
    lvc <- log(45.57); label("Apparent volume of distribution (V/F, L) at LBW 48 kg")                            # Table 2 final model V/F 45.57 (RSE 6.21%); Abstract

    # Covariate effects. Monolix 2023R1 parameterizes every covariate effect on
    # the log-parameter scale, so a "power function" on a log-transformed
    # covariate and a "proportional shift" on an untransformed one are the same
    # log-linear machinery:
    #   log(CL_i) = log(CL_pop) + e_crcl_cl * log(CrCL_i / 54)
    #                           + e_..._inh_cl * INH_i          + eta_CL
    #   log(V_i)  = log(V_pop)  + e_ffm_vc * (LBW_i - 48)       + eta_V
    # Results 3.3.2 names the two CL/F forms ("CrCL with a power function and the
    # use of a CYP3A4/P-gp inhibitor with a proportional shift function") and
    # Results 3.3.1 names the V/F form ("the model incorporating LBW using a
    # proportional-shift functional form"). The exponential reading of the
    # categorical shift is confirmed arithmetically by the Discussion: exp(-0.38)
    # = 0.684, i.e. a 31.6% reduction, and the Discussion states "Concurrent use
    # of CYP3A4/P-gp inhibitors reduced CL/F by 32%". A linear (1 + beta) reading
    # would give 38% and does not match.
    e_crcl_cl                  <-  0.64; label("Creatinine-clearance power exponent on CL/F (unitless)")            # Table 2 final model "CrCl on CL/F 0.64 (RSE 19.78%)"; Results 3.3.2
    e_conmed_cyp3a4_pgp_inh_cl <- -0.38; label("Log-scale shift in CL/F when a CYP3A4 or P-gp inhibitor is coadministered (unitless)")  # Table 2 final model "CYP3A4/P-gp inhibitors on CL/F -0.38 (RSE 5.44%)"; Results 3.3.2
    e_ffm_vc                   <-  0.01; label("Log-scale shift in V/F per kg of lean body weight above 48 kg (1/kg)")                  # Table 2 final model "LBW on V/F 0.01 (RSE 41.37%)"; Results 3.3.1

    # IIV, log-normal on both CL/F and V/F ("IIV included the apparent volume of
    # distribution (V/F) and CL/F, both of which were best described by an
    # exponential model", Results 3.2). Table 2 prints each omega to two decimals
    # (0.59 for both in the final model) AND the derived CV% to four significant
    # figures. The CV% column is the more precise of the two, so omega^2 is taken
    # from it via the log-normal identity omega^2 = log(CV^2 + 1):
    #   CL/F CV 64.20% -> omega^2 = log(0.6420^2 + 1) = 0.345123 -> omega 0.5875
    #   V/F  CV 63.93% -> omega^2 = log(0.6393^2 + 1) = 0.342670 -> omega 0.5854
    # Both back-solved omegas round to the printed 0.59, so this reading is
    # consistent with BOTH columns of Table 2; using 0.59^2 = 0.3481 directly
    # would instead imply CV 64.53% for both and contradict the printed CV%.
    # The base model's numbers reproduce the same way (0.71 / 81.42% -> 0.7131).
    # A correlation between the two etas is not reported and was not modelled.
    etalcl ~ 0.345123  # Table 2 final model omega_CL 0.59 (RSE 7.67%), CV 64.20%
    etalvc ~ 0.342670  # Table 2 final model omega_V  0.59 (RSE 6.27%), CV 63.93%

    # Residual error. Monolix's "constant" error model is a pure ADDITIVE error
    # on the observation scale: "The RV can be best described using a constant
    # error model" (Results 3.2). Table 2's residual-variability row is labelled
    # "a" and glossed in the table footnote as "residual unexplained variability
    # (constant error model)", so 14.36 is an additive SD in ng/mL.
    addSd <- 14.36; label("Additive residual error (ng/mL)")  # Table 2 final model a = 14.36 (RSE 17.72%)
  })
  model({
    # Individual parameters. See the ini() comment block for the Monolix
    # log-linear covariate parameterization and for why the categorical shift is
    # exponential rather than linear.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (CRCL / 54)^e_crcl_cl *
      exp(e_conmed_cyp3a4_pgp_inh_cl * CONMED_CYP3A4_PGP_INH)
    vc <- exp(lvc + etalvc) * exp(e_ffm_vc * (FFM - 48))

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x1000 converts to the ng/mL used
    # throughout the paper.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}

Kolowrat_2025_apixaban <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model with first-order ",
    "absorption and linear elimination for oral apixaban 2.5 or 5 mg twice ",
    "daily in hospitalized adults with nonvalvular atrial fibrillation ",
    "(Kolowrat 2025), quantifying the real-world amiodarone drug-drug ",
    "interaction from salvaged clinical plasma samples. Apparent oral ",
    "clearance CL/F = 1.5 L/h at the cohort median age of 77 years without ",
    "amiodarone, scaled by a power function of age (AGE/77)^-1.52 and ",
    "multiplied by exp(-0.4) = 0.670 (a 33% reduction, 95% CI 12% to 48%) ",
    "during concomitant amiodarone 200 mg. Apparent volume of distribution ",
    "V/F = 45.57 L carries no covariate. The absorption rate constant ka was ",
    "fixed to 0.82 1/h from the Gaspar 2023 OptimAT apixaban model because ",
    "its interindividual variability was imprecisely estimated and fixing it ",
    "improved the corrected Bayesian information criterion; see ",
    "modellib('Gaspar_2023_apixaban'). Interindividual variability is ",
    "supported on CL/F (62.23% CV) and V/F (52.66% CV) but not on ka. ",
    "Residual variability is combined proportional (15%) and additive ",
    "(21.28 ng/mL). Two-compartment models were explored but the peripheral ",
    "volume and intercompartmental clearance were not estimable from the ",
    "sparse real-world data. Renal function (estimated glomerular filtration ",
    "rate, creatinine clearance), body weight and body mass index were ",
    "screened and not retained."
  )
  reference <- paste0(
    "Kolowrat S, Riley C, Lam K, Thomson L, Stickle DF, Kraft WK. ",
    "Real-world impact of amiodarone on apixaban population ",
    "pharmacokinetics in hospitalized patients. ",
    "Clin Transl Sci. 2025;18(11):e70392. doi:10.1111/cts.70392. ",
    "The fixed absorption rate constant ka = 0.82 1/h is inherited from ",
    "Gaspar F, Terrier J, Favre S, et al. Population pharmacokinetics of ",
    "apixaban in a real-life hospitalized population from the OptimAT ",
    "study. CPT Pharmacometrics Syst Pharmacol. 2023;12(10):1541-1552. ",
    "doi:10.1002/psp4.13032; see modellib('Gaspar_2023_apixaban')."
  )
  vignette <- "Kolowrat_2025_apixaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Kolowrat 2025 Methods (apixaban
  # measured in plasma by LC-MS/MS, lower limit of quantification 5 ng/mL)
  # and Results 3.2 (one-compartment model with first-order absorption of
  # oral apixaban tablets, no lag time).
  compartmentData <- list(
    depot   = list(analyte = "apixaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as the median-centered power term (AGE / 77)^-1.52 per ",
        "the Kolowrat 2025 Table 2 note 'CL/F_i = 1.5 x (Age/77)^beta_age x ",
        "e^(beta_amio x 1_{CAT=1}) x e^(IIV_CL/F_i)', which instantiates the ",
        "general continuous-covariate form of Methods 2.3, ",
        "log(theta_i) = log(theta_pop) + beta_theta x log(COV_i / ",
        "COV_median) + eta_theta,i. The reference value 77 years is the ",
        "cohort median, stated in Results 3.2 ('the addition of age ",
        "(normalized to a median of 77) on CL resulted in a better model ",
        "(dBICc = -9.72)'). Note that 77 is the OVERALL median and is not ",
        "printed in Table 1, which reports the two treatment groups ",
        "separately: 79 years (IQR 71-86) for apixaban alone and 74 years ",
        "(IQR 64-83) for apixaban plus amiodarone. The exponent is strongly ",
        "NEGATIVE (-1.52), so apparent clearance FALLS with increasing age ",
        "-- the authors attribute this to age-related decline in hepatic ",
        "activity and renal elimination in a cohort whose median age is ",
        "close to 80 (Discussion). Baseline value only: Limitations note ",
        "that covariate data pulled from the electronic medical record ",
        "'only reflected the first reported value per encounter'. The ",
        "fitted age range is not printed; the reported interquartile ranges ",
        "span 64-86 years, so the term is not supported in younger adults ",
        "and extrapolating it downward inflates clearance steeply (at age ",
        "40 the multiplier is 2.6-fold the value at 77)."
      ),
      source_name        = "Age"
    ),
    CONMED_AMIO = list(
      description        = "Concomitant amiodarone therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant amiodarone)",
      notes              = paste0(
        "1 = subject is receiving a stable dose of amiodarone 200 mg taken ",
        "for at least 30 days prior to admission, 0 = apixaban alone. ",
        "Cohort prevalence 51/106 (48.1%) (Table 1). Entered as a binary ",
        "variable, 'i.e., 0 for no and 1 for yes' (Methods 2.3), under the ",
        "categorical-covariate form log(theta_i) = log(theta_pop) + ",
        "beta_theta x 1_{CAT_i=1} + eta_theta,i, so the effect is the ",
        "multiplicative factor exp(-0.4) = 0.670 on CL/F. Time-fixed at the ",
        "analysis baseline: the 30-day pre-admission stability requirement ",
        "and amiodarone's ~50-day half-life (Introduction) mean inhibition ",
        "is already fully established at the first modelled dose, and all ",
        "patients were assumed to be at steady state with both apixaban and ",
        "amiodarone (Methods 2.1). Per the register convention the ",
        "per-paper definition is 'stable amiodarone 200 mg for >= 30 days', ",
        "which is narrower than 'any concurrent amiodarone use' -- the ",
        "amiodarone dose itself was not modelled and Limitations state that ",
        "'the effect of amiodarone dose nor its metabolites on apixaban PK ",
        "were not assessed'. The mechanism the authors propose is combined ",
        "inhibition of CYP3A4 by the minor metabolite ",
        "N-monodesethylamiodarone and of P-glycoprotein by parent ",
        "amiodarone (Introduction, Discussion). Patients prescribed a ",
        "STRONG CYP3A4 or P-gp inducer or inhibitor during admission were ",
        "excluded, so this covariate is not confounded by strong ",
        "perpetrators; mild and moderate perpetrators were present in 31 ",
        "patients and were deliberately NOT modelled 'given the number of ",
        "individuals per category' (Results 3.1, 3.2)."
      ),
      source_name        = "concomitant amiodarone"
    )
  )

  # Covariates screened by Kolowrat 2025 but NOT retained in the final model.
  # Documentation only: these are not referenced in model(). Recorded here so
  # the provenance of the paper's covariate screen is preserved. Results 3.2:
  # 'All other covariates had weak associations or were not statistically
  # significant for inclusion.' No point estimate is printed for any of them,
  # so no effect can be encoded even as an option.
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Renal function (CKD-EPI 2021 eGFR and Cockcroft-Gault creatinine clearance)",
      units       = "mL/min/1.73 m^2 (eGFR) and mL/min (Cockcroft-Gault)",
      type        = "continuous",
      notes       = paste0(
        "TWO renal-function columns were constructed and screened, and ",
        "neither was retained; they are documented under one register entry ",
        "because both are members of this canonical. (1) Estimated ",
        "glomerular filtration rate by the CKD-EPI 2021 equation, ",
        "BSA-normalised mL/min/1.73 m^2; median 48 in both groups (IQR ",
        "23-70 apixaban alone, 24.5-62 plus amiodarone; Table 1). Tested on ",
        "CL/F centered at its median of 48 and rejected: 'adding eGFR ",
        "(normalized to a median of 48) did not [improve the model] (dBICc ",
        "= -0.83)' against a retention threshold of 2 BICc units (Results ",
        "3.2, Methods 2.4). (2) Creatinine clearance by the Cockcroft-Gault ",
        "equation, raw mL/min and NOT BSA-normalised; median 38.4 mL/min ",
        "(IQR 26.6-52.5) apixaban alone and 38.2 (24.0-55.1) plus ",
        "amiodarone. Methods 2.3 specifies a non-standard weight input for ",
        "this column: ideal body weight (IBW) was used rather than actual ",
        "body weight 'to avoid extreme fluctuations associated with this ",
        "method due to body weight', unless actual weight was below IBW, ",
        "and where actual weight was 30% or more above IBW an adjusted body ",
        "weight AdjBW_i = 0.4 x (WT_i - IBW_i) + IBW_i was substituted. A ",
        "user supplying a conventional actual-body-weight Cockcroft-Gault ",
        "value is therefore NOT supplying the same column. The cohort is ",
        "markedly renally impaired (both group medians near 38 mL/min), ",
        "which the Discussion invokes to explain why this model predicts ",
        "higher exposures than the Morath 2025 comparator; the authors also ",
        "note that eGFR WAS significantly associated with CL during model ",
        "building even though it did not survive the BICc criterion in the ",
        "presence of age, so age and renal function are collinear here and ",
        "the retained age term is partly a proxy for renal decline."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Actual reported body weight was collected (median 86.4 kg, IQR ",
        "62.4-93.3 apixaban alone; 83.4 kg, IQR 65.8-96.7 plus amiodarone; ",
        "Table 1) and used to derive ideal and adjusted body weight for the ",
        "Cockcroft-Gault creatinine clearance column (Methods 2.3). Not ",
        "retained on any PK parameter. Patients weighing 120 kg or more ",
        "were excluded by design (Methods), so the model carries no ",
        "information about the obese extreme where an apparent-volume ",
        "effect would be most likely. No point estimate is reported, so no ",
        "effect can be encoded."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste0(
        "Named in the Introduction as one of the three prespecified ",
        "secondary-objective covariates ('potential covariates such as age, ",
        "renal function, and body mass index'). Not retained, and unlike ",
        "age and eGFR no dBICc is reported for it; Table 1 does not ",
        "tabulate BMI, reporting weight and height separately (height 168 ",
        "+/- 11.5 cm apixaban alone, 170 +/- 10.1 cm plus amiodarone). No ",
        "point estimate is reported, so no effect can be encoded."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 106L,
    n_studies        = 1L,
    n_observations   = 360L,
    age_median       = "77 years (overall cohort median, the centering value for the CL/F age term)",
    age_range        = "IQR 71-86 years apixaban alone (median 79); IQR 64-83 years apixaban plus amiodarone (median 74). Full range not reported.",
    weight_median    = "86.4 kg apixaban alone; 83.4 kg apixaban plus amiodarone",
    weight_range     = "IQR 62.4-93.3 kg apixaban alone; IQR 65.8-96.7 kg apixaban plus amiodarone. Weight >= 120 kg excluded by design.",
    sex_female_pct   = 48.1,
    race_ethnicity   = c(
      White = 58.5, Black = 28.3, Hispanic = 2.8, Asian = 3.8, Other = 0.9
    ),
    disease_state    = paste0(
      "Hospitalized adults (age >= 18 years) with a history of nonvalvular ",
      "atrial fibrillation on a stable dose of apixaban 2.5 mg or 5 mg ",
      "twice daily, either alone (n = 55) or with a stable dose of ",
      "amiodarone 200 mg taken for at least 30 days prior to admission (n = ",
      "51). Markedly renally impaired: median estimated glomerular ",
      "filtration rate 48 mL/min/1.73 m^2 and median Cockcroft-Gault ",
      "creatinine clearance about 38 mL/min in both groups. Median serum ",
      "creatinine 1.31 mg/dL (apixaban alone) and 1.37 mg/dL (plus ",
      "amiodarone). Exclusions: concurrent prescription of a STRONG CYP3A4 ",
      "or P-gp inducer or inhibitor during admission, body weight >= 120 ",
      "kg, or apixaban prescribed for venous thromboembolism rather than ",
      "atrial fibrillation. Mild / moderate CYP3A4 and P-gp perpetrators ",
      "were permitted and present in 31 patients (CYP3A4 inhibitor 27.3% / ",
      "29.4%, CYP3A4 inducer 9.1% / 11.8%, P-gp inducer 12.7% / 15.7%, P-gp ",
      "inhibitor 12.7% / 13.7% in the two groups); they were not modelled ",
      "because of the small number of individuals per category."
    ),
    dose_range       = paste0(
      "Oral apixaban 2.5 mg or 5 mg twice daily (every 12 h) at steady ",
      "state per the drug label. All patients were assumed to be at steady ",
      "state with apixaban and amiodarone at the time of the first modelled ",
      "dose (Methods 2.1). Concentrations came from SALVAGED plasma -- ",
      "samples drawn for routine clinical care that had passed a 5-day ",
      "discard threshold -- so sampling times were driven by clinical care ",
      "rather than by a PK schedule and are not restricted to troughs. ",
      "Median 5 samples per patient (range 1-19). Apixaban by LC-MS/MS, ",
      "lower limit of quantification 5 ng/mL; 7 of 360 samples were below ",
      "it and were censored."
    ),
    regions          = "United States (single-centre retrospective observational study, Thomas Jefferson University Hospital, Philadelphia, Pennsylvania; IRB iRISID-2023-2228, approved 6 December 2024)",
    co_medication    = "Concomitant amiodarone 200 mg in 51/106 (48.1%); mild / moderate CYP3A4 or P-gp perpetrators other than amiodarone in 31/106 (29.2%)",
    renal_function   = "Median eGFR 48 mL/min/1.73 m^2 (CKD-EPI 2021) and median Cockcroft-Gault creatinine clearance about 38 mL/min in both groups; a renally impaired cohort",
    notes            = paste0(
      "Baseline demographics in Table 1; base- and final-model parameter ",
      "estimates in Table 2; simulated exposure metrics in Table 3 and ",
      "Figure 2. Patients were screened from electronic medical records ",
      "from December 2023 to July 2024. Software: R 4.4.2 for data ",
      "assembly, Monolix 2024R1 (SAEM) for estimation, Sycomore 2024R1 for ",
      "model management, Simulx for the exposure simulations. Model ",
      "selection used the corrected Bayesian information criterion with a ",
      "retention threshold of 2 BICc units, forward addition then backward ",
      "deletion, supported by diagnostic plots, Pearson correlation and ",
      "Wald tests. No formal power analysis was conducted. Evaluation used ",
      "a dose-normalized visual predictive check stratified by receipt of ",
      "amiodarone (Figure S4) plus goodness-of-fit, IWRES and NPDE plots ",
      "(Figures S1-S3). Structural-model development: two-compartment ",
      "models were explored but 'did not allow PK parameter estimation with ",
      "adequate precision despite removing IIV on each PK parameter' -- ",
      "specifically the peripheral volume and intercompartmental clearance ",
      "were not supported -- so a one-compartment model was selected. IIV ",
      "on ka was likewise imprecise; removing it lowered BICc by 7.37 units ",
      "and additionally fixing ka to the Gaspar 2023 value of 0.82 1/h ",
      "lowered it by a further 9.48. Exposure metrics were computed for 104 ",
      "of the 106 patients (2 excluded for an in-hospital apixaban dose ",
      "reduction) and cross-checked against 1000 simulated individuals ",
      "sampled from the individual covariate parameters, with 90% ",
      "confidence intervals from nonparametric bootstrap (1000 replicates). ",
      "Limitations noted by the authors: covariates reflect only the first ",
      "reported value per encounter and so do not track the in-hospital ",
      "time course; apixaban was temporarily held in several patients for ",
      "procedures or nonadherence, and only the encounter with the most ",
      "samples was analysed per patient; metabolizer phenotypes and ",
      "CYP / P-gp genotypes were unavailable; sampling times were clinical ",
      "rather than protocol-driven; and neither the amiodarone dose nor its ",
      "metabolites were modelled. The authors' clinical conclusion is that ",
      "the 1.1- to 1.7-fold exposure increase does not warrant empiric ",
      "apixaban dose reduction with concomitant amiodarone."
    )
  )

  ini({
    # Structural parameters -- Kolowrat 2025 Table 2, FINAL model column.
    # Table 2 also prints a base model (no covariates: V/F 43.45 L, CL/F 1.32
    # L/h); per the replicate-author-structure policy only the final model is
    # extracted.

    # Absorption rate constant, FIXED. Table 2 prints '0.82 (Fixed)' with no
    # RSE and no confidence interval in both the base and final columns.
    # Results 3.2 gives the provenance and the reason: 'there was suboptimal
    # precision around the estimation of the IIV on K a (absorption rate
    # constant). Removing IIV on K a decreased BICc by 7.37 units, which was
    # further improved by fixing K a to 0.82 (dBICc = -9.48) [22]', where
    # reference [22] is Gaspar 2023 (OptimAT). That value is independently
    # confirmed against the sibling extraction: Gaspar_2023_apixaban.R carries
    # lka <- log(0.82) from its own Table 2 (ka = 0.82 1/h, RSE 11%).
    lka <- fixed(log(0.82)); label("Absorption rate constant ka (1/h), from Gaspar 2023")  # Kolowrat 2025 Table 2: Ka = 0.82 1/h (Fixed), inherited from Gaspar 2023 Table 2

    # Apparent oral clearance at the reference covariate values, i.e. age 77
    # years (the cohort median) and no amiodarone. This is the 1.5 L/h
    # intercept of the published final equation in the Table 2 note.
    lcl <- log(1.5); label("Apparent oral clearance CL/F (L/h) at AGE = 77 years without amiodarone")  # Kolowrat 2025 Table 2 final: CL/F = 1.5 L/h (RSE 8.99%, 95% CI 1.26-1.78)

    # Apparent volume of distribution. No covariate and no reference value --
    # V/F was tested and no covariate was retained on it. The Discussion notes
    # this estimate 'is within the range established by prior studies
    # (~25-53 L)'.
    lvc <- log(45.57); label("Apparent volume of distribution V/F (L)")  # Kolowrat 2025 Table 2 final: V/F = 45.57 L (RSE 10.6%, 95% CI 37.11-55.96)

    # Power exponent on (AGE / 77) for CL/F. Applied as (AGE / 77)^e_age_cl
    # per the Table 2 note, which is the log-linear continuous-covariate form
    # of Methods 2.3 written multiplicatively. NEGATIVE: apparent clearance
    # decreases with age.
    e_age_cl <- -1.52; label("Power exponent on (AGE / 77) for CL/F (unitless)")  # Kolowrat 2025 Table 2 final: beta_CL,age = -1.52 (RSE 25%, 95% CI -2.27 to -0.78)

    # Log-scale effect of concomitant amiodarone on CL/F, applied as
    # exp(e_amio_cl * CONMED_AMIO): the exponent collapses to 0 (factor 1)
    # without amiodarone and to -0.4 (factor exp(-0.4) = 0.670) with it. The
    # paper reports this coefficient already on the log scale, per the
    # categorical form of Methods 2.3 -- unlike the sibling
    # Morath_2025_apixaban.R, which reports the multiplicative factor 0.679
    # and therefore encodes log(0.679). Both files end up with the same
    # exp(e_amio_cl * CONMED_AMIO) shape in model().
    #
    # The encoding is confirmed arithmetically by the paper's own prose:
    # exp(-0.4) = 0.670 is a 33.0% reduction, and the Abstract / Results /
    # Discussion all state 'a 33% decrease (ranging from 12% to 48%) in
    # clearance'. The stated range is the confidence interval transformed the
    # same way: exp(-0.66) = 0.517 (48.3% reduction) and exp(-0.13) = 0.878
    # (12.2% reduction). A multiplicative (1 + beta) reading would instead
    # give a 40% reduction and would not reproduce the 12-48% range.
    e_amio_cl <- -0.4; label("Log of the multiplicative factor on CL/F for concomitant amiodarone")  # Kolowrat 2025 Table 2 final: beta_CL,amio = -0.4 (RSE 33.8%, 95% CI -0.66 to -0.13)

    # Interindividual variability -- Kolowrat 2025 Table 2 final, IIV rows.
    # IIV was supported on CL/F and V/F but NOT on ka (Results 3.2); the eta
    # on ka was removed before ka itself was fixed. The two etas are entered
    # as independent diagonal elements: the paper reports no CL/F-V/F
    # correlation or covariance, so none is imposed.
    #
    # SCALE CONVENTION. The Table 2 rows are headed 'IIV V (CV%)' and 'IIV CL
    # (CV%)' with values 52.66 and 62.23, while the adjacent '95% CI' cells
    # hold 0.35-0.69 and 0.47-0.69 -- i.e. the estimate is a percent CV and
    # the interval is on the omega (SD) scale, two different scales in one
    # row. The internal log-scale variance is therefore recovered by
    # inverting the log-normal CV transform, omega^2 = log(CV^2 + 1):
    #   V/F : omega = sqrt(log(0.5266^2 + 1)) = 0.49473, omega^2 = 0.244754
    #   CL/F: omega = sqrt(log(0.6223^2 + 1)) = 0.57213, omega^2 = 0.327329
    # This reading is confirmed against all four printed intervals (base and
    # final model, both parameters) using the log-scale Wald form
    # omega * exp(+/- 1.96 * RSE), which Monolix uses for variance
    # parameters. It reproduces every bound after rounding -- V/F final
    # 0.352-0.696 vs printed 0.35-0.69; CL/F final 0.474-0.691 vs printed
    # 0.47-0.69; V/F base 0.361-0.676 vs printed 0.36-0.67; CL/F base
    # 0.528-0.754 vs printed 0.53-0.75. Reading the CV% column as omega
    # directly instead misses 6 of those 8 bounds (e.g. CL/F base would give
    # 0.585-0.835 against a printed 0.53-0.75), so the convention is pinned.
    etalvc ~ 0.244754  # Kolowrat 2025 Table 2 final: IIV V = 52.66% CV (RSE 17.4%, omega CI 0.35-0.69); omega^2 = log(0.5266^2 + 1)
    etalcl ~ 0.327329  # Kolowrat 2025 Table 2 final: IIV CL = 62.23% CV (RSE 9.65%, omega CI 0.47-0.69); omega^2 = log(0.6223^2 + 1)

    # Residual unexplained variability -- combined proportional and additive,
    # selected over additive-only and proportional-only alternatives
    # (Methods 2.2, Results 3.2).
    #
    # The proportional row is headed '(CV%)' like the IIV rows but its value
    # 0.15 and interval 0.11-0.21 are plainly FRACTIONS, not percents: a
    # 0.15% proportional error alongside a 21 ng/mL additive term is
    # implausible, and the log-scale Wald interval on 0.15 at RSE 18.3% is
    # 0.105-0.215, matching the printed 0.11-0.21. So propSd = 0.15 (15%).
    propSd <- 0.15; label("Proportional residual error (fraction)")  # Kolowrat 2025 Table 2 final: proportional error = 0.15 (RSE 18.3%, 95% CI 0.11-0.21)

    # The additive row carries no unit in Table 2; it is in the units of the
    # observation, ng/mL. Its magnitude (21.28) sits well above the 5 ng/mL
    # assay lower limit of quantification, consistent with real-world
    # salvage-sample data whose recorded dose and sample times are uncertain.
    addSd <- 21.28; label("Additive residual error (ng/mL)")  # Kolowrat 2025 Table 2 final: additive error = 21.28 (RSE 25.4%, 95% CI 13.21-34.26)
  })

  model({
    # 1. Individual parameters.
    #    CL/F carries both retained covariate effects and the CL eta. The
    #    published final equation is 'CL/F_i = 1.5 x (Age/77)^beta_age x
    #    e^(beta_amio x 1_{CAT=1}) x e^(IIV_CL/F_i)' (Kolowrat 2025 Table 2
    #    note), with log-normal IIV applied multiplicatively via exp(etalcl).
    cl <- exp(lcl + etalcl) * (AGE / 77)^e_age_cl * exp(e_amio_cl * CONMED_AMIO)

    # No covariate was retained on V/F, but IIV was (Results 3.2).
    vc <- exp(lvc + etalvc)

    # ka is fixed and carries no IIV -- its eta was removed before the value
    # itself was fixed to the Gaspar 2023 estimate (Results 3.2).
    ka <- exp(lka)

    # 2. Micro-constants.
    kel <- cl / vc

    # 3. One-compartment ODE system with first-order oral absorption and no
    #    lag time (Methods 2.2: 'first-order absorption and linear
    #    elimination without lag time'; Results 3.2 confirms no lag time in
    #    the final model). Bioavailability F is not identifiable from
    #    oral-only data and is absorbed into the apparent parameters CL/F and
    #    V/F, so no f(depot) term is applied.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Observation and error.
    #    central is in mg and vc in L, so central / vc is mg/L; multiply by
    #    1000 to report ng/mL (1 mg/L = 1 ug/mL = 1000 ng/mL), matching the
    #    paper's reporting units for concentrations, the 5 ng/mL assay lower
    #    limit of quantification, and the ng/mL Cmax / Cmin and ng*h/mL AUC
    #    of Table 3. Same convention as the sibling apixaban models
    #    Gaspar_2023_apixaban.R and Morath_2025_apixaban.R.
    Cc <- 1000 * central / vc

    # Combined proportional plus additive residual error. The source used
    # Monolix, which offers two combined parameterisations -- combined1
    # (SD = a + b*f) and combined2 (SD = sqrt(a^2 + (b*f)^2)) -- and the
    # paper names only 'combined proportional and additive error' without
    # saying which. nlmixr2's default add() + prop() is the combined2 form;
    # it is used here rather than combined1() because the source makes no
    # explicit combined1 declaration. See the vignette Errata.
    Cc ~ add(addSd) + prop(propSd)
  })
}

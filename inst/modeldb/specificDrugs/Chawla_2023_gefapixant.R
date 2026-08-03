Chawla_2023_gefapixant <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for the P2X3-receptor antagonist gefapixant in healthy volunteers and adults with refractory or unexplained chronic cough (Chawla 2023)"
  reference <- "Chawla A, Largajolli A, Hussain A, et al. Population pharmacokinetic analysis of the P2X3-receptor antagonist gefapixant. CPT Pharmacometrics Syst Pharmacol. 2023;12(8):1107-1118. doi:10.1002/psp4.12978"
  vignette <- "Chawla_2023_gefapixant"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate (eGFR)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Creatinine-based estimated GFR, BSA-normalized. Power effect on CL/F centered on the",
        "data median of 87.2 mL/min/1.73 m^2 (Chawla 2023 Table 1 footnote b). Added to CL/F as",
        "part of base-model development (not via the stepwise covariate method) because gefapixant",
        "is primarily renally eliminated, and retained through the final model. A sensitivity",
        "analysis substituting a power relationship of creatinine clearance for eGFR slightly",
        "worsened the objective function value and was not retained. Cohort median 86.9",
        "(range 13-243) mL/min/1.73 m^2 per Table S3. Individuals with end-stage renal disease or",
        "on hemodialysis were excluded, so the eGFR-CL/F relationship is not validated below the",
        "severe-RI range."
      ),
      source_name        = "eGFR"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effects on both CL/F and Vc/F centered on the data median of 74 kg",
        "(Chawla 2023 Table 1 footnote b). Cohort median 74 kg (range 35-159) per Table S3.",
        "Body mass index was deliberately not evaluated because of its high correlation with",
        "body weight. Fixed allometric scaling on CL/F and the volume parameters was tested as a",
        "sensitivity analysis and gave no improvement over these estimated exponents."
      ),
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effects on both CL/F and Vc/F centered on the data median of 59 years",
        "(Chawla 2023 Table 1 footnote b). Cohort median 59 years (range 18-89) per Table S3.",
        "The age effect on Vc/F was the least precisely estimated fixed effect in the final",
        "model (RSE 37.1%)."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) -- inverted relative to the canonical SEXF reference",
      notes              = paste(
        "Chawla 2023 uses a MALE indicator with FEMALE as the reference category (Table 1",
        "footnote b: 'Categorical covariate references: sex = female'), which is the inverse of",
        "the canonical SEXF orientation (1 = female, reference 0 = male). To store the column",
        "under the canonical SEXF name while reproducing the paper's published CL/F = 10.3 L/h",
        "and Vc/F = 101 L verbatim (both of which are female-reference typical values), the",
        "effects are applied in model() through sex_male <- 1 - SEXF, so SEXF = 1 gives a factor",
        "of exactly 1 and SEXF = 0 gives the paper's male-vs-female fractional change. Same",
        "pattern as Bajaj_2017_nivolumab.R. Cohort 71.3% female per Table S3."
      ),
      source_name        = "SEX (male indicator)"
    ),
    FED = list(
      description        = "Fed-vs-fasted state at the time of dosing (1 = fed, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "General fed-vs-fasted dose-record flag, not a high-fat-meal challenge, so FED applies",
        "rather than FED_HIGHFAT. Fractional-change effect on Ka, and applied ONLY to records",
        "carrying the F02 formulation (see FORM_GEF_F02): dedicated phase I relative-bioavailability",
        "studies had already shown that fed status affects gefapixant exposure for F02 but not for",
        "F04, so food effects were tested for F02 only. Food effects on absorption lag time and on",
        "bioavailability were also tested during base-model development but were not retained",
        "(nonsignificant objective-function drop or poor precision)."
      ),
      source_name        = "FED"
    ),
    FORM_GEF_F02 = list(
      description        = "Gefapixant F02 (wet-granulation, citric-acid-containing, film-coated immediate-release tablet) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (F04 / F04A gefapixant-citrate formulations)",
      notes              = paste(
        "1 = F02, the earlier wet-granulation immediate-release tablet containing citric acid as",
        "an acidulant (7.5, 20, and 50 mg strengths), used in several phase I studies and the",
        "phase IIb chronic-cough study. 0 = F04 (20A film coating, two phase I studies) or F04A",
        "(03K film coating, COUGH-1 and COUGH-2), both of which contain gefapixant citrate as the",
        "active ingredient. This indicator carries NO main effect on any PK parameter in the final",
        "model -- no relative-bioavailability difference between F02 and F04/F04A was retained. It",
        "exists solely to gate the FED effect on Ka to F02 records. Formulation assignment by study",
        "is tabulated in Chawla 2023 Table S2. The marketed F04B formulation (F04A without citric",
        "acid) is bioequivalent to F04A, so F04B records take FORM_GEF_F02 = 0."
      ),
      source_name        = "FORM"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor (omeprazole) use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested during base-model development on absorption rate constant, absorption lag time,",
        "and relative bioavailability for the F02 formulation only (dedicated phase I studies had",
        "shown no PPI effect on F04). No PPI relationship reached the retention criteria, so PPI",
        "use appears nowhere in the final model. Chawla 2023 Results, 'Model development'."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was screened on CL/F and Vc/F in the stepwise covariate method and was statistically",
        "significant on CL/F in the integrated phase I-III analysis, but the signal was driven",
        "entirely by the 'multiple' race category (effect size 0.19 vs ~0.05 for every other",
        "category) which comprised only ~5% of the target population. 'Multiple' was therefore",
        "merged into 'other' and the race effect was removed from the final model as not clinically",
        "relevant. Chawla 2023 Results and Discussion; no retained point estimate exists."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened with the other race indicators; race removed from the final model. See RACE_WHITE note."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened with the other race indicators; race removed from the final model. See RACE_WHITE note."
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator (includes American Indian or Alaskan Native, multiple or other, and Native Hawaiian or Pacific Islander)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened with the other race indicators; the 'multiple' category was merged into 'other'",
        "before race was dropped from the final model. See RACE_WHITE note."
      )
    ),
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Ethnicity was screened on CL/F and Vc/F in the stepwise covariate method and was not",
        "retained in the final model. Cohort 17.0% Hispanic per Chawla 2023 Table S3."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1618,
    n_studies      = 9,
    age_range      = "18-89 years",
    age_median     = "59 years",
    weight_range   = "35-159 kg",
    weight_median  = "74 kg",
    sex_female_pct = 71.3,
    race_ethnicity = c(White = 78.5, Asian = 8.3, Black = 3.6, Other = 9.5),
    disease_state  = "Healthy volunteers (n = 121) and adults with refractory or unexplained chronic cough (n = 1497)",
    dose_range     = "7.5-150 mg oral b.i.d. (phase I, single and multiple dose); up to 50 mg b.i.d. (phase IIb); 15 and 45 mg b.i.d. (phase III COUGH-1 and COUGH-2)",
    regions        = "Global; United Kingdom and United States (phase IIb), Japan (one phase I study), multinational (phase III)",
    renal_function = paste(
      "eGFR median 86.9 mL/min/1.73 m^2 (range 13-243). Among the 1555 chronic-cough participants,",
      "664 had normal renal function, 817 mild renal impairment, and 74 moderate renal impairment;",
      "no phase II/III participant had severe renal impairment, which is represented only by a",
      "dedicated phase I renal-impairment study. Individuals with end-stage renal disease and",
      "individuals on hemodialysis were excluded."
    ),
    formulations   = "F02, F04, and F04A; the earlier F01 formulation was excluded",
    notes          = paste(
      "Baseline demographics from Chawla 2023 Table S3; per-study participant counts and",
      "formulation / food / PPI assignments from Table S2; study designs from Table S1.",
      "1677 participants were included in the data set and 1618 had evaluable PK data",
      "(8886 measurable concentrations) after exclusions for inconsistent dosing times, an extra",
      "dose within 3 days before sampling, measurable pre-first-dose concentrations, missing dose",
      "times, and |CWRES| > 5 outliers on the newly added phase II/III data."
    )
  )

  ini({
    # Structural parameters. CL/F and Vc/F are typical values at the reference
    # covariate set: eGFR 87.2 mL/min/1.73 m^2, age 59 y, body weight 74 kg, female.
    # Ka is the typical value under the fasted-state reference.
    lka   <- log(2.25);  label("Absorption rate constant (Ka, 1/h), fasted reference")                              # Chawla 2023 Table 1: Ka = 2.25 (95% CI 1.86, 2.85), RSE 9.9%
    lcl   <- log(10.3);  label("Apparent clearance (CL/F, L/h) at reference eGFR, age, weight, and female sex")     # Chawla 2023 Table 1: CL/F = 10.3 (10.1, 10.5), RSE 1.1%
    lvc   <- log(101);   label("Apparent central volume of distribution (Vc/F, L) at reference age, weight, sex")   # Chawla 2023 Table 1: Vc/F = 101 (96.9, 104), RSE 1.8%
    lq    <- log(3.51);  label("Apparent intercompartmental clearance (Q/F, L/h)")                                  # Chawla 2023 Table 1: Q/F = 3.51 (2.7, 4.42), RSE 12.8%
    lvp   <- log(32.8);  label("Apparent peripheral volume of distribution (Vp/F, L)")                              # Chawla 2023 Table 1: Vp/F = 32.8 (26.9, 46.6), RSE 12.7%
    ltlag <- log(0.432); label("Absorption lag time (ALAG, h)")                                                     # Chawla 2023 Table 1: ALAG = 0.432 (0.415, 0.445), RSE 1.7%

    # Covariate effects on CL/F.
    # Continuous covariates enter as median-centered power functions and categorical
    # covariates as fractional changes, per the two Chawla 2023 Methods equations:
    #   P_i = P_TV * (cov_i / median(cov))^theta * exp(eta_i)
    #   P_i = P_TV * (1 + theta_cov,c * I_cov,c,i) * exp(eta_i)
    e_crcl_cl <- 0.375;  label("Power exponent of eGFR on CL/F, centered at 87.2 mL/min/1.73 m^2 (unitless)")       # Chawla 2023 Table 1: Cl_eGFR = 0.375 (0.317, 0.429), RSE 8.4%
    e_age_cl  <- -0.229; label("Power exponent of age on CL/F, centered at 59 years (unitless)")                    # Chawla 2023 Table 1: CL_AGE = -0.229 (-0.284, -0.171), RSE 13.3%
    e_wt_cl   <- 0.35;   label("Power exponent of body weight on CL/F, centered at 74 kg (unitless)")               # Chawla 2023 Table 1: CL_BW = 0.35 (0.27, 0.43), RSE 11.3%
    e_sex_cl  <- 0.0931; label("Fractional change in CL/F for male sex vs female reference (applied as (1 - SEXF))") # Chawla 2023 Table 1: CL_SEX = 0.0931 (0.0545, 0.133), RSE 22.7%

    # Covariate effects on Vc/F.
    e_age_vc  <- 0.0911; label("Power exponent of age on Vc/F, centered at 59 years (unitless)")                    # Chawla 2023 Table 1: Vc_AGE = 0.0911 (0.0239, 0.162), RSE 37.1%
    e_wt_vc   <- 0.541;  label("Power exponent of body weight on Vc/F, centered at 74 kg (unitless)")               # Chawla 2023 Table 1: Vc_BW = 0.541 (0.452, 0.627), RSE 8.4%
    e_sex_vc  <- 0.181;  label("Fractional change in Vc/F for male sex vs female reference (applied as (1 - SEXF))") # Chawla 2023 Table 1: Vc_SEX = 0.181 (0.138, 0.234), RSE 13.6%

    # Covariate effect on Ka. Applies only to F02-formulation records; dedicated
    # phase I studies showed no food effect for the F04 formulations.
    e_fed_ka  <- -0.594; label("Fractional change in Ka for the fed state, F02 formulation only (vs fasted reference)") # Chawla 2023 Table 1: Ka_FEED = -0.594 (-0.675, -0.496), RSE 7.7%

    # Interindividual variability. Chawla 2023 Table 1 reports these as variances
    # (omega^2) with footnote d giving %CV = sqrt(exp(omega^2) - 1) * 100.
    # CL/F and Vc/F IIV are correlated; the reported CL/F-Vc/F term is a covariance.
    etalcl + etalvc ~ c(0.0708,
                        0.0176, 0.0161)  # Chawla 2023 Table 1: omega^2_CL = 0.0708 (27.1% CV, RSE 6.7%); CL-Vc covariance = 0.0176 (RSE 22.9%); omega^2_Vc = 0.0161 (12.7% CV, RSE 41.8%)
    etalka          ~ 0.551              # Chawla 2023 Table 1: omega^2_Ka = 0.551 (85.7% CV, RSE 18.7%); informed by phase I data only

    # Residual error. Chawla 2023 Table 1 reports these as standard deviations
    # (sigma), not variances: the additive term carries ng/mL units and the
    # proportional term corresponds to 30.3% CV.
    propSd <- 0.303; label("Proportional residual error (fraction)")  # Chawla 2023 Table 1: sigma_PROP = 0.303 (0.29, 0.316), RSE 2.3%
    addSd  <- 3.04;  label("Additive residual error (ng/mL)")         # Chawla 2023 Table 1: sigma_ADD = 3.04 (2.3, 3.81) ng/mL, RSE 15.4%
  })
  model({
    # 1. Derived covariate terms.
    # Chawla 2023 uses a male indicator with female as the reference category, the
    # inverse of the canonical SEXF orientation; (1 - SEXF) reproduces the paper's
    # male = 1 column while keeping the published female-reference CL/F and Vc/F.
    sex_male <- 1 - SEXF

    # Food slows gefapixant absorption for the F02 formulation only.
    fed_ka <- 1 + e_fed_ka * FED * FORM_GEF_F02

    # 2. Individual parameters.
    ka   <- exp(lka + etalka) * fed_ka
    tlag <- exp(ltlag)
    cl   <- exp(lcl + etalcl) *
      (CRCL / 87.2)^e_crcl_cl *
      (AGE / 59)^e_age_cl *
      (WT / 74)^e_wt_cl *
      (1 + e_sex_cl * sex_male)
    vc   <- exp(lvc + etalvc) *
      (AGE / 59)^e_age_vc *
      (WT / 74)^e_wt_vc *
      (1 + e_sex_vc * sex_male)
    q    <- exp(lq)
    vp   <- exp(lvp)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Absorption lag time.
    alag(depot) <- tlag

    # 6. Observation and error.
    # Dose in mg and volume in L give mg/L; * 1000 converts to ng/mL to match the
    # additive residual error units and the reported concentration scale.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

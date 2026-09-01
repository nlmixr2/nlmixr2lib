Nielsen_2015_vigabatrin <- function() {
  description <- paste(
    "Population dose-response model for daily seizure count in adult and pediatric patients with",
    "refractory complex partial seizures (rCPS) receiving adjunctive vigabatrin, from Nielsen 2015",
    "(pooled analysis of two adult and three pediatric randomized controlled trials; 621 patients,",
    "112,168 daily seizure records). Daily counts follow a negative-binomial distribution with mean",
    "seizure rate lambda and overdispersion ovdp, so that variance = lambda * (1 + ovdp * lambda).",
    "The mean rate is multiplicative in three terms: a baseline rate that rises with decreasing age;",
    "an asymptotic non-drug time effect switched on at the first randomized dose (astime, approached",
    "with rate constant exp(lktime), half-life about 103 days); and a quadratic drug effect in",
    "normalized dosage. Normalized dosage (dosenorm) is the actual total daily vigabatrin dosage",
    "rescaled to a 60-kg-equivalent exposure surrogate via (WT/60)^-0.608 -- an in-model stand-in for",
    "systemic exposure, estimated because pediatric PK data were unavailable. There is no PK ODE and",
    "no ODE state at all: exposure enters as the covariate DOSE_VGB_MGD. Overdispersion and the",
    "baseline-rate IIV magnitude are study-specific; the Box-Cox transform of the baseline-rate random",
    "effect applies to the adult cohorts only. rxode2 exposes no negative-binomial endpoint, so the",
    "observation is declared as a Poisson likelihood and the overdispersion is exposed as the model",
    "variable ovdp; see the vignette for negative-binomial sampling via rnbinom(mu, size = 1/ovdp).",
    sep = " "
  )
  reference <- paste(
    "Nielsen JC, Hutmacher MM, Wesche DL, Tolbert D, Patel M, Kowalski KG.",
    "Population dose-response analysis of daily seizure count following vigabatrin therapy",
    "in adult and pediatric patients with refractory complex partial seizures.",
    "J Clin Pharmacol. 2015;55(1):81-92.",
    "doi:10.1002/jcph.378.",
    sep = " "
  )
  vignette <- "Nielsen_2015_vigabatrin"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "seizures/day (the model output is a mean daily seizure rate, not a drug concentration)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Patient age at study entry.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters the log baseline seizure rate as",
        "e_age_rbase * (log(AGE) - log(24)), i.e. a power model on age centred at the pooled median",
        "of 24 years (Nielsen 2015 Table 2; Equation 10b). The coefficient is negative, so younger",
        "patients have a HIGHER baseline seizure rate -- the paper's 'Baseline seizure rate increased",
        "with decreasing age'. Age replaced the three per-study shifts on lambda during covariate",
        "analysis because pediatric-only and adult-only enrolment made study and age collinear",
        "(Nielsen 2015 Results, Covariate Model). Pooled range 3-63 years."
      ),
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline weight). Sole covariate on normalized dosage:",
        "dosenorm = DOSE_VGB_MGD * (WT/60)^e_wt_dosenorm with e_wt_dosenorm = -0.608 and a 60 kg",
        "reference (Nielsen 2015 Equation 10e). Creatinine clearance was also tested on normalized",
        "dosage but was small and imprecisely estimated and was dropped (Results, Covariate Model);",
        "the paper attributes this to body size being the dominant determinant of glomerular",
        "filtration in a cohort without significant renal impairment. Pooled range 12-136 kg,",
        "median 62 kg (Table 2)."
      ),
      source_name        = "WT"
    ),
    DOSE_VGB_MGD = list(
      description        = "Total daily vigabatrin dosage.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-record (time-varying) covariate: dosage was titrated in all five studies and given as a",
        "twice-daily regimen, so this column carries the SUM across the day and is updated as the",
        "patient escalates. Set to 0 mg/day during the baseline run-in and for placebo subjects, which",
        "makes dosenorm = 0 and the drug-effect term exp(0) = 1. Studies 118/192/221 dosed pediatric",
        "patients by weight band; studies 24 and 25 used flat adult dosages of 1, 3 or 6 g/day",
        "(Nielsen 2015 Table 1). NOTE the Study 118 dosing-weight cap: participants weighing more than",
        "60 kg received their mg/kg dosage as if they weighed 60 kg, a correction the authors applied",
        "to the analysis dataset before the final model run (Results, Covariate Model)."
      ),
      source_name        = "DOSE"
    ),
    STUDY_118 = list(
      description        = "1 = subject enrolled in pediatric Study 118; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled adult studies 24 and 25 are the reference)",
      notes              = paste(
        "Time-fixed per subject. Shifts BOTH the log overdispersion (e_study_118_ovdp) and the",
        "asymptotic time effect (e_study_118_astime), and scales the baseline-rate IIV SD",
        "(e_study_118_sd_rbase). n = 125; placebo or vigabatrin 20/60/100 mg/kg/day (Table 1).",
        "Study 118 is the study whose AS shift moved when active-treatment data were added and whose",
        "maximum time effect is SMALLER than adults' -- the paper flags this as unexplained",
        "(Results, Vigabatrin Treatment Effects; Discussion)."
      ),
      source_name        = "ST118"
    ),
    STUDY_192 = list(
      description        = "1 = subject enrolled in pediatric Study 192; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled adult studies 24 and 25 are the reference)",
      notes              = paste(
        "Time-fixed per subject. Shifts the log overdispersion (e_study_192_ovdp) and the asymptotic",
        "time effect (e_study_192_astime), and scales the baseline-rate IIV SD (e_study_192_sd_rbase).",
        "n = 55; placebo or weight-banded vigabatrin 0.5-4 g/day (Table 1). Study 192's placebo arm is",
        "the one group whose observed seizure rate fell outside the final model's 90% prediction",
        "interval in the VPC (Results, Model Diagnostics; Figure 2)."
      ),
      source_name        = "ST192"
    ),
    STUDY_221 = list(
      description        = "1 = subject enrolled in pediatric Study 221; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled adult studies 24 and 25 are the reference)",
      notes              = paste(
        "Time-fixed per subject. Shifts the log overdispersion (e_study_221_ovdp) and the asymptotic",
        "time effect (e_study_221_astime), and scales the baseline-rate IIV SD (e_study_221_sd_rbase).",
        "n = 85; placebo or weight-banded vigabatrin 0.5-4 g/day (Table 1)."
      ),
      source_name        = "ST221"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Baseline creatinine clearance (mL/min).",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened as a covariate on normalized dosage in the full model but NOT retained: the effect",
        "was 'small and imprecisely estimated' and removing it (with age on AS and age on b) raised",
        "the objective function by only 0.435 units (Nielsen 2015 Results, Covariate Model). No point",
        "estimate is published, so the effect cannot be reconstructed. Pooled median 99 mL/min",
        "(range 21-273; Table 2)."
      )
    ),
    CONMED_AED = list(
      description = "Concomitant antiepileptic drug use (carbamazepine, valproic acid, lamotrigine, hydantoins, gabapentin, barbiturates, benzodiazepines, methsuximide).",
      units       = "(binary, one indicator per drug class)",
      type        = "binary",
      notes       = paste(
        "Screened in the second stage of covariate analysis as a shift in ln(lambda), as eight",
        "additional parameters. Not retained: inclusion reduced the objective function by only 5.4",
        "units for 8 parameters, so concomitant AEDs were not predictive of baseline seizure rate",
        "(Nielsen 2015 Results, Covariate Model). Per-class prevalences are in Table 2; no effect",
        "estimates are published."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 621L,
    n_studies      = 5L,
    age_range      = "3-63 years (median 24; pediatric median 11, adult median 33)",
    weight_range   = "12-136 kg (median 62; pediatric median 44, adult median 70)",
    sex_female_pct = 53.6,
    race_ethnicity = c(White = 91.5, Black = 5.5, Asian = 0.6, Other = 2.4),
    disease_state  = paste(
      "Refractory complex partial seizures (rCPS), with or without secondary generalization, on a",
      "stable regimen of one or two background antiepileptic drugs. Patients with generalized",
      "epilepsy, progressive neurological disorders, treatable causes of seizures or non-epileptic",
      "seizures were excluded (Nielsen 2015 Table 1)."
    ),
    dose_range     = paste(
      "Adults (studies 24, 25): placebo or vigabatrin 1, 3 or 6 g/day. Pediatrics (study 118):",
      "placebo or 20, 60 or 100 mg/kg/day, with the dosing weight capped at 60 kg. Pediatrics",
      "(studies 192, 221): placebo or weight-banded 0.5-1.5 g/day (10-15 kg), 0.5-2.0 g/day",
      "(16-30 kg), 1.0-3.0 g/day (31-50 kg), 1.0-4.0 g/day (over 50 kg). Twice-daily and titrated",
      "in all five studies."
    ),
    regions        = "not reported",
    biomarkers     = paste(
      "Daily seizure count from patient diaries. 41,282 daily records from 356 adults and 70,886",
      "daily records from 265 pediatric patients (112,168 records total). Study phases: baseline",
      "run-in 6-10 weeks (pediatric) or 10 weeks (adult); dosage titration 6-10 weeks (pediatric) or",
      "4-6 weeks (adult); maintenance 7-8 weeks (pediatric) or 12 weeks (adult)."
    ),
    notes          = paste(
      "Demographics from Nielsen 2015 Table 2 (pooled column); study designs and enrolment from",
      "Table 1. The three pediatric studies were suspended before completing planned enrolment for",
      "administrative reasons and were individually underpowered for dose response, which is the",
      "stated motivation for the pooled analysis. Model fitted in NONMEM 7.2 with the Laplace method;",
      "final objective function 195069.298, condition number 149 (Table 3)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # All values are final estimates from Nielsen 2015 Table 3 'Final Model
    # Parameter Estimates' (page 87), which reflects the re-run on the
    # dataset corrected for the Study 118 60-kg dosing-weight cap (Results,
    # Covariate Model: 'The parameter estimates and model diagnostics
    # reported here reflect the updated parameter estimates in the revised
    # data set'). Standard errors are quoted for provenance only.
    #
    # Table 3 rows prefixed 'LN' are on the natural-log scale in the source
    # and are carried on the log scale here unchanged.
    # ------------------------------------------------------------------

    # -- Baseline daily seizure rate (Equation 10b) --------------------
    lrbase <- -1.070
    label("Typical log baseline daily seizure rate for a 24-year-old adult (log(seizures/day))")
    # Table 3 'LN lambda - Adults' = -1.070 (SE 0.0394); exp(-1.07) = 0.343 seizures/day,
    # matching the paper's 'typical predicted mean daily baseline seizure rate ... was 0.343
    # for a 24-year-old patient' (Results, Covariate Model).

    e_age_rbase <- -0.4180
    label("Coefficient on (log(AGE) - log(24)) for the log baseline seizure rate (unitless)")
    # Table 3 'Age on LN lambda' = -0.4180 (SE 0.0592). Negative, so baseline seizure rate
    # increases with DECREASING age (Abstract; Discussion).

    # -- Negative-binomial overdispersion (Equation 1) -----------------
    # variance = lambda * (1 + ovdp * lambda); ovdp -> 0 collapses to Poisson.
    lovdp <- 0.0924
    label("Typical log negative-binomial overdispersion for the adult studies (log unitless)")
    # Table 3 'LN OVDP - Adults' = 0.0924 (SE 0.0195)

    e_study_118_ovdp <- -0.6710
    label("Shift on log overdispersion for pediatric Study 118 (unitless)")
    # Table 3 'D LN OVDP - Study 118' = -0.6710 (SE 0.0412)

    e_study_192_ovdp <- -0.9400
    label("Shift on log overdispersion for pediatric Study 192 (unitless)")
    # Table 3 'D LN OVDP - Study 192' = -0.9400 (SE 0.0743)

    e_study_221_ovdp <- -1.4200
    label("Shift on log overdispersion for pediatric Study 221 (unitless)")
    # Table 3 'D LN OVDP - Study 221' = -1.4200 (SE 0.0730)

    # -- Asymptotic non-drug time effect (Equations 10c, 10d) ----------
    astime <- -0.1560
    label("Asymptotic maximum time effect on log seizure rate, adult reference (unitless, log fractional change)")
    # Table 3 'AS - Adults' = -0.1560 (SE 0.2430)

    e_study_118_astime <- 0.1200
    label("Shift on the asymptotic time effect for pediatric Study 118 (unitless)")
    # Table 3 'D AS - Study 118' = 0.1200 (SE 0.4060). Positive: Study 118 has a SMALLER
    # maximum decrease in seizure rate from time effects than adults (Discussion).

    e_study_192_astime <- -0.5840
    label("Shift on the asymptotic time effect for pediatric Study 192 (unitless)")
    # Table 3 'D AS - Study 192' = -0.5840 (SE 0.3780)

    e_study_221_astime <- -0.3740
    label("Shift on the asymptotic time effect for pediatric Study 221 (unitless)")
    # Table 3 'D AS - Study 221' = -0.3740 (SE 0.3640)

    lktime <- -5.000
    label("Log first-order rate constant for approach of the time effect to astime (log(1/day))")
    # Table 3 'LN k' = -5.000 (SE 0.2780). exp(-5) = 0.006738 /day, i.e. a half-life of
    # log(2)/exp(-5) = 102.9 days, matching the paper's 'half-life of approximately 100 days'
    # (Results, Covariate Model).

    # -- Quadratic drug effect on normalized dosage (Equation 10f) -----
    # fdrug = exp(aquad * x^2 - blin * x) with x = dosenorm / 3000. Both
    # coefficients are estimated on the log scale so they stay positive
    # ('exponential error models were assumed for parameters that were
    # constrained to be positive', Methods, Base Model).
    laquad <- -1.2900
    label("Log quadratic coefficient of the drug-effect polynomial in (dosenorm/3000) (log unitless)")
    # Table 3 'LN a' = -1.2900 (SE 0.2280); exp(-1.29) = 0.2753

    lblin <- -0.1220
    label("Log linear coefficient of the drug-effect polynomial in (dosenorm/3000) (log unitless)")
    # Table 3 'LN b' = -0.1220 (SE 0.1180); exp(-0.122) = 0.8851

    # -- Normalized dosage (Equation 10e) ------------------------------
    e_wt_dosenorm <- -0.6080
    label("Power exponent on (WT/60) for normalized dosage (unitless)")
    # Table 3 'WT on DNORM' = -0.6080 (SE 0.1740). The printed Equation 10e drops the minus
    # sign; the NEGATIVE exponent from Table 3 is the one that reproduces the paper's own
    # worked equivalence that a 20-kg patient on 1.54 g/day, a 60-kg patient on 3 g/day and a
    # 100-kg patient on 4.09 g/day share the same normalized dosage (Results, Covariate Model)
    # -- giving 3003, 3000 and 2998 mg respectively.

    # -- Box-Cox shape on the baseline-rate random effect (Equation 9) --
    bc_shape <- 0.7450
    label("Box-Cox shape parameter for the baseline-rate random effect, adult cohorts only (unitless)")
    # Table 3 'SHAPE - Adults' = 0.7450 (SE 0.1040). Applied to the ADULT data only: 'The
    # Box-Cox transformation in the random effects for lambda was not required for the
    # pediatric data. Therefore, the shape parameter was not estimated' (Results, Time Effects).

    # -- Baseline-rate IIV magnitude, log-SD scale ---------------------
    # Table 3 footnote: 'LNv = interindividual variability expressed as the
    # natural log of the standard deviation'. The source makes this SD
    # study-dependent, which ini() cannot express as a covariate-dependent
    # OMEGA, so etalrbase below is a standard normal that model() rescales.
    lsd_rbase <- -0.3190
    label("Log SD of the baseline log-rate random effect, adult reference (log unitless)")
    # Table 3 'LN v - lambda - Adults' = -0.3190 (SE 0.0552); SD = exp(-0.319) = 0.7269

    e_study_118_sd_rbase <- 0.5860
    label("Shift on the log SD of the baseline log-rate random effect for Study 118 (unitless)")
    # Table 3 'D LN v - lambda - Study 118' = 0.5860 (SE 0.0867)

    e_study_192_sd_rbase <- 0.8860
    label("Shift on the log SD of the baseline log-rate random effect for Study 192 (unitless)")
    # Table 3 'D LN v - lambda - Study 192' = 0.8860 (SE 0.1220)

    e_study_221_sd_rbase <- 0.3580
    label("Shift on the log SD of the baseline log-rate random effect for Study 221 (unitless)")
    # Table 3 'D LN v - lambda - Study 221' = 0.3580 (SE 0.0989)

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 3 reports each IIV as the natural
    # log of the STANDARD DEVIATION, so variances below are exp(2 * LNv).
    # ------------------------------------------------------------------

    # Standard normal by construction: the study-dependent SD above is
    # applied inside model(). Fixed at 1 so it is not re-estimated.
    etalrbase ~ fix(1)

    etaastime ~ 4.0068284
    # Table 3 'LN v - AS' = 0.6940 (SE 0.2130); SD = exp(0.694) = 2.0017, var = exp(2*0.694)

    etalblin ~ 0.7603321
    # Table 3 'LN v - b' = -0.1370 (SE 0.0699); SD = exp(-0.137) = 0.8720, var = exp(2*-0.137)
  })

  model({
    # ------------------------------------------------------------------
    # Nielsen 2015 Equations 10a-10f. There is no ODE state: the paper's
    # model is a NONMEM $PRED-style likelihood on daily seizure counts, and
    # drug exposure enters through the covariate DOSE_VGB_MGD rather than a
    # PK model (pediatric PK data were unavailable; Introduction).
    # ------------------------------------------------------------------

    # -- Cohort stratification -----------------------------------------
    # Studies 118, 192 and 221 enrolled only children; studies 24 and 25
    # enrolled only adults (Table 1), so the three indicators sum to the
    # pediatric flag and its complement is the adult flag.
    peds <- STUDY_118 + STUDY_192 + STUDY_221

    # -- Baseline daily seizure rate (Equations 9, 10b) ----------------
    # Study-dependent SD of the baseline-rate random effect.
    sd_rbase <- exp(lsd_rbase +
                      e_study_118_sd_rbase * STUDY_118 +
                      e_study_192_sd_rbase * STUDY_192 +
                      e_study_221_sd_rbase * STUDY_221)
    eta_rbase <- etalrbase * sd_rbase

    # Box-Cox transform of the baseline-rate random effect (Equation 9),
    # applied to the ADULT cohorts only; the pediatric branch keeps the
    # untransformed eta (Results, Time Effects).
    eta_rbase_bc <- (exp(eta_rbase)^bc_shape - 1) / bc_shape
    eta_rbase_eff <- (1 - peds) * eta_rbase_bc + peds * eta_rbase

    rbase <- exp(lrbase + e_age_rbase * (log(AGE) - log(24)) + eta_rbase_eff)

    # -- Non-drug time effect (Equations 10c, 10d) ---------------------
    astime_i <- astime +
      e_study_118_astime * STUDY_118 +
      e_study_192_astime * STUDY_192 +
      e_study_221_astime * STUDY_221 +
      etaastime
    ktime <- exp(lktime)

    # I(DAY >= 1) from Equations 3, 4 and 10d: the time effect switches on
    # at the first randomized dose. 'It was assumed that seizure rate was
    # constant during the baseline period, and changes in seizure rate due
    # to time effects were initiated on Day 1' (Methods, Base Model). `t`
    # is time in days relative to that first randomized dose, so baseline
    # run-in records carry t < 1.
    iday1 <- (t >= 1)
    ftime <- exp(astime_i * (1 - exp(-ktime * t)) * iday1)

    # -- Normalized dosage and quadratic drug effect (Eq. 10e, 10f) ----
    # dosenorm is the paper's DNORM: the actual total daily dosage rescaled
    # to a 60-kg-equivalent exposure surrogate. On placebo and during the
    # baseline run-in DOSE_VGB_MGD is 0, so dosenorm is 0 and fdrug is 1.
    dosenorm <- DOSE_VGB_MGD * (WT / 60)^e_wt_dosenorm
    dscaled <- dosenorm / 3000
    aquad <- exp(laquad)
    blin <- exp(lblin + etalblin)
    fdrug <- exp(aquad * dscaled^2 - blin * dscaled)

    # -- Mean seizure rate and overdispersion (Equations 1, 10a) -------
    lambda <- rbase * ftime * fdrug
    ovdp <- exp(lovdp +
                  e_study_118_ovdp * STUDY_118 +
                  e_study_192_ovdp * STUDY_192 +
                  e_study_221_ovdp * STUDY_221)

    # -- Observation ----------------------------------------------------
    # The source likelihood is negative-binomial with mean lambda and
    # overdispersion ovdp (Equation 1). rxode2 5.1.7 / nlmixr2 7.0.1 expose
    # llikNbinomMu() and rxnbinomMu() but accept no negative-binomial
    # ENDPOINT, so -- following the ddmore/Schoemaker_2018_levetiracetam.R
    # and ddmore/Plan_2012_pain.R precedents in this library -- a Poisson
    # likelihood with the same mean is declared and ovdp is exposed as a
    # model variable. The deterministic mean-rate trajectory is unaffected;
    # only the dispersion of stochastic draws differs. The vignette draws
    # true negative-binomial counts with rnbinom(mu = lambda, size = 1/ovdp),
    # the size mapping implied by variance = lambda * (1 + ovdp * lambda).
    seizure_count <- lambda
    seizure_count ~ pois(lambda)
  })
}

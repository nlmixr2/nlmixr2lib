Burns_2015_tenofovir <- function() {
  description <- paste(
    "Two-compartment population PK model for plasma tenofovir after oral",
    "tenofovir disoproxil fumarate in healthy women, parameterised with",
    "micro-rate constants and first-order absorption with a lag time, linked",
    "by a first-order rate constant to a peripheral-blood-mononuclear-cell",
    "tenofovir-diphosphate compartment with first-order elimination; body",
    "weight on central volume, and a fixed adherence-adjustment",
    "bioavailability on self-administered (unobserved) doses."
  )
  reference <- paste(
    "Burns RN, Hendrix CW, Chaturvedula A. Population pharmacokinetics of",
    "tenofovir and tenofovir-diphosphate in healthy women.",
    "J Clin Pharmacol. 2015;55(6):629-638. doi:10.1002/jcph.461"
  )
  vignette <- "Burns_2015_tenofovir"

  # Methods, "Population Pharmacokinetic Model Development": "Dosing was in TFV
  # equivalents (136 mg TFV/300 mg TDF) and in micromoles whereas TFV
  # concentrations were converted to nanomoles/mL (TFV molecular weight of
  # 288.1 g/mole)." Supplementary Figure S2 labels the dose box "micromoles TFV
  # equivalents" and the plasma box "Plasma, TFV nm/ml". One 300 mg tenofovir
  # disoproxil fumarate tablet therefore enters the depot as
  # 136 / 288.1 * 1000 = 472.058 umol.
  units <- list(time = "h", dosing = "umol", concentration = "nmol/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Supplementary Figure S2 (structural
  # model overview) and the Methods unit-conversion paragraph.
  #
  # NOTE on the `pbmc_tfvdp` state's units. The authors fitted the system with
  # NONMEM ADVAN5, whose rate constants act on compartment AMOUNTS, so the
  # state is an amount in umol carried in the same mass balance as `central`
  # ("By using the linear ADVAN5 subroutine we included TFV-DP in the mass
  # balance equations", Discussion). No scaling parameter S4 / V4 is reported
  # anywhere in the paper or in Figure S2 -- the figure labels the PBMC box
  # directly as "PBMC, TFV-DP nm/L" with no volume -- so the NONMEM default
  # S4 = 1 applies and the compartment amount IS the predicted TFV-DP
  # concentration numerically. The whole amount-to-concentration scaling is
  # absorbed into K24, which the authors describe as apparent ("we modeled only
  # 1 million PBMCs and thus the parameters were apparent in theory",
  # Discussion). See the `kmet_tfvdp` note in ini() for the numeric check that
  # settles this reading.
  compartmentData <- list(
    depot       = list(analyte = "tenofovir", units = "umol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tenofovir", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tenofovir", units = "umol", specimen = "plasma", verified = TRUE),
    pbmc_tfvdp  = list(analyte = "tenofovir diphosphate", units = "umol", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Enters Vc/F as a",
        "LINEAR, NOT allometric, term centred on the cohort median weight of",
        "73 kg: Results gives the relationship verbatim as",
        "'TVV = 385.71 - 2.16 * (73 - weight (kg))', so Vc/F RISES with weight",
        "(equivalently 385.71 + 2.16 * (WT - 73)). Adding it dropped the",
        "objective function by 11 points and cut between-subject variability on",
        "Vc from 24.28% to 19.3% CV. Cohort weight mean 79.6 kg, range 43-145 kg",
        "(Supplementary Table S1); over that range Vc/F spans about 321-541 L.",
        "The published equation's minus-sign form is preserved verbatim in",
        "model(), with e_wt_vc = -2.16 as tabulated in Table 1, rather than",
        "being algebraically re-centred, so that the coefficient in ini()",
        "matches the coefficient in the paper's own table."
      ),
      source_name        = "WT"
    ),
    SELFADMIN = list(
      description        = paste(
        "1 = the dose was self-administered at home without study-staff",
        "supervision (the 'preclinic' dose of the end-of-period visit);",
        "0 = the dose was taken under observation in clinic."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (observed in-clinic dose; relative bioavailability 1)",
      notes              = paste(
        "Carries the paper's adherence adjustment, which is the methodological",
        "point of the analysis. MTN-001 adherence was suspected to be",
        "suboptimal, so the authors applied Gibiansky's method: a",
        "bioavailability parameter on the previous day's (self-reported,",
        "unobserved) dose, FIXED to 0.5 with its between-subject variance FIXED",
        "to a deliberately large value so that each subject's empirical Bayes",
        "estimate can absorb either missed doses (low F1) or extra doses (high",
        "F1). Methods: 'The bioavailability for the in-clinic dose was set to 1;",
        "thus all estimated parameters are apparent.' The estimation dataset",
        "therefore has exactly one SELFADMIN = 1 record (the preclinic dose,",
        "mean 12.9 h, range 1.95-35 h before the clinic dose) and one",
        "SELFADMIN = 0 record per subject-period.",
        "FOR SIMULATION, set SELFADMIN = 0 on every dose to reproduce the",
        "full-compliance scenario the authors simulated for their steady-state",
        "and single-dose qualification runs; the F1 machinery then collapses to",
        "F = 1 and plays no part.",
        "Polarity matches the canonical (self-administered = 1) and matches the",
        "paper, in which the unobserved dose is the one carrying the",
        "adjustment."
      ),
      source_name        = "F1 (preclinic dose)"
    )
  )

  # Covariates the authors screened but did NOT retain in the final model.
  # Documented here for provenance; neither is referenced in model().
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance.",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened on K20 with a linear relationship and significant in forward",
        "addition (-4 objective-function points), but Results states that once",
        "weight on Vc was in the model 'there was no significant drop in the",
        "objective function with the addition of CrCl on K20 and the model",
        "failed to converge', so only weight on Vc was retained. Cohort mean",
        "130.1 mL/min, range 64-257 mL/min (Supplementary Table S1). The paper",
        "does not report the estimating equation used for creatinine clearance."
      )
    ),
    RACE_BLACK = list(
      description = "Race indicator, 1 = Black, 0 = non-Black.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on clearance as a two-level black / non-black split and not",
        "retained: 'Race ... on clearance was not shown to be a significant",
        "covariate and this finding was in agreement with other literature",
        "reports.' Results further notes the variable was confounded with trial",
        "site, since the African sites contributed the Black participants.",
        "Cohort composition White 32, Black 60, Other 9 (Supplementary",
        "Table S1)."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 101,
    n_observations   = 875,
    n_studies        = 1,
    age_range        = "18-45 years",
    age_mean         = "31.3 years",
    weight_range     = "43-145 kg",
    weight_mean      = "79.6 kg",
    weight_median    = "73 kg",
    sex_female_pct   = 100,
    race_ethnicity   = c(White = 32, Black = 60, Other = 9),
    disease_state    = "healthy HIV-negative women (pre-exposure prophylaxis target population)",
    renal_function   = "creatinine clearance mean 130.1 (range 64-257) mL/min",
    dose_range       = "tenofovir disoproxil fumarate 300 mg orally once daily (136 mg tenofovir equivalents = 472.058 umol)",
    regions          = "United States, South Africa, Uganda, Zimbabwe (MTN-001 sites)",
    notes            = paste(
      "MTN-001, a 21-week Phase II open-label three-period crossover study of",
      "daily oral tenofovir disoproxil fumarate 300 mg and/or 1% vaginal",
      "tenofovir gel in healthy women. 168 enrolled, 144 completed at least one",
      "follow-up visit in each period, 141 had PK measurements, and 101 entered",
      "this analysis contributing 476 plasma tenofovir and 399 PBMC",
      "tenofovir-diphosphate concentrations (Results; demographics in",
      "Supplementary Table S1, which reports means rather than medians, except",
      "for the 73 kg weight median used to centre the weight covariate).",
      "ONLY the oral period's end-of-visit data were used for model building;",
      "the dual (oral + vaginal gel) period was held back for internal",
      "qualification, and the vaginal-only period was not modelled at all.",
      "Below-quantification-limit data were excluded (3.4% of tenofovir and",
      "12.9% of tenofovir-diphosphate observations); the M3 method was",
      "attempted but did not converge. Eight further points were excluded as",
      "outliers. Sampling was intensive at some sites (1, 2, 4, 6 and 8 h",
      "postdose) and a single postdose sample between 1 and 8 h at others.",
      "The race counts are subject counts, not percentages, and sum to 101."
    )
  )

  ini({
    # ---- Structural parameters: Table 1, "Final Model / Value (%RSE)" column.
    # The authors parameterised in NONMEM micro-rate constants rather than
    # clearances: "Micro-rate constants were used as typical Cl, V
    # parametrization encountered numerical problems and the TFV-DP portion of
    # the model must be described in micro-rate constants" (Results). Their
    # ADVAN5 compartment numbering is 1 = depot, 2 = plasma tenofovir,
    # 3 = peripheral tenofovir, 4 = PBMC tenofovir-diphosphate
    # (Supplementary Figure S2), which maps onto the canonical names as
    # K23 -> k12, K32 -> k21, K20 -> kel, K24 -> kmet_tfvdp, K40 -> kel_tfvdp.
    lka <- log(9.79); label("First-order absorption rate constant KA (1/h)")                             # Table 1 final model (KA = 9.79 1/h, 65.18% RSE)
    lvc <- log(385.71); label("Apparent central volume Vc/F at the median weight of 73 kg (L)")          # Table 1 final model (Vc/F = 385.71 L, 14.84% RSE)
    lkel <- log(0.13); label("Elimination rate constant K20 out of plasma tenofovir (1/h)")              # Table 1 final model (K20 = 0.13 1/h, 17.81% RSE)
    lk12 <- log(0.631); label("Plasma-to-peripheral tenofovir rate constant K23 (1/h)")                  # Table 1 final model (K23 = 0.631 1/h, 24.7% RSE)
    lk21 <- log(0.396); label("Peripheral-to-plasma tenofovir rate constant K32 (1/h)")                  # Table 1 final model (K32 = 0.396 1/h, 23.24% RSE)
    ltlag <- log(0.5); label("Absorption lag time (h)")                                                  # Table 1 final model (Absorption lag = 0.5 h, 35.49% RSE); adding it dropped the objective function by 192

    # K24 is simultaneously the tenofovir-diphosphate formation rate constant
    # and a second apparent elimination pathway for plasma tenofovir --
    # "It is important to note that both K20 and K24 are apparent clearance
    # terms for TFV. Therefore, when assessing TFV clearance both terms were
    # included" (Results) -- which is what makes the published total apparent
    # clearance (K20 + K24) * Vc = 56.7 L/h an exact closed-form check on this
    # parameterisation.
    #
    # K24 also carries the entire amount-to-concentration scaling of the PBMC
    # compartment, because no S4 / V4 is reported (see the compartmentData
    # note). The reading is settled numerically and non-circularly against the
    # paper's OWN steady-state simulation: solving this model at its typical
    # values with 300 mg tenofovir disoproxil fumarate daily and full
    # compliance gives a tenofovir-diphosphate trough of 167 (compartment
    # units read as nmol/L) = 47.0 fmol per million cells after dividing by the
    # paper's 282 fL/cell conversion, against the published simulated median
    # trough of 49.9 fmol per million cells (Results, "Simulation of Steady
    # State TFV-DP Levels"). Any reading that introduced a genuine PBMC volume
    # would miss that target by many orders of magnitude.
    lkmet_tfvdp <- log(0.017); label("Apparent tenofovir-diphosphate formation rate constant K24 out of plasma tenofovir (1/h)")  # Table 1 final model (K24 = 0.017 1/h, 72.48% RSE)
    lkel_tfvdp <- log(0.013); label("Tenofovir-diphosphate elimination rate constant K40 (1/h)")         # Table 1 final model (K40 = 0.013 1/h, 16.63% RSE); log(2)/0.013 = 53.3 h, the paper's quoted half-life

    # ---- Covariate effect. Results: "The covariate relationship between
    # weight and volume was linear and was centered on the median weight
    # (73 kg) as follows: TVV = 385.71 - 2.16 * (73 - weight (kg))". Table 1
    # tabulates the coefficient as -2.16 (34.52% RSE) with a bootstrap median
    # of -1.78 and a 95% CI of -3.37 to -0.16, i.e. wholly negative, so the
    # sign is settled by the bootstrap column as well as by the abstract's
    # rendering of the equation.
    e_wt_vc <- -2.16; label("Linear weight effect on Vc/F, entering as Vc = Vc_73kg + e_wt_vc * (73 - WT) (L/kg)")  # Table 1 final model (cov WT on Vc = -2.16, 34.52% RSE)

    # ---- Adherence-adjustment bioavailability on the self-administered
    # (preclinic) dose. Both the point value and its variance were FIXED by
    # design, not estimated: Methods, "fixing the bioavailability parameter to
    # 0.5 and its omega distribution to a high value, in our case 10". Table 1
    # reports for F1 only the empirical Bayes estimates (mean 0.98, range
    # 0.015-4.17), flagged by the table footnote 'F1 values are the empirical
    # Bayes estimates mean (min-max)', which confirms F1 itself was not an
    # estimated typical value.
    lfdepot <- fixed(log(0.5)); label("Relative bioavailability F1 of a self-administered (unobserved) dose, carrying the adherence adjustment (unitless)")  # Methods, Gibiansky's method (F1 = 0.5 FIX); Methods equation 'F1 = 0.5 * exp(eta_i)'

    # ---- Between-subject variability. Modelled exponentially (Methods:
    # 'theta_i = theta_typical * exp(eta_i)'), supported on KA, Vc, K20 and K24
    # (Results). Table 1 reports these as %CV, so the variance on the log scale
    # is omega^2 = log(1 + CV^2).
    etalka ~ 1.2712043      # Table 1 final model 'BSV KA (%CV)' = 160.2; Results text gives 160.16 to one more digit; log(1 + 1.6016^2)
    etalvc ~ 0.0365720      # Table 1 final model 'BSV Vc (%CV)' = 19.3; log(1 + 0.193^2)
    etalkel ~ 0.1232692     # Table 1 final model 'BSV K20 (%CV)' = 36.22; log(1 + 0.3622^2)
    etalkmet_tfvdp ~ 1.2651731  # Table 1 final model 'BSV K24 (%CV)' = 159.49; log(1 + 1.5949^2)

    # Fixed, deliberately enormous variance on the adherence bioavailability.
    # Read as the NONMEM $OMEGA entry, which holds a VARIANCE, so omega^2 = 10
    # (omega = 3.16 on the log scale). This is not a fitted dispersion: its
    # whole purpose is to make the prior on F1 nearly uninformative so the
    # individual empirical Bayes estimate is driven by the data and can span
    # missed doses through to double doses. The authors note it makes ordinary
    # Monte-Carlo simulation useless -- 'the Monte-Carlo method that NONMEM
    # uses for simulation generated a large distribution of F1 values that
    # resulted in exaggerated prediction intervals' -- and imputed empirical
    # Bayes estimates instead. Simulate with SELFADMIN = 0 unless you are
    # deliberately reproducing the estimation dataset's dosing structure.
    etalfdepot ~ fixed(10)  # Methods, Gibiansky's method ('its omega distribution to a high value, in our case 10')

    # ---- Residual error. Methods specified a combined additive-plus-
    # proportional form and states 'Reduction of the residual error model was
    # tested during model development'; Table 1 reports a proportional term
    # only for each analyte, so the final model is proportional-only with a
    # separate term per analyte.
    propSd <- 0.2748; label("Proportional residual error for plasma tenofovir (fraction)")                        # Table 1 final model 'Proportional, TFV (%CV)' = 27.48
    propSd_Cpbmc_tfvdp <- 0.3118; label("Proportional residual error for PBMC tenofovir-diphosphate (fraction)")        # Table 1 final model 'Proportional, TFV-DP (PBMC) (%CV)' = 31.18
  })

  model({
    # 1. Individual parameters. The weight effect is additive in litres on the
    #    typical central volume, then the exponential between-subject term is
    #    applied to the covariate-adjusted typical value.
    ka <- exp(lka + etalka)
    vc <- (exp(lvc) + e_wt_vc * (73 - WT)) * exp(etalvc)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kmet_tfvdp <- exp(lkmet_tfvdp + etalkmet_tfvdp)
    kel_tfvdp <- exp(lkel_tfvdp)

    # 2. ODE system, Supplementary Figure S2. Tenofovir leaves `central` by
    #    three routes -- kel (K20) out of the system, k12 (K23) to the
    #    peripheral compartment, and kmet_tfvdp (K24) into the PBMC
    #    tenofovir-diphosphate pool -- and the last of these is part of the
    #    mass balance, which is why the paper's total apparent clearance is
    #    (K20 + K24) * Vc rather than K20 * Vc.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot -
      (kel + k12 + kmet_tfvdp) * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(pbmc_tfvdp) <- kmet_tfvdp * central - kel_tfvdp * pbmc_tfvdp

    alag(depot) <- exp(ltlag)

    # 3. Adherence adjustment. A self-administered (unobserved) dose carries
    #    the fixed F1 = 0.5 and its very wide between-subject term; an
    #    observed in-clinic dose has bioavailability 1 by construction, which
    #    is what makes every other parameter an apparent one.
    f(depot) <- ifelse(SELFADMIN == 1, exp(lfdepot + etalfdepot), 1)

    # 4. Observations. Plasma tenofovir is an amount in umol over a volume in
    #    L, giving umol/L = nmol/mL as declared in `units` and drawn in
    #    Figure S2. The tenofovir-diphosphate state is read directly, with no
    #    volume: see the compartmentData and kmet_tfvdp notes. Its numeric
    #    value is nmol/L, so divide by 1e-6 / (282e-15 * 1e6) = 3.5461 to
    #    recover the fmol per million cells the assay reported.
    Cc <- central / vc
    Cpbmc_tfvdp <- pbmc_tfvdp

    Cc ~ prop(propSd)
    Cpbmc_tfvdp ~ prop(propSd_Cpbmc_tfvdp)
  })
}

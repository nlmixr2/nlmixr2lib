Nyangwa_2026_pretomanid <- function() {
  description <- paste(
    "One-compartment population PK model for oral pretomanid 200 mg once",
    "daily in adults with rifampicin-resistant tuberculosis (RR-TB) treated",
    "with the BPaL, BPaLM, or BPaLC regimen in the TB-PRACTECAL trial PKPD",
    "sub-study (Nyang'wa 2026). First-order absorption (ka 0.316 1/h, no",
    "between-subject variability) into a single central compartment with",
    "first-order elimination; apparent clearance CL/F 3.10 L/h and apparent",
    "central volume V/F 102 L, both at the cohort median fat-free mass of",
    "45.5 kg. Fat-free mass is the only retained covariate and enters as",
    "a priori allometric scaling with fixed exponents 0.75 on clearance and",
    "1 on volume; FFM was selected over total body weight and body mass",
    "index on base-model fit. Between-subject variability is a diagonal",
    "pair on clearance (32.9% CV) and central volume (33.6% CV), and",
    "residual error is combined proportional (32.2%) plus additive",
    "(0.368 mg/L). Female sex, Black race, and the BPaL regimen were",
    "significant on volume in forward selection but none survived backward",
    "elimination, so the final model carries no covariate other than FFM."
  )
  reference <- paste(
    "Nyang'wa BT, Motta I, Moodliar R, Solodovnikova V, Rajaram S, Rasool M,",
    "Berry C, Huang Z, Davies G, Moore DAJ, Kloprogge F (2026).",
    "Population pharmacokinetics and target attainment of pretomanid in",
    "rifampicin-resistant tuberculosis patients.",
    "Sci Rep 16:46217. doi:10.1038/s41598-026-46217-2",
    sep = " "
  )
  vignette <- "Nyangwa_2026_pretomanid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "pretomanid",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "pretomanid",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate retained in the final model. Enters as a priori",
        "allometric body-size scaling on both apparent clearance (fixed",
        "exponent 0.75) and apparent central volume (fixed exponent 1), per",
        "Nyang'wa 2026 Materials and methods: 'Allometric body size scaling",
        "was applied a-priori to clearance and volume parameters, with fixed",
        "exponents of 0.75 for clearance and 1 for volume'. FFM was chosen",
        "over baseline body weight and body mass index because it gave the",
        "largest drop in objective function value relative to a model with",
        "no body-size scaling (Results paragraph 4).",
        "",
        "REFERENCE VALUE. The supplementary Appendix 1 nlmixr2 code scales",
        "on a pre-computed data column named logFFM",
        "(cl <- exp(lcl + eta.cl + logFFM * covffmPow1)) and never states",
        "the constant that column is centred on. The centring is recovered",
        "arithmetically and is not an assumption about the model form: the",
        "code's exp(lcl) = 3.096 L/h and exp(lvc) = 101.7 L reproduce the",
        "Table 2 typical values of 3.10 L/h and 102 L exactly, so logFFM",
        "must be zero at the reference, i.e. logFFM = log(FFM / FFM_ref).",
        "The value of FFM_ref is then pinned by the paper's own reported",
        "median AUC(0-24) of 64,000 ug*h/L (Table 3): at steady state",
        "AUC(0-24) = Dose / CL, so the median subject's CL is",
        "200 mg / 64 mg*h/L = 3.13 L/h, which equals exp(lcl) only when the",
        "reference FFM is the cohort median. Taking the Table 1 pretomanid",
        "PK-cohort median FFM of 45.5 kg gives a predicted median AUC(0-24)",
        "of 64,600 ug*h/L (0.9% from the reported 64,000); the conventional",
        "70 kg adult reference would instead give 89,200 ug*h/L (39% high)",
        "and an uncentred log(FFM) would give 3,690 ug*h/L, so both are",
        "excluded. See the vignette Errata for the full adjudication.",
        "",
        "The paper does not report which equation produced its FFM column.",
        "A downstream user must supply FFM directly, or derive it with a",
        "documented formula such as Janmahasatian et al. (Clin Pharmacokinet",
        "2005;44:1051-1065); the cohort's median total body weight of",
        "56.8 kg against a median FFM of 45.5 kg is the only in-paper",
        "calibration point available."
      ),
      source_name = "FFM"
    )
  )

  # Covariates the paper screened but did NOT retain in the final model.
  # Documentation only: none is referenced in model(), and the paper reports
  # no usable point estimate for any of them. Appendix 2 is a covariate-
  # versus-eta correlation matrix over the continuous entries below; the
  # stepwise result is Results paragraph 4.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened head-to-head against FFM and BMI as the allometric",
        "body-size descriptor and lost: FFM 'improved the base model fit",
        "better than weight or BMI' (Results paragraph 4). Cohort median",
        "56.8 kg (range 39.2-144.4), Table 1."
      )
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = paste(
        "Screened as an alternative allometric body-size descriptor and",
        "lost to FFM (Results paragraph 4). Cohort median 19.7 kg/m^2",
        "(range 14.3-47.1), Table 1."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Significant on volume of distribution at the forward-inclusion",
        "threshold (p < 0.05, dOFV > 3.84) but eliminated at the backward",
        "threshold (p < 0.001, dOFV > 10.83), so it carries no estimate in",
        "the final model (Results paragraph 4, Appendix 2). 34 of 94",
        "participants (36.2%) were female, Table 1."
      )
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Significant on volume of distribution in forward inclusion but",
        "eliminated in backward elimination (Results paragraph 4). 52 of 94",
        "participants (55.3%) were Black, 40 (42.6%) Caucasian, 1 Asian and",
        "1 other, Table 1."
      )
    ),
    CONMED_MOXIFLOXACIN = list(
      description = "Concomitant moxifloxacin (the BPaLM regimen arm)",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "The regimen arm was screened on clearance and did not improve the",
        "fit: 'including BPaLM and BPaLC as covariates on clearance did not",
        "improve the model fit significantly, suggesting none or limited",
        "impact of the accompanying anti-TB drugs in the regimen on",
        "pretomanid exposure' (Discussion paragraph 4). 38 of 94",
        "participants (40.4%) received BPaLM, Table 1."
      )
    ),
    CONMED_CLOFAZIMINE = list(
      description = "Concomitant clofazimine (the BPaLC regimen arm)",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Screened on clearance with BPaLM and not retained (Discussion",
        "paragraph 4). Separately, 'BPaL regimen' was significant on volume",
        "of distribution in forward inclusion but eliminated in backward",
        "elimination (Results paragraph 4). 30 of 94 participants (31.9%)",
        "received BPaLC and 26 (27.7%) received BPaL alone, Table 1."
      )
    ),
    HIV_POS = list(
      description = "HIV-positive comorbidity indicator",
      units = "(binary)",
      type = "categorical",
      notes = paste(
        "Not a significant covariate; no further exploration of individual",
        "antiretroviral effects was performed (Discussion paragraph 4). All",
        "39 participants living with HIV (41.5%, Table 1) were on integrase",
        "inhibitor plus nucleoside/nucleotide reverse transcriptase",
        "inhibitor regimens. Potent CYP450 inducers such as efavirenz and",
        "lopinavir/ritonavir, which are known to reduce pretomanid exposure,",
        "were contraindicated in the trial."
      )
    ),
    CRCL = list(
      description = "Estimated creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Included in the Appendix 2 covariate-versus-eta correlation matrix",
        "and not retained. Cohort median 105.4 mL/min (range 43.4-243.8),",
        "Table 1. The trial excluded patients with moderate renal function",
        "abnormality (Discussion, limitations), so the model is not",
        "informed about renal impairment."
      )
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = paste(
        "Included in the Appendix 2 correlation matrix and not retained.",
        "Cohort median 3.6 mmol/L (range 1.7-8.5), Table 1."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = paste(
        "Included in the Appendix 2 correlation matrix and not retained.",
        "Cohort median 19.5 IU/L (range 4-113), Table 1. The trial excluded",
        "patients with moderate liver function abnormality (Discussion,",
        "limitations)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = paste(
        "Included in the Appendix 2 correlation matrix and not retained.",
        "Cohort median 22 IU/L (range 4-82), Table 1."
      )
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Included in the Appendix 2 correlation matrix and not retained.",
        "Cohort median 77 g/L (range 61-118), Table 1. Note this is a",
        "screened PK covariate and is unrelated to the 85-95% plasma",
        "protein binding assumed in the paper's PTA analysis, which is an",
        "assumption carried from the literature rather than a measurement:",
        "'pretomanid protein binding in vivo has not yet been determined'",
        "(Introduction paragraph 3)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 1L,
    age_range = "19-71 years (median 36)",
    age_median = "36 years",
    weight_range = "39.2-144.4 kg (median 56.8)",
    weight_median = "56.8 kg",
    ffm_range = "28.6-75.5 kg (median 45.5); the allometric reference used here",
    bmi_range = "14.3-47.1 kg/m2 (median 19.7)",
    sex_female_pct = 36.2,
    race_ethnicity = "Black 52 (55.3%), Caucasian 40 (42.6%), Asian 1 (1.1%), other 1 (1.1%)",
    disease_state = paste(
      "Rifampicin-resistant pulmonary tuberculosis. 39 participants (41.5%)",
      "were living with HIV, all on integrase-inhibitor plus nucleoside /",
      "nucleotide reverse transcriptase inhibitor antiretroviral therapy.",
      "Patients with moderate liver or renal function abnormality were",
      "excluded from the trial. Median estimated creatinine clearance",
      "105.4 mL/min, median ALT 19.5 IU/L, median AST 22 IU/L."
    ),
    dose_range = paste(
      "Pretomanid 200 mg orally once daily for 24 weeks in every arm.",
      "Co-administered with bedaquiline (400 mg daily for 2 weeks then",
      "200 mg three times weekly for 22 weeks) and linezolid (600 mg daily",
      "for 16 weeks then 300 mg daily for 8 weeks), plus moxifloxacin",
      "400 mg daily in the BPaLM arm (38 participants, 40.4%) or",
      "clofazimine 100 mg daily in the BPaLC arm (30, 31.9%); 26 (27.7%)",
      "received BPaL alone. Participants were encouraged to eat before",
      "dosing but meals were neither standardised nor recorded, so the",
      "estimates reflect real-world mixed fed / fasted absorption."
    ),
    regions = "South Africa and Belarus",
    notes = paste(
      "PRACTECAL-PKPD sub-study of the TB-PRACTECAL randomised controlled",
      "trial (ClinicalTrials.gov NCT04081077). 952 timed plasma samples",
      "(86 pre-first-dose, 866 post-dose) spanning the full 24-week",
      "treatment course and follow-up visits to week 72. Sampling was on",
      "day 1 (0, 2, 23 h), week 8 (predose, 6.5, 23 h), and weeks 12, 16,",
      "20, 24, 32 and 72. Observed concentrations ranged 19.1-11,566 ng/mL",
      "with a median trough of 1,789 ng/mL (IQR 1,126-2,689). The assay",
      "lower limit of quantification was 7 ng/mL; 234 samples were below",
      "it, of which 151 were collected after treatment completion, leaving",
      "9.5% of on-treatment samples BLQ, handled by the M1 method",
      "(discarded as missing). Estimation used FOCE-I in nlmixr2 under",
      "R 4.1.2. Baseline characteristics are Table 1; the final parameter",
      "estimates are Table 2 and the nlmixr2 model code is",
      "supplementary Appendix 1."
    )
  )

  ini({
    # --- Structural parameters.
    #
    # Every value below is the supplementary Appendix 1 nlmixr2 ini() block
    # transcribed at full precision. Appendix 1 is the authors' own final
    # model code, so these are final estimates and not initial values; each
    # back-transforms to the Table 2 point estimate, which is quoted in the
    # trailing comment as the independent cross-check.
    lka <- -1.15259732127294
    label("First-order absorption rate constant from depot to central (1/h)")   # Appendix 1 lka; exp() = 0.3158, Table 2 k a 0.316 (RSE 19.6%), bootstrap 95% CI 0.203-0.492
    lcl <- 1.13008124732802
    label("Apparent clearance CL/F at FFM = 45.5 kg (L/h)")                     # Appendix 1 lcl; exp() = 3.0959, Table 2 CL/F 3.10 (RSE 3.35%), bootstrap 95% CI 2.87-3.33
    lvc <- 4.62170931536226
    label("Apparent central volume of distribution V/F at FFM = 45.5 kg (L)")   # Appendix 1 lvc; exp() = 101.67, Table 2 V/F 102 (RSE 2.45%), bootstrap 95% CI 81.4-127

    # --- Covariate effects. Both exponents were fixed a priori to the
    # theory-based allometric values rather than estimated (Materials and
    # methods paragraph 4), and Appendix 1 wraps both in fix().
    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on apparent clearance (unitless)")        # Appendix 1 covffmPow1 <- fix(0.75); Materials and methods 'fixed exponents of 0.75 for clearance'
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on apparent central volume (unitless)")   # Appendix 1 covffmPow2 <- fix(1); Materials and methods 'and 1 for volume'

    # --- Between-subject variability, a diagonal pair on clearance and
    # central volume with no random effect on absorption (Appendix 1 has
    # only eta.cl and eta.vc). The Appendix 1 values are log-scale
    # VARIANCES, confirmed against the Table 2 %CV column via the
    # log-normal identity CV = sqrt(exp(omega^2) - 1):
    # sqrt(exp(0.1026746) - 1) = 32.88% and sqrt(exp(0.1070584) - 1) =
    # 33.62%, reproducing the printed 32.9% and 33.6%. Reading the Table 2
    # percentages instead as omega standard deviations would give variances
    # of 0.1082 and 0.1129, which are 5.4% high; the Appendix 1 variances
    # are used verbatim.
    etalcl ~ 0.102674627530168
    label("Between-subject variability in apparent clearance (log-scale variance)")       # Appendix 1 eta.cl; Table 2 CL/F %CV 32.9 [shrinkage 9.65%]
    etalvc ~ 0.10705836808838
    label("Between-subject variability in apparent central volume (log-scale variance)")  # Appendix 1 eta.vc; Table 2 V/F %CV 33.6 [shrinkage 36.7%]

    # --- Residual error, combined proportional plus additive (Results
    # paragraph 3: 'a combined residual error model was used to characterise
    # the unexplained variability'). Appendix 1 declares each as
    # c(lower_bound, estimate), so the estimates are the second elements.
    #
    # UNITS. Appendix 1 gives the additive term as 367.83, which is on the
    # ng/mL scale the paper reports its concentrations in; Table 2 prints the
    # same quantity as 0.368 with an explicit '(mg/L)' unit tag. This model
    # declares concentration in mg/L (units$concentration above), so the
    # Appendix 1 value is carried at full precision divided by 1000. The two
    # sources agree to every printed digit -- 367.831.../1000 = 0.36783...
    # rounds to the tabulated 0.368 -- so this is a unit expression, not a
    # rescaling choice.
    propSd <- 0.321731821549103
    label("Proportional residual error (fraction)")                                       # Appendix 1 prop.err <- c(0, 0.321731821549103); Table 2 Proportional 0.322
    addSd <- 0.367831309683669
    label("Additive residual error (mg/L)")                                               # Appendix 1 add.err <- c(0, 367.831309683669) ng/mL = 0.3678 mg/L; Table 2 Additive (mg/L) 0.368
  })

  model({
    # 1. Individual parameters. Appendix 1 writes the allometry additively on
    #    the log scale against a pre-computed, centred data column,
    #      cl <- exp(lcl + eta.cl + logFFM * covffmPow1)
    #      vc <- exp(lvc + eta.vc + logFFM * covffmPow2)
    #    which is algebraically the power-of-ratio form below once logFFM is
    #    resolved to log(FFM / 45.5). The reference of 45.5 kg is the
    #    Table 1 pretomanid-cohort median fat-free mass; it is recovered
    #    from the paper rather than assumed, because exp(lcl) and exp(lvc)
    #    equal the Table 2 typical values (so the column is zero at the
    #    reference) and because only a cohort-median reference reproduces
    #    the Table 3 median AUC(0-24) of 64,000 ug*h/L. The full
    #    adjudication, including the two rejected alternatives, is in
    #    covariateData[[FFM]]$notes and the vignette Errata.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (FFM / 45.5)^e_ffm_cl
    vc <- exp(lvc + etalvc) * (FFM / 45.5)^e_ffm_vc

    kel <- cl / vc

    # 2. One-compartment disposition with first-order absorption (Results
    #    paragraph 3: 'A one-compartment first order absorption and
    #    elimination model best described the observed pretomanid PK time
    #    series data'). Appendix 1 expresses this as linCmt(); the explicit
    #    ODE form below is the same system and makes the depot and central
    #    states addressable for simulation.
    #
    #    Bioavailability is not identifiable from oral-only data and the
    #    authors did not estimate it: CL/F and V/F are apparent parameters,
    #    so no f(depot) term is applied and the full dose enters the depot.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 3. Observation. Dose in mg and vc in L give central/vc in mg/L, which
    #    matches units$concentration and the mg/L scale of addSd.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}

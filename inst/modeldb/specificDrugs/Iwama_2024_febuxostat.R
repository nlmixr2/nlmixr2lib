Iwama_2024_febuxostat <- function() {
  description <- "Integrated two-compartment population PK model for febuxostat spanning Japanese pediatric patients with hyperuricemia including gout and a Japanese adult population of healthy subjects and patients with renal dysfunction (Iwama 2024). First-order absorption from a depot with an absorption lag time, first-order elimination, and apparent (oral) disposition parameters. Apparent clearance is a power function of body weight and of BSA-normalized estimated glomerular filtration rate, each centered on the pooled-analysis median (60.6 kg, 98.6 mL/min/1.73 m^2). Relative bioavailability is reduced to 0.838 when the dose is taken fed, with the fasted state as the structural anchor. Age was screened but not retained, which is what makes the single model applicable across the pediatric and adult populations."
  reference <- paste(
    "Iwama R, Nishida K, Ishii D, Iijima T.",
    "An integrated population pharmacokinetic model of febuxostat in",
    "pediatric patients with hyperuricemia including gout and adult",
    "population of healthy subjects and patients with renal dysfunction.",
    "Pharmacol Res Perspect. 2024;12(6):e70032.",
    "doi:10.1002/prp2.70032.",
    sep = " "
  )
  vignette <- "Iwama_2024_febuxostat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "febuxostat", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "febuxostat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "febuxostat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL/F centered on the pooled-analysis median of",
        "60.6 kg: (WT / 60.6)^0.584 per the Section 3.3 CL/F equation.",
        "The 60.6 kg centering value is stated in Results Section 3.3",
        "(equation) and again in Section 3.5 ('subjects having median body",
        "weight (60.6 kg)'); it is the median of the pooled pediatric +",
        "adult analysis data set, so it does not equal either the pediatric",
        "median (46.70 kg) or the adult median (61.30 kg) in Table 1.",
        "Table 1 cohort range 26.7-94.3 kg. Integrating the pediatric data",
        "is what let the authors resolve the weight effect at all - the",
        "Discussion notes that earlier adult-only febuxostat popPK analyses",
        "failed to detect it because their weight distributions were",
        "clustered above 60-70 kg."
      ),
      source_name        = "WGT"
    ),
    CRCL = list(
      description        = "BSA-normalized estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL/F centered on the pooled-analysis median of",
        "98.6 mL/min/1.73 m^2: (CRCL / 98.6)^0.324 per the Section 3.3 CL/F",
        "equation; the centering value is restated in Section 3.5 and in the",
        "Figure 3 and Figure 5 captions. Assay form: eGFR computed from the",
        "equation recommended by the Japanese Society of Nephrology and the",
        "Japanese Society for Pediatric Nephrology (Methods Section 2.1,",
        "reference 46; the equation itself is in the Supporting method,",
        "which is not on disk). This is NOT the Cockcroft-Gault creatinine",
        "clearance, which the paper also tabulates separately as CLCr in",
        "mL/min (Table 1) - the authors tested both and retained eGFR",
        "because it gave the greater OFV reduction (Discussion). Table 1",
        "cohort range 20.9-145.3 mL/min/1.73 m^2, spanning normal through",
        "severe renal dysfunction."
      ),
      source_name        = "EGFR"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator (1 = fed, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (fasted), which the paper fixes as the bioavailability anchor:",
        "'F (FASTED) = 1' in Section 3.3"
      ),
      notes              = paste(
        "Per-dose-record indicator. Methods Section 2.1 gives the",
        "operational definition: 'fed: orally administered within 0.5 h",
        "after a meal; fasted: orally administered under fasting conditions",
        "or more than 0.5 h after a meal'. That is a general",
        "fed-vs-fasted flag rather than a high-fat-meal challenge, so the",
        "canonical FED column applies and not FED_HIGHFAT (same reading as",
        "Chen_2023_nemonoxacin.R). Single multiplicative power-form effect",
        "on relative bioavailability: F1 = 1 x 0.838^FED (Section 3.3).",
        "The Discussion states the resulting exposure change directly:",
        "'the present PopPK analysis suggests that relative bioavailability,",
        "AUCtau,ss, and Cmax,ss were each reduced by 16.2% in the fed",
        "group' - and 1 - 0.838 = 0.162 exactly. Note that fasted/fed status",
        "was NOT retained on ka; the Discussion flags that omission as a",
        "likely reason for the very poor precision of the ka IIV."
      ),
      source_name        = "FED"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Methods Section 2.2.2) but not",
        "retained: 'Age was not included as a significant covariate'",
        "(Section 3.3). This is the load-bearing negative result of the",
        "paper - because age drops out, the single integrated model applies",
        "to pediatric patients and adults alike (Section 3.7)."
      )
    ),
    CLCR_RAW = list(
      description = "Creatinine clearance by the Cockcroft-Gault formula (not BSA-normalized)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate on CL/F and rejected in favour of",
        "eGFR: 'the impact of renal function on CL/F was represented by",
        "incorporating eGFR, which showed a greater reduction in the OFV",
        "than when using CLCr' (Discussion). Tabulated in Table 1 (cohort",
        "median 109.83 mL/min pediatric, 121.02 mL/min adult). Retained here",
        "only to record that the two renal-function metrics were tested",
        "head-to-head; the model uses CRCL (the BSA-normalized eGFR)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Not considered as a candidate covariate because the cohort was",
        "overwhelmingly male: 'there were only 9 females (6 children and 3",
        "adults) compared to 133 males' (Discussion). Post hoc CL/F adjusted",
        "by body weight showed no clear male-female difference (Figure S2)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 142,
    n_studies      = 6,
    age_range      = "8-72 years (pediatric 8-18; adult 20-72)",
    age_median     = "13.0 years (pediatric); 24.0 years (adult)",
    weight_range   = "26.7-94.3 kg",
    weight_median  = "46.70 kg (pediatric); 61.30 kg (adult); 60.6 kg (pooled analysis median used for covariate centering)",
    sex_female_pct = 6.3,
    race_ethnicity = c(Japanese = 100),
    disease_state  = paste(
      "29 pediatric patients with hyperuricemia including gout;",
      "113 adults who were healthy or had renal dysfunction.",
      "Renal function by eGFR category across the whole analysis set:",
      "87 normal, 34 mild, 19 moderate, 2 severe."
    ),
    dose_range     = paste(
      "Pediatric: 5, 10, 20 or 30 mg once daily orally for body weight",
      "< 40 kg and 10, 20, 40 or 60 mg once daily for body weight >= 40 kg,",
      "for 52 weeks with up-titration as needed.",
      "Adult: 10-160 mg as single or repeated oral doses."
    ),
    regions        = "Japan",
    n_observations = 2611,
    renal_function = paste(
      "Deliberately enriched for renal impairment. Pediatric eGFR median",
      "76.65 (range 33.9-145.3) mL/min/1.73 m^2; adult eGFR median 99.47",
      "(range 20.9-143.6). The proportion with reduced renal function was",
      "higher among the pediatric patients than the adults (Section 3.1)."
    ),
    notes          = paste(
      "Pooled analysis of six Japanese studies: two Phase 2 pediatric",
      "studies (Study 1 evaluation phase, Study 2 continuation phase),",
      "three Phase 1 studies in healthy adult males (Studies 3, 4, 5), and",
      "one Phase 1 repeated-dose study in adults with normal renal function",
      "or renal dysfunction (Study 6). 2611 plasma concentration records",
      "from 142 subjects (110 from pediatric patients < 40 kg, 190 from",
      "pediatric patients >= 40 kg, 2311 from adults). Baseline demographics",
      "in Table 1; per-study designs in Table S1 (supplement not on disk).",
      "30 pediatric patients were enrolled in the Phase 2 studies but 29",
      "contributed to the PopPK analysis set."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Iwama 2024 Table 2 'Fixed Effect' block; the
    # same six typical values are restated in the Section 3.3 prose. All
    # disposition parameters are APPARENT (oral) quantities, i.e. CL/F,
    # V2/F, Q/F, V3/F, because no intravenous data were available - the
    # fasted bioavailability is the structural anchor (see lfdepot).
    # ------------------------------------------------------------------
    lcl <- log(6.53); label("Apparent clearance CL/F at WT = 60.6 kg and eGFR = 98.6 mL/min/1.73 m^2 (L/h)") # Table 2 row 'CL/F' = 6.53 L/h (RSE 3.0%, bootstrap 6.52, 95% CI 6.14-6.92); Section 3.3 equation intercept
    lvc <- log(19.4); label("Apparent central volume of distribution V2/F (L)")                              # Table 2 row 'V2/F' = 19.4 L (RSE 4.3%, bootstrap 19.4, 95% CI 17.8-21.0)
    lq  <- log(1.80); label("Apparent intercompartmental clearance Q/F (L/h)")                               # Table 2 row 'Q/F' = 1.80 L/h (RSE 7.2%, bootstrap 1.78, 95% CI 1.55-2.05)
    lvp <- log(15.6); label("Apparent peripheral volume of distribution V3/F (L)")                           # Table 2 row 'V3/F' = 15.6 L (RSE 5.3%, bootstrap 15.5, 95% CI 14.0-17.2)

    # ka is estimated with poor precision (RSE 24.6%, bootstrap 95% CI
    # 2.47-14.5 1/h) and carries an extreme IIV (CV 306.3%, see etalka).
    # The Discussion attributes this to fasted/fed status not being carried
    # onto ka, compounded by febuxostat being a weakly acidic, poorly
    # soluble BCS class II compound whose absorption tracks gastric
    # emptying. The point estimate is used as reported.
    lka   <- log(3.43);  label("First-order absorption rate constant ka (1/h)") # Table 2 row 'Ka' = 3.43 1/h (RSE 24.6%, bootstrap 4.24, 95% CI 1.78-5.08)
    ltlag <- log(0.437); label("Absorption lag time ALAG1 (h)")                 # Table 2 row 'ALAG1' = 0.437 h (RSE 4.6%, bootstrap 0.441, 95% CI 0.398-0.476)

    # Bioavailability. Section 3.3 states the food effect as a ratio
    # against a fixed fasted reference: 'F (FASTED) = 1' and
    # 'F1 (FED) = 0.838'. The leading 1 is a structural anchor - Table 2
    # reports F1 (FED) but no estimated typical F - so it is encoded as
    # fixed, matching the Chen_2023_nemonoxacin.R treatment of the same
    # construct.
    lfdepot <- fixed(log(1)); label("Relative bioavailability anchor F, fasted (fraction)") # Section 3.3 'F (FASTED) = 1'; no typical-value row in Table 2

    # ------------------------------------------------------------------
    # Covariate effects. Section 3.3 gives the CL/F equation verbatim:
    #   CL/F_i = 6.53 x (WGT_i / 60.6)^0.584 x (EGFR_i / 98.6)^0.324
    # Both effects are power functions centered on the pooled-analysis
    # medians, as the surrounding prose states ('modeled using power
    # functions that were centered around their respective medians').
    # ------------------------------------------------------------------
    e_wt_cl   <- 0.584; label("Power exponent of body weight on CL/F, centered at 60.6 kg (unitless)")                 # Table 2 row 'CLWGT' = 0.584 (RSE 20.4%, bootstrap 0.589, 95% CI 0.351-0.817); Section 3.3 equation
    e_crcl_cl <- 0.324; label("Power exponent of eGFR on CL/F, centered at 98.6 mL/min/1.73 m^2 (unitless)")           # Table 2 row 'CLEGFR' = 0.324 (RSE 17.3%, bootstrap 0.326, 95% CI 0.214-0.434); Section 3.3 equation
    e_fed_fdepot <- 0.838; label("Fed-state multiplicative factor on relative bioavailability (power-form base)")      # Table 2 row 'F1 (FED)' = 0.838 (RSE 4.0%, bootstrap 0.838, 95% CI 0.773-0.903); Section 3.3 'F1 (FED) = 0.838'

    # ------------------------------------------------------------------
    # Between-subject variability. Table 2 reports omega^2 VARIANCES on the
    # log-normal exponential-IIV scale; the table note pins the scale
    # unambiguously: 'CV (%) for between-subject variability = SQRT (EXP
    # (omega^2)-1) x 100'. Each variance below reproduces the paper's own
    # CV% column to the printed precision, e.g. sqrt(exp(0.0662)-1) =
    # 26.2% and sqrt(exp(2.34)-1) = 306.3%.
    #
    # Two 2x2 correlation blocks were estimated (Section 3.3: 'Covariances
    # were set for the between-subject variability of CL/F-V2/F and
    # Q/F-V3/F'). Implied correlations are 0.0483/sqrt(0.0662 x 0.0638) =
    # 0.743 and 0.307/sqrt(0.411 x 0.259) = 0.941; both blocks are positive
    # definite as reported.
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(0.0662,
                        0.0483, 0.0638)  # Table 2 'omega CL/F^2' = 0.0662 (CV 26.2%, shrinkage 6.0%), 'omega CL/F-V2/F^2' = 0.0483, 'omega V2/F^2' = 0.0638 (CV 25.7%, shrinkage 22.0%)
    etalq  + etalvp ~ c(0.411,
                        0.307,  0.259)   # Table 2 'omega Q/F^2'  = 0.411  (CV 71.3%, shrinkage 14.5%), 'omega Q/F-V3/F^2'  = 0.307,  'omega V3/F^2'  = 0.259  (CV 54.4%, shrinkage 13.5%)

    etalka ~ 2.34  # Table 2 'omega Ka^2' = 2.34 (CV 306.3%, shrinkage 13.5%) -> log-scale SD 1.53

    # No eta on the lag time or on bioavailability: Table 2 reports
    # 'omega ALAG1^2 = 0, FIX' and 'omega F1 (FED)^2 = 0, FIX', i.e. both
    # IIV terms were fixed to zero rather than estimated. A zero-variance
    # eta is not representable in nlmixr2, so those etas are omitted, which
    # is the exact equivalent of the NONMEM encoding.

    # ------------------------------------------------------------------
    # Residual error. Section 3.3: 'Residual variability was described by
    # an exponential error model'. NONMEM's exponential residual
    # Y = F x EXP(EPS(1)) is nlmixr2's lnorm() residual, whose argument is
    # the log-scale SD. Table 2 reports sigma^2 = 0.136 as a VARIANCE - the
    # table note 'CV (%) for residual variability = sigma x 100' confirms
    # the scale, since sqrt(0.136) = 0.369 reproduces the printed 36.9%.
    # ------------------------------------------------------------------
    expSd <- sqrt(0.136); label("Exponential (log-normal) residual error, log-scale SD") # Table 2 row 'sigma^2 (exponential error)' = 0.136 (CV 36.9%, RSE 6.2%, shrinkage 6.5%) -> expSd = 0.369
  })

  model({
    # Individual PK parameters. The CL/F covariate model is Section 3.3
    # verbatim; V2/F, Q/F and V3/F carry no covariates (Section 3.3
    # identifies body weight and eGFR on CL/F and fasted/fed status on F as
    # the only retained covariates).
    cl <- exp(lcl + etalcl) * (WT / 60.6)^e_wt_cl * (CRCL / 98.6)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    ka     <- exp(lka + etalka)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot) * e_fed_fdepot^FED

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with first-order absorption from a depot
    # and a lag time on the absorption input.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Plasma febuxostat concentration. Doses are in mg and volumes in L, so
    # central / vc is mg/L; the factor 1000 converts to the ng/mL used
    # throughout the paper (Table 3, Figure 1, Figure 5).
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}

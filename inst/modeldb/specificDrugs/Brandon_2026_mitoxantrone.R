Brandon_2026_mitoxantrone <- function() {
  description <- "Two-compartment intravenous population PK model for mitoxantrone in children with acute myeloid leukaemia (Brandon 2026), developed from 282 plasma concentrations (45.7% below the 5 ng/mL assay LLOQ, handled by the likelihood-based M3 method) in 44 patients aged 0.9-17 years receiving 1-h infusions of 12 mg/m2/day, or 0.4 mg/kg/day for infants below 12 months, at or under 10 kg, or below 0.5 m2 body surface area. CL and Q carry a body-weight allometric exponent of 0.75 and Vc and Vp an exponent of 1, both referenced to the cohort median weight of 27.5 kg. Vc is fixed to the paediatric literature value of 23.2 L from O'Brien 2010 to stabilise the fit; interindividual variability is a correlated block on CL and Vp. No other covariate was retained in the final model."
  reference <- "Brandon AM, Huisman-Siebinga H, Barnett S, Wetherell P, Kearns P, Gibson B, Heaney N, Smith O, Baruchel A, Petit A, Moore A, Ogungbenro K, Huitema ADR, Veal GJ. Population pharmacokinetics and dose-response relationships of mitoxantrone in children with acute myeloid leukaemia. Br J Clin Pharmacol. 2026;92(6):1760-1770. doi:10.1002/bcp.70436"
  vignette <- "Brandon_2026_mitoxantrone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "mitoxantrone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mitoxantrone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "The only covariate retained in the final model. Drives fixed allometric scaling on CL and Q (exponent 0.75) and on Vc and Vp (exponent 1), referenced to the dataset median body weight of 27.5 kg (Brandon 2026 equation 1 and Methods section 2.3: 'allometric scaling of body weight was implemented in the base model, with parameters scaled to the dataset median body weight of 27.5 kg'). Cohort mean 33.3 kg, median 27.5 kg, range 9.5-69.5 kg (Table 1). The supplement's final NONMEM control stream (Supporting Information S2) writes the same scaling as LTVCL = LOG(THETA(1)) + THETA(5)*LOG(WT/27.5) and TVV1 = THETA(2)*(WT/27.5)**THETA(6).",
      source_name = "WT"
    )
  )

  # Screened during covariate model building but NOT retained in the final
  # model. Brandon 2026 Results section 3.2: "No other clinically relevant
  # covariate effect was found to be significant. Including BSA as a covariate
  # instead of allometric scaling did not improve the model, nor did estimating
  # the allometric scaling exponents." Documented here to preserve the
  # provenance of the covariate screen without carrying convention warnings.
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units = "m2",
      type = "continuous",
      notes = "Tested as a power model scaled to the population median, as an alternative to body-weight allometry (Methods 2.3). Not retained: 'Including BSA as a covariate instead of allometric scaling did not improve the model' (Results 3.2). BSA nevertheless determines the administered dose for the 12 mg/m2/day regimen. Cohort mean 1.09 m2, median 0.98, range 0.42-1.84 (Table 1)."
    ),
    AGE = list(
      description = "Age at trial entry",
      units = "years",
      type = "continuous",
      notes = "Tested as linear, power and maturation (sigmoidal Emax / Hill) functions (Methods 2.3). Not retained; no point estimate for an age effect is reported anywhere in the paper. Cohort mean 9.6 years, median 9.7, range 0.9-17.0 (Table 1). Age enters the study only through the dosing-regimen eligibility rule (below 12 months qualifies for mg/kg dosing), not through the PK model."
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Tested as linear and power models on clearance under the 'liver function' heading (Methods 2.3), motivated by the reported three-fold AUC increase in adults with bilirubin above 58.1 umol/L (Introduction). Not retained. Source column BIL in the supplement's $INPUT record. Cohort mean 7.4 umol/L, median 6.0, range 3-19 (Table 1) -- entirely within the normal range, so the cohort had no power to detect the adult hepatic-impairment effect."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested as linear and power models on clearance under the 'liver function' heading (Methods 2.3). Not retained. Source column ALT in the supplement's $INPUT record. Cohort mean 53.0 U/L, median 42.5, range 6-184 (Table 1); imputed as the population median 42.5 U/L for the two patients with no value (Results 3.1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Carried in the analysis dataset (supplement $INPUT column SEX, coded 0 = male / 1 = female per the Figure S1 legend, which matches the canonical SEXF orientation without value transformation) and plotted against the other covariates in Supporting Information Figure S1, but not reported as tested on any PK parameter and not retained. 17 of 44 patients (38.6%) were female (Table 1)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Supplement $INPUT column HT; plotted in Figure S1 but not tested as a PK covariate. Missing for 21 of 44 patients and estimated from UK-WHO growth charts using sex, age and weight (Results 3.1). Used to derive BSA and eGFR rather than entering the model directly. Cohort mean 134.8 cm, median 134.2, range 71.0-179.6 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Supplement $INPUT column SCR (a registered alias of the canonical CREAT); plotted in Figure S1 but not reported as tested on any PK parameter. Imputed as the population median 36.0 umol/L for the two patients with no value (Results 3.1). Cohort mean 37.7 umol/L, median 36.0, range 14-74 (Table 1). Mitoxantrone is primarily eliminated by biliary excretion (Introduction), so renal markers were not the physiologically motivated covariates here."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (revised Schwartz equation), BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Supplement $INPUT column GFR (a registered alias of the canonical CRCL, which covers BSA-normalized creatinine-based eGFR); plotted in Figure S1 but not reported as tested on any PK parameter. Missing for 19 patients and estimated from serum creatinine with known or estimated height using the revised Schwartz equation (Results 3.1). Cohort mean 138.1, median 133.7, range 88.7-229.7 mL/min/1.73 m^2 (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Supplement $INPUT column ALB; plotted in Figure S1 but not reported as tested on any PK parameter. Imputed as the population median 35.5 g/L for the two patients with no value (Results 3.1). Cohort mean 33.8 g/L, median 35.5, range 4-43 (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 1L,
    age_range = "0.9-17.0 years",
    age_median = "9.7 years",
    weight_range = "9.5-69.5 kg",
    weight_median = "27.5 kg",
    sex_female_pct = 38.6,
    race_ethnicity = "Not reported",
    disease_state = "Newly diagnosed acute myeloid leukaemia (40 patients, 90.9%), isolated myeloid sarcoma (3, 6.8%) and high-risk myelodysplastic syndrome (1, 2.3%); all under 18 years at trial entry",
    dose_range = "Mitoxantrone 12 mg/m2/day by 1-h intravenous infusion, once daily for up to 4 days. Infants below 12 months old, weighing 10 kg or less, or with a body surface area below 0.5 m2 instead received 0.4 mg/kg/day (3 of 44 patients). Administered dose 3.6-23 mg (median 12 mg); infusion duration 0.5-2.5 h (median 1.08 h).",
    regions = "United Kingdom, Ireland, France and Australia",
    n_observations = "282 plasma concentrations from 145 dosing events, 45.7% below the 5 ng/mL LLOQ (313 samples were originally collected; artificially high end-of-infusion samples were excluded)",
    samples_plasma = "Pre-dose; immediately after end of infusion on day 1; 0.5, 1, 2 and 6 h post-infusion on day 1; immediately before the day 2 infusion; and 48 and 72 h after the end of the final day's infusion. Mean 7.1 samples per patient (range 3-8).",
    co_medication = "All patients received non-chemotherapeutic concomitant medications for symptom and side-effect management; 16 of 44 (36.4%) received concomitant cytarabine (30-180 mg IV) within 7 days of starting mitoxantrone.",
    notes = "Demographics from Brandon 2026 Table 1. Data came from a single phase III paediatric AML trial (ISRCTN12389567; EudraCT 2014-005066-30). Assay was HPLC with photo-diode-array detection, calibration range 5-1000 ng/mL and LLOQ 5 ng/mL. Estimation was by SAEM with a separate importance-sampling evaluation step, using the likelihood-based M3 method for the below-LLOQ observations. Serum creatinine, ALT and serum albumin were missing for two patients each and imputed as the population medians 36.0 umol/L, 42.5 U/L and 35.5 g/L."
  )

  ini({
    # -------- Structural disposition (Brandon 2026 Table 2; the table
    # footnote states "Parameter values are based on the median body weight
    # of 27.5 kg"). The supplement's final NONMEM control stream (Supporting
    # Information S2) carries the identical $THETA record, confirming these
    # are the final estimates and not initial values:
    #   39.1     ; CL ; L/h ;(1)
    #   23.2 FIX ; V1 ; L   ;(2)
    #   27.6     ; Q  ; L/h ;(3)
    #   85.9     ; V2 ; L   ;(4)
    # --------
    lcl <- log(39.1); label("Clearance at the reference body weight of 27.5 kg (L/h)")                        # Table 2: CL = 39.1 L/h, RSE 9.61%
    lq  <- log(27.6); label("Intercompartmental clearance at the reference body weight of 27.5 kg (L/h)")     # Table 2: Q = 27.6 L/h, RSE 25.0%
    lvp <- log(85.9); label("Peripheral volume of distribution at the reference body weight of 27.5 kg (L)")  # Table 2: V2 = 85.9 L, RSE 22.6%

    # Vc was not estimated. Brandon 2026 Results 3.2: "V1 was fixed to
    # 23.2 L to improve model stability", the value taken from the
    # comparable paediatric cohort of O'Brien et al. 2010 (Discussion:
    # "In the presented model, V1 was fixed to the value reported by
    # O'Brien et al, which was based on a comparable patient population
    # to ours"). Marked FIX in the supplement's $THETA record.
    lvc <- fixed(log(23.2)); label("Central volume of distribution at the reference body weight of 27.5 kg, literature value taken from O'Brien 2010 (L)")  # Table 2: V1 = 23.2 L fixed

    # -------- Allometric exponents, fixed at the standard 3/4 and 1 power
    # model. Brandon 2026 Results 3.2: "allometric scaling of body weight
    # was included in the base model, using fixed exponents of 0.75 and 1.0
    # for clearances and volumes, respectively", and Methods 2.3 records
    # that estimating them was tested and rejected. Both appear as FIX in
    # the supplement's $THETA record ((5) WTCLQ and (6) WTV), and the
    # control stream applies THETA(5) to both CL and Q and THETA(6) to both
    # V1 and V2 -- hence the shared-exponent naming. --------
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless)")   # Table 2: 'WT effect on CL, Q' = 0.75 fixed
    e_wt_vc_vp <- fixed(1);    label("Shared allometric exponent on Vc and Vp (unitless)")  # Table 2: 'WT effect on V1, V2' = 1 fixed

    # -------- Interindividual variability, an estimated 2x2 block on CL and
    # Vp. Brandon 2026 Results 3.2: "The final model comprised two
    # compartments, parameterized as CL, V1, Q and V2, with IIV on CL and
    # V2." Equation 2 defines the exponential (log-normal) IIV model.
    #
    # Table 2 prints these as CV% with the covariance on the variance scale;
    # the supplement's $OMEGA BLOCK(2) record gives the estimated variances
    # directly and is used here verbatim:
    #   $OMEGA BLOCK(2)
    #   0.344        ; IIV CL ;;(ETA1)
    #   -0.426 1.02  ; IIV V2 ;;(ETA2)
    #
    # The two settle each other on the exact log-normal relation
    # omega^2 = log(CV^2 + 1): log(1 + 0.640^2) = 0.343 and
    # log(1 + 1.33^2) = 1.018, reproducing the printed 64.0% and 133% to
    # the precision the variances are given to. (The commonly-used
    # approximation CV% = 100*sqrt(omega^2) would instead imply variances
    # of 0.410 and 1.769, which the control stream contradicts.) The
    # implied correlation is -0.426 / sqrt(0.344 * 1.02) = -0.72.
    # --------
    etalcl + etalvp ~ c(0.344,
                        -0.426, 1.02)  # Supplement S2 $OMEGA BLOCK(2); Table 2 'IIV CL' 64.0% CV (RSE 26.2%, shrinkage 10.65%), 'Covariance of IIV CL, V2' -0.426 (RSE 43.0%), 'IIV V2' 133% CV (RSE 64.4%, shrinkage 34.11%)

    # -------- Residual unexplained variability. Brandon 2026 Results 3.2:
    # "Residual unexplained variability was included as a combined
    # (additive + proportional) error model, with the additive error fixed
    # to 2.5 ng/mL, reflecting half the assay LLOQ." The supplement's
    # $ERROR block builds exactly the nlmixr2 add()+prop() variance sum:
    #   W = SQRT((THETA(7)*IPRED)**2+THETA(8)**2)
    #   Y = IPRED + W*EPS(1)                       with $SIGMA 1 FIX
    # so THETA(7) is the proportional SD and THETA(8) the additive SD, both
    # on the ng/mL observation scale. --------
    propSd <- 0.382;      label("Proportional residual error (fraction)")  # Table 2: proportional error = 0.382, RSE 8.91%
    addSd  <- fixed(2.5); label("Additive residual error, set to half the 5 ng/mL assay LLOQ (ng/mL)")  # Table 2: additive error = 2.5 fixed
  })

  model({
    # 1. Individual PK parameters. Fixed allometric size scaling to the
    # 27.5 kg cohort-median reference (Brandon 2026 equation 1:
    # TVP = theta * (WT_i / WT_med)^k). IIV is on CL and Vp only; Vc and Q
    # carry no random effect (Results 3.2).
    cl <- exp(lcl + etalcl) * (WT / 27.5)^e_wt_cl_q
    vc <- exp(lvc)          * (WT / 27.5)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 27.5)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 27.5)^e_wt_vc_vp

    # 2. Micro-constants for the two-compartment linear disposition
    # (NONMEM ADVAN3 TRANS4 in the supplement's $SUBROUTINES record).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system. Mitoxantrone is given as a 1-h intravenous infusion, so
    # the dose enters the central compartment directly and there is no
    # absorption or bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Observation. Amounts are in mg and volumes in L, giving mg/L
    # (= ug/mL); multiply by 1000 to report ng/mL as the paper and the assay
    # do. This reproduces the supplement's scaling record S1 = V1/1000,
    # which divides the central amount by V1/1000 for the same effect.
    # Unit check against the paper's own exposures: a typical 27.5 kg
    # patient receiving the median 12 mg dose has AUCinf = Dose / CL =
    # 12 / 39.1 = 0.307 mg*h/L = 307 ug*h/L, against the reported
    # mg/m2-group mean of 317 +/- 184 ug*h/L (Results 3.3). Note 1 ug/L is
    # 1 ng/mL, so the paper's ug*h/L exposures are on this same scale.
    Cc <- (central / vc) * 1000

    # Combined additive plus proportional residual error; nlmixr2 sums the
    # two on the variance scale, matching the supplement's
    # W = sqrt((prop*IPRED)^2 + add^2).
    Cc ~ add(addSd) + prop(propSd)
  })
}

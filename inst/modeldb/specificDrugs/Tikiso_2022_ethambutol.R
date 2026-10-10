Tikiso_2022_ethambutol <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ethambutol in plasma in",
    "African children with tuberculosis, with and without HIV (pooled DATiC,",
    "DNDi super-boosted lopinavir and SHINE studies; Tikiso 2022).",
    "Absorption is a Savic transit chain (MTT 40.4 min, estimated chain",
    "length 4.82) feeding a first-order depot (ka 1.43 1/h); elimination is",
    "first order. Clearance (15.9 L/h), central volume (44.3 L),",
    "intercompartmental clearance (11.5 L/h) and peripheral volume (86.2 L)",
    "are allometrically scaled on paediatric fat-free mass (computed in the",
    "model with the Al-Sallami / Janmahasatian formula from weight, height,",
    "age and sex) against a 7.7 kg reference with fixed 0.75 / 1 exponents.",
    "Clearance matures with postmenstrual age through a sigmoid Emax",
    "function (PMA50 10.8 months, Hill 3.25). Bioavailability is 32% lower",
    "with lopinavir/ritonavir and declines linearly by 8.53% per year below",
    "an estimated age hinge of 3.16 years; absorption is 23.6% slower in the",
    "DNDi and SHINE studies than in DATiC. Between-subject variability on",
    "clearance and two-occasion between-occasion variability on",
    "bioavailability (inflated 1.37-fold for the unobserved pre-dose",
    "occasion), ka and MTT; combined proportional (17.7%) plus fixed",
    "additive (0.0168 mg/L) residual error. Doses are mg of ethambutol base",
    "(0.737 x ethambutol dihydrochloride)."
  )
  reference <- paste(
    "Tikiso T, McIlleron H, Abdelwahab MT, Bekker A, Hesseling A,",
    "Chabala C, Davies G, Zar HJ, Rabie H, Andrieux-Meyer I, Lee J,",
    "Wiesner L, Cotton MF, Denti P (2022).",
    "Population pharmacokinetics of ethambutol in African children: a",
    "pooled analysis. J Antimicrob Chemother 77:1949-1959.",
    "doi:10.1093/jac/dkac127",
    sep = " "
  )
  vignette <- "Tikiso_2022_ethambutol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ethambutol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight, an input to the paediatric fat-free mass formula.",
      units = "kg",
      type = "continuous",
      source_name = "WT",
      notes = paste(
        "Not used directly for scaling: FFM beat total body weight as the",
        "size descriptor (dOFV 12.1 vs 6.91, Results). Cohort median",
        "(range) 9.6 (3.9-34.5) kg, Table 1."
      )
    ),
    HT = list(
      description = "Body height, an input to the paediatric fat-free mass formula.",
      units = "cm",
      type = "continuous",
      source_name = "HT",
      notes = paste(
        "The supplementary control stream converts HT to metres",
        "(HTM = HT/100) before the Janmahasatian term."
      )
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male; selects the sex-specific constants of the fat-free mass formula.",
      units = "(binary)",
      type = "binary",
      source_name = "SEX",
      reference_category = "0 (male)",
      notes = paste(
        "The control stream codes SEX = 0 for female ('IF (SEX.EQ.0) THEN ;",
        "female'), so SEXF = 1 - SEX. No sex effect on any PK parameter",
        "other than through FFM."
      )
    ),
    AGE = list(
      description = "Postnatal age in years.",
      units = "years",
      type = "continuous",
      source_name = "AGE",
      notes = paste(
        "Used three ways, as in the control stream: (1) the age term of",
        "the Al-Sallami paediatric FFM multiplier, (2) postmenstrual age",
        "PMAGE = AGE + GESTATION/52 (years) for clearance maturation, and",
        "(3) the hockey-stick effect on bioavailability below the",
        "estimated 3.16-year hinge. Cohort median (range) 1.9",
        "(0.3-12.6) years, Table 1."
      )
    ),
    GA = list(
      description = "Gestational age at birth in weeks; added to postnatal age to form postmenstrual age for clearance maturation.",
      units = "weeks",
      type = "continuous",
      source_name = "GESTATION",
      notes = paste(
        "Methods: 'Gestational age was used to calculate PMAGE; when",
        "unavailable, a gestational age of 39 weeks was used.' Supply",
        "GA = 39 when unknown."
      )
    ),
    CONMED_LOPINAVIR = list(
      description = "1 = child receiving lopinavir/ritonavir-based ART, 0 = not.",
      units = "(binary)",
      type = "binary",
      source_name = "LPVRTV",
      reference_category = "0 (no lopinavir/ritonavir)",
      notes = paste(
        "Control stream: LPVR = 1 if LPVRTV > 0, then BIO multiplied by",
        "THETA(14). 95 of 103 HIV+ children were on LPV/r + ABC + 3TC",
        "(Table 1); 88% of them on super-boosted lopinavir (4:4",
        "lopinavir:ritonavir), 12% on adjusted-dose lopinavir/ritonavir",
        "three times daily. The authors could not separate a drug-drug",
        "interaction from HIV-related malabsorption or formulation."
      )
    ),
    STUDY_DNDI_SHINE = list(
      description = "1 = record from the DNDi super-boosted lopinavir study or the SHINE trial, 0 = the DATiC study (reference).",
      units = "(binary)",
      type = "binary",
      source_name = "STUDY",
      reference_category = "0 (DATiC)",
      notes = paste(
        "Control stream: STUDY_ABS = THETA(13) if STUDY > 1 (STUDY = 1",
        "DATiC, 2 DNDi, 3 SHINE), applied as KA * STUDY_ABS and",
        "MTT / STUDY_ABS. Derive STUDY_DNDI_SHINE = as.integer(STUDY > 1).",
        "Captures assay, sampling-schedule and formulation differences",
        "between studies (Discussion); set to 0 for a DATiC-like child."
      )
    ),
    OCC = list(
      description = paste(
        "Occasion index for between-occasion variability. 1 = the",
        "unobserved dose taken the day before the PK visit (producing the",
        "pre-dose sample), 2 = the observed dose on the PK visit day."
      ),
      units = "(count)",
      type = "categorical",
      source_name = "OCC",
      notes = paste(
        "The control stream selects ETA(9/11/13/15) on OCC = 1 and",
        "ETA(10/12/14/16) on OCC = 2 with $OMEGA BLOCK(1) SAME, and scales",
        "the OCC = 1 bioavailability deviate by THETA(18) = 1.37 (Table 2",
        "'Scaling of BOV in F for unobserved doses'). Records with OCC = 0",
        "carry no between-occasion variability."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 188L,
    n_studies = 3L,
    age_range = "0.3-12.6 years",
    age_median = "1.9 years",
    weight_range = "3.9-34.5 kg",
    weight_median = "9.6 kg",
    ffm_median = "7.7 kg (range 3.3-29.2)",
    sex_female_pct = 54.8,
    race_ethnicity = "Black African children (South Africa, Malawi, Zambia); race not tabulated.",
    disease_state = paste(
      "Children with tuberculosis on first-line treatment; 103 (54.8%)",
      "HIV-positive on ART, 92% of them on lopinavir/ritonavir."
    ),
    dose_range = paste(
      "Ethambutol dihydrochloride 15-25 mg/kg orally once daily per WHO",
      "guidelines (median 20.2, range 12.5-25.0 mg/kg); whole, crushed,",
      "syringe or nasogastric-tube administration."
    ),
    regions = "South Africa, Malawi (DATiC), Zambia (SHINE)",
    notes = paste(
      "Table 1: DATiC n = 79 (368 samples), DNDi n = 84 (471), SHINE",
      "n = 25 (173); 85 of 188 male. 1012 plasma concentrations, 121 (12%)",
      "below the 0.0844 mg/L LLOQ handled by the M6 method. Steady-state",
      "intensive sampling up to 8-12 h post-dose."
    )
  )

  ini({
    # Disposition. Table 2 typical values refer to 'a child weighing 10 kg
    # (at fully mature clearance) not co-treated with lopinavir/ritonavir'
    # (footnote b); the control stream scales on FFM with TVFFM = 7.7 kg,
    # the cohort median FFM (Table 1).
    # Values are apparent and relative to ethambutol BASE dose: the dose of
    # ethambutol dihydrochloride was multiplied by 0.737 (Methods).
    lcl <- log(15.9)
    label("Apparent clearance CL/F at FFM 7.7 kg, fully mature (L/h)") # Table 2 CL 15.9 L/h (14.8-17.3)
    lvc <- log(44.3)
    label("Apparent central volume Vc/F at FFM 7.7 kg (L)") # Table 2 Central Vd 44.3 L (37.2-51.1)
    lq <- log(11.5)
    label("Apparent intercompartmental clearance Q/F at FFM 7.7 kg (L/h)") # Table 2 Inter-compartmental clearance 11.5 L/h (10.1-13.2)
    lvp <- log(86.2)
    label("Apparent peripheral volume Vp/F at FFM 7.7 kg (L)") # Table 2 Peripheral Vd 86.2 L (73.6-102)

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent on fat-free mass for CL/F and Q/F (unitless)") # Methods: exponents fixed to 3/4 for clearance; control stream ALLMCL_FFM_CH = (FFM_CH/TVFFM)**0.75
    e_ffm_vc <- fixed(1)
    label("Allometric exponent on fat-free mass for Vc/F and Vp/F (unitless)") # Methods: exponents fixed to 1 for volume; control stream ALLMV_FFM_CH = (FFM_CH/TVFFM)

    # Clearance maturation (sigmoid Emax on postmenstrual age, not normalised).
    ltm50_cl <- log(10.8)
    label("Postmenstrual age at 50% clearance maturation PMAGE50 (months)") # Table 2 PMAGE50 10.8 months (9.66-11.7)
    lhill_cl <- log(3.25)
    label("Shape (Hill) coefficient of the clearance maturation function (unitless)") # Table 2 gamma-shape of maturation function 3.25 (2.76-3.77)

    # Absorption: Savic transit chain feeding a first-order depot.
    lka <- log(1.43)
    label("First-order absorption rate constant ka in DATiC (1/h)") # Table 2 ka 1.43 1/h (1.10-1.84)
    lmtt <- log(40.4 / 60)
    label("Mean absorption transit time MTT in DATiC (h)") # Table 2 MTT 40.4 min (35.3-46.1) = 0.6733 h
    lntr <- log(4.82)
    label("Number of absorption transit compartments (unitless)") # Table 2 n 4.82 (3.70-6.32)
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability F (unitless)") # Table 2 F 1 Fixed

    # Covariate effects.
    e_study_dndi_shine_ka <- -0.236
    label("Fractional change in speed of absorption (ka up, MTT down) in DNDi and SHINE vs DATiC (unitless)") # Table 2 Change in speed of absorption in DNDi and SHINE -23.6% (-31.6 to -14.4); footnote c Ka = theta_Ka x theta_change, MTT = theta_MTT / theta_change
    e_lopinavir_fdepot <- -0.320
    label("Fractional change in F with lopinavir/ritonavir (unitless)") # Table 2 Change in F when on LPV/r -32.0% (-38.9 to -23.8)
    lage_hinge <- log(3.16)
    label("Age hinge above which F no longer depends on age (years)") # Table 2 Breakpoint for age effect on F 3.16 years (2.18-4.14)
    e_age_fdepot <- 0.0853
    label("Fractional change in F per year of age below the age hinge (1/year)") # Table 2 Age on F, fractional change +0.0853 per year (+0.0463 to +0.130)

    # Random effects. Table 2 reports each variability as a percentage that
    # is the omega standard deviation on the log scale: the supplementary
    # control stream's $OMEGA initials reproduce the tabulated percentages
    # as sqrt(omega^2) (BOV ka 0.403 -> 63.5% vs 63.7%; BOV MTT 0.233 ->
    # 48.3% vs 48.5%), whereas sqrt(exp(omega^2) - 1) gives 70.4% and 51.0%.
    # BSV on V, ka, F, Q and Vp and BOV on CL are '0 FIX' in the control
    # stream, so those etas are absent.
    etalcl ~ 0.013456 # Table 2 CL BSV 11.6% (7.12-15.3); 0.116^2 = 0.013456

    # BOV on F. The control stream scales the OCC = 1 deviate,
    # 'IF (OCC==1) SCALE_BOVBIO=THETA(18)', so the unobserved-dose occasion
    # has variance 1.37^2 * 0.331^2 = 0.205635. This is an exact
    # re-expression of the published parameterisation.
    etaiov_fdepot_1 ~ fixed(0.205635) # Table 2 Scaling of BOV in F for unobserved doses 1.37 (1.12-1.65); 1.37^2 * 0.331^2
    etaiov_fdepot_2 ~ 0.109561 # Table 2 F BOV 33.1% (30.1-36.1); 0.331^2 = 0.109561
    etaiov_ka_1 ~ 0.405769 # Table 2 ka BOV 63.7% (51.3-77.5); 0.637^2 = 0.405769
    etaiov_ka_2 ~ fixed(0.405769) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 0.235225 # Table 2 MTT BOV 48.5% (41.2-56.9); 0.485^2 = 0.235225
    etaiov_mtt_2 ~ fixed(0.235225) # equal to occasion 1 per control stream $OMEGA BLOCK(1) SAME

    # Residual error: W = SQRT(ADD**2 + PROP**2) with ADD = LLOQ/5 + THETA(8),
    # THETA(8) FIX 0 and LLOQ = 0.084 mg/L.
    propSd <- 0.177
    label("Proportional residual error (fraction)") # Table 2 Proportional error 17.7% (16.5-19.1)
    addSd <- fixed(0.0168)
    label("Additive residual error (mg/L)") # Table 2 Additive error 16.8 ug/L Fixed (footnote d: 20% of LLOQ)
  })

  model({
    # 1. Occasion indicators for the two-occasion BOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2

    # 2. Paediatric fat-free mass (control stream $PK, Al-Sallami 2015 age
    # multiplier on the Janmahasatian 2005 adult formula), WT in kg, height
    # in metres, AGE in years. The constants select on sex.
    htm <- HT / 100
    ffm_alpha <- 0.88 + (1.11 - 0.88) * SEXF
    ffm_a50 <- 13.4 + (7.1 - 13.4) * SEXF
    ffm_gamma <- 12.7 + (1.1 - 12.7) * SEXF
    ffm_whsmax <- 42.92 + (37.99 - 42.92) * SEXF
    ffm_whs50 <- 30.93 + (35.98 - 30.93) * SEXF
    agegam <- AGE^ffm_gamma
    a50gam <- ffm_a50^ffm_gamma
    ffm <- ((agegam + ffm_alpha * a50gam) / (agegam + a50gam)) *
      ((ffm_whsmax * htm^2 * WT) / (ffm_whs50 * htm^2 + WT))

    # 3. Clearance maturation on postmenstrual age. The control stream works
    # in years (PMAGE = AGE + GESTATION/52); Table 2 reports PMAGE50 in
    # months, so PMA is converted to months here. The function is not
    # normalised: CL reaches its Table 2 value only at full maturation.
    pma <- (AGE + GA / 52) * 12
    tm50_cl <- exp(ltm50_cl)
    hill_cl <- exp(lhill_cl)
    mat_cl <- pma^hill_cl / (pma^hill_cl + tm50_cl^hill_cl)

    # 4. Bioavailability covariates. Hockey-stick on age: linear below the
    # hinge, flat above it (control stream THETA(17) FIX 0).
    age_hinge <- exp(lage_hinge)
    f_age <- 1 + e_age_fdepot * (min(AGE, age_hinge) - age_hinge)
    f_lpv <- 1 + e_lopinavir_fdepot * CONMED_LOPINAVIR
    f_study_abs <- 1 + e_study_dndi_shine_ka * STUDY_DNDI_SHINE

    # 5. Individual parameters.
    cl <- exp(lcl + etalcl) * (ffm / 7.7)^e_ffm_cl * mat_cl
    vc <- exp(lvc) * (ffm / 7.7)^e_ffm_vc
    q <- exp(lq) * (ffm / 7.7)^e_ffm_cl
    vp <- exp(lvp) * (ffm / 7.7)^e_ffm_vc
    ka <- exp(lka + iov_ka) * f_study_abs
    mtt <- exp(lmtt + iov_mtt) / f_study_abs
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot) * f_age * f_lpv

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 6. Savic transit-compartment input, written out as in the control
    # stream $DES:
    #   KTR = (NN+1)/MTT; PIZZA = LOG(BIO*PD*KTR) - GAMLN(NN+1)
    #   DADT(1) = EXP(PIZZA + NN*LOG(KTR*TEMPO) - KTR*TEMPO) - KA*A(1)
    # with PD the latest dose amount and TEMPO the time after it. podo() and
    # tad() supply PD and TEMPO; rxode2's transit() built-in is not used
    # because combined with f(depot) <- 0 it yields a zero input rate.
    tdos <- tad(depot)
    ktr <- (ntr + 1) / mtt
    ktt <- ktr * tdos
    trin <- 0
    if (ktt > 0) {
      trin <- exp(log(fdepot * podo(depot) * ktr) - lgamma(ntr + 1) +
                    ntr * log(ktt) - ktt)
    }

    # 7. ODE system.
    d/dt(depot) <- trin - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # The dose is delivered entirely through the transit density, so the
    # ordinary bolus into the depot is suppressed (control stream F1 = 0).
    f(depot) <- 0

    # 8. Observation and residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

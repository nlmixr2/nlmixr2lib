Willmann_2019_moxifloxacin <- function() {
  description <- paste(
    "Three-compartment population PK model for moxifloxacin in children and",
    "adolescents aged 3 months to <18 years (Willmann 2019), with first-order",
    "absorption from an oral depot, logit-scale absolute bioavailability, and",
    "a-priori (not estimated) allometric body-weight scaling of all clearance",
    "and volume parameters with the canonical exponents 0.75 and 1.0. The",
    "structural parameters are reported per-kg (CL, Q in L/h/kg^0.75; V in",
    "L/kg) rather than normalized to a reference weight, so no centering",
    "weight enters the model. Interindividual variability was retained on CL",
    "and on the central volume only; the source held the absorption-rate and",
    "bioavailability etas at zero. Residual variability is proportional with",
    "three separate magnitudes, selected by development phase and route:",
    "phase I intravenous, phase III intravenous and phase III oral. The",
    "covariate screen (age, serum creatinine, estimated glomerular filtration",
    "rate, study and sex) retained nothing beyond the a-priori weight scaling.",
    "The companion whole-body PBPK model of the same paper is a PK-Sim / MoBi",
    "platform model whose physiological parameters come from the vendor's",
    "internal databases and are not printed in the paper; it is therefore not",
    "reproducible as an rxode2 model and is not encoded here.",
    sep = " "
  )
  reference <- paste(
    "Willmann S, Frei M, Sutter G, Coboeken K, Wendl T, Eissing T,",
    "Lippert J, Stass H (2019).",
    "Application of physiologically-based and population pharmacokinetic",
    "modeling for dose finding and confirmation during the pediatric",
    "development of moxifloxacin.",
    "CPT: Pharmacometrics & Systems Pharmacology 8(9):654-663.",
    "doi:10.1002/psp4.12446.",
    "Final parameter estimates are from Supplementary Table S2 (final",
    "population PK model O6/run021); the model structure, the a-priori",
    "allometric exponents, the logit bioavailability parameterization, the",
    "zero-fixed absorption-rate and bioavailability etas and the three-way",
    "residual-error stratification are from the NONMEM control stream",
    "reproduced verbatim in Supplementary Material S1 (Model Code).",
    sep = " "
  )
  vignette <- "Willmann_2019_moxifloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "A-priori (not estimated) allometric scaling of every disposition",
        "parameter, applied with the commonly accepted exponents 0.75 for all",
        "clearances and 1.0 for all volumes (Willmann 2019 Methods,",
        "'Pediatric popPK model'; Supplementary Methods 'Model selection",
        "criteria'). The control stream writes the scaling without a",
        "reference weight -- TCL = THETA(2)*WGHT**0.75, TV2 = THETA(3)*WGHT,",
        "and likewise for Q3/V3/Q4/V4 -- so the reported typical values are",
        "per-kg quantities (L/h/kg^0.75 and L/kg) and there is no centering",
        "weight to supply. Source column WGHT. The covariate analysis found",
        "no further significant effect on top of this weight relationship.",
        sep = " "
      ),
      source_name = "WGHT"
    ),
    STUDY_PHASE3 = list(
      description = "Phase III study indicator (pooled phase I / phase III analysis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = phase I study 11826 (single-dose intravenous, N = 31)",
      notes = paste(
        "1 = record from the phase III complicated intra-abdominal infection",
        "study 11643, 0 = record from the phase I study 11826. Selects the",
        "residual-error magnitude only; it touches no structural parameter, so",
        "the typical-value prediction is identical for either setting",
        "(Willmann 2019 Supplementary Table S2 residual-error rows, and the",
        "control stream $ERROR branch IF(STUD.EQ.11826)). When simulating,",
        "use STUDY_PHASE3 = 0 to reproduce the rich phase I profiles and",
        "STUDY_PHASE3 = 1 to reproduce the phase III scatter. Source column",
        "STUD, an integer study number (11826 / 11643); derive",
        "STUDY_PHASE3 = as.integer(STUD == 11643).",
        sep = " "
      ),
      source_name = "STUD"
    ),
    ROUTE_ORAL = list(
      description = "Oral administration indicator (reference = intravenous infusion)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = intravenous infusion",
      notes = paste(
        "1 = observation record following oral moxifloxacin, 0 = intravenous",
        "infusion. Like STUDY_PHASE3 this selects the residual-error",
        "magnitude only; the structural route difference is carried by the",
        "dose record's target compartment (depot for oral, central for",
        "intravenous), so the covariate is read by the error model alone. In",
        "the source $ERROR block the oral branch is written last and",
        "therefore overrides the study branch, which is reproduced here by",
        "nesting the STUDY_PHASE3 selection inside (1 - ROUTE_ORAL). Only",
        "the phase III study contributed oral records (28 subjects switched",
        "to oral treatment), so the phase I oral cell is empty in the source",
        "data. Source column ROUT, coded 1 = oral and 2 = intravenous;",
        "derive ROUTE_ORAL = as.integer(ROUT == 1).",
        sep = " "
      ),
      source_name = "ROUT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened by stepwise forward inclusion (P < 0.01) / backward",
        "elimination (P < 0.001) and not retained: 'the covariate analysis",
        "did not yield any significant effect on top of the relation between",
        "body weight and the disposition parameters' (Willmann 2019 Results,",
        "'PopPK model in children'). No point estimate is reported, so the",
        "effect cannot be carried.",
        sep = " "
      ),
      source_name = "AGE"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "not reported",
      type = "continuous",
      notes = "Screened in the covariate analysis and not retained (Willmann 2019 Methods, 'Pediatric popPK model'; Results, 'PopPK model in children'). Source column SCRE; the unit is not printed in the paper or the supplement.",
      source_name = "SCRE"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate",
      units = "not reported",
      type = "continuous",
      notes = "Screened in the covariate analysis and not retained (Willmann 2019 Methods, 'Pediatric popPK model'). Source column EGFR; neither the estimating equation nor the unit is printed in the paper or the supplement.",
      source_name = "EGFR"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the covariate analysis and not retained (Willmann 2019 Methods, 'Pediatric popPK model'). Source column SEX.",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 186,
    n_studies = 2,
    age_range = "3 months to <18 years",
    disease_state = paste(
      "Phase III: complicated intra-abdominal infection.",
      "Phase I: pediatric patients receiving a single intravenous dose.",
      sep = " "
    ),
    dose_range = paste(
      "Phase I: single intravenous dose 5-10 mg/kg.",
      "Phase III: multiple intravenous doses followed by oral doses --",
      "400 mg once daily for 12 to <18 years and >=45 kg;",
      "4 mg/kg twice daily for 12 to <18 years and <45 kg and for",
      "6 to <12 years; 5 mg/kg twice daily for 2 to <6 years;",
      "6 mg/kg twice daily for 3 months to <2 years.",
      sep = " "
    ),
    notes = paste(
      "186 pediatric subjects contributing 1,562 moxifloxacin plasma",
      "concentrations were available for model development; 33 concentrations",
      "were excluded as influential outliers and 14 concentrations from 12",
      "subjects were retained as noninfluential outliers (Willmann 2019",
      "Results, 'PopPK model in children'). Pooled from a phase I study",
      "(study 11826, N = 31; 2 subjects >=12 years, 29 subjects <12 years)",
      "and a phase III study in complicated intra-abdominal infection (study",
      "11643, N = 451 randomized of whom 301 were exposed to moxifloxacin and",
      "155 contributed PK samples; 98 subjects >=12 years, 57 subjects",
      "<12 years); 28 phase III subjects switched to oral treatment",
      "(Willmann 2019 Methods, 'Data sources'; Results, 'Demographic",
      "analysis'). Young children were underrepresented: only 11 subjects",
      "younger than 3 years contributed PK information (10 from phase I, 1",
      "from phase III), i.e. 5.9% of the study population. The lower limit of",
      "quantification was 10 ug/L for moxifloxacin and values below it were",
      "excluded (Supplementary Methods, 'Bioanalytical methods'). Baseline",
      "weight, sex and race distributions are not tabulated in the paper or",
      "the supplement.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Willmann 2019 Supplementary Table S2, final
    # population PK model O6/run021. Values are per-kg because the control
    # stream scales by WGHT**0.75 / WGHT without a reference weight.
    lka <- log(0.537)
    label("Log of the first-order absorption rate constant from depot (1/h)")
    # Table S2 TKA = 0.537 1/h (RSE 26.6%, 95% CI 0.257-0.817)

    lcl <- log(0.45)
    label("Log of the weight-scaled clearance from central (L/h/kg^0.75)")
    # Table S2 TCL = 0.45 L/h/kg^0.75 (RSE 2.56%, 95% CI 0.427-0.473)

    lvc <- log(0.906)
    label("Log of the weight-scaled central volume of distribution (L/kg)")
    # Table S2 TV2 = 0.906 L/kg (RSE 5.61%, 95% CI 0.806-1.01)

    lq <- log(1.74)
    label("Log of the weight-scaled intercompartmental clearance to peripheral1 (L/h/kg^0.75)")
    # Table S2 TQ3 = 1.74 L/h/kg^0.75 (RSE 15.1%, 95% CI 1.22-2.26)

    lvp <- log(0.732)
    label("Log of the weight-scaled first peripheral volume of distribution (L/kg)")
    # Table S2 TV3 = 0.732 L/kg (RSE 6.15%, 95% CI 0.644-0.82)

    lq2 <- log(0.0889)
    label("Log of the weight-scaled intercompartmental clearance to peripheral2 (L/h/kg^0.75)")
    # Table S2 TQ4 = 0.0889 L/h/kg^0.75 (RSE 13.6%, 95% CI 0.0652-0.113)

    lvp2 <- log(0.615)
    label("Log of the weight-scaled second peripheral volume of distribution (L/kg)")
    # Table S2 TV4 = 0.615 L/kg (RSE 16.0%, 95% CI 0.422-0.808)

    logitfdepot <- log(0.866 / (1 - 0.866))
    label("Logit of the absolute oral bioavailability (unitless; F = 0.866)")
    # Table S2 TF = 0.866 (RSE 4.26%, 95% CI 0.794-0.938). The control stream
    # estimates F on the logit scale:
    # PHI = LOG(THETA(7)/(1-THETA(7))); F1 = EXP(PHI+ETA(4))/(1+EXP(PHI+ETA(4))).
    # An absorption lag time ALAG1 = THETA(6) was carried in the control
    # stream but FIXED at 0 and is therefore structurally absent here; it is
    # not listed in Table S2.

    # A-priori allometric exponents. Willmann 2019 Methods, 'Pediatric popPK
    # model': 'all CL and volume parameters were a priori scaled by body
    # weight using an allometric model with commonly accepted scaling
    # coefficients of 0.75 (for CL) and 1.0 (for volumes)'. They were never
    # estimated (they are hardcoded in the control stream $PK block as
    # WGHT**0.75 and WGHT), hence fixed().
    e_wt_cl_q_q2 <- fixed(0.75)
    label("Allometric exponent of body weight shared by CL, Q and Q2 (unitless)")
    # Control stream $PK: TCL = THETA(2)*WGHT**0.75, TQ3 = THETA(4)*WGHT**0.75, TQ4 = THETA(8)*WGHT**0.75

    e_wt_vc_vp_vp2 <- fixed(1)
    label("Allometric exponent of body weight shared by Vc, Vp and Vp2 (unitless)")
    # Control stream $PK: TV2 = THETA(3)*WGHT, TV3 = THETA(5)*WGHT, TV4 = THETA(9)*WGHT

    # Interindividual variability -- Table S2 'Random effects', reported as
    # OMEGA variances on the log scale with the CV% given by
    # SQRT(EXP(OMEGA2)-1)*100. IIV on ka (OMEGA 1) and on the logit
    # bioavailability (OMEGA 4) was FIXED at 0 in the control stream and is
    # therefore omitted rather than carried as a zero-variance eta.
    etalcl ~ 0.113 # Table S2 CL omega^2 = 0.113 (RSE 16.7%, 95% CI 0.0760-0.150), CV 34.6%; eta-shrinkage 4.43%
    etalvc ~ 0.257 # Table S2 V2 omega^2 = 0.257 (RSE 20.4%, 95% CI 0.154-0.360), CV 54.1%; eta-shrinkage 16.5%

    # Proportional residual error, three magnitudes selected by development
    # phase and route. Table S2 reports variances; the SD entered here is
    # SQRT(SIGMA2), which is the CV% the table prints alongside each variance.
    propSdPh1Iv <- 0.152643
    label("Proportional residual error SD, phase I intravenous (fraction)")
    # Table S2 '11826, IV' sigma^2 = 0.0233 (RSE 20.8%, 95% CI 0.0138-0.0328), CV 15.3%; sqrt(0.0233) = 0.152643

    propSdPh3Iv <- 0.337639
    label("Proportional residual error SD, phase III intravenous (fraction)")
    # Table S2 '11643, IV' sigma^2 = 0.114 (RSE 10.7%, 95% CI 0.0901-0.138), CV 33.8%; sqrt(0.114) = 0.337639

    propSdPh3Oral <- 0.475395
    label("Proportional residual error SD, phase III oral (fraction)")
    # Table S2 '11643, PO' sigma^2 = 0.226 (RSE 16.9%, 95% CI 0.151-0.301), CV 47.5%; sqrt(0.226) = 0.475395
  })

  model({
    # A-priori allometric scaling on body weight, written exactly as the
    # control stream does it: a per-kg typical value multiplied by WT raised
    # to the fixed exponent, with no reference weight.
    cl <- exp(lcl + etalcl) * WT^e_wt_cl_q_q2
    vc <- exp(lvc + etalvc) * WT^e_wt_vc_vp_vp2
    q <- exp(lq) * WT^e_wt_cl_q_q2
    vp <- exp(lvp) * WT^e_wt_vc_vp_vp2
    q2 <- exp(lq2) * WT^e_wt_cl_q_q2
    vp2 <- exp(lvp2) * WT^e_wt_vc_vp_vp2
    ka <- exp(lka)
    fdepot <- expit(logitfdepot)

    # Micro-constants for the explicit three-compartment ODE system
    # (NONMEM ADVAN12 TRANS4: depot + central + two peripherals).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Absolute bioavailability applies to the oral depot only; intravenous
    # doses enter central directly and never read it.
    f(depot) <- fdepot

    Cc <- central / vc

    # Stratum-specific proportional residual error. The source $ERROR block
    # writes the oral branch last so that it overrides the study branch:
    #   Y = IPRED + W*EPS(1)                      ; phase III, intravenous
    #   IF(STUD.EQ.11826) Y = IPRED + W*EPS(2)    ; phase I, intravenous
    #   IF(ROUT.EQ.1)     Y = IPRED + W*EPS(3)    ; oral
    # which is reproduced by nesting the study selection inside the
    # non-oral branch.
    propSdCc <- ROUTE_ORAL * propSdPh3Oral +
      (1 - ROUTE_ORAL) *
        (STUDY_PHASE3 * propSdPh3Iv + (1 - STUDY_PHASE3) * propSdPh1Iv)

    Cc ~ prop(propSdCc)
  })
}

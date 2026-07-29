Abrantes_2017_moroctocog <- function() {
  description <- "Two-compartment population PK model for factor VIII activity (IU/dL) following intravenous administration of moroctocog alfa (B-domain-deleted recombinant FVIII, marketed as ReFacto, ReFacto AF and Xyntha) in patients with moderate to severe hemophilia A; pooled analysis of 754 patients across 13 clinical trials over 20 years (Abrantes 2017). The exogenous-drug component is added to a constant endogenous-baseline FVIII activity (severe-subpopulation typical value, 0.474 IU/dL; the paper's full model is a two-class mixture, see vignette deviations). Clearance and inter-compartmental clearance scale allometrically with body weight at theory-based exponent 0.75; central and peripheral volumes share an estimated allometric exponent 0.812. Clearance has a piecewise-linear age effect (increasing from birth to 1 year of age, then decreasing into adulthood; centered at 20 years), a +166% inhibitor (ADA_POS) effect, and a -34.7% study B1831090 effect. The peripheral volume is +88.4% larger in Black subjects. Bioavailability F carries multiplicative covariate effects for product (1.38x for Xyntha vs ReFacto), assay (-39.0% for OSA central, -14.6% for OSA local laboratory) following Abrantes 2017 Table 2 footnote g. Proportional residual error is 19.2% (CSA reference) and switches to 26.9% (+40.3%) for OSA-assayed samples."
  reference <- "Abrantes JA, Nielsen EI, Korth-Bradley J, Harnisch L, Jonsson S. Elucidation of Factor VIII Activity Pharmacokinetics: A Pooled Population Analysis in Patients With Hemophilia A Treated With Moroctocog Alfa. Clin Pharmacol Ther. 2017 Jul;102(1):113-121. doi:10.1002/cpt.716. PMID:28437834."
  vignette <- "Abrantes_2017_moroctocog"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling at reference 70 kg: fixed theory-based exponent 0.75 on CL and Q (Abrantes 2017 Table 2 footnotes b, d), shared estimated exponent 0.812 on V1 and V2 (Table 2 footnotes c, e). Cohort median 69 kg (range 3.0-134 kg; Table 1). The pooled analysis used total body weight; lean body weight and total body water were tested but weight was retained because predictions did not differ substantially and r^2 between WT and TBW was 0.963 (Results). Missing weights imputed by carrying values backward/forward within subject (Methods).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Piecewise-linear effect on CL with breakpoint at 1 year of age (Abrantes 2017 Table 2 footnote b). Reference 20 years (where age_eff_cl = 1). Cohort median 23 years, range 0.0027 years (1 day) to 73 years; 234 of 754 subjects (31%) were < 17 years at study entry (Table 1). Below 1 year: age_eff_cl = 1 + 0.149 * (AGE - 1) - (20 - 1) * (-0.00678) = 1.1288 + 0.149 * (AGE - 1); above 1 year: age_eff_cl = 1 + (-0.00678) * (AGE - 20). The two pieces meet at age_eff_cl = 1.1288 when AGE = 1. Mechanism hypothesised by Abrantes 2017 Discussion to reflect age-dependent von Willebrand factor levels (high at birth, decreasing through year 1, gradually increasing through childhood).",
      source_name        = "AGE"
    ),
    ADA_POS = list(
      description        = "FVIII inhibitor status (1 = inhibitor-positive at >= 0.6 Bethesda units/mL, 0 = negative). FVIII inhibitors are neutralizing alloantibodies to administered FVIII measured by the Bethesda functional clotting-inhibition assay; mapped onto the canonical ADA_POS column per the NAB-subset alias documented in inst/references/covariate-columns.md.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (FVIII inhibitor negative)",
      notes              = "Abrantes 2017 Table 1: 174 of 7363 observations (2.4%) from inhibitor-positive patients; 6138 (83%) negative; 1051 (14%) missing inhibitor status (imputed forward/backward). 92.5% of inhibitor-positive subjects were low-titer (< 5 BU/mL). The dichotomous status was statistically more informative than the continuous Bethesda titer (Discussion); the paper attributes this to between-laboratory variability of the inhibitor assay and to different inhibitor kinetics across types. Inhibitor status is per-occasion (titer measured per visit), although it tends to be stable over the study window for most subjects. Source column INH; canonical ADA_POS used because FVIII inhibitors are neutralizing anti-drug alloantibodies (the NAB alias of ADA_POS).",
      source_name        = "INH"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black race: White, Asian, or Other)",
      notes              = "Abrantes 2017 Table 1: 10 of 754 subjects (1.3%) self-identified as Black; 657 (87%) White; 58 (7.7%) Asian; 29 (3.8%) Other. The +88.4% V2 increase in Black subjects has wide uncertainty (95% CI 36.7-151%) reflecting the small Black-race n. The paper hypothesises higher von Willebrand factor levels in Black subjects as the mechanism (Discussion); clinical consequence considered limited. Source column RACE (multi-level: Asian / Black / Other / White); derive RACE_BLACK = as.integer(RACE == 'Black') with all non-Black races collapsed to the reference.",
      source_name        = "RACE"
    ),
    STUDY_B1831090 = list(
      description        = "Study B1831090 cohort indicator within the Abrantes 2017 13-study pooled analysis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other study in the 13-study pool: B1831003, B1831004, B1831015, B1831053, B1831054, B1831061, B1831066, B1831067, B1831068, B1831070, B1831071, B1831077)",
      notes              = "Abrantes 2017 Table 3: study B1831090 is a phase I bioequivalence study of ReFacto (n = 18, severe and moderately severe hemophilia A, rich sampling, CSA assay). The -34.7% CL reduction for B1831090 subjects relative to the other studies (Table 2 footnote b) is retained in the model so that the typical PK parameters describe the remaining 12 studies. Subject-level / time-fixed; derive once per subject from the trial identifier (STUDY_B1831090 = 1 if trial == 'B1831090', else 0). Source column STUD.",
      source_name        = "STUD"
    ),
    FORM_XYNTHA = list(
      description        = "Xyntha vs ReFacto / ReFacto AF drug-product indicator (moroctocog alfa is the same active moiety in all three; the products differ in potency calibration)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ReFacto or ReFacto AF; CSA-calibrated potency)",
      notes              = "Abrantes 2017 Methods + Table 2 footnote g. Xyntha is OSA-calibrated and ReFacto/ReFacto AF are CSA-calibrated; the potency calibration difference is encoded as a +38.0% bioavailability multiplier (1.38^FORM_XYNTHA, equivalent to 1 + 0.38 * FORM_XYNTHA when binary) on the dose so that the model-predicted FVIII activity is in CSA-reference units (Methods). Per-observation when subjects could switch between products across trials (e.g., study B1831066 enrolled ReFacto and ReFacto AF in the same protocol); subject-level otherwise. Source column PROD.",
      source_name        = "PROD"
    ),
    ASSAY_OSA = list(
      description        = "One-stage activated partial thromboplastin time clotting assay (OSA) indicator (vs chromogenic substrate assay, CSA, as the reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CSA reference)",
      notes              = "Abrantes 2017 Table 2 footnote g + Methods. Switches both the bioavailability multiplier (F decreases by 39.0% for OSA-measured observations) and the proportional residual error magnitude (+40.3% for OSA). The OSA-CSA bias is well documented in the literature for B-domain-deleted FVIII products (Discussion). Per-observation indicator: an FVIII activity sample is assayed by one method; pooled datasets across studies may mix methods. Source column METH1.",
      source_name        = "METH1"
    ),
    ASSAY_OSA_LOCAL = list(
      description        = "OSA performed at a local (vs central) laboratory indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (central laboratory OSA or any non-OSA-local sample)",
      notes              = "Abrantes 2017 Table 3 + Table 2 footnote g. Only study B1831003 used a local laboratory for FVIII activity measurement; the -14.6% OSA-local effect on F is interpreted as inter-laboratory variability of the OSA assay (Discussion). Sub-variant of ASSAY_OSA: ASSAY_OSA_LOCAL = 1 implies ASSAY_OSA = 1 (the local-lab samples are a subset of OSA-assayed samples). Per-observation indicator. Source column METH2.",
      source_name        = "METH2"
    )
  )

  population <- list(
    species               = "human",
    n_subjects            = 754L,
    n_studies             = 13L,
    age_range             = "0.0027-73 years (1 day to 73 years)",
    age_median            = "23 years",
    weight_range          = "3.0-134 kg",
    weight_median         = "69 kg",
    sex_female_pct        = 0.1,
    race_ethnicity        = c(White = 87.0, Asian = 7.7, Black = 1.3, Other = 3.8),
    ethnicity_hispanic_pct = 5.2,
    disease_state         = "Moderate to severe hemophilia A; per-study severity criteria reported in Abrantes 2017 Table 3 (most studies enrolled severe FVIII activity < 1 IU/dL or severe + moderately severe FVIII activity < 2 IU/dL). The pooled cohort contains 91 previously untreated patients (PUPs; study B1831054).",
    inhibitor_status      = "Per-observation FVIII inhibitor status (Bethesda assay > 0.6 BU/mL = positive): 174 of 7363 observations (2.4%) positive, 6138 (83%) negative, 1051 (14%) missing; 92.5% of inhibitor-positive samples were low-titer (< 5 BU/mL). Patients with a history of inhibitors were not usually eligible for enrollment, so the inhibitor prevalence in the analysis cohort is lower than the population-level estimate (Discussion).",
    pediatric_cohorts     = "234 of 754 subjects (31%) < 17 years at study entry. Age-cohort weights (Table 1): 0 to < 1 yr n = 62, weight median 8 kg (range 3-12); 1 to < 2 yr n = 21, 11 kg (9-14); 2 to < 6 yr n = 8, 17 kg (11-20); 6 to < 12 yr n = 25, 30 kg (21-57); 12 to < 17 yr n = 118, 56 kg (34-109). Two subjects counted in both the 6-<12 yr and 12-<17 yr cohorts (different studies).",
    dose_range            = "Intravenous moroctocog alfa (ReFacto, ReFacto AF, or Xyntha): median dose 2000 IU (range 100-40000 IU), corresponding to median 32 IU/kg (range 2-1300 IU/kg). Per-study dosing protocols ranged from on-demand to prophylactic to perioperative (Abrantes 2017 Table 3 and Supplement Text S1).",
    regions               = "25 countries (multinational); studies conducted between 1993 and 2013.",
    assay_distribution    = "Per-study analytical assay: 9 studies used CSA (chromogenic substrate assay) at a central laboratory, 3 studies used OSA (one-stage clotting assay) at a central laboratory, and 1 study (B1831003) used OSA at a local laboratory (Abrantes 2017 Table 3).",
    n_observations        = "7363 FVIII activity observations; 910 (12.4%) below the lower limit of quantification (LLOQ: 1 IU/dL in 11 studies, 2 IU/dL in 2 studies); 312 (4.2%) sampling times had 2 or more replicate measurements.",
    notes                 = "Pooled analysis of 13 Pfizer/Wyeth-sponsored clinical trials (Abrantes 2017 Table 3). Single female enrolled (homozygous hemophilia A; B1831090). The full population PK model includes (a) a two-class mixture on endogenous FVIII activity (severe subpopulation 0.474 IU/dL at P = 0.803; moderately severe 1.59 IU/dL at P = 0.197 in 12 studies, or P = 0.110 within studies B1831015 and B1831053), (b) an additional residual-FVIII-activity component capturing incomplete-washout pre-dose activity for subjects whose first-dose pre-dose FVIII exceeded the per-study severity cutoff, and (c) inter-occasion variability of 34.7% on CL and 41.0% on V2; none of (a)-(c) are encoded in the model-library file (see vignette deviations)."
  )

  ini({
    # Structural PK at the reference covariate values: WT = 70 kg, AGE = 20 yr,
    # ADA_POS = 0 (inhibitor negative), RACE_BLACK = 0 (non-Black race),
    # STUDY_B1831090 = 0 (the other 12 studies), FORM_XYNTHA = 0 (ReFacto or
    # ReFacto AF), ASSAY_OSA = 0 (CSA reference), ASSAY_OSA_LOCAL = 0 (central
    # laboratory). Volumes encoded in dL (10x the L value reported in Table 2)
    # so that central / V1 yields the paper's IU/dL observable directly while
    # clearances stay in dL/h as reported. V_ss = V1 + V2 = 24.5 + 9.23 dL
    # = 3.37 L, matching the paper's stated V_ss = 3.38 L (Discussion).
    lcl <- log(2.76)   ; label("Typical clearance CL (dL/h)")                           # Abrantes 2017 Table 2: CL = 2.76 dL/h
    lvc <- log(24.5)   ; label("Typical central volume V1 (dL; = 2.45 L)")              # Abrantes 2017 Table 2: V1 = 2.45 L (encoded as 24.5 dL)
    lq  <- log(25.1)   ; label("Typical inter-compartmental clearance Q (dL/h)")        # Abrantes 2017 Table 2: Q = 25.1 dL/h
    lvp <- log(9.23)   ; label("Typical peripheral volume V2 (dL; = 0.923 L)")          # Abrantes 2017 Table 2: V2 = 0.923 L (encoded as 9.23 dL)

    # Bioavailability anchor (F = 1 for ReFacto/ReFacto AF + CSA central). The
    # covariate-driven F multiplier overrides this anchor for Xyntha-product
    # and OSA-assayed observations (see f_typical in model() below).
    lfdepot <- fixed(log(1)) ; label("Bioavailability anchor (F = 1 for ReFacto + CSA central reference)")  # Abrantes 2017 Methods: F set to 1 for ReFacto and 1.38 for Xyntha (relative to ReFacto)

    # Endogenous baseline FVIII activity. The paper's full model is a two-class
    # mixture on baseline severity (severe 0.474 IU/dL at P = 0.803;
    # moderately severe 1.59 IU/dL at P = 0.197). The severe-subpopulation
    # typical value is encoded here as the default baseline because severe
    # hemophilia A is the majority indication and primary prophylactic-
    # simulation use case. Setting rbase to the moderately severe value (1.59
    # IU/dL) or sampling from the mixture is documented in the vignette.
    lrbase <- log(0.474) ; label("Typical baseline endogenous FVIII activity for the severe hemophilia A subpopulation (IU/dL)")  # Abrantes 2017 Table 2: Endogenous FVIII activity 1 = 0.474 IU/dL (P = 0.803 subpopulation)

    # Allometric exponents. Theory-based 0.75 on CL and Q are held fixed; the
    # shared exponent 0.812 on V1 and V2 was estimated jointly.
    e_wt_cl    <- fixed(0.75) ; label("Allometric exponent of (WT/70) on CL (theory-based; fixed)")   # Abrantes 2017 Table 2 footnote b: "(WT / 70)^3/4"
    e_wt_q     <- fixed(0.75) ; label("Allometric exponent of (WT/70) on Q (theory-based; fixed)")    # Abrantes 2017 Table 2 footnote d: "(WT / 70)^3/4"
    e_wt_vc_vp <- 0.812       ; label("Shared estimated allometric exponent of (WT/70) on V1 and V2 (unitless)")  # Abrantes 2017 Table 2: Allometric exponent for V1 and V2 = 0.812

    # Piecewise age-on-CL slopes (Abrantes 2017 Table 2 footnote b). Two
    # pieces meet at AGE = 1 yr; AGE reference (where age_eff_cl = 1) is
    # AGE = 20 yr.
    e_age_below1_cl <- 0.149    ; label("Slope of (AGE - 1) on age-on-CL factor below 1 year of age (1/year)")  # Abrantes 2017 Table 2: Age on CL up to 1 year old = 0.149
    e_age_above1_cl <- -0.00678 ; label("Slope of (AGE - 20) on age-on-CL factor above 1 year of age (1/year)") # Abrantes 2017 Table 2: Age on CL above 1 year old = -0.00678

    # Inhibitor + study effects on CL (Abrantes 2017 Table 2 footnote b).
    e_ada_pos_cl        <- 1.66   ; label("Fractional change in CL when FVIII-inhibitor positive (unitless)")     # Abrantes 2017 Table 2: Inhibitor status on CL = +166% if positive
    e_study_b1831090_cl <- -0.347 ; label("Fractional change in CL for study B1831090 cohort (unitless)")         # Abrantes 2017 Table 2: Study effect on CL = -34.7% if B1831090

    # Race-on-V2 effect (Abrantes 2017 Table 2 footnote e).
    e_race_black_vp <- 0.884 ; label("Fractional change in V2 in Black-race subjects (unitless)")  # Abrantes 2017 Table 2: Race on V2 = +88.4% if Black

    # F covariate effects (Abrantes 2017 Table 2 footnote g).
    e_form_xyntha_f     <- 0.380  ; label("Fractional change in F for Xyntha (vs ReFacto/ReFacto AF) (unitless)")          # Abrantes 2017 Table 2 footnote g: F = 1 * 1.38^PROD => +38.0% if Xyntha
    e_assay_osa_f       <- -0.390 ; label("Fractional change in F for OSA central-laboratory measurement (unitless)")      # Abrantes 2017 Table 2: Activity bias = -39.0% if OSA
    e_assay_osa_local_f <- -0.146 ; label("Fractional change in F for OSA local-laboratory measurement (unitless)")        # Abrantes 2017 Table 2: Activity bias = -14.6% if OSA local

    # Residual-error magnitude effect (Abrantes 2017 Table 2).
    e_assay_osa_propsd <- 0.403 ; label("Fractional increase in proportional residual SD when OSA-assayed (unitless)")  # Abrantes 2017 Table 2: Assay on residual error = +40.3% CV% if OSA

    # IIV (log-normal). omega^2 = log(1 + CV^2) for each estimated component
    # (Abrantes 2017 Table 2 IIV %CV values).
    #   CL    IIV 30.5%  CV  -> log(1 + 0.305^2)  = 0.08894
    #   F     IIV 13.0%  CV  -> log(1 + 0.130^2)  = 0.01676
    #   rbase IIV 7.19%  CV  -> log(1 + 0.0719^2) = 0.005154
    # Inter-occasion variability of 34.7% on CL and 41.0% on V2 (Table 2) and
    # IIV on the residual-activity component are NOT encoded structurally
    # here; see vignette deviations.
    etalcl     ~ 0.08894    # Abrantes 2017 Table 2: IIV CL = 30.5% CV
    etalfdepot ~ 0.01676    # Abrantes 2017 Table 2: IIV F = 13.0% CV
    etalrbase  ~ 0.005154   # Abrantes 2017 Table 2: IIV endogenous FVIII activity = 7.19% CV

    # Proportional residual error (CSA reference value; OSA magnitude switched
    # in the model() body via the per-observation ASSAY_OSA indicator).
    propSd <- 0.192 ; label("Proportional residual SD for CSA-assayed observations (fraction)")  # Abrantes 2017 Table 2: Proportional residual error = 19.2% CV
  })

  model({
    # Allometric weight scaling at WT = 70 kg reference.
    wt_cl    <- (WT / 70) ^ e_wt_cl
    wt_q     <- (WT / 70) ^ e_wt_q
    wt_vc_vp <- (WT / 70) ^ e_wt_vc_vp

    # Piecewise age effect on CL (Abrantes 2017 Table 2 footnote b). Pieces
    # meet at AGE = 1 yr; the AGE = 20-year reference defines age_eff_cl = 1.
    #   AGE <= 1: age_eff = 1 + 0.149 * (AGE - 1) - (20 - 1) * (-0.00678)
    #                     = 1.1288 + 0.149 * (AGE - 1)
    #   AGE >  1: age_eff = 1 + (-0.00678) * (AGE - 20)
    age_eff_below1 <- 1 + e_age_below1_cl * (AGE - 1) - (20 - 1) * e_age_above1_cl
    age_eff_above1 <- 1 + e_age_above1_cl * (AGE - 20)
    age_eff_cl <- (AGE < 1) * age_eff_below1 + (AGE >= 1) * age_eff_above1

    # Categorical / binary covariate multipliers.
    inh_eff_cl   <- 1 + e_ada_pos_cl        * ADA_POS
    study_eff_cl <- 1 + e_study_b1831090_cl * STUDY_B1831090
    race_eff_vp  <- 1 + e_race_black_vp     * RACE_BLACK

    # F multiplier (covariate-driven). Abrantes 2017 Table 2 footnote g:
    #   F = 1 * 1.38^PROD * (1 - 0.390 * METH1) * (1 - 0.146 * METH2)
    # Since FORM_XYNTHA / ASSAY_OSA / ASSAY_OSA_LOCAL are binary, the
    # 1.38^FORM_XYNTHA factor is equivalent to (1 + 0.38 * FORM_XYNTHA).
    f_typical <- (1 + e_form_xyntha_f     * FORM_XYNTHA)    *
                 (1 + e_assay_osa_f       * ASSAY_OSA)      *
                 (1 + e_assay_osa_local_f * ASSAY_OSA_LOCAL)

    # Individual PK parameters with the Abrantes 2017 covariate equations.
    cl     <- exp(lcl + etalcl) * wt_cl * age_eff_cl * inh_eff_cl * study_eff_cl
    vc     <- exp(lvc)          * wt_vc_vp
    q      <- exp(lq)           * wt_q
    vp     <- exp(lvp)          * wt_vc_vp * race_eff_vp
    fdepot <- exp(lfdepot + etalfdepot) * f_typical
    rbase  <- exp(lrbase + etalrbase)

    # Micro-constants for explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV disposition. Dose enters the central compartment
    # directly; the bioavailability multiplier is applied to the central
    # compartment via f(central) <- fdepot to mirror the NONMEM F treatment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1
    f(central) <- fdepot

    # Observation: FVIII activity (IU/dL) is the sum of the exogenous-drug
    # component and the endogenous baseline. With dose in IU and V1 in dL,
    # central / vc yields IU/dL directly.
    Cc <- central / vc + rbase

    # Per-observation proportional residual SD: OSA-assayed observations have
    # +40.3% magnitude relative to the CSA reference. Switched via the
    # ASSAY_OSA indicator (Andrews 2017 tacrolimus precedent).
    propSd_eff <- propSd * (1 + e_assay_osa_propsd * ASSAY_OSA)
    Cc ~ prop(propSd_eff)
  })
}

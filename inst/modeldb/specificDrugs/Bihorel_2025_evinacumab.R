Bihorel_2025_evinacumab <- function() {
  description <- "Population PK/PD model for evinacumab in children (5 to <12 years), adolescents, and adults with homozygous familial hypercholesterolemia and in phase 1 participants (Bihorel 2025): two-compartment PK with first-order SC absorption (lag time, partial bioavailability) and parallel linear plus Michaelis-Menten elimination, with linear disposition parameters allometrically scaled by time-varying body weight, linked to a type 1 indirect-response model in which evinacumab inhibits LDL-C production and a second, time-varying LDL-C elimination process quantifies the lipoprotein-apheresis effect."
  reference   <- "Bihorel S, Dingman R, Mendell J, Wang Y, Banerjee P, Pordy R, Davis JD, DiCioccio AT, Harnisch L. Population pharmacokinetics and exposure-response modeling for evinacumab in children, adolescents, and adults with homozygous familial hypercholesterolemia. CPT Pharmacometrics Syst Pharmacol. 2025;14(11):1823-1834. doi:10.1002/psp4.70016. Refines the earlier model of Pu X et al. CPT Pharmacometrics Syst Pharmacol. 2021;10(11):1412-1421 (doi:10.1002/psp4.12711); see modellib('Pu_2021_evinacumab')."
  vignette    <- "Bihorel_2025_evinacumab"

  units       <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TIME-VARYING, not baseline. Bihorel 2025 Methods section 2.5 states that a refinement over the predecessor model was defining allometric scaling 'as a function of time-varying body weight as opposed to a fixed baseline body weight'; the supplement control streams confirm this by centring on WGT (CWEIGHT = WGT/72), not WGTBL. Enters PK as power scaling on CL, Vc, Q and Vp with reference 72 kg (Table 2 typical-value equations), and enters PD as an ADDITIVE linear term on the logit of Imax centred by subtraction at the same 72 kg (supplement Code 2 CLWEIGHT = WGT-72), so larger patients have a smaller maximum inhibitory effect. The 72 kg reference is the rounded model-centring constant and is distinct from the 74.1 kg median baseline weight used as the reference patient in the Figure 2 tornado plots.",
      source_name        = "WGT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (time-fixed). Power effect on the estimated baseline LDL-C concentration with reference 43 years (Table 3 typical-value equation, supplement Code 2 CAGE = AGE/43); the negative exponent means younger patients have a higher baseline LDL-C. Age was newly investigated in this refinement and was retained on baseline LDL-C only; it was not retained on any PK parameter.",
      source_name        = "AGE"
    ),
    ANGPTL3 = list(
      description        = "Baseline total serum angiopoietin-like protein 3 concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline only (per-subject time-fixed). Power-form effect on the maximum saturable elimination rate Vmax with reference 0.0908 mg/L (supplement Code 1 CANGBL = ANGBL/0.0908, printed in the Table 2 typical-value equation). Note the reference differs from the 0.08 mg/L used by the predecessor Pu 2021 model. The assay detects both free ANGPTL3 and ANGPTL3 bound to evinacumab (Methods section 2.2; LLOQ 0.0195 mg/L in neat human serum). Higher baseline target predicts a faster saturable elimination.",
      source_name        = "ANGBL"
    ),
    DIS_HOFH = list(
      description        = "Homozygous familial hypercholesterolemia disease indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 1 participant reference)",
      notes              = "Time-fixed per subject. The paper states explicitly that 'HoFH is 0 for phase 1 adult participants and 1 for patients with HoFH' (Table 2 footnote), so the reference cohort here is the pooled phase 1 population rather than a healthy-volunteer cohort per se. Multiplicative effect on Vmax: the printed equation is Vmax = 3.03 * (ANGPTL3/0.0908)^0.395 * 0.75^HoFH, i.e. HoFH patients have a 25% lower target-mediated Vmax. Encoded here in the log domain as exp(e_dis_hofh_vmax * DIS_HOFH), which reproduces supplement Code 1 verbatim: EXP(DISTYPN*LOG(THETA(12))).",
      source_name        = "DISTYPN"
    ),
    APHERESIS_LIPO_ACTIVE = list(
      description        = "Lipoprotein-apheresis-active indicator (time-varying per-session gate)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no apheresis session running)",
      notes              = "WITHIN-SUBJECT TIME-VARYING gate, 1 only while a lipoprotein-apheresis (LA) session is physically running and 0 otherwise (including in the 81.7% of the analysis population never treated with LA). Switches on the second, LA-driven first-order elimination arm of the LDL-C turnover pool; reproduces supplement Code 2 $DES term -KAPH*APHON*A(4). Quantifying this effect separately from the drug effect is the central novelty of this refinement. Session duration is carried by the event table, not by this column: the supplement's imputation rule assumes a 2.5-h session where start or stop time was missing, and the paper reports that the estimated rate of 8.36 1/day gives 'approximately a 60% LDL-C reduction for a typical 2.5-h LA session' (1 - exp(-8.36 * 2.5/24) = 0.581, matching). LA frequency in the analysis population was weekly (8.07%), bi-weekly (9.63%) or monthly (0.62%) per Table 1.",
      source_name        = "APHON"
    )
  )

  covariatesDataExcluded <- list(
    RACE_WHITE = list(
      description = "White (Caucasian) race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Retained on Imax in the predecessor Pu 2021 model but formally dropped here: Results section 3.3 states 'the effect of racial classification on Imax included in the original PK/PD model was not found to be statistically significant in the refined model and was, thus, removed' in a limited backward-elimination step. The supplement Code 2 $PK block still derives RAC1 from RACEN but never uses it in any parameter expression. Documented rather than implemented so the provenance of the removal is preserved."
    ),
    LDLC = list(
      description = "Baseline serum low-density lipoprotein cholesterol concentration",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Not a covariate in this model. The predecessor Pu 2021 model read the observed baseline LDL-C from the data both to initialise the LDL-C state and to scale IC50; this refinement instead ESTIMATES the baseline as a model parameter (lrbase, Table 3 LDLC0 = 214 mg/dL) with IIV and an age effect, and IC50 carries no covariate. Supplement Code 2 uses LDLBL for a single subject only (ID 244) as a data-driven override, which is a data-handling exception rather than a covariate relationship and is not reproduced here."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "evinacumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "evinacumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "evinacumab", units = "mg", specimen = "serum", verified = TRUE),
    # The LDL-C turnover state holds a CONCENTRATION, not an amount: supplement Code 2
    # defines no scaling factor for compartment 4 and kin is reported in mg/dL/day
    # (Table 3), so A(4) is directly the serum LDL-C concentration in mg/dL.
    ldl         = list(analyte = "low-density lipoprotein cholesterol", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 322,
    n_studies      = 7,
    age_range      = "5-75 years",
    age_median     = "mean 40.0 years (SD 15.0)",
    weight_range   = "19.7-152 kg",
    weight_median  = "mean 73.7 kg (SD 20.1); median baseline 74.1 kg",
    sex_female_pct = 46.3,
    race_ethnicity = c(White = 73.0, `Black/African American` = 4.04, Asian = 16.1,
                       `American Indian/Alaska Native` = 0.621, Other = 2.8, Unknown = 3.42),
    disease_state  = "homozygous familial hypercholesterolemia (n = 139) pooled with phase 1 participants having elevated triglycerides and/or LDL-C but otherwise healthy (n = 183)",
    dose_range     = "5-20 mg/kg IV single or repeated; 75-450 mg SC single or repeated. The regimen carried through the paper's simulations is 15 mg/kg IV every 4 weeks.",
    notes          = "Baseline demographics from Table 1; study designs from Table S1. Age strata: 20 children (5 to <12 years), 14 adolescents (12 to <18 years), 288 adults (>=18 years). The PK analysis used 5085 quantifiable evinacumab concentrations (plus 613 below the 0.078 mg/L LLOQ, handled by the M3 likelihood method in supplement Code 1). The PK/PD analysis is a SUBSET: 3316 LDL-C measurements in the 139 patients with HoFH aged 5-75 years only, since phase 1 participants were excluded from the PD dataset. Lipoprotein apheresis at baseline: 81.7% none, 8.07% weekly, 9.63% bi-weekly, 0.62% monthly. Anti-drug antibodies were positive in only 2 of 322 participants (0.621%) and were not tested as a covariate."
  )

  ini({
    # ---- Structural PK: reference weight 72 kg, reference baseline ANGPTL3 0.0908 mg/L ----
    # Values are the full-precision final estimates that supplement Code 2 hardcodes from the
    # PK run (run234) as PKTH1-PKTH15; each rounds to the value printed in Table 2.
    lcl    <- log(0.0993107); label("Linear elimination clearance, CL (L/day)")                     # Suppl Code 2 PKTH1 (Table 2: 0.0993)
    lvc    <- log(2.77705);   label("Central volume, Vc (L)")                                       # Suppl Code 2 PKTH2 (Table 2: 2.78)
    lq     <- log(0.199605);  label("Distribution clearance, Q (L/day)")                            # Suppl Code 2 PKTH3 (Table 2: 0.200)
    lvp    <- log(2.06583);   label("Peripheral volume, Vp (L)")                                    # Suppl Code 2 PKTH4 (Table 2: 2.07)
    lvmax  <- log(3.02687);   label("Maximum saturable elimination rate, Vmax (mg/day/L)")          # Suppl Code 2 PKTH5 (Table 2: 3.03)
    lkm    <- log(2.61118);   label("Michaelis-Menten constant, km (mg/L)")                         # Suppl Code 2 PKTH6 (Table 2: 2.61)
    lka    <- log(0.187704);  label("First-order SC absorption rate constant, ka (1/day)")          # Suppl Code 2 PKTH7 (Table 2: 0.188)
    ltlag  <- log(0.136084);  label("SC absorption lag time, ALAG1 (day)")                          # Suppl Code 2 PKTH8 (Table 2: 0.136)
    lfdepot <- log(0.709165); label("Bioavailability after SC dosing, F1 (fraction)")               # Suppl Code 2 PKTH9 (Table 2: 0.709)

    # ---- PK covariate effects ----
    e_wt_cl  <- 0.775851;  label("Allometric exponent on CL for WT/72 (unitless)")                  # Suppl Code 2 PKTH13 (Table 2: 0.776)
    e_wt_vc  <- 0.668333;  label("Allometric exponent on Vc for WT/72 (unitless)")                  # Suppl Code 2 PKTH10 (Table 2: 0.668)
    e_wt_q   <- 1.08198;   label("Allometric exponent on Q for WT/72 (unitless)")                   # Suppl Code 2 PKTH14 (Table 2: 1.08)
    e_wt_vp  <- 0.986423;  label("Allometric exponent on Vp for WT/72 (unitless)")                  # Suppl Code 2 PKTH15 (Table 2: 0.986)
    e_angptl3_vmax   <- 0.394701;  label("Power exponent on Vmax for ANGPTL3/0.0908 (unitless)")    # Suppl Code 2 PKTH11 (Table 2: 0.395)
    # Table 2 reports the HoFH effect as the PROPORTIONAL multiplier 0.750 (a 25% lower Vmax).
    # Supplement Code 1/2 apply it as EXP(DISTYPN*LOG(THETA(12))), so the log is taken here.
    e_dis_hofh_vmax  <- log(0.749573); label("Log proportional effect of HoFH on Vmax (unitless)")  # Suppl Code 2 PKTH12 (Table 2: 0.750)

    # ---- Structural PD (patients with HoFH only) ----
    lrbase <- log(214);   label("Baseline LDL-C concentration, LDLC0 (mg/dL)")                      # Table 3
    lkin   <- log(34.7);  label("LDL-C production rate, kin (mg/dL/day)")                           # Table 3
    lic50  <- log(32.7);  label("Evinacumab concentration at half-maximal inhibition, IC50 (mg/L)") # Table 3
    # Imax is parameterised on the LOGIT scale (Methods 2.6: 're-parameterisation of the maximum
    # inhibitory effect of evinacumab (Imax) on the logit scale'). Table 3 prints the natural-scale
    # typical value 0.574; supplement Code 2 forms LIMAX = LOG(THETA(3)/(1-THETA(3))) + ... .
    logitimax <- log(0.574 / (1 - 0.574)); label("Logit of maximum inhibition of LDL-C production, Imax (unitless)")  # Table 3 (Imax = 0.574)
    # Second elimination arm of the LDL-C pool, active only while an apheresis session runs.
    lkout_apheresis <- log(8.36); label("First-order LDL-C elimination rate constant via lipoprotein apheresis (1/day)")  # Table 3

    # ---- PD covariate effects ----
    e_age_rbase <- -0.320;   label("Power exponent on baseline LDL-C for AGE/43 (unitless)")        # Table 3
    # ADDITIVE on the logit scale, per kg, centred by SUBTRACTION at 72 kg (Code 2: CLWEIGHT = WGT-72).
    e_wt_imax   <- -0.0186;  label("Linear effect of (WT - 72) on logit(Imax) (1/kg)")              # Table 3

    # ---- IIV ----
    # Tables 2 and 3 report IIV as "% CV". The back-transform is omega^2 = log(1 + CV^2):
    # the alternative reading omega = CV is excluded by the Table 5 stochastic simulation, where
    # the simulated baseline LDL-C SD/mean is 0.555-0.577 across age strata. Because age varies
    # within each stratum and can only ADD spread, the pure-eta CV must lie at or below those
    # values; omega = CV would force it to 0.578, above two of the three strata, whereas
    # omega^2 = log(1 + CV^2) gives 0.537, consistent with all three.
    etalcl + etalvc ~ c(0.072391,
                        0.032, 0.131576)  # Table 2: IIV in CL 27.4% CV, cov(CL,Vc) 0.032 (reported directly on the omega scale), IIV in Vc 37.5% CV
    etalka   ~ 0.553655   # Table 2: IIV in ka 86.0% CV
    etaltlag ~ 1.401799   # Table 2: IIV in ALAG1 175% CV
    etalrbase ~ 0.253377  # Table 3: IIV in LDLC0 53.7% CV
    etalkin   ~ 0.260073  # Table 3: IIV in kin 54.5% CV
    # IIV on Imax is additive on the LOGIT scale. Supplement Code 2 tags this eta with its own
    # reporting formula, ;--eta2- IIV in IMAX [cv=100*(1-th3)*eta2], i.e. the delta-method
    # approximation CV = (1 - Imax) * omega. Inverting it: omega = 0.300 / (1 - 0.574) = 0.704225.
    etalogitimax ~ 0.495933  # Table 3: IIV in Imax 30.0% CV, via the Code 2 [cv=...] tag

    # ---- Residual error ----
    # NONMEM $SIGMA holds VARIANCES, and supplement Code 1/2 confirm it by forming
    # W = SQRT(IPRED**2*SIGMA(1,1) + SIGMA(2,2)); the tabulated values are therefore squared SDs.
    propSd     <- 0.2582634; label("Proportional residual error for evinacumab (fraction)")        # Table 2: 0.0667 variance
    addSd      <- 0.2839014; label("Additive residual error for evinacumab (mg/L)")                # Table 2: 0.0806 variance
    propSd_ldl <- 0.2511971; label("Proportional residual error for LDL-C (fraction)")             # Table 3: 0.0631 variance
    addSd_ldl  <- 7.7653075; label("Additive residual error for LDL-C (mg/dL)")                    # Table 3: 60.3 variance
  })

  model({
    # ---- Individual PK parameters (reference weight 72 kg; WT is time-varying) ----
    cl   <- exp(lcl + etalcl) * (WT / 72)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT / 72)^e_wt_vc
    q    <- exp(lq)           * (WT / 72)^e_wt_q
    vp   <- exp(lvp)          * (WT / 72)^e_wt_vp
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    km   <- exp(lkm)
    # Vmax carries no IIV in the source model.
    vmax <- exp(lvmax + e_dis_hofh_vmax * DIS_HOFH) * (ANGPTL3 / 0.0908)^e_angptl3_vmax

    Cc <- central / vc

    # ---- PK ODEs (Figure 1; supplement Code 1 $DES) ----
    # The Michaelis-Menten term is Vmax * A(2) / (km + Cc): an AMOUNT in the numerator against a
    # CONCENTRATION in the denominator, which is why Table 2 gives Vmax the units mg/day/L rather
    # than mg/day. It is equivalent to a conventional Vmax * Cc / (km + Cc) whose maximum rate is
    # Vmax * Vc, so the saturable capacity scales with the individual central volume.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          (q  / vc) * central +
                          (q  / vp) * peripheral1 -
                          vmax * central / (km + Cc)
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # SC doses enter the depot and carry the lag time and partial bioavailability;
    # IV doses bypass the depot by dosing central directly.
    f(depot)    <- exp(lfdepot)
    alag(depot) <- tlag

    # ---- Individual PD parameters (patients with HoFH; reference age 43 y, weight 72 kg) ----
    rbase <- exp(lrbase + etalrbase) * (AGE / 43)^e_age_rbase
    kin   <- exp(lkin + etalkin)
    ic50  <- exp(lic50)
    # Kept on two lines, mirroring supplement Code 2 (LIMAX, then the inverse-logit): the
    # individual logit must be a simple additive expression for mu-referencing to be recognised.
    limax <- logitimax + e_wt_imax * (WT - 72) + etalogitimax
    imax  <- expit(limax)
    kout_apheresis <- exp(lkout_apheresis)
    # kout is derived so that the LDL-C pool sits at its individual baseline before treatment
    # (supplement Code 2: KOUT = KIN/LDL0), using the INDIVIDUAL kin and baseline.
    kout  <- kin / rbase

    # ---- LDL-C indirect response (type 1: inhibition of production) ----
    # The apheresis arm is a second first-order loss that is gated on and off by the
    # time-varying APHERESIS_LIPO_ACTIVE indicator (supplement Code 2 $DES: -KAPH*APHON*A(4)).
    ldl(0)    <- rbase
    d/dt(ldl) <- kin * (1 - imax * Cc / (ic50 + Cc)) -
                 kout * ldl -
                 kout_apheresis * APHERESIS_LIPO_ACTIVE * ldl

    # ---- Observation and error models ----
    # Cc in mg/L (= ug/mL); ldl in mg/dL.
    Cc  ~ add(addSd)     + prop(propSd)
    ldl ~ add(addSd_ldl) + prop(propSd_ldl)
  })
}

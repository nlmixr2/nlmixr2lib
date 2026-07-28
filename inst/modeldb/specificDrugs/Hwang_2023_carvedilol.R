Hwang_2023_carvedilol <- function() {
  description <- "Sequential population PK-PD model for oral carvedilol in 21 healthy Korean male volunteers genotyped for CYP2D6 (Hwang 2023). PK is a two-compartment model (central + peripheral1) with zero-order absorption of duration D1 into the central compartment preceded by an absorption lag time Tlag, and first-order elimination; volumes and clearances are apparent (CL/F, Vc/F, Vp/F, Q/F). The only retained covariate is the CYP2D6*10/*10 (intermediate-metabolizer-2) genotype, which reduces CL/F by 32.8 percent relative to the pooled CYP2D6*1/*1, *1/*2, *1/*10, *2/*10 reference. PD is a direct-effect competitive-antagonism Emax model for the heart-rate response to an isoproterenol sensitivity test: HR = E0 + Emax * D / (ED50 * (1 + C / IC50) + D), where D is the isoproterenol challenge dose in ug (covariate DOSE_ISOPROTERENOL_UG), C is the carvedilol plasma concentration in ng/mL, and the Hill exponent was fixed at 1. Carvedilol competitively shifts the isoproterenol dose-response curve to the right without changing Emax; no covariate, including CYP2D6 phenotype, was retained on any PD parameter."
  reference <- paste(
    "Hwang S, Lee S, Yoon J, Chung JY.",
    "Population Pharmacokinetic-Pharmacodynamic Modeling of Carvedilol to Evaluate the Effect of Cytochrome P450 2D6 Genotype on the Heart Rate Reduction.",
    "J Korean Med Sci. 2023 Jun 5;38(22):e173.",
    "doi:10.3346/jkms.2023.38.e173",
    sep = " "
  )
  vignette <- "Hwang_2023_carvedilol"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    CYP2D6_STAR10_HOM = list(
      description        = "CYP2D6*10 (rs1065852) homozygous-mutant indicator; 1 = CYP2D6*10/*10 (the paper's intermediate-metabolizer-2, IM-2, group), 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled CYP2D6*1/*1, *1/*2 extensive metabolizers and CYP2D6*1/*10, *2/*10 intermediate-metabolizer-1 subjects)",
      notes              = paste(
        "Time-fixed germline genotype. Hwang 2023 stratified the 21 subjects into three",
        "phenotype groups -- EM (*1/*1, *1/*2; n = 6), IM-1 (*1/*10, *2/*10; n = 7) and IM-2",
        "(*10/*10; n = 8) -- but the final covariate model distinguished only IM-2 from the",
        "pooled EM + IM-1 reference, so a single binary indicator is sufficient. Table 2",
        "footnote a gives the retained relationship explicitly as CL/F = 153 * (1 - 0.328 * GT)",
        "with GT = 1 for IM-2 and GT = 0 for EM and IM-1. Adding this covariate reduced the",
        "CL/F inter-individual variability from 26.8 percent to 14.6 percent CV",
        "(delta OFV = -27.50; Hwang 2023 Results 'PK modeling'). Because the paper pooled the",
        "heterozygous *10 carriers (IM-1) with the wild-type EM group, the paired canonical",
        "CYP2D6_STAR10_HET indicator carries no effect in this model and is documented under",
        "covariatesDataExcluded instead.",
        sep = " "
      ),
      source_name        = "GT (Hwang 2023 Table 2 footnote a); CYP2D6 phenotype group IM-2 in Table 1"
    ),
    DOSE_ISOPROTERENOL_UG = list(
      description        = "Per-observation intravenous isoproterenol challenge dose administered during the isoproterenol sensitivity test (IST), ug",
      units              = "ug",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-observation challenge-dose column that drives the algebraic direct-effect Emax",
        "term; it is NOT an rxode2 dosing event (the PD model is algebraic in D, with no",
        "isoproterenol PK compartment and no isoproterenol elimination rate constant reported).",
        "Set DOSE_ISOPROTERENOL_UG = 0 on rows with no isoproterenol challenge, which recovers",
        "HR = E0. Hwang 2023 Methods 'Data': the IST escalation was 0.25, 0.5, 1 and 2 ug at",
        "baseline and 5, 10, 20 and 40 ug after single and multiple carvedilol doses;",
        "isoproterenol was injected until heart rate exceeded 140 bpm or rose 30 bpm above the",
        "pre-dose rate. Follows the DOSE_AGT_UG precedent established by",
        "BuchwalderCsajka_1999_angiotensin.R.",
        sep = " "
      ),
      source_name        = "D (Hwang 2023 equations 3-4 and Fig. 1; 'isoproterenol dose')"
    )
  )

  covariatesDataExcluded <- list(
    CYP2D6_STAR10_HET = list(
      description        = "CYP2D6*10 (rs1065852) heterozygote indicator; 1 = CYP2D6*1/*10 or *2/*10 (the paper's IM-1 group)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2D6*1/*1, *1/*2)",
      notes              = paste(
        "Hwang 2023 genotyped and analysed the heterozygous *10 group (IM-1, n = 7) as a",
        "distinct phenotype stratum, but the final PK covariate model pooled IM-1 with the EM",
        "group as the reference: only the *10/*10 homozygotes had a detectably lower CL/F",
        "(Results 'PK modeling'; Table 2 footnote a assigns GT = 0 to both EM and IM-1).",
        "Documented here to preserve the three-level genotype screen without carrying an",
        "unused covariateData entry.",
        sep = " "
      ),
      source_name        = "CYP2D6 phenotype group IM-1 (Hwang 2023 Table 1)"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on both the PK and the PK-PD parameters via the median-normalised power model of Hwang 2023 equation (6); not significant and not retained (Results 'PK modeling' and 'PK-PD modeling'). Cohort mean (SD) 27.3 (4.7) years (Table 1).",
      source_name        = "Age (Hwang 2023 Table 1)"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on both the PK and the PK-PD parameters; not significant and not retained (Results 'PK modeling' and 'PK-PD modeling'). Cohort mean (SD) 174.2 (6.6) cm (Table 1).",
      source_name        = "Height (Hwang 2023 Table 1)"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on both the PK and the PK-PD parameters; not significant and not retained (Results 'PK modeling' and 'PK-PD modeling'). Cohort mean (SD) 70.4 (9.6) kg (Table 1). Hwang 2023 Discussion notes that other carvedilol popPK analyses in congestive-heart-failure patients did retain total body weight, and attributes the difference to the narrow demographic range of this healthy-volunteer cohort.",
      source_name        = "Weight (Hwang 2023 Table 1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    n_observations = "450 carvedilol plasma concentrations and 1003 heart-rate observations (Hwang 2023 Results 'Study population')",
    age_range      = "mean (SD) 27.3 (4.7) years; group means 24.4-29.6 years",
    weight_range   = "mean (SD) 70.4 (9.6) kg; group means 64.5-79.3 kg",
    height_range   = "mean (SD) 174.2 (6.6) cm",
    bmi_range      = "mean (SD) 23.1 (2.4) kg/m^2",
    sex_female_pct = 0,
    race_ethnicity = "Korean",
    disease_state  = "Healthy adult male volunteers",
    dose_range     = "Carvedilol 12.5 mg oral single dose followed by 25 mg oral once daily multiple dosing; isoproterenol sensitivity test boluses of 0.25, 0.5, 1 and 2 ug at baseline and 5, 10, 20 and 40 ug after single and multiple carvedilol doses",
    regions        = "Republic of Korea (Seoul National University Bundang Hospital, single centre)",
    genotype       = "CYP2D6 extensive metabolizer (*1/*1, *1/*2) n = 6; intermediate metabolizer-1 (*1/*10, *2/*10) n = 7; intermediate metabolizer-2 (*10/*10) n = 8",
    notes          = paste(
      "Open-label, one-sequence, multiple-dosing study; all 21 subjects were male",
      "(Hwang 2023 Table 1). PK sampling at 0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12 and 24 h",
      "after the single 12.5 mg dose and after multiple 25 mg once-daily doses. There were",
      "no significant demographic differences between the three CYP2D6 phenotype groups",
      "(Table 1). Only the parent drug was modelled: metabolic ratios for the three",
      "carvedilol metabolites (4'-hydroxyphenyl carvedilol, 5'-hydroxyphenyl carvedilol and",
      "O-desmethyl carvedilol) were below 0.1 in most subjects (Hwang 2023 Discussion).",
      "Carvedilol was administered as the racemate; the model describes total (racemic)",
      "carvedilol plasma concentrations, and the IC50 of 16.5 ng/mL is therefore a",
      "racemic-total value even though only S(-)-carvedilol carries the beta-blocking",
      "activity (Hwang 2023 Discussion).",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # PK: two-compartment model, zero-order absorption of duration D1 into the
    # central compartment delayed by Tlag, first-order elimination.
    # Hwang 2023 Results 'PK modeling' and Fig. 1; estimates from Table 2.
    # All disposition parameters are apparent (divided by the unknown oral
    # bioavailability F), so no separate bioavailability term is estimated.
    # ========================================================================
    ld1   <- log(0.383); label("Zero-order absorption duration D1 (h)")                          # Hwang 2023 Table 2: D1 = 0.383 h (RSE 70%; bootstrap median 0.390, 95% CI 0.330-0.515)
    ltlag <- log(0.215); label("Absorption lag time Tlag (h)")                                   # Hwang 2023 Table 2: Tlag = 0.215 h (RSE 0.6%; bootstrap median 0.205, 95% CI 0.166-0.246)
    lcl   <- log(153);   label("Apparent clearance CL/F (L/h)")                                  # Hwang 2023 Table 2: CL/F = 153 L/h (RSE 0.7%; bootstrap median 147.6, 95% CI 131-170)
    lvc   <- log(440);   label("Apparent central volume of distribution Vc/F (L)")               # Hwang 2023 Table 2: Vc/F = 440 L (RSE 35%; bootstrap median 414, 95% CI 365-501)
    lvp   <- log(754);   label("Apparent peripheral volume of distribution Vp/F (L)")            # Hwang 2023 Table 2: Vp/F = 754 L (RSE 2.3%; bootstrap median 732, 95% CI 553-992)
    lq    <- log(41.3);  label("Apparent inter-compartmental clearance Q/F (L/h)")               # Hwang 2023 Table 2: Q/F = 41.3 L/h (RSE 5.1%; bootstrap median 41.3, 95% CI 35.6-50.9)

    # ---- Retained covariate effect on CL/F ---------------------------------
    # Hwang 2023 equation (7), P_i = P_TV * (1 - Pdiff) * exp(eta_i), made
    # explicit in Table 2 footnote a as
    #     CL/F = 153 * (1 - 0.328 * GT),  GT = 1 (IM-2) or 0 (EM, IM-1)
    # so the fractional reduction is 32.8% and CL/F(IM-2) = 102.8 L/h.
    e_cyp2d6_star10_hom_cl <- 0.328; label("Fractional reduction in CL/F for CYP2D6*10/*10 (IM-2) subjects (unitless)")  # Hwang 2023 Table 2: IM-2 effect on CL/F = 0.328 (RSE 14.6%; bootstrap median 0.314, 95% CI 0.203-0.417)

    # ========================================================================
    # PD: direct-effect Emax model for the isoproterenol sensitivity test,
    # modified for competitive antagonism by carvedilol. Hwang 2023 Fig. 1
    # prints the final-model form
    #     HR = E0 + Emax * D / (ED50 * (1 + C / IC50) + D)
    # which is the general equation (4) with the Hill exponent gamma fixed at
    # 1 (Results 'PK-PD modeling': "A simple Emax model was selected with
    # gamma fixed as 1"). Estimates from Table 3.
    # ========================================================================
    le0   <- log(60.4);  label("Baseline heart rate E0 (bpm)")                                                    # Hwang 2023 Table 3: E0 = 60.4 bpm (RSE 3.1%; bootstrap median 60.4, 95% CI 56.9-64.0)
    lemax <- log(30.7);  label("Maximal isoproterenol-induced heart-rate increase Emax (bpm)")                    # Hwang 2023 Table 3: Emax = 30.7 bpm (RSE 21.9%; bootstrap median 31.7, 95% CI 20.2-52.1)
    led50 <- log(0.685); label("Isoproterenol dose producing half-maximal effect ED50 (ug)")                      # Hwang 2023 Table 3: ED50 = 0.685 ug (RSE 30.9%; bootstrap median 0.709, 95% CI 0.387-1.512); Results text rounds this to 0.69 ug
    lic50 <- log(16.5);  label("Carvedilol concentration producing half-maximal inhibition IC50 (ng/mL)")         # Hwang 2023 Table 3: IC50 = 16.5 ng/mL (RSE 34.4%; bootstrap median 16.6, 95% CI 9.3-39.4)

    # ========================================================================
    # Inter-individual variability
    # Hwang 2023 equation (1) specifies an exponential (log-normal) IIV model,
    # P_i = P_TV * exp(eta_i). Tables 2 and 3 report the IIV as a coefficient
    # of variation in percent, so the internal variance is
    #     omega^2 = log(1 + CV^2).
    # Bracketed values below are the published eta shrinkages.
    # ========================================================================
    etaltlag ~ 0.2796176   # Hwang 2023 Table 2: omega_Tlag CV 56.8% (RSE 12.9%) [shrinkage 7%];  log(1 + 0.568^2) = 0.2796176
    etalcl   ~ 0.0210920   # Hwang 2023 Table 2: omega_CL/F CV 14.6% (RSE 9.8%)  [shrinkage 11%]; log(1 + 0.146^2) = 0.0210920
    etalvc   ~ 0.0659010   # Hwang 2023 Table 2: omega_Vc   CV 26.1% (RSE 9.2%)  [shrinkage 31%]; log(1 + 0.261^2) = 0.0659010
    etale0   ~ 0.0180609   # Hwang 2023 Table 3: omega_E0   CV 13.5% (RSE 12.7%) [shrinkage 0%];  log(1 + 0.135^2) = 0.0180609
    etalemax ~ 0.3560760   # Hwang 2023 Table 3: omega_Emax CV 65.4% (RSE 27%)   [shrinkage 8%];  log(1 + 0.654^2) = 0.3560760

    # ========================================================================
    # Residual error
    # Hwang 2023 equation (2) is proportional on the concentration and
    # equation (5) is additive on the heart rate. Tables 2 and 3 report the
    # NONMEM $SIGMA estimates on the variance scale, so both are converted to
    # standard deviations with sqrt(): an additive HR error of 64.2 is only
    # interpretable as a variance (sqrt = 8.01 bpm; 64.2 bpm as a standard
    # deviation would exceed the 60.4 bpm baseline), and the proportional
    # concentration error is reported in the same style in the sibling table.
    # See the vignette Assumptions and deviations section.
    # ========================================================================
    propSd   <- sqrt(0.379); label("Proportional residual SD on carvedilol Cc (fraction)")  # Hwang 2023 Table 2: proportional error 0.379 (RSE 2.2%; bootstrap median 0.376, 95% CI 0.339-0.404); sqrt(0.379) = 0.6156
    addSd_HR <- sqrt(64.2);  label("Additive residual SD on heart rate (bpm)")              # Hwang 2023 Table 3: additive error 64.2 (RSE 5.9%; bootstrap median 63.1, 95% CI 50.6-79.1); sqrt(64.2) = 8.012 bpm
  })

  model({
    # ---- Individual PK parameters ------------------------------------------
    # Hwang 2023 Table 2 footnote a: CL/F = 153 * (1 - 0.328 * GT).
    cl   <- exp(lcl + etalcl) * (1 - e_cyp2d6_star10_hom_cl * CYP2D6_STAR10_HOM)
    vc   <- exp(lvc + etalvc)
    vp   <- exp(lvp)
    q    <- exp(lq)
    d1   <- exp(ld1)
    tlag <- exp(ltlag + etaltlag)

    # ---- Micro-constants ---------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- PK ODE system ----------------------------------------------------
    # Oral dose enters the central compartment as a zero-order input of
    # duration d1 starting after the lag time tlag; because CL/F and Vc/F are
    # apparent, no separate depot compartment or bioavailability term is
    # needed (Hwang 2023 Fig. 1). Dose records must carry rate = -2 so that
    # rxode2 uses the modelled duration dur(central) = d1; without it the dose
    # collapses to an instantaneous bolus. Mean absorption time is
    # tlag + d1 / 2 = 0.407 h.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    dur(central)  <- d1
    alag(central) <- tlag

    # Central amount is in mg and vc in L, so central/vc is mg/L; multiply by
    # 1000 to report ng/mL, the unit in which IC50 was estimated.
    Cc <- (central / vc) * 1000

    # ---- Individual PD parameters ------------------------------------------
    e0   <- exp(le0 + etale0)
    emax <- exp(lemax + etalemax)
    ed50 <- exp(led50)
    ic50 <- exp(lic50)

    # ---- PD: competitive-antagonism direct-effect Emax ---------------------
    # Hwang 2023 Fig. 1 final-model equation (equation (4) with gamma = 1):
    #   HR = E0 + Emax * D / (ED50 * (1 + C / IC50) + D)
    # Carvedilol acts as a competitive antagonist, inflating the apparent
    # isoproterenol ED50 by the factor (1 + Cc / IC50) while leaving Emax
    # unchanged. With DOSE_ISOPROTERENOL_UG = 0 this reduces to HR = E0.
    HR <- e0 + emax * DOSE_ISOPROTERENOL_UG /
      (ed50 * (1 + Cc / ic50) + DOSE_ISOPROTERENOL_UG)

    # ---- Observations and residual error -----------------------------------
    Cc ~ prop(propSd)
    HR ~ add(addSd_HR)
  })
}

Biesdorf_2026_cofetuzumabPelidotin <- function() {
  description <- paste0(
    "Population pharmacokinetic model of the cofetuzumab pelidotin ",
    "(Cofe-P, ABBV-647, PF-06647020) CONJUGATE, the anti-PTK7 ",
    "antibody-drug conjugate carrying the Aur0101 auristatin payload, in ",
    "adults with advanced solid tumours (Biesdorf 2026, n = 202 pooled ",
    "from the first-in-human phase 1 study B7661001 and the phase 1b ",
    "NSCLC study M19-611). Linear two-compartment disposition after ",
    "intravenous infusion, with interindividual variability on ",
    "clearance and central volume and a combined proportional plus ",
    "additive residual error. Body weight enters as a fixed allometric ",
    "exponent of 0.75 on CL and Q and 1 on Vc and Vp; the phase 1b study ",
    "has a lower CL and Vc than the first-in-human study. The companion ",
    "unconjugated Aur0101 payload model, which Biesdorf 2026 developed ",
    "separately, is Biesdorf_2026_cofetuzumabPelidotin_aur0101; the five ",
    "companion exposure-response logistic models are the ",
    "Biesdorf_2026_cofetuzumabPelidotin_orr / _clinical_benefit_6mo / ",
    "_neutropenia_grade3 / _neuropathy_grade2 / _rash_grade2 family."
  )
  reference <- paste(
    "Biesdorf C, Rinas M, Engelhardt B, Saab R, Guo C, Ferlini C, Freise KJ,",
    "Menon RM, Polepally AR. Population Pharmacokinetics and Exposure-Response",
    "Analyses to Support Cofetuzumab Pelidotin (ABBV-647) Dose Optimization in",
    "Non-Small Cell Lung Cancer Patients. Clin Pharmacol Ther. 2026.",
    "doi:10.1002/cpt.70432. PMCID PMC13473755.",
    "All parameter values are from Table 1, 'Parameter estimates for Cofe-P",
    "(ADC) population pharmacokinetic model'.",
    sep = " "
  )
  vignette <- "Biesdorf_2026_cofetuzumabPelidotin"

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  # CL and Q are in L/day and Vc and Vp in L (Table 1), so an amount in mg
  # divided by a volume in L gives mg/L == ug/mL, the scale of the
  # Cofe-P Cavg axis in Biesdorf 2026 Figures 1-3 and Table 3 and of the
  # 'ADC Conjugate' panel of Figure S4.
  compartmentData <- list(
    central = list(analyte = "cofetuzumab pelidotin conjugate", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "cofetuzumab pelidotin conjugate", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Normalised to 70 kg inside model(). Biesdorf 2026 Methods: 'Continuous covariates were normalized to a reference value (median value of the covariate from the dataset) and incorporated in the model using a power function', but the dataset median weight is not printed anywhere; Table S3 gives only the mean (SD), 70.3 (18.2) kg overall. 70 kg is the rounded dataset mean and the standard allometric reference; see the vignette Assumptions and deviations section. The exponents are FIXED (Table 1, %RSE 0): 0.75 shared by CL and Q, 1 shared by Vc and Vp. Because dosing is weight-based (mg/kg), the fixed exponent of 1 on Vc makes the end-of-infusion concentration weight-independent.",
      source_name = "Weight"
    ),
    STUDY_M19611 = list(
      description = "1 = patient from the phase 1b study M19-611 (NCT04189614); 0 = patient from the first-in-human phase 1 study B7661001 (NCT02222922).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (first-in-human study B7661001)",
      notes = "Table 1 tabulates 'Study on CL' 0.794 and 'Study on Vc' 0.832 without naming the reference study. The phase 1b study is the one carrying the factors: (1) Table 3's median steady-state Cavg of 13.5 ug/mL at 2.8 mg/kg Q3W for the phase 1b target population is reproduced only when CL is 0.911 * 0.794 L/day (a reference-study CL of 0.911 L/day gives about 10.2 ug/mL); (2) the population-prediction panel of Figure S2 shows end-of-infusion bands near 80 and 65 ug/mL at 2.8 mg/kg, i.e. 2.8 mg/kg over 2.95 L/70 kg and over 0.832 * 2.95 L/70 kg, which confirms the size of the Vc factor (though not which study carries it). Set 0 to simulate the first-in-human population and 1 for the phase 1b population; Biesdorf 2026 Table 3 predictions are for phase 1b patients.",
      source_name = "Study"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 202L,
    n_studies = 2L,
    n_observations = "3391 Cofe-P (ADC) serum concentrations",
    age_range = "mean 61.1 years (SD 11.0) (Biesdorf 2026 Table S3)",
    weight_range = "mean 70.3 kg (SD 18.2) (Biesdorf 2026 Table S3)",
    sex_female_pct = 69.8,
    disease_state = "Advanced solid tumours. First-in-human B7661001 (n = 137): PTK7-selected triple-negative breast cancer and NSCLC and non-PTK7-selected platinum-resistant or refractory ovarian cancer. Phase 1b M19-611 (n = 65): recurrent NSCLC with PTK7 expression, enrolment narrowing to PTK7-high (>= 90% of cells at >= 2+ intensity), non-squamous, EGFR wild-type disease.",
    dose_range = "B7661001: 0.2-3.7 mg/kg every 3 weeks and 2.1-3.2 mg/kg every 2 weeks (days 1 and 15 of 28-day cycles). M19-611: 2.8 mg/kg every 3 weeks. All doses as a 60-minute intravenous infusion; 161 of 202 patients (79.7%) received 2.8 mg/kg Q3W.",
    ecog = "ECOG 0: 70 (34.7%), 1: 131 (64.9%), 2: 1 (0.5%) (Table S3)",
    liver_metastasis_pct = 23.8,
    albumin = "mean 37.2 g/L (SD 5.16) (Table S3)",
    notes = paste0(
      "Pooled data cut-off 24 May 2023. Concentrations below the ",
      "lower limit of quantitation (90 ng/mL for the conjugate) were ",
      "excluded. NONMEM 7.5, FOCE with interaction. A total-antibody ",
      "model was not developed because total-antibody exposure was ",
      "highly correlated with conjugate exposure (Figure S1)."
    )
  )

  ini({
    # ==================================================================
    # Structural disposition, Biesdorf 2026 Table 1, reported on the
    # natural scale with %RSE and 95% CI and therefore wrapped in log()
    # here. Typical values are for a 70 kg patient from the
    # first-in-human study (STUDY_M19611 = 0).
    # ==================================================================
    lcl <- log(0.911);  label("Clearance of the conjugate (L/day)")                   # Table 1 'CL (L/day)' 0.911, %RSE 4.14, 95% CI (0.840, 0.988)
    lvc <- log(2.95);   label("Central volume of distribution of the conjugate (L)")  # Table 1 'Vc (L)' 2.95, %RSE 2.83, 95% CI (2.79, 3.12)
    lq  <- log(0.0974); label("Intercompartmental clearance of the conjugate (L/day)") # Table 1 'Q (L/day)' 0.0974, %RSE 10.8, 95% CI (0.0789, 0.120)
    lvp <- log(1.09);   label("Peripheral volume of distribution of the conjugate (L)") # Table 1 'Vp (L)' 1.09, %RSE 6.56, 95% CI (0.955, 1.23)

    # ==================================================================
    # Covariate effects. Categorical covariates enter multiplicatively
    # and continuous ones as power functions of covariate / reference
    # (Biesdorf 2026 Methods). The two allometric exponents carry a
    # %RSE of 0 and the row label '(fixed)'.
    # ==================================================================
    e_study_m19611_cl <- 0.794;         label("Multiplicative factor on CL for the phase 1b study M19-611 versus the first-in-human study (unitless)") # Table 1 'Study on CL' 0.794, %RSE 7.11, 95% CI (0.691, 0.912)
    e_study_m19611_vc <- 0.832;         label("Multiplicative factor on Vc for the phase 1b study M19-611 versus the first-in-human study (unitless)") # Table 1 'Study on Vc' 0.832, %RSE 4.95, 95% CI (0.755, 0.916)
    e_wt_cl_q         <- fixed(0.75);   label("Allometric exponent on (WT/70 kg) for CL and Q (unitless)")  # Table 1 'Weight effect on CL and Q (fixed)' 0.750, %RSE 0
    e_wt_vc_vp        <- fixed(1);      label("Allometric exponent on (WT/70 kg) for Vc and Vp (unitless)") # Table 1 'Weight effect on Vc and Vp (fixed)' 1.00, %RSE 0

    # ==================================================================
    # Interindividual variability. The first column of the IIV rows is
    # the VARIANCE omega^2: the Table 1 footnote defines %CV as
    # SQRT[exp(omega^2) - 1] x 100, and sqrt(exp(0.211) - 1) = 48.4%
    # and sqrt(exp(0.0913) - 1) = 30.9% reproduce the printed %CV
    # column. No CL-Vc covariance is tabulated.
    # ==================================================================
    etalcl ~ 0.211   # Table 1 'IIV on CL' 0.211, 48.4 %CV, 2.17% shrinkage
    etalvc ~ 0.0913  # Table 1 'IIV on Vc' 0.0913, 30.9 %CV, 7.59% shrinkage

    # ==================================================================
    # Residual error. Table 1 prints the proportional and additive terms
    # as bare numbers. They are read as NONMEM $SIGMA VARIANCES and
    # entered as square roots, because (a) for N = 3391 observations the
    # expected %RSE of a variance is about 100*sqrt(2/N) = 2.4% and of
    # an SD about 1.2%, and the printed 3.20% matches the variance
    # reading; (b) the additive SD then comes to sqrt(0.0108) = 0.104
    # ug/mL, at the 0.090 ug/mL LLOQ, as an additive term usually does,
    # whereas the SD reading (0.0108 ug/mL) would sit at one-eighth of
    # the LLOQ. See the vignette Assumptions and deviations section.
    # ==================================================================
    propSd <- sqrt(0.118);  label("Proportional residual error on Cc (fraction)") # Table 1 'Proportional error' 0.118 (variance), %RSE 3.20, 95% CI (0.110, 0.125)
    addSd  <- sqrt(0.0108); label("Additive residual error on Cc (ug/mL)")        # Table 1 'Additive Error' 0.0108 (variance), %RSE 29.2, 95% CI (0.00463, 0.0170)
  })

  model({
    # ----- Individual parameters -----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * e_study_m19611_cl^STUDY_M19611
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp * e_study_m19611_vc^STUDY_M19611
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- Two-compartment disposition, intravenous infusion into central -----
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----- Observation -----
    # Cc is the conjugate serum concentration in ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}

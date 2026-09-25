Zhang_2019b_nivolumab <- function() {
  description <- "Two-compartment population PK model with sigmoidal time-varying clearance for intravenous nivolumab (anti-PD-1 IgG4) in Chinese and global patients with previously treated advanced solid tumors, including NSCLC and nasopharyngeal carcinoma (Zhang 2019, J Clin Pharmacol)"
  reference <- "Zhang J, Cai J, Bello A, Roy A, Sheng J. Model-Based Population Pharmacokinetic Analysis of Nivolumab in Chinese Patients With Previously Treated Advanced Solid Tumors, Including Non-Small Cell Lung Cancer. J Clin Pharmacol. 2019;59(10):1415-1424. doi:10.1002/jcph.1432"
  vignette <- "Zhang_2019b_nivolumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: the analyte and the specimen were both
  # confirmed against the source (Methods, 'PK Samples and Quantification of
  # Nivolumab Concentrations': samples were 'analyzed for nivolumab serum
  # concentrations' by validated ELISA / ECL ligand-binding assays).
  compartmentData <- list(
    central = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "nivolumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL (CLBBWT) and on Vc (VCBBWT) with reference weight 80 kg, taken from the Zhang 2019 Figure 1 caption ('The reference patient is a white/other male weighing 80 kg ...'). Q and Vp carry no body-weight term in this analysis, matching the Bajaj 2017 parent model that was re-estimated. Cohort body weight: mean 73.5 kg (SD 17.3), median 71.0 kg (range 34.9-157.9), 1 missing of 1200 (Table S6).",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL (CLGFR) with reference 90 mL/min/1.73 m^2 from the Zhang 2019 Figure 1 caption. Stored under the canonical CRCL; the source column is eGFR. Zhang 2019 does not restate which creatinine equation was used -- the covariate and its reference value are inherited from the Bajaj 2017 model whose parameters this analysis re-estimated, where eGFR was CKD-EPI. Cohort eGFR: mean 85.4 (SD 20.0), median 88.5 (range 31.1-135.4), 6 missing of 1200 (Table S6).",
      source_name = "eGFR"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Exponential effect on CL (CLSEX = -0.182) and on Vc (VCSEX = -0.132); both are negative, i.e. female patients have lower CL and lower Vc than the male reference. The polarity is fixed by the Figure 1 caption, which names the reference patient as male, and is corroborated by the Results text ('Female patients appeared to have a lower clearance and volume of the central compartment versus male patients'). Cohort: 388 female of 1200 (32.3%; Table S6).",
      source_name = "SEX"
    ),
    ECOG_GE1 = list(
      description = "Baseline Eastern Cooperative Oncology Group performance-status indicator (1 if ECOG PS > 0, else 0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status = 0, i.e. fully active)",
      notes = "Exponential effect on CL (CLPS = 0.138, i.e. exp(0.138) = 1.148 fold higher CL). Zhang 2019 collapses the ECOG scale to PS = 0 versus PS > 0; the Results text reports this as 'patients with Eastern Cooperative Oncology Group PS > 0 versus PS = 0 exhibited a 15% increase in nivolumab clearance', which is the back-transform of this coefficient. Cohort: PS 0 334 (27.8%), PS 1 858 (71.5%), PS 2 8 (0.7%) of 1200 (Table S6). Renamed from the source column PS to the canonical ECOG_GE1.",
      source_name = "PS"
    ),
    RACE_BLACK = list(
      description = "Race indicator: Black / African American",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (white or other; the Zhang 2019 reference patient is 'white/other')",
      notes = "Exponential effect on CL (CLRAAA = -0.00409). Effectively null and by far the least precisely estimated parameter in the model (RSE 1430%, 95% CI -0.120 to 0.105, spanning zero). Retained here because the published equation retains it. Cohort: 48 of 1200 (4.0%; Table S6). Renamed from the source column RAAA to the canonical RACE_BLACK.",
      source_name = "RAAA"
    ),
    RACE_ASIAN = list(
      description = "Race indicator: Asian (inclusive of the Chinese cohort)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (white or other; same composite reference as RACE_BLACK)",
      notes = "Exponential effect on CL (CLRAAS = -0.0891, i.e. exp(-0.0891) = 0.915, an 8.5% reduction reported in the Results text as 'a 9% lower baseline clearance for Asian versus non-Asian patients'). IMPORTANT: this indicator is 1 for the Chinese cohort as well as for non-Chinese Asian patients. Table S6 tabulates 'Chinese' (314) and 'Asian' (21) as separate REPORTING rows, but the modelled covariate RAAS is an Asian-race indicator covering both; setting RACE_ASIAN = 0 for Chinese subjects fails to reproduce the paper's own predicted Chinese baseline CL of 10.2 mL/h (Figure 2B). Renamed from the source column RAAS to the canonical RACE_ASIAN.",
      source_name = "RAAS"
    ),
    TUMTP_NPC = list(
      description = "Tumor-type indicator: nasopharyngeal carcinoma",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NSCLC, the reference tumor type)",
      notes = "Exponential effect on CL (CLNPC = 0.0889, i.e. exp(0.0889) = 1.093, reported in the Results text as 'a 9% higher baseline clearance for patients with NPC versus NSCLC'). Not statistically significant (RSE 90.9%, 95% CI -0.068 to 0.254, spanning zero) and interpreted with caution by the authors because of the small NPC sample. NPC was newly added as a categorical covariate in this analysis relative to the Bajaj 2017 parent model. Cohort: 23 of 1200 (1.9%; Table S6), all of them Chinese.",
      source_name = "NPC"
    ),
    TUMTP_OTHER = list(
      description = "Tumor-type indicator: tumor types other than NSCLC and NPC",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NSCLC, the reference tumor type)",
      notes = "Exponential effect on CL (CLOTH = 0.0718, i.e. exp(0.0718) = 1.074). Not statistically significant (RSE 57.5%, 95% CI -0.0129 to 0.153, spanning zero). Per-paper composition of the 'other tumor types' pool, enumerated in the Zhang 2019 Results text: renal cell carcinoma 35, colorectal cancer 33, prostate cancer 25, melanoma 116, hepatocellular carcinoma 2 -- i.e. 211 of 1200 (17.6%) once the 959 NSCLC (including 7 of unknown histology) and 23 NPC subjects are removed (Table S6). Renamed from the source column OTH to the canonical TUMTP_OTHER.",
      source_name = "OTH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1200L,
    n_studies = 7L,
    n_observations = 6945L,
    age_range = "Adults. Age is tabulated only for the Chinese cohort (Table S5): CheckMate 077 240 mg median 45 y (range 29-60), CheckMate 077 3 mg/kg median 52 y (range 27-59), CheckMate 078 3 mg/kg median 60 y (range 27-78). Not tabulated for the pooled global cohort.",
    weight_range = "34.9-157.9 kg (min-max, pooled cohort); mean 73.5 kg (SD 17.3), median 71.0 kg. Chinese cohort is lighter: mean 59.8 kg (SD 10.6) / 59.3 (8.2) / 64.0 (10.2) in the three Chinese study arms (Table S5).",
    sex_female_pct = 32.3,
    race_ethnicity = "Chinese 314 (26.2%), White 797 (66.4%), Black/African American 48 (4.0%), Asian (non-Chinese) 21 (1.8%), Other 12 (1.0%), Unknown 6 (0.5%), Missing 2 (0.2%) (Table S6). The analysis groups subjects as Chinese (n = 314), non-Chinese Asian (n = 21) and non-Asian (n = 865); the 'global population' of the paper is the 886 non-Chinese subjects.",
    renal_function = "Baseline eGFR mean 85.4 mL/min/1.73 m^2 (SD 20.0), median 88.5 (range 31.1-135.4), 6 missing (Table S6).",
    performance_status = "ECOG PS 0 334 (27.8%), PS 1 858 (71.5%), PS 2 8 (0.7%) (Table S6).",
    disease_state = "Previously treated advanced / recurrent solid tumors. Non-squamous NSCLC 544 (45.3%), squamous NSCLC 415 (34.6%), melanoma 116 (9.7%), renal cell carcinoma 35 (2.9%), colorectal cancer 33 (2.8%), prostate cancer 25 (2.1%), NPC 23 (1.9%), NSCLC of unknown histology 7 (0.6%), hepatocellular carcinoma 2 (0.2%) (Table S6).",
    dose_range = "Nivolumab 0.1, 0.3, 1, 3 or 10 mg/kg IV Q2W, plus a 240 mg Q2W flat-dose arm (Chinese only, n = 20). Infusions were administered over 60 minutes in the phase 1 studies (Table S1).",
    regions = "Two predominantly Chinese studies (CheckMate 077, CheckMate 078) and five global studies (MDX1106-01, CA209-003, CheckMate 017, CheckMate 057, CheckMate 063).",
    notes = "Baseline demographics per Tables S5 and S6 of the supporting information. This model was developed by RE-ESTIMATING the parameters of the Bajaj 2017 nivolumab monotherapy popPK model (reference 18 of Zhang 2019; shipped here as Bajaj_2017_nivolumab.R) on a pooled data set that adds the Chinese cohort, and by revising the tumor-type covariate to carry both NSCLC (reference) and NPC. The THETA indices in Table S7 are non-contiguous (theta1-4, 6, 7, 9, 12-15, 17, 18, 24-28) because they are inherited from the parent model's numbering, where the intervening THETAs were dropped during backward elimination."
  )

  ini({
    # Structural parameters. Reference patient (Zhang 2019 Figure 1 caption):
    # white/other male, 80 kg, baseline ECOG PS 0, eGFR 90 mL/min/1.73 m^2,
    # second-line-or-later NSCLC.
    # CL and Q are reported in mL/h in Table S7; converted to L/day (x 24 / 1000)
    # because this model keeps time in days, matching the rest of the nlmixr2lib
    # mAb library and its Bajaj 2017 / Zhang 2019 nivolumab siblings.
    lcl <- log(11.6 * 24 / 1000)
    label("Baseline clearance CLTV at the reference covariates (L/day)") # Zhang 2019 Table S7: CLTV (theta1) = 11.6 mL/h
    lvc <- log(4.19)
    label("Central volume of distribution VCTV (L)") # Zhang 2019 Table S7: VCTV (theta2) = 4.19 L
    lq <- log(29.3 * 24 / 1000)
    label("Intercompartmental clearance QTV (L/day)") # Zhang 2019 Table S7: QTV (theta3) = 29.3 mL/h
    lvp <- log(2.64)
    label("Peripheral volume of distribution VPTV (L)") # Zhang 2019 Table S7: VPTV (theta4) = 2.64 L

    # Covariate effects on CL (Zhang 2019 model equation, p. 1418):
    #   CL_i = CLTV * (BW_i/BW_TV)^CLBBWT * (eGFR_i/eGFR_TV)^CLGFR
    #          * e^CLSEX * e^CLPS * e^CLRAAA * e^CLRAAS * e^CLNPC * e^CLOTH
    # The exponential terms are written without their indicators in the printed
    # equation; each applies only to the subjects in that category, with the
    # reference patient of the Figure 1 caption carrying none of them.
    e_wt_cl <- 0.529
    label("Power exponent of WT on CL (unitless)") # Zhang 2019 Table S7: CLBBWT (theta7) = 0.529
    e_crcl_cl <- 0.132
    label("Power exponent of CRCL (eGFR) on CL (unitless)") # Zhang 2019 Table S7: CLGFR (theta9) = 0.132
    e_sexf_cl <- -0.182
    label("Exponential coefficient of female sex on CL (unitless)") # Zhang 2019 Table S7: CLSEX (theta12) = -0.182
    e_ecog_ge1_cl <- 0.138
    label("Exponential coefficient of ECOG PS > 0 on CL (unitless)") # Zhang 2019 Table S7: CLPS (theta13) = 0.138
    e_tumtp_npc_cl <- 0.0889
    label("Exponential coefficient of NPC tumor type on CL, versus NSCLC (unitless)") # Zhang 2019 Table S7: CLNPC (theta14) = 0.0889
    e_tumtp_other_cl <- 0.0718
    label("Exponential coefficient of other tumor types on CL, versus NSCLC (unitless)") # Zhang 2019 Table S7: CLOTH (theta15) = 0.0718
    e_race_black_cl <- -0.00409
    label("Exponential coefficient of Black race on CL (unitless)") # Zhang 2019 Table S7: CLRAAA (theta27) = -0.00409
    e_race_asian_cl <- -0.0891
    label("Exponential coefficient of Asian race on CL (unitless)") # Zhang 2019 Table S7: CLRAAS (theta28) = -0.0891

    # Covariate effects on Vc. Only body weight and sex were retained; Q and Vp
    # carry no covariates in this analysis (Table S7 has no QBBWT / VPBBWT row).
    e_wt_vc <- 0.740
    label("Power exponent of WT on VC (unitless)") # Zhang 2019 Table S7: VCBBWT (theta17) = 0.740
    e_sexf_vc <- -0.132
    label("Exponential coefficient of female sex on VC (unitless)") # Zhang 2019 Table S7: VCSEX (theta18) = -0.132

    # Time-varying clearance: sigmoidal Emax function of time since the start of
    # treatment (Methods, 'the current PPK model involved ... time-varying
    # clearance (sigmoidal function of the estimate of the maximal change in
    # clearance [Emax])'), in the form inherited from the Bajaj 2017 parent:
    #   CL(t) = CL_baseline * exp(cl_time_max * t^HILL / (T50^HILL + t^HILL))
    # cl_time_max is negative, so CL FALLS with time on treatment. The paper's
    # own back-transforms of this coefficient are the tightest available check:
    # 'Steady-state clearance (% of baseline) = exp(Emax) x 100%' = 68.5% and
    # 'Maximal decrease (% of baseline) = (1 - exp[Emax]) x 100%' = 31.5%, which
    # the Results report as a maximal decrease of 'approximately 32%'.
    # T50 is reported in hours in Table S7 and converted to days here.
    cl_time_max <- -0.378
    label("Maximal change in log-clearance at full time-effect saturation (unitless)") # Zhang 2019 Table S7: CLEMAX (theta24) = -0.378
    lcl_t50 <- log(1.38e3 / 24)
    label("Time at which half of the maximal change in CL is reached (log days)") # Zhang 2019 Table S7: CLt50 (theta25) = 1.38e3 h
    lcl_time_hill <- log(1.92)
    label("Hill / sigmoidicity exponent of time on CL (log unitless)") # Zhang 2019 Table S7: CLHILL (theta26) = 1.92

    # Inter-individual variability (Zhang 2019 Table S7, 'Random effects').
    # Estimates are VARIANCES; the value in parentheses in the source table is
    # the corresponding standard deviation (sqrt(0.119) = 0.345, sqrt(0.101) =
    # 0.318, sqrt(0.283) = 0.532, sqrt(0.0951) = 0.308), except on the ZCL:ZVC
    # row where the parenthesised 0.503 is the CORRELATION
    # (0.0551 / (0.345 * 0.318) = 0.502).
    # CL and VC form a 2x2 log-normal block; VP is an independent log-normal
    # eta; Emax carries an independent ADDITIVE eta on the linear scale, as in
    # the Bajaj 2017 parent model.
    etalcl + etalvc ~ c(
      0.119,
      0.0551, 0.101
    ) # Zhang 2019 Table S7: omega1,1 ZCL = 0.119; omega1,2 ZCL:ZVC = 0.0551; omega2,2 ZVC = 0.101
    etalvp ~ 0.283 # Zhang 2019 Table S7: omega3,3 ZVP = 0.283
    etacl_time_max ~ 0.0951 # Zhang 2019 Table S7: omega4,4 ZEMAX = 0.0951

    # Residual error. Zhang 2019 Methods states a proportional residual error
    # model, and Table S7 reports a single residual-error row, PEER (theta6) =
    # 0.224, carried as a THETA rather than a SIGMA. Read as the proportional
    # coefficient on the standard-deviation scale (22.4% CV), consistent with
    # the 0.215 reported the same way by the Bajaj 2017 parent model.
    propSd <- 0.224
    label("Proportional residual error (fraction)") # Zhang 2019 Table S7: PEER (theta6) = 0.224
  })
  model({
    # Individual baseline clearance, i.e. CL at t = 0, before any time effect.
    # Continuous covariates enter as (cov / reference)^exponent; categorical
    # covariates as exp(coefficient * indicator). Reference covariates:
    # WT 80 kg, eGFR 90 mL/min/1.73 m^2, male, ECOG PS 0, white/other race,
    # NSCLC (Zhang 2019 Figure 1 caption).
    cl_base <- exp(lcl + etalcl) *
      (WT / 80)^e_wt_cl *
      (CRCL / 90)^e_crcl_cl *
      exp(e_sexf_cl * SEXF) *
      exp(e_ecog_ge1_cl * ECOG_GE1) *
      exp(e_race_black_cl * RACE_BLACK) *
      exp(e_race_asian_cl * RACE_ASIAN) *
      exp(e_tumtp_npc_cl * TUMTP_NPC) *
      exp(e_tumtp_other_cl * TUMTP_OTHER)

    vc <- exp(lvc + etalvc) *
      (WT / 80)^e_wt_vc *
      exp(e_sexf_vc * SEXF)

    # Q and Vp are covariate-free in this analysis.
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Time-varying clearance. The driver t is time since the start of treatment,
    # so event tables must place the first nivolumab dose at t = 0.
    cl_t50 <- exp(lcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L, so central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

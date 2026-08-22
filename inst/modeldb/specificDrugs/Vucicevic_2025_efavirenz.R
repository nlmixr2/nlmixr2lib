Vucicevic_2025_efavirenz <- function() {
  description <- "One-compartment population pharmacokinetic-pharmacogenetic (PopPK-PGx) model for oral efavirenz in Caucasian HIV-1-infected adults (Vucicevic 2025), with binary CYP2B6 c.516G>T (rs3745274) carrier and CYP2B6 c.485-18C>T (rs4803419) TT-homozygote effects on apparent oral clearance CL/F. Absorption rate ka and apparent volume V/F were not estimable from the one-sample-per-patient sparse design and are fixed to literature values (0.3 1/h and 237 L)."
  reference <- paste(
    "Vucicevic K, Michalickova D, Obradovic B, Ranin J, Jevtovic D,",
    "Lukic R, Owen A, Dragovic G.",
    "Population pharmacokinetic-pharmacogenetic (PopPK-PGx) model of",
    "efavirenz in HIV-1-infected patients.",
    "Cureus. 2025;17(7):e88533.",
    "doi:10.7759/cureus.88533."
  )
  vignette <- "Vucicevic_2025_efavirenz"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant. Vucicevic 2025 encodes this as the binary indicator RSA ('whether the patient is a carrier of CYP2B6 516G>T (1) or not (0)', Table 2 footnote), which the canonical count column reproduces deterministically as `SNP_CYP2B6_RS3745274_T_COUNT >= 1`.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (germline genotype). Cohort distribution: 30/89 (33.71%) carriers (Table 1 row 'Carriers of CYP2B6 c.516G>T (rs3745274)'); Results paragraph 2 identifies the carriers as GT heterozygotes (n = 30) versus GG wild-type (n = 60), i.e. no TT homozygotes were observed, so in this dataset the carrier indicator is effectively a heterozygote indicator. (Results reports 30 + 60 = 90 genotyped subjects against the stated n = 89; see vignette Errata.) Because the published parameterization is a carrier indicator rather than a per-allele or HET/HOM decomposition, applying the model to a TT homozygote assigns the same -36.4% CL/F shift as to a heterozygote -- an extrapolation beyond the fitted cohort.",
      source_name        = "RSA (binary indicator: carrier of CYP2B6 516G>T)"
    ),
    SNP_CYP2B6_RS4803419_T_COUNT = list(
      description        = "Count of CYP2B6 c.485-18C>T (rs4803419, intron 3) T-alleles per subject (0/1/2). 0 = CC homozygous wild-type, 1 = CT heterozygous, 2 = TT homozygous variant. Vucicevic 2025 encodes this as the binary indicator RSB, which the canonical count column reproduces deterministically as `SNP_CYP2B6_RS4803419_T_COUNT == 2` (a recessive / TT-versus-rest encoding).",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (germline genotype). Cohort distribution: 12/89 (13.48%) flagged (Table 1). Although the Table 2 footnote calls RSB a 'carrier' indicator, three independent statements in the paper establish that the flagged group is the TT homozygous stratum, not any-T carriage: Results paragraph 2 ('a 26.8% decrease in CL/F in patients with TT genotype (n=12)') gives the same n = 12; Table 3 contrasts the simulated genotype strata as 'CYP2B6 c.485-18TT' versus 'CYP2B6 c.485-18CC'; and 13.48% is far below the ~45% any-T carrier frequency expected in a European-ancestry cohort while matching the expected TT-homozygote frequency. Heterozygotes (CT) are therefore pooled with CC wild-type in the reference group (recessive encoding, same shape as SNP_ABCG2_RS2231142_HOM).",
      source_name        = "RSB (binary indicator: CYP2B6 c.485-18C>T TT genotype)"
    )
  )

  # Covariates screened during stepwise covariate model building but NOT
  # retained in the final model (Vucicevic 2025 Methods 'PopPK-PGx model
  # development' paragraph 3 and Results paragraph 2). Documented so the
  # provenance of the covariate screen is preserved; none is referenced in
  # model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested on CL/F in linear and power functions during stepwise covariate model building; not retained in the final model. Cohort median 76 kg (IQR 66-85), Table 1. Note that V/F was fixed to a literature value quoted as '237 L/70 kg'; the 70 kg is the reference subject of the cited literature estimate, not a weight-scaling term in this model -- no weight term appears in the final-model equation of Table 2, and the covariate screen tested weight on CL/F only."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Tested on CL/F in linear and power functions; not retained in the final model. No cohort distribution reported."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F as a categorical covariate; not retained in the final model. Table 1 reports 'Sex, male; n (%) 21 (23.60)', i.e. 76.4% female; see vignette Errata."
    ),
    SMOKE = list(
      description = "Current-smoker indicator (1 = current smoker, 0 = non-smoker)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F and entered the model during forward selection, but was removed in the backward-elimination step (Results paragraph 2: 'Smoking status was excluded from the model in the backward step of covariate model building'). Cohort prevalence 60/89 (67.41%), Table 1. Smoking WAS associated with efavirenz concentrations in the source cohort's earlier single-timepoint analysis (Olagunju 2014, reference 17)."
    ),
    SNP_NR1I3_RS2307424_T_COUNT = list(
      description = "Count of NR1I3 (CAR) c.540C>T (rs2307424) T-alleles per subject (0/1/2)",
      units       = "(count, 0/1/2)",
      type        = "continuous",
      notes       = "Screened as a categorical covariate on CL/F; no significant effect detected (Discussion paragraph 4: 'No impact of these polymorphisms were observed in the current analysis'). Cohort: 53/89 (59.60%) carriers, Table 1. No point estimate is published, so no effect can be encoded."
    ),
    SNP_NR1I3_RS3003596_C_COUNT = list(
      description = "Count of NR1I3 (CAR) c.152-1089T>C (rs3003596) C-alleles per subject (0/1/2)",
      units       = "(count, 0/1/2)",
      type        = "continuous",
      notes       = "Screened as a categorical covariate on CL/F; no significant effect detected (Discussion paragraph 4). Cohort: 58/89 (65.17%) carriers, Table 1. No point estimate is published, so no effect can be encoded."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 89,
    n_studies      = 1,
    age_range      = "18 years and older (inclusion criterion); interquartile range 32-48 years",
    age_median     = "40 years (IQR 32-48)",
    weight_range   = "interquartile range 66-85 kg (full range not reported)",
    weight_median  = "76 kg (IQR 66-85)",
    sex_female_pct = 76.40,
    race_ethnicity = c(Caucasian = 100),
    disease_state  = "Adults with confirmed HIV-1 infection receiving efavirenz-based antiretroviral therapy for at least three months. Excluded: concomitant antituberculosis medication or other drugs known to interact significantly with efavirenz metabolism, age under 18 years, pregnancy, and incomplete pharmacokinetic or clinical data.",
    dose_range     = "600 mg orally once daily (all subjects)",
    regions        = "Serbia (outpatient HIV/AIDS Center, Clinic for Infectious and Tropical Diseases, University of Belgrade Teaching Hospital)",
    smoking        = "60/89 (67.41%) smokers (Table 1)",
    cyp2b6_freq    = "CYP2B6 c.516G>T (rs3745274): 30/89 (33.71%) carriers, identified in Results paragraph 2 as GT heterozygotes (n = 30) versus GG wild-type (n = 60), with no TT homozygotes. CYP2B6 c.485-18C>T (rs4803419): 12/89 (13.48%) flagged, identified in Results paragraph 2 as the TT genotype stratum (Table 1, Table 3).",
    nr1i3_freq     = "NR1I3 c.540C>T (rs2307424): 53/89 (59.60%) carriers. NR1I3 c.152-1089T>C (rs3003596): 58/89 (65.17%) carriers. Neither polymorphism had a detectable effect on efavirenz CL/F in this analysis (Table 1, Discussion paragraph 4).",
    notes          = "Sparse design: exactly ONE steady-state efavirenz plasma sample per patient, drawn 2.83-13.42 h after the last dose (Results paragraph 1; Appendix 1 Figure 4). Data are a re-analysis of the cohort first reported in Olagunju 2014 (reference 17). NONMEM 7.4 / PsN 3.4.2 under Pirana 2.9.0, FOCE-I estimation, ADVAN2 TRANS2. Internal validation by 1000-replicate bootstrap and NPDE (mean 0.058, variance 1.009; both not significantly different from 0 and 1)."
  )

  ini({
    # ---- Structural PK parameters ----
    # Vucicevic 2025 Methods 'PopPK-PGx model development' paragraph 3 and
    # Discussion paragraph 2: with only one sample per patient, ka and V/F
    # "could not be estimated" and were fixed to literature-reported values.
    lka <- fixed(log(0.3)) ; label("First-order absorption rate constant ka (1/h)")  # Vucicevic 2025 Methods para 3 / Discussion para 2: ka fixed to 0.3 1/h, cited to Csajka 2003 (reference 20)
    lvc <- fixed(log(237)) ; label("Apparent volume of distribution V/F (L)")        # Vucicevic 2025 Methods para 3 / Discussion para 2: V/F fixed to "237 L/70 kg" (references 9, 11, 18, 19), i.e. the literature estimate for a ~70 kg subject. No weight term appears in the Table 2 final-model equation and weight was screened on CL/F only, so V/F is encoded as a constant (see covariatesDataExcluded$WT and the vignette Errata).
    lcl <- log(13.9)       ; label("Population apparent oral clearance CL/F (L/h) -- CYP2B6 516GG, c.485-18 non-TT reference") # Vucicevic 2025 Table 2 CLp = 13.9 L/h (RSE 4%; bootstrap median 13.9, 95% CI 12.73-15.21)

    # ---- Covariate effects on CL/F (Vucicevic 2025 Table 2 final-model equation) ----
    #   CL/F = CLp * (1 + theta_RSA)^RSA * (1 + theta_RSB)^RSB
    # RSA and RSB are 0/1 indicators, so each factor collapses to 1 when the
    # indicator is 0 and to (1 + theta) when it is 1.
    e_cyp2b6_516t_cl  <- -0.364 ; label("Fractional change in CL/F for CYP2B6 516G>T carriers (unitless)")   # Vucicevic 2025 Table 2 theta_RSA = -0.364 (RSE 14%; bootstrap median -0.368, 95% CI -0.457 to -0.255); Results para 2: CL/F 36.4% lower in GT versus GG
    e_cyp2b6_485tt_cl <- -0.268 ; label("Fractional change in CL/F for CYP2B6 c.485-18C>T TT homozygotes (unitless)") # Vucicevic 2025 Table 2 theta_RSB = -0.268 (RSE 16%; no bootstrap value reported); Results para 2: 26.8% decrease in CL/F in the TT genotype

    # ---- IIV (log-normal on CL/F; the only random effect in the model) ----
    # Vucicevic 2025 Table 2 "CV CL (%) 13.1 (35)" and Results paragraph 2:
    # "IIV in CL/F was best described by an exponential model, and the
    # estimated value was 13.1% (35% RSE)".
    # CV-to-variance conversion: omega^2 = log(CV^2 + 1) = log(0.131^2 + 1).
    etalcl ~ 0.01702  # log(0.131^2 + 1) = 0.017016; bootstrap median CV 13.4%, 95% CI 3.24-20.9%

    # ---- Residual error (proportional) ----
    propSd <- 0.25 ; label("Proportional residual error (fraction)")  # Vucicevic 2025 Table 2 "Proportional error 0.25 (11)"; Results para 2: "residual variability was explained by a standard deviation of 0.25 (11%)"; bootstrap median 0.244, 95% CI 0.182-0.299
  })

  model({
    # 1. Recover the paper's two binary genotype indicators from the canonical
    #    0/1/2 allele-count covariate columns.
    #    RSA = carrier of CYP2B6 516G>T (Table 2 footnote: "carrier ... (1) or not (0)").
    #    RSB = CYP2B6 c.485-18C>T TT homozygote (Results paragraph 2; Table 3
    #          contrasts "c.485-18TT" against "c.485-18CC").
    rsa <- (SNP_CYP2B6_RS3745274_T_COUNT >= 1)
    rsb <- (SNP_CYP2B6_RS4803419_T_COUNT == 2)

    # 2. Typical CL/F, written exactly as the Table 2 final-model equation:
    #    CL/F = CLp * (1 + theta_RSA)^RSA * (1 + theta_RSB)^RSB
    tvcl <- exp(lcl) *
      (1 + e_cyp2b6_516t_cl)^rsa *
      (1 + e_cyp2b6_485tt_cl)^rsb

    # 3. Individual PK parameters
    ka <- exp(lka)
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc)

    # 4. Micro-constants
    kel <- cl / vc

    # 5. ODE system: one-compartment with first-order absorption (ADVAN2 TRANS2)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

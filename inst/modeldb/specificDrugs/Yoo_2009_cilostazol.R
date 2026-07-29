Yoo_2009_cilostazol <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral cilostazol with ",
    "first-order absorption from the depot and an absorption lag time, ",
    "estimated in 104 healthy Korean male volunteers receiving a single ",
    "50- or 100-mg dose (Yoo 2009). Apparent oral clearance CL/F is ",
    "modulated by two pharmacogenetic covariates entered in linear-",
    "fractional form: a three-level CYP3A5 genotype (CYP3A5*1/*1 = ",
    "reference, *1/*3 = -22.3%, *3/*3 = -40.7%) and a three-level ",
    "CYP2C19 metabolizer phenotype (extensive metabolizer = reference, ",
    "intermediate = -14.7%, poor = -27.2%). The final NONMEM ADVAN4/",
    "TRANS4 model places exponential IIV on CL/F, Vc/F, Q/F and Vp/F ",
    "with a partial OMEGA BLOCK retaining the (Vp/F, CL/F), (Q/F, Vc/F) ",
    "and (Vp/F, Q/F) covariances; the remaining off-diagonals are held ",
    "at zero. Residual error is combined additive plus proportional."
  )
  reference   <- paste0(
    "Yoo HD, Cho HY, Lee YB. Population pharmacokinetic analysis of ",
    "cilostazol in healthy subjects with genetic polymorphisms of ",
    "CYP3A5, CYP2C19 and ABCB1. Br J Clin Pharmacol. 2010;69(1):27-37. ",
    "doi:10.1111/j.1365-2125.2009.03558.x. Online: 2009-12-04."
  )
  vignette    <- "Yoo_2009_cilostazol"
  units       <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CYP3A5_STAR1_HET = list(
      description        = "CYP3A5*1/*3 heterozygote indicator (one functional CYP3A5*1 allele).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*1/*1 homozygote-expresser in Yoo 2009, paired with CYP3A5_STAR1_HOM = 1; or CYP3A5*3/*3 nonexpresser when both indicators are 0).",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = CYP3A5*1/*3 ",
        "heterozygote (one *1 allele at rs776746); 0 = otherwise. Yoo ",
        "2009 cohort (Table 1, n = 104): *1/*1 6/104 (5.8%), *1/*3 ",
        "42/104 (40.4%), *3/*3 56/104 (53.8%). Paired with ",
        "CYP3A5_STAR1_HOM. Yoo 2009 places the reference genotype at ",
        "*1/*1 (highest CL/F), so the linear-fractional coefficients in ",
        "the final-model equation are negative for both the *1/*3 ",
        "(CYP3A5_STAR1_HET = 1, CYP3A5_STAR1_HOM = 0) and *3/*3 ",
        "(both indicators 0) strata. The *3/*3 indicator is derived ",
        "inside model() as 1 - CYP3A5_STAR1_HET - CYP3A5_STAR1_HOM."
      ),
      source_name        = "CYP3A5*1/*3"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5*1/*1 homozygote indicator (two functional CYP3A5*1 alleles).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (CYP3A5*1/*1 homozygote-expresser is the Yoo 2009 typical-value reference for CL/F).",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = CYP3A5*1/*1 ",
        "homozygote (two *1 alleles at rs776746); 0 = otherwise. Yoo ",
        "2009 cohort: 6/104 (5.8%). In Yoo 2009 the *1/*1 stratum is ",
        "the reference for the CYP3A5 covariate effect on CL/F (Table ",
        "2 Model 4, see source-trace table in the vignette). Paired ",
        "with CYP3A5_STAR1_HET; together they uniquely identify the ",
        "three-level CYP3A5 genotype."
      ),
      source_name        = "CYP3A5*1/*1"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C19 extensive metabolizer, *1/*1, when paired with CYP2C19_PM = 0).",
      notes              = paste0(
        "Time-fixed per subject (germline genotype-derived phenotype). ",
        "1 = subject is a CYP2C19 intermediate metabolizer (one ",
        "functional and one loss-of-function allele; Yoo 2009 ",
        "definition: CYP2C19*1/*2 or *1/*3); 0 = otherwise. Yoo 2009 ",
        "cohort (Methods 'Genotype analysis' and Results): EM (*1/*1) ",
        "52/104 (50.0%), IM 38/104 (36.5%), PM 14/104 (13.5%). ",
        "Paired with CYP2C19_PM to encode the three-level EM ",
        "(reference) / IM / PM phenotype with two binary indicators ",
        "following the Zhao 2018 omeprazole precedent."
      ),
      source_name        = "CYP2C19 IM"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C19 extensive metabolizer, *1/*1, when paired with CYP2C19_IM = 0).",
      notes              = paste0(
        "Time-fixed per subject (germline genotype-derived phenotype). ",
        "1 = subject is a CYP2C19 poor metabolizer (two loss-of-",
        "function alleles; Yoo 2009 definition: CYP2C19*2/*2, *2/*3 ",
        "or *3/*3); 0 = otherwise. Yoo 2009 cohort: PM 14/104 (13.5%). ",
        "Paired with CYP2C19_IM."
      ),
      source_name        = "CYP2C19 PM"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Adult age (years).",
      units       = "years",
      type        = "continuous",
      notes       = paste0(
        "Tested with stepwise covariate search in PsN; not retained in ",
        "the Yoo 2009 final model (Results: 'Other covariates including ",
        "age, body weight, BSA and the ABCB1 genotype did not appear to ",
        "affect any of the PK parameters'). Cohort range 19-28 years, ",
        "mean 23.7 +/- 1.6 (Methods 'Subjects')."
      )
    ),
    WT = list(
      description = "Body weight (kg).",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Tested via stepwise covariate search; not retained. Cohort ",
        "range 44-83.5 kg, mean 65.9 +/- 7.7 (Methods 'Subjects')."
      )
    ),
    BSA = list(
      description = "Body surface area (m^2).",
      units       = "m^2",
      type        = "continuous",
      notes       = paste0(
        "Tested via stepwise covariate search; not retained. Cohort ",
        "range 1.420-2.017 m^2, mean 1.783 +/- 0.124 (Methods ",
        "'Subjects')."
      )
    ),
    CYP3A4_STAR1B = list(
      description = "CYP3A4*1B allele indicator (rs2740574, -392A>G).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Genotyped but not detected in any of the 104 Yoo 2009 subjects ",
        "(Results 'Genetic analysis'); allele is rare in Asians. Not ",
        "carried as a model covariate."
      )
    ),
    ABCB1_C1236T_HET = list(
      description = "ABCB1 exon 12 C1236T heterozygote indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Tested; not retained. Yoo 2009 cohort: CC 28/104 (26.9%), CT ",
        "50/104 (48.1%), TT 26/104 (25.0%)."
      )
    ),
    ABCB1_G2677T_HET = list(
      description = "ABCB1 exon 21 G2677T/A heterozygote indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Tested; not retained. Yoo 2009 cohort: GG 32/104 (30.8%), ",
        "GT/GA 56/104 (53.9%), TT/TA/AA 16/104 (15.3%)."
      )
    ),
    ABCB1_C3435T_HET = list(
      description = "ABCB1 exon 26 C3435T heterozygote indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Tested; not retained. Yoo 2009 cohort: CC 50/104 (48.1%), CT ",
        "42/104 (40.4%), TT 12/104 (11.5%). Discussion: 'no significant ",
        "differences were observed in the cilostazol PK parameters ",
        "among the ABCB1 genotype groups'."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 104L,
    n_studies       = 4L,
    age_range       = "19-28 years (Methods 'Subjects'; mean 23.7 +/- 1.6).",
    age_median      = "23.7 years (mean reported in Methods 'Subjects').",
    weight_range    = "44-83.5 kg (mean 65.9 +/- 7.7).",
    weight_median   = "65.9 kg (mean reported in Methods 'Subjects').",
    sex_female_pct  = 0,
    race_ethnicity  = c(Asian = 100),
    disease_state   = paste0(
      "Healthy adult volunteers screened with physical exam, blood ",
      "chemistry, complete blood count and urinalysis; no history of ",
      "illness or drug hypersensitivity; refrained from medications, ",
      "alcohol and other drugs for at least 1 week prior to and ",
      "throughout the study (Methods 'Subjects')."
    ),
    dose_range      = paste0(
      "Single oral dose of 50 mg (n = 52) or 100 mg (n = 52) cilostazol ",
      "(Pletaal tablet, Otsuka Pharmaceuticals) with 240 mL water after ",
      "an overnight fast (Methods 'Study design')."
    ),
    regions         = "Korea (single-centre, Chonnam National University, Gwangju).",
    sampling_design = paste0(
      "Serum samples collected pre-dose and at 1, 2, 2.5, 3, 3.5, 4, 6, ",
      "8, 12, 24 and 48 h post-dose (12 samples per subject across ",
      "0-48 h; Methods 'Study design'). Pooled from four separate ",
      "bioequivalence (BE) studies sharing the same protocol; only ",
      "reference-formulation data were retained for this analysis."
    ),
    cyp3a5_genotype = c(`*1/*1` = 5.8, `*1/*3` = 40.4, `*3/*3` = 53.8),
    cyp2c19_phenotype = c(EM = 50.0, IM = 36.5, PM = 13.5),
    notes           = paste0(
      "Demographics from Yoo 2009 Methods 'Subjects' and Table 1 ",
      "(genotype frequencies). Software: NONMEM v6 level 1.1 (FOCE-I ",
      "with eta-epsilon interaction); Wings for NONMEM v614 for ",
      "bootstrap; PsN v2.3.1 for stepwise covariate selection; Xpose ",
      "v4.0.4 for diagnostics. Stepwise significance levels: P < 0.05 ",
      "forward, P < 0.01 backward."
    )
  )

  ini({
    # ===== Structural parameters (Yoo 2009 Table 3, "Final model Estimate") =====
    # Apparent oral clearance CL/F at the CYP3A5*1/*1 + CYP2C19 EM
    # combined reference (the typical-value subject); 12.8 L/h (95% CI
    # 8.2, 17.4; %RSE 13.9). Bootstrap median 12.7 (CI 8.9, 18.5).
    lcl   <- log(12.8);    label("Apparent oral clearance CL/F at the CYP3A5*1/*1 / CYP2C19 EM reference (L/h)")  # Yoo 2009 Table 3, TV(CL/F)
    lvc   <- log(20.5);    label("Apparent central volume of distribution Vc/F (L)")                              # Yoo 2009 Table 3, TV(Vc/F)
    lq    <- log(5.64);    label("Apparent inter-compartmental clearance Q/F (L/h)")                              # Yoo 2009 Table 3, TV(Q/F)
    lvp   <- log(73.1);    label("Apparent peripheral volume of distribution Vp/F (L)")                            # Yoo 2009 Table 3, TV(Vp/F)
    lka   <- log(0.244);   label("First-order absorption rate constant Ka (1/h)")                                  # Yoo 2009 Table 3, TV(Ka)
    ltlag <- log(0.565);   label("Absorption lag time Tlag (h)")                                                   # Yoo 2009 Table 3, TV(TLag)

    # ===== CYP3A5 covariate effects on CL/F (Yoo 2009 final-model equation) =====
    # Linear-fractional form: CL/F = TVCL * (1 + e_*1/*3 * G_HET +
    # e_*3/*3 * G_MUT + e_IM * G_IM + e_PM * G_PM) where the
    # G_HET / G_MUT / G_IM / G_PM are paired binary indicators
    # (Methods 'Covariates analysis' equation for categorical
    # covariates; final equation reproduced in Results just below
    # Table 3). Reference combination: CYP3A5*1/*1 + CYP2C19 EM
    # (typical CL/F = 12.8 L/h).
    e_cyp3a5_star1_het_cl <- -0.223; label("Fractional change in CL/F for CYP3A5*1/*3 heterozygote (unitless)")           # Yoo 2009 Table 3, CL/F theta *1/*3
    e_cyp3a5_star3_hom_cl <- -0.407; label("Fractional change in CL/F for CYP3A5*3/*3 homozygous nonexpresser (unitless)")  # Yoo 2009 Table 3, CL/F theta *3/*3

    # ===== CYP2C19 covariate effects on CL/F =====
    e_cyp2c19_im_cl       <- -0.147; label("Fractional change in CL/F for CYP2C19 intermediate metabolizer (unitless)")    # Yoo 2009 Table 3, CL/F theta IMs
    e_cyp2c19_pm_cl       <- -0.272; label("Fractional change in CL/F for CYP2C19 poor metabolizer (unitless)")            # Yoo 2009 Table 3, CL/F theta PMs

    # ===== Inter-individual variability (Yoo 2009 Table 3 OMEGA BLOCK) =====
    # Exponential (log-normal) model: X_i = TV(X) * exp(eta_i). The
    # OMEGA matrix is reported as variances on the eta scale; the
    # paper expresses CV% as CV = sqrt(omega^2) (Discussion: "108%
    # for Q/F and 82.9% for V3/F" matches sqrt(1.16) = 1.077 and
    # sqrt(0.688) = 0.829). Partial OMEGA BLOCK on CL, Vc, Q, Vp;
    # off-diagonals for (CL, Vc), (CL, Q), (Vc, Vp) are held at zero
    # (Methods: correlations were 'explored using the OMEGA BLOCK
    # option'; Table 3 reports only the three non-zero off-diagonals).
    # No IIV on Ka or Tlag (Results paragraph 1: 'the effect of
    # including a random effect on KA and ALAG1 was not significant').
    #
    # Order in the c(...) lower-triangular: (var_cl, cov_cl_vc, var_vc,
    # cov_cl_q, cov_vc_q, var_q, cov_cl_vp, cov_vc_vp, cov_q_vp, var_vp).
    etalcl + etalvc + etalq + etalvp ~ c(
      0.075,
      fixed(0), 0.481,
      fixed(0), -0.355, 1.16,
      0.0685, fixed(0), 0.460, 0.688)                                                                          # Yoo 2009 Table 3, omega^2 / off-diagonals

    # ===== Residual error (Yoo 2009 Table 3; combined additive + proportional) =====
    # The paper reports variances sigma^2_pro = 0.0131 (SD = 0.114
    # fraction = 11.4% proportional) and sigma^2_add = 0.00261 (SD =
    # 0.0511 ug/mL).
    propSd <- 0.114;   label("Proportional residual error (fraction)")        # Yoo 2009 Table 3, sigma^2_pro = 0.0131 -> SD = 0.114
    addSd  <- 0.0511;  label("Additive residual error (ug/mL)")               # Yoo 2009 Table 3, sigma^2_add = 0.00261 -> SD = 0.0511
  })

  model({
    # ----- 1. Derived CYP3A5 indicator for the *3/*3 stratum -----
    # The dataset carries two registered binary indicators
    # CYP3A5_STAR1_HET (= 1 if *1/*3) and CYP3A5_STAR1_HOM (= 1 if
    # *1/*1). The *3/*3 stratum is implicit when both are 0; derive
    # it here so the linear-fractional CL/F equation reproduces Yoo
    # 2009 Methods 'Covariates analysis' exactly.
    cyp3a5_star3_hom <- 1 - CYP3A5_STAR1_HET - CYP3A5_STAR1_HOM

    # ----- 2. Individual parameters -----
    ka   <- exp(lka)
    tlag <- exp(ltlag)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq  + etalq)
    vp   <- exp(lvp + etalvp)
    cl   <- exp(lcl + etalcl) * (1 +
              e_cyp3a5_star1_het_cl * CYP3A5_STAR1_HET +
              e_cyp3a5_star3_hom_cl * cyp3a5_star3_hom +
              e_cyp2c19_im_cl       * CYP2C19_IM +
              e_cyp2c19_pm_cl       * CYP2C19_PM)

    # ----- 3. Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- 4. ODE system (NONMEM ADVAN4 / TRANS4 equivalent) -----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- 5. Absorption lag time on the depot -----
    alag(depot) <- tlag

    # ----- 6. Observation and combined additive + proportional error -----
    # Cc has units of mg/L = ug/mL (depot amount in mg, vc in L), so
    # propSd is dimensionless and addSd has units ug/mL matching the
    # reported assay scale.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

Khwarg_2024_proguanil <- function() {
  description <- paste(
    "Semi-physiologic joint parent-metabolite population PK model for oral",
    "proguanil and its active metabolite cycloguanil in healthy Korean male",
    "adults, quantifying the effect of the SLC22A1 (OCT1) 1022C>T",
    "(rs2282143) polymorphism on hepatic uptake (Khwarg 2024). Proguanil is",
    "absorbed first-order from the `depot` directly into a permeability-",
    "considered well-stirred `liver` compartment, so hepatic first-pass",
    "metabolism is structural rather than folded into a bioavailability",
    "term. The liver exchanges with `central` at the hepatic plasma flow QH",
    "(inflow QH, return flow QH * (1 - EH)), and `central` additionally",
    "exchanges with `peripheral1` at Q/F and is cleared at CL/F, the",
    "apparent clearance via the non-cycloguanil pathway. The hepatic",
    "extraction ratio EH follows the well-stirred model on a blood basis",
    "from the blood unbound fraction fub, the intrinsic clearance CLint and",
    "the hepatic blood flow QBH; the resulting hepatic plasma clearance",
    "CLH = EH * QBH * (CB/CP) carries proguanil out of the liver and into",
    "the cycloguanil chain. Cycloguanil is formed in `liver_cycloguanil`",
    "and effluxes to `central_cycloguanil` through two transit compartments",
    "with a shared rate constant mKT, then is cleared at CLM/F; its central",
    "volume was set equal to the proguanil central volume because the two",
    "were not separately identifiable. Liver volume (1 L) and hepatic blood",
    "flow (90 L/h) are fixed physiologic constants for a 70 kg adult,",
    "allometrically scaled on body weight with exponents 1 and 0.75; the",
    "plasma unbound fraction (0.25), blood-to-plasma ratio (2.7) and",
    "hematocrit (0.45) are fixed literature values. The SLC22A1 1022C>T",
    "heterozygote (CT) genotype enters as an exponent on FUP = 0.416, the",
    "relative fraction of OCT1-mediated hepatocyte uptake versus the CC",
    "wild type, reducing EH and thereby raising proguanil and lowering",
    "cycloguanil exposure."
  )
  reference <- paste(
    "Khwarg J, Yang E, Park CS, Ji SC, Yu K-S, Lee S. Effect of SLC22A1",
    "polymorphism on the pharmacokinetics of proguanil in Korean: A",
    "semi-physiologic population pharmacokinetic approach. Clin Transl Sci.",
    "2024;17:e70103. doi:10.1111/cts.70103"
  )
  vignette <- "Khwarg_2024_proguanil"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")
  compartmentData <- list(
    depot                 = list(analyte = "proguanil",   units = "mg", specimen = "administration site", verified = TRUE),
    liver                 = list(analyte = "proguanil",   units = "mg", specimen = "tissue",              verified = TRUE),
    central               = list(analyte = "proguanil",   units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1           = list(analyte = "proguanil",   units = "mg", specimen = "plasma",              verified = TRUE),
    liver_cycloguanil     = list(analyte = "cycloguanil", units = "mg", specimen = "tissue",              verified = TRUE),
    transit1_cycloguanil  = list(analyte = "cycloguanil", units = "mg", specimen = "tissue",              verified = TRUE),
    transit2_cycloguanil  = list(analyte = "cycloguanil", units = "mg", specimen = "tissue",              verified = TRUE),
    central_cycloguanil   = list(analyte = "cycloguanil", units = "mg", specimen = "plasma",              verified = TRUE)
  )
  covariateData <- list(
    WT = list(
      description = paste(
        "Total body weight, used only as the allometric size descriptor for",
        "the fixed physiologic liver volume and hepatic blood flow"
      ),
      units = "kg",
      type = "continuous",
      notes = paste(
        "Khwarg 2024 Methods 'Structural model' Equations 5 and 6 scale the",
        "70 kg reference liver volume (1 L, exponent 1) and hepatic blood",
        "flow (90 L/h, exponent 0.75) on body weight. Body weight was also",
        "screened by stepwise covariate modeling on CL/F, CLint and the",
        "central volume and was NOT retained (Results 'Covariate analysis':",
        "'no significant covariates were identified'), so WT acts here only",
        "through the physiologic sub-model. Cohort mean 76.12 kg (SD 9.82).",
        "Time-fixed per subject."
      ),
      source_name = "Body weight (kg)"
    ),
    SNP_SLC22A1_RS2282143 = list(
      description = paste(
        "SLC22A1 (OCT1) 1022C>T (rs2282143) variant carrier indicator;",
        "1 = CT heterozygote, 0 = CC homozygous wild type"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CC homozygous wild type)",
      notes = paste(
        "Khwarg 2024 Results 'Effect of SLC22A1 genetic polymorhism on OCT1",
        "activity': 'The SLC22A1 genotype (CC genotype, GENO = 0; CT",
        "genotype, GENO = 1) was applied as an exponent of FUP.' No TT",
        "homozygotes were observed in the cohort (Table 1: 9 CC, 4 CT of 13",
        "genotyped), so the indicator is heterozygote-vs-wild-type rather",
        "than a general variant-allele count. Time-fixed per subject",
        "(germline genotype)."
      ),
      source_name = "GENO"
    )
  )
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate (eGFR)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = paste(
        "Screened by stepwise covariate modeling on CL/F but not retained",
        "(Khwarg 2024 Results 'Covariate analysis'). No point estimate is",
        "reported, so the effect cannot be encoded."
      )
    ),
    AST = list(
      description = "Aspartate transaminase",
      units = "U/L",
      type = "continuous",
      notes = paste(
        "Screened by stepwise covariate modeling on CLint but not retained",
        "(Khwarg 2024 Results 'Covariate analysis'). No point estimate is",
        "reported, so the effect cannot be encoded."
      )
    ),
    ALT = list(
      description = "Alanine transaminase",
      units = "U/L",
      type = "continuous",
      notes = paste(
        "Screened by stepwise covariate modeling on CLint but not retained",
        "(Khwarg 2024 Results 'Covariate analysis'). No point estimate is",
        "reported, so the effect cannot be encoded."
      )
    )
  )
  population <- list(
    species = "human",
    n_subjects = 16,
    n_studies = 1,
    age_range = "19-50 years (eligibility); cohort mean 34.31 y (SD 7.38)",
    weight_range = "cohort mean 76.12 kg (SD 9.82); BMI 19.0-30.0 kg/m^2 eligibility",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "single oral atovaquone/proguanil 250/100 mg (100 mg proguanil hydrochloride)",
    regions = "Korea (Seoul National University Hospital)",
    genotype = paste(
      "All CYP2C19 normal metabolizers (*1/*1). SLC22A1 1022C>T",
      "(rs2282143): 9 CC and 4 CT of the 13 subjects genotyped; the 3",
      "ungenotyped subjects were assigned by a mixture model with fixed",
      "Hardy-Weinberg probabilities 0.7 (CC) / 0.3 (CT) and were all",
      "classified CC. The five other assayed SNPs (rs202220802,",
      "rs12208357, rs34059508, rs55918055, rs4646277) were homozygous",
      "wild type in all 13 subjects."
    ),
    n_observations = paste(
      "160 proguanil and 159 cycloguanil plasma concentrations over 48 h"
    ),
    notes = paste(
      "Post hoc analysis of the PK data from the CYP2C19-mediated",
      "tegoprazan/proguanil drug-drug-interaction study NCT04568772. Only",
      "the proguanil-alone and proguanil + tegoprazan 50 mg arms were used;",
      "esomeprazole and vonoprazan arms were excluded because of observed",
      "DDIs. Tegoprazan co-administration is not a covariate in the final",
      "model."
    )
  )
  ini({
    # ==================================================================
    # PROGUANIL (PARENT) STRUCTURAL PARAMETERS
    # ==================================================================
    # All typical values are from Khwarg 2024 Table 3 ("Parameter
    # estimates of the final population pharmacokinetic model of proguanil
    # and cycloguanil"), Final model / Estimate column. The second column
    # of that table is %RSE, not an IIV. Clearances are in L/h, volumes in
    # L, rate constants in 1/h.

    lka <- log(0.101)
    label("Absorption rate constant, depot -> liver (KA, 1/h)")                      # Khwarg 2024 Table 3: K_A = 0.101 1/h (RSE 6.0%)

    lcl <- log(33.5)
    label("Apparent clearance of proguanil via the non-cycloguanil pathway (CL/F, L/h)")  # Khwarg 2024 Table 3: CL/F = 33.5 L/h (RSE 16.4%). Figure 2 draws CL/F leaving the central compartment; the table abbreviation list defines it as "apparent clearance of proguanil via non-cycloguanil metabolic pathway".

    lvc <- log(58.4)
    label("Apparent central volume of distribution of proguanil (Vc/F, L)")          # Khwarg 2024 Table 3: V_c/F = 58.4 L (RSE 24.1%). Also used as the cycloguanil central volume -- Results 'Structural model': "The central volume of distribution of cycloguanil was assumed to be the same as the central volume of distribution of proguanil, due to identifiability issues."

    lq <- log(41.4)
    label("Apparent inter-compartmental clearance of proguanil (Q/F, L/h)")          # Khwarg 2024 Table 3: Q/F = 41.4 L/h (RSE 7.2%)

    lvp <- log(1400)
    label("Apparent peripheral volume of distribution of proguanil (Vp/F, L)")       # Khwarg 2024 Table 3: V_p/F = 1400 L (RSE 26.0%)

    lclint <- log(78.3)
    label("Hepatic intrinsic clearance of proguanil (CLint, L/h)")                   # Khwarg 2024 Table 3: CL_int = 78.3 L/h (RSE 18.6%). Feeds the well-stirred extraction ratio of Equation 8.

    # ==================================================================
    # SLC22A1 (OCT1) 1022C>T EFFECT ON OCT1-MEDIATED HEPATIC UPTAKE
    # ==================================================================
    # Khwarg 2024 Equation 8 multiplies the fub * CLint product inside the
    # well-stirred extraction ratio by FUP^GENO:
    #     EH = fub * CLint * FUP^GENO / (QBH + fub * CLint * FUP^GENO)
    # With the binary GENO indicator this is a plain multiplier of 1 for
    # the CC wild type and of FUP for the CT heterozygote, but it is coded
    # as the printed exponent form so the equation matches the source.
    # Results: "The hepatocyte uptake of proguanil was estimated to be
    # 0.42-fold lower in the CT genotype, compared to the CC genotype."
    # (0.416 rounds to the 0.42 quoted in the Abstract and Results.)

    e_snp_oct1_clint <- 0.416
    label("Relative fraction of OCT1-mediated hepatocyte uptake, SLC22A1 1022C>T CT vs CC (FUP, unitless)")  # Khwarg 2024 Table 3: FUP = 0.416 (RSE 54.1%)

    # ==================================================================
    # CYCLOGUANIL (METABOLITE) STRUCTURAL PARAMETERS
    # ==================================================================

    lktr_cycloguanil <- log(1.03)
    label("Cycloguanil efflux transit rate constant (mKT, 1/h)")                     # Khwarg 2024 Table 3: mK_T = 1.03 1/h (RSE 5.8%). Shared by all three transfers of the liver -> transit1 -> transit2 -> central chain (Figure 2 labels each arrow mK_T).

    lcl_cycloguanil <- log(97.0)
    label("Apparent clearance of cycloguanil (CLM/F, L/h)")                          # Khwarg 2024 Table 3: CL_M/F = 97.0 L/h (RSE 13.8%)

    # ==================================================================
    # FIXED PHYSIOLOGIC CONSTANTS OF THE WELL-STIRRED LIVER SUB-MODEL
    # ==================================================================
    # None of these were estimated. Khwarg 2024 Methods 'Structural model':
    # "Fixed values of fraction unbound in plasma (fu) of 0.25 and a
    # blood-to-plasma concentration ratio (CB/CP) of 2.7 was used to
    # calculate fub (Equation 3)" (citing refs 7 and 24); "The liver volume
    # (VH) of 1 L and hepatic blood flow (QBH) of 90 L/h of a typical 70 kg
    # adult was assumed and allometrically scaled based on body weight
    # (Equations 5 and 6). The hepatic plasma flow (QH) was calculated from
    # QBH using hematocrit of 0.45 (Equation 7)."

    fu_fix <- fixed(0.25)
    label("Fraction unbound in plasma (fu, unitless); literature value")             # Khwarg 2024 Methods 'Structural model', from refs 7 and 24

    bpr_fix <- fixed(2.7)
    label("Blood-to-plasma concentration ratio (CB/CP, unitless); literature value") # Khwarg 2024 Methods 'Structural model', from refs 7 and 24

    hct_fix <- fixed(0.45)
    label("Hematocrit used to convert hepatic blood flow to plasma flow (unitless); assumed") # Khwarg 2024 Methods 'Structural model' / Equation 7

    lvh_70 <- fixed(log(1))
    label("Liver volume of a typical 70 kg adult (VH, L); assumed")                  # Khwarg 2024 Methods 'Structural model' / Equation 5: V_H = 1 L * (body weight / 70 kg)

    lqbh_70 <- fixed(log(90))
    label("Hepatic blood flow of a typical 70 kg adult (QBH, L/h); assumed")         # Khwarg 2024 Methods 'Structural model' / Equation 6: Q_BH = 90 L/h * (body weight / 70 kg)^0.75

    e_wt_vh <- fixed(1)
    label("Allometric exponent on body weight for liver volume (unitless)")          # Khwarg 2024 Equation 5 prints (body weight / 70 kg) with no exponent, i.e. exponent 1

    e_wt_qbh <- fixed(0.75)
    label("Allometric exponent on body weight for hepatic blood flow (unitless)")    # Khwarg 2024 Equation 6 prints (body weight / 70 kg)^0.75

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY
    # ==================================================================
    # Khwarg 2024 Results 'Structural model': "The IIV of proguanil
    # parameters was incorporated in the apparent clearance of proguanil
    # via non-cycloguanil metabolic pathway (CL/F), central Vd, and CLint.
    # ... The IIV of cycloguanil was applied on the transit rate constant
    # (mKT)." Methods 'Structural model': "The inter-individual
    # variability (IIV) of PK parameters was evaluated using the
    # exponential error model on the assumption of being log-normally
    # distributed."
    #
    # Table 3 reports IIV as %CV, and its footnote gives the exact
    # back-transformation used: "CV, coefficient of variation calculated
    # as sqrt(e^omega^2 - 1) * 100". Inverting, omega^2 = log(CV^2 + 1):
    #   CL/F   31.7% -> log(0.317^2 + 1) = 0.095755   (shrinkage 5.9%)
    #   Vc/F   87.2% -> log(0.872^2 + 1) = 0.565532   (shrinkage 12.2%)
    #   CLint  53.7% -> log(0.537^2 + 1) = 0.253377   (shrinkage 4.2%)
    #   mKT    15.5% -> log(0.155^2 + 1) = 0.023741   (shrinkage 18.8%)
    # Each round-trips back to the printed %CV to one decimal place.
    # Table 3 reports no off-diagonal covariances, so the omega matrix is
    # diagonal.

    etalcl ~ 0.095755
    etalvc ~ 0.565532
    etalclint ~ 0.253377
    etalktr_cycloguanil ~ 0.023741

    # ==================================================================
    # RESIDUAL UNEXPLAINED VARIABILITY
    # ==================================================================
    # Khwarg 2024 Results 'Structural model': "Residual variability for
    # proguanil was modeled using a proportional error model. ... Residual
    # variability for cycloguanil was modeled using a combined error
    # model."
    #
    # Table 3 reports these on the STANDARD-DEVIATION scale, not as
    # variances. The discriminator is the unit label on the additive term:
    # Table 3 names the row "Additive error (ug/L)" -- a variance would
    # carry (ug/L)^2. Read as SDs the three values are also the ordinary
    # magnitudes for a Phase-1 LC-MS/MS assay (16% and 12.4% proportional);
    # read as variances they would imply 40% and 35% proportional error,
    # which is implausible for this data. See the vignette Assumptions and
    # deviations section.

    propSd <- 0.16
    label("Proportional residual SD for proguanil (unitless)")                       # Khwarg 2024 Table 3, Proguanil: Proportional error = 0.16 (RSE 9.1%)

    propSd_cycloguanil <- 0.124
    label("Proportional residual SD for cycloguanil (unitless)")                     # Khwarg 2024 Table 3, Cycloguanil: Proportional error = 0.124 (RSE 12.7%)

    addSd_cycloguanil <- 0.372
    label("Additive residual SD for cycloguanil (ug/L)")                             # Khwarg 2024 Table 3, Cycloguanil: Additive error (ug/L) = 0.372 (RSE 18.5%)
  })
  model({
    # ------------------------------------------------------------------
    # Unit-conversion and stoichiometric constants
    # ------------------------------------------------------------------
    # Doses are in mg and volumes in L, so amount / volume is mg/L; the
    # observed concentrations (and the additive residual SD) are in ug/L,
    # hence the factor of 1000 on both observations.
    mgPerL_to_ugPerL <- 1000

    # Molar-mass ratio converting a proguanil-mass metabolic flux into the
    # equivalent cycloguanil mass at 1:1 stoichiometry. Khwarg 2024
    # Results 'Structural model': "The ratio of molecular weight between
    # proguanil hydrochloride (290.19 g/mol) and cycloguanil
    # hydrochloride (288.17 g/mol) was used for molecular weight
    # conversion of the proguainil to cycloguanil."
    # Check: 288.17 / 290.19 = 0.993004.
    mwr_cycloguanil <- 288.17 / 290.19

    # ------------------------------------------------------------------
    # Individual PK parameters
    # ------------------------------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    clint <- exp(lclint + etalclint)
    ktr_cycloguanil <- exp(lktr_cycloguanil + etalktr_cycloguanil)
    cl_cycloguanil <- exp(lcl_cycloguanil)

    # ------------------------------------------------------------------
    # Fixed physiologic sub-model, Khwarg 2024 Equations 3, 5, 6, 7
    # ------------------------------------------------------------------
    #   Eq 3:  fub = fu / (CB/CP)              = 0.25 / 2.7 = 0.092593
    #   Eq 5:  VH  = 1 L    * (WT / 70 kg)
    #   Eq 6:  QBH = 90 L/h * (WT / 70 kg)^0.75
    #   Eq 7:  QH  = QBH * (1 - Hematocrit)
    fub <- fu_fix / bpr_fix
    vh <- exp(lvh_70) * (WT / 70)^e_wt_vh
    qbh <- exp(lqbh_70) * (WT / 70)^e_wt_qbh
    qh <- qbh * (1 - hct_fix)

    # ------------------------------------------------------------------
    # Well-stirred hepatic extraction with the OCT1 genotype effect
    # ------------------------------------------------------------------
    #   Eq 8:  EH  = fub * CLint * FUP^GENO / (QBH + fub * CLint * FUP^GENO)
    #   Eq 2:  CLBH = EH * QBH                       (hepatic BLOOD clearance)
    #   Eq 4:  CLH  = CLBH * (CB/CP)                 (hepatic PLASMA clearance)
    # Equation 1 is the GENO = 0 special case of Equation 8 and is not
    # coded separately. EH is built on the blood-side flow QBH and blood
    # unbound fraction fub; CLH is then referenced back to plasma so that
    # it can multiply a plasma concentration in the ODEs.
    uptake_clint <- fub * clint * e_snp_oct1_clint^SNP_SLC22A1_RS2282143
    eh <- uptake_clint / (qbh + uptake_clint)
    cl_h <- eh * qbh * bpr_fix

    # ------------------------------------------------------------------
    # ODE system, transcribed from Khwarg 2024 Figure 2
    # ------------------------------------------------------------------
    # Figure 2 arrow by arrow:
    #   Dose            -> Proguanil Depot
    #   Proguanil Depot -> Proguanil Liver            at KA
    #   Proguanil Liver -> Cycloguanil Liver          at CLH
    #   Proguanil Central -> Proguanil Liver          at QH
    #   Proguanil Liver -> Proguanil Central          at QH * (1 - EH)
    #   Proguanil Central <-> Proguanil Peripheral    at Q/F
    #   Proguanil Central -> eliminated               at CL/F
    #   Cycloguanil Liver -> transit1 -> transit2 -> Cycloguanil Central,
    #                                                 each at mKT
    #   Cycloguanil Central -> eliminated             at CLM/F
    #
    # The dose enters the liver rather than the central compartment, so
    # hepatic first-pass extraction is structural: the fraction of an
    # absorbed dose surviving the liver is QH * (1 - EH) / (QH * (1 - EH)
    # + CLH), and the complement is converted to cycloguanil before ever
    # reaching systemic plasma.
    #
    # Independent check of this topology against the paper's own reported
    # exposures, at the cohort mean 76.12 kg (steady-state mass balance,
    # no IIV): CC genotype gives AUCinf 1526 h*ug/L for proguanil and 500
    # h*ug/L for cycloguanil, against the Table 2 AUClast values of 1205
    # and 445 plus a log-linear tail extrapolation from the 48 h Figure 1
    # profiles of roughly 1434 and 477 -- both within ~6%. The CT/CC
    # cycloguanil exposure ratio comes out at 0.58, matching the Abstract's
    # "0.6-fold lower exposure of cycloguanil".
    d/dt(depot) <- -ka * depot
    d/dt(liver) <-
      ka * depot +
      qh * (central / vc) -
      qh * (1 - eh) * (liver / vh) -
      cl_h * (liver / vh)
    d/dt(central) <-
      qh * (1 - eh) * (liver / vh) -
      qh * (central / vc) -
      cl * (central / vc) -
      q * (central / vc) +
      q * (peripheral1 / vp)
    d/dt(peripheral1) <-
      q * (central / vc) -
      q * (peripheral1 / vp)

    d/dt(liver_cycloguanil) <-
      mwr_cycloguanil * cl_h * (liver / vh) -
      ktr_cycloguanil * liver_cycloguanil
    d/dt(transit1_cycloguanil) <-
      ktr_cycloguanil * liver_cycloguanil -
      ktr_cycloguanil * transit1_cycloguanil
    d/dt(transit2_cycloguanil) <-
      ktr_cycloguanil * transit1_cycloguanil -
      ktr_cycloguanil * transit2_cycloguanil
    d/dt(central_cycloguanil) <-
      ktr_cycloguanil * transit2_cycloguanil -
      cl_cycloguanil * (central_cycloguanil / vc)

    # ------------------------------------------------------------------
    # Observations. Both analytes are sampled from systemic plasma, and
    # cycloguanil shares the proguanil central volume (Results
    # 'Structural model').
    # ------------------------------------------------------------------
    Cc <- mgPerL_to_ugPerL * central / vc
    Cc_cycloguanil <- mgPerL_to_ugPerL * central_cycloguanil / vc

    Cc ~ prop(propSd)
    Cc_cycloguanil ~ add(addSd_cycloguanil) + prop(propSd_cycloguanil)
  })
}

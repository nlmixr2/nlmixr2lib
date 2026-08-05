MohammedAli_2023_tacrolimus <- function() {
  description <- paste(
    "Two-compartment population PK model for once-daily extended-release",
    "LCP-Tac tacrolimus (MeltDose technology, Envarsus) in stable adult renal",
    "transplant recipients (Mohammed Ali 2023), parameterized in apparent",
    "elimination clearance CL/F, apparent distributional clearance CLD/F and",
    "apparent central and peripheral volumes Vc/F and Vp/F. Delayed absorption",
    "is described by a Savic 2007 transit-compartment chain with the number of",
    "transit compartments fixed at NN = 2 (mean transit time MTT = 2.91 h)",
    "feeding a first-order absorption step ka into the central compartment. A",
    "combined CYP3A4/CYP3A5 cluster phenotype (high, intermediate and poor",
    "metabolizer, reconstructed inside model() from the recipient CYP3A5",
    "expresser status and the CYP3A4*22 rs35599367 carrier indicator) gives",
    "three distinct typical CL/F values. Inter-individual variability is",
    "carried on CL/F, Vc/F, Vp/F and MTT, with inter-occasion variability on",
    "CL/F and a proportional residual error on whole-blood concentrations.",
    sep = " "
  )
  reference <- paste(
    "Mohammed Ali Z, Meertens M, Fernandez B, Fontova P, Vidal-Alabro A,",
    "Rigo-Bonnin R, Melilli E, Cruzado JM, Grinyo JM, Colom H, Lloberas N.",
    "CYP3A5*3 and CYP3A4*22 Cluster Polymorphism Effects on LCP-Tac Tacrolimus",
    "Exposure: Population Pharmacokinetic Approach.",
    "Pharmaceutics. 2023;15(12):2699. doi:10.3390/pharmaceutics15122699.",
    sep = " "
  )
  vignette <- "MohammedAli_2023_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = FALSE),
    transit1    = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = FALSE),
    transit2    = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = FALSE),
    transit3    = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "tacrolimus", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CYP3A5_EXPR = list(
      description        = "Recipient CYP3A5 expresser status (rs776746)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes              = paste(
        "1 = at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3);",
        "0 = CYP3A5*3/*3. Genotyped by TaqMan assay (Mohammed Ali 2023 Methods",
        "section 2.3). Not used on its own in the final model: it is one of the",
        "two inputs from which the three-level CYP3A4/CYP3A5 cluster phenotype",
        "(HM / IM / PM) is reconstructed inside model(), following the",
        "two-binary-input convention documented in the CYP3A5_EXPR_DONOR entry",
        "of inst/references/covariate-columns.md. Cohort distribution",
        "(Table 1, all patients): *1/*1 4 (4.1%), *1/*3 17 (17.3%),",
        "*3/*3 77 (78.6%), so CYP3A5_EXPR = 1 in 21 of 98 subjects.",
        sep = " "
      ),
      source_name        = "CYP3A5 genotype"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description        = "CYP3A4*22 (rs35599367, c.522-191C>T) carrier indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A4*1/*1 non-carrier)",
      notes              = paste(
        "1 = carrier of at least one CYP3A4*22 (T) allele (genotype *1/*22;",
        "no *22/*22 homozygotes were observed); 0 = CYP3A4*1/*1. Genotyped by",
        "TaqMan assay (Mohammed Ali 2023 Methods section 2.3). Second input to",
        "the cluster-phenotype reconstruction inside model(). Cohort",
        "distribution (Table 1, all patients): *1/*1 86 (86.7%),",
        "*1/*22 12 (13.3%).",
        sep = " "
      ),
      source_name        = "CYP3A4 genotype"
    ),
    OCC = list(
      description        = "Sampling-occasion index for inter-occasion variability on CL/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Integer occasion index 1..5 decomposed inside model() into",
        "mutually-exclusive binary indicators multiplexing per-occasion etas",
        "on log CL/F. Mohammed Ali 2023 reports a single IOV magnitude",
        "(Table 3 'IOV (%) 44.8') without stating the number of occasions",
        "used; the upper bound of 5 is taken from Methods section 2.1,",
        "'from one to five C0 samples per patient were obtained'. Occasions",
        "2-5 have their variance fixed equal to occasion 1 (the NONMEM",
        "$OMEGA BLOCK(1) SAME pattern; nlmixr2 has no SAME shortcut).",
        "Records with OCC outside 1..5 carry no IOV.",
        sep = " "
      ),
      source_name        = "occasion"
    )
  )

  # Covariates screened by Mohammed Ali 2023 but NOT retained in the final
  # model. Documented for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested univariately on CL/F, CLD/F, Vc/F and Vp/F; no significant MOFV",
        "drop (p > 0.05), and allometric inclusion of body weight on these",
        "parameters worsened the model (Results section 3.2). Cohort:",
        "74.73 kg mean (IQR 65-81.13).",
        sep = " "
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Tested univariately on the disposition parameters; no significant",
        "MOFV drop (Results section 3.2). Cohort: 26.33 kg/m^2 mean",
        "(IQR 22.94-28.94).",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Tested univariately; no significant MOFV drop (Results section 3.2).",
        "Cohort: 56 years mean (IQR 46-68).",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Listed among the covariates investigated (Methods 'Covariate Model');",
        "not retained in the final model. Cohort: 30 of 98 female (30.6%).",
        sep = " "
      )
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = paste(
        "Tested as a covariate and also as a standardizing factor on the",
        "residual error (concentrations standardized to a hematocrit of 45%);",
        "neither improved the model (Results section 3.2). The Discussion",
        "attributes the absence of a hematocrit effect to the near-normal",
        "hematocrit of these stable recipients (IQR 37.4-44.0%) compared with",
        "'de novo' cohorts.",
        sep = " "
      )
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 c.3435C>T (rs1045642) variant carrier indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Genotyped and tested on CL/F; the ABCB1 *C/*C SNP failed to influence",
        "the model significantly (Results section 3.2) and no differences in",
        "dose-normalized C0 or AUC were found between *C carriers and",
        "non-carriers (Table 2). Cohort: *T/*T 21%, *C/*T 47%, *C/*C 32%.",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 98,
    n_studies      = 1,
    n_observations = 655,
    age_range      = "56 years mean (IQR 46-68)",
    weight_range   = "74.73 kg mean (IQR 65-81.13)",
    sex_female_pct = 30.6,
    disease_state  = paste(
      "stable adult renal transplant recipients at least 6 months",
      "post-transplant, converted from twice-daily IR-Tac to once-daily",
      "LCP-Tac at a 0.7 dose-conversion ratio, on triple immunosuppression",
      "(tacrolimus + mycophenolate mofetil + prednisone)",
      sep = " "
    ),
    dose_range     = "0.5-12 mg once daily (median total daily dose 2 mg)",
    regions        = "Spain (single centre, Hospital Universitari de Bellvitge, Barcelona)",
    renal_function = "eGFR (CKD-EPI) 47.72 mL/min mean (IQR 36-58); serum creatinine 146.9 umol/L mean (IQR 116-163)",
    hematocrit     = "40.53% mean (IQR 37.4-44.0)",
    genotypes      = paste(
      "CYP3A5 *1/*1 4 (4.1%), *1/*3 17 (17.3%), *3/*3 77 (78.6%);",
      "CYP3A4 *1/*1 86 (86.7%), *1/*22 12 (13.3%);",
      "CYP3A4/A5 cluster HM 19 (19.4%), IM 68 (69.4%), PM 11 (11.2%);",
      "ABCB1 *T/*T 21 (21%), *C/*T 46 (47%), *C/*C 31 (32%)",
      sep = " "
    ),
    notes          = paste(
      "Mohammed Ali 2023 Table 1. Pooled dataset: 30 patients from the",
      "open-label clinical trial NCT02961608 contributed rich steady-state",
      "profiles (480 of 655 observations, 10-18 samples over 24 h), and 68",
      "patients from routine hospital follow-up contributed 1-5 trough (C0)",
      "samples each (175 observations). Whole-blood tacrolimus measured by",
      "LC-MS/MS with a 1.0 ng/mL limit of quantitation.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Apparent elimination clearance CL/F -- three typical values, one per
    # CYP3A4/CYP3A5 cluster phenotype (Mohammed Ali 2023 Table 3). Each was
    # estimated as its own THETA (each carries its own RSE), so all three
    # are encoded directly rather than as ratios to a reference group.
    # =====================================================================
    lcl_hm <- log(19.6)
    label("Apparent elimination clearance CL/F in high metabolizers (L/h)")
    # Table 3 row 'CL/F HM (L.h-1) = 19.6 (RSE 10%)'

    lcl_im <- log(10.6)
    label("Apparent elimination clearance CL/F in intermediate metabolizers (L/h)")
    # Table 3 row 'CL/F IM (L.h-1) = 10.6 (RSE 5.2%)'

    lcl_pm <- log(7.37)
    label("Apparent elimination clearance CL/F in poor metabolizers (L/h)")
    # Table 3 row 'CL/F PM (L.h-1) = 7.37 (RSE 11.9%)'

    # =====================================================================
    # Disposition (Mohammed Ali 2023 Table 3)
    # =====================================================================
    lvc <- log(169)
    label("Apparent central volume of distribution Vc/F (L)")
    # Table 3 row 'Vc/F (L) = 169 (RSE 17.2%)'

    lq <- log(37.6)
    label("Apparent distributional clearance CLD/F (L/h)")
    # Table 3 row 'CL D /F (L.h-1) = 37.6 (RSE 13.5%)'

    lvp <- log(460)
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Table 3 row 'Vp (L) = 460 (RSE 27.8%)'; Vc/F + Vp/F = 629 L matches the Discussion

    # =====================================================================
    # Delayed absorption: Savic 2007 transit chain (reference [39]) with the
    # number of transit compartments fixed at 2, feeding a first-order
    # absorption step ka into the central compartment.
    # =====================================================================
    lmtt <- log(2.91)
    label("Mean absorption transit time MTT through the Savic chain (h)")
    # Table 3 row 'MTT (h) = 2.91 (RSE 15.5%)'; Results section 3.2 'mean absorption transit time of 2 h 55 min'

    nn_fix <- fixed(2)
    label("Number of Savic transit compartments NN (integer, unitless)")
    # Table 3 row 'NN = 2 FIX'; Results section 3.2 'the number of absorption compartments were fixed to 2'

    lka <- log(0.72)
    label("First-order absorption rate constant ka out of the transit chain (1/h)")
    # Table 3 row 'Ka (h-1) = 0.72 (RSE 33.2%)'

    # =====================================================================
    # Inter-individual variability. Table 3 reports IIV as a CV percentage;
    # converted to the log-scale variance with omega^2 = log(CV^2 + 1).
    # =====================================================================
    etalcl ~ log(1 + 0.379^2)
    # Table 3 'CL/F HM ... IIV % 37.9 (RSE 17.9%)' -> omega^2 = log(1 + 0.379^2) = 0.1342

    etalvc ~ log(1 + 0.70^2)
    # Table 3 'Vc/F (L) ... IIV % 70 (RSE 41.4%)' -> omega^2 = log(1 + 0.70^2) = 0.3988

    etalvp ~ log(1 + 0.75^2)
    # Table 3 'Vp (L) ... IIV % 75 (RSE 44.3%)' -> omega^2 = log(1 + 0.75^2) = 0.4463

    etalmtt ~ log(1 + 0.546^2)
    # Table 3 'MTT (h) ... IIV % 54.6 (RSE 37.1%)' -> omega^2 = log(1 + 0.546^2) = 0.2609

    # =====================================================================
    # Inter-occasion variability on CL/F. A single IOV magnitude is reported
    # (Table 3 'IOV (%) 44.8'), so occasions 2-5 are fixed equal to occasion
    # 1 (NONMEM $OMEGA BLOCK(1) SAME).
    # =====================================================================
    etaiov_cl_1 ~ log(1 + 0.448^2)
    # Table 3 row 'IOV (%) = 44.8 (RSE 27.4%)' -> omega^2 = log(1 + 0.448^2) = 0.1829
    etaiov_cl_2 ~ fixed(log(1 + 0.448^2))
    etaiov_cl_3 ~ fixed(log(1 + 0.448^2))
    etaiov_cl_4 ~ fixed(log(1 + 0.448^2))
    etaiov_cl_5 ~ fixed(log(1 + 0.448^2))

    # =====================================================================
    # Residual error -- proportional (Results section 3.2 'A proportional
    # error model best described the RE distribution').
    # =====================================================================
    propSd <- 0.0967
    label("Proportional residual error (fraction)")
    # Table 3 row 'RE (%) = 9.67 (RSE 8%)'
  })

  model({
    # ------------------------------------------------------------------
    # 1. CYP3A4/CYP3A5 cluster phenotype, reconstructed from the two
    #    genotype inputs exactly as defined in Methods section 2.3:
    #      PM = CYP3A4*22 carriers + CYP3A5*3/*3
    #      IM = CYP3A4*22 non-carriers + CYP3A5*3/*3,
    #           OR CYP3A4*22 carriers + CYP3A5*1 carriers
    #      HM = CYP3A4*22 non-carriers + CYP3A5*1 carriers
    #    The three indicators are mutually exclusive and sum to 1.
    # ------------------------------------------------------------------
    is_hm <- (1 - SNP_CYP3A4_RS35599367) * CYP3A5_EXPR
    is_pm <- SNP_CYP3A4_RS35599367 * (1 - CYP3A5_EXPR)
    is_im <- 1 - is_hm - is_pm

    # ------------------------------------------------------------------
    # 2. Occasion indicators multiplexing the per-occasion CL/F etas.
    # ------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    # ------------------------------------------------------------------
    # 3. Individual PK parameters. The cluster phenotype selects one of the
    #    three estimated typical log CL/F values; IIV and IOV are additive
    #    on the log scale.
    # ------------------------------------------------------------------
    lcl_typ <- lcl_hm * is_hm + lcl_im * is_im + lcl_pm * is_pm
    cl  <- exp(lcl_typ + etalcl + iov_cl)
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq)
    ka  <- exp(lka)
    mtt <- exp(lmtt + etalmtt)
    nn  <- nn_fix

    # ------------------------------------------------------------------
    # 4. Micro-constants. Savic 2007 parameterizes the transit chain by the
    #    mean transit time: ktr = (NN + 1) / MTT, i.e. NN + 1 = 3 first-order
    #    transfers at rate ktr carry the dose from the dosing compartment to
    #    the absorption compartment, so the mean time to reach that
    #    compartment is exactly MTT. The absorption compartment then empties
    #    into central at the separately estimated rate ka.
    # ------------------------------------------------------------------
    ktr <- (nn + 1) / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 5. ODE system. Chain: depot -> transit1 -> transit2 -> transit3, all
    #    at ktr (NN + 1 = 3 transfers), then transit3 -> central at ka.
    #    transit3 is the Savic absorption compartment; the two Savic transit
    #    compartments proper are transit1 and transit2 (NN = 2).
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ka  * transit3
    d/dt(central)     <-  ka  * transit3 - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # ------------------------------------------------------------------
    # 6. Observation. Doses are in mg and volumes in L, so central / vc is
    #    mg/L; the 1000 factor converts to the ng/mL units in which whole-
    #    blood tacrolimus concentrations are reported.
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}

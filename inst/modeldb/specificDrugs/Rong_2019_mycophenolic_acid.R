Rong_2019_mycophenolic_acid <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) in stable adult kidney transplant recipients co-treated with immediate-release tacrolimus and free of corticosteroids (Rong 2019). Apparent clearance CL/F carries two metabolite-exposure covariates: a power effect of the acyl mycophenolic acid glucuronide trough concentration (ACMPAG_CC, exponent -0.09) and a power effect of the dose-normalised MPAG-to-MPA AUC ratio (formed inside model() as AUC_MPAG / AUC_MPA, exponent 0.68). Inter-individual variability is log-normal on all six structural parameters. Transcribed from the executable mrgsolve source distributed as S1 File of Maizaud 2025, which re-implements the Rong 2019 model; the Rong 2019 primary publication was not available at extraction time (see vignette Errata)."
  reference <- paste(
    "Rong Y, Mayo P, Ensom MHH, Kiang TKL. Population pharmacokinetics of",
    "mycophenolic acid co-administered with tacrolimus in corticosteroid-free",
    "adult kidney transplant patients. Clin Pharmacokinet. 2019;58(11):1483-1495.",
    "doi:10.1007/s40262-019-00771-3.",
    "Parameter values transcribed from the mrgsolve implementation in S1 File of",
    "Maizaud F, Arraki-Zava S, Sayadi H, Fromage Y, Marquet P, Woillard J-B,",
    "Monchaud C. Population pharmacokinetic modeling of missed mycophenolate",
    "mofetil doses: impact on exposure and exploration of mitigation strategies.",
    "PLoS One. 2025;20(8):e0330854. doi:10.1371/journal.pone.0330854.",
    sep = " "
  )
  vignette <- "Maizaud_2025_missed_mycophenolate_doses"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "mycophenolic acid (MPA)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "mycophenolic acid (MPA)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    ACMPAG_CC = list(
      description = paste(
        "Acyl mycophenolic acid glucuronide (AcMPAG) plasma concentration used as a",
        "covariate on MPA apparent clearance. AcMPAG is the pharmacologically active",
        "acyl-linked glucuronide conjugate of MPA formed by UGT2B7, distinct from the",
        "inactive 7-O-glucuronide MPAG formed by UGT1A9."
      ),
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as an uncentered power term ACMPAG_CC^(-0.09), so the typical",
        "value lka/lcl parameter lcl = log(2.87) is an INTERCEPT rather than a typical",
        "clearance: with the source default covariates the realised CL/F is",
        "2.87 * 0.54^(-0.09) * (588.8/53)^0.68 = 15.60 L/h, which reproduces the",
        "steady-state AUC0-12h values of Maizaud 2025 Table 2 to within 1-3% across",
        "all three Rong dose arms. The default value 0.54 mg/L is the mrgsolve",
        "'@covariates' default in Maizaud 2025 S1 File; its provenance (presumably a",
        "Rong 2019 cohort median trough) is not stated in any on-disk source. The S1",
        "annotation text reads 'Mycophenolic acid glucuronide C0' while both the",
        "parameter name (AcMPAG) and Maizaud 2025 Table 1 identify the ACYL",
        "glucuronide; the parameter name and Table 1 are followed here. See vignette",
        "Errata."
      ),
      source_name = "AcMPAG"
    ),
    AUC_MPAG = list(
      description = paste(
        "Dose-normalised area under the plasma concentration-time curve of",
        "mycophenolic acid 7-O-glucuronide (MPAG), the inactive phenolic glucuronide",
        "metabolite of MPA, over the dosing interval."
      ),
      units = "mg*h/L per g MMF",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Rong 2019 prints a single covariate on CL/F -- the MPAG:MPA AUC RATIO -- with",
        "power exponent 0.68. Per the covariate-register decision for this extraction,",
        "the ratio is NOT registered as its own canonical: the two component AUCs are",
        "registered as members of the existing AUC_<DRUG> family and the ratio is",
        "formed inside model() as (AUC_MPAG / AUC_MPA). Both AUCs are normalised by",
        "the same MMF dose, so the normalisation cancels in the ratio and the",
        "internally-computed value reproduces the paper's covariate exactly",
        "(588.8 / 53 = 11.11). Default 588.8 mg*h/L/g is the mrgsolve '@covariates'",
        "default in Maizaud 2025 S1 File; provenance not stated on disk."
      ),
      source_name = "AUCMPAG"
    ),
    AUC_MPA = list(
      description = paste(
        "Dose-normalised area under the plasma concentration-time curve of",
        "mycophenolic acid (MPA) over the dosing interval; the denominator of the",
        "MPAG:MPA AUC ratio covariate on CL/F."
      ),
      units = "mg*h/L per g MMF",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Denominator of the dose-normalised MPAG:MPA AUC ratio covariate; see the",
        "AUC_MPAG notes for why the ratio is formed inside model() from the two",
        "family members rather than registered as a single ratio canonical. Default",
        "53 mg*h/L/g is the mrgsolve '@covariates' default in Maizaud 2025 S1 File;",
        "provenance not stated on disk. Note this is a COVARIATE describing the",
        "subject's own historical MPA exposure, not a model output -- it is a",
        "measured per-subject quantity in the Rong 2019 dataset."
      ),
      source_name = "AUCMPA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    disease_state = "stable adult kidney transplant recipients, corticosteroid-free",
    co_medication = "immediate-release tacrolimus as the concomitant calcineurin inhibitor",
    dose_range = "500-1000 mg MMF BID in the Maizaud 2025 simulations",
    notes = paste(
      "Population described in Maizaud 2025 Methods 2.1: 'The first model, published",
      "by Rong et al. [27], was developed from 27 stable adult kidney transplant",
      "recipients co-treated with immediate-release tacrolimus.' The",
      "corticosteroid-free status comes from the Rong 2019 title as reproduced in the",
      "Maizaud 2025 S1 File section heading. Detailed baseline demographics",
      "(age, weight, sex, race) are in the Rong 2019 primary publication, which was",
      "not available at extraction time."
    )
  )

  ini({
    # Structural parameters. Maizaud 2025 S1 File, Rong '[PARAM] @annotated' block;
    # cross-checked against Maizaud 2025 Table 1 column 'Rong et al. [27]'.
    ltlag <- log(0.162); label("Absorption lag time (h)")                       # S1 TVTlag = 0.162; Table 1 'Tlag=0.162 h'
    lka   <- log(1.98);  label("Absorption rate constant (1/h)")                # S1 TVKA = 1.98; Table 1 'Ka=1.98 h-1'
    lcl   <- log(2.87);  label("Apparent clearance intercept (CL/F, L/h)")      # S1 TVCL = 2.87; Table 1 'CL/F=2.87 L/h'
    lvc   <- log(25);    label("Apparent central volume (V1/F, L)")             # S1 TVVC = 25; Table 1 'V1/F=25 L'
    lq    <- log(36.7);  label("Apparent intercompartmental clearance (Q/F, L/h)") # S1 TVQ = 36.7; Table 1 'Q/F=36.7 L/h'
    lvp   <- log(607);   label("Apparent peripheral volume (V2/F, L)")          # S1 TVVP = 607; Table 1 'V2/F=607 L'

    # Covariate effects on CL/F. Both are uncentered power terms, so lcl above is an
    # intercept and not a typical clearance -- see the ACMPAG_CC covariateData notes.
    e_acmpag_cc_cl     <- -0.09; label("Power exponent of AcMPAG concentration on CL/F (unitless)")     # S1 'B1 : - 0.09 : covariate parameter estimate mycophenolic acid acyl- glucuronide'
    e_auc_mpag_mpa_cl  <-  0.68; label("Power exponent of the MPAG:MPA AUC ratio on CL/F (unitless)")   # S1 'B2 : 0.68 : covariate parameter estimate'

    # IIV. S1 '[OMEGA] @annotated' entries are VARIANCES; Maizaud 2025 Table 1
    # reports the same quantities as SDs (omega). All six cross-check: Table 1
    # omega^2 = 0.99^2 = 0.980, 0.18^2 = 0.032, 1.08^2 = 1.166, 0.23^2 = 0.053,
    # 0.27^2 = 0.073, 1.08^2 = 1.166 against S1 0.98 / 0.03 / 1.16 / 0.05 / 0.07 / 1.16.
    etaltlag ~ 1.16  # S1 'ETATlag : 1.16 : IIV on Tlag'; Table 1 'omega Tlag=1.08'
    etalka   ~ 0.98  # S1 'ETAKA : 0.98 : IIV on Ka absorbance'; Table 1 'omega Ka=0.99'
    etalcl   ~ 0.05  # S1 'ETACL : 0.05 : IIV on apparent clearance'; Table 1 'omega CL/F=0.23'
    etalvc   ~ 0.03  # S1 'ETAVC : 0.03 : IVV on Central compartment volume'; Table 1 'omega V1/F=0.18'
    etalq    ~ 0.07  # S1 'ETAQ : 0.07 : IVV on apparent intercompartmental clearance'; Table 1 'omega Q/F=0.27'
    etalvp   ~ 1.16  # S1 'ETAVP : 1.16 : IVV on apparent peripheral compartment volume'; Table 1 'omega V2/F=1.08'

    # Residual error. These are the ORIGINAL Rong 2019 values as reported in Maizaud
    # 2025 Table 1, not the deliberately-decreased sigma = 0.001 that Maizaud 2025
    # substituted for its own simulations (Methods 2.1: 'residual variability was
    # fixed at a minimal value (sigma = 0.001) to isolate the influence of
    # inter-individual variability and covariates on exposure').
    propSd <- 0.32; label("Proportional residual error (fraction)") # Table 1 'Original proportional error sigma =0.32'
    addSd  <- 0.08; label("Additive residual error (mg/L)")         # Table 1 'Original additive error sigma =0.08 mg/L'
  })

  model({
    # Individual parameters. The two CL/F covariates are UNCENTERED power terms, i.e.
    # exactly as written in the Maizaud 2025 S1 File '[MAIN]' block:
    #   double CL = TVCL*pow(AcMPAG,B1)*pow(AUCMPAG/AUCMPA,B2)*exp(ETACL);
    # Rong 2019 prints a single covariate for the second term -- the dose-normalised
    # MPAG:MPA AUC ratio -- which is formed here from the two registered AUC_<DRUG>
    # family members. Both AUCs are normalised by the same MMF dose so the
    # normalisation cancels and the ratio reproduces the paper's covariate exactly.
    auc_mpag_mpa_ratio <- AUC_MPAG / AUC_MPA

    cl   <- exp(lcl + etalcl) * ACMPAG_CC^e_acmpag_cc_cl * auc_mpag_mpa_ratio^e_auc_mpag_mpa_cl
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka + etalka)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp)
    tlag <- exp(ltlag + etaltlag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODE system. MMF is dosed into the depot; MPA appears in the central
    # compartment. The MMF-to-MPA molar conversion is absorbed into the apparent
    # parameters CL/F and V/F, so dose amounts stay in mg of MMF.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

Saqr_2025_mycophenolic_acid <- function() {
  description <- paste(
    "Semi-mechanistic population PK model for mycophenolic acid (MPA) and its",
    "glucuronide metabolites MPAG and AcMPAG after intravenous mycophenolate",
    "mofetil (MMF) in adult hematopoietic cell transplant recipients",
    "(Saqr 2025). One-compartment MPA and MPAG and two-compartment AcMPAG",
    "dispositions linked through fixed metabolic fractions (0.87 to MPAG,",
    "0.13 to AcMPAG). Enterohepatic circulation (EHC) runs MPAG -> bile, with",
    "25% draining continuously into the gut and 75% stored in a gallbladder",
    "that empties into the gut for 30 min starting at an estimated time after",
    "each dose; gut MPAG is deglucuronidated and returned to the MPA central",
    "compartment. The extent of EHC (FrEHC, logit scale) depends on the stool",
    "relative abundances of two bacterial communities (Bacteroides",
    "uniformis-dominant and Bacteroides vulgatus-dominant). Fixed allometric",
    "body-weight scaling (reference 70 kg). Doses are mg of MMF; states hold",
    "MMF-molar-equivalent mg and outputs are mg/L of MPA, MPAG and AcMPAG."
  )
  reference <- paste(
    "Saqr A, Al-Kofahi M, Staley C, Rashidi A, Carlson B, Holtan SG, Kaiser T,",
    "Jacobson PA. (2025). Microbiome-Informed Dosing: Exploring Gut Microbial",
    "Communities Impact on Mycophenolate Enterohepatic Circulation and",
    "Therapeutic Target Achievement. Clin Pharmacol Ther 118(6):1477-.",
    "doi:10.1002/cpt.3740",
    sep = " "
  )
  vignette <- "Saqr_2025_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The source analysed doses and concentrations in molar equivalents
  # (Supplementary File 3, Structural Model) and assumed MMF is converted
  # completely and instantaneously to MPA. Every state here holds a molar
  # amount expressed as mg of MMF (mg MMF / 433.49 g/mol = mmol), so a dose
  # in mg MMF enters `central` unchanged, the 0.87 / 0.13 metabolic split and
  # the MPAG -> MPA reconversion conserve moles, and each concentration is
  # converted to mg/L of its own analyte with its molecular weight in model().
  compartmentData <- list(
    central = list(
      analyte = "mycophenolic acid (MPA)",
      units = "mg (MMF molar equivalent)",
      specimen = "plasma",
      verified = TRUE
    ),
    central_mpag = list(
      analyte = "mycophenolic acid 7-O-glucuronide (MPAG)",
      units = "mg (MMF molar equivalent)",
      specimen = "plasma",
      verified = TRUE
    ),
    central_acmpag = list(
      analyte = "mycophenolic acid acyl-glucuronide (AcMPAG)",
      units = "mg (MMF molar equivalent)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_acmpag = list(
      analyte = "mycophenolic acid acyl-glucuronide (AcMPAG)",
      units = "mg (MMF molar equivalent)",
      specimen = "tissue",
      verified = TRUE
    ),
    gallbladder_mpag = list(
      analyte = "mycophenolic acid 7-O-glucuronide (MPAG)",
      units = "mg (MMF molar equivalent)",
      specimen = "bile",
      verified = TRUE
    ),
    gut_mpag = list(
      analyte = "mycophenolic acid 7-O-glucuronide (MPAG)",
      units = "mg (MMF molar equivalent)",
      specimen = "faeces",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight measured on the pharmacokinetic sampling day.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling, reference 70 kg (Supplementary File 3",
        "Eq.3): exponent 0.75 on clearances, 1 on volumes and -0.25 on rate",
        "constants. The source does not list which rate constants; both the",
        "estimated kicv (K51) and the fixed gallbladder emptying constant kehc",
        "(KGG) are scaled here."
      ),
      source_name = "WT"
    ),
    STOOL_COMMUNITY_BACTEROIDES_UNIFORMIS_PCT = list(
      description = paste(
        "Stool relative abundance (%) of bacterial community 1, the",
        "Bacteroides uniformis-dominant community: the summed shotgun",
        "metagenomic relative abundance of B. uniformis, B. thetaiotaomicron,",
        "B. faecis, Akkermansia muciniphila, Clostridium clostridioforme and",
        "Ruminococcus gnavus (Supplementary File 7 Figure S3)."
      ),
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters logit(FrEHC) as e * log(abundance / 19.97), where 19.97% is",
        "the cohort median (main text Eq.3; Supplementary Table S4; observed",
        "range 0-73.44%). The authors replaced a 0% abundance by one tenth of",
        "the smallest non-zero relative abundance in their analysis dataset",
        "before log-transformation (Supplementary File 5); that floor is not",
        "printed. Supply a small positive value for a zero abundance: a value",
        "of exactly 0 gives log(0) = -Inf and therefore FrEHC = 0 (no EHC)."
      ),
      source_name = "Community 1"
    ),
    STOOL_COMMUNITY_BACTEROIDES_VULGATUS_PCT = list(
      description = paste(
        "Stool relative abundance (%) of bacterial community 2, the",
        "Bacteroides vulgatus-dominant community: the summed shotgun",
        "metagenomic relative abundance of B. vulgatus, B. ovatus,",
        "B. massiliensis, B. stercoris and B. finegoldii (Supplementary File 7",
        "Figure S3)."
      ),
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters logit(FrEHC) as e * log(abundance / 13.23), where 13.23% is",
        "the cohort median (main text Eq.3; Supplementary Table S4; observed",
        "range 0-89.23%, at least 25% of subjects at 0%). Zero handling as",
        "for STOOL_COMMUNITY_BACTEROIDES_UNIFORMIS_PCT."
      ),
      source_name = "Community 2"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Screened by GAM on the EBEs; not significant (Results, Clinical covariates)."
    ),
    AGE = list(
      description = "Age at transplant.",
      units = "years",
      type = "continuous",
      notes = "Screened by GAM on the EBEs; not significant (Results, Clinical covariates)."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened by GAM on the EBEs; not significant (Results, Clinical covariates)."
    ),
    ALT = list(
      description = "Alanine transaminase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened by GAM on the EBEs; not significant (Results, Clinical covariates)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate.",
      units = "mL/min/1.73m^2",
      type = "continuous",
      notes = "Screened by GAM on the EBEs; not significant (Results, Clinical covariates)."
    ),
    N_ANTIBIOTICS = list(
      description = "Number of antibiotics given within 72 h of pharmacokinetic sampling (2, 3 or >= 4).",
      units = "(count)",
      type = "categorical",
      notes = "Screened by GAM on the EBEs; not significant (Results, Clinical covariates)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = "480 plasma concentrations (160 each of MPA, MPAG and AcMPAG; Results).",
    age_median = "46 years (IQR 36-55; Table 1)",
    age_range = "18-75 years (eligibility; Supplementary File 1)",
    weight_median = "87 kg (IQR 76-103; Table 1)",
    sex_female_pct = 40,
    race_ethnicity = c(White = 90, `American Indian` = 10),
    disease_state = paste(
      "Adult hematopoietic cell transplant recipients with hematologic",
      "malignancies (acute leukemia 75%, acute myeloid leukemia 5%, lymphoma",
      "10%, myelodysplastic syndrome 10%), studied on day 6-8 after",
      "transplant with post-transplant cyclophosphamide GVHD prophylaxis",
      "followed by MMF and tacrolimus."
    ),
    dose_range = "MMF 1000 mg IV every 8 h infused over 2 h, after at least 3 unchanged doses.",
    regions = "United States (University of Minnesota).",
    co_medication = paste(
      "Tacrolimus in all subjects; levofloxacin prophylaxis (85%), cefepime",
      "(55-70%), vancomycin (40%) and other antibiotics (Table 1)."
    ),
    sampling_window = "0 (pre-infusion), 2, 3, 4, 5, 6, 7 and 8 h after the start of the infusion (Supplementary File 1).",
    notes = "Prospective observational pharmacomicrobiome study; stool for shotgun metagenomics collected within 48 h of PK sampling."
  )

  ini({
    # MPA (reference 70 kg)
    lcl <- log(35.19); label("MPA metabolic clearance CL_MPA (L/h/70 kg)") # Table 3 'CL MPA' = 35.19
    lvc <- log(26.88); label("MPA central volume Vc_MPA (L/70 kg)") # Table 3 'Vc MPA' = 26.88

    # MPAG
    lcl_mpag <- log(3.69); label("MPAG systemic clearance CL_MPAG (L/h/70 kg)") # Table 3 'CL MPAG' = 3.69
    lvc_mpag <- log(18.83); label("MPAG central volume Vc_MPAG (L/70 kg)") # Table 3 'Vc MPAG' = 18.83

    # AcMPAG
    lcl_acmpag <- log(73.78); label("AcMPAG systemic clearance CL_AcMPAG (L/h/70 kg)") # Table 3 'CL AcMPAG' = 73.78
    lvc_acmpag <- log(19.2); label("AcMPAG central volume Vc_AcMPAG (L/70 kg)") # Table 3 'Vc AcMPAG' = 19.2
    lq_acmpag <- log(37.11); label("AcMPAG inter-compartmental clearance Q (L/h/70 kg)") # Table 3 'Q AcMPAG' = 37.11
    lvp_acmpag <- log(77.26); label("AcMPAG peripheral volume Vp_AcMPAG (L/70 kg)") # Table 3 'Vp AcMPAG' = 77.26

    # Fractions of MPA metabolic clearance forming each glucuronide
    fm_mpag <- fixed(0.87); label("Fraction of MPA metabolised to MPAG (molar)") # Suppl. File 3 Structural Model 'fraction of metabolism of 0.87'; Figure 1a
    fm_acmpag <- fixed(0.13); label("Fraction of MPA metabolised to AcMPAG (molar)") # Suppl. File 3 Structural Model '0.13'; Figure 1a

    # Enterohepatic circulation
    logitfrehc <- qlogis(0.39); label("Fraction of MPAG recirculated to MPA, FrEHC, at median community abundances (logit scale)") # Table 3 'FrEHC' = 0.39; main text Eq.1-2
    lmtime1 <- log(4.35); label("Gallbladder emptying start time after each dose, MTIME(1) (h)") # Table 3 'Gallbladder emptying time (h)' = 4.35
    dge <- fixed(0.5); label("Duration of gallbladder emptying (h)") # Suppl. File 3: MTIME(2) is 30 minutes after MTIME(1); Figure 1a
    lkehc <- fixed(log(3.71)); label("Gallbladder-to-gut emptying rate constant KGG (1/h, 70 kg)") # Table 3 'K45 (h-1) [fixed]' = 3.71; Figure 1 legend 'KGG ... fixed at 3.71 h-1'
    lkicv <- log(0.18); label("Gut MPAG-to-MPA central reactivation rate constant K51 (1/h, 70 kg)") # Table 3 'K51 (h-1)' = 0.18

    # Microbiome community effects on logit(FrEHC)
    e_stool_community_bacteroides_uniformis_pct_frehc <- 0.19; label("Community 1 log-abundance effect on logit(FrEHC) (unitless)") # Table 3 'Bacterial community 1' = 0.19
    e_stool_community_bacteroides_vulgatus_pct_frehc <- 0.10; label("Community 2 log-abundance effect on logit(FrEHC) (unitless)") # Table 3 'Bacterial community 2' = 0.10

    # IIV: omega^2 = log(1 + CV^2) from the Table 3 CVs (the tabled omega^2 are
    # rounded to two decimals; Table 3 footnote d defines CV = sqrt(exp(omega^2) - 1)).
    # The CL_MPA / CL_MPAG covariance is the tabled correlation 0.536 times
    # sqrt(omega2_CL_MPA * omega2_CL_MPAG) (tabled covariance 0.04).
    etalcl + etalcl_mpag ~ c(0.0170154, 0.0367442, 0.276188) # Table 3: CL MPA CV 13.1%; covariance 0.04 (correlation 53.6%); CL MPAG CV 56.4%
    etalvc_mpag ~ 0.103964 # Table 3 'Vc MPAG' CV 33.1%
    etalmtime1 ~ 1.38971 # Table 3 'Gallbladder emptying time' CV 173.6%
    etalogitfrehc ~ 0.584381 # Table 3 'FrEHC' CV 89.1% (eta added on the logit scale, main text Eq.3)
    etalcl_acmpag ~ 0.222344 # Table 3 'CL AcMPAG' CV 49.9%

    # Residual error: proportional, propSd = tabled CV = sqrt(sigma^2) (Table 3 footnote e)
    propSd <- 0.254; label("Proportional residual error, MPA (fraction)") # Table 3 'Proportional MPA plasma' CV 25.4%
    propSd_mpag <- 0.108; label("Proportional residual error, MPAG (fraction)") # Table 3 'Proportional MPAG plasma' CV 10.8%
    propSd_acmpag <- 0.222; label("Proportional residual error, AcMPAG (fraction)") # Table 3 'Proportional AcMPAG plasma' CV 22.2%
  })
  model({
    # Molecular weights (g/mol) for the molar bookkeeping; the paper states
    # the molar conversion but not the weights. MMF C23H31NO7, MPA C17H20O6,
    # MPAG and AcMPAG C23H28O12 (both glucuronides share one formula).
    mw_mmf <- 433.49
    mw_mpa <- 320.34
    mw_mpag <- 496.46
    mw_acmpag <- 496.46

    # Fixed allometric scaling to 70 kg (Supplementary File 3 Eq.3). The
    # source scales "rate constants" by WT^-0.25 without listing them; both
    # kicv (K51) and kehc (KGG) are scaled.
    wt_cl <- (WT / 70)^0.75
    wt_v <- WT / 70
    wt_k <- (WT / 70)^(-0.25)

    cl <- exp(lcl + etalcl) * wt_cl
    vc <- exp(lvc) * wt_v
    cl_mpag <- exp(lcl_mpag + etalcl_mpag) * wt_cl
    vc_mpag <- exp(lvc_mpag + etalvc_mpag) * wt_v
    cl_acmpag <- exp(lcl_acmpag + etalcl_acmpag) * wt_cl
    vc_acmpag <- exp(lvc_acmpag) * wt_v
    q_acmpag <- exp(lq_acmpag) * wt_cl
    vp_acmpag <- exp(lvp_acmpag) * wt_v
    kicv <- exp(lkicv) * wt_k
    kehc <- exp(lkehc) * wt_k
    mtime1 <- exp(lmtime1 + etalmtime1)

    # FrEHC on the logit scale with the log-abundance covariate terms
    # centred at the cohort medians 19.97% and 13.23% (main text Eq.1-4;
    # Supplementary Table S4).
    logitfrehc_ind <- logitfrehc + etalogitfrehc +
      e_stool_community_bacteroides_uniformis_pct_frehc *
        log(STOOL_COMMUNITY_BACTEROIDES_UNIFORMIS_PCT / 19.97) +
      e_stool_community_bacteroides_vulgatus_pct_frehc *
        log(STOOL_COMMUNITY_BACTEROIDES_VULGATUS_PCT / 13.23)
    frehc <- expit(logitfrehc_ind)

    # Micro-constants (Figure 1a). The MPAG -> bile rate KEHC (Suppl. File 3
    # Eq.2) is `kbm` here: FrEHC = kbm / (kbm + kel_mpag).
    kel <- cl / vc
    kel_mpag <- cl_mpag / vc_mpag
    kbm <- frehc * kel_mpag / (1 - frehc)
    kel_acmpag <- cl_acmpag / vc_acmpag
    k12_acmpag <- q_acmpag / vc_acmpag
    k21_acmpag <- q_acmpag / vp_acmpag

    # Gallbladder emptying gate (NONMEM MTIME): open from mtime1 to
    # mtime1 + dge after each dose (Figure 1a inset, 'time after dose').
    # tad() is measured from the start of the most recent infusion.
    tafterdose <- tad()
    gb_open <- (tafterdose >= mtime1) * (tafterdose < (mtime1 + dge))

    d/dt(central) <- -kel * central + kicv * gut_mpag
    d/dt(central_mpag) <- fm_mpag * kel * central - kel_mpag * central_mpag -
      kbm * central_mpag
    d/dt(central_acmpag) <- fm_acmpag * kel * central -
      kel_acmpag * central_acmpag - k12_acmpag * central_acmpag +
      k21_acmpag * peripheral1_acmpag
    d/dt(peripheral1_acmpag) <- k12_acmpag * central_acmpag -
      k21_acmpag * peripheral1_acmpag
    # 75% of the biliary MPAG flux is stored in the gallbladder and 25% drains
    # continuously into the gut (Supplementary File 3; Figure 1a).
    d/dt(gallbladder_mpag) <- 0.75 * kbm * central_mpag -
      kehc * gb_open * gallbladder_mpag
    d/dt(gut_mpag) <- 0.25 * kbm * central_mpag +
      kehc * gb_open * gallbladder_mpag - kicv * gut_mpag

    # Concentrations in mg/L (= ug/mL) of each analyte
    Cc <- central / vc * mw_mpa / mw_mmf
    Cc_mpag <- central_mpag / vc_mpag * mw_mpag / mw_mmf
    Cc_acmpag <- central_acmpag / vc_acmpag * mw_acmpag / mw_mmf

    Cc ~ prop(propSd)
    Cc_mpag ~ prop(propSd_mpag)
    Cc_acmpag ~ prop(propSd_acmpag)
  })
}

Gidal_2018_eslicarbazepine <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for ESLICARBAZEPINE, the primary active ",
    "metabolite of the once-daily antiepileptic prodrug eslicarbazepine ",
    "acetate (ESL), in adults with focal-onset seizures and in healthy ",
    "subjects (Gidal 2018, n = 1,039, 5,965 plasma concentrations pooled ",
    "across 11 phase 1 studies and the three phase 3 adjunctive-therapy ",
    "trials 2093-301 / 2093-302 / 2093-304). ESL is hydrolysed to ",
    "eslicarbazepine on first pass, so only the metabolite is modelled and ",
    "CL/F and V/F are apparent quantities referenced to the ESL dose. ",
    "Apparent oral clearance is ADDITIVE in its covariate terms rather ",
    "than the usual multiplicative-exponential form: a 2.43 L/h base is ",
    "incremented by a power function of the concomitant daily ",
    "carbamazepine dose, by a fixed 1.24 L/h for phenobarbital-like ",
    "enzyme inducers, and by 0.0132 L/h per kg above 70 kg, with the whole ",
    "bracket then scaled by a power function of creatinine clearance ",
    "(Gidal 2018 Appendix S1 Eq. E-1). Apparent volume is additive in sex ",
    "and inducer status and then scaled allometrically by weight ",
    "(Eq. E-2). Residual variability is stratified by development phase: ",
    "the richly sampled phase 1 studies and the trough-only phase 3 ",
    "studies carry separate proportional and additive components. This is ",
    "the PK backbone of the Gidal_2018_eslicarbazepine_* family; its ",
    "empirical-Bayes exposure metrics (AUC0-24, Cmax, Cav-ss) drive seven ",
    "companion exposure-response models for safety and efficacy."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Model equations and parameter tables are in Appendix S1",
    "(supporting information), Table S-2 and Equations E-1 and E-2.",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"

  # The paper estimates four residual-error components -- a proportional and
  # an additive term in each of two development-phase strata -- and selects
  # between the strata with STUDY_PHASE3 inside model(). rxode2 attaches one
  # error model to Cc, so the selected SD is assembled as a derived variable
  # (the Friberg_2012_voriconazole.R pattern).
  paper_specific_residual_sds <- c(
    "propSdPhase1",
    "addSdPhase1",
    "propSdPhase3",
    "addSdPhase3"
  )

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Acts on BOTH apparent clearance and apparent volume, but in two",
        "different functional forms, and the difference is load-bearing.",
        "On CL/F the effect is a LINEAR additive slope of 0.0132 L/h per kg",
        "centred at 70 kg, inside the additive bracket of Eq. E-1; on V/F it",
        "is a POWER term (WT / 70)^0.617 multiplying the whole additive",
        "bracket of Eq. E-2. Gidal 2018 Discussion notes that these",
        "exponents were estimated rather than fixed at allometric values,",
        "which is the change from the earlier three-study model. Reference",
        "70 kg is approximately the population median (mean 72.69 kg,",
        "SD 16.01; Table S1). Observed range 34-140 kg, stated in the",
        "Discussion as the span of virtual patients used to test the",
        "covariate. At 34 and 140 kg the model predicts AUCss 24.3% higher",
        "and 27.5% lower than at 70 kg, reproducing Appendix S1 exactly."
      ),
      source_name = "wt (Eq. E-1, E-2)"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance. RAW mL/min, NOT normalised to body surface",
        "area -- the register's Cockcroft-Gault alias branch."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as the power term (CRCL / 115.7)^0.195 scaling the ENTIRE",
        "additive clearance bracket of Eq. E-1, so it multiplies the",
        "carbamazepine and inducer increments as well as the base",
        "clearance. Reference 115.7 mL/min is the population median",
        "(mean 116.95, SD 26.18; Table S1). Gidal 2018 does not name the",
        "estimating equation; the cohort distribution (median 115.7,",
        "essentially normal renal function in 83.7% of subjects) is",
        "consistent with an absolute Cockcroft-Gault clearance rather than",
        "a BSA-normalised eGFR, and the renal-function categories in",
        "Table S1 are defined on the unnormalised mL/min scale",
        "(normal >= 90, mild 60-89, moderate 30-59). The positive exponent",
        "means CL/F FALLS with worsening renal function: 2.45, 2.26 and",
        "2.06 L/h at 120, 80 and 50 mL/min, reproducing Appendix S1.",
        "Only 0.9% of the analysis population had moderate impairment and",
        "none had severe, so the exponent is extrapolated below about",
        "30 mL/min."
      ),
      source_name = "CrCL (Eq. E-1)"
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Gidal 2018 writes the indicator as flag_sexf with 0 = male and",
        "1 = female, the same orientation as the canonical, so no value",
        "transformation is needed. Acts ONLY on V/F, as an additive shift",
        "of -9.9 L applied BEFORE the weight power term, so a woman and a",
        "man of the same body weight differ by 9.9 L * (WT / 70)^0.617",
        "rather than by a constant. Not retained on CL/F. 44.5% of the PK",
        "analysis population was female (Table S1, 577 of 1,039 male)."
      ),
      source_name = "flag_sexf (Eq. E-2)"
    ),
    DOSE_CBZ_MGD = list(
      description = "Total daily dose of concomitant carbamazepine.",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Set to 0 for subjects not taking carbamazepine; the power term",
        "(DOSE_CBZ_MGD / 800)^0.411 is then exactly 0 and the increment",
        "vanishes, so no separate carbamazepine indicator is needed on the",
        "PK model. Reference 800 mg/day is carbamazepine 400 mg twice",
        "daily, the modal regimen. The paper characterises the effect over",
        "400 mg/day (200 mg BID) to 1,200 mg/day (400 mg TID), where",
        "predicted CL/F is 33.4% and 52.5% above the no-carbamazepine",
        "value; both reproduce from Eq. E-1 exactly. 40.6% of the PK",
        "analysis population and 49.6% of the phase 3 subset took",
        "carbamazepine (Table S1). Mechanism per the Discussion:",
        "carbamazepine induces the UGT enzymes that glucuronidate",
        "eslicarbazepine, roughly one third of its elimination."
      ),
      source_name = "dose_carbamazepine (Eq. E-1)"
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no phenobarbital)",
      notes = paste(
        "One of the three drugs that set Gidal 2018's composite",
        "flag_phenobarbital-like indicator. The paper carries a single",
        "flag for 'phenobarbital or phenobarbital-like metabolic inducers",
        "(phenytoin or primidone)' (Table S-2 footnote a); this extraction",
        "decomposes it into the three constituent CONMED_<INN> indicators",
        "and rebuilds the flag inside model() as",
        "1 - (1 - CONMED_PB) * (1 - CONMED_PHT) * (1 - CONMED_PRM),",
        "which is exactly a logical OR of the three. Decomposing rather",
        "than minting a group canonical keeps the column meaning",
        "unambiguous and lets a user encode a real regimen. NOTE that",
        "carbamazepine is NOT part of this flag even though it is also an",
        "enzyme inducer -- it has its own dose-dependent term in the same",
        "equation -- so CONMED_EIAED would be the wrong canonical here.",
        "14.8% of the PK analysis population took phenobarbital",
        "(Table S1)."
      ),
      source_name = "flag_phenobarbital-like (component)"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no phenytoin)",
      notes = paste(
        "Second constituent of Gidal 2018's composite phenobarbital-like",
        "inducer flag; see the CONMED_PB notes for the reconstruction.",
        "The flag is an all-or-nothing indicator: the paper explicitly",
        "states the 51.0% clearance increase is 'without considering the",
        "effect of phenytoin/primidone dose', so unlike carbamazepine no",
        "dose-response was estimated for these three drugs. 1.0% of the PK",
        "analysis population took phenytoin (Table S1)."
      ),
      source_name = "flag_phenobarbital-like (component)"
    ),
    CONMED_PRM = list(
      description = "Concomitant primidone; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no primidone)",
      notes = paste(
        "Third constituent of Gidal 2018's composite phenobarbital-like",
        "inducer flag; see the CONMED_PB notes for the reconstruction.",
        "Primidone is metabolised to phenobarbital, which is why the paper",
        "groups it here. 0.8% of the PK analysis population took primidone",
        "(Table S1), so this indicator is the sparsest of the three and",
        "the grouping is what makes the effect estimable."
      ),
      source_name = "flag_phenobarbital-like (component)"
    ),
    STUDY_PHASE3 = list(
      description = paste(
        "Development-phase stratum of the record; 1 = phase 3 study",
        "(2093-301 / 2093-302 / 2093-304), 0 = one of the 11 phase 1",
        "studies."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the pooled phase 1 studies)",
      notes = paste(
        "Selects the residual-error magnitudes ONLY; it touches no",
        "structural or covariate parameter, so the typical-value",
        "prediction is identical either way. The stratification exists",
        "because the sampling designs differ by an order of magnitude: the",
        "phase 1 studies drew 8-28 samples per subject from 30 min to 24 h",
        "post-dose, while the phase 3 studies contributed mostly a single",
        "pre-dose trough at randomisation plus one maintenance-period",
        "sample (Appendix S1, Study data and analysis). Set to 1 when",
        "simulating a clinical-trial-like sparse design and 0 when",
        "simulating a dense phase 1 profile. See the ini() block for how",
        "the four printed residual rows were back-solved."
      ),
      source_name = "Phase 1 / Phase 3 (Table S-2, RV block)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "eslicarbazepine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "eslicarbazepine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1039L,
    n_studies = 14L,
    n_observations = "5,965 plasma eslicarbazepine concentrations",
    age_range = "16-80 years (median 36); 97.7% were under 65 years",
    age_median = "36 years",
    weight_range = "34-140 kg (Gidal 2018 Discussion, observed span)",
    weight_median = "mean 72.69 kg, SD 16.01 (Table S1); reference 70 kg",
    sex_female_pct = 44.5,
    race_ethnicity = c(
      Caucasian = 81.7,
      Black = 5.1,
      `Asian/Pacific Islander` = 6.3,
      Other = 4.9,
      Unknown = 2.0
    ),
    disease_state = paste(
      "224 healthy subjects or special populations from 11 phase 1 studies",
      "(PK data only) plus 815 adults with focal-onset seizures from three",
      "randomised, double-blind, placebo-controlled adjunctive-therapy",
      "trials; phase 3 patients had at least four focal-onset seizures in",
      "the 4 weeks before screening despite 1-3 concomitant antiepileptic",
      "drugs"
    ),
    dose_range = paste(
      "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily",
      "(400 mg in studies 301 and 302 only); phase 1 studies spanned the",
      "same 400-1,200 mg range"
    ),
    regions = paste(
      "Western Europe, North America, Latin America and Rest of World;",
      "studies 301 and 302 were European, study 304 was North American"
    ),
    renal_function = paste(
      "creatinine clearance mean 116.95 mL/min (SD 26.18), median 115.7;",
      "83.7% normal (>= 90 mL/min), 15.4% mild impairment (60-89), 0.9%",
      "moderate impairment (30-59), none severe (Table S1)"
    ),
    co_medication = paste(
      "carbamazepine 40.6%, gabapentin 20.8%, topiramate 18.2%, lamotrigine",
      "14.0%, phenobarbital 14.8%, clobazam 11.3%, valproate 11.4%,",
      "levetiracetam 7.0%, clonazepam 4.9%, zonisamide 1.6%, phenytoin",
      "1.0%, primidone 0.8% (Table S1, overall column). Oxcarbazepine was",
      "excluded by protocol because of its metabolic similarity to ESL, and",
      "felbamate was excluded from studies 301 and 302"
    ),
    notes = paste(
      "500 of the 1,039 subjects were also in the earlier three-study model",
      "of Nunes 2013. Study 2093-303 was excluded for Good Clinical",
      "Practice deficiencies found in a sponsor audit and was replaced by",
      "study 304. Estimation was FOCE / FOCE-with-interaction in NONMEM",
      "version 6 level 2.0. All parameters were estimated with less than",
      "36% SEM."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-2 and Equations E-1 / E-2:
    #
    #   CL/F = [2.43 + 1.08*(Dcbz/800)^0.411 + 1.24*flag_pb
    #                + 0.0132*(WT - 70)] * (CrCL/115.7)^0.195
    #   V/F  = (61.3 - 9.9*flag_sexf + 12.0*flag_pb) * (WT/70)^0.617
    #
    # Every published derived number in Appendix S1 reproduces from these
    # two equations to the printed precision; the checks are recorded on
    # the individual parameter lines and in the vignette source trace.
    # ==================================================================

    # ----- Structural -----
    lka <- log(2.34) ; label("Absorption rate constant (1/h)")  # Table S-2, k_a population mean 2.34 1/h, 9.6% SEM
    lcl <- log(2.43) ; label("Apparent oral clearance of eslicarbazepine for a 70 kg subject with creatinine clearance 115.7 mL/min and no concomitant antiepileptic drug (L/h)")  # Table S-2, 'CL/F for no carbamazepine use' 2.43 L/h, 1.3% SEM; Eq. E-1
    lvc <- log(61.3) ; label("Apparent volume of distribution for a 70 kg male not taking phenobarbital-like inducers (L)")  # Table S-2, 'V/F' 61.3 L, 2.0% SEM; Eq. E-2. Implied terminal half-life log(2)*61.3/2.43 = 17.5 h, inside the 13-20 h plasma half-life quoted in the Introduction

    # ----- Covariate effects on CL/F (all ADDITIVE, in L/h, inside the
    # ----- bracket; the CrCL power term multiplies the whole bracket) -----
    e_dose_cbz_mgd_cl <- 1.08 ; label("Additional apparent oral clearance at a concomitant carbamazepine dose of 800 mg/day (L/h)")  # Table S-2, 'Additional CL/F when carbamazepine dose = 800 mg' 1.08 L/h, 5.4% SEM; Eq. E-1
    e_dose_cbz_mgd_cl_pow <- 0.411 ; label("Power exponent on (carbamazepine daily dose / 800 mg/day) in the clearance increment (unitless)")  # Table S-2, 'Power term for effect of carbamazepine dose on CL/F' 0.411, 35.8% SEM. Check: at 400 and 1,200 mg/day the bracket rises to 3.242 and 3.706 L/h, i.e. 33.4% and 52.5% above 2.43, matching Appendix S1 exactly
    e_pblike_cl <- 1.24 ; label("Additional apparent oral clearance with concomitant phenobarbital, phenytoin or primidone (L/h)")  # Table S-2, 'Additive shift of concomitant phenobarbital or phenobarbital-like EIAED on CL/F' 1.24 L/h, 6.7% SEM. Check: (2.43+1.24)/2.43 = 1.510, the 51.0% increase stated in Appendix S1, and the reciprocal 33.8% AUCss reduction
    e_wt_cl <- 0.0132 ; label("Slope of body weight on apparent oral clearance, centred at 70 kg (L/h/kg)")  # Table S-2, 'Slope term for effect of body weight on CL/F' 0.0132 L/h/kg, 24.0% SEM. Check: 62, 70, 81 kg give 2.32, 2.43, 2.58 L/h as printed
    e_crcl_cl <- 0.195 ; label("Power exponent on (creatinine clearance / 115.7 mL/min) applied to the whole clearance bracket (unitless)")  # Table S-2, 'Power term for effect of CrCL on CL/F' 0.195, 33.9% SEM. Check: 120, 80, 50 mL/min give 2.45, 2.26, 2.06 L/h as printed

    # ----- Covariate effects on V/F (additive shifts in L applied BEFORE
    # ----- the weight power term) -----
    e_sexf_vc <- -9.9 ; label("Additive shift in apparent volume of distribution for females versus males of the same body weight (L)")  # Table S-2, 'Additive shift of female gender on V/F' -9.9 L, 18.2% SEM; Eq. E-2
    e_pblike_vc <- 12.0 ; label("Additive shift in apparent volume of distribution with concomitant phenobarbital, phenytoin or primidone (L)")  # Table S-2, 'Additive shift of concomitant phenobarbital or phenobarbital-like EIAED on V/F' 12.0 L, 30.3% SEM. Check: 12.0/61.3 = 19.6%, the increase stated in Appendix S1
    e_wt_vc <- 0.617 ; label("Power exponent on (body weight / 70 kg) applied to the whole volume bracket (unitless)")  # Table S-2, 'Power term for effect of body weight on V/F' 0.617, 15.0% SEM. Check: 62 and 81 kg give 56.9 and 67.1 L, the two printed values -- note Appendix S1 mislabels those weights as 61 and 79 kg; see vignette Errata item 1

    # ----- Interindividual variability -----
    # Table S-2 reports IIV as %CV; the variance on the log scale is
    # omega^2 = log(CV^2 + 1). Values below are that conversion.
    etalka ~ 0.9555   # Table S-2, k_a IIV 126.49 %CV, 18.4% SEM: log(1.2649^2 + 1) = 0.9555
    etalcl ~ 0.070573 # Table S-2, CL/F IIV 27.04 %CV, 10.5% SEM: log(0.2704^2 + 1) = 0.070573
    etalvc ~ 0.030814 # Table S-2, V/F IIV 17.69 %CV, 15.6% SEM: log(0.1769^2 + 1) = 0.030814

    # ----- Residual variability, stratified by development phase -----
    # Table S-2 prints the RV block in an unusual parameterisation: a
    # proportional VARIANCE component and a RATIO of the additive to the
    # proportional component on the SD scale. The additive SD is therefore
    # ratio * sqrt(proportional variance), and the total variance is
    #   Var = sigma_prop^2 * Cp^2 + sigma_add^2.
    # Both strata were confirmed against the table's own %CV footnotes,
    # which pin the curve at two concentrations each and over-determine
    # the pair:
    #   Phase 1, footnote b: 23.74 %CV at 2,400 ng/mL and 11.24 %CV at
    #     33,000 ng/mL solve to prop variance 0.012402 (printed 0.0124)
    #     and additive variance 253,214 ng^2/mL^2, i.e. SD 503.2 ng/mL.
    #     The printed ratio 4,520 recovers as 503.2/sqrt(0.0124) = 4,519.
    #   Phase 3, footnote c: 311.15 %CV at 740 ng/mL and 15.68 %CV at
    #     39,200 ng/mL solve to prop variance 0.021144 and additive
    #     variance 5,289,922 ng^2/mL^2 (printed 5,290,000), i.e. SD
    #     2,300 ng/mL. The printed ratio 0.0000632 recovers as
    #     sqrt(0.021144)/2,300 = 0.00006322.
    # See the vignette source trace for the full arithmetic.
    propSdPhase1 <- 0.111355 ; label("Proportional residual SD in the phase 1 studies (fraction)")  # Table S-2, 'Proportional RV component (sigma_1), Phase 1' 0.0124 as a variance, 7.8% SEM; sqrt(0.0124) = 0.111355
    addSdPhase1 <- 503.2 ; label("Additive residual SD in the phase 1 studies (ng/mL)")  # Table S-2, 'Ratio of additive/proportional RV components, Phase 1' 4,520, 18.4% SEM; 4520 * 0.111355 = 503.3, and footnote b back-solves 503.2
    propSdPhase3 <- 0.145410 ; label("Proportional residual SD in the phase 3 studies (fraction)")  # Table S-2 footnote c back-solve: proportional variance 0.021144, sqrt = 0.145410. The proportional component is not printed directly for phase 3; only the ratio and the additive term are
    addSdPhase3 <- 2300.0 ; label("Additive residual SD in the phase 3 studies (ng/mL)")  # Table S-2, 'Additive RV component (sigma_2), Phase 3' 5,290,000 as a variance, 17.1% SEM; sqrt(5,290,000) = 2,300.0
  })

  model({
    # ----- Composite phenobarbital-like inducer flag (Table S-2 footnote a)
    # Logical OR of the three constituent comedications, written as a
    # product of complements so it stays exactly 0 or 1.
    flag_pblike <- 1 - (1 - CONMED_PB) * (1 - CONMED_PHT) * (1 - CONMED_PRM)

    # ----- Individual parameters -----
    ka <- exp(lka + etalka)

    # Eq. E-1. The carbamazepine term is 0 when DOSE_CBZ_MGD is 0, so no
    # separate carbamazepine indicator is needed. The creatinine-clearance
    # power term scales the whole bracket, inducer increments included.
    cl <- (exp(lcl) +
      e_dose_cbz_mgd_cl * (DOSE_CBZ_MGD / 800)^e_dose_cbz_mgd_cl_pow +
      e_pblike_cl * flag_pblike +
      e_wt_cl * (WT - 70)) *
      (CRCL / 115.7)^e_crcl_cl *
      exp(etalcl)

    # Eq. E-2. The sex and inducer shifts are absolute litres applied
    # before the allometric term, not multipliers.
    vc <- (exp(lvc) +
      e_sexf_vc * SEXF +
      e_pblike_vc * flag_pblike) *
      (WT / 70)^e_wt_vc *
      exp(etalvc)

    kel <- cl / vc

    # ----- ODE system -----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ----- Observation -----
    # Dose is in mg and vc in L, so central/vc is ug/mL; the factor 1000
    # puts Cc on the ng/mL scale of Table S-2 and of every exposure metric
    # consumed by the companion exposure-response models.
    Cc <- 1000 * central / vc

    # Phase-stratified residual. Var = (prop * Cp)^2 + add^2 is the NONMEM
    # form the table parameterises; assembling one SD keeps it in a single
    # error model.
    propSdPhase <- propSdPhase1 * (1 - STUDY_PHASE3) + propSdPhase3 * STUDY_PHASE3
    addSdPhase <- addSdPhase1 * (1 - STUDY_PHASE3) + addSdPhase3 * STUDY_PHASE3
    sdCc <- sqrt((propSdPhase * Cc)^2 + addSdPhase^2)

    Cc ~ add(sdCc)
  })
}

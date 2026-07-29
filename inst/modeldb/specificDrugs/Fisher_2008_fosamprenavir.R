Fisher_2008_fosamprenavir <- function() {
  description <- "Two-compartment population PK model with first-order absorption for orally administered fosamprenavir (FPV), measured as the active amprenavir (APV) metabolite, in HIV-1-infected pediatric patients aged 4 weeks to 18 years (Fisher 2008). Allometric scaling on apparent clearance (CL/F, Q) at a fixed exponent of 0.75 and on apparent volumes (V2/F, V3) at a fixed exponent of 1.0 (reference 70 kg). Apparent CL/F is reduced ~60% by concomitant ritonavir (RTV) co-administration (maximal CYP3A4 inhibition assumed at the RTV doses used), and is further modified by a piecewise age-maturation factor (linearly declining additive offset for AGE <= 2*AG50, zero above), by sex (lower in females), by race (separate multipliers for Black and for the non-Black non-White composite vs the White reference), and by a power effect of serum alpha-1-acid glycoprotein (AAG, centred at 0.77 g/L). Apparent V2/F also carries a power effect of AAG. Bioavailability is anchored on suspension-under-fed conditions (F=1), with a separate relative bioavailability for the tablet formulation (F_tab) and a separate relative bioavailability for the suspension administered fasted (F_food,sus). Inter-occasion variability on CL/F (~34% CV) reported by the source poster is NOT structurally encoded here (no operational occasion column is defined for the model-library use case); downstream users who want to simulate IOV can add an OCC indicator and a per-occasion eta in rxode2."
  reference <- "Fisher J, Gastonguay MR, Knebel W, Gibiansky L, Wire MB. Population Pharmacokinetic Modeling of Fosamprenavir in Pediatric HIV-Infected Patients. American Conference on Pharmacometrics (ACOP) poster, Tucson, AZ, 2008. Available at https://metrumrg.com/wp-content/uploads/2018/08/acop_2008_fosamprenavir.pdf"
  vignette <- "Fisher_2008_fosamprenavir"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject. Allometric power scaling with reference 70 kg and fixed theory-based exponents 0.75 on CL/F and Q, 1.0 on V2/F and V3 (Fisher 2008 Methods, 'Model and Modeling Assumptions'). Cohort baseline mean 34.96 kg, median 32.9 kg, range 5.9-102.8 kg (Table 1).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Postnatal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject (baseline age 0.72-18 years in the source cohort; Table 1). Drives a piecewise age-maturation factor on CL/F that is linearly decreasing in AGE for AGE <= 2*AG50 (AG50 = 2.05 years) and is identically 1 for AGE > 2*AG50. The reference covariate-effect category is AGE > ~4 years (where the age factor = 1).",
      source_name        = "AGE"
    ),
    AAG = list(
      description        = "Serum alpha-1-acid glycoprotein (orosomucoid) concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (median value per patient over occasions was used in Table 1's summary). Cohort baseline mean 0.87 g/L, median 0.80 g/L, range 0.41-2.69 g/L (Table 1). Enters CL/F and V2/F via two independent power-law forms centred at the cohort median 0.77 g/L (Figure 3 caption: 'AAG = 0.77 g/L (population median)'): (AAG/0.77)^(-0.626) on CL/F and (AAG/0.77)^(-0.369) on V2/F (Table 3 AAG_CL and AAG_V2).",
      source_name        = "AAG"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male (Fisher 2008 Table 3 reference category)",
      notes              = "Time-fixed per subject. Fisher 2008 Table 2 encodes SEX as 0 = males / 1 = females, which matches the canonical SEXF orientation (1 = female). Female patients have ~15% lower CL/F than males at the same WT, AAG, race, age, and RTV status (Table 3 theta_12 = 0.846).",
      source_name        = "SEX"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-Black (White, Asian, Hispanic, American Indian, or Other in this cohort; reference category is White)",
      notes              = "Time-fixed per subject. Fisher 2008 Table 2 encodes race as a 6-level categorical: 1 = White (57.7%), 2 = Black (27%), 3 = Asian (1.5%), 4 = Hispanic (5.1%), 5 = American Indian (5.8%), 6 = Other (2.9%). The model collapses the 6 levels into 3 effective categories: White (reference), Black (RACE_BLACK = 1), and 'non-Black ethnic origin' (RACE_NONBLACK_NONWHITE = 1). RACE_BLACK on CL/F: Table 3 theta_13 = 0.940 (6% lower than White reference).",
      source_name        = "RACE (level 2)"
    ),
    RACE_NONBLACK_NONWHITE = list(
      description        = "Non-Black, non-White composite race indicator (Asian + Hispanic + American Indian + Other in the source cohort)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White or Black (the larger-N pooled reference in the source paper)",
      notes              = "Time-fixed per subject. Fisher 2008 'non-Black ethnic origin' composite of Race = 3 (Asian), 4 (Hispanic), 5 (American Indian), and 6 (Other) collapsed into a single indicator vs the White reference. Patients in this composite have ~6% higher CL/F than the White reference (Table 3 theta_14 = 1.06). Mutually exclusive with RACE_BLACK (a subject cannot be both Black and non-Black-non-White).",
      source_name        = "RACE (levels 3, 4, 5, 6)"
    ),
    CONMED_RTV = list(
      description        = "Concomitant low-dose ritonavir (RTV) co-administration indicator (CYP3A4 boost)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no concomitant ritonavir (unboosted FPV)",
      notes              = "Time-varying per dosing record. Fisher 2008 encodes RTV as 0 = administered without ritonavir (n=18 patients, 13.1%) / 1 = administered with ritonavir (n=119, 86.9%; Table 2). The source poster's Methods explicitly state 'Maximal inhibition of FPV CL/F was assumed at the RTV doses included in the model' (RTV mg ranged 0-200 in the dataset; Table 1), so the model treats RTV as a binary switch rather than a dose-response. Reduces typical CL/F from 84.4 L/h (without RTV) to 34.1 L/h (with RTV), a 60% reduction (Table 3 theta_1 / theta_6).",
      source_name        = "RTV"
    ),
    FORM_TABLET = list(
      description        = "Tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = suspension formulation (the typical-value bioavailability reference when combined with FED = 1)",
      notes              = "Time-varying per dosing record. Fisher 2008 encodes formulation as 1 = suspension (n=89 records, 65%) / 2 = tablets (n=48 records, 35%) in Table 2; this is recoded to FORM_TABLET = 1 if tablet, 0 if suspension. Tablet formulation has F_tab = 1.09 relative to the fed-suspension reference (Table 3 theta_7).",
      source_name        = "Formulation (level 2)"
    ),
    FED = list(
      description        = "Fed-state indicator at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Time-varying per dosing record. Fisher 2008 encodes 'Food Intake' as -1 = Missing (43 records, 31.4%) / 0 = Administered without food (20 records, 14.6%) / 1 = Administered with food (74 records, 54%) in Table 2. This is recoded to FED = 1 if administered with food, 0 if fasted; the missing food category is treated as the fed reference (FED = 1) in this model file (downstream users with a different imputation should set FED accordingly). The food effect only applies to suspension administration: F_food,sus = 0.87 (Table 3 theta_8) is the relative bioavailability of the fasted suspension vs the fed suspension; tablet bioavailability F_tab does NOT additionally depend on food in the source paper's structural form.",
      source_name        = "Food Intake"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 137L,
    n_studies        = 3L,
    n_observations   = 1322L,
    age_range        = "0.72-18 years (baseline; eligibility 4 weeks to 18 years)",
    age_median       = "10 years",
    weight_range     = "5.9-102.8 kg (baseline)",
    weight_median    = "32.9 kg",
    sex_female_pct   = 54.7,
    race_ethnicity   = c(White = 57.7, Black = 27.0, Asian = 1.5, Hispanic = 5.1, AmericanIndian = 5.8, Other = 2.9),
    disease_state    = "HIV-1-infected pediatric patients, protease-inhibitor (PI) naive or PI-experienced, on antiretroviral therapy that included fosamprenavir (FPV) +/- ritonavir (RTV) booster. 119 of 137 patients (86.9%) received FPV with RTV; 18 (13.1%, all aged 2-6 years) received FPV alone.",
    dose_range       = "Oral FPV (suspension or tablet) at per-kilogram doses based on age group: FPV/RTV BID 18-36 mg/kg (max 700 mg), FPV/RTV QD 36-72 mg/kg (max 1400 mg), or FPV BID 17-38 mg/kg (max 1400 mg). RTV dose (when given) 0-200 mg (cohort mean 109.5 mg, median 100 mg; Table 1).",
    regions          = "Multinational: Argentina, Belgium, Canada, Italy, Mexico, Netherlands, Portugal, Romania, Russia, Spain, USA (per Acknowledgments).",
    aag_distribution = "Baseline AAG mean 0.87 g/L, median 0.80 g/L, range 0.41-2.69 g/L (Table 1); reference (cohort median used in the model centring) 0.77 g/L (Figure 3 caption).",
    iov_structure    = "Inter-occasion variability on CL/F was identified at omega^2 = 0.114 (~34% CV) in the final model (Table 3 theta_9). This model file does NOT encode IOV structurally -- the source poster does not define an operational occasion column for downstream simulation use, and the nlmixr2lib convention (Brooks 2021 / Andrews 2017 precedent) is to omit IOV when no occasion mapping is defined; see vignette Assumptions and deviations. Downstream users who want to simulate IOV can add an OCC indicator and a per-occasion eta in rxode2.",
    studies          = "APV20002 (Phase II, 4 weeks to <2 years, FPV +/- RTV, n=9 patients, 6.6%); APV20003 (Phase II, 2-18 years, FPV/RTV QD or BID, n=59 patients, 43.1%); APV29005 (Phase II, FPV/RTV BID 2-18 years and FPV BID 2-<6 years PI-naive, n=69 patients, 50.4%). Each study contributed one extensive PK sampling day after at least 10 days of multiple-dose therapy plus up to 9 trough samples spanning up to 18 months.",
    notes            = "Pooled analysis from three multinational pediatric studies sponsored by GlaxoSmithKline. Fosamprenavir (FPV) is a prodrug of amprenavir (APV); APV is the measured species (plasma assay, units ug/mL). APV is ~90% protein-bound, primarily to alpha-1-acid glycoprotein (AAG). APV is metabolized primarily by CYP3A4. RTV is a potent CYP3A4 inhibitor and is co-administered at low 'booster' doses to increase APV exposure. The poster was used to support FDA approval of pediatric FPV dosing recommendations for children 2-18 years (FPV) and 6-18 years (FPV/RTV)."
  )

  ini({
    # Structural parameters (Fisher 2008 Table 3 final-model column).
    # Reference subject for the typical-value structural parameters:
    #   WT = 70 kg, AGE > 4 years (so the age factor = 1), male (SEXF = 0),
    #   White (RACE_BLACK = 0, RACE_NONBLACK_NONWHITE = 0), AAG = 0.77 g/L,
    #   FORM_TABLET = 0 + FED = 1 (suspension under fed conditions; F = 1
    #   reference), and CONMED_RTV = 0 (unboosted; CL/F = theta_6 = 84.4 L/h).
    # When CONMED_RTV = 1 (boosted, the more common cohort state), CL/F
    # reduces to theta_1 = 34.1 L/h via the e_conmed_rtv_cl effect below.
    lka <- log(1.13);  label("First-order absorption rate constant ka (1/h)")               # Table 3 theta_5 = 1.13 1/h (RSE 30%)
    lcl <- log(84.4);  label("Apparent oral clearance at CONMED_RTV = 0 (CL/F, L/h)")        # Table 3 theta_6 = 84.4 L/h (RSE 11%; CL/F without RTV)
    lvc <- log(288);   label("Apparent central volume of distribution at WT = 70 kg (V2/F, L)")  # Table 3 theta_2 = 288 L (RSE 23%)
    lq  <- log(63.5);  label("Apparent inter-compartmental clearance at WT = 70 kg (Q, L/h)")    # Table 3 theta_3 = 63.5 L/h (RSE 15%)
    lvp <- log(1630);  label("Apparent peripheral volume of distribution at WT = 70 kg (V3, L)") # Table 3 theta_4 = 1630 L (RSE 28%)

    # Allometric exponents fixed at theory-based values (Methods, 'Model and
    # Modeling Assumptions': "An allometric model was assumed to describe the
    # relationship between body weight and all clearance and volume parameters
    # using a normalized reference weight of 70 kg and fixed allometric power
    # values of 0.75 (CL and Q) and 1 (V2/F and V3)."). ka is NOT allometrically
    # scaled in Fisher 2008.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F with WT (unitless; fixed at theory)")  # Methods, Model and Modeling Assumptions
    e_wt_q  <- fixed(0.75); label("Allometric exponent on Q with WT (unitless; fixed at theory)")     # Methods, Model and Modeling Assumptions
    e_wt_vc <- fixed(1);    label("Allometric exponent on V2/F with WT (unitless; fixed at theory)")  # Methods, Model and Modeling Assumptions
    e_wt_vp <- fixed(1);    label("Allometric exponent on V3 with WT (unitless; fixed at theory)")    # Methods, Model and Modeling Assumptions

    # Concomitant ritonavir effect on CL/F. Fisher 2008 reports two CL/F values
    # in Table 3 (theta_1 = 34.1 L/h with RTV; theta_6 = 84.4 L/h without RTV).
    # Encoded here as a fractional-deviation form analogous to Colombo 2006
    # atazanavir: cl = exp(lcl + etalcl) * (1 + e_conmed_rtv_cl * CONMED_RTV).
    # The fractional change is (34.1 - 84.4) / 84.4 = -0.5961, i.e. CL/F is
    # ~60% lower when RTV is co-administered (consistent with the Conclusions
    # paragraph "co-administration of RTV was estimated to decrease plasma APV
    # CL/F by approximately 60%").
    e_conmed_rtv_cl <- -0.5961; label("Fractional change in CL/F when concomitant ritonavir is co-administered (unitless)")  # Derived from Table 3 theta_1 / theta_6: (34.1 - 84.4)/84.4 = -0.5961

    # Age maturation parameters for the piecewise age effect on CL/F. The
    # source poster's Model section gives the equation explicitly:
    #   For AGE <= 2*AG50:  CL = theta_TVCL * (WT/70)^0.75 * [1 + AMAX*(1 - 0.5*AGE/AG50)] * exp(eta_CL)
    #   For AGE >  2*AG50:  CL = theta_TVCL * (WT/70)^0.75 * exp(eta_CL)
    # The two pieces are continuous at AGE = 2*AG50 (the bracket = 1). AMAX is
    # the maximum age-effect multiplicative shift at birth (AGE = 0), AG50 is
    # the age at which the effect has decayed to half its birth value. The
    # piecewise structure is implemented in model() via max(0, ...).
    amax_cl <- 0.789; label("Maximum age-effect shift on CL/F at birth (unitless)")                                       # Table 3 theta_10 = 0.789 (RSE 65%)
    ag50_cl <- 2.05;  label("Age at which the age effect on CL/F has decayed to half its birth value (years)")            # Table 3 theta_11 = 2.05 (RSE 26%)

    # Sex effect on CL/F (multiplicative; CL/F multiplier for females vs males).
    e_sexf_cl <- 0.846; label("CL/F multiplier for females (vs male reference)")  # Table 3 theta_12 = 0.846 (RSE 7%)

    # Race effects on CL/F (multiplicative; CL/F multiplier vs the White
    # reference). The non-Black non-White composite pools Asian + Hispanic +
    # American Indian + Other into a single indicator.
    e_race_black_cl <- 0.940; label("CL/F multiplier for Black race (vs White reference)")                                # Table 3 theta_13 = 0.940 (RSE 8%)
    e_race_nbnw_cl  <- 1.06;  label("CL/F multiplier for non-Black non-White composite race (vs White reference)")        # Table 3 theta_14 = 1.06 (RSE 13%)

    # AAG power effects on CL/F and V2/F. Fisher 2008 Methods caption ("Other
    # continuous covariate effects were generally modeled using a normalized
    # power model"); centring value 0.77 g/L from Figure 3 caption (population
    # median used as the reference). Negative exponents reflect lower CL/F
    # and lower V2/F at higher AAG (which sequesters protein-bound APV).
    e_aag_cl <- -0.626; label("AAG power exponent on CL/F, centred at 0.77 g/L (unitless)")     # Table 3 theta_15 = -0.626 (RSE 9%)
    e_aag_vc <- -0.369; label("AAG power exponent on V2/F, centred at 0.77 g/L (unitless)")     # Table 3 theta_16 = -0.369 (RSE 92%)

    # Bioavailability terms. Reference state: suspension under fed conditions
    # (F = 1, anchor for the relative bioavailability terms). Tablet vs
    # suspension: F_tab = 1.09 (Table 3 theta_7). Fasted suspension vs fed
    # suspension: F_food,sus = 0.87 (Table 3 theta_8). Per the source poster's
    # Model paragraph, the food effect only applies to the suspension form;
    # tablet bioavailability does not additionally depend on food in the final
    # structural model.
    lftab     <- log(1.09); label("Relative bioavailability of tablet vs fed suspension (F_tab)")             # Table 3 theta_7 = 1.09 (RSE 8%)
    lffoodsus <- log(0.87); label("Relative bioavailability of fasted vs fed suspension (F_food,sus)")        # Table 3 theta_8 = 0.87 (RSE 8%)

    # Inter-individual variability. Block covariance matrix on (etalcl, etalvc)
    # per Methods "Model and Modeling Assumptions" ("A block covariance matrix
    # for the inter-individual random effects (Omega) for CL/F and V2/F was
    # included"). Diagonal omega on etalq. Inter-occasion variability on CL/F
    # (omega^2 = 0.114, Table 3 theta_9, ~34% CV) is NOT encoded structurally
    # here (see vignette Assumptions and deviations).
    #
    # Variance interpretation: Fisher 2008 reports the diagonal omegas with
    # both an omega^2 value and a CV-style label. For CL/F the label "30% CV"
    # is consistent with log-normal: sqrt(exp(0.0901) - 1) = 30.7%. For V2/F
    # the label "66% CV" corresponds to sqrt(0.438) = 66.2%, i.e. the paper
    # reports the omega^2-square-root as the CV (not log-normal); for Q the
    # label "73%" is sqrt(0.536) = 73.2%. The variance values 0.0901, 0.438,
    # 0.536 are taken verbatim from Table 3 as the omega^2 entries on the
    # internal log-scale; the (slightly inconsistent) CV labels are
    # informational and do not modify the encoded variances.
    etalcl + etalvc ~ c(0.0901,
                        0.0945, 0.438)   # Table 3: Omega_11 = 0.0901; Omega_12 = 0.0945; Omega_22 = 0.438
    etalq ~ 0.536                         # Table 3: Omega_33 = 0.536 (RSE 36%)

    # Residual error. Fisher 2008 Table 3 reports sigma^2 values for the
    # combined proportional + additive error model: sigma_1^2 = 0.0828
    # (proportional, on the fractional scale) and sigma_2^2 = 0.0760
    # (additive, in ug/mL). nlmixr2's propSd / addSd are SDs, so
    # propSd = sqrt(0.0828) = 0.2877 and addSd = sqrt(0.0760) = 0.2757.
    propSd <- 0.2877; label("Proportional residual SD (fraction)")  # Table 3 sigma_1^2 = 0.0828 (RSE 7%); propSd = sqrt(0.0828)
    addSd  <- 0.2757; label("Additive residual SD (ug/mL)")          # Table 3 sigma_2^2 = 0.0760 (RSE 12%); addSd  = sqrt(0.0760)
  })

  model({
    # Body-weight scaling reference: 70 kg (Methods, Model and Modeling Assumptions).
    wt70 <- WT / 70

    # Piecewise age maturation factor on CL/F.
    #   For AGE <= 2*AG50:  f_age = 1 + AMAX * (1 - 0.5 * AGE / AG50)
    #   For AGE >  2*AG50:  f_age = 1
    # Implemented via max(0, ...) so the bracketed term is clamped to zero
    # once AGE exceeds 2*AG50; the two pieces match at the boundary (the
    # bracket equals 0 at AGE = 2*AG50, and the source equations both give
    # f_age = 1 there).
    f_age <- 1 + amax_cl * max(0, 1 - 0.5 * AGE / ag50_cl)

    # Sex multiplicative factor on CL/F (1 for males, e_sexf_cl for females).
    f_sexf <- 1 + (e_sexf_cl - 1) * SEXF

    # Race multiplicative factor on CL/F. White reference: both indicators 0,
    # f_race = 1. Black: f_race = e_race_black_cl. Non-Black non-White
    # composite: f_race = e_race_nbnw_cl. The two indicators are mutually
    # exclusive (a subject is at most one of {Black, non-Black non-White}).
    f_race <- 1 + (e_race_black_cl - 1) * RACE_BLACK + (e_race_nbnw_cl - 1) * RACE_NONBLACK_NONWHITE

    # AAG power-law factors on CL/F and V2/F, centred at 0.77 g/L.
    f_aag_cl <- (AAG / 0.77) ^ e_aag_cl
    f_aag_vc <- (AAG / 0.77) ^ e_aag_vc

    # Concomitant ritonavir fractional-deviation factor on CL/F.
    # f_rtv = 1 when CONMED_RTV = 0 (off RTV reference); f_rtv = 1 + e_conmed_rtv_cl
    # = 0.404 when CONMED_RTV = 1 (on RTV).
    f_rtv <- 1 + e_conmed_rtv_cl * CONMED_RTV

    # Individual PK parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * wt70 ^ e_wt_cl * f_rtv * f_age * f_sexf * f_race * f_aag_cl
    vc <- exp(lvc + etalvc) * wt70 ^ e_wt_vc * f_aag_vc
    q  <- exp(lq  + etalq)  * wt70 ^ e_wt_q
    vp <- exp(lvp)          * wt70 ^ e_wt_vp

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral disposition (apparent CL/F, V2/F, Q, V3 -- F
    # bioavailability terms enter via f(depot) below).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability. Anchor F = 1 at suspension-under-fed; tablet has
    # F = F_tab regardless of FED; fasted suspension has F = F_food,sus.
    # log(F) = log(F_tab) * FORM_TABLET + log(F_food,sus) * (1 - FORM_TABLET) * (1 - FED)
    # Reproduces the three structural F levels exactly: 1.0 (suspension, fed),
    # 1.09 (tablet, any food state), 0.87 (suspension, fasted).
    f(depot) <- exp(lftab * FORM_TABLET + lffoodsus * (1 - FORM_TABLET) * (1 - FED))

    # Observed amprenavir (APV) concentration. Dose in mg, vc in L -> central/vc
    # is mg/L = ug/mL, matching the source poster's reported concentration units.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

Hazendonk_2016_factor_viii <- function() {
  description <- "Two-compartment population PK model for factor VIII (FVIII) concentrates in severe and moderate hemophilia A patients (adults and children, FVIII plasma concentration < 0.05 IU/mL) undergoing elective / minor / major surgery (Hazendonk 2016). PK parameters are allometrically scaled to a 68 kg reference body weight with fixed exponents of 0.75 on clearances and 1.0 on volumes; typical CL, V1, Q, V2 at the reference body weight are 150 mL/h, 2810 mL, 160 mL/h, 1900 mL. Clearance carries three covariate effects (age with power exponent -0.17 centered at 40 years; +26% for blood group O; -7% for a major surgical procedure) and central volume carries an age effect (power exponent -0.09 centered at 40 years). B-domain-deleted recombinant FVIII products (Refacto AF) are under-detected by the one-stage clotting assay by a fixed 34%, encoded as a multiplicative correction on the predicted concentration. IIV on CL and V1 is 37% and 27% (exponential); IIV on Q and V2 was not estimable. Residual error is combined (proportional 18% + additive 0.15 IU/mL) using the values reported for the majority center cluster (centers 1, 2, 3)."
  reference   <- "Hazendonk H, Fijnvandraat K, Lock J, Driessens M, van der Meer F, Meijer K, Kruip M, Laros-van Gorkom B, Peters M, de Wildt S, Leebeek F, Cnossen M, Mathot R; OPTI-CLOT study group. A population pharmacokinetic model for perioperative dosing of factor VIII in hemophilia A patients. Haematologica. 2016 Oct;101(10):1159-1169. doi:10.3324/haematol.2015.136275. PMID:27390359."
  vignette    <- "Hazendonk_2016_factor_viii"
  units       <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling at reference 68 kg (the total-cohort median): fixed theory-based exponent 0.75 on CL and Q, fixed exponent 1.0 on V1 and V2 (Hazendonk 2016 Structural model development and Table 5). Total-cohort weight median 75 kg (range 5-111 kg); adult median 80 kg (45-111 kg); pediatric median 18.5 kg (5-85 kg). Treated as time-fixed baseline body weight per subject.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL (exponent -0.17) and V1 (exponent -0.09) centered at 40 years (Hazendonk 2016 Table 5). Older patients have lower FVIII CL and V1; the paper reports typical CL of 214, 169, 150, 142 mL/h/68 kg for a non-O, minor-surgery patient aged 5, 20, 40, 55 years. Total-cohort age median 40 years (range 0.2-78); adult median 48 years (19-78); pediatric median 4.3 years (0.2-17.3).",
      source_name        = "AGE"
    ),
    BLOOD_GROUP_O = list(
      description        = "ABO blood group O indicator (1 = O, 0 = non-O: A, B, or AB)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-O ABO blood group: A, B, or AB)",
      notes              = "Hazendonk 2016 Table 5: `1.26 ^ blood_group` on CL, so blood group O subjects show 26% higher FVIII clearance than non-O. Mechanism per Hazendonk 2016 Discussion: blood group O individuals have ~25% lower baseline VWF (the FVIII-protective carrier protein), leading to accelerated FVIII proteolytic degradation. Cohort prevalence 50% (51 of 101 with recorded blood group; 34/68 adults, 17/33 pediatric). Time-fixed per subject.",
      source_name        = "BLOOD_GROUP_O"
    ),
    SURG_SEV_MAJOR = list(
      description        = "Major (vs minor) surgical procedure severity indicator (Koshy 1995 classification: major and high-risk collapsed to 1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (minor surgical procedure)",
      notes              = "Hazendonk 2016 Table 5: `0.93 ^ severity` on CL, so major-surgery cases show 7% lower typical FVIII clearance than minor cases. The severity classification follows Koshy et al. 1995 (Hazendonk 2016 Methods reference 19; major and high-risk categories collapsed together). Cohort prevalence 49% of surgical procedures (97 of 198; 61.4% adults, 19% pediatric). The paper's Discussion notes that the negative sign is confounded by age -- older patients underwent more major procedures and had lower CL -- so the univariate association was retained in the multivariate model but interpreted with caution. Time-fixed per surgical procedure (a subject can have different values across multiple procedures).",
      source_name        = "SURG_SEV_MAJOR"
    ),
    FORM_FVIII_BDD = list(
      description        = "B-domain-deleted (BDD) recombinant FVIII product indicator (1 = Refacto AF or other BDD recombinant FVIII; 0 = full-length recombinant or plasma-derived FVIII)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (full-length recombinant or plasma-derived FVIII)",
      notes              = "Hazendonk 2016 Methods: `C_pred,bdp = C_pred * (1 - theta_bdp)` with theta_bdp = 0.34, applied to the predicted plasma concentration to correct for a well-documented ~34% under-detection of B-domain-deleted FVIII by the one-stage clotting assay (Refacto AF was the only BDD product administered in this cohort). Cohort prevalence 14% of surgical procedures received a BDD product (Discussion: 'of which 14% were a B-domain deleted FVIII concentrate'). Per-observation (per-row) indicator; a subject may have received different products across multiple surgical procedures.",
      source_name        = "FORM_FVIII_BDD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 119L,
    n_studies      = 1L,
    n_surgeries    = 198L,
    age_range      = "0.2-78 years (adult 19-78, pediatric 0.2-17.3)",
    age_median     = "40 years total cohort (adult median 48; pediatric median 4.3)",
    weight_range   = "5-111 kg (adult 45-111, pediatric 5-85)",
    weight_median  = "75 kg total cohort (adult median 80; pediatric median 18.5)",
    sex_female_pct = 0,
    race_ethnicity = "not reported in source; Netherlands multi-center cohort",
    disease_state  = "Severe or moderate hemophilia A (FVIII plasma concentration < 0.05 IU/mL). 70% severe (< 0.01 IU/mL); 70% on prophylaxis. History of FVIII inhibitors recorded but inhibitor status at study entry not analysed as a covariate.",
    dose_range     = "Perioperative FVIII replacement: pre-operative bolus ~50 IU/kg followed by continuous or bolus infusion targeting FVIII plasma concentration 0.80-1.00 IU/mL for the first 24 h, 0.50-0.80 for 24-120 h, 0.30-0.50 for >120 h per National Hemophilia Consensus (Hazendonk 2016 Methods). Approximately 3-4 mL/kg/h continuous infusion rate; 58% of procedures used continuous infusion, 42% bolus.",
    regions        = "Netherlands (5 Academic Hemophilia Treatment Centers)",
    products       = "Recombinant FVIII (77% of procedures): Kogenate FS, Helixate FS, Advate, Recombinate, Refacto AF (the B-domain-deleted product, 14% of recombinant procedures); plasma-derived (23%): Aafact, Hemofil M",
    surgical_mix   = "Major or high-risk 49%; primary types orthopedic (47.5%), central-venous-catheter placement (16.2%, mostly pediatric), miscellaneous (15.7%), urology (6.1%), ENT (5.6%)",
    notes          = "Retrospective study of hemophilia A patients undergoing 198 elective, minor, or major surgical procedures between 2000 and 2013. 1389 total FVIII plasma concentration measurements (median ~7 per patient) collected as trough / peak / steady-state samples from pre-operative through post-day-6+. FVIII plasma concentrations were measured by one-stage clotting assay at all centers. The population is Netherlands-only; ethnicity / race not reported. Hemophilia A is X-linked recessive, so the cohort is essentially all-male (sex_female_pct = 0 assumed; the paper does not tabulate sex explicitly). Center-specific residual error was reported (0.15 IU/mL additive + 18% proportional for centers 1-3; 0.05 IU/mL additive + 23% proportional for centers 4-5); only the majority center-cluster values are encoded here (see vignette Assumptions and deviations)."
  )

  ini({
    # Structural PK parameters at the reference body weight of 68 kg
    # (total-cohort median). Values are the Hazendonk 2016 Table 4 "Final model"
    # column, reported in mL and mL/h; convert to L and L/h for internal
    # consistency with other FVIII models in nlmixr2lib.
    lcl <- log(0.150); label("Clearance CL at the reference 68 kg, 40-year non-O minor-surgery patient (L/h)")             # Hazendonk 2016 Table 4: CL = 150 mL/h/68 kg (RSE 8%)
    lvc <- log(2.810); label("Central volume V1 at the reference 68 kg, 40-year patient (L)")                              # Hazendonk 2016 Table 4: V1 = 2810 mL/68 kg (RSE 4%)
    lq  <- log(0.160); label("Inter-compartmental clearance Q at the reference 68 kg patient (L/h)")                        # Hazendonk 2016 Table 4: Q  = 160 mL/h/68 kg (RSE 20%)
    lvp <- log(1.900); label("Peripheral volume V2 at the reference 68 kg patient (L)")                                    # Hazendonk 2016 Table 4: V2 = 1900 mL/68 kg (RSE 11%)

    # B-domain-deleted product multiplicative correction on the predicted
    # concentration (paper equation: C_pred,bdp = C_pred * (1 - theta_bdp)).
    theta_bdp <- 0.34; label("Fractional under-prediction of plasma concentration for B-domain-deleted FVIII products (dimensionless)") # Hazendonk 2016 Table 4: theta B-domain deleted = 0.34 (RSE 13%)

    # Covariate coefficients (estimated). Values encode the paper's Table 5
    # equation form; see model() for how they combine with the covariates.
    e_age_cl   <- -0.17;  label("Power exponent of (AGE / 40) on CL (dimensionless)")                              # Hazendonk 2016 Table 4: theta5 CL-Age = -0.17 (RSE 22%)
    e_age_vc   <- -0.09;  label("Power exponent of (AGE / 40) on V1 (dimensionless)")                              # Hazendonk 2016 Table 4: theta8 V1-Age = -0.09 (RSE 28%)
    e_blood_cl <-  0.26;  label("Fractional effect of blood group O on CL (dimensionless; 1 + f * BLOOD_GROUP_O)") # Hazendonk 2016 Table 4: theta6 CL-Blood group O = 26% (RSE 7%)
    e_surg_cl  <- -0.07;  label("Fractional effect of major surgery on CL (dimensionless; 1 + f * SURG_SEV_MAJOR)") # Hazendonk 2016 Table 4: theta7 CL-Major surgery = -7% (RSE 6%)

    # Fixed allometric exponents (paper explicitly fixed these at theory-based
    # values; Hazendonk 2016 Structural model development and Table 3 footnote).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT / 68) on CL and Q (fixed)")  # Hazendonk 2016 Structural model development: "power exponents fixed at 0.75 for clearances"
    e_wt_vc <- fixed(1.00); label("Allometric exponent of (WT / 68) on V1 and V2 (fixed)") # Hazendonk 2016 Structural model development: "power exponents fixed at ... 1.0 for volumes of distribution"

    # Inter-individual variability (exponential): omega^2 = CV^2 following the
    # convention of related FVIII popPK models (Chelle 2019, Nestorov 2014).
    # IIV on Q and V2 was tested but not retained because shrinkage exceeded
    # 40% (Hazendonk 2016 Results, "Structural model development").
    etalcl ~ 0.1369  # Hazendonk 2016 Table 4: IIV CL = 37% (RSE 14%); omega^2 = 0.37^2 = 0.1369
    etalvc ~ 0.0729  # Hazendonk 2016 Table 4: IIV V1 = 27% (RSE 14%); omega^2 = 0.27^2 = 0.0729

    # Residual error (combined). Values are the Hazendonk 2016 Table 4 "Center
    # 1, 2, 3" row -- the majority center cluster (3 of 5 centers). See vignette
    # Assumptions and deviations for the alternative "Center 4, 5" residual
    # error (add 0.05 IU/mL + prop 23%).
    addSd  <- 0.15; label("Additive residual error at centers 1, 2, 3 (IU/mL)")                          # Hazendonk 2016 Table 4: Additive residual error Center 1,2,3 = 0.15 IU/mL (RSE 12%)
    propSd <- 0.18; label("Proportional residual error at centers 1, 2, 3 (fraction of predicted)")     # Hazendonk 2016 Table 4: Proportional residual error Center 1,2,3 = 0.18 (RSE 15%)
  })
  model({
    # Covariate multipliers on CL (paper Table 5):
    #   CL = 150 * (WT/68)^0.75 * (AGE/40)^-0.17 * 1.26^blood_group * 0.93^severity
    # The two categorical multipliers 1.26^b and 0.93^s are equivalent to
    # (1 + 0.26 * b) and (1 + (-0.07) * s) for binary indicators; the linear
    # form is used here for parameter identifiability and interpretability.
    age_eff_cl   <- (AGE / 40) ^ e_age_cl
    age_eff_vc   <- (AGE / 40) ^ e_age_vc
    blood_eff_cl <- 1 + e_blood_cl * BLOOD_GROUP_O
    surg_eff_cl  <- 1 + e_surg_cl  * SURG_SEV_MAJOR

    # Allometric scaling on body weight (reference 68 kg).
    ws <- WT / 68

    # Individual PK parameters (Hazendonk 2016 Table 5). CL, Q scale with WT^0.75;
    # V1, V2 scale with WT^1.0. Only CL and V1 carry IIV etas.
    cl <- exp(lcl + etalcl) * ws ^ e_wt_cl * age_eff_cl * blood_eff_cl * surg_eff_cl
    vc <- exp(lvc + etalvc) * ws ^ e_wt_vc * age_eff_vc
    q  <- exp(lq)           * ws ^ e_wt_cl
    vp <- exp(lvp)          * ws ^ e_wt_vc

    # Micro-constants for the explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV model. Doses (bolus or continuous infusion) enter the
    # central compartment directly; the event table supplies rate / duration.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation: FVIII activity in IU/mL. Dose in IU, vc in L give
    # central / vc in IU/L; divide by 1000 to reach IU/mL, the paper's unit.
    # (1 - theta_bdp * FORM_FVIII_BDD) applies the B-domain-deleted-product
    # under-detection correction; equals 1 for full-length / plasma-derived
    # products and (1 - 0.34) = 0.66 for BDD products (Hazendonk 2016 Methods).
    Cc <- (central / vc / 1000) * (1 - theta_bdp * FORM_FVIII_BDD)
    Cc ~ add(addSd) + prop(propSd)
  })
}

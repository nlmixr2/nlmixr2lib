Rozman_2017_rituximab <- function() {
  description <- "Two-compartment population PK model of rituximab in adults with diffuse large B-cell lymphoma (DLBCL) receiving R-CHOP; total CL is the sum of a time-stationary non-specific (IgG-catabolic) component cl_ss and a mono-exponentially decaying target-mediated component cl_time * exp(-kdes * time), with age and body weight on cl_ss, sex on V1 (central volume), and a post-hoc progression-free-survival event indicator (PFS_EVENT) on kdes (Rozman 2017)."
  reference <- "Rozman S, Grabnar I, Novakovic S, Mrhar A, Jezersek Novakovic B. Population pharmacokinetics of rituximab in patients with diffuse large B-cell lymphoma and association with clinical outcome. Br J Clin Pharmacol. 2017;83(8):1782-1790. doi:10.1111/bcp.13271"
  vignette <- "Rozman_2017_rituximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on the time-stationary clearance cl_ss with estimated exponent 1.23 (95% CI 0.70-1.73, which brackets the theoretical allometric 0.75 per Rozman 2017 Results). Reference weight 70 kg matches the typical-male subject Rozman 2017 used for the Figure 3 progressor / non-progressor simulation.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at treatment start",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Linear effect on the time-stationary clearance cl_ss, centered at 60 years: cl_ss is reduced by 0.82% per year above the 60-year reference (Rozman 2017 Results, Table 2 theta = -0.00820).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Rozman 2017 anchors the typical-value central volume V1 = 4.62 L to the male reference; women have 21.4% lower V1 (Rozman 2017 Results, Table 2 theta = -0.214). Applied via (1 + e_sexf_vc * SEXF) = (1 - 0.214) for females, (1) for males.",
      source_name        = "SEX"
    ),
    PFS_EVENT = list(
      description        = "Progression-free-survival event indicator (1 = patient experienced disease progression during follow-up, 0 = no progression observed through end of follow-up)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no progression event observed during follow-up)",
      notes              = "Post-hoc outcome-derived covariate: Rozman 2017 assessed disease progression during follow-up (median observation 52.9 months) per revised response criteria for malignant lymphoma (International Harmonization Project). 6 / 29 patients (20.7%) experienced progression. Applied as a proportional multiplier on the time-varying CL decay-rate kdes: kdes in progressors is 82.2% lower than in non-progressors (0.0254/day vs 0.143/day; Rozman 2017 Table 2 theta = -0.822, 95% CI -0.950 to -0.334). Biological interpretation (paper Discussion): slow decay of target-mediated CL reflects sustained CD20 target burden, which is associated with poorer response. Using an outcome variable as a covariate on PK carries selection-bias implications documented in the vignette Assumptions and deviations section.",
      source_name        = "progression status (paper does not disclose the NONMEM $INPUT column name)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 29L,
    n_studies      = 1L,
    n_observations = 512L,
    age_range      = "48-84 years",
    age_median     = "62 years",
    weight_range   = "54-100 kg",
    weight_median  = "74 kg",
    sex_female_pct = 44.8,
    race_ethnicity = "Caucasian (100%)",
    disease_state  = "Newly diagnosed diffuse large B-cell lymphoma (DLBCL); Ann Arbour clinical stage I-II 37.9% / III-IV 62.1%; International Prognostic Index (IPI) 0-2 65.5% / IPI 3-5 34.5%; bulky disease 34.5%.",
    dose_range     = "Rituximab 375 mg/m^2 IV every 3 weeks for 8 cycles of R-CHOP (rituximab + cyclophosphamide + doxorubicin + vincristine + methylprednisolone). Median dose per cycle 700 mg (range 500-800 mg). First-cycle slow infusion 50 mg/h ramped to a maximum 400 mg/h; subsequent cycles fast infusion (20% of dose over 30 min then 80% over 60 min).",
    regimen        = "Q3W IV infusion x 8 cycles. Sixteen of 29 patients received additional radiotherapy after R-CHOP.",
    regions        = "Slovenia (single-centre Institute of Oncology Ljubljana study).",
    followup       = "Median observation 52.9 months (range 9.7-66.3 months). Six patients (20.7%) experienced disease progression during follow-up (three with IPI = 4, three with IPI = 1); four progressed at 3 months post-therapy, one at 4 months, one at 5 months. Projected 5-year progression-free and overall survival both 79%.",
    sampling       = "18 samples per patient (16-18 available). Cycles 1-7: peak (15 min to 3 h post-infusion) and trough (immediately before the next cycle's infusion) per cycle. Cycle 8: peak plus 1, 3, and 6 months post-final-infusion.",
    assay          = "Enzyme-linked immunosorbent assay (ELISA) with rat anti-rituximab IgG2a capture (MB2A4) and goat peroxidase-conjugated anti-human IgG detection; five-parameter logistic calibration; between-run precision <= 13.8% CV, within-run precision <= 9.8% CV, accuracy <= 13.7% deviation. Calibration range 10-2000 mg/L (samples diluted 1/20000).",
    notes          = "First prospective popPK study of rituximab in newly diagnosed DLBCL patients on standard 3-week R-CHOP. Covariates tested and not retained in the final model: IPI, irradiation-following-R-CHOP indicator, and treatment response category (complete / partial / stable / progressive; the binary disease-progression indicator on KD was retained instead). Comorbidities (arterial hypertension, benign prostatic hyperplasia, hypercholesterolaemia) and concomitant medications (ACE inhibitors / ARB antagonists, proton-pump inhibitors, NSAIDs) were not expected to affect rituximab PK and were not tested."
  )

  ini({
    # ---- Structural PK parameters (Rozman 2017 Table 2, final covariate model) ----
    # Reference subject: male, 70 kg, 60 years, no progression during follow-up
    # (matches the 'typical male patient' Figure 3 simulation baseline).
    lcl_ss   <- log(0.252); label("Time-stationary non-specific clearance CL1 (L/day)")             # Rozman 2017 Table 2 CL1 = 0.252 L/day (95% CI 0.227-0.279)
    lcl_time <- log(0.278); label("Time-varying target-mediated clearance at t=0 CL2,0 (L/day)")    # Rozman 2017 Table 2 CL2,0 = 0.278 L/day (95% CI 0.181-0.390)
    lkdes    <- log(0.143); label("Decay rate constant of the time-varying CL arm KD (1/day)")      # Rozman 2017 Table 2 KD = 0.143 /day (95% CI 0.0478-0.418) at PFS_EVENT = 0 reference
    lvc      <- log(4.62);  label("Central volume of distribution V1 (L) for males")                 # Rozman 2017 Table 2 V1 = 4.62 L (95% CI 4.34-4.93) at male reference
    lvp      <- log(8.61);  label("Peripheral volume of distribution V2 (L)")                        # Rozman 2017 Table 2 V2 = 8.61 L (95% CI 7.45-9.81)
    lq       <- log(1.02);  label("Inter-compartmental clearance Q (L/day)")                         # Rozman 2017 Table 2 Q = 1.02 L/day (95% CI 0.664-1.95)

    # ---- Covariate effects (Rozman 2017 Table 2) ----
    # AGE on CL1 -- linear centered at 60 years; -0.82% per year above 60.
    e_age_cl_ss    <- -0.00820; label("Linear age effect on CL1 (per year above 60; cl_ss = cl_ss_typ * (1 + e_age_cl_ss * (AGE - 60)))")  # Rozman 2017 Table 2 age effect on CL1
    # WT on CL1 -- power / allometric-like, reference 70 kg.
    e_wt_cl_ss     <- 1.23;     label("Power exponent of (WT/70) on CL1 (unitless; cl_ss = cl_ss_typ * (WT/70)^e_wt_cl_ss)")               # Rozman 2017 Table 2 weight effect on CL1 (95% CI 0.70-1.73, brackets allometric 0.75)
    # SEXF on V1 -- proportional multiplier for females; male reference.
    e_sexf_vc      <- -0.214;   label("Proportional effect of SEXF on V1 (fraction; vc = vc_typ_male * (1 + e_sexf_vc * SEXF))")           # Rozman 2017 Table 2 sex effect on V1 (women 21.4% lower V1)
    # PFS_EVENT on KD -- proportional multiplier for progressors; no-progression reference.
    e_pfs_event_kdes <- -0.822; label("Proportional effect of PFS_EVENT on KD (fraction; kdes = kdes_typ * (1 + e_pfs_event_kdes * PFS_EVENT))")  # Rozman 2017 Table 2 disease-progression effect on KD (progressors 82.2% lower; 95% CI -0.950 to -0.334)

    # ---- Inter-individual variability (Rozman 2017 Table 2, IIV as CV%) ----
    # Log-normal variance transform: omega^2 = log(1 + CV^2).
    etalcl_ss ~ 0.03365  # CV = 18.5% -> log(1 + 0.185^2) = 0.03365; Rozman 2017 Table 2 IIV CL1
    etalkdes  ~ 1.27874  # CV = 161%  -> log(1 + 1.61^2)  = 1.27874; Rozman 2017 Table 2 IIV KD (very large; 22.7% shrinkage)
    etalvc    ~ 0.01337  # CV = 11.6% -> log(1 + 0.116^2) = 0.01337; Rozman 2017 Table 2 IIV V1

    # ---- Residual error (Rozman 2017 Table 2, combined additive + proportional) ----
    addSd  <- 2.46;  label("Additive residual error (mg/L)")           # Rozman 2017 Table 2 additive = 2.46 mg/L (95% CI 1.05-4.36)
    propSd <- 0.159; label("Proportional residual error (fraction)")   # Rozman 2017 Table 2 proportional = 15.9% (95% CI 14.3-17.5)
  })

  model({
    # ---- Individual structural parameters ---------------------------------------
    # Reference subject: male (SEXF = 0), 70 kg, 60 years, no follow-up progression (PFS_EVENT = 0).
    # Age enters as a linear proportional shift centred at 60 years (Rozman 2017 Methods
    # "Covariate model", Results paragraph reporting "0.82% decrease per year above 60 years").
    # Weight enters as a power / allometric-like factor with reference 70 kg (Equation 3 in
    # Rozman 2017; the 95% CI on the exponent brackets the theoretical allometric 0.75).
    sf_wt_cl_ss  <- (WT / 70)^e_wt_cl_ss
    sf_age_cl_ss <- 1 + e_age_cl_ss * (AGE - 60)
    cl_ss   <- exp(lcl_ss   + etalcl_ss) * sf_wt_cl_ss * sf_age_cl_ss
    cl_time <- exp(lcl_time)
    kdes    <- exp(lkdes    + etalkdes)  * (1 + e_pfs_event_kdes * PFS_EVENT)
    vc      <- exp(lvc      + etalvc)    * (1 + e_sexf_vc * SEXF)
    vp      <- exp(lvp)
    q       <- exp(lq)

    # ---- Time-varying total clearance --------------------------------------------
    # Rozman 2017 Methods "Structural model development", Equation (1):
    #   CL2(t) = CL2,0 * exp(-KD * t)
    # Total clearance from the central compartment is the sum of the time-stationary
    # non-specific arm (CL1, IgG catabolism) and the target-mediated time-varying arm
    # (CL2). The paper reports (Results): "Rituximab elimination is a sum of linear
    # (CL1) and time varying (CL2) clearance." The integration variable `time` starts
    # at the first-dose event (matching NONMEM TIME).
    cl <- cl_ss + cl_time * exp(-kdes * time)

    # ---- Micro-constants ---------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system --------------------------------------------------------------
    # Two-compartment linear disposition; IV infusion into the central compartment.
    # kel(t) = cl(t) / vc is time-varying via the exp(-kdes * time) term above.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- Observation and residual error -----------------------------------------
    # Dose in mg, volume in L -> Cc in mg/L (numerically equivalent to ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

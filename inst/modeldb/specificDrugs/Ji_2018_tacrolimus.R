Ji_2018_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model for oral tacrolimus in ",
    "Korean adult living-donor liver-transplant recipients during the first 14 ",
    "days post-transplantation (Ji 2018). First-order absorption with ka fixed ",
    "at 4.48 1/h from prior reports; CL/F = 6.33 * POD^0.257 multiplied by a ",
    "combinational CYP3A5 recipient-and-donor categorical factor (2.314 if both ",
    "recipient and donor are CYP3A5 expressers; 1.523 if the recipient is a ",
    "CYP3A5 expresser and the donor is a nonexpresser; 1.0 otherwise); V/F = ",
    "465 * POD^0.322; exponential IIV on CL/F and V/F; combined proportional + ",
    "additive residual error on whole-blood tacrolimus concentration."
  )
  reference <- paste0(
    "Ji E, Kim MG, Oh JM. CYP3A5 genotype-based model to predict tacrolimus ",
    "dosage in the early postoperative period after living donor liver ",
    "transplantation. Ther Clin Risk Manag. 2018;14:2119-2126. ",
    "doi:10.2147/TCRM.S184376"
  )
  vignette <- "Ji_2018_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    POD = list(
      description        = "Days post-transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-varying within subject. Enters CL/F and V/F as a direct power ",
        "covariate (POD^0.257 on CL/F, POD^0.322 on V/F; Ji 2018 final-model ",
        "equation in Results), not as a centred-deviation effect. The model ",
        "was developed using daily trough samples on POD 1 through POD 14, ",
        "so POD is expected to take integer-valued days >= 1 over that range. ",
        "Because POD enters as a power, values of POD = 0 collapse CL/F and ",
        "V/F to zero; downstream users simulating predose / day-of-surgery ",
        "events should restrict simulations to POD >= 1 or clamp POD = max(POD, 1)."
      ),
      source_name        = "POD"
    ),
    CYP3A5_EXPR = list(
      description        = "Recipient CYP3A5 expresser status (rs776746 / CYP3A5*3 polymorphism)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 non-expresser recipient)",
      notes              = paste0(
        "Time-fixed germline genotype of the graft recipient. 1 = at least one ",
        "functional CYP3A5*1 allele (genotype *1/*1 or *1/*3); 0 = homozygous ",
        "*3/*3. Derived from the recipient half of the combinational CYP3A5 ",
        "group label used by Ji 2018 (REDE/REDN: recipient expresser; ",
        "RNDE/RNDN: recipient nonexpresser). The combinational effect on CL/F ",
        "depends on both CYP3A5_EXPR and CYP3A5_EXPR_DONOR: the donor-side ",
        "genotype only modifies CL/F when the recipient is itself an ",
        "expresser. Cohort distribution: 23/58 (39.7%) recipient expressers ",
        "(REDE n=10, REDN n=13), 35/58 (60.3%) recipient nonexpressers (RNDE ",
        "n=8, RNDN n=27); Hardy-Weinberg equilibrium (Ji 2018 Methods)."
      ),
      source_name        = "CYP3A5 recipient (REDE/REDN -> 1; RNDE/RNDN -> 0)"
    ),
    CYP3A5_EXPR_DONOR = list(
      description        = "Donor CYP3A5 expresser status (rs776746 / CYP3A5*3 polymorphism)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 non-expresser donor)",
      notes              = paste0(
        "Time-fixed germline genotype of the liver-graft donor. 1 = donor ",
        "carries at least one functional CYP3A5*1 allele (genotype *1/*1 or ",
        "*1/*3); 0 = donor is homozygous *3/*3. Derived from the donor half ",
        "of the combinational CYP3A5 group label used by Ji 2018 (REDE/RNDE: ",
        "donor expresser; REDN/RNDN: donor nonexpresser). The combinational ",
        "effect on CL/F is applied only when the recipient is an expresser ",
        "(CYP3A5_EXPR = 1): when recipient and donor are both expressers ",
        "(REDE), CL/F is multiplied by 2.314; when the recipient is an ",
        "expresser and the donor is not (REDN), CL/F is multiplied by 1.523; ",
        "when the recipient is a nonexpresser (RNDE or RNDN) the donor's ",
        "genotype has no effect on CL/F and the factor is 1. Cohort ",
        "distribution: 18/58 (31.0%) donor expressers (REDE n=10, RNDE n=8), ",
        "40/58 (69.0%) donor nonexpressers (REDN n=13, RNDN n=27); ",
        "Hardy-Weinberg equilibrium (Ji 2018 Methods)."
      ),
      source_name        = "CYP3A5 donor (REDE/RNDE -> 1; REDN/RNDN -> 0)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 58L,
    n_studies        = 1L,
    n_observations   = 605L,
    age_range        = "19-65 years",
    age_median       = "49.2 years (mean +/- SD 49.2 +/- 8.7)",
    weight_range     = "40.1-85.5 kg",
    weight_median    = "61.4 kg (mean +/- SD 61.4 +/- 10.1)",
    sex_female_pct   = 20.7,
    race_ethnicity   = "Korean (single-centre Seoul, Republic of Korea)",
    disease_state    = paste0(
      "Korean adult patients receiving de novo living-donor liver ",
      "transplantation (LDLT). Patients were on a triple- or double-drug ",
      "immunosuppression regimen including tacrolimus and corticosteroids, ",
      "with or without mycophenolate mofetil."
    ),
    dose_range       = paste0(
      "Oral tacrolimus (Astellas Prograf), starting on postoperative day 1, ",
      "administered twice daily at 10:00 and 22:00 on an empty stomach. ",
      "Per-dose range 0.1-6 mg (mean +/- SD 1.9 +/- 1.2 mg). Doses were ",
      "adjusted by therapeutic drug monitoring to achieve target whole-blood ",
      "trough concentrations of 8-13 ng/mL (triple regimen) or 13-17 ng/mL ",
      "(double regimen) over the first 14 days post-transplant."
    ),
    regions          = "Republic of Korea (Seoul National University Hospital, single-centre).",
    cyp3a5_distribution = paste0(
      "Four combinational recipient/donor groups: REDE (recipient expresser ",
      "+ donor expresser) n = 10; REDN (recipient expresser + donor ",
      "nonexpresser) n = 13; RNDE (recipient nonexpresser + donor expresser) ",
      "n = 8; RNDN (recipient nonexpresser + donor nonexpresser) n = 27. ",
      "Hardy-Weinberg equilibrium for both recipient and donor genotype ",
      "distributions (all P > 0.05; Ji 2018 Methods)."
    ),
    pod_range        = "1-14 days (first 14 days post-transplant).",
    sampling_design  = paste0(
      "Routine therapeutic drug monitoring trough samples drawn daily around ",
      "09:00 (~11 hours after the prior 22:00 dose), starting the day after ",
      "the first dose and continuing until discharge. Total 605 trough ",
      "concentrations (median 10.4 samples per patient, range 3-13). ",
      "Whole-blood tacrolimus measured by enzyme immunoassay. Observed ",
      "concentration mean +/- SD 9.7 +/- 3.9 ng/mL (range 1.6-21.4)."
    ),
    notes            = paste0(
      "Single-centre retrospective analysis of 58 patients from a prior ",
      "study (Ji 2018 reference 3). Population pharmacokinetic modelling ",
      "performed in NONMEM v7.4 + PsN with FOCE-I; internal validation by ",
      "bootstrap (500 resamples) and visual predictive check. Baseline ",
      "demographics from Ji 2018 Table 1; covariate-stratified clearances ",
      "in Figure 1. Other patient baseline values (mean +/- SD, range): ",
      "graft-to-recipient body weight ratio 1.07 +/- 0.24 (0.59-1.60); ",
      "hematocrit 26.7 +/- 4.5%% (15.8-40.9); total bilirubin 2.4 +/- 1.9 ",
      "mg/dL (0.4-12.0); ALT 141 +/- 113 IU/L (7-570); albumin 3.0 +/- 0.3 ",
      "g/dL (1.7-3.9). All of these were evaluated as candidate covariates ",
      "during stepwise forward/backward modelling but only POD and the ",
      "combinational CYP3A5 group were retained."
    )
  )

  ini({
    # ----- Structural PK (Ji 2018 Table 2 final-model column) -----
    # One-compartment oral disposition. Reference values: typical subject at
    # POD = 1 with a recipient CYP3A5 nonexpresser (i.e. RNDE or RNDN; the
    # combinational genotype factor is 1). ka was held fixed at 4.48 1/h
    # because the observed concentrations were all troughs and the absorption
    # phase was not identifiable from the data; the literature value was
    # taken from Ji 2018 references 14 and 15 (Methods).
    lka <- fixed(log(4.48)) ; label("Absorption rate constant ka (1/h)")                                                       # Ji 2018 Methods + Table 2: ka fixed at 4.48 1/h (references 14, 15)
    lcl <- log(6.33)        ; label("Apparent oral clearance CL/F at POD = 1, recipient CYP3A5 nonexpresser (L/h)")            # Ji 2018 Table 2 final CL/F = 6.33 L/h (RSE 16%, bootstrap 4.48-8.06)
    lvc <- log(465)         ; label("Apparent volume of distribution V/F at POD = 1 (L)")                                       # Ji 2018 Table 2 final V/F = 465 L (RSE 11%, bootstrap 384-560)

    # ----- Covariate effects -----
    # POD enters CL/F and V/F as a direct power covariate POD^e (Ji 2018
    # Results final-model equation):
    #   CL/F = 6.33 * POD^0.257 * CYP3A5_factor
    #   V/F  = 465  * POD^0.322
    # The CYP3A5 combinational categorical factor takes one of three values:
    # 2.314 for REDE, 1.523 for REDN, 1.0 for RNDE / RNDN (reference; the
    # paper merged RNDE and RNDN because their estimated effects were similar
    # and the combination had a lower AIC).
    e_pod_cl  <- 0.257 ; label("Power exponent of POD on CL/F (unitless)")  # Ji 2018 Table 2 final POD on CL/F = 0.257 (RSE 30%)
    e_pod_vc  <- 0.322 ; label("Power exponent of POD on V/F (unitless)")   # Ji 2018 Table 2 final POD on V/F  = 0.322 (RSE 16%)
    e_rede_cl <- 2.314 ; label("Multiplicative CL/F factor for REDE (recipient + donor both CYP3A5 expressers; vs RNDE/RNDN reference)")  # Ji 2018 Table 2 final CYP3A5 REDE = 2.314 (RSE 12%, bootstrap 1.907-2.742)
    e_redn_cl <- 1.523 ; label("Multiplicative CL/F factor for REDN (recipient CYP3A5 expresser + donor nonexpresser; vs RNDE/RNDN reference)")  # Ji 2018 Table 2 final CYP3A5 REDN = 1.523 (RSE 32%, bootstrap 1.201-1.923)

    # ----- Inter-individual variability (Ji 2018 Table 2) -----
    # Exponential IIV. Ji 2018 reports omega^2 (variance on the log-scale eta)
    # together with the approximation IIV (%) = sqrt(omega^2) * 100 in the
    # Table 2 footnote: 0.117 -> 34.2%, 0.207 -> 45.5%. No CL/F-V/F covariance
    # is reported, so the two etas are treated as independent.
    etalcl ~ 0.117  # Ji 2018 Table 2 omega^2 CL/F = 0.117 (IIV 34.2%; RSE 23%, bootstrap 0.071-0.160)
    etalvc ~ 0.207  # Ji 2018 Table 2 omega^2 V/F  = 0.207 (IIV 45.5%; RSE 23%, bootstrap 0.139-0.289)

    # ----- Residual error (Ji 2018 Table 2) -----
    # Combined proportional + additive on the linear whole-blood
    # concentration scale (ng/mL). Ji 2018 reports the variances
    # sigma^2_prop = 0.182 and sigma^2_add = 0.838; nlmixr2 takes standard
    # deviations, so propSd = sqrt(0.182) = 0.4266 and addSd = sqrt(0.838)
    # = 0.9154 ng/mL.
    propSd <- 0.4266 ; label("Proportional residual error (fraction)")             # Ji 2018 Table 2 sigma^2 proportional = 0.182 (paper-reported 42.7%; sqrt(0.182) = 0.4266; RSE 9%, bootstrap 0.151-0.203)
    addSd  <- 0.9154 ; label("Additive residual error on whole-blood Cc (ng/mL)") # Ji 2018 Table 2 sigma^2 additive = 0.838 ng^2/mL^2; sqrt(0.838) = 0.9154 ng/mL (RSE 22%, bootstrap 0.509-1.156)
  })

  model({
    # ----- 1. Derived covariate terms -----
    # The Ji 2018 combinational CYP3A5 categorical effect on CL/F applies
    # different factors to the REDE and REDN combinations and treats the
    # remaining cells (RNDE, RNDN) as the unit reference. The two binary
    # canonical inputs (CYP3A5_EXPR for the recipient, CYP3A5_EXPR_DONOR for
    # the donor) compose to give mutually exclusive indicators for REDE and
    # REDN; both are zero when the recipient is a nonexpresser (RNDE / RNDN
    # reference). The factor on CL/F is then
    #   e_rede_cl^is_rede * e_redn_cl^is_redn,
    # which evaluates to 2.314 / 1.523 / 1.0 in the three identified groups.
    is_rede <- CYP3A5_EXPR * CYP3A5_EXPR_DONOR
    is_redn <- CYP3A5_EXPR * (1 - CYP3A5_EXPR_DONOR)

    # ----- 2. Individual PK parameters -----
    # POD enters as a direct power covariate (Ji 2018 final-model equation in
    # Results: CL/F = 6.33 * POD^0.257 * ...; V/F = 465 * POD^0.322).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * POD^e_pod_cl *
          e_rede_cl^is_rede * e_redn_cl^is_redn
    vc <- exp(lvc + etalvc) * POD^e_pod_vc

    # ----- 3. Micro-rate constants -----
    kel <- cl / vc

    # ----- 4. ODE system (one-compartment oral) -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 5. Observation and error -----
    # Dose in mg, central amount in mg, vc in L -> mg/L; multiply by 1000 to
    # convert to ng/mL, the unit used throughout Ji 2018 (Table 1: mean
    # observed tacrolimus concentration 9.7 +/- 3.9 ng/mL; Table 2 reports
    # the additive residual variance on that scale).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

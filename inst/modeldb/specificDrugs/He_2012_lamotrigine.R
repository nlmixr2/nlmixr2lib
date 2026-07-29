He_2012_lamotrigine <- function() {
  description <- "One-compartment population PK model for oral lamotrigine in Chinese paediatric patients with epilepsy aged 0.5-17 years (He 2012). First-order absorption with Ka fixed at 1.0 1/h and bioavailability fixed at 1 (lamotrigine steady-state trough therapeutic-drug-monitoring data, which do not identify Ka or F), and first-order elimination from a single central compartment. Apparent oral clearance is scaled by an estimated power of total body weight (exponent 0.635) and modified exponentially by concomitant antiepileptic comedication: valproate (CONMED_VPA) reduces CL, while the enzyme-inducers carbamazepine (CONMED_CBZ) and phenobarbital (CONMED_PB) increase CL. Apparent central volume is fixed at 16.7 L at the 27.87 kg reference weight, scaled linearly with total body weight (allometric exponent fixed at 1.0)."
  reference   <- "He DK, Wang L, Lu W, Qin J, Zhang S, Li L, Zhang JM, Bao WQ, Song XQ, Liu HT. Population pharmacokinetics of lamotrigine in Chinese children with epilepsy. Acta Pharmacol Sin. 2012 Nov;33(11):1417-1423. doi:10.1038/aps.2012.118"
  vignette    <- "He_2012_lamotrigine"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (paper notation TBW).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate used for allometric scaling of both apparent oral clearance (estimated power exponent 0.635) and apparent central volume (linear, exponent fixed at 1.0). Reference weight 27.87 kg is the mean total body weight of the PPK model group (n=116; He 2012 Table 1).",
      source_name        = "TBW"
    ),
    CONMED_VPA = list(
      description        = "Concomitant valproate (VPA) coadministration indicator: 1 = patient is taking valproate as a concomitant antiepileptic drug, 0 = not taking valproate.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant valproate)",
      notes              = "Exponential multiplicative effect on apparent oral clearance: cl *= exp(-0.753 * CONMED_VPA), corresponding to a 53% reduction in CL relative to no-VPA reference (factor 0.471). VPA inhibits UGT enzymes that mediate lamotrigine glucuronidation (He 2012 Discussion paragraph 4; Table 3 theta3).",
      source_name        = "VPA"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine (CBZ) coadministration indicator: 1 = patient is taking carbamazepine, 0 = not taking carbamazepine.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = "Exponential multiplicative effect on apparent oral clearance: cl *= exp(0.868 * CONMED_CBZ), corresponding to a ~2.4-fold increase in CL relative to no-CBZ reference. CBZ is an inducer of hepatic UGT and CYP enzymes responsible for lamotrigine elimination (He 2012 Discussion paragraph 3; Table 3 theta4).",
      source_name        = "CBZ"
    ),
    CONMED_PB = list(
      description        = "Concomitant phenobarbital (PB) coadministration indicator: 1 = patient is taking phenobarbital, 0 = not taking phenobarbital.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes              = "Exponential multiplicative effect on apparent oral clearance: cl *= exp(0.633 * CONMED_PB), corresponding to a ~1.9-fold increase in CL relative to no-PB reference. PB is a broad-spectrum CYP/UGT inducer (He 2012 Discussion paragraph 3; Table 3 theta5).",
      source_name        = "PB"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 116,
    n_studies       = 1,
    age_range       = "0.5-17 years",
    age_mean        = "6.91 years (PPK model group)",
    weight_range    = "8-85 kg (PPK model group)",
    weight_mean     = "27.87 kg (PPK model group)",
    sex_female_pct  = round(48 / 116 * 100, 1),
    race_ethnicity  = "Chinese (Han majority; not formally tested as covariate).",
    disease_state   = "Paediatric epilepsy (diagnosis by clinician based on seizures and electroencephalogram; LTG therapy duration > 1 month so concentrations reflect steady-state trough; normal liver and kidney function).",
    dose_range      = "Oral lamotrigine 12.5-525 mg/day in the PPK model group (mean 135 mg/day); reported as steady-state trough concentrations (no dosing-history detail in the publication).",
    regions         = "China (multi-province referrals to the Department of Paediatrics, Peking University First Hospital, Beijing).",
    n_observations  = "191 steady-state trough lamotrigine plasma concentrations across 116 patients in the PPK model group (1-10 samples per patient, mean 1.65); a separate PPK valid group of 168 patients with 213 concentrations was used for external validation.",
    co_medication   = "Concomitant AEDs analysed: valproate (VPA, 63.4% of patients), carbamazepine (CBZ, 27.2%), phenobarbital (PB, 3.7%), oxcarbazepine (OXC, 7.9%), clonazepam (CZP, 8.9%), levetiracetam (LEV, 3.1%), topiramate (TPM, 5.2%). Only VPA, CBZ, and PB were retained in the final CL model; OXC, CZP, LEV, TPM had non-significant effects in univariate screening (He 2012 Table 2 models 8-11; SE > 100%).",
    notes           = "Sparse therapeutic-drug-monitoring data: blood drawn before breakfast and before the morning dose at steady state, i.e. trough concentrations only. The model was developed in NONMEM V (ADVAN2 TRANS2) using a one-compartment open-kinetic model with first-order absorption and elimination. Because the trough data do not inform Ka or F, both were fixed (Ka = 1.0 1/h, F = 1.0) per the upstream reference paper cited as [7]. Age, sex, and concomitant OXC/CZP/LEV/TPM were tested but not retained in the final model (He 2012 Results paragraph 1; Table 2 univariate screening summarised in Table 2 rows 7-11)."
  )

  ini({
    # Structural parameters - final-model NONMEM estimates from He 2012
    # Table 3 ("The final model" column). All clearance and volume terms
    # are apparent (X/F) because lamotrigine was administered orally and
    # bioavailability was fixed at 1.0; absolute F is not separately
    # identifiable from steady-state trough data.

    lka <- fixed(log(1.0)); label("First-order absorption rate constant (Ka, 1/h, fixed)")          # He 2012 Methods, page 1418: "absorption rate (Ka) was fixed at 1.0 h-1"
    lcl <- log(1.01);       label("Apparent oral clearance at 27.87 kg reference TBW (CL/F, L/h)")  # He 2012 Table 3 final model: theta1 = 1.01 (RSE 4.48%)
    lvc <- fixed(log(16.7)); label("Apparent central volume at 27.87 kg reference TBW (V/F, L, fixed)") # He 2012 Table 3 final model: theta6 = 16.7 (RSE "-", i.e. fixed); derived by solving 45 = theta6 * (75/27.87) so the model reproduces 45 L at a 75 kg adult (He 2012 Results, page 1420)

    # Allometric exponents on total body weight. CL exponent estimated;
    # Vc exponent fixed at 1.0 (linear scaling per the paper's stated
    # Vd formula Vd = 16.7 * (TBW/27.87)).
    e_wt_cl <- 0.635;       label("Allometric TBW exponent on CL/F (unitless)")              # He 2012 Table 3 final model: theta2 = 0.635 (RSE 7.91%)
    e_wt_vc <- fixed(1.0);  label("Allometric TBW exponent on V/F (unitless, fixed)")        # He 2012 Results, page 1420: "It was intended to be a fixed allometric covariate model such as Vd(L) = 16.7 * (WT/27.87)"

    # AED-coadministration effects on apparent oral clearance.
    # Exponential categorical form: cl *= exp(theta * indicator) for each
    # concomitant AED, matching the NONMEM final regression model
    # CL = theta1 * (TBW/27.87)^theta2 * exp(theta3*VPA) * exp(theta4*CBZ) * exp(theta5*PB)
    # printed by He 2012 in the Results paragraph 2.
    e_conmed_vpa_cl <- -0.753; label("Effect of CONMED_VPA on CL/F (exponential)")           # He 2012 Table 3 final model: theta3 = -0.753 (RSE 6.28%)
    e_conmed_cbz_cl <-  0.868; label("Effect of CONMED_CBZ on CL/F (exponential)")           # He 2012 Table 3 final model: theta4 =  0.868 (RSE 5.76%)
    e_conmed_pb_cl  <-  0.633; label("Effect of CONMED_PB on CL/F (exponential)")            # He 2012 Table 3 final model: theta5 =  0.633 (RSE 13.1%)

    # Interindividual variability. He 2012 Methods (page 1418) describes
    # the IIV model as CL_j = CL_pop * exp(eta_j), where eta has mean
    # zero and variance omega^2. The Table 3 row omega^2_CL therefore
    # reports the variance on the internal log scale directly; the value
    # is pasted into nlmixr2's variance-scale ~ syntax as-is.
    # He 2012 Table 3 final model: omega^2_CL = 0.067 (RSE 25.0%);
    # reported CV = sqrt(exp(0.067) - 1) ~= 26%, matching the
    # 25.8% CL CV reported in the Results paragraph 3 narrative.
    etalcl ~ 0.067

    # Residual error. He 2012 Methods (page 1418) writes the residual
    # model as C_obs = C_pred * (1 + epsilon) with epsilon ~ N(0, sigma^2);
    # sigma^2 in Table 3 is the variance of epsilon on the linear scale,
    # so the proportional SD passed to nlmixr2 is sqrt(sigma^2).
    propSd <- sqrt(0.045); label("Proportional residual error (fraction)")                   # He 2012 Table 3 final model: sigma^2 = 0.045 (RSE 37.6%); propSd = sqrt(0.045) = 0.2121
  })

  model({
    # Reference TBW for allometric scaling: the mean TBW of the PPK model
    # group (He 2012 Table 1: 27.87 kg).
    ref_wt <- 27.87

    # AED-coadministration multiplier on CL/F. Exponential categorical
    # form per He 2012 Results paragraph 2 final regression equation.
    # Patients on combination therapy (e.g., LTG + VPA + CBZ) get the
    # product of indicator-specific exponentials; the paper's Table 5
    # "LTG + VPA + (PB or CBZ)" row reflects this multiplicative
    # cancellation (the inducer and VPA effects approximately balance).
    aed_cl <- exp(e_conmed_vpa_cl * CONMED_VPA) *
              exp(e_conmed_cbz_cl * CONMED_CBZ) *
              exp(e_conmed_pb_cl  * CONMED_PB)

    # Individual parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl * aed_cl
    vc <- exp(lvc) * (WT / ref_wt)^e_wt_vc

    # Micro-constant.
    kel <- cl / vc

    # ODE system: one-compartment with first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. Dose in mg, volume in L -> mg/L (numerically equal to
    # the paper's ug/mL units).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

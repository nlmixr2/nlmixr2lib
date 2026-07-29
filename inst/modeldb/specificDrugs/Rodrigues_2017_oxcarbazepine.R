Rodrigues_2017_oxcarbazepine <- function() {
  description <- "Parent-metabolite population PK model for oral oxcarbazepine (OXC) and its active monohydroxy derivative (MHD) in epileptic children aged 2-12 years (Rodrigues 2017). Two-compartment OXC + one-compartment MHD with first-order absorption, complete metabolic conversion (Fm fixed to 1), reversible MHD-to-OXC back-transformation (KBT), empirical allometric weight scaling on CL_OXC/F, Vc_OXC/F, CL_MHD/F, and Vc_MHD/F (no scaling on Q_OXC/F or Vp_OXC/F), and a 29.3% increase in MHD clearance under concomitant enzyme-inducing antiepileptic drugs."
  reference   <- "Rodrigues C, Chiron C, Rey E, Dulac O, Comets E, Pons G, Jullien V. Population pharmacokinetics of oxcarbazepine and its monohydroxy derivative in epileptic children. Br J Clin Pharmacol. 2017 Dec;83(12):2695-2708. doi:10.1111/bcp.13392"
  vignette    <- "Rodrigues_2017_oxcarbazepine"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; constant within an individual in the source dataset).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Empirical allometric exponents on CL_OXC/F (0.798), Vc_OXC/F (2.4), CL_MHD/F (0.549), and Vc_MHD/F (1.09); Q_OXC/F and Vp_OXC/F are not weight-scaled in the final empirical model (Rodrigues 2017 Table 3 and final-model equation block, page 2699). Reference weight 70 kg.",
      source_name        = "WT"
    ),
    CONMED_EIAED = list(
      description        = "Concomitant enzyme-inducing antiepileptic drug indicator: 1 = patient is taking carbamazepine, phenobarbital, or phenytoin; 0 = none of these EIAEDs.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no EIAED coadministration; in the source paper this is the perturbation arm via MED = 1 -> CL_MHD reduced 22.7%).",
      notes              = "Source paper uses MED with the inverted value convention (MED = 1 if EIAEDs are ABSENT, 0 if present). The canonical CONMED_EIAED column inverts this so 1 = on EIAED (matching the CONMED_* convention). The model() block applies the source coefficient via the absence indicator (1 - CONMED_EIAED), so CL_MHD = 4.11 L/h/70 kg with EIAEDs and 4.11 * exp(-0.257) = 3.18 L/h/70 kg without. Rodrigues 2017 Table 1 enumerates the per-subject AED comedication; vigabatrin, clobazam, valproic acid, clonazepam, lamotrigine, diazepam, ethosuccimide, and progabide are NOT counted as EIAEDs.",
      source_name        = "MED"
    )
  )

  population <- list(
    n_subjects     = 31,
    n_studies      = 1,
    age_range      = "2.25-12.5 years",
    age_median     = "8.08 years",
    weight_range   = "12.7-56 kg",
    weight_median  = "23 kg",
    sex_female_pct = 41.9,
    race_ethnicity = "Not reported (single-center French paediatric epilepsy cohort).",
    disease_state  = "Children with inadequately controlled partial-onset and/or generalised atonic, tonic, or tonic-clonic seizures (>= 1 seizure/week despite 1-3 background AEDs unchanged for >= 1 month). 24 of 31 patients were comedicated with at least one enzyme-inducing AED (carbamazepine, phenobarbital, or phenytoin).",
    dose_range     = "Single oral OXC dose of 5 or 15 mg/kg as oral suspension after an overnight fast (14 patients received 5 mg/kg, 17 patients received 15 mg/kg).",
    regions        = "France (Cochin, Saint-Vincent de Paul, and Saint-Anne hospitals).",
    n_observations = "277 OXC and 279 MHD plasma samples (sampling at baseline, ~1, 2, 4, 6, 8, 12, 24, 36, 48 h post-dose); LLOQ 0.05 mg/L for OXC and 0.10 mg/L for MHD. After keeping only the first BLQ per patient (M3 method), 13.7% of OXC and 6.5% of MHD observations remained BLQ.",
    co_medication  = "Background AEDs reported in Table 1: carbamazepine 61.3%, vigabatrin 45.2%, clobazam 25.8%, phenytoin 16.1%, valproic acid 12.9%, clonazepam 9.7%, lamotrigine 9.7%, diazepam 6.5%, phenobarbital 6.5%, ethosuccimide 3.2%, progabide 3.2%. Six patients on one AED, 19 on two AEDs, six on three AEDs.",
    notes          = "Patients with renal or hepatic failure, untreated hypothyroidism, congenital metabolic disease, or body weight outside +/- 2 SD of normal were excluded. The model is reported by the authors to apply only to 2-12-year-old patients within the inclusion-criteria weight range. Demographics from Rodrigues 2017 Patient characteristics paragraph and Table 1."
  )

  ini({
    # Structural parameters - final empirical-allometry model from Rodrigues 2017
    # equation block (page 2699) and Table 3 column 2 ("Model with estimated
    # allometric exponents"). All clearances and volumes are apparent (X/F);
    # bioavailability F was fixed to 1 per the cited Flesch 2011 absolute-F
    # estimate of 0.99 (Rodrigues 2017 Methods, page 2697).
    lka     <- log(1.83);  label("First-order absorption rate constant (1/h)")                                  # Rodrigues 2017 Table 3, empirical model: Ka = 1.83 h^-1, RSE 4%
    lcl     <- log(140);   label("OXC apparent clearance CL_OXC/F at 70 kg (L/h)")                              # Rodrigues 2017 Table 3 + final-model equation: CL_OXC/F = 140 L/h/70 kg, RSE 24%
    lvc     <- log(337);   label("OXC apparent central volume Vc_OXC/F at 70 kg (L)")                           # Rodrigues 2017 Table 3 + final-model equation: Vc_OXC/F = 337 L/70 kg, RSE 41%
    lq      <- log(62.5);  label("OXC apparent inter-compartmental clearance Q_OXC/F (L/h, not weight-scaled)") # Rodrigues 2017 Table 3 + final-model equation: Q_OXC/F = 62.5 L/h, RSE 21%; no allometric scaling in empirical model
    lvp     <- log(60.7);  label("OXC apparent peripheral volume Vp_OXC/F (L, not weight-scaled)")              # Rodrigues 2017 Table 3 + final-model equation: Vp_OXC/F = 60.7 L, RSE 25%; no allometric scaling in empirical model

    # Metabolite parameters (suffix "_mhd" for the 10-monohydroxy derivative;
    # registered in conventions.R::registeredMetabolites). Reference value of
    # CL_MHD/F = 4.11 L/h/70 kg corresponds to a child WITH EIAEDs; the
    # absence-of-EIAED multiplier exp(-0.257) reduces it to 3.18 L/h/70 kg
    # (Rodrigues 2017 Discussion, page 2701).
    lcl_mhd <- log(4.11);  label("MHD apparent clearance CL_MHD/F at 70 kg, EIAED+ reference (L/h)")            # Rodrigues 2017 Table 3 + final-model equation: CL_MHD/F = 4.11 L/h/70 kg with EIAEDs, RSE 14%
    lvc_mhd <- log(54.8);  label("MHD apparent volume Vc_MHD/F at 70 kg (L)")                                   # Rodrigues 2017 Table 3 + final-model equation: Vc_MHD/F = 54.8 L/70 kg, RSE 16%

    # Back-transformation rate constant (MHD -> OXC). First-order on the MHD
    # central compartment amount; mass-balanced so the same flux is added to
    # the OXC central compartment.
    lkbt    <- log(0.0622); label("MHD -> OXC back-transformation rate constant (1/h)")                         # Rodrigues 2017 Table 3 + final-model equation: KBT = 0.0622 h^-1, RSE 15%

    # Bioavailability anchor: F fixed to 1 per the source paper.
    lfdepot <- fixed(log(1)); label("OXC depot bioavailability (fixed to 1)")                                   # Rodrigues 2017 Methods, page 2697: "Based on previous results evidencing a bioavailability of OXC of 0.99, this parameter was fixed to 1"

    # Empirical allometric exponents on the four weight-scaled parameters.
    # Reported as estimated point values (with RSEs), not "FIX" -- these are
    # NOT wrapped in fixed().
    e_wt_cl     <- 0.798;  label("Allometric exponent on CL_OXC/F (unitless)")                                  # Rodrigues 2017 Table 3, empirical model: theta_WT_CLOXC = 0.798, RSE 26%
    e_wt_vc     <- 2.4;    label("Allometric exponent on Vc_OXC/F (unitless)")                                  # Rodrigues 2017 Table 3, empirical model: theta_WT_VcOXC = 2.4, RSE 17%
    e_wt_cl_mhd <- 0.549;  label("Allometric exponent on CL_MHD/F (unitless)")                                  # Rodrigues 2017 Table 3, empirical model: theta_WT_CLMHD = 0.549, RSE 21%
    e_wt_vc_mhd <- 1.09;   label("Allometric exponent on Vc_MHD/F (unitless)")                                  # Rodrigues 2017 Table 3, empirical model: theta_WT_VcMHD = 1.09, RSE 13%

    # EIAED comedication effect on CL_MHD/F. Source coefficient is on the
    # absence-of-EIAED indicator (MED in the paper), so the model() block
    # applies it via (1 - CONMED_EIAED). Negative sign reduces CL_MHD when
    # EIAEDs are absent (Rodrigues 2017 Discussion: "EIAEDs increased MHD
    # clearance by 29.3%", i.e. 4.11/3.18 = 1.293).
    e_eiaed_cl_mhd <- -0.257; label("Exponential coefficient on absence-of-EIAED for CL_MHD/F (unitless)")      # Rodrigues 2017 Table 3, empirical model: theta_nEIAEDs_CLMHD = -0.257, RSE 42%

    # IIV. Rodrigues 2017 reports omega values on the log-scale (consistent
    # with the exponential IIV model in Eq. 1: theta_i = theta_TV * exp(eta_i)
    # with eta_i ~ N(0, omega^2)). The reported numeric values therefore are
    # interpreted as the lognormal SDs; ini() takes variances, so each line
    # below is omega^2.
    etalcl     ~ 0.154         # 0.393^2; Rodrigues 2017 Table 3, empirical model: omega_CLOXC = 0.393, RSE 15%
    etalvc     ~ 0.361         # 0.601^2; Rodrigues 2017 Table 3, empirical model: omega_VcOXC = 0.601, RSE 22%
    etalq      ~ 0.844         # 0.919^2; Rodrigues 2017 Table 3, empirical model: omega_QOXC = 0.919, RSE 18%
    etalvp     ~ 1.588         # 1.26^2;  Rodrigues 2017 Table 3, empirical model: omega_VpOXC = 1.26, RSE 15%
    etalcl_mhd ~ 0.0552        # 0.235^2; Rodrigues 2017 Table 3, empirical model: omega_CLMHD = 0.235, RSE 14%
    etalvc_mhd ~ 0.0445        # 0.211^2; Rodrigues 2017 Table 3, empirical model: omega_VcMHD = 0.211, RSE 25%
    etalkbt    ~ 0.397         # 0.63^2;  Rodrigues 2017 Table 3, empirical model: omega_KBT = 0.63, RSE 16%

    # Residual error. Rodrigues 2017 Methods page 2697: "The residual error
    # model used was proportional for OXC and combined for MHD." The Monolix
    # "combined" form maps to nlmixr2's add() + prop() composition; SDs are
    # taken at the values reported in Table 3 footnote (a = additive in mg/L,
    # b = proportional fraction).
    propSd     <- 0.32;    label("OXC proportional residual SD (fraction)")                                    # Rodrigues 2017 Table 3, empirical model: sigma_OXC = 0.32, RSE 7%
    addSd_mhd  <- 0.993;   label("MHD additive residual SD (mg/L)")                                             # Rodrigues 2017 Table 3, empirical model: sigma_MHD a (additive) = 0.993, RSE 13%
    propSd_mhd <- 0.0398;  label("MHD proportional residual SD (fraction)")                                     # Rodrigues 2017 Table 3, empirical model: sigma_MHD b (proportional) = 0.0398, RSE 21%
  })

  model({
    # Reference body weight for empirical allometric scaling (Rodrigues 2017
    # Methods Eq. 2: cov_median fixed to standard adult value of 70 kg).
    ref_wt <- 70

    # Individual parameters - parent OXC. Q_OXC/F and Vp_OXC/F are NOT
    # weight-scaled in the final empirical model.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Individual parameters - MHD metabolite. The EIAED effect from Table 3
    # is on the ABSENCE-of-EIAED indicator; convert from canonical
    # CONMED_EIAED (1 = on EIAED) via (1 - CONMED_EIAED).
    cl_mhd <- exp(lcl_mhd + etalcl_mhd + e_eiaed_cl_mhd * (1 - CONMED_EIAED)) *
              (WT / ref_wt)^e_wt_cl_mhd
    vc_mhd <- exp(lvc_mhd + etalvc_mhd) * (WT / ref_wt)^e_wt_vc_mhd

    # Back-transformation rate (1/h).
    kbt <- exp(lkbt + etalkbt)

    # Micro-constants. With Fm fixed to 1 (Rodrigues 2017 Results page 2699:
    # "The fraction of OXC metabolized to MHD (Fm) was fixed to 1"), the
    # parent's elimination flux enters the metabolite compartment 1:1 in
    # mass-equivalent units (no MW correction in the source paper).
    kel  <- cl / vc
    kelm <- cl_mhd / vc_mhd
    k12  <- q / vc
    k21  <- q / vp

    # ODE system (parent depot + 2-compartment OXC, 1-compartment MHD with
    # bidirectional MHD <-> OXC exchange via KBT).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (kel + k12) * central + k21 * peripheral1 + kbt * central_mhd
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_mhd) <-  kel * central - kelm * central_mhd - kbt * central_mhd

    # Bioavailability on depot (anchored to 1 per the paper).
    f(depot) <- exp(lfdepot)

    # Observations.
    Cc     <- central     / vc
    Cc_mhd <- central_mhd / vc_mhd

    Cc     ~ prop(propSd)
    Cc_mhd ~ add(addSd_mhd) + prop(propSd_mhd)
  })
}

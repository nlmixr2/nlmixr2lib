Wang_2012_levetiracetam <- function() {
  description <- "One-compartment population PK model for levetiracetam (LEV) in Chinese pediatric epilepsy patients (Wang 2012). First-order oral absorption and linear elimination (NONMEM ADVAN2 TRANS2). Body weight is the only retained covariate; it enters CL/F as a power-style allometric term with reference weight 25 kg (cohort median)."
  reference   <- "Wang YH, Wang L, Lu W, Shang DW, Wei MJ, Wu Y. Population pharmacokinetics modeling of levetiracetam in Chinese children with epilepsy. Acta Pharmacol Sin. 2012 Jun;33(6):845-851. doi:10.1038/aps.2012.57"
  vignette    <- "Wang_2012_levetiracetam"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 25 kg = cohort median (Wang 2012 Discussion paragraph 1; Table 1). Power-style scaling on CL/F only via (WT/25)^e_wt_cl. Wang 2012 treated weight as baseline-only; no time-varying schedule was reported.",
      source_name        = "WEIG"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 311,
    n_studies       = 1,
    n_observations  = 368,
    age_range       = "0.5-14 years (Wang 2012 Table 1).",
    age_median      = "Mean 6.34 years (Wang 2012 Table 1).",
    weight_range    = "5-70 kg (Wang 2012 Table 1).",
    weight_median   = "25 kg (Wang 2012 Discussion paragraph 1); mean 25.17 kg (Table 1).",
    sex_female_pct  = 48.6,
    race_ethnicity  = c(Asian_Chinese = 100),
    disease_state   = "Children with epilepsy: partial, generalized, and undetermined syndromes (Wang 2012 Methods - Patients). Normal renal and hepatic function in all subjects.",
    dose_range      = "20-60 mg/kg/day (Methods - Patients); total daily dose 250-2000 mg, mean 655 mg/d (Table 1).",
    regions         = "China (single centre: Peking University First Hospital, Beijing).",
    co_medication   = "40% on LEV monotherapy; 60% on combination AED therapy. Most common concomitant AEDs: valproic acid (VPA), lamotrigine (LTG), carbamazepine (CBZ), oxcarbazepine (OXC), topiramate (TPM). No significant interaction with concomitant AEDs detected (Wang 2012 Results paragraph 2; Table 2).",
    notes           = "Retrospective single-centre population PK analysis. The PPK model group used 311 patients (368 concentration-time points); a separate 50-patient PPK validation group (50 points) was used for external prediction-error checks. Sampling intervals 0.1-13 hours after last dose; trough-style with limited absorption-phase coverage."
  )

  ini({
    # Structural parameters - final-model estimates from Wang 2012 Table 4.
    lka <- log(1.56);  label("Absorption rate constant (Ka, 1/h)")               # Wang 2012 Table 4: Ka = 1.56 /h (RSE 14.3%, 95% CI 1.230-1.997)
    lvc <- log(12.1);  label("Apparent central volume of distribution (V/F, L)") # Wang 2012 Table 4: V/F = 12.1 L (RSE 5.6%, 95% CI 10.767-13.433)
    lcl <- log(1.04);  label("Apparent oral clearance (CL/F, L/h)")              # Wang 2012 Table 4: CL/F = 1.04 L/h (RSE 1.4%, 95% CI 1.011-1.069)

    # Covariate effect: body weight on CL/F, power form
    # CL/F = 1.04 * (WT/25)^e_wt_cl  (Wang 2012 final-model equation, page 848)
    # The in-text equation on page 848 reports the exponent as 0.563; Table 4
    # reports 0.567 with RSE 6.7% and 95% CI 0.489-0.637. The table value (0.567)
    # is used here as the formal point estimate; the in-text 0.563 is a rounding
    # or typo of the same quantity.
    e_wt_cl <- 0.567;  label("Power exponent: WT on CL/F (unitless)")            # Wang 2012 Table 4: 0.567 (RSE 6.7%, 95% CI 0.489-0.637)

    # Inter-individual variability: log-normal P_i = P_TV * exp(eta) (Wang 2012
    # Statistical model, page 847). Table 4 lists the omega values as variances.
    # No eta on Ka: Wang 2012 fixed omega_Ka at 0 because the sparse sampling did
    # not characterise the absorption phase (Discussion, "Patient data" paragraph,
    # page 848). Ka therefore has no eta term in model() below.
    etalcl ~ 0.195                                                                # Wang 2012 Table 4: omega^2 of CL/F = 0.195 (RSE 7.0%)
    etalvc ~ 0.163                                                                # Wang 2012 Table 4: omega^2 of V/F = 0.163 (RSE 45.5%)

    # Residual error. Wang 2012 writes the residual model in linear-concentration
    # space as E_ij^0 = E_ij + eps_ij with "variance is sigma_E^2" (Statistical
    # model, page 847), i.e. additive. Table 4 reports the variance estimate
    # epsilon = 0.028 (RSE 18.6%, 95% CI 0.018-0.035), so the additive SD is
    # sqrt(0.028) ~= 0.1673 mg/L. The magnitude is unusually small for an HPLC
    # assay over a 4.85-116.11 mg/L concentration range and likely reflects a
    # large fraction of residual variability absorbed into the IIV terms (omega
    # V/F shrinkage 44.9%); see vignette "Assumptions and deviations".
    addSd <- sqrt(0.028); label("Additive residual error SD (mg/L)")             # Wang 2012 Table 4: sigma_E^2 = 0.028 -> SD = sqrt(0.028)
  })

  model({
    # No IIV on Ka (omega_Ka fixed at 0 by the source paper).
    ka <- exp(lka)

    # Power-style covariate model on CL/F with reference weight 25 kg.
    cl <- exp(lcl + etalcl) * (WT / 25)^e_wt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment with first-order absorption (NONMEM ADVAN2 TRANS2).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L (the source concentration unit).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}

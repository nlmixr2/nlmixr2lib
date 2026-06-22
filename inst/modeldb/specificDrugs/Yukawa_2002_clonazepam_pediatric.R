Yukawa_2002_clonazepam_pediatric <- function() {
  description <- "Steady-state population PK model for clonazepam relative clearance (CL/F) in 137 Japanese pediatric and adult epileptic patients (Yukawa 2002 Table III row 4). CL/F is a body-weight power function with a 3-tier drug-interaction factor for concomitant antiepileptic drugs (monotherapy, +1 AED (CBZ or VPA), +>=2 AEDs)."
  reference   <- "Yukawa E, Satou M, Nonaka T, Yukawa M, Ohdo S, Higuchi S, Kuroda T, Goto Y. Pharmacoepidemiologic investigation of clonazepam relative clearance by mixed-effect modeling using routine clinical pharmacokinetic data in Japanese patients. J Clin Pharmacol. 2002;42(1):81-88."
  vignette    <- "Yukawa_2002_clonazepam_pediatric"
  units       <- list(time = "hr", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yukawa 2002 Methods Eq. 'CL = theta1 * TBW^theta2' (Table III, page 84 row 'CL = theta1 * TBW^theta2'). The paper does not state an explicit reference weight; the canonical (WT/1 kg) parameterisation is used inside model() to preserve the paper's reported theta1 = 152 ml/kg/h and theta2 = -0.181 verbatim. Cohort total-body-weight range 5.5-75 kg (Table II, page 83).",
      source_name        = "TBW"
    ),
    CONMED_AED = list(
      description        = "Indicator for any concomitant antiepileptic drug coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (clonazepam monotherapy)",
      notes              = "Yukawa 2002 Results page 84: DIF = 1 (CONMED_AED = 0; monotherapy reference), DIF = 1.18 (CONMED_AED = 1 and CONMED_AED_GE2 = 0; concomitant carbamazepine OR valproate alone), DIF = 2.12 * TBW^(-0.119) (CONMED_AED = 1 and CONMED_AED_GE2 = 1; >=2 AEDs alongside clonazepam). Paired with CONMED_AED_GE2 inside model() to recover the paper's 3-tier comed factor. Qualifying AEDs in the source cohort (Table I): carbamazepine (CBZ), valproate (VPA), phenobarbital (PB), phenytoin (PHT), zonisamide (ZSM), ethosuximide (ETOX).",
      source_name        = "CONMED_AED"
    ),
    CONMED_AED_GE2 = list(
      description        = "Indicator for >=2 concomitant antiepileptic drugs",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (zero or one concomitant AED)",
      notes              = "Yukawa 2002 Results page 84: the >=2-AEDs tier carries a weight-dependent DIF = 2.12 * TBW^(-0.119). CONMED_AED_GE2 = 1 must imply CONMED_AED = 1 (data assemblers must enforce this consistency: a record with CONMED_AED_GE2 = 1 and CONMED_AED = 0 is malformed). The paper's prose 'more than two antiepileptic drugs' is a translation artifact; Table I unambiguously bins Polytherapy (B) as >=2 AEDs alongside clonazepam (e.g. CZP+VPA+CBZ has exactly 2 AEDs and is grouped in Polytherapy B).",
      source_name        = "CONMED_AED_GE2"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 137L,
    n_observations   = 259L,
    n_studies        = 1L,
    age_range        = "0.3-28 years",
    age_median       = "monotherapy mean 10.4 (SD 4.9) yr; polytherapy A mean 8.1 (SD 4.6) yr; polytherapy B mean 12.6 (SD 6.3) yr",
    weight_range     = "5.5-75 kg",
    weight_median    = "monotherapy mean 31.3 (SD 16.9) kg; polytherapy A mean 24.5 (SD 13.8) kg; polytherapy B mean 34.7 (SD 17) kg",
    sex_female_pct   = 41.3,
    race_ethnicity   = "Japanese",
    disease_state    = "Epileptic patients (pediatric + adult) on chronic oral clonazepam maintenance therapy. 31 patient-periods on clonazepam monotherapy; 62 on clonazepam + 1 of {carbamazepine, valproate}; 67 on clonazepam + >=2 AEDs (Table II, page 83). 160 patient-periods total across 137 unique patients (some patients contributed to multiple comed tiers as regimens changed).",
    dose_range       = "Daily clonazepam dose 5.5-128.2 ug/kg/day across the three tiers (Table II); two or three divided doses per day as tablet or fine-granule preparation. Sampling at 2-6 hours after the morning dose at steady state (Methods page 82).",
    regions          = "Japan",
    notes            = "Yukawa 2002 Tables I-II baseline demographics (page 83). 66 female and 94 male patient-periods (parentheses in Table II indicate female counts: 16 + 23 + 27). All patients had normal renal and hepatic function and had been on clonazepam for >1 month before the analysis window. Serum concentrations measured by HPLC (Nonaka et al. method, CV <10%) as part of routine TDM care. The model is a steady-state CL/F regression: Css_ij = D_ij / (CL_ij * tau_ij). Bioavailability (F) and clearance are not separable; the estimated CL is a relative clearance because Css_ij is sampled at 2-6 h post-morning-dose rather than as a time-averaged Css (Methods page 83)."
  )

  ini({
    # Structural CL/F - Yukawa 2002 Results page 84 / Table III row 4 (model #4 in
    # the hypothesis-test sequence). The paper reports clearance in ml/kg/h (per-kg
    # body weight). The conversion to total CL in L/h is applied inside model() so
    # the ODE elimination term uses canonical L/h units while ini() values map 1:1
    # to the paper-reported thetas.
    lcl     <- log(152);    label("CL/F per kg body weight at TBW = 1 kg, clonazepam monotherapy (ml/kg/h)")  # Table III row 4: theta1 = 152 (95% CI 119-185); 27.6 SE/RSE not reported separately
    e_wt_cl <- fixed(-0.181); label("Allometric exponent on TBW for body-weight-normalised CL/F (unitless)")  # Table III row 4: theta2 = -0.181 (95% CI -0.248 to -0.114)

    # Drug-interaction factor (DIF) covariate coefficients. Yukawa 2002 Results
    # page 84 / Table III row 4. All fixed: the paper provides point estimates
    # and 95 percent CIs but does not separately estimate IIV on the DIF
    # coefficients themselves.
    e_conmed_aed_cl         <- fixed(1.18);   label("Multiplicative factor on CL/F for CONMED_AED = 1 in the +1-AED tier (unitless)")  # Table III row 4: theta3 = 1.18 (95% CI 1.09-1.27)
    e_conmed_aed_ge2_cl     <- fixed(2.12);   label("Multiplicative intercept on CL/F for CONMED_AED_GE2 = 1 (>=2-AEDs tier) at TBW = 1 kg (unitless)")  # Table III row 4: theta4 = 2.12 (95% CI 1.05-3.35)
    e_wt_cl_conmed_aed_ge2  <- fixed(-0.119); label("Additional allometric exponent on TBW for CL/F in the >=2-AEDs tier (unitless)")  # Table III row 4: theta5 = -0.119 (95% CI -0.254 to 0.016; CI includes zero -> marginal but retained in the final model)

    # Structural ODE constants - NOT from Yukawa 2002. The paper's published model
    # is purely a steady-state CL/F regression (Css = D / (CL * tau)) and does not
    # require Vc or ka; both are needed only as ODE structural constants for
    # time-resolved simulation. Steady-state Cavg is independent of Vc and ka.
    # Values are clonazepam-typical literature constants for total drug; the
    # vignette Assumptions section documents the choice.
    lvc <- fixed(log(3.0 * 70)); label("Volume of distribution Vc (L; not from Yukawa 2002; clonazepam literature ~3 L/kg, 70 kg reference adult)")
    lka <- fixed(log(1.5));      label("First-order oral absorption rate ka (1/h; not from Yukawa 2002; clonazepam Tmax 1-4 h per paper page 7)")

    # IIV - Yukawa 2002 Results page 84: 'estimate of coefficient of variation
    # for interpatient variability in clearance was 20.1 percent (95% CI
    # 15.0-24.2)'. Following the squared-fractional-CV convention used by sister
    # Japanese popPK models in this registry (Yukawa 1990 phenytoin, Hashimoto
    # 1994 zonisamide): omega^2 = (CV/100)^2.
    etalcl ~ 0.0404  # 20.1 percent CV; 0.201^2 = 0.0404

    # Residual error - Yukawa 2002 Results page 84: 'coefficient of variation
    # for residual variability yielded 18.6 percent (95% CI 15.4-21.4)'.
    # Eq. for residual error (page 83): Css_ij = Css_ij_predicted * (1 + epsilon_ij).
    propSd <- 0.186; label("Proportional residual SD on Cc (fraction)")  # Results page 84: 18.6 percent CV
  })

  model({
    # Tier indicators for the 3-tier drug-interaction factor (Yukawa 2002 Results
    # page 84 / Table III row 4). Operator-confirmed encoding via two binary
    # canonicals (CONMED_AED and CONMED_AED_GE2) per the 2026-06-21 sidecar:
    #   monotherapy:        CONMED_AED = 0, CONMED_AED_GE2 = 0  -> ind_mono = 1
    #   1 AED (CBZ or VPA): CONMED_AED = 1, CONMED_AED_GE2 = 0  -> ind_one  = 1
    #   >=2 AEDs:           CONMED_AED = 1, CONMED_AED_GE2 = 1  -> ind_many = 1
    ind_mono <- 1 - CONMED_AED
    ind_one  <- CONMED_AED * (1 - CONMED_AED_GE2)
    ind_many <- CONMED_AED_GE2

    # DIF (drug interaction factor) per Yukawa 2002 Results page 84:
    #   monotherapy:        DIF = 1
    #   1 AED:              DIF = 1.18                       (= e_conmed_aed_cl)
    #   >=2 AEDs:           DIF = 2.12 * TBW^(-0.119)        (= e_conmed_aed_ge2_cl * WT^e_wt_cl_conmed_aed_ge2)
    dif <- ind_mono * 1 +
           ind_one  * e_conmed_aed_cl +
           ind_many * e_conmed_aed_ge2_cl * WT^e_wt_cl_conmed_aed_ge2

    # Body-weight-normalised CL in ml/kg/h (paper's Results page 84 expression):
    cl_per_kg <- exp(lcl + etalcl) * WT^e_wt_cl * dif
    # Convert to total CL in L/h for the ODE elimination term.
    cl <- cl_per_kg * WT / 1000

    # Volume of distribution and absorption rate (literature constants; not from
    # Yukawa 2002).
    vc <- exp(lvc)
    ka <- exp(lka)

    # 1-compartment oral PK. The paper's CL absorbs bioavailability F (CL is
    # really CL/F per Methods page 83); F = 1 in the ODE so the ratio of dose to
    # absorption-compartment is preserved at the paper's reported magnitude.
    # Units: dose in mg, central in mg, vc in L; central/vc is in mg/L = 1000
    # ng/mL. Multiply by 1000 so Cc reads in ng/mL, matching the paper's reported
    # Css magnitudes (Table II monotherapy mean 10.4 ng/mL, etc.).
    kel <- cl / vc
    Cc  <- 1000 * central / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc ~ prop(propSd)
  })
}

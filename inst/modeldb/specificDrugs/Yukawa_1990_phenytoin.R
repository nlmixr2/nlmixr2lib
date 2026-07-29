Yukawa_1990_phenytoin <- function() {
  description <- "Steady-state Michaelis-Menten population PK model for phenytoin in 334 Japanese epilepsy outpatients on chronic oral phenytoin (Yukawa 1990 Model 2). Covariate effects on Vmax (allometric body weight, co-anticonvulsants) and Km (age <15 yr, co-anticonvulsants); dose-dependent powder bioavailability."
  reference   <- "Yukawa E, Higuchi S, Aoyama T. Population pharmacokinetics of phenytoin from routine clinical data in Japan: an update. Chem Pharm Bull (Tokyo). 1990;38(7):1973-1976. doi:10.1248/cpb.38.1973"
  vignette    <- "Yukawa_1990_phenytoin"
  units       <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 60 kg per Yukawa 1990 Methods (page 1974, Eq. 2: 'Vm and Km are the parameter values for the standard patient (adult, weight 60 kg, PHT alone, tablet)'). Power-form effect on Vmax with exponent 0.737 (Model 2 estimate, Table III).",
      source_name        = "WT"
    ),
    CHILD = list(
      description        = "Indicator for pediatric subject (age < 15 years)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult, age >= 15 years)",
      notes              = "Yukawa 1990 Methods page 1974 (Eq. 3) defines an AGE indicator equal to 1 if the patient is more than 15 years old and theta_AGE otherwise; the canonical CHILD indicator inverts this orientation (CHILD = 1 - AGE_indicator) so 0 is the adult reference. When CHILD = 1, Km is multiplied by theta_AGE = 0.752 (Model 2 estimate, Table III), i.e., Km in <15 yr olds is 24.8 percent less than adults.",
      source_name        = "AGE_LT15"
    ),
    CONMED_AED = list(
      description        = "Indicator for any concomitant antiepileptic drug coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (PHT alone)",
      notes              = "Yukawa 1990 Methods page 1974 (Eqs. 2 and 3) defines a CO indicator equal to 1 if the patient takes PHT alone and theta_co otherwise; the canonical CONMED_AED indicator inverts this orientation (CONMED_AED = 1 - CO_indicator) so 0 is the PHT-monotherapy reference. When CONMED_AED = 1, Vmax is multiplied by theta_coVm = 1.08 and Km by theta_coKm = 1.32 (Model 2 estimates, Table III). Co-anticonvulsants in the source dataset (Table I) include phenobarbital, carbamazepine, valproate, primidone, clonazepam, sultiame, ethotoin, ethosuximide, acetazolamide, and diazepam.",
      source_name        = "CO_AED"
    ),
    FORM_POWDER = list(
      description        = "Indicator for phenytoin powder oral formulation (1 = powder, 0 = tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet)",
      notes              = "Yukawa 1990 Methods page 1974 (Eq. 4) defines a BA indicator equal to 1 if PHT is prescribed as a tablet and 0 if as a powder; the canonical FORM_POWDER indicator inverts this (FORM_POWDER = 1 - BA_indicator) so 0 is the tablet reference. For tablet (FORM_POWDER = 0) bioavailability F = 1; for powder (FORM_POWDER = 1) Model 2's F = 1 - exp(-9.92 / DOSE_PHT_MGKGD). Both formulations are Aleviatin brand from Dainippon Pharmaceutical Co., Ltd.",
      source_name        = "FORM_POWDER"
    ),
    DOSE_PHT_MGKGD = list(
      description        = "Patient's own total daily dose of phenytoin per kg body weight",
      units              = "mg/kg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yukawa 1990 Methods page 1974 (Eq. 4): Dij is the daily dose of PHT for the ith Cpss in the jth patient, expressed in mg/kg/d. Used only in the powder-form bioavailability F = 1 - exp(-9.92 / DOSE_PHT_MGKGD); for tablet doses (FORM_POWDER = 0) the value is multiplied by 0 in the F expression and has no effect, but a non-NA non-zero placeholder must still be supplied. Compute as total daily dose (mg/d, summed across the 2-3 daily doses) divided by current body weight (kg).",
      source_name        = "Dij"
    )
  )

  population <- list(
    n_subjects       = 334L,
    n_observations   = 756L,
    n_studies        = 1L,
    age_range        = "0.6-71.1 years",
    age_median       = "mean 24.3 (SD 14.1) years",
    weight_range     = "9.0-115.0 kg",
    weight_median    = "mean 49.1 (SD 15.5) kg",
    sex_female_pct   = 49.1,
    race_ethnicity   = "Japanese (single-centre cohort at Kyushu University Hospital, Fukuoka)",
    disease_state    = "Epileptic outpatients on chronic oral phenytoin maintenance therapy. 101 patients on PHT monotherapy; 233 on PHT combined with phenobarbital, carbamazepine, valproate, primidone, clonazepam, sultiame, ethotoin, ethosuximide, acetazolamide, and/or diazepam in various combinations (Yukawa 1990 Table I).",
    dose_range       = "Daily dose mean 225.8 (SD 73.1) mg/d. Aleviatin brand tablets and powders (Dainippon Pharmaceutical Co., Ltd., Osaka, Japan) prescribed two to three times daily.",
    regions          = "Japan (Kyushu University Hospital, Fukuoka, single centre).",
    notes            = "Yukawa 1990 Tables I and II baseline demographics. Steady-state PHT serum concentration mean 9.78 (SD 7.77) ug/mL (target therapeutic range 10-20 ug/mL). 170 male, 164 female. Sample timing: 2-5 hours post-dose at steady state, with concentration determined at least 30 days after any dose change. Concurrent therapy was not altered during the analysis window. All patients had normal renal and hepatic function. Two to four steady-state R-Cpss observations per patient at different daily doses (756 paired records total)."
  )

  ini({
    # Structural Michaelis-Menten parameters - Yukawa 1990 Table III Model 2 column.
    # Reference covariates: 60 kg adult on PHT alone with tablet formulation.
    lvmax <- log(325);  label("Maximum elimination rate Vmax at reference covariates (mg/d)")  # Table III Model 2: Vm = 325 mg/d (SEM 9.18)
    lkm   <- log(2.41); label("Michaelis-Menten constant Km at reference covariates (mg/L)")  # Table III Model 2: Km = 2.41 ug/mL (SEM 0.25); 1 ug/mL = 1 mg/L

    # Structural ODE constants - NOT estimated by Yukawa 1990. The paper's published
    # model is purely a steady-state R-Cpss regression (R*F = Vmax*Cpss/(Km+Cpss))
    # which does not require Vc or ka; both are needed only as ODE structural
    # constants for time-resolved simulation. Steady-state Cpss is independent of
    # Vc and ka. Values are phenytoin-typical literature constants for total drug.
    lvc <- fixed(log(36));      label("Volume of distribution Vc (L; not from Yukawa 1990; phenytoin literature 0.6 L/kg x 60 kg)")
    lka <- fixed(log(1.5 * 24)); label("First-order oral absorption rate ka (1/day; not from Yukawa 1990; phenytoin tablet literature 1.5 1/h)")

    # Covariate effects - Yukawa 1990 Table III Model 2; all fixed since the paper
    # reports point estimates that fully specify the covariate model without IIV
    # on the effect coefficients themselves.
    e_wt_vmax         <- fixed(0.737); label("Power exponent of (WT/60) on Vmax (unitless)")  # Table III Model 2: theta_pw = 0.737 (SEM 0.036)
    e_child_km        <- fixed(0.752); label("Multiplicative factor on Km when CHILD = 1 (unitless)")  # Table III Model 2: theta_AGE = 0.752 (SEM 0.073)
    e_conmed_aed_vmax <- fixed(1.08);  label("Multiplicative factor on Vmax when CONMED_AED = 1 (unitless)")  # Table III Model 2: theta_coVm = 1.08 (SEM 0.039)
    e_conmed_aed_km   <- fixed(1.32);  label("Multiplicative factor on Km when CONMED_AED = 1 (unitless)")  # Table III Model 2: theta_coKm = 1.32 (SEM 0.165)
    e_form_powder_f   <- fixed(9.92);  label("Powder bioavailability parameter (F_powder = 1 - exp(-theta/DOSE_PHT_MGKGD)) (mg/kg/d)")  # Table III Model 2: theta_BA2 = 9.92 (SEM 0.52)

    # IIV - Yukawa 1990 Table III Model 2 reports CV percent for inter-individual
    # variability on Vm and Km. Following the squared-fractional-CV approximation
    # omega^2 ~ (CV/100)^2 used by older NONMEM popPK papers for exponential-IIV
    # models (matches the Hennig 2015 phenytoin convention in this registry).
    etalvmax ~ 0.0372  # Table III Model 2: 19.3 percent CV (SEM 2.2); omega^2 = 0.193^2 = 0.0372
    etalkm   ~ 0.398   # Table III Model 2: 63.1 percent CV (SEM 5.1); omega^2 = 0.631^2 = 0.398

    # Residual error - Yukawa 1990 Table III Model 2: 10.3 percent CV intra-individual.
    # Yukawa 1990 Eq. 6 applies the residual to the predicted daily dose R, not to
    # observed Cpss. At steady state R is monotonic in Cpss but the elasticity
    # Km/(Km+Cpss) makes the equivalent Cc proportional CV larger than 10.3 percent
    # in the saturated regime. The vignette Assumptions section discusses this
    # re-parameterization.
    propSd <- 0.103; label("Proportional residual SD on Cc (fraction)")  # Table III Model 2: theta_E = 10.3 percent (SEM 0.5)
  })

  model({
    # Individual M-M parameters with covariate effects (reference 60 kg adult, PHT alone, tablet).
    # Multiplicative ^CONMED_AED form: 1 when CONMED_AED = 0 (reference, PHT alone), 1.08 or 1.32
    # when CONMED_AED = 1; mirrors the indicator-as-multiplier formulation in Yukawa 1990 Eqs. 2-3.
    vmax <- exp(lvmax + etalvmax) * (WT / 60)^e_wt_vmax * e_conmed_aed_vmax^CONMED_AED
    km   <- exp(lkm + etalkm) * e_child_km^CHILD * e_conmed_aed_km^CONMED_AED

    # Bioavailability: F = 1 for tablet (FORM_POWDER = 0), F = 1 - exp(-9.92 / DOSE_PHT_MGKGD)
    # for powder (FORM_POWDER = 1). Yukawa 1990 Methods page 1974 Eq. 4 (Model 2).
    fdepot_powder <- 1 - exp(-e_form_powder_f / DOSE_PHT_MGKGD)
    fdepot        <- (1 - FORM_POWDER) + FORM_POWDER * fdepot_powder

    # Volume of distribution and absorption rate (literature constants; not from Yukawa 1990).
    vc <- exp(lvc)
    ka <- exp(lka)

    # 1-compartment oral PK with Michaelis-Menten elimination (no IV route used in source).
    # vmax in mg/day, km in mg/L, Cc in mg/L gives a M-M term in mg/day matching d/dt(central).
    Cc <- central / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - vmax * Cc / (km + Cc)

    f(depot) <- fdepot

    Cc ~ prop(propSd)
  })
}

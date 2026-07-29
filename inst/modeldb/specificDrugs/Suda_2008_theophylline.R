Suda_2008_theophylline <- function() {
  description <- "Steady-state population PK model for oral theophylline in 52 Japanese premature neonates and infants with apnea (Suda 2008). One-compartment first-order absorption structure; oral clearance CL/F is the only structural parameter the paper estimates (steady-state trough analysis Css = R / CL/F). Body-weight allometric scaling and a binary indicator for the Apnecut formulation (vs the in-house theophylline-alcohol comparator) on CL/F."
  reference   <- "Suda Y, Hanada K, Tsuchiwata S, Saito M, Nakamura T, Ito Y, Ishikawa Y, Kushida K, Ogata H. Population pharmacokinetic analysis of two theophylline formulations in premature neonates and infants with apnea. Yakugaku Zasshi. 2008;128(4):635-641. doi:10.1248/yakushi.128.635"
  vignette    <- "Suda_2008_theophylline"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (current, time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Suda 2008 final model (page 638 / Table 3) scales CL/F as (BW(g)/1000)^1.08, i.e. body weight in grams divided by 1000. With WT expressed in kg in nlmixr2lib, this is equivalent to (WT/1)^1.08, so the implicit reference body weight is 1 kg (1000 g). The studied range is 841-2548 g (Table 1); simulations at weights well outside this range are an extrapolation.",
      source_name        = "BW"
    ),
    FORM_THEO_APNECUT = list(
      description        = "Binary indicator of the Apnecut oral theophylline product (Kowa Co., Ltd.); 1 = Apnecut (APC), 0 = theophylline alcohol (TA) in-house preparation reference",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (TA in-house preparation, the registry reference)",
      notes              = "Suda 2008 final model (page 638): CL/F is multiplied by (1 - 0.282 * AP) with AP = 1 for Apnecut and 0 for theophylline alcohol. Canonical FORM_THEO_APNECUT preserves the source orientation. The two formulations are both oral liquids: TA is an in-house theophylline-alcohol preparation (5 mg/mL, with ethanol final concentration 10 percent) compounded at the National Center for Child Health and Development; APC is the commercial Apnecut internal-use solution (4 mg/mL aqueous) from Kowa Co., Ltd. (Methods Section 3, page 637).",
      source_name        = "AP"
    )
  )

  population <- list(
    species          = "human (Japanese premature neonates and infants with apnea)",
    n_subjects       = 52L,
    n_observations   = 90L,
    n_studies        = 1L,
    age_range        = "9-107 postnatal days (mean 34, SD 21)",
    weight_range     = "841-2548 g (mean 1465, SD 365)",
    birth_weight_range = "528-2158 g (mean 1236, SD 380)",
    pca_range        = "30-42 weeks postconceptional age (mean 34, SD 2)",
    sex_female_pct   = 51.9,
    race_ethnicity   = "Japanese (single-centre cohort at National Center for Child Health and Development, Tokyo)",
    disease_state    = "Premature neonatal apnea (apnea of prematurity) - apnea episodes >20 seconds, or shorter episodes accompanied by bradycardia or cyanosis. Frequent in infants born <34 weeks gestational age.",
    dose_range       = "Approximately 4 mg/kg/d oral theophylline divided every 12 hours (q12h). Mean dose 2.14 mg per administration. TA (n=36): 2.49 (SD 0.78) mg; APC (n=16): 1.53 (SD 0.51) mg per dose (Table 2).",
    regions          = "Japan (National Center for Child Health and Development, Setagaya-ku, Tokyo - single centre).",
    notes            = "Retrospective chart review of patients hospitalised between June 2005 and January 2007. Plasma theophylline measured as part of routine therapeutic drug monitoring on the JCA-BM1650 automated analyser (JEOL Ltd.) by homogeneous enzyme immunoassay. Only steady-state trough samples (>= 4 days on the oral product, >= 7 hours since last dose) were included (Methods Section 2, page 637). Variables collected but not retained in the final model: sex (SEXF), postnatal age (PNA, days), corrected postconceptional age (PCA, weeks; significantly correlated with body weight per Figure 2, so dropped to avoid collinearity), gestational age (GA, weeks), birth weight (WT_BIRTH, g), Apgar scores at 1 and 5 minutes (AS1, AS5), and oxygen supply status (OXY). Sex and postnatal age were significant in univariate models but were excluded from the full model (Table 3). Oxygen supply was tested and rejected; the authors note (Discussion, page 640) that prior reports of an oxygen-supply effect were in much younger infants and that their steady-state cohort is older. The therapeutic theophylline window is 5-13 ug/mL; toxicity (tachycardia, abdominal distension, emesis) onsets quickly above 13 ug/mL (Methods Section 1, page 636)."
  )

  ini({
    # Structural CL/F - Suda 2008 final model (page 638, equation reproduced in Table 3 final row).
    # CL/F (L/h) = 0.0201 * (BW(g)/1000)^1.08 * (1 - 0.282 * AP)
    # with AP = 0 for TA and AP = 1 for APC. Intercept 0.0201 L/h is the CL/F at the implicit
    # reference of 1000 g (1 kg) with the TA formulation (the AP = 0 reference).
    lcl <- log(0.0201); label("Apparent oral clearance CL/F at reference 1 kg on TA formulation (L/h)")  # Suda 2008 final model, page 638

    # Allometric exponent on WT - Suda 2008 final model.
    e_wt_cl <- fixed(1.08); label("Allometric exponent of (WT/1 kg) on CL/F (unitless)")  # Suda 2008 final model, page 638; reported as the u_2 exponent in Table 3 row 'u_1 * (BW/1000)^u_2'

    # Formulation effect on CL/F - Suda 2008 final model.
    e_form_apnecut_cl <- fixed(-0.282); label("Multiplicative shift on CL/F for APC formulation vs TA reference (unitless)")  # Suda 2008 final model, page 638; reported as the u_3 coefficient in Table 3 row 'Final model ... * (1 + u_3 * AP)'

    # Structural parameters NOT estimated by Suda 2008. The published model is a steady-state
    # regression of dose rate against trough concentration (Methods Section 5, page 637:
    # 'steady-state blood drug concentration = dose rate / oral clearance') and does not require
    # Vc or ka. Both are needed only as ODE structural constants for time-resolved simulation.
    # Vc/F is derived from the paper-cited half-life t1/2 = 22.3 h (Methods Section 5, page 637:
    # 'mean half-life = 22.3 hours: Apnecut 10 mg attached document'); ka is a literature-typical
    # theophylline oral absorption rate.
    lvc <- fixed(log(0.647)); label("Apparent central volume Vc/F at reference 1 kg (L; not from Suda 2008; derived as CL/F_TA(ref) / kel with kel = ln(2)/22.3 h from paper-cited package insert)")  # not estimated by Suda 2008; derived from t1/2 = 22.3 h (paper-cited, Methods Section 5)
    lka <- fixed(log(1.5));   label("First-order oral absorption rate ka (1/h; not from Suda 2008; literature-typical theophylline oral absorption)")                                                       # not estimated by Suda 2008; literature value

    # IIV on CL/F - Suda 2008 final model reports an inter-individual CV of 15.0 percent
    # (page 638 and Table 3 final-model row). The paper specifies a log-normal IIV model
    # (Methods Section 5, page 637, Eq. 1: P_j = P_typ * exp(eta_Pj)), so the natural
    # variance conversion is omega^2 = log(1 + CV^2).
    etalcl ~ 0.02226  # Suda 2008 final model: CV 15.0 percent; omega^2 = log(1 + 0.15^2) = 0.02226

    # Residual (intra-individual) error - Suda 2008 final model reports an intra-individual
    # CV of 15.3 percent (page 638, Table 3 final-model row). The paper specifies a
    # proportional (relative) residual model (Methods Section 5, page 637, Eq. 2:
    # Cp_ij_obs = Cp_ij_pred + Cp_ij_pred * eps_ij1).
    propSd <- 0.153; label("Proportional residual SD on Cc (fraction)")  # Suda 2008 final model: 15.3 percent CV intra-individual
  })

  model({
    # Individual oral clearance with allometric WT scaling and Apnecut-vs-TA formulation shift.
    # Suda 2008 final model: CL/F = 0.0201 * (BW(g)/1000)^1.08 * (1 - 0.282 * AP), reproduced
    # in nlmixr2lib's (WT(kg)/1)^exponent form because (BW(g)/1000) numerically equals WT(kg).
    cl <- exp(lcl + etalcl) * (WT / 1)^e_wt_cl * (1 + e_form_apnecut_cl * FORM_THEO_APNECUT)

    # Apparent central volume scales with body weight allometrically (literature default
    # exponent 1 for theophylline distribution in neonates); ka is body-size-independent.
    vc <- exp(lvc) * (WT / 1)
    ka <- exp(lka)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

Bergmann_2014_tacrolimus <- function() {
  description <- "Two-compartment population PK model for oral tacrolimus in adult kidney transplant recipients (Bergmann 2014), with first-order absorption after a lag time, allometric (WT/70 kg)^0.75 scaling on apparent clearance, multiplicative CYP3A5*1-carrier effect on CL/F, linear hematocrit and post-transplant-day effects on CL/F, linear free prednisolone Cmax effect on V1/F, correlated inter-individual variability across V1/F, ka, and V2/F, and proportional residual error."
  reference <- "Bergmann TK, Hennig S, Barraclough KA, Isbel NM, Staatz CE. Population Pharmacokinetics of Tacrolimus in Adult Kidney Transplant Patients: Impact of CYP3A5 Genotype on Starting Dose. Ther Drug Monit. 2014;36(1):62-70. doi:10.1097/FTD.0b013e31829f1ab8"
  vignette <- "Bergmann_2014_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Bergmann 2014. Allometric power scaling on CL/F with reference 70 kg and the theory-based exponent 0.75 fixed (Bergmann 2014 Table 2 footnote). Study median 79 kg (10-90 percentile 59-101).",
      source_name        = "WT"
    ),
    HCT = list(
      description        = "Hematocrit, expressed as a fraction of total blood volume",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Centred at 0.33 (Bergmann 2014 study-population median; 10-90 percentile 0.25-0.40). Source paper reports HCT as a fraction (0-1), not as percent (0-100); the canonical-register HCT entry's units (%) are explicitly overridden here so the centring value 0.33 and the linear-deviation coefficient -1.01 reproduce the paper's equation directly. To use a dataset that records HCT in percent, multiply the column by 0.01 before passing it to this model.",
      source_name        = "HEM"
    ),
    POD = list(
      description        = "Days post-transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject. Centred at 22.7 days (Bergmann 2014 Table 2 footnote; close to the 23-day median in Table 1). Values greater than 180 days are capped at 180 inside model() per the source paper's covariate equation (linear effect plateaus beyond 180 days post-transplant). 83% of dataset observations were within the first 90 days post-transplant.",
      source_name        = "POD"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser)",
      notes              = "Time-fixed germline genotype derived from rs776746 (CYP3A5 6986A>G): the *1 (A) allele encodes functional CYP3A5 protein; the *3 (G) allele creates a cryptic splice site and yields nonfunctional protein. In the Bergmann 2014 cohort the genotype distribution was *1/*1 = 3 (1.7%), *1/*3 = 23 (13.3%), *3/*3 = 146 (84.4%); 1 patient with failed genotyping was assigned to *3/*3. CYP3A5_EXPR = 1 for the 26 *1 carriers, 0 for the 147 nonexpressers. Multiplicative effect on CL/F as `theta_CYP3A5 ^ CYP3A5_EXPR` with `theta_CYP3A5 = 1.60` (60% higher CL/F in expressers).",
      source_name        = "X"
    ),
    PRED_CMAX_FREE = list(
      description        = "Maximum free (ultrafiltrable) plasma prednisolone concentration over a tacrolimus dosing interval",
      units              = "nmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject co-medication-exposure covariate (one Cmax per subject). Centred at 155.5 nmol/L in Bergmann 2014 Table 2 footnote (study median 162 nmol/L per Table 1). Values derived from limited-sampling free prednisolone concentrations at 1, 2, and 4 hours postdose per Bergmann 2014 Methods. Linear deviation effect on V1/F: every 1 nmol/L increase from 155.5 nmol/L decreases apparent central volume by 0.28%.",
      source_name        = "PredCmax,free"
    )
  )

  population <- list(
    n_subjects        = 173L,
    n_studies         = 2L,
    n_observations    = 1554L,
    age_range         = "33-64 years (10-90 percentile)",
    age_median        = "50 years",
    weight_range      = "59-101 kg (10-90 percentile)",
    weight_median     = "79 kg",
    sex_female_pct    = 33.5,
    race_ethnicity    = "Not reported in source paper (single-centre Brisbane, Australia cohort).",
    disease_state     = "Adult kidney transplant recipients on protocol immunosuppression (basiliximab + methylprednisolone induction, oral tacrolimus + oral prednisolone + mycophenolate mofetil maintenance).",
    dose_range        = "Oral tacrolimus 3-8 mg twice daily (median 5 mg per dose). Initial dose 0.075 mg/kg twice daily started preoperatively, adjusted by therapeutic drug monitoring to a target trough of 6-10 ug/L during the first 3 months posttransplant.",
    regions           = "Australia (Princess Alexandra Hospital, Brisbane).",
    cyp3a5_distribution = "*1/*1 (rs776746 AA) n = 3 (1.7%); *1/*3 (AG) n = 23 (13.3%); *3/*3 (GG) n = 146 (84.4%); 1 patient with failed genotyping assigned to *3/*3. Population in Hardy-Weinberg equilibrium.",
    pod_range         = "2-104 days (10-90 percentile); 83% of measurements within the first 90 days posttransplant.",
    hct_median        = "0.33 (10-90 percentile 0.25-0.40)",
    pred_cmax_free_median = "162 nmol/L (10-90 percentile 85-260)",
    notes             = "Pooled analysis of two consecutive prospective studies at the same centre. Study 1 (n = 20): full 13-point concentration-time profile (predose plus 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 6, 9, and 12 h postdose); half sampled in their first posttransplant week, half >90 days posttransplant. Study 2 (n = 153): limited concentration-time profile (predose plus 1, 2, and 4 h postdose) on at least one occasion 4-12 months posttransplant. Tacrolimus measured by HPLC-MS/MS in whole blood (assay range 0.5-50 ug/L)."
  )

  ini({
    # Structural PK -- Bergmann 2014 Table 2 final-model estimates. Time in
    # hours; apparent clearances (CL/F, Q/F) in L/h; apparent volumes (V1/F,
    # V2/F) in L. The Table 2 reference is a 70-kg adult who is a CYP3A5
    # nonexpresser (*3/*3) with HCT = 0.33, POD = 22.7 days, and free
    # prednisolone Cmax = 155.5 nmol/L; the table values are reproduced
    # verbatim and the covariate-equation factors below recreate
    # individual-specific values inside model().
    lka   <- log(0.35)  ; label("Absorption rate constant ka (1/h)")                                                 # Bergmann 2014 Table 2 final ka = 0.35 1/h
    ltlag <- log(0.44)  ; label("Absorption lag time (h)")                                                            # Bergmann 2014 Table 2 final Lag time = 0.44 h
    lcl   <- log(25.5)  ; label("Apparent oral clearance CL/F at WT = 70 kg, CYP3A5 nonexpresser, HCT 0.33, POD 22.7 d (L/h)") # Bergmann 2014 Table 2 final CL/F = 25.5 L/h
    lvc   <- log(113.0) ; label("Apparent central volume V1/F at free prednisolone Cmax = 155.5 nmol/L (L)")          # Bergmann 2014 Table 2 final V1/F = 113 L
    lq    <- log(67.9)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                                   # Bergmann 2014 Table 2 final Q/F = 67.9 L/h
    lvp   <- log(1060)  ; label("Apparent peripheral volume V2/F (L)")                                                # Bergmann 2014 Table 2 final V2/F = 1060 L

    # Covariate effects on CL/F -- Bergmann 2014 Table 2 footnote covariate
    # equation: CL/F_i = theta_CL/F * theta_CYP3A5^X * (1 + theta_HEM * (HEM -
    # 0.33)) * (WT/70)^0.75 * (1 + theta_POD * (POD - 22.7)).
    e_cyp3a5_expr_cl <- 1.60     ; label("CYP3A5*1-carrier multiplicative factor on CL/F (theta_CYP3A5; expressers have 60% higher CL/F)") # Bergmann 2014 Table 2 theta_CYP3A5 = 1.60
    e_hct_cl         <- -1.01    ; label("Hematocrit linear-deviation coefficient on CL/F (per unit fraction, centred at 0.33)")          # Bergmann 2014 Table 2 theta_HEM = -1.01
    e_pod_cl         <- -0.0021  ; label("Post-transplant-day linear-deviation coefficient on CL/F (per day, centred at 22.7 d, capped at 180 d)") # Bergmann 2014 Table 2 theta_POD = -0.21% per day
    e_wt_cl          <- 0.75     ; label("Allometric exponent of (WT/70 kg) on CL/F (unitless; fixed)")                                    # Bergmann 2014 Table 2 footnote: power 0.75 fixed (allometric theory)

    # Covariate effect on V1/F -- Bergmann 2014 Table 2 footnote: V1/F_i =
    # theta_V1/F * (1 + theta_PRED * (PredCmax,free - 155.5)).
    e_pred_cmax_free_vc <- -0.0028 ; label("Free prednisolone Cmax linear-deviation coefficient on V1/F (per nmol/L, centred at 155.5)")  # Bergmann 2014 Table 2 theta_PRED = -0.28% per nmol/L

    # Inter-individual variability -- Bergmann 2014 Table 2 reports IIV as
    # %CV. Convert to log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F  CV 29.5% -> log(0.295^2 + 1) = 0.0834
    #   V1/F  CV 46.8% -> log(0.468^2 + 1) = 0.1981
    #   V2/F  CV 89.4% -> log(0.894^2 + 1) = 0.5874
    #   ka    CV 47.6% -> log(0.476^2 + 1) = 0.2043
    # CL/F is independent. V1/F, ka, and V2/F form a correlated block; the
    # paper's Table 2 final-model correlations are r(V1/F, ka) = 0.677,
    # r(V1/F, V2/F) = -0.049, r(ka, V2/F) = -0.013. Off-diagonal covariances
    # = r * sqrt(omega^2_i) * sqrt(omega^2_j):
    #   cov(V1/F, ka)   =  0.677  * sqrt(0.1981) * sqrt(0.2043) =  0.1362
    #   cov(V1/F, V2/F) = -0.049  * sqrt(0.1981) * sqrt(0.5874) = -0.01672
    #   cov(ka,   V2/F) = -0.013  * sqrt(0.2043) * sqrt(0.5874) = -0.00451
    etalcl ~ 0.0834                                               # Bergmann 2014 Table 2 IIV CL/F = 29.5%
    etalvc + etalka + etalvp ~ c(0.1981,
                                  0.1362,   0.2043,
                                 -0.01672, -0.00451, 0.5874)      # Bergmann 2014 Table 2 IIV/correlations on V1/F, ka, V2/F

    # Residual unexplained variability -- Bergmann 2014 Table 2 reports a
    # proportional error model on the linear concentration scale.
    propSd <- 0.183 ; label("Proportional residual error (fraction)") # Bergmann 2014 Table 2 final Proportional RUV = 18.3%
  })

  model({
    # Cap POD at 180 days per Bergmann 2014 Table 2 footnote ("POD, capped at
    # 180 days. Final model: ...").
    pod_capped <- min(POD, 180)

    # Individual PK parameters with Bergmann 2014 covariate equations.
    # Reference subject: WT 70 kg, CYP3A5 nonexpresser (CYP3A5_EXPR = 0),
    # HCT 0.33, POD <= 22.7 d (capped at 180 d), free prednisolone Cmax
    # 155.5 nmol/L.
    cl <- exp(lcl + etalcl) *
          e_cyp3a5_expr_cl ^ CYP3A5_EXPR *
          (1 + e_hct_cl * (HCT - 0.33)) *
          (WT / 70) ^ e_wt_cl *
          (1 + e_pod_cl * (pod_capped - 22.7))
    vc <- exp(lvc + etalvc) * (1 + e_pred_cmax_free_vc * (PRED_CMAX_FREE - 155.5))
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)

    # Two-compartment oral PK with first-order absorption + lag time. Dose
    # lands in `depot`; bioavailability is fixed at 1 (Bergmann 2014 Methods:
    # "bioavailability fixed to 1 during model building"), so f(depot) is
    # not assigned. The micro-constants are spelled out for clarity.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Tacrolimus assay reports concentrations in ug/L (= ng/mL). Convert
    # internal mg / L (dose mg, vc L) to ng/mL by multiplying by 1000.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}

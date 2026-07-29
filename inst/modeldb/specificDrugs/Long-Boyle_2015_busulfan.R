"Long-Boyle_2015_busulfan" <- function() {
  description <- "One-compartment IV PK model with Michaelis-Menten elimination for busulfan in pediatric and young adult patients (0.1-24 yrs) undergoing hematopoietic cell transplant. Allometric body-weight scaling on intrinsic clearance (CLin, exponent fixed 0.75) and central volume (Vc, exponent fixed 1) with reference weight 22 kg; hockey-stick age effect on CLin (linear increase below the 12-yr breakpoint applied to AGE directly, multiplicative linear decrease above). Correlated IIV on CLin and Vc; combined proportional + additive residual error (Long-Boyle 2015)."
  reference   <- "Long-Boyle JR, Savic R, Yan S, Bartelink I, Musick L, French D, Law J, Horn B, Cowan MJ, Dvorak CC. Population Pharmacokinetics of Busulfan in Pediatric and Young Adult Patients Undergoing Hematopoietic Cell Transplant: A Model-Based Dosing Algorithm for Personalized Therapy and Implementation Into Routine Clinical Use. Ther Drug Monit. 2015;37(2):236-245. doi:10.1097/FTD.0000000000000131"
  vignette    <- "Long-Boyle_2015_busulfan"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CLin (exponent 0.75, fixed) and Vc (exponent 1, fixed) with reference weight 22 kg (cohort median). Treated as baseline (4-day busulfan course). Source range 3-101 kg (Long-Boyle 2015 Table 2).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age in years",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hockey-stick effect on CLin with breakpoint fixed at 12 yrs (Long-Boyle 2015 p. 240 CLin equations). Below the breakpoint the multiplier is (1 + SL,bp * AGE) with SL,bp = 0.032 per yr applied to AGE directly (structural reference AGE = 0; the peak multiplier 1 + 0.032 * 12 = 1.384 occurs at the breakpoint). Above the breakpoint the multiplier is (1 + SL,bp * 12) * (1 + SL>bp * (AGE - 12)) with SL>bp = -0.0138 per yr, giving a slow decline back toward the structural baseline at higher ages. Source range 0.1-24 yrs (Long-Boyle 2015 Table 2).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 90L,
    n_studies      = 1L,
    age_range      = "0.1-24 years",
    age_median     = "7 years",
    weight_range   = "3-101 kg",
    weight_median  = "22 kg",
    sex_female_pct = 41,
    disease_state  = "Pediatric and young adult patients undergoing autologous or allogeneic hematopoietic cell transplantation (HCT) for malignant and nonmalignant pediatric disorders. Conditioning chemotherapy included busulfan plus one of: fludarabine + serotherapy (ATG or alemtuzumab); fludarabine + thiotepa + serotherapy; fludarabine + clofarabine + serotherapy; or melphalan + serotherapy. Seizure prophylaxis with lorazepam or levetiracetam.",
    dose_range     = "Intravenous busulfan over a 2-hour infusion every 6 hours for 16 doses. Initial doses for 79 of 90 patients used the conventional weight-band nomogram: 1.1 mg/kg/dose for patients <=12 kg and 0.8 mg/kg/dose for patients >12 kg. In 11 patients a 0.5 mg/kg test dose 3-4 days before conditioning was used to estimate individual CL and select the first dose.",
    regions        = "USA (single center: UCSF Benioff Children's Hospital, San Francisco)",
    notes          = "Retrospective routine-TDM data collected at UCSF Benioff between January 2007 and April 2013 (Long-Boyle 2015 Methods and Table 2). 1165 quantifiable plasma busulfan concentrations in 90 subjects analyzed with NONMEM v7 FOCE-I. Baseline laboratory values (medians, Long-Boyle 2015 Table 2): serum creatinine 0.3 mg/dL (0.3-0.95), creatinine clearance 169 mL/min/m^2 (70-286), alkaline phosphatase 159 IU/L (46-1760), AST 30 IU/L (10-265), ALT 28 IU/L (5-525), total bilirubin 0.6 mg/dL (0.1-3.7). Therapeutic Css target 600-900 ng/mL (midpoint 750 ng/mL; AUC 4.5 mg.h/L over the 6-h interval)."
  )

  ini({
    # Structural CLin reference value. The paper Table 3 footnote labels this
    # 4.32 L/h as the typical value "for a child weighing 22 kg and 7 years of
    # age", but the explicit CLin equations on p. 240 multiply this THETA by
    # the AGE factor (1 + SL,bp * AGE) -- so 4.32 is the model THETA at the
    # structural reference (WT = 22 kg, AGE = 0); the equation predicts
    # CLin = 4.32 * 1.224 = 5.29 L/h at AGE = 7 yrs. Table 4 model-based doses
    # confirm this interpretation: e.g., 20 kg / 6 yr -> 22.1 mg implies CLi =
    # 22.1 / 4.5 = 4.91 L/h, which the AGE-applied-to-AGE formula reproduces
    # within 1-2 percent.
    lcl <- log(4.32)  ; label("Intrinsic clearance CLin THETA at 22 kg (L/h)")   # Long-Boyle 2015 Table 3 (RSE 8 percent)
    lvc <- log(15.7)  ; label("Central volume of distribution at 22 kg (L)")     # Long-Boyle 2015 Table 3 (RSE 3 percent)

    # Michaelis-Menten constant. Paper reports Km = 6704 ng/mL = 6.704 mg/L.
    # MM elimination is parameterised as rate = CLin * Cp / (1 + Cp/Km); at low
    # Cp the effective clearance is CLin (intrinsic clearance), at high Cp it
    # saturates at Vmax = CLin * Km.
    lkm <- log(6.704) ; label("Michaelis-Menten constant Km (mg/L)")              # Long-Boyle 2015 Table 3 (6704 ng/mL; RSE 43 percent)

    # Allometric exponents on WT, both fixed at canonical values per the paper
    # ("Exponent for effect of weight on CL in: 0.75 (fixed)" and "Exponent for
    # effect of weight on V c: 1 (fixed)", Long-Boyle 2015 Table 3).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CLin (unitless, fixed)")  # Long-Boyle 2015 Table 3
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on Vc (unitless, fixed)")    # Long-Boyle 2015 Table 3

    # Hockey-stick AGE effect on CLin (Long-Boyle 2015 p. 240 CLin equations).
    # The breakpoint is fixed at 12 yrs (Table 3: "BP ... 12 (fixed)"). Below
    # the breakpoint the multiplier is (1 + SL,bp * AGE); at AGE = 0 this
    # equals 1 (the structural reference) and at AGE = 12 it equals 1.384
    # (peak). Above the breakpoint the multiplier is the peak times
    # (1 + SL>bp * (AGE - 12)), with SL>bp = -0.0138 per yr producing a slow
    # decline back toward the structural baseline.
    e_age_le12 <-  0.032   ; label("Slope of AGE effect on CLin for AGE <= 12 yrs (1/yr)")  # Long-Boyle 2015 Table 3 (SL,bp; RSE 32 percent)
    e_age_gt12 <- -0.0138  ; label("Slope of AGE effect on CLin for AGE >  12 yrs (1/yr)")  # Long-Boyle 2015 Table 3 (SL>bp; RSE 46 percent)

    # Correlated IIV on CLin and Vc. Paper reports IIV as %CV (22 percent on
    # CLin, 29 percent on Vc) and correlation coefficient 0.42 between CLin
    # and Vc (Long-Boyle 2015 Table 3). Converted to the log-normal internal
    # variance scale via omega^2 = log(CV^2 + 1):
    #   var_lcl = log(0.22^2 + 1) = 0.047251
    #   var_lvc = log(0.29^2 + 1) = 0.080744
    # Covariance from the reported correlation:
    #   cov = 0.42 * sqrt(var_lcl * var_lvc) = 0.42 * sqrt(0.003816) = 0.025952
    etalcl + etalvc ~ c(0.047251,
                        0.025952, 0.080744)  # Long-Boyle 2015 Table 3 (IIV CLin 22 percent, IIV Vc 29 percent, corr 0.42)

    # Combined proportional + additive residual error. Paper reports
    # proportional 14.8 percent and additive 47 ng/mL (Long-Boyle 2015 Table 3).
    # Additive converted to mg/L: 47 ng/mL = 0.047 mg/L.
    propSd <- 0.148 ; label("Proportional residual error (fraction)")  # Long-Boyle 2015 Table 3 (14.8 percent; RSE 14.6 percent)
    addSd  <- 0.047 ; label("Additive residual error (mg/L)")          # Long-Boyle 2015 Table 3 (47 ng/mL; RSE 25 percent)
  })

  model({
    # Hockey-stick AGE effect on CLin (Long-Boyle 2015 p. 240 explicit
    # equations). The piecewise multiplicative form combines into a single
    # smooth expression using min(AGE, 12) and max(AGE - 12, 0):
    #   AGE <= 12: factor = 1 + e_age_le12 * AGE
    #              = (1 + e_age_le12 * min(AGE, 12)) * (1 + e_age_gt12 * 0)
    #   AGE >  12: factor = (1 + e_age_le12 * 12) * (1 + e_age_gt12 * (AGE - 12))
    #              = (1 + e_age_le12 * min(AGE, 12)) * (1 + e_age_gt12 * max(AGE - 12, 0))
    # Structural reference is AGE = 0 (factor = 1); peak at AGE = 12 is
    # factor = 1.384; AGE = 7 gives factor = 1.224 (i.e., the CLin THETA of
    # 4.32 L/h scaled to 5.29 L/h for the 22 kg / 7 yr median child).
    age_factor_cl <- (1 + e_age_le12 * min(AGE, 12)) *
                     (1 + e_age_gt12 * max(AGE - 12, 0))

    # Size scaling on actual body weight, reference 22 kg.
    size_cl <- (WT / 22)^e_wt_cl
    size_vc <- (WT / 22)^e_wt_vc

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * size_cl * age_factor_cl
    vc <- exp(lvc + etalvc) * size_vc
    km <- exp(lkm)

    # One-compartment IV with Michaelis-Menten elimination (Long-Boyle 2015
    # Results, p. 240). Dose targets `central` directly (busulfan is given as
    # a 2-h IV infusion); no depot. Equivalent forms:
    #   d/dt(central) = -CLin * Cp * Km / (Km + Cp)
    #                 = -CLin * Cp / (1 + Cp/Km)
    # so at Cp << Km the elimination is linear with CL = CLin, and at Cp >> Km
    # it saturates at Vmax = CLin * Km.
    conc <- central / vc
    d/dt(central) <- -cl * conc / (1 + conc / km)

    # Concentration in mg/L (= ug/mL). The paper reports Cc in ng/mL; multiply
    # the simulated Cc by 1000 in post-processing for direct comparison.
    Cc <- conc
    Cc ~ add(addSd) + prop(propSd)
  })
}

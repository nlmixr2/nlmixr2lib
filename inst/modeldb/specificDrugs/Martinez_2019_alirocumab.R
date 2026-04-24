Martinez_2019_alirocumab <- function() {
  description <- "Two-compartment population PK model for alirocumab in healthy volunteers and adults with hypercholesterolemia (Martinez 2019, Part I), with first-order SC absorption (with lag time), linear plus Michaelis-Menten (target-mediated) elimination from the central compartment, and logit-transformed bioavailability."
  reference   <- "Martinez JM, Brunet A, Hurbin F, DiCioccio AT, Rauch C, Fabre D. Population Pharmacokinetic Analysis of Alirocumab in Healthy Volunteers or Hypercholesterolemic Subjects Using a Michaelis-Menten Approximation of a Target-Mediated Drug Disposition Model - Support for a Biologics License Application Submission: Part I. Clin Pharmacokinet. 2019;58(1):101-113. doi:10.1007/s40262-018-0669-y"
  vignette    <- "Martinez_2019_alirocumab"
  units       <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying per the LOCF imputation in Martinez 2019)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on linear clearance CLL (additive): CLL = TVCLL + 2.92e-4 L/h/kg * (WT - 82.9) + 6.44e-3 L/h * STATIN. Reference 82.9 kg is the median body weight in the pooled Martinez 2019 dataset. Converted to day-units in the model file (x24).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age (time-varying per LOCF; practically time-fixed since trial durations <= 104 weeks)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on the peripheral volume of distribution V3: V3 = TVV3 * (AGE/60)^0.310. Reference age 60 years is the median of the Martinez 2019 dataset.",
      source_name        = "AGE"
    ),
    STATIN = list(
      description        = "Concomitant statin administration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no statin coadministration)",
      notes              = "Martinez 2019 codes STATIN = 1 for coadministration of rosuvastatin (< 20 mg/day), atorvastatin (< 40 mg/day), or simvastatin (any dose), and 0 otherwise. Additive effect on linear clearance CLL: +0.00644 L/h = +0.15456 L/day when STATIN = 1. Other lipid-lowering therapies (ezetimibe, fibrates) are not captured by STATIN.",
      source_name        = "STATIN"
    ),
    FPCSK9 = list(
      description        = "Free (unbound) serum proprotein convertase subtilisin/kexin type 9 concentration (time-varying, LOCF for missing)",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying additive effect on Km: Km = TVKM + (-0.541) * (FPCSK9/72.9). Reference FPCSK9 = 72.9 ng/mL is the median time-varying value in the Martinez 2019 dataset. Observed 5th-95th percentile range 0-392 ng/mL. Distinct from total PCSK9; only the free (drug-unbound) PCSK9 is the pharmacologically active target fraction.",
      source_name        = "FPCSK9"
    )
  )

  population <- list(
    n_subjects          = 2799L,
    n_observations      = 13717L,
    n_studies           = 13L,
    phases              = "Phase I, II, and III",
    age_median          = "60 years",
    weight_median       = "82.9 kg",
    disease_state       = "Healthy volunteers (HV) and adult patients with hypercholesterolemia (including familial and non-familial hypercholesterolemia; subset with established coronary heart disease) not adequately controlled on a maximally tolerated statin regimen or with statin intolerance.",
    dose_range          = "IV: 0.3-12 mg/kg single dose (n=30, phase I only). SC: 50-300 mg single or repeated (Q2W or Q4W over up to 104 weeks; current analysis included data up to 24 weeks). Marketed SC regimens are 75 mg or 150 mg Q2W.",
    regions             = "Multi-regional pool of 13 Sanofi/Regeneron trials including two Japanese cohorts (NCT01448317 phase I, NCT01812707 phase II); ODYSSEY phase III programme (MONO, COMBO II, FH I, LONG TERM).",
    fpcsk9_median       = "72.9 ng/mL (time-varying); 283 ng/mL baseline",
    fpcsk9_range_5_95   = "0-392 ng/mL (time-varying 5th-95th percentiles)",
    concomitant         = "STATIN = 1 for rosuvastatin (< 20 mg/day), atorvastatin (< 40 mg/day), or simvastatin (any dose). STATIN = 0 for monotherapy, other statins, or statin-intolerant subjects.",
    notes               = "Median body weight, age, and free-PCSK9 reference values are taken from Table 2 footnotes a-c of Martinez 2019. The covariate-screening list (sex, race, renal function, BMI, ADA status, injection site, injection device) was non-significant and is not carried as a model covariate. Trials included in the pooled dataset are listed in Table 1 of Martinez 2019 with NCT identifiers."
  )

  ini({
    # Structural PK parameters - Martinez 2019 Table 2 final-model typical values.
    # Reference covariates: WT = 82.9 kg, AGE = 60 years, STATIN = 0, FPCSK9 = 72.9 ng/mL.
    # Paper reports rates in h^-1; converted to day^-1 for the nlmixr2lib convention (x24).
    lka          <- log(7.68e-3 * 24); label("First-order SC absorption rate Ka (1/day)")                        # Martinez 2019 Table 2 (7.68e-3 /h)
    lcl          <- log(0.0124 * 24);  label("Linear clearance CLL at reference covariates (L/day)")             # Martinez 2019 Table 2 (0.0124 L/h)
    lvc          <- log(3.19);         label("Central (depot-to-central) volume of distribution V2 (L)")        # Martinez 2019 Table 2
    lvp          <- log(2.79);         label("Peripheral volume of distribution V3 at reference age (L)")       # Martinez 2019 Table 2 (age reference 60 y)
    lq           <- log(0.0185 * 24);  label("Intercompartmental clearance Q (L/day)")                          # Martinez 2019 Table 2 (0.0185 L/h)
    lvm          <- log(0.183 * 24);   label("Maximum Michaelis-Menten elimination rate Vm (mg/day)")           # Martinez 2019 Table 2 (0.183 mg/h; table footer 'mg.h/L' is a typo - text and dimensional analysis confirm mg/h)
    lkm          <- log(7.73);         label("Michaelis-Menten constant Km at reference FPCSK9 (mg/L)")          # Martinez 2019 Table 2 (FPCSK9 reference 72.9 ng/mL)
    llag         <- log(0.641 / 24);   label("SC absorption lag time (day)")                                    # Martinez 2019 Table 2 (0.641 h)
    logitfdepot  <- log(0.862 / (1 - 0.862)); label("Logit of SC bioavailability F (unitless; F_pop = 0.862)")  # Martinez 2019 Table 2 (typical F = 0.862)

    # Covariate effects (additive on CLL, additive on Km, power on V3; Martinez 2019 Table 2 and equations).
    # Additive slopes on CL are in L/h in the paper and converted to L/day (x24).
    e_wt_cl       <-  2.92e-4 * 24; label("Additive slope of (WT - 82.9 kg) on CLL (L/day per kg)")                 # Martinez 2019 Table 2 (2.92e-4 L/h/kg)
    e_statin_cl   <-  6.44e-3 * 24; label("Additive effect of STATIN on CLL (L/day when STATIN=1)")                 # Martinez 2019 Table 2 (6.44e-3 L/h)
    e_age_vp      <-  0.310;        label("Power exponent of AGE/60 on peripheral volume V3 (unitless)")            # Martinez 2019 Table 2
    e_fpcsk9_km   <- -0.541;        label("Additive slope of (FPCSK9/72.9) on Km (mg/L per unit of FPCSK9/72.9)")   # Martinez 2019 Table 2

    # Inter-individual variability. Martinez 2019 Table 2 reports omega^2 directly (NONMEM variance for
    # the exponential error model on CL, V2, V3, Km). The paper notes that LAG, Q, Ka, and Vm had no IIV.
    # The block V3/Km was estimated with a covariance parameter; Table 2 footnote d reports the correlation
    # coefficient r = -0.793, so cov(etalvp, etalkm) = -0.793 * sqrt(0.0735 * 0.298) = -0.11738.
    # F IIV is estimated in logit space (Table 2 footnote e) with omega^2_F = 1.060 (103%).
    etalcl                 ~ 0.232                                # Martinez 2019 Table 2 (CV 48.2%)
    etalvc                 ~ 0.589                                # Martinez 2019 Table 2 (CV 76.7%)
    etalvp + etalkm        ~ c(0.0735, -0.11738, 0.298)           # Martinez 2019 Table 2 (V3 CV 27.1%, Km CV 54.6%, r = -0.793)
    etalogitfdepot         ~ 1.060                                # Martinez 2019 Table 2 (logit-space IIV, 103%)

    # Residual error - Martinez 2019 Table 2 reports proportional SD = 0.259 (25.9%) and additive SD = 0.0465 mg/L.
    propSd <- 0.259;  label("Proportional residual error (SD, fraction)")   # Martinez 2019 Table 2 (theta8, 25.9%)
    addSd  <- 0.0465; label("Additive residual error (SD, mg/L)")            # Martinez 2019 Table 2 (theta9)
  })
  model({
    # Individual PK parameters. Martinez 2019 parameterizes the covariate effects on the
    # population *typical value* (additively on CLL and Km, as a power on V3) and then applies
    # exponential inter-individual variability to the typical value.

    # Linear clearance: additive covariate model for WT (relative to 82.9 kg) and STATIN.
    cl_tv <- exp(lcl) + e_wt_cl * (WT - 82.9) + e_statin_cl * STATIN
    cl    <- cl_tv * exp(etalcl)

    # Central volume (V2): no covariates.
    vc <- exp(lvc + etalvc)

    # Peripheral volume (V3): power-form age effect relative to 60 y.
    vp_tv <- exp(lvp) * (AGE / 60)^e_age_vp
    vp    <- vp_tv * exp(etalvp)

    # Intercompartmental clearance, absorption rate, and Vm: no IIV per the paper.
    q  <- exp(lq)
    ka <- exp(lka)
    vm <- exp(lvm)

    # Michaelis-Menten Km: additive linear-deviation effect of free PCSK9 (reference 72.9 ng/mL).
    km_tv <- exp(lkm) + e_fpcsk9_km * (FPCSK9 / 72.9)
    km    <- km_tv * exp(etalkm)

    # Bioavailability on logit scale with IIV (logit-space).
    logit_f <- logitfdepot + etalogitfdepot
    fdepot  <- 1 / (1 + exp(-logit_f))

    # Absorption lag time (no IIV).
    lag <- exp(llag)

    # Observation and state abbreviations.
    Cc <- central / vc

    # Two-compartment ODE with parallel linear (CLL) and Michaelis-Menten elimination from the
    # central compartment. MM term vm * Cc / (km + Cc) yields mg/day when vm is mg/day.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          vm * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # SC bioavailability and absorption lag apply to the depot. IV doses (phase I only in the
    # pooled dataset) bypass the depot via cmt=central on the dose event.
    f(depot)    <- fdepot
    alag(depot) <- lag

    Cc ~ add(addSd) + prop(propSd)
  })
}

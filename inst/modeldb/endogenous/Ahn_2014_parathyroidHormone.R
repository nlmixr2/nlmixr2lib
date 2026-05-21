Ahn_2014_parathyroidHormone <- function() {
  description <- "Semi-mechanistic indirect-response PD model of parathyroid hormone (PTH) suppression after oral calcium intake in healthy adults; absorbed but unobserved ionized calcium inhibits PTH secretion via an Emax (fixed to 1) negative-feedback term, with a parallel homeostatic indirect-response model for observed plasma ionized calcium."
  reference   <- "Ahn JE, Jeon S, Lee J, Han S, Yim DS. Modeling of the Parathyroid Hormone Response after Calcium Intake in Healthy Subjects. Korean J Physiol Pharmacol. 2014 Jun;18(3):217-223. doi:10.4196/kjpp.2014.18.3.217"
  vignette    <- "Ahn_2014_parathyroidHormone"
  units       <- list(
    time          = "hour",
    dosing        = "unit (normalized: depot at t=0 set to 1 for the thermal-water reference per Ahn 2014 p. 218)",
    concentration = "mmol/L (ionized Ca, observation Cc) and pg/mL (PTH, observation PTH)"
  )

  covariateData <- list(
    FORM_CACO3 = list(
      description        = "1 = calcium carbonate tablet (500 mg CaCO3 = 200 mg elemental Ca x 2 tablets, with 240 mL normal saline or 340 mL purified water); 0 = Geumjin thermal spring water (240 mL containing 400 mg elemental calcium followed by 100 mL purified water, the Ahn 2014 reference arm).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Geumjin thermal spring water; bioavailability anchor F = 1).",
      notes              = "Per-subject (treatment-arm) categorical indicator. The two CaCO3 arms in Ahn 2014 Methods (CaCO3 + saline and CaCO3 + water) were combined into a single CaCO3 treatment per the paper's Methods (p. 218): 'the 12 subjects that received calcium carbonate tablets were combined and treated as one treatment'. Encoded as a multiplicative effect on depot bioavailability per Ahn 2014 Results (Table 2 'Relative F1' = 1.98, 95% CI 1.06-2.90).",
      source_name        = "treatment (paper-narrative categorical: thermal spring water vs calcium carbonate)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24,
    n_studies      = 1,
    age_range      = "21-39 years",
    age_median     = "26 years",
    weight_range   = "55.1-79.3 kg",
    weight_median  = "68.5 kg (mean reported)",
    sex_female_pct = 8.33,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Healthy adults",
    dose_range     = "400 mg elemental calcium / day (240 mL Geumjin thermal spring water, n=12) or 500 mg calcium carbonate / day (two 200-mg calcium tablets with 240 mL normal saline or 340 mL purified water, n=12; the two CaCO3 sub-arms were pooled in the analysis per Methods p. 218). Days 1 and 7 received the full morning dose; Days 2-6 split the daily dose into twice-daily administrations.",
    regions        = "Korea (Seoul St. Mary's Hospital, The Catholic University of Korea)",
    notes          = "Randomized parallel clinical trial of Geumjin thermal spring water vs CaCO3 tablets (2:1:1 ratio across the three treatment arms; pooled to thermal water vs CaCO3 in the final analysis). 24 subjects produced 423 PTH and 423 ionized Ca observations across the 8-hour post-dose sampling window (0, 0.5, 1, 1.5, 2, 3, 4, 6, 8 hours) on Days 1 and 7 (Ahn 2014 Methods p. 218 and Results p. 219). Inclusion criteria required baseline 25(OH) vitamin D3 within 4.8-52.8 ng/mL and baseline calcium within 8.0-10.0 mg/dL; weight within +/- 20% of ideal body weight calculated as (height_cm - 100) * 0.9. Lifestyle covariates (smoking 7/17, alcohol 14/10, caffeine 14/10) were reported in Table 1 but not tested as model covariates -- only treatment formulation was tested (Methods p. 218)."
  )

  ini({
    # Structural parameters -- Ahn 2014 Table 2 final estimates.

    # Absorption rate of calcium from the depot (oral intake) into the unobserved Ca pool.
    lka       <- log(0.796)        ; label("First-order absorption rate constant for oral calcium (1/hr)")  # Ahn 2014 Table 2 (ka = 0.796 1/hr; 95% CI 0.438-1.15)

    # Indirect-response model for observed plasma ionized Ca: kin / kout turnover.
    lkin_ca   <- log(3.41)         ; label("Zero-order input rate for observed Ca, kin_ca (mmol/L/hr)")     # Ahn 2014 Table 2 (kin_ca = 3.41 mmol/L/hr; 95% CI 3.39-3.43)
    lkout_ca  <- fixed(log(2.86))  ; label("First-order output rate constant for Ca, kout_ca (1/hr)")        # Ahn 2014 Table 2 (kout_ca fixed at literature value 2.862 1/hr, Abraham et al. 2009 -- reference [6] of Ahn 2014)

    # Indirect-response model for PTH: kin_pth production with inhibition by unobserved Ca,
    # first-order kout_pth elimination.
    lkin_pth  <- log(21.6)         ; label("Zero-order input rate for PTH, kin_pth (pg/mL/hr)")               # Ahn 2014 Table 2 (kin_pth = 21.6 pg/mL/hr; 95% CI 12.4-30.8)
    lkout_pth <- log(0.849)        ; label("First-order output rate constant for PTH, kout_pth (1/hr)")       # Ahn 2014 Table 2 (kout_pth = 0.849 1/hr; 95% CI 0.513-1.19)

    # Inhibitory Emax model on PTH secretion (Ahn 2014 Methods p. 218; Eq. on the same page).
    lec50     <- log(0.158)        ; label("Half-maximal inhibitory concentration of unobserved Ca, EC50 (mmol/L)")  # Ahn 2014 Table 2 (EC50 = 0.158 mmol/L; 95% CI 0.0924-0.224)
    emax      <- fixed(1)          ; label("Maximum fractional PTH-secretion suppression, Emax (dimensionless)")     # Ahn 2014 Methods p. 218 + Results p. 220: 'It was also necessary to fix the Emax to 1, as the maximal calcium effect on PTH was unlikely to be observed with the calcium amount evaluated in this study.'

    # Relative bioavailability: thermal water = reference (F = 1 fixed), CaCO3 = 1.98.
    lfdepot              <- fixed(log(1))    ; label("Bioavailability anchor for thermal-water reference (log scale)")  # Ahn 2014 Methods p. 218 (paper: 'The bioavailability of calcium carbonate was estimated relative to that of thermal spring water', i.e. thermal water F fixed at 1)
    e_form_caco3_fdepot  <- log(1.98)        ; label("Log relative bioavailability of CaCO3 vs thermal water")          # Ahn 2014 Table 2 ('Relative F1' = 1.98; 95% CI 1.06-2.90; RSE 24%)

    # Inter-individual variability (log-normal; variance on the internal scale).
    # CV% relates via omega^2 = log(CV^2 + 1).
    etalka       ~ 0.406       # Ahn 2014 Table 2 (IIV variance in ka; sqrt(exp(0.406)-1) ~= 70.8% CV)
    etalkin_ca   ~ 0.000229    # Ahn 2014 Table 2 (IIV variance in kin_ca; CV ~= 1.5%. Published 95% CI in Table 2 is reported as '0.000814, 0.000377' which is order-inverted in the paper; bootstrap CI 0.0001-0.0004 confirms the magnitude.)
    etalkin_pth  ~ 0.0453      # Ahn 2014 Table 2 (IIV variance in kin_pth; sqrt(exp(0.0453)-1) ~= 21.4% CV, matching Ahn 2014 Results p. 220 narrative '21.4%')

    # NOTE: Ahn 2014 also reports inter-occasion variability (IOV) on kin_pth of 0.00969
    # (~9.84% CV between Day 1 and Day 7 occasions). IOV is not encoded as a separate
    # construct in the packaged model file; for occasion-aware simulation the user can add an
    # occasion-indexed eta on lkin_pth following the standard nlmixr2 IOV pattern.

    # Residual error -- Ahn 2014 Table 2.
    addSd      <- 0.0258   ; label("Additive residual SD on observed Ca, Cc (mmol/L)")           # Ahn 2014 Table 2 (SD_ca = 0.0258 mmol/L; 95% CI 0.0212-0.0304)
    propSd_PTH <- 0.21     ; label("Proportional residual SD on observed PTH (fraction)")        # Ahn 2014 Table 2 (CV_pth = 21.0%; 95% CI 18.8-23.2%)
  })

  model({
    # Individual structural parameters.
    ka       <- exp(lka + etalka)
    kin_ca   <- exp(lkin_ca + etalkin_ca)
    kout_ca  <- exp(lkout_ca)
    kin_pth  <- exp(lkin_pth + etalkin_pth)
    kout_pth <- exp(lkout_pth)
    ec50     <- exp(lec50)
    fdepot   <- exp(lfdepot + e_form_caco3_fdepot * FORM_CACO3)

    # Inhibitory effect of unobserved Ca on PTH secretion (Ahn 2014 Methods p. 218).
    inhib <- emax * ca_unobs / (ec50 + ca_unobs)

    # ODE system -- Ahn 2014 Methods p. 218 equations.
    # depot:    oral calcium dose, first-order absorption out of depot.
    # ca_unobs: net absorbed but unobserved Ca, driven by absorption from depot and
    #           eliminated at the same kout_ca as observed Ca (Ahn 2014 Results p. 220).
    # ca:       observed Ca, independent indirect-response (homeostasis); not coupled to
    #           the calcium dose because the dose does not noticeably perturb observed Ca.
    # pth:      observed PTH, indirect response with inhibition by ca_unobs.
    d/dt(depot)    <- -ka * depot
    d/dt(ca_unobs) <-  ka * depot - kout_ca * ca_unobs
    d/dt(ca)       <-  kin_ca - kout_ca * ca
    d/dt(pth)      <-  kin_pth * (1 - inhib) - kout_pth * pth

    # Bioavailability of the depot per treatment arm (Ahn 2014 Table 2 Relative F1).
    f(depot) <- fdepot

    # Steady-state initial conditions (Ahn 2014 Methods assumption that no circadian rhythm
    # or systemic difference in baseline Ca / PTH levels was present on Days 1 and 7 beyond
    # random variability; baselines are the drug-free steady states of the two homeostatic
    # indirect-response models).
    ca_unobs(0) <- 0
    ca(0)       <- kin_ca / kout_ca
    pth(0)      <- kin_pth / kout_pth

    # Observations.
    Cc  <- ca
    PTH <- pth

    Cc  ~ add(addSd)
    PTH ~ prop(propSd_PTH)
  })
}

Livio_2014_tobramycin <- function() {
  description <- "One-compartment population PK model with first-order absorption for systemic tobramycin released from an implanted calcium-sulfate bone-graft substitute (Osteoset T) in adults undergoing orthopedic surgery (Livio 2014); clearance equated to Cockcroft-Gault creatinine clearance under the assumption that absorbed tobramycin is exclusively eliminated by glomerular filtration, and absolute bioavailability differing between the 10 g (262 mg tobramycin) and 20 g (524 mg tobramycin) Osteoset T cast cohorts."
  reference <- "Livio F, Wahl P, Csajka C, Gautier E, Buclin T. Tobramycin exposure from active calcium sulfate bone graft substitute. BMC Pharmacol Toxicol. 2014;15:12. doi:10.1186/2050-6511-15-12"
  vignette <- "Livio_2014_tobramycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Raw Cockcroft-Gault CrCL in mL/min (NOT BSA-normalized to mL/min/1.73 m^2); computed from each patient's serum creatinine, body weight, and sex per Livio 2014 Methods. Linear interpolation was used for days without a fresh creatinine measurement. CL was equated to CrCL in the structural model under the paper's assumption that tobramycin absorbed from Osteoset T is exclusively eliminated by glomerular filtration. The typical clearance of 7.14 L/h reported in Table 1 corresponds to the population-mean CrCL of 119 mL/min via CL = CrCL x 60/1000. Encoded here as a power-1 effect on (CRCL/119); the reference 119 mL/min is the cohort mean (Livio 2014 Results / Table 1). Stored under canonical CRCL with `CLCR` (raw Cockcroft-Gault mL/min) as the documented source-name alias.",
      source_name        = "CLcr"
    ),
    DOSE = list(
      description        = "Subject's assigned single Osteoset T tobramycin dose at implantation (262 mg for 10 g Osteoset T cast, 524 mg for 20 g cast)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a binary indicator for the 20 g cast cohort: f(depot) <- exp(lfdepot + etalfdepot) * (1 + e_dose_524mg_fdepot * (DOSE == 524)) reproduces Livio 2014 Table 1's distinct bioavailability estimates of 0.63 (10 g cast, reference) and 0.32 (20 g cast). Each subject received a single Osteoset T implantation at time 0 so DOSE is constant per subject. The 4 % w/w tobramycin sulfate loading of Osteoset T combined with the tobramycin/tobramycin-sulfate salt factor of 0.655 (Methods) gives 262 mg tobramycin for a 10 g cast and 524 mg for a 20 g cast. Use case (a) of the canonical DOSE register (per-subject assigned dose level used as a stratified covariate); precedent: Castro-Surez 2020 nimotuzumab uses (DOSE == 50) on V1.",
      source_name        = "(derived from AMT; 262 mg = 10 g cast, 524 mg = 20 g cast)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 12L,
    n_studies       = 1L,
    age_range       = "19-82 years",
    age_mean        = "52 years (SD 20)",
    weight_range    = "53-116 kg",
    weight_mean     = "73 kg (SD 17)",
    crcl_range      = "34-288 mL/min (Cockcroft-Gault, raw mL/min)",
    crcl_mean       = "119 mL/min (SD 55)",
    sex_female_pct  = 41.7,
    race_ethnicity  = "Not reported (single-centre Swiss orthopaedic cohort, Fribourg Cantonal Hospital)",
    disease_state   = "Adults undergoing orthopaedic surgery for established bone / soft-tissue / prosthetic-joint infection (9 patients) or surgical-site infection prophylaxis (3 patients), with Osteoset T (4 % tobramycin sulfate calcium-sulfate bone-graft substitute) implanted as adjunct local antimicrobial therapy.",
    dose_range      = "Single Osteoset T implantation: 10 g cast containing 262 mg tobramycin (8 patients) or 20 g cast containing 524 mg tobramycin (4 patients).",
    implantation_sites = c(`tibia/fibula` = 6L, hip = 2L, calcaneum = 2L, femur = 1L, `lumbar spine` = 1L),
    regions         = "Single centre: Department of Orthopedic Surgery, Cantonal Hospital, Fribourg, Switzerland; data collected October 2006 - March 2008.",
    clinicaltrials  = "NCT01938417",
    n_observations  = "9 blood samples per patient at 3, 6, 12, 24, 48 hours and on days 3, 5, 7, 10 post-implantation (concentrations below the LOQ of 0.1 mg/L: first sample below LOQ set to LOQ/2 = 0.05 mg/L, subsequent BLQ samples dropped).",
    notes           = "Baseline demographics from Livio 2014 Results section first paragraph (paper does not present a dedicated Table 1 of demographics). Tourniquet release-time replaced implantation time as dose-event time when a tourniquet was used during the operation. No intravenous aminoglycoside co-medication. No wound drains."
  )

  ini({
    # Structural parameters - Livio 2014 Table 1 (final-model column).
    lka     <- log(0.0603); label("Absorption rate constant ka (1/h)")                                       # Livio 2014 Table 1: ka = 0.0603 h^-1 (s.e. 19 %)
    lvc     <- log(16.6);   label("Apparent volume of distribution V (L)")                                   # Livio 2014 Table 1: V = 16.6 L (s.e. 35 %)
    lfdepot <- log(0.63);   label("Bioavailability for the 10 g Osteoset T cast (reference cohort) (fraction)") # Livio 2014 Table 1: F (10 g) = 0.63 (s.e. 19 %)

    # CL equated to Cockcroft-Gault CrCL per Livio 2014 Methods ("As
    # aminoglycosides are known to be completely and exclusively eliminated
    # by glomerular filtration, tobramycin CL was equated to the creatinine
    # clearance (CLcr) value"). Encoded as a fixed power-1 effect on
    # (CRCL/119) anchored at the cohort-mean CL of 7.14 L/h (Table 1),
    # which corresponds to CRCL = 119 mL/min via CL [L/h] = CRCL [mL/min] *
    # 60 / 1000. Both anchors are wrapped in fixed() because the paper did
    # not estimate either parameter independently of the covariate.
    lcl       <- fixed(log(7.14)); label("Typical CL at the population-mean CrCL of 119 mL/min (L/h)")  # Livio 2014 Table 1: CL = 7.14 L/h (equated to CLcr)
    e_crcl_cl <- fixed(1.0);       label("Power exponent on (CRCL/119) for CL (unitless); fixed at 1 (CL equated to CrCL)") # Livio 2014 Methods: CL equated to CLcr

    # Bioavailability change for the 20 g cast cohort vs the 10 g reference.
    # Encoded as a multiplicative fractional change on F via the canonical
    # DOSE indicator (Castro-Surez 2020 nimotuzumab precedent: (DOSE == 50)
    # for cohort effect). Derived as F_20g / F_10g - 1 = 0.32 / 0.63 - 1 =
    # -0.4921 (the s.e. 19 % in Table 1 is for F_20g itself, not for the
    # relative change; the underlying NONMEM parameterization is not made
    # explicit in the publication).
    e_dose_524mg_fdepot <- -0.49206; label("Relative change in bioavailability for the 524 mg / 20 g cast cohort vs the 262 mg / 10 g reference (fraction)") # Livio 2014 Table 1: F (20 g) = 0.32 (0.32/0.63 - 1 = -0.49206)

    # Inter-individual variability (Livio 2014 Table 1, expressed as CV%).
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    #   V CV 89 % -> omega^2 = log(1 + 0.89^2) = 0.58339
    #   F CV 74 % -> omega^2 = log(1 + 0.74^2) = 0.43671
    etalvc      ~ 0.58339  # Livio 2014 Table 1: V CV 89 % (s.e. 74 %)
    etalfdepot  ~ 0.43671  # Livio 2014 Table 1: F CV 74 % (s.e. 72 %); paper Discussion notes the large V and F IIV likely incorporates some variability in CL and ka
    # No IIV on CL (CL is fixed-equated to CrCL per individual; Table 1 IIV column shows '-' for CL).
    # No IIV on ka (Table 1 IIV column shows '-' for ka).

    # Residual error - combined additive + proportional (Livio 2014 Methods
    # and Table 1).
    propSd <- 0.29;  label("Proportional residual error (fraction)")   # Livio 2014 Table 1: sigma_prop CV = 29 % (s.e. 50 %)
    addSd  <- 0.062; label("Additive residual error SD (mg/L)")        # Livio 2014 Table 1: sigma_add SD = 0.062 mg/L (s.e. 22 %)
  })

  model({
    # 1. Individual PK parameters
    ka <- exp(lka)                                                          # ka has no IIV (Table 1)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl) * (CRCL / 119)^e_crcl_cl                                 # CL equated to CrCL per Methods (no estimated IIV)

    # 2. Micro-constant
    kel <- cl / vc

    # 3. ODE system - one-compartment with first-order absorption from the
    #    implant depot. The "depot" represents the Osteoset T cast acting
    #    as a slow-release reservoir of tobramycin; first-order release
    #    rate ka was selected over a two-compartment alternative because
    #    the two-compartment model did not improve the fit (Methods /
    #    Results).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Bioavailability for the depot. 10 g cast (DOSE == 262) is the
    #    reference (F = 0.63); the 20 g cast cohort (DOSE == 524) shifts F
    #    multiplicatively via e_dose_524mg_fdepot. Per-subject IIV
    #    etalfdepot is shared across both cohorts (Table 1 reports a single
    #    F CV 74 %).
    f(depot) <- exp(lfdepot + etalfdepot) * (1 + e_dose_524mg_fdepot * (DOSE == 524))

    # 5. Observation and error. Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

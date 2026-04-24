Timmermann_2019_brodalumab <- function() {
  description <- "Two-compartment population PK model for brodalumab in adults with moderate-to-severe plaque psoriasis (Timmermann 2019), with first-order SC absorption, fixed bioavailability, and combined linear plus Michaelis-Menten (target-mediated) elimination from the central compartment."
  reference <- "Timmermann S, Hall A. Population pharmacokinetics of brodalumab in patients with moderate to severe plaque psoriasis. Basic Clin Pharmacol Toxicol. 2019;125(1):16-25. doi:10.1111/bcpt.13202"
  vignette <- "Timmermann_2019_brodalumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, V1 (central volume) and Vmax, each normalized as WT/90 kg per Timmermann 2019 Table 3 and final-model narrative (reference patient body-weight 90 kg).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 622L,
    n_observations = 7725L,
    n_studies      = 6L,
    age_range      = "18-75 years",
    age_median     = "46 years",
    weight_range   = "43-186 kg",
    weight_median  = "87.8 kg",
    sex_female_pct = 33,
    race_ethnicity = c(White = 93, Other = 7),
    disease_state  = "Moderate-to-severe plaque psoriasis (adults). Median baseline PASI 17.7 (range 8.8-60.6); phase III trials required baseline PASI >= 12.",
    dose_range     = "70-700 mg; SC at 70, 140, 210 or 280 mg (Q2W, Q2W+1, Q4W, Q8W) and IV at 700 mg (phase I single dose). Approved regimen is 210 mg SC at weeks 0, 1, 2 followed by 210 mg Q2W.",
    regions        = "Multi-regional (six pooled phase I/II/III trials: NCT00867100, NCT01937260, NCT00975637 with long-term extension NCT01101100, NCT01708590 AMAGINE-1, NCT01708603 AMAGINE-2, NCT01708629 AMAGINE-3).",
    notes          = "Baseline demographics from Timmermann 2019 Table 2 (PK analysis set of 622 psoriasis patients with rich non-trough sampling; 7725 quantifiable and 2508 BLQ concentrations). BMI median 29.6 kg/m^2 (16.7-66.1). Healthy volunteers were deliberately excluded from this analysis because IL-17RA expression (and therefore target-mediated disposition) differs between diseased and healthy populations. Race not tested as a covariate because 93% of patients were Caucasian. Binding ADA not tested due to low ADA incidence and transient nature of positive samples."
  )

  ini({
    # Structural PK parameters - Timmermann 2019 Table 3 final-model estimates.
    # Reference patient: body-weight 90 kg.
    lka <- log(0.300); label("Absorption rate Ka (1/day)")                              # Timmermann 2019 Table 3, Ka row
    lvc <- log(4.68);  label("Central volume V1 at 90 kg (L)")                          # Timmermann 2019 Table 3, V1 row
    lcl <- log(0.155); label("Linear clearance CL at 90 kg (L/day)")                    # Timmermann 2019 Table 3, CL row
    lvm <- log(6.07);  label("Maximum Michaelis-Menten elimination rate Vmax at 90 kg (mg/day)")  # Timmermann 2019 Table 3, Vmax row
    lq  <- log(0.328); label("Intercompartmental clearance Q (L/day)")                  # Timmermann 2019 Table 3, Q row
    lvp <- log(2.41);  label("Peripheral volume V2 (L)")                                # Timmermann 2019 Table 3, V2 row

    # Fixed parameters (Km carried from Endres 2014 previous analysis; F estimated in
    # the initial structural-model run then fixed before the covariate analysis per
    # Timmermann 2019 Methods).
    km      <- fixed(0.02); label("Michaelis-Menten constant Km (ug/mL)")               # Timmermann 2019 Table 3, Km row (fixed, same value as previous analysis)
    lfdepot <- log(0.548);  label("SC bioavailability F (fraction)")                    # Timmermann 2019 Table 3, F row (fixed at 54.8%)

    # Covariate exponents on reference 90 kg body weight - Timmermann 2019 Table 3.
    e_wt_vc <- 0.938; label("Power exponent of WT/90 on V1 (unitless)")                 # Timmermann 2019 Table 3: Power of weight on V1
    e_wt_cl <- 0.767; label("Power exponent of WT/90 on CL (unitless)")                 # Timmermann 2019 Table 3: Power of weight on CL
    e_wt_vm <- 0.769; label("Power exponent of WT/90 on Vmax (unitless)")               # Timmermann 2019 Table 3: Power of weight on Vmax

    # Inter-individual variability. Timmermann 2019 Table 3 reports IIV as %CV on the
    # linear-parameter scale, with log-normal (multiplicative exponential) random
    # effects per the Methods. Convert CV to NONMEM-style log-normal variance via
    # omega^2 = log(CV^2 + 1):
    #   Ka    CV 62.6% -> omega^2 = log(0.626^2 + 1) = 0.33065
    #   V1    CV 25.5% -> omega^2 = log(0.255^2 + 1) = 0.06300
    #   CL    CV 57.5% -> omega^2 = log(0.575^2 + 1) = 0.28565
    #   Vmax  CV 24.7% -> omega^2 = log(0.247^2 + 1) = 0.05922
    #   Q     CV 91.0% -> omega^2 = log(0.910^2 + 1) = 0.60328
    #   V2    CV 189%  -> omega^2 = log(1.890^2 + 1) = 1.51997
    # CL-V1 correlation 0.75 -> cov = 0.75 * sqrt(0.28565 * 0.06300) = 0.10061.
    etalcl + etalvc ~ c(0.28565, 0.10061, 0.06300)   # Timmermann 2019 Table 3: CL IIV 57.5% CV, V1 IIV 25.5% CV, CL-V1 correlation 0.75
    etalka ~ 0.33065                                  # Timmermann 2019 Table 3: Ka IIV 62.6% CV
    etalvm ~ 0.05922                                  # Timmermann 2019 Table 3: Vmax IIV 24.7% CV
    etalq  ~ 0.60328                                  # Timmermann 2019 Table 3: Q IIV 91% CV
    etalvp ~ 1.51997                                  # Timmermann 2019 Table 3: V2 IIV 189% CV

    # Combined additive + proportional residual error on the linear concentration
    # scale. Table 3 footnote: "Proportional and additive errors are given as %CV
    # and standard deviation."
    CcpropSd <- 0.355; label("Proportional residual error (fraction)")   # Timmermann 2019 Table 3: Proportional residual error 35.5% CV
    CcaddSd  <- 3.00;  label("Additive residual error (ug/mL)")          # Timmermann 2019 Table 3: Additive residual error 3.00 ug/mL (SD)
  })
  model({
    # Individual PK parameters with allometric weight scaling on CL, V1, Vmax.
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / 90)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 90)^e_wt_cl
    vm <- exp(lvm + etalvm) * (WT / 90)^e_wt_vm
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment PK with parallel linear and Michaelis-Menten elimination from
    # the central compartment. SC doses enter `depot` with bioavailability F; IV
    # doses bolus directly into `central` (user sets `cmt = central` on IV rows).
    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * central -
                          vm * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    f(depot) <- exp(lfdepot)

    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}

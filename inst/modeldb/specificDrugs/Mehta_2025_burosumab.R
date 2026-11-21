Mehta_2025_burosumab <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "linear elimination for subcutaneous burosumab (a human anti-FGF23 IgG1",
    "monoclonal antibody) in pediatric and adult patients with X-linked",
    "hypophosphatemia, coupled to a direct Emax model for absolute serum",
    "phosphorus. CL/F and V/F carry estimated allometric exponents on body",
    "weight (70 kg reference), and EC50 is a power function of each",
    "patient's own baseline serum phosphorus. This model generates the",
    "per-cycle average serum phosphorus that drives the companion fracture",
    "count model Mehta_2025_burosumab_fracture."
  )
  reference <- paste(
    "Mehta K, Storopoli J, Ramwani N, Quattrocchi E, Gobburu J, Weber T,",
    "Hruska MW, Marsteller D. Pharmacodynamic Exposure-Response Analysis of",
    "Fracture Count Data Following Treatment with Burosumab in Patients with",
    "XLH. J Clin Pharmacol. 2025;65(2):253. doi:10.1002/jcph.6140"
  )
  vignette <- "Mehta_2025_burosumab"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects on CL/F (exponent 0.87) and V/F (0.96), Table 2",
        "'CLWT' and 'VcWT'. The reference weight is not printed; 70 kg is",
        "used, as in the predecessor model (Lee 2022, shipped as",
        "Lee_2022_burosumab), whose time-varying weight was centred on 70 kg."
      ),
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "burosumab",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "burosumab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 278,
    n_studies = 9,
    age_range = "pediatric and adult patients (the paper does not tabulate ages for the PK/PD dataset)",
    weight_range = "not reported",
    disease_state = "X-linked hypophosphatemia (XLH)",
    dose_range = "0.05-1.2 mg/kg subcutaneous, single dose or every 2 or 4 weeks",
    regions = "multinational (NCT00830674, NCT01340482, NCT01571596, NCT02163577, NCT02750618, NCT02312687, NCT02526160, NCT02915705, NCT02537431)",
    notes = paste(
      "Methods 'Clinical Studies and Data' and Supplement S1: phase 1",
      "(KRN23-US-02, KRN23-INT-001, KRN23-INT-002), phase 2 (UX023-CL201,",
      "UX023-CL205, UX023-CL203) and phase 3 (UX023-CL303, UX023-CL301,",
      "UX023-CL304) studies. Concentrations below the LLOQ (36 or 50 ng/mL",
      "by assay) were set to missing. Fitted in Pumas 2.3.1. Demographics",
      "are not tabulated in this paper; the dataset updates the 277-patient",
      "dataset of Lee 2022."
    )
  )

  ini({
    # ---- Structural PK, Mehta 2025 Table 2 'PPK model parameters' (typical values at 70 kg) ----
    lka <- log(0.016)
    label("First-order absorption rate constant (1/h)") # Table 2: Ka = 0.016 1/h (RSE 8.8%); no IIV on Ka
    lcl <- log(0.012)
    label("Apparent clearance CL/F at 70 kg (L/h)") # Table 2: CL = 0.012 L/h (RSE 2.8%)
    lvc <- log(7.45)
    label("Apparent central volume V/F at 70 kg (L)") # Table 2: Vc = 7.45 L (RSE 3.0%)
    e_wt_cl <- 0.87
    label("Power exponent of (WT/70) on CL/F (unitless)") # Table 2: CLWT = 0.87 (RSE 3.5%)
    e_wt_vc <- 0.96
    label("Power exponent of (WT/70) on V/F (unitless)") # Table 2: VcWT = 0.96 (RSE 3.0%)

    # ---- PK IIV: Table 2 CV% on log-normal parameters, omega^2 = log(1 + CV^2) ----
    etalcl + etalvc ~ c(
      log(1 + 0.321^2), # CL
      0.46 * sqrt(log(1 + 0.321^2) * log(1 + 0.303^2)), # CL-Vc
      log(1 + 0.303^2) # Vc
    ) # Table 2: IIV CL 32.1%, IIV Vc 30.3%, correlation (CL, Vc) 0.46

    # ---- PK residual error ----
    addSd <- 39.32
    label("Additive residual error on serum burosumab (ng/mL)") # Table 2: sigma add = 39.32, unit printed 'ng/L' (read as ng/mL, the assay unit)
    propSd <- 0.233
    label("Proportional residual error on serum burosumab (fraction)") # Table 2: sigma prop = 23.3% (RSE 2.4%)

    # ---- PK-PD direct Emax model, Table 2 'PPKPD model parameters' ----
    le0 <- log(2.14)
    label("Baseline serum phosphorus BSL (mg/dL)") # Table 2: BSL = 2.14 mg/dL (RSE 1.0%)
    lec50 <- log(8270)
    label("EC50 coefficient: EC50 at a baseline serum phosphorus of 1 mg/dL (ng/mL)") # Table 2: EC50 = 8270 ng/mL (RSE 18.1%); uncentred power form, see model()
    lemax <- log(1.55)
    label("Maximal increase in serum phosphorus from baseline, Emax (mg/dL)") # Table 2: Emax = 1.55 mg/dL (RSE 4.7%)
    e_e0_ec50 <- -0.64
    label("Power exponent of the individual baseline serum phosphorus on EC50 (unitless)") # Table 2: dEC50dBPHOS = -0.64 (RSE 40%)

    # ---- PK-PD IIV: Table 2 'Between-patient variability (CV%)' ----
    # No IIV on EC50 (Discussion: 'IIV on EC50 was negligible and not in the final model').
    etale0 ~ log(1 + 0.140^2) # Table 2: IIV BSL = 14.0%
    etalemax ~ log(1 + 0.326^2) # Table 2: IIV Emax = 32.6%

    # ---- PK-PD residual error ----
    propSd_serum_phosphorus <- 0.133
    label("Proportional residual error on serum phosphorus (fraction)") # Table 2: residual variability sigma prop = 13.3% (RSE 3.5%)
  })

  model({
    # ---- 1. Individual PK parameters ----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    # ---- 2. One-compartment disposition with first-order SC absorption ----
    # Only subcutaneous data were modelled, so CL and Vc are apparent
    # (CL/F, V/F) and no bioavailability is estimated.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- 3. Serum burosumab concentration ----
    # central (mg) / vc (L) = mg/L = ug/mL; x1000 gives ng/mL, the unit of
    # the assay, the additive error and EC50.
    Cc <- 1000 * central / vc

    # ---- 4. Direct Emax model for serum phosphorus (Supplement S1 equation) ----
    # Phosphate_ij = Phosphate_i(0) + Emax * Conc_ij / (EC50 + Conc_ij).
    # EC50 = 8270 * BSL_i^-0.64: an uncentred power function of the
    # individual baseline e0. The paper does not print the form of
    # 'dEC50dBPHOS'; the uncentred power on the model's own baseline is the
    # form that reproduces Supplement Figure S1.1 (see the vignette).
    e0 <- exp(le0 + etale0)
    ec50 <- exp(lec50) * e0^e_e0_ec50
    emax <- exp(lemax + etalemax)
    serum_phosphorus <- e0 + emax * Cc / (ec50 + Cc)

    # ---- 5. Observations ----
    # Pumas combined error: SD = sqrt(add^2 + (prop * Cc)^2).
    Cc ~ add(addSd) + prop(propSd) + combined2()
    serum_phosphorus ~ prop(propSd_serum_phosphorus)
  })
}

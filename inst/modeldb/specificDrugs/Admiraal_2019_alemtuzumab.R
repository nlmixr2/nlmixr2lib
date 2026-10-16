Admiraal_2019_alemtuzumab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten elimination for intravenous alemtuzumab in children undergoing allogeneic hematopoietic cell transplantation, with a body-weight-dependent exponent on linear clearance (Admiraal 2019)"
  reference <- "Admiraal R, Jol-van der Zijde CM, Furtado Silva JM, Knibbe CAJ, Lankester AC, Boelens JJ, Hale G, Etuk A, Wilson M, Adams S, Veys P, van Kesteren C, Bredius RGM. Population Pharmacokinetics of Alemtuzumab (Campath) in Pediatric Hematopoietic Cell Transplantation: Towards Individualized Dosing to Improve Outcome. Clin Pharmacokinet. 2019;58(12):1609-1620. doi:10.1007/s40262-019-00782-0"
  vignette <- "Admiraal_2019_alemtuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline actual body weight. Reference (median) weight 17.25 kg as coded in the ESM control stream (MWT = WT/17.25); Table 2 footnote prints 17.3 kg and Table 1 17.2 kg. Enters V1 and Q as fixed power functions and linear CL as a body-weight-dependent exponent, (WT/17.25)^((e_wt_cl_bdecoef * WT)^e_wt_cl_bdeexp).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "alemtuzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "alemtuzumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 206,
    n_studies = 2,
    age_range = "0.2-19 years",
    age_median = "4.8 years",
    weight_range = "IQR 11-32 kg",
    weight_median = "17.2 kg",
    sex_female_pct = 33,
    race_ethnicity = "Not reported",
    disease_state = "Children receiving alemtuzumab serotherapy as part of conditioning for allogeneic hematopoietic cell transplantation (immune deficiency 52%, hematologic malignancy 25%, bone marrow failure 18%, metabolic disease 4%, benign hematology 1%)",
    dose_range = "Cumulative IV dose 0.2-1.5 mg/kg, usually 0.2 mg/kg/day for 4-5 consecutive days starting a median of 8 days before graft infusion; a few patients received 20 mg added to the graft bag",
    regions = "Netherlands (Leiden University Medical Center) and United Kingdom (Great Ormond Street Hospital)",
    n_transplants = 212,
    n_observations = 1146,
    notes = "Table 1 of Admiraal 2019: 139 patients from London and 67 (73 HCTs) from Leiden, treated 2003-2015; 28 infants under 1 year. Concentrations measured by an anti-idiotype ELISA (84%) or a quantitative flow-cytometry assay (16%)."
  )

  ini({
    # Structural parameters: Admiraal 2019 Table 2 final estimates; model form
    # from the ESM NONMEM control stream (run28, $PK and $DES).
    lcl <- log(0.25); label("Linear clearance at the 17.25 kg reference weight (L/day)") # Table 2 'CL pop (L/day)' = 0.25; ESM THETA(1)
    lvc <- log(2.13); label("Central volume of distribution at the 17.25 kg reference weight (L)") # Table 2 'V1,pop (L)' = 2.13; ESM THETA(2)
    lvp_ratio <- log(0.7); label("Log ratio of peripheral to central volume of distribution (unitless)") # Table 2 'V2,pop (factor of V1)' = 0.7; ESM 'V2 = V1 * THETA(3)'
    lq <- log(0.18); label("Intercompartmental clearance at the 17.25 kg reference weight (L/day)") # Table 2 'Q (L/day)' = 0.18; ESM THETA(4)
    lvmax <- log(0.42); label("Maximum elimination rate of the saturable pathway (mg/day)") # Table 2 'Vmax,pop (AU/day)' = 0.42; ESM THETA(5)
    lkm <- log(1.38); label("Michaelis-Menten constant of the saturable pathway (mg/L)") # Table 2 'Km,pop (AU/L)' = 1.38; ESM THETA(6)

    # Body-weight effects. Linear CL uses a body-weight-dependent exponent in
    # power-function form (Methods Eq. 7): exponent = (a * WT)^b.
    e_wt_cl_bdecoef <- 0.038; label("Coefficient a of the body-weight-dependent exponent on linear CL (1/kg)") # Table 2 'a' = 0.038; ESM THETA(7)
    e_wt_cl_bdeexp <- -0.79; label("Power b of the body-weight-dependent exponent on linear CL (unitless)") # Table 2 'b' = -0.79; ESM THETA(8)
    e_wt_vc <- 0.58; label("Power exponent of body weight on V1 (unitless)") # Table 2 'c' = 0.58; ESM THETA(9)
    e_wt_q <- 0.74; label("Power exponent of body weight on Q (unitless)") # Table 2 'd' = 0.74; ESM THETA(10)

    # IIV reported as %CV; omega^2 = log(CV^2 + 1). No IIV on Q or Vmax (ESM $PK).
    etalcl ~ 0.73314 # Table 2 'Interindividual variability on CL (%)' = 104
    etalvc ~ 0.33426 # Table 2 'Interindividual variability on V1 (%)' = 63
    etalkm ~ 1.06623 # Table 2 'Interindividual variability on Km (%)' = 138

    # Additive error on log-transformed concentrations (ESM $ERROR
    # 'Y = IPRED + ERR(1)', IPRED = LOG(F)), reported as proportional 34%.
    expSd <- 0.34; label("Residual error on the log scale (SD)") # Table 2 'Proportional residual error (%)' = 34
  })
  model({
    # Reference weight 17.25 kg from the ESM control stream (MWT = WT/17.25).
    wt_ratio <- WT / 17.25
    e_wt_cl_bde <- (e_wt_cl_bdecoef * WT)^e_wt_cl_bdeexp

    cl <- exp(lcl + etalcl) * wt_ratio^e_wt_cl_bde
    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc
    vp <- vc * exp(lvp_ratio)
    q <- exp(lq) * wt_ratio^e_wt_q
    vmax <- exp(lvmax)
    km <- exp(lkm + etalkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    d/dt(central) <- -kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc ~ lnorm(expSd)
  })
}

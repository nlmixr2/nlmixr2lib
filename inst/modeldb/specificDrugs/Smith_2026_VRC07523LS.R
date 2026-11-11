Smith_2026_VRC07523LS <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for the broadly neutralizing HIV-1 monoclonal antibody VRC07-523LS in adults with and without HIV, with clearance increased by time-varying plasma HIV-1 RNA through a hyperbolic (Emax, hill = 1) function (Smith 2026)."
  reference <- paste(
    "Smith NM, Ho BM, Bar KJ, Gama L, Dziubla G, Koup RA, Pensiero M,",
    "Serebryannyy L, Caskey M, Wilkin TJ, Wood TD, Morse GD, Venuto CS,",
    "Cha R, Ma Q. Enhanced Clearance of HIV-1 Broadly Neutralizing Antibody",
    "VRC07-523-LS During Viremia: Influences on Trial Design and Analysis.",
    "Clin Transl Sci. 2026;19(2):e70488. doi:10.1111/cts.70488.",
    sep = " "
  )
  vignette <- "Smith_2026_VRC07523LS"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    HIV_VLOAD = list(
      description = "Plasma HIV-1 RNA concentration at time t (time-varying)",
      units = "copies/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING, on the LINEAR scale. Smith 2026 Methods: 'HIV viral load was included as a time-dependent covariate to allow for the model to incorporate viral load effects dynamically'; Results equation CL_i(t) = CL_TV * (1 + Emax_i * HIV(t)^gamma / (EC50^gamma + HIV(t)^gamma)) * exp(eta_CL,i). Set to 0 for people without HIV and for virally suppressed people (the model then reduces to linear clearance CL_TV). Values below the assay limit of quantitation (< 50 copies/mL) have a negligible effect (50 / 16,596 < 0.3% of Emax), so recording either 0 or the LOQ gives the same prediction. The paper does not state its BLQ convention.",
      source_name = "HIV(t)"
    )
  )

  covariatesDataExcluded <- list(
    HIV_POS = list(
      description = "HIV serostatus (1 = person with HIV)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (without HIV)",
      notes = "Collected and screened (Smith 2026 Methods, 'HIV status, sex, height, weight, age, baseline serum creatinine, and baseline hematocrit were collected to test for covariate effects'); not in the final model (Table 1). The disease effect enters only through HIV_VLOAD."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened; not retained in the final model (Smith 2026 Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained in the final model (Smith 2026 Table 1). The clinical-trial simulations use flat mg doses with no weight effect."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained in the final model (Smith 2026 Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained in the final model (Smith 2026 Table 1)."
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained in the final model (Smith 2026 Table 1)."
    ),
    HCT = list(
      description = "Baseline hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained in the final model (Smith 2026 Table 1)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "VRC07523LS", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "VRC07523LS", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "VRC07523LS", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 2L,
    age_range = "Median 26 years (IQR 23.25-31.75)",
    weight_range = "Median 75.1 kg (IQR 61.5-85.05)",
    height_range = "Median 170.25 cm (IQR 164.25-176.75)",
    sex_female_pct = 44.1,
    disease_state = "Healthy adults without HIV (VRC 605 phase 1 dose escalation, NCT02840474; n = 25) and adults with chronic HIV-1 viremia not on antiretroviral therapy (NCT03015181; n = 9; viral load >= 500 copies/mL, CD4 >= 350 cells/mm3), who started antiretroviral therapy about 14 days after the antibody infusion.",
    dose_range = "Without HIV: 1, 5, 20 or 40 mg/kg IV single dose (n = 3 each); 5 mg/kg SC single dose (n = 3); 5 mg/kg SC or 20 mg/kg IV at day 0 and weeks 12 and 24 (n = 5 each). With HIV: 40 mg/kg IV single dose.",
    route = "IV (n = 26, 76.5%) and SC (n = 8, 23.5%)",
    viral_load = "Mean 35,264 copies/mL at antibody initiation in the 9 viremic participants; nadir at week 1, rebound to a peak of 4,410 copies/mL at week 3, suppression (mean 22.2 copies/mL) by weeks 7-8 after ART start.",
    labs = "Serum creatinine median 0.875 mg/dL (IQR 0.725-0.945); hematocrit median 42.1% (IQR 38.9-43.6)",
    sampling = "Median 17 (IQR 17-23) PK samples per participant; serum VRC07-523LS by anti-idiotype (clone 5C9) MSD electrochemiluminescence assay, limit of detection 0.5 ug/mL.",
    regions = "United States",
    notes = "Smith 2026 Supplementary Table S1 and Methods. Table S1 lists the dose groups as '1 mg / 2 mg / 20 mg / 40 mg' with n = 3 / 11 / 8 / 12; the Methods and the group sizes show these are 1 / 5 / 20 / 40 mg/kg."
  )

  ini({
    # Structural parameters: Smith 2026 Table 1 final estimates. The paper
    # reports volumes in mL and clearances in mL/day; they are carried here
    # in L and L/day so that mg / L gives mg/L.
    lka <- log(log(2) / 2.02); label("First-order SC absorption rate constant ka = ln(2) / Tabs (1/day)")          # Table 1 'T abs' = 2.02 day (first-order absorption half-life)
    logitfdepot <- logit(0.486); label("SC bioavailability F (fraction; logit scale)")                             # Table 1 'F bio' = 0.486
    lcl <- log(99.3 / 1000); label("Clearance without viral load CL_TV (L/day)")                                    # Table 1 'CL' = 99.3 mL/day
    lvc <- log(2080 / 1000); label("Central volume V1 (L)")                                                         # Table 1 'V1' = 2080 mL
    lvp <- log(2760 / 1000); label("Peripheral volume V2 (L)")                                                      # Table 1 'V2' = 2760 mL
    lq <- log(518 / 1000); label("Intercompartmental clearance CLD (L/day)")                                        # Table 1 'CLD' = 518 mL/day

    # Viral-load effect on clearance:
    # CL(t) = CL_TV * (1 + Emax * HIV(t)^gamma / (EC50^gamma + HIV(t)^gamma)).
    le_vload_cl <- log(5.46); label("Maximal fractional increase in CL from HIV-1 viral load, Emax (unitless)")      # Table 1 'CLmax' = 5.46; Results 'increase to 5.46-fold higher (i.e., CL = 641 mL/day)'
    lvload50_cl <- log(10^4.22); label("Plasma HIV-1 RNA at half-maximal CL increase, EC50 (copies/mL)")           # Table 1 'EC50' = 4.22 log10(copy/mL); Results 'EC50 ... 16,595 copies/mL'
    hill <- fixed(1); label("Hill coefficient of the viral-load effect on CL, gamma (unitless)")                   # Table 1 'gamma' = 1 (FIXED)

    # Inter-individual variability: Table 1 reports omega as SD, so the
    # variance is omega^2. CL and V2 are correlated (rho = 0.819):
    # cov = 0.819 * 0.296 * 0.390 = 0.0945454.
    etalka ~ 0.231361                       # Table 1 'omega Tabs' = 0.481 (SD); 0.481^2. Same magnitude on log(ka) as on log(Tabs).
    etalogitfdepot ~ 1.1236                 # Table 1 'omega F' = 1.06 (SD); 1.06^2, applied on the logit scale
    etalcl + etalvp ~ c(0.087616,
                        0.0945454, 0.1521)  # Table 1 'omega CL' = 0.296, 'omega V2' = 0.390, 'rho CL ~ V2' = 0.819
    etalvc ~ 0.145924                       # Table 1 'omega V1' = 0.382 (SD); 0.382^2
    etalq ~ 0.373321                        # Table 1 'omega CLD' = 0.611 (SD); 0.611^2
    etale_vload_cl ~ 0.200704               # Table 1 'omega CLMAX' = 0.448 (SD); 0.448^2

    # Residual error: combined additive + proportional (Table 1).
    addSd <- 0.523; label("Additive residual error (mg/L)")                                                         # Table 1 'sigma a' = 0.523 mg/L
    propSd <- 0.204; label("Proportional residual error (fraction)")                                                # Table 1 'sigma p' = 0.204
  })

  model({
    ka <- exp(lka + etalka)
    logitf_ind <- logitfdepot + etalogitfdepot
    fdepot <- expit(logitf_ind)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    # Viral-load-driven increase in clearance (Smith 2026 Results equation).
    e_vload_cl <- exp(le_vload_cl + etale_vload_cl)
    vload50_cl <- exp(lvload50_cl)
    f_vload_cl <- 1 + e_vload_cl * HIV_VLOAD^hill / (vload50_cl^hill + HIV_VLOAD^hill)
    cl <- exp(lcl + etalcl) * f_vload_cl

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # mg / L = mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

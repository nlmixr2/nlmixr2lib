Vet_2020_salbutamol <- function() {
  description <- "Enantiomer-resolved population PK model for continuous IV racemic salbutamol in children with refractory status asthmaticus in the PICU (Vet 2020). The racemate enters a dosing compartment that splits 50/50 into R- and S-salbutamol central compartments (k = 10 1/h, fixed); each enantiomer has its own two-compartment system with shared central volume, peripheral volume and intercompartmental clearance and an enantiomer-specific clearance. Allometric weight scaling to 70 kg with fixed exponents (0.75 on CL and Q, 1 on volumes). The fitted pre-infusion 'baseline dose' from prior inhalation is not encoded (see vignette)."
  reference <- paste(
    "Vet NJ, de Winter BCM, Koninckx M, Boeschoten SA, Boehmer ALM,",
    "Verhallen JT, Plotz FB, Vaessen-Verberne AA, van der Nagel BCH,",
    "Knibbe CAJ, Buysse CMP, de Wildt SN, Koch BCP, de Hoog M. (2020).",
    "Population Pharmacokinetics of Intravenous Salbutamol in Children with",
    "Refractory Status Asthmaticus. Clin Pharmacokinet 59(2):257-264.",
    "doi:10.1007/s40262-019-00811-y.",
    sep = " "
  )
  vignette <- "Vet_2020_salbutamol"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "racemic salbutamol", units = "ug", specimen = "plasma", verified = TRUE),
    central_r_enant = list(analyte = "R-salbutamol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_r_enant = list(analyte = "R-salbutamol", units = "ug", specimen = "plasma", verified = TRUE),
    central_s_enant = list(analyte = "S-salbutamol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_s_enant = list(analyte = "S-salbutamol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling to a 70-kg reference with exponents fixed at 0.75 on CL (both enantiomers) and Q, and 1 on Vc and Vp (Vet 2020 Section 2.3 and Section 3.2; estimating the exponents 'resulted in similar values and did not improve the model').",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the covariate analysis (Vet 2020 Section 2.3) and not retained ('No other covariates resulted in a significant improvement of the model', Section 3.2)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the covariate analysis (Vet 2020 Section 2.3) and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Tested in the covariate analysis (Vet 2020 Section 2.3) and not retained; the Discussion attributes this to the small cohort and the absence of renal dysfunction."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Tested in the covariate analysis (Vet 2020 Section 2.3) and not retained."
    ),
    BUN = list(
      description = "Urea concentration",
      units = "mmol/L",
      type = "continuous",
      notes = "The paper's 'urea' (Table 1, mmol/L). Tested in the covariate analysis (Vet 2020 Section 2.3) and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19L,
    n_studies = 1L,
    age_range = "0.8-15.3 years (median 4.9)",
    weight_range = "7.8-70 kg (median 18)",
    sex_female_pct = 26,
    race_ethnicity = NULL,
    disease_state = "Refractory status asthmaticus requiring continuous IV salbutamol in the paediatric intensive care unit; patients needing ECMO were excluded.",
    dose_range = "Continuous IV racemic salbutamol started at 0.1 ug/kg/min and titrated every 10 min to clinical response, up to 10 ug/kg/min (median 0.96, range 0.07-10 ug/kg/min); median infusion duration 22.5 h (range 7.6-192.5 h). All patients had received inhaled salbutamol before the IV infusion.",
    regions = "Netherlands (Erasmus MC-Sophia Children's Hospital PICU, Rotterdam)",
    baseline_dose = "Salbutamol already in the body from inhalation before the IV infusion was fitted as a per-patient 'baseline dose' of 448 ug (RSE 25%, IPV 99%, shrinkage 15%; Vet 2020 Table 2), split equally between R- and S-salbutamol. The paper does not state its timing, target compartment or weight scaling, so it is not encoded in ini(); supply it as a user dose record (see vignette).",
    notes = "Demographics from Vet 2020 Table 1 (n = 19; 14 male, 5 female; creatinine 29 (17-70) umol/L, ALT 15 (5-26) IU/L, urea 3.2 (1.6-4.8) mmol/L). 111 plasma samples for R- and S-salbutamol (median 3 per patient, range 1-22) after removing 6 samples contaminated by nebulised drug. Pre-infusion (inhalation only) plasma levels in 8 patients: R 24.5 (7-60) and S 49 (31-94) ug/L."
  )

  ini({
    # Split of the racemic dose from the dosing compartment into the two
    # enantiomer centrals (Vet 2020 Section 3.2 and Fig. 1).
    lka <- fixed(log(10)); label("First-order rate constant from the dosing compartment to the R- and S-salbutamol central compartments (1/h)") # Vet 2020 Section 3.2: 'k = 10 h-1 fixed'

    # Structural parameters for a 70-kg patient (Vet 2020 Table 2).
    lcl_r_enant <- log(16.3); label("R-salbutamol clearance for a 70-kg patient (L/h)") # Vet 2020 Table 2: CL R-salbutamol 16.3 L/h/70 kg (RSE 12.4%)
    lcl_s_enant <- log(8.82); label("S-salbutamol clearance for a 70-kg patient (L/h)") # Vet 2020 Table 2: CL S-salbutamol 8.82 L/h/70 kg (RSE 10.4%)
    lvc <- log(12.9); label("Central volume of each enantiomer for a 70-kg patient (L)") # Vet 2020 Table 2: Vc 12.9 L/70 kg (RSE 65%); identical for R and S (Table 2 footnote)
    lq <- log(20.4); label("Intercompartmental clearance of each enantiomer for a 70-kg patient (L/h)") # Vet 2020 Table 2: Q 20.4 L/h/70 kg (RSE 41.5%); identical for R and S
    lvp <- log(45.5); label("Peripheral volume of each enantiomer for a 70-kg patient (L)") # Vet 2020 Table 2: Vp 45.5 L/70 kg (RSE 9.6%); identical for R and S

    # Allometric exponents (Vet 2020 Section 2.3 and Section 3.2).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (both enantiomers) and Q (unitless)") # Vet 2020 Section 3.2: 'fixed at 0.75 for clearance and Q'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of body weight on Vc and Vp (unitless)") # Vet 2020 Section 3.2: 'at 1 for Vc and Vp'

    # IPV, exponential model; Table 2 'IPV (%)' read as CV, omega^2 = log(1 + CV^2).
    etalcl_r_enant ~ 0.16246 # Vet 2020 Table 2: IPV CL R-salbutamol 42%; log(1 + 0.42^2)
    etalcl_s_enant ~ 0.12831 # Vet 2020 Table 2: IPV CL S-salbutamol 37%; log(1 + 0.37^2)
    etalvc ~ 2.1856 # Vet 2020 Table 2: IPV Vc 281%; log(1 + 2.81^2); one eta shared by both enantiomers' Vc

    # Residual error: additive on log-transformed concentrations (Vet 2020
    # Section 2.3), one value for both enantiomers (Table 2). nlmixr2 needs a
    # separate residual parameter per endpoint, so the shared value is carried
    # twice.
    expSd_r_enant <- 0.351; label("Residual error SD on the log scale, R-salbutamol (shared value)") # Vet 2020 Table 2: residual error 0.351 (RSE 21.9%)
    expSd_s_enant <- 0.351; label("Residual error SD on the log scale, S-salbutamol (shared value)") # Vet 2020 Table 2: residual error 0.351 (RSE 21.9%)
  })

  model({
    # Allometric scaling to 70 kg (Vet 2020 Section 2.3).
    wt_ref <- 70

    ka <- exp(lka)
    cl_r_enant <- exp(lcl_r_enant + etalcl_r_enant) * (WT / wt_ref)^e_wt_cl
    cl_s_enant <- exp(lcl_s_enant + etalcl_s_enant) * (WT / wt_ref)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc_vp
    q <- exp(lq) * (WT / wt_ref)^e_wt_cl
    vp <- exp(lvp) * (WT / wt_ref)^e_wt_vc_vp

    kel_r_enant <- cl_r_enant / vc
    kel_s_enant <- cl_s_enant / vc
    k12 <- q / vc
    k21 <- q / vp

    # Racemic IV salbutamol is dosed into depot; half of the outflow goes to
    # each enantiomer (Vet 2020 Fig. 1: 'rapidly split into 50% R-salbutamol
    # and 50% S-salbutamol').
    d/dt(depot) <- -ka * depot
    d/dt(central_r_enant) <- 0.5 * ka * depot - kel_r_enant * central_r_enant - k12 * central_r_enant + k21 * peripheral1_r_enant
    d/dt(peripheral1_r_enant) <- k12 * central_r_enant - k21 * peripheral1_r_enant
    d/dt(central_s_enant) <- 0.5 * ka * depot - kel_s_enant * central_s_enant - k12 * central_s_enant + k21 * peripheral1_s_enant
    d/dt(peripheral1_s_enant) <- k12 * central_s_enant - k21 * peripheral1_s_enant

    Cc_r_enant <- central_r_enant / vc
    Cc_s_enant <- central_s_enant / vc

    Cc_r_enant ~ lnorm(expSd_r_enant)
    Cc_s_enant ~ lnorm(expSd_s_enant)
  })
}

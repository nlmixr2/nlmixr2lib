Ivaturi_2017_fludarabine <- function() {
  description <- paste(
    "Two-compartment population PK model for plasma f-ara-a (the circulating",
    "nucleoside of intravenous fludarabine phosphate) in children and young",
    "adults undergoing hematopoietic cell transplantation, with allometric",
    "weight scaling (reference 15 kg) and a linear creatinine-clearance effect",
    "on clearance, linked to an intracellular f-ara-ATP compartment in",
    "peripheral blood mononuclear cells whose uptake rate constant declines as",
    "a power of time after the first dose."
  )
  reference <- paste(
    "Ivaturi V, Dvorak CC, Chan D, Liu T, Cowan MJ, Wahlstrom J, Stricherz M,",
    "Jennissen C, Orchard PJ, Tolar J, Pai SY, Huang L, Aweeka F, Long-Boyle J.",
    "Pharmacokinetics and Model-Based Dosing to Optimize Fludarabine Therapy in",
    "Pediatric Hematopoietic Cell Transplant Recipients.",
    "Biol Blood Marrow Transplant. 2017;23(10):1701-1713.",
    "doi:10.1016/j.bbmt.2017.06.021"
  )
  vignette <- "Ivaturi_2017_fludarabine"

  # Doses are in mg of f-ara-a. Methods ('PopPK Analysis'): "An equivalent dose
  # of f-ara-a (molecular weight 285 g/mol) to that of the administered
  # monophosphate form (f-ara-AMP, molecular weight 365 g/mol) was calculated
  # and used for model building", so multiply a fludarabine phosphate dose by
  # 285/365 before dosing this model (the paper's inverse factor is 1.28).
  # Plasma f-ara-a is reported in ng/mL (Figure 1A); the intracellular
  # f-ara-ATP output `Cpbmc_faraatp` is in pmol per million PBMCs (Methods
  # 'PK Sampling': "Results for f-ara-ATP were normalized and presented as
  # pmoles per million cells").
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "fludarabine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fludarabine", units = "mg", specimen = "plasma", verified = TRUE),
    pbmc_faraatp = list(
      analyte = "fludarabine triphosphate",
      units = "pmol per million cells",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Actual (total) body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on every clearance and volume, normalised to 15 kg",
        "(Results equation: '(Weight / 15 kg)^0.75' on CL; Table 2 expresses CL",
        "and Vp per 15 kg). The paper states 'A body weight-based allometric",
        "model was added to all clearance and volume parameters' but prints",
        "only the CL exponent (0.75); the standard 0.75 is used for Q and 1 for",
        "Vc and Vp. Collected on each day of PK sampling (Methods 'Covariate",
        "Modeling'). Cohort median 15.7 kg, range 3-111 kg (Table 1)."
      ),
      source_name = "Weight"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance normalised to 1.73 m^2 body surface area,",
        "estimated by the Schwartz equation in children and by the",
        "Cockcroft-Gault equation with ideal body weight in young adults."
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on clearance centred at 100 mL/min/1.73 m^2: Results",
        "equation 'CL = CLpopulation x (Weight/15 kg)^0.75 x [1 + (CRCL -",
        "100 mL/min/1.73 m^2) x CRCLEFF]'. Serum creatinine reported as",
        "< 0.3 mg/dL was imputed as 0.15 mg/dL before applying the Schwartz",
        "equation (Methods 'Covariate Modeling'). The estimation data were",
        "capped at 150 mL/min/1.73 m^2 (Table 1 median 150, range 35-150; the",
        "follow-up refinement of this model, Brooks 2022 Pharmaceutics",
        "14:2462, states the cap explicitly), so cap CRCL at 150 before",
        "simulating; the cap is not applied inside model(). Table 1 prints the",
        "unit as 'mL/min/m2'; the Results equation and the Brooks 2022",
        "refinement both use mL/min/1.73 m^2."
      ),
      source_name = "CRCL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Evaluated as a continuous covariate on f-ara-a CL and not retained (Methods 'Covariate Modeling'; Results 'PopPK Analysis'). A renal-function maturation factor also did not improve the fit."
    ),
    HT = list(
      description = "Height.",
      units = "cm",
      type = "continuous",
      notes = "Evaluated on f-ara-a CL and not retained (Methods 'Covariate Modeling')."
    ),
    BSA = list(
      description = "Body surface area.",
      units = "m^2",
      type = "continuous",
      notes = "Evaluated on f-ara-a CL and not retained (Methods 'Covariate Modeling'). Cohort median 0.68 m^2, range 0.2-2.3 (Table 1)."
    ),
    BUN = list(
      description = "Blood urea nitrogen.",
      units = "mg/dL",
      type = "continuous",
      notes = "Evaluated on f-ara-a CL and not retained (Methods 'Covariate Modeling')."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/dL",
      type = "continuous",
      notes = "Evaluated on f-ara-a CL and not retained (Methods 'Covariate Modeling')."
    ),
    WBC = list(
      description = "White blood cell count.",
      units = "10^9/L",
      type = "continuous",
      notes = "Evaluated on f-ara-a CL and not retained (Methods 'Covariate Modeling')."
    ),
    NEUT = list(
      description = "Absolute neutrophil count.",
      units = "10^9/L",
      type = "continuous",
      notes = "Evaluated on f-ara-a CL and not retained (Methods 'Covariate Modeling')."
    ),
    SEXF = list(
      description = "Sex, 1 = female.",
      units = "(binary)",
      type = "binary",
      notes = "Gender was evaluated as a categorical covariate and not retained (Methods 'Covariate Modeling'). Conditioning regimen, indication for transplantation and study center were likewise tested and not retained (Results 'PopPK Analysis')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 140,
    n_studies = 1,
    age_range = "0.2-17.9 years",
    age_median = "5.0 years",
    weight_range = "3-111 kg",
    weight_median = "15.7 kg",
    sex_female_pct = 46,
    disease_state = paste(
      "children undergoing allogeneic hematopoietic cell transplantation (or",
      "autologous HCT with gene therapy) with fludarabine-containing",
      "conditioning, for hematologic malignancy (44%), primary immune",
      "deficiency (14%), hemoglobinopathy (6%), inherited metabolic disorder",
      "(16%), bone marrow failure (16%) or epidermolysis bullosa (4%)"
    ),
    renal_function = "creatinine clearance median 150 (range 35-150) mL/min/1.73 m^2; serum creatinine median 0.3 (0.14-1.89) mg/dL",
    dose_range = "fludarabine phosphate 12.5-40 mg/m^2 or 0.9-1.33 mg/kg IV over 30-60 min once daily for 3-5 days",
    regions = "United States (UCSF Benioff Children's Hospital, University of Minnesota Masonic Children's Hospital, Boston Children's Hospital)",
    notes = paste(
      "Prospective multicenter study NCT01316549, 2010-2015. 140 subjects",
      "completed PK assessments and were used for PK model building; Table 1",
      "demographics describe the 133 also included in the PK-PD analysis.",
      "840 plasma f-ara-a and 239 PBMC f-ara-ATP concentrations; f-ara-ATP",
      "was collected only at the UCSF site (n = 66). Sampling on two",
      "occasions: f-ara-a at 2, 4, 8 and 24 h after the start of the first",
      "sampled infusion and at 2 and 24 h after a later one; f-ara-ATP at 2",
      "and 24 h on both occasions. Below-quantification values were entered",
      "as reported."
    )
  )

  ini({
    # ---- f-ara-a structural parameters: Table 2, 'Final Model Results'.
    # All four are typical values at 15 kg body weight. Table 2 prints Vc as
    # 'L/kg' and Q as 'L/h/kg' but CL as 'L/h/15 kg' and Vp as 'L/15'; read
    # literally, Vc = 13.4 L/kg would give a 15 kg child a 201 L central
    # volume and f-ara-a peaks about 15-fold below the ~1000 ng/mL of
    # Figure 1A, so all four are taken per 15 kg.
    lcl <- log(3.1); label("Typical f-ara-a clearance at 15 kg and CRCL 100 mL/min/1.73 m^2 (L/h)") # Table 2 'Typical value for f-ara-a CL, L/h/15 kg' = 3.1 (RSE 3.3%); Results: '3.1L/hour is the typical value of f-ara-a CL'
    e_crcl_cl <- 0.006; label("Linear creatinine-clearance effect on f-ara-a clearance per mL/min/1.73 m^2 above 100 (1/(mL/min/1.73 m^2))") # Table 2 'Effect of creatinine clearance on f-ara-a CL' = .006 (RSE 9.5%); Results equation CRCLEFF
    lvc <- log(13.4); label("Central volume of distribution at 15 kg (L)") # Table 2 'Vc, L/kg' = 13.4 (RSE 6.3%); read per 15 kg
    lq <- log(2.2); label("Intercompartmental clearance at 15 kg (L/h)") # Table 2 'Intercompartmental CL, L/h/kg' = 2.2 (RSE 10.4%); read per 15 kg
    lvp <- log(13.4); label("Peripheral volume of distribution at 15 kg (L)") # Table 2 'Vp, L/15' = 13.4 (RSE 4.8%)

    # ---- Intracellular f-ara-ATP: Table 2. Kin multiplies the plasma f-ara-a
    # concentration in ng/mL and yields pmol per million cells per hour; that
    # pairing is the only one consistent with the observed f-ara-ATP
    # magnitudes in Table 3. The functional form of the 'Time effect on kin'
    # is not printed; it is encoded as a power of time after the first dose in
    # hours, Kin(t) = Kin * TAFD^-0.41, the only form among those tested that
    # reproduces the Table 3 observed medians and the within-day decline of
    # Figure 1B at the printed exponent (see the vignette).
    lkmet_faraatp <- log(0.005); label("f-ara-ATP formation rate constant Kin at 1 h after the first dose ((pmol/10^6 cells)/h per ng/mL)") # Table 2 'K in' = .005 (RSE 15.4%)
    e_time_kmet_faraatp <- -0.41; label("Power exponent of time after first dose on Kin (unitless)") # Table 2 'Time effect on kin' = -.41 (RSE 12.7%; bootstrap 95% CI -.52 to -.30)
    lkel_faraatp <- log(0.09); label("f-ara-ATP elimination rate constant Kout from PBMCs (1/h)") # Table 2 'K out' = .09 (RSE 7.2%)

    # ---- Between-subject variability: Table 2 labels the percentages %CV, but
    # they are omega x 100 (the SD of eta), so omega^2 = (P/100)^2. The
    # printed RSEs decide it: the RSE of omega implied by each row's bootstrap
    # 95% CI matches the printed RSE under this reading (Kin 11.3 vs 10.8
    # printed, Vc 22.2 vs 21.9, CL 7.6 vs 7.8), and misses under
    # omega^2 = log(1 + CV^2) (Kin 9.0, Vc 20.1, CL 7.3). Covariances from the
    # printed correlations as r * omega_i * omega_j.
    etalcl + etalvc + etalkmet_faraatp ~ c(
      0.100489,
      0.1044135, 0.146689,
      0.07529067, 0.1850694, 0.670761
    ) # Table 2 IIV CL 31.7, Vc 38.3, Kin 81.9 (omega x 100); correlations CL-Vc .86, CL-Kin .29, Vc-Kin .59

    # ---- Residual error: 'A separate additive error on the log scale for
    # f-ara-a and f-ara-ATP' (Results). Read on the same scale as the IIV rows,
    # the printed percentage is the log-scale SD x 100.
    expSd <- 0.248; label("Log-scale residual SD for plasma f-ara-a (unitless)") # Table 2 'f-ara-A residual unexplained variability' = 24.8
    expSd_Cpbmc_faraatp <- 0.715; label("Log-scale residual SD for PBMC f-ara-ATP (unitless)") # Table 2 'f-ara-ATP residual unexplained variability' = 71.5
  })

  model({
    # 1. Individual f-ara-a parameters. Results equation:
    #    CL = CLpop * (WT/15)^0.75 * [1 + (CRCL - 100) * CRCLEFF].
    cl <- exp(lcl + etalcl) * (WT / 15)^0.75 * (1 + (CRCL - 100) * e_crcl_cl)
    vc <- exp(lvc + etalvc) * (WT / 15)
    q <- exp(lq) * (WT / 15)^0.75
    vp <- exp(lvp) * (WT / 15)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 2. f-ara-ATP uptake. tafd(central) is the time since the first dose in
    #    hours; the floor keeps TAFD^-0.41 finite at the instant of the first
    #    dose, where the plasma concentration (and so the uptake flux) is zero.
    tfd <- tafd(central)
    tfdPos <- max(tfd, 1e-4)
    kmet_faraatp <- exp(lkmet_faraatp + etalkmet_faraatp) * tfdPos^e_time_kmet_faraatp
    kel_faraatp <- exp(lkel_faraatp)

    # 3. ODEs. Plasma f-ara-a in mg; the PBMC state is f-ara-ATP in pmol per
    #    million cells, driven by the plasma concentration. The uptake removes
    #    no mass from `central`: the paper models the PBMC pool as a separate
    #    per-cell quantity, and the flux is negligible on the plasma scale.
    Cc <- central / vc * 1000
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(pbmc_faraatp) <- kmet_faraatp * Cc - kel_faraatp * pbmc_faraatp

    # 4. Observations: f-ara-a in ng/mL; f-ara-ATP in pmol per million cells.
    Cpbmc_faraatp <- pbmc_faraatp

    Cc ~ lnorm(expSd)
    Cpbmc_faraatp ~ lnorm(expSd_Cpbmc_faraatp)
  })
}

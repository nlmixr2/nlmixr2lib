Cuffari_2016_mesalamine <- function() {
  description <- "Joint parent-metabolite population PK model for oral multimatrix mesalamine (5-aminosalicylic acid, 5-ASA) and its N-acetyl metabolite Ac-5-ASA in 52 children and adolescents (5-17 years, 18-82 kg) with ulcerative colitis (Cuffari 2016). Each dose is split between two parallel first-order absorption depots: a fraction F1 enters depot 1 (slow ka1, lag ALAG1) and the remainder depot 2 (faster ka2, absorbed after a further lag ALAG3 beyond ALAG1). 5-ASA is a one-compartment central pool cleared renally (CLR/F, to a cumulative urine compartment) and by acetylation (CLM/F) to Ac-5-ASA, which is itself one-compartment with renal (CLRM/F, to its own urine compartment) and non-renal (CLNRM/F) clearance. All clearances and volumes are apparent (divided by the unknown oral bioavailability) and allometrically scaled by body weight with fixed exponents 0.75 and 1 at a 70 kg reference."
  reference <- paste(
    "Cuffari C, Pierce D, Korczowski B, Fyderek K, Van Heusen H, Hossack S,",
    "Wan H, Edwards AYZ, Martin P.",
    "Randomized clinical trial: pharmacokinetics and safety of multimatrix",
    "mesalamine for treatment of pediatric ulcerative colitis.",
    "Drug Des Devel Ther. 2016;10:593-607.",
    "doi:10.2147/DDDT.S95316.",
    sep = " "
  )
  vignette <- "Cuffari_2016_mesalamine"
  dosing <- c("depot", "depot2")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "mesalamine", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "mesalamine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mesalamine", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "mesalamine", units = "mg", specimen = "urine", verified = TRUE),
    central_ac5asa = list(
      analyte = "N-acetyl-5-aminosalicylic acid",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    urine_ac5asa = list(analyte = "N-acetyl-5-aminosalicylic acid", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of every apparent clearance (exponent 0.75, fixed) and",
        "both apparent central volumes (exponent 1, fixed) at a 70 kg reference",
        "(Results 'Population pharmacokinetics'; Table S2 rows 'CLR/F~WT' ...",
        "'VcM/F~WT' and footnote a). Randomisation was stratified by body weight",
        "18-24, 25-49 and 50-82 kg.",
        sep = " "
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Named as a potential explanatory covariate (Methods 'Population pharmacokinetic analysis'); the final model (Table S2) carries no age effect."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Named as a potential explanatory covariate (Methods 'Population pharmacokinetic analysis'); the final model (Table S2) carries no sex effect."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 52L,
    n_studies = 1L,
    age_range = "5-17 years (children 5-12, adolescents 13-17)",
    age_mean = "13.3 years (SD 3.06); 14.0 / 13.9 / 10.6 years in the 30 / 60 / 100 mg/kg/day arms (Table 1)",
    weight_range = "18-82 kg",
    bmi_mean = "19.44 kg/m^2 (SD 3.308) (Table 1)",
    sex_female_pct = 57.7,
    race_ethnicity = "White 98.1%, Asian 1.9%; Hispanic or Latino 1.9% (Table 1)",
    disease_state = "Ulcerative colitis diagnosed at least 3 months before the first dose (not confined to the rectum); median 1.83 years since diagnosis (Table 1).",
    dose_range = "Multimatrix mesalamine 30, 60 or 100 mg/kg once daily each morning for 7 days, capped at 4800 mg/day; total daily doses 900-4800 mg using 300, 600 and 1200 mg tablets. Arms: 21 / 22 / 9 subjects.",
    regions = "Multicentre: United States, Poland and Slovakia (the ethics-committee list also names sites in Australia and the United Kingdom).",
    notes = paste(
      "Phase I randomised open-label three-arm study (NCT01130844). Plasma 5-ASA and",
      "Ac-5-ASA sampled pre-dose on days 5 and 6 and at 0, 2, 4, 6, 9, 12, 16 and 24 h",
      "on day 7; 0-24 h urine collected on day 7. LC-MS/MS assay 5-5000 ng/mL in",
      "plasma. Model fitted in NONMEM 7.2 (importance-sampling EM with mu-referencing);",
      "covariate model built on an interim cut of 40 subjects, final model refitted to",
      "all 52.",
      sep = " "
    )
  )

  ini({
    # Absorption (Table S2; Results 'Population pharmacokinetics'). The paper's
    # depot 1 and depot 3 (NONMEM compartments 1 and 3, Figure S1) are depot and
    # depot2 here.
    lka <- log(0.0334); label("Absorption rate constant from depot 1, ka1 (1/h)") # Table S2: Ka1 = 0.0334 1/h (RSE 7.26%)
    lka2 <- log(0.273); label("Absorption rate constant from depot 2 (paper depot 3), ka3 (1/h)") # Table S2: Ka3 = 0.273 1/h (RSE 19.5%)
    ltlag <- log(5.10); label("Absorption lag time of depot 1, ALAG1 (h)") # Table S2: ALAG1 = 5.10 h (RSE 5.34%)
    ldtlag2 <- log(15.0); label("Additional lag of depot 2 beyond ALAG1, ALAG3 (h)") # Table S2: ALAG3 = 15.0 h (RSE 1.23%); footnote 'absorption lag time from depot 3 in addition to the lag time from depot 1'
    logitfrel <- log(0.734 / (1 - 0.734)); label("Logit of the fraction of the dose absorbed from depot 1, F1 (unitless)") # Table S2: F1 = 0.734 (RSE 22.3%); logit(0.734) = 1.0157

    # 5-ASA disposition, 70 kg reference (Table S2 footnote a)
    lcl_renal <- log(1.15); label("Apparent renal clearance of 5-ASA, CLR/F (L/h)") # Table S2: CLR/F = 1.15 L/h (RSE 45.9%)
    lcl_met <- log(85.6); label("Apparent metabolic clearance of 5-ASA to Ac-5-ASA, CLM/F (L/h)") # Table S2: CLM/F = 85.6 L/h (RSE 1.38%)
    lvc <- log(109); label("Apparent central volume of 5-ASA, Vc/F (L)") # Table S2: Vc/F = 109 L (RSE 4.78%)

    # Ac-5-ASA disposition, 70 kg reference
    lcl_renal_ac5asa <- log(2.54); label("Apparent renal clearance of Ac-5-ASA, CLRM/F (L/h)") # Table S2: CLRM/F = 2.54 L/h (RSE 6.34%)
    lcl_nonren_ac5asa <- log(74.4); label("Apparent non-renal clearance of Ac-5-ASA, CLNRM/F (L/h)") # Table S2: CLNRM/F = 74.4 L/h (RSE 1.31%)
    lvc_ac5asa <- log(7.10); label("Apparent central volume of Ac-5-ASA, VcM/F (L)") # Table S2: VcM/F = 7.10 L (RSE 7.55%)

    # Allometric exponents (Table S2: every clearance ~WT 0.75 FIX, both
    # volumes ~WT 1 FIX)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on every apparent clearance (unitless)") # Table S2: CLR/F, CLM/F, CLRM/F, CLNRM/F ~WT = 0.75 FIX
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on both apparent central volumes (unitless)") # Table S2: Vc/F, VcM/F ~WT = 1 FIX

    # Inter-individual variability: Table S2 omega^2 column (the CV% column is
    # sqrt(omega^2) x 100 for every row).
    etalcl_renal ~ 0.181 # Table S2: omega2 CLR = 0.181 (CV 42.5%)
    etalcl_met ~ 0.115 # Table S2: omega2 CLM = 0.115 (CV 33.9%)
    etalvc ~ 0.412 # Table S2: omega2 Vc = 0.412 (CV 64.2%)
    etalcl_renal_ac5asa ~ 0.133 # Table S2: omega2 CLRM = 0.133 (CV 36.5%)
    etalcl_nonren_ac5asa ~ 0.0823 # Table S2: omega2 CLNRM = 0.0823 (CV 28.7%)
    etalvc_ac5asa ~ fixed(0.0225) # Table S2: omega2 VcM = 0.0225 FIX (CV 15.0%)
    etalka ~ 2.72 # Table S2: omega2 Ka1 = 2.72 (CV 165%)
    etalka2 ~ 1.87 # Table S2: omega2 Ka3 = 1.87 (CV 137%)
    etaltlag ~ 0.256 # Table S2: omega2 ALAG1 = 0.256 (CV 50.6%)
    etaldtlag2 ~ fixed(0.0225) # Table S2: omega2 ALAG3 = 0.0225 FIX (CV 15.0%)
    etalogitfrel ~ 0.662 # Table S2: omega2 F1 = 0.662 (CV 81.4%); applied on the logit scale, see vignette Assumptions

    # Residual error: proportional, Table S2 sigma^2 -> SD = sqrt(sigma^2)
    propSd <- 0.36878; label("Proportional residual SD, plasma 5-ASA (fraction)") # Table S2: sigma2 prop plasma 5-ASA = 0.136 (CV 36.9%); sqrt(0.136) = 0.36878
    propSd_ac5asa <- 0.28792; label("Proportional residual SD, plasma Ac-5-ASA (fraction)") # Table S2: sigma2 prop plasma Ac-5-ASA = 0.0829 (CV 28.8%); sqrt(0.0829) = 0.28792
    propSd_urine <- fixed(0.06); label("Proportional residual SD, urine 5-ASA amount (fraction)") # Table S2: sigma2 prop urine 5-ASA = 0.00360 FIX (CV 6.00%); sqrt(0.0036) = 0.06
    propSd_urine_ac5asa <- fixed(0.06); label("Proportional residual SD, urine Ac-5-ASA amount (fraction)") # Table S2: sigma2 prop urine Ac-5-ASA = 0.00360 FIX (CV 6.00%); sqrt(0.0036) = 0.06
  })

  model({
    # ---- Absorption ------------------------------------------------
    ka <- exp(lka + etalka)
    ka2 <- exp(lka2 + etalka2)
    tlag <- exp(ltlag + etaltlag)
    dtlag2 <- exp(ldtlag2 + etaldtlag2)
    tlag2 <- tlag + dtlag2
    logitfrel_ind <- logitfrel + etalogitfrel
    frel <- expit(logitfrel_ind)

    # ---- Disposition (allometric, 70 kg reference) -----------------
    cl_renal <- exp(lcl_renal + etalcl_renal) * (WT / 70)^e_wt_cl
    cl_met <- exp(lcl_met + etalcl_met) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    cl_renal_ac5asa <- exp(lcl_renal_ac5asa + etalcl_renal_ac5asa) * (WT / 70)^e_wt_cl
    cl_nonren_ac5asa <- exp(lcl_nonren_ac5asa + etalcl_nonren_ac5asa) * (WT / 70)^e_wt_cl
    vc_ac5asa <- exp(lvc_ac5asa + etalvc_ac5asa) * (WT / 70)^e_wt_vc

    # Ac-5-ASA states are carried in mg of Ac-5-ASA. Each mg of 5-ASA
    # acetylated forms 195.15 / 153.14 mg of Ac-5-ASA (molecular weights from
    # Methods 'Non-compartmental pharmacokinetic analysis').
    mwr_ac5asa <- 195.15 / 153.14

    kel_renal <- cl_renal / vc
    kel_met <- cl_met / vc
    kel_renal_ac5asa <- cl_renal_ac5asa / vc_ac5asa
    kel_nonren_ac5asa <- cl_nonren_ac5asa / vc_ac5asa

    # ---- ODE system (Figure S1) ------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - (kel_renal + kel_met) * central
    d/dt(urine) <- kel_renal * central
    d/dt(central_ac5asa) <- mwr_ac5asa * kel_met * central -
      (kel_renal_ac5asa + kel_nonren_ac5asa) * central_ac5asa
    d/dt(urine_ac5asa) <- kel_renal_ac5asa * central_ac5asa

    # Every dose is given to both depots; F1 and 1 - F1 partition it, so the
    # total bioavailability is 1 and all clearances and volumes are apparent.
    f(depot) <- frel
    f(depot2) <- 1 - frel
    alag(depot) <- tlag
    alag(depot2) <- tlag2

    # ---- Observations ----------------------------------------------
    # mg / L x 1000 = ng/mL (= ug/L, the units of the paper's Table 2 and Figure 2)
    Cc <- 1000 * central / vc
    Cc_ac5asa <- 1000 * central_ac5asa / vc_ac5asa

    Cc ~ prop(propSd)
    Cc_ac5asa ~ prop(propSd_ac5asa)
    urine ~ prop(propSd_urine)
    urine_ac5asa ~ prop(propSd_urine_ac5asa)
  })
}

Lim_2024_dihydroxyflavone_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Two-compartment plasma pharmacokinetic model for",
    "7,8-dihydroxyflavone (7,8-DHF, tropoflavin) in postnatal day 9 C57BL/6",
    "mice after a single 5 mg/kg intraperitoneal injection given 10 min after",
    "an experimentally induced unilateral (left-sided) hypoxic-ischemic (HI)",
    "brain injury or a sham operation, with the intraperitoneal dose treated as",
    "an immediate (intravenous-like) input. HI injury lowers clearance",
    "(CL = theta1 * exp(theta2 * HIE_POST), a 21% reduction); sex was tested",
    "and not retained. Right and left brain-hemisphere concentrations follow a",
    "sequentially fitted perfusion-limited distribution model driven by the",
    "plasma concentration, with the left-hemisphere partition coefficient",
    "expressed as the right-hemisphere Kp times a left-to-right ratio (Lim",
    "2024). All parameters are apparent and weight-normalised (per kg body",
    "weight) because intraperitoneal bioavailability is unknown.",
    sep = " "
  )
  reference <- paste(
    "Lim SY, Scarlett C, Yapici S, Ferrazzano P, Cengiz P. Pharmacokinetics of",
    "7,8-dihydroxyflavone in neonatal mice with hypoxia-ischemia related brain",
    "injury. Front Pharmacol. 2024;15:1508696. doi:10.3389/fphar.2024.1508696.",
    "Final parameter estimates from Table 1; model equations from Section 2.4;",
    "base-model and covariate-search tables from Supplementary Tables S2-S3.",
    sep = " "
  )
  vignette <- "Lim_2024_dihydroxyflavone_mouse"
  units <- list(time = "min", dosing = "mg/kg", concentration = "ng/mL")

  # Brain-hemisphere states are not canonical compartments: the paper resolves
  # the right (contralateral) and left (ipsilateral to the carotid ligation)
  # hemispheres as two separately sampled tissues.
  paper_specific_compartments <- c("brain_right", "brain_left")

  # Amounts are per kilogram of body weight because the dose and every fitted
  # volume and clearance are weight-normalised (mg/kg, mL/kg, mL/min/kg).
  compartmentData <- list(
    central = list(
      analyte = "7,8-dihydroxyflavone",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "7,8-dihydroxyflavone",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    brain_right = list(
      analyte = "7,8-dihydroxyflavone",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    brain_left = list(
      analyte = "7,8-dihydroxyflavone",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    HIE_POST = list(
      description = "Post hypoxic-ischemic event indicator (1 = mouse underwent the left-sided hypoxia-ischemia procedure; 0 = sham-operated)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sham-operated)",
      notes = "Lim 2024 Section 2.2: left common carotid artery electro-cauterisation followed by 50 min of 10% oxygen at 36 degC; sham mice had the same procedure without cauterisation or hypoxia. 7,8-DHF was injected 10 min after hypoxia ended, so HIE_POST is constant over the whole PK observation window (1 for HI mice, 0 for sham). Enters as CL = theta1 * exp(theta2 * HIE_POST) (Table 1 footnote b).",
      source_name = "surgical treatment (HI versus sham)"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested on CL in the stepwise covariate search and not added (Supplementary Table S3: delta -2LL = -1.64 alone, -1.83 on top of surgery, below the 6.64 forward threshold); also tested on Kp,R and Rkp and not significant (Section 3.2)."
    )
  )

  population <- list(
    species = "mouse (C57BL/6, postnatal day 9 neonates)",
    n_subjects = 80L,
    n_studies = 1L,
    age_range = "postnatal day 9",
    sex_female_pct = 50,
    disease_state = paste(
      "Unilateral left-sided hypoxic-ischemic brain injury (Vannucci model:",
      "left common carotid artery cauterisation then 50 min at 10% oxygen) or",
      "sham operation; 40 HI and 40 sham mice.",
      sep = " "
    ),
    dose_range = "7,8-DHF 5 mg/kg single intraperitoneal injection, 10 min after hypoxia",
    notes = paste(
      "Section 2.2.1: destructive sampling at 10, 30, 60 and 180 min after",
      "dosing, with 5 female and 5 male HI mice and 5 female and 5 male sham",
      "mice per time point (one plasma sample plus right and left brain",
      "hemispheres per mouse). Body weights were not reported. Plasma LLOQ",
      "0.75 ng/mL; brain LLOQ 2.25 ng/g (Section 3.1); BLQ values were",
      "censored in Phoenix NLME 8.3 (Laplacian).",
      sep = " "
    )
  )

  ini({
    # Plasma PK (Table 1, 'Plasma PK'; apparent, weight-normalised values)
    lcl <- log(270.1)
    label("Clearance in sham-operated mice, theta1 (mL/min/kg)") # Table 1 'theta1' = 270.1 (RSE 17.7%)
    lvc <- log(1591)
    label("Central volume of distribution V1 (mL/kg)") # Table 1 'V1 (mL/kg)' = 1,591 (RSE 48.3%)
    lq <- log(90.50)
    label("Distributional clearance CLd (mL/min/kg)") # Table 1 'CLd (mL/min/kg)' = 90.50 (RSE 29.0%)
    lvp <- log(5802)
    label("Peripheral volume of distribution V2 (mL/kg)") # Table 1 'V2 (mL/kg)' = 5,802 (RSE 34.5%)

    e_hie_post_cl <- -0.236
    label("Log-scale effect of HI injury on CL, theta2 (unitless)") # Table 1 'theta2' = -0.236 (RSE 18.8%); footnote b: HI CL = 270.1 * exp(-0.236) = 213.3

    # Brain PK (Table 1, 'Brain PK'; sequential fit driven by plasma C1)
    lv_brain_hemi <- fixed(log(25))
    label("Volume of each brain hemisphere Vb,R = Vb,L (mL/kg)") # Table 1 'Vb,R and Vb,L' = 2.5% of body weight, held constant (25 g/kg, density 1 g/mL)
    lqspec_brain <- fixed(log(0.54))
    label("Cerebral blood flow per gram of brain (mL/min/g)") # Table 1 'Qb,R and Qb,L' = 54 mL/min/100 g brain, held constant (Ek 2015)
    lkp_brain_right <- log(0.086)
    label("Right brain hemisphere-to-plasma partition coefficient Kp,R (unitless)") # Table 1 'Kp,R' = 0.086 (RSE 5.37%)
    lrkp <- log(1.158)
    label("Left-to-right brain partition coefficient ratio Rkp (unitless)") # Table 1 'Rkp' = 1.158 (RSE 3.17%); footnote c: Kp,L = 0.099

    # Interindividual variability (exponential model; Table 1 reports %CV,
    # converted with omega^2 = log(1 + CV^2))
    etalcl ~ 0.04039 # Table 1 'IIV CL (%)' = 20.3 -> log(1 + 0.203^2)
    etalkp_brain_right ~ 0.19271 # Table 1 'IIV Kp,R' = 46.1 -> log(1 + 0.461^2)
    etalrkp ~ 0.04685 # Table 1 'IIV Rkp' = 21.9 -> log(1 + 0.219^2)

    # Residual error, held at the bioanalytical QC variability (Section 2.4)
    propSd <- fixed(0.05)
    label("Proportional residual error on plasma concentration (fraction)") # Table 1 plasma 'Proportional error (%)' = 5, held constant
    propSd_Cbrain_right <- fixed(0.10)
    label("Proportional residual error on right-hemisphere concentration (fraction)") # Table 1 brain 'Proportional error (%)' = 10, held constant
    propSd_Cbrain_left <- fixed(0.10)
    label("Proportional residual error on left-hemisphere concentration (fraction)") # Table 1 brain 'Proportional error (%)' = 10, held constant
  })

  model({
    # Individual parameters (Section 2.4; Table 1 footnote b)
    cl <- exp(lcl + e_hie_post_cl * HIE_POST + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    v_brain <- exp(lv_brain_hemi)
    q_brain <- exp(lqspec_brain) * v_brain
    kp_brain_right <- exp(lkp_brain_right + etalkp_brain_right)
    kp_brain_left <- kp_brain_right * exp(lrkp + etalrkp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Plasma: V1 dC1/dt = input - CLd (C1 - C2) - CL C1; V2 dC2/dt = CLd (C1 - C2).
    # The intraperitoneal dose goes straight into central (immediate
    # absorption, Section 2.4).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Brain: Vb dCb/dt = Qb (C1 - Cb / Kp) per hemisphere. The brain model was
    # fitted sequentially with plasma C1 as a forcing function, so brain uptake
    # is not subtracted from central.
    cp_mass <- central / vc
    d/dt(brain_right) <- q_brain *
      (cp_mass - brain_right / v_brain / kp_brain_right)
    d/dt(brain_left) <- q_brain *
      (cp_mass - brain_left / v_brain / kp_brain_left)

    # mg/mL -> ng/mL (plasma) and mg/g -> ng/g (brain tissue, 1 mL = 1 g)
    Cc <- 1e6 * cp_mass
    Cbrain_right <- 1e6 * brain_right / v_brain
    Cbrain_left <- 1e6 * brain_left / v_brain

    Cc ~ prop(propSd)
    Cbrain_right ~ prop(propSd_Cbrain_right)
    Cbrain_left ~ prop(propSd_Cbrain_left)
  })
}

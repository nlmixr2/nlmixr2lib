Yusuff_2026_voriconazole <- function() {
  description <- "One-compartment population PK model of intravenous voriconazole in critically ill adults receiving extracorporeal membrane oxygenation (ECMO), with time-varying total clearance: an intrinsic clearance that rises logistically in time from an early to a late level, plus a parallel apparent ECMO-circuit sequestration clearance that decays mono-exponentially from treatment start. Volume scales linearly with body weight; CYP2C19 intermediate/poor metabolizer status lowers the late intrinsic clearance (Yusuff 2026)"
  reference <- "Yusuff H, Gadsby J, Isgro G, Zochios V, Jenkins D, Cooke S, Mulla H. Time-varying voriconazole clearance during extracorporeal membrane oxygenation. Antimicrobial Agents and Chemotherapy. 2026;70(6):e00098-26. doi:10.1128/aac.00098-26"
  vignette <- "Yusuff_2026_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # `central` is the single disposition compartment of the one-compartment
  # model (Yusuff 2026 supplement Monolix code: pkmodel(V, Cl)); all
  # observations are arterial-line plasma concentrations (Yusuff 2026 Methods
  # 'Voriconazole administration and sampling strategy').
  compartmentData <- list(
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear scaling of the central volume, V_i = V * (WT / 70), taken from the authors' Monolix code in the supplement ('V_WT = V * (WT/70)', WT declared with use=regressor). The main text does not mention the weight scaling; Table 2 reports V = 144.79 L, which is therefore the volume of a 70 kg patient. Cohort mean 86.8 kg (range 49-134) per Yusuff 2026 Table 1.",
      source_name = "WT"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with CYP2C19_PM = 0, i.e. the pooled normal (*1/*1), rapid (*1/*17) and ultrarapid (*17/*17) metabolizer 'wild-type' group of Yusuff 2026 Figure 2 and Figure S5.",
      notes = "Yusuff 2026 estimates ONE shared effect for intermediate and poor metabolizers together ('patients carrying reduced-function alleles (intermediate/poor metabolizers) had 36% lower CL2 compared with normal metabolizers (beta = -0.44, RSE 60%)', Results and Table 2). It is encoded with the two register indicators summed, exp(e_cyp2c19impm_cl_late * (CYP2C19_IM + CYP2C19_PM)); the indicators are mutually exclusive, so the sum is 0 or 1. Genotype-to-phenotype map from Yusuff 2026 Results: IM = *1/*2 and *2/*17 (n = 12). The effect acts on the LATE intrinsic clearance only.",
      source_name = "CYP2C19 genotype (IM/PM vs WT)"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with CYP2C19_IM = 0 (pooled normal / rapid / ultrarapid metabolizers).",
      notes = "PM = *2/*2 (n = 1 in Yusuff 2026). Shares the single IM/PM coefficient with CYP2C19_IM; see that entry.",
      source_name = "CYP2C19 genotype (IM/PM vs WT)"
    )
  )

  # Covariates Yusuff 2026 screened but did not retain (Results 'Final model
  # parameters and covariate effects'). Documentation only.
  covariatesDataExcluded <- list(
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened; showed a suggestive association with reduced clearance but did not meet retention criteria. Cohort mean 137.1 mg/L (median 105, range 5-409) per Yusuff 2026 Table S1."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort mean 27.1 g/L (range 16-41) per Yusuff 2026 Table S1."
    ),
    BILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as a hepatic-function marker; not retained. Cohort mean 27.3 umol/L (range 3-128) per Yusuff 2026 Table S1."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened as a hepatic-function marker; not retained. Cohort mean 83.2 U/L (range 8-526) per Yusuff 2026 Table S1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31,
    n_studies = 1,
    n_observations = 131,
    age_range = "19-58 years (median 40, mean 39.9)",
    weight_range = "49-134 kg (mean 86.8)",
    sex_female_pct = 54.8,
    race_ethnicity = c(WhiteBritish = 51.6, Asian = 25.8, MixedBritish = 16.1, Other = 6.5),
    disease_state = "Critically ill adults on ECMO for acute severe respiratory and/or cardiac failure (influenza 64.5%, COVID-19 19.4%, other 16.1%) receiving intravenous voriconazole for suspected or confirmed invasive aspergillosis.",
    ecmo_mode = "Venovenous 93.5%, venoarterial 6.5%; CentriMag pump with Paragon oxygenators",
    cyp2c19_phenotype = c(NM = 12, IM = 12, RM = 5, PM = 1, UM = 1),
    co_medication = "All patients received a proton pump inhibitor (omeprazole or lansoprazole); 87.2% received corticosteroids; 61.3% received CVVHDF for all or part of the study.",
    dose_range = "IV 6 mg/kg q12h for two loading doses then 4 mg/kg q12h, infused over 1.5-2 h; adjusted body weight for obese patients.",
    regions = "Single centre: University Hospitals Leicester NHS Trust, United Kingdom.",
    notes = "Prospective observational study (NCT04868188), August 2021 - May 2025. Up to five protocol samples per patient plus routine TDM samples (median 4 per patient, range 1-7). Estimation in Monolix 2024 (SAEM). Baseline characteristics: Yusuff 2026 Table 1; laboratory values: Table S1; parameter estimates: Table 2."
  )

  ini({
    # Volume of distribution for a 70 kg patient (scaled by WT/70 in model()).
    lvc <- log(144.79); label("Central volume of distribution at 70 kg (L)") # Yusuff 2026 Table 2: V = 144.79 L (RSE 8.04%)

    # Intrinsic clearance: logistic transition in time from the early level
    # (CL1) to the late level (CL2), with midpoint Tswitch.
    lcl <- log(6.22); label("Early intrinsic clearance CL1 (L/h)") # Yusuff 2026 Table 2: CL1 = 6.22 L/h (RSE 24.8%)
    lcl_late <- log(22.26); label("Late intrinsic clearance CL2 in CYP2C19 normal/rapid/ultrarapid metabolizers (L/h)") # Yusuff 2026 Table 2: CL2 = 22.26 L/h (RSE 24.4%)
    lcl_t50 <- log(87.96); label("Midpoint time Tswitch of the logistic CL1-to-CL2 transition (h)") # Yusuff 2026 Table 2: Tswitch = 87.96 h (RSE 29.4%)

    # Apparent ECMO-circuit sequestration clearance, decaying mono-exponentially
    # from its value at t = 0. The paper reports the decay as a half-life
    # (T1/2,seq = 4.17 h), so the rate is log(2) / 4.17 = 0.1662 1/h.
    lcl_exp_component <- log(17.68); label("Initial sequestration clearance CLseq0 at t = 0 (L/h)") # Yusuff 2026 Table 2: CLseq0 = 17.68 L/h (RSE 117.3%)
    lcl_exp_kdes <- log(log(2) / 4.17); label("Decay rate of the sequestration clearance, log(2)/T1/2,seq (1/h)") # Yusuff 2026 Table 2: T1/2,seq = 4.17 h (RSE 83.7%)

    # CYP2C19 intermediate/poor metabolizer effect on the late clearance,
    # natural-log scale: exp(-0.44) = 0.644, i.e. 36% lower CL2.
    e_cyp2c19impm_cl_late <- -0.44; label("Log-scale effect of CYP2C19 IM/PM status on late clearance CL2 (unitless)") # Yusuff 2026 Table 2: genotype effect on CL2 = -0.44 (RSE 60.0%)

    # IIV: log-normal; Table 2 prints CV% = sqrt(exp(omega^2) - 1) * 100,
    # so omega^2 = log(1 + CV^2). No IIV on CLseq0 or T1/2,seq (Table 2 'NC').
    etalvc ~ 0.07341 # Yusuff 2026 Table 2: IIV V = 27.6% CV -> log(1 + 0.276^2)
    etalcl ~ 0.46556 # Yusuff 2026 Table 2: IIV CL1 = 77.0% CV -> log(1 + 0.770^2)
    etalcl_late ~ 0.56058 # Yusuff 2026 Table 2: IIV CL2 = 86.7% CV -> log(1 + 0.867^2)
    etalcl_t50 ~ 0.51328 # Yusuff 2026 Table 2: IIV Tswitch = 81.9% CV -> log(1 + 0.819^2)

    # Combined additive + proportional residual error.
    addSd <- 0.091; label("Additive residual error SD (mg/L)") # Yusuff 2026 Table 2: residual additive error SD = 0.091 mg/L
    propSd <- 0.12; label("Proportional residual error (fraction)") # Yusuff 2026 Table 2: residual proportional error = 0.12
  })

  model({
    # Individual parameters (log-normal IIV).
    vc <- exp(lvc + etalvc) * (WT / 70)
    cl_early <- exp(lcl + etalcl)
    cl_late <- exp(lcl_late + e_cyp2c19impm_cl_late * (CYP2C19_IM + CYP2C19_PM) + etalcl_late)
    cl_t50 <- exp(lcl_t50 + etalcl_t50)
    cl_exp_component <- exp(lcl_exp_component)
    cl_exp_kdes <- exp(lcl_exp_kdes)

    # Yusuff 2026 equation 1 and supplement code: logistic transition with the
    # steepness k fixed at 0.1 1/h ('k is the slope parameter (fixed at 0.1)';
    # written as the literal 0.1 in the authors' Monolix code). t is time since
    # the first voriconazole dose.
    cl_intr <- cl_early + (cl_late - cl_early) / (1 + exp(-0.1 * (t - cl_t50)))

    # Yusuff 2026 equation 2: sequestration clearance decays from CLseq0.
    cl_seq <- cl_exp_component * exp(-cl_exp_kdes * t)

    # Yusuff 2026 equation 3: total clearance is the sum of both pathways.
    cl <- cl_intr + cl_seq
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose mg / volume L = mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}

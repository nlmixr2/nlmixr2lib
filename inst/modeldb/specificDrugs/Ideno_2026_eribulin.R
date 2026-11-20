Ideno_2026_eribulin <- function() {
  description <- "Three-compartment IV PK model for eribulin (typical values and covariate equations fixed from Majid 2014) driving an effect compartment and a Friberg myelosuppression model with three transit compartments (Ktr = Kprol, separate Kout, feedback (BNEU/Neu)^gamma, linear drug effect Slope * Ce) for neutrophil counts in Japanese patients with metastatic breast cancer (Ideno 2026, SCOUT Part 3). PD parameters were re-estimated by Ideno 2026 without covariates on pooled phase 2 and post-marketing-survey data, with exponential IIV on Keo, Kprol, Kout, Gamma and Slope and a proportional residual error on the neutrophil count. Doses are in mg of eribulin free base."
  reference <- paste(
    "Ideno Y, Kasai H, Tanigawara Y.",
    "SCOUT: An Exploratory Approach to Scouting Dose-Relevant Covariates.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(4):e70235.",
    "doi:10.1002/psp4.70235.",
    "PD estimates from Data S3; ODE structure from the Part 3 control stream in Data S1.",
    "PK component: Majid O, Gupta A, Reyderman L, Olivo M, Hussein Z.",
    "Population pharmacometric analyses of eribulin in patients with locally",
    "advanced or metastatic breast cancer previously treated with anthracyclines",
    "and taxanes. J Clin Pharmacol. 2014;54(10):1134-1143. doi:10.1002/jcph.315;",
    "typical-value equations as reproduced in Kawamura T, Kasai H, Fermanelli V,",
    "et al. Cancer Sci. 2018;109(9):2822-2829. doi:10.1111/cas.13708 (section 2.3).",
    "PD structure: Friberg LE, et al. J Clin Oncol. 2002;20(24):4713-4721;",
    "effect-compartment variant from Kasai H, et al. Ther Drug Monit. 2023;45:318-326.",
    sep = " "
  )
  vignette <- "Ideno_2026_scout"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    anc = "cells/uL"
  )

  # What each ODE state holds. The SCOUT control stream (Data S1, Part 3)
  # doses eribulin free base into CENTRAL in ug (F1 = DOSE * 1E3) and divides
  # by V1 (L), so its concentration is in ug/L = ng/mL. Here the dose is in mg
  # and Cc is scaled by 1000 to the same ng/mL. The effect compartment holds a
  # concentration (DADT(4) = KEO*(C1 - CE)), and the five Friberg states hold
  # neutrophil counts (cells/uL here; 10^3/uL in the source).
  compartmentData <- list(
    central = list(analyte = "eribulin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "eribulin", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "eribulin", units = "mg", specimen = "tissue", verified = TRUE),
    effect = list(analyte = "eribulin", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    precursor1 = list(analyte = "neutrophil progenitors", units = "cells/uL", specimen = "tissue", verified = TRUE),
    precursor2 = list(analyte = "maturing neutrophils", units = "cells/uL", specimen = "tissue", verified = TRUE),
    precursor3 = list(analyte = "maturing neutrophils", units = "cells/uL", specimen = "tissue", verified = TRUE),
    precursor4 = list(analyte = "maturing neutrophils", units = "cells/uL", specimen = "tissue", verified = TRUE),
    circ = list(analyte = "neutrophils", units = "cells/uL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Majid 2014 PK covariate (Kawamura 2018 section 2.3): CL, Q2, Q3 scale as (WT/68.7)^0.75 and V1, V2, V3 as (WT/68.7). Cohort weight approximately 33-100 kg (Ideno 2026 Figure 4).",
      source_name = "WT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Majid 2014 PK covariate on CL, (ALB/4.0 g/dL)^0.946 (Kawamura 2018 section 2.3). The source column is in g/dL; model() converts g/L to g/dL (x 0.1). Cohort approximately 1.3-5.0 g/dL (Ideno 2026 Figure 4). Not a covariate on the PD layer in Ideno 2026.",
      source_name = "ALB"
    ),
    ALP = list(
      description = "Serum alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Majid 2014 PK covariate on CL, (ALP/132)^-0.209 (Kawamura 2018 section 2.3). Cohort approximately 110-3000 U/L (Ideno 2026 Figure 4, log axis).",
      source_name = "ALP"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Majid 2014 PK covariate on CL, (BILI/0.5 mg/dL)^-0.180 (Kawamura 2018 section 2.3). The source column BILI is in mg/dL; model() converts umol/L to mg/dL (/ 17.1). Cohort approximately 0.1-2.0 mg/dL (Ideno 2026 Figure 4).",
      source_name = "BILI"
    ),
    NEUT = list(
      description = "Baseline absolute neutrophil count (per subject, before the first dose)",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "BNEU in the Data S1 Part 3 control stream (10^3/uL there; supply cells/uL here, i.e. BNEU * 1000). Sets the Friberg initial conditions (progenitor and transit states at BNEU*Kout/Ktr, circulating state at BNEU) and the feedback ratio (BNEU/Neu)^gamma. Cohort approximately 1-15 x 10^3/uL (Ideno 2026 Figure 4).",
      source_name = "BNEU"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 2L,
    weight_range = "approximately 33-100 kg (Ideno 2026 Figure 4)",
    sex_female_pct = 100,
    race_ethnicity = c(Asian = 100),
    disease_state = "Japanese patients with (heavily pretreated) metastatic or recurrent breast cancer receiving eribulin mesilate",
    dose_range = "Approved regimen 1.4 mg/m^2 eribulin mesilate (1.23 mg/m^2 free base) on days 1 and 8 of a 21-day cycle; SCOUT-estimated maximum tolerated doses averaged 1.26 (SD 0.59) mg/m^2 free base for that regimen and 1.97 (SD 1.00) mg/m^2 for the every-other-week regimen",
    regions = "Japan",
    notes = "Pooled anonymised data from a Japanese phase 2 study (Aogi 2012) and the Japanese post-marketing survey analysed by Kawamura 2018 (Ideno 2026 Methods, Part 3). The number of subjects is not stated in the text; patient IDs on the Figure 4 x-axis run to roughly 490."
  )

  ini({
    # ---- PK: Majid 2014 typical values, fixed (not estimated in Ideno 2026) ----
    # Ideno 2026 Methods Part 3: 'the three-compartment model developed by
    # Majid et al. [20] was used for the PK component'; Data S3 reports only
    # PD parameters. Values below are Majid 2014 as printed in Kawamura 2018
    # section 2.3 (same research group, same post-marketing data).
    lcl  <- fixed(log(3.11)); label("Clearance CL at reference covariates (L/h)")    # Kawamura 2018 section 2.3: CL = 3.11*(WGT/68.7)^0.75*(ALB/4.0)^0.946*(ALP/132)^-0.209*(BILI/0.5)^-0.180
    lvc  <- fixed(log(4.06)); label("Central volume V1 at 68.7 kg (L)")               # Kawamura 2018 section 2.3: V1 = 4.06*(WGT/68.7)
    lq   <- fixed(log(2.64)); label("Intercompartmental clearance Q2 at 68.7 kg (L/h)") # Kawamura 2018 section 2.3: Q2 = 2.64*(WGT/68.7)^0.75
    lvp  <- fixed(log(2.42)); label("Peripheral volume V2 at 68.7 kg (L)")            # Kawamura 2018 section 2.3: V2 = 2.42*(WGT/68.7)
    lq2  <- fixed(log(6.60)); label("Intercompartmental clearance Q3 at 68.7 kg (L/h)") # Kawamura 2018 section 2.3: Q3 = 6.60*(WGT/68.7)^0.75
    lvp2 <- fixed(log(121)); label("Peripheral volume V3 at 68.7 kg (L)")             # Kawamura 2018 section 2.3: V3 = 121*(WGT/68.7)

    e_wt_cl    <- fixed(0.75);   label("Allometric exponent of WT on CL, Q2, Q3 (unitless)") # Kawamura 2018 section 2.3: (WGT/68.7)^0.75
    e_alb_cl   <- fixed(0.946);  label("Power exponent of ALB on CL (unitless)")            # Kawamura 2018 section 2.3: (ALB/4.0)^0.946
    e_alp_cl   <- fixed(-0.209); label("Power exponent of ALP on CL (unitless)")            # Kawamura 2018 section 2.3: (ALP/132)^(-0.209)
    e_tbili_cl <- fixed(-0.180); label("Power exponent of BILI on CL (unitless)")           # Kawamura 2018 section 2.3: (BILI/0.5)^(-0.180)

    # ---- PD: Ideno 2026 Data S3 final estimates ----
    lke0   <- log(0.0203); label("Effect-compartment equilibration rate constant Keo (1/h)")  # Data S3 table: Keo (/h) = 0.0203 (RSE 0.110)
    lkprol <- log(0.0570); label("Neutrophil proliferation rate constant Kprol = Ktr (1/h)")  # Data S3 table: Kprol (/h) = 0.0570 (RSE 0.0596)
    lkout  <- log(0.0257); label("Circulating neutrophil elimination rate constant Kout (1/h)") # Data S3 table: Kout (/h) = 0.0257 (RSE 0.0883)
    lgamma <- log(0.237);  label("Feedback exponent Gamma (unitless)")                       # Data S3 table: Gamma = 0.237 (RSE 0.0364)
    lslope <- log(0.0443); label("Linear drug-effect coefficient Slope (mL/ng)")            # Data S3 table: Slope (L/ug) = 0.0443 (RSE 0.0485)

    # Data S3 reports Omega values as SDs of an exponential IIV model;
    # nlmixr2 takes the variance, so each value below is SD^2.
    etalke0   ~ 0.732736  # Data S3 table: Omega-Keo SD = 0.856 (shrinkage 0.422); variance 0.856^2
    etalkprol ~ 0.139876  # Data S3 table: Omega-Kprol SD = 0.374 (shrinkage 0.278); variance 0.374^2
    etalkout  ~ 0.231361  # Data S3 table: Omega-Kout SD = 0.481 (shrinkage 0.431); variance 0.481^2
    etalgamma ~ 0.145924  # Data S3 table: Omega-Gamma SD = 0.382 (shrinkage 0.223); variance 0.382^2
    etalslope ~ 0.180625  # Data S3 table: Omega-Slope SD = 0.425 (shrinkage 0.253); variance 0.425^2

    propSd <- 0.360; label("Proportional residual error on neutrophil count (fraction)") # Data S3 table: RUV prop = 0.360 (RSE 0.0214)
  })

  model({
    # Unit conversions: the Majid 2014 equations use ALB in g/dL and BILI in mg/dL.
    alb_gdl   <- ALB * 0.1     # g/L -> g/dL
    tbili_mgdl <- TBILI / 17.1 # umol/L -> mg/dL
    bneu <- NEUT               # cells/uL

    # ---- PK (Majid 2014; typical values, no IIV) ----
    cl  <- exp(lcl) * (WT / 68.7)^e_wt_cl * (alb_gdl / 4.0)^e_alb_cl * (ALP / 132)^e_alp_cl * (tbili_mgdl / 0.5)^e_tbili_cl
    vc  <- exp(lvc) * (WT / 68.7)
    q   <- exp(lq) * (WT / 68.7)^e_wt_cl
    vp  <- exp(lvp) * (WT / 68.7)
    q2  <- exp(lq2) * (WT / 68.7)^e_wt_cl
    vp2 <- exp(lvp2) * (WT / 68.7)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- PD individual parameters ----
    ke0   <- exp(lke0 + etalke0)
    kprol <- exp(lkprol + etalkprol)
    kout  <- exp(lkout + etalkout)
    gamma <- exp(lgamma + etalgamma)
    slope <- exp(lslope + etalslope)
    ktr   <- kprol # Data S1 Part 3: KTR = KPROL

    # Data S1 Part 3: C1 = A(1)/V1 with A(1) in ug, so C1 is in ng/mL.
    Cc <- 1000 * central / vc

    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(effect)      <- ke0 * (Cc - effect)

    # Data S1 Part 3: EDRUG = SLOPE*CE; NEU floored at 0.01 x 10^3/uL (10 cells/uL)
    # in both the feedback ratio and the circulating-cell elimination term.
    edrug <- slope * effect
    neu <- circ
    if (neu < 10) neu <- 10

    d/dt(precursor1) <- kprol * precursor1 * (1 - edrug) * (bneu / neu)^gamma - ktr * precursor1
    d/dt(precursor2) <- ktr * (precursor1 - precursor2)
    d/dt(precursor3) <- ktr * (precursor2 - precursor3)
    d/dt(precursor4) <- ktr * (precursor3 - precursor4)
    d/dt(circ)       <- ktr * precursor4 - kout * neu

    # Data S1 Part 3: A_0(5..8) = BNEU*KOUT/KTR; A_0(9) = BNEU (system at steady state)
    precursor1(0) <- bneu * kout / ktr
    precursor2(0) <- bneu * kout / ktr
    precursor3(0) <- bneu * kout / ktr
    precursor4(0) <- bneu * kout / ktr
    circ(0)       <- bneu

    ANC <- circ
    ANC ~ prop(propSd)
  })
}

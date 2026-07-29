Kim_2018_tacrolimus <- function() {
  description <- "Integrated population PK model of the tacrolimus (TAC) - mycophenolate mofetil (MMF) drug-drug interaction in healthy Korean male volunteers (Kim 2018, final integrated model). TAC follows a two-compartment model with first-order absorption and a lag time; apparent oral clearance (CL/F) is increased 1.48-fold in CYP3A5 expressers and is suppressed by co-administered mycophenolic acid (MPA) through an inverse-exponential interaction (CL/F = 13.8 / exp(0.0294*Cmpa) * 1.48^CYP3A5). MPA (the active moiety of MMF) follows a two-compartment model with first-order absorption; MPA is metabolised to MPAG (7-O-glucuronide; 85% of metabolism) and AcMPAG (acyl glucuronide; 15%). MPAG undergoes enterohepatic recirculation via a gallbladder compartment that empties into the MPA absorption compartment during a meal-triggered window. Tacrolimus concentrations are in ng/mL; MPA, MPAG and AcMPAG are in ug/mL."
  reference <- paste(
    "Kim JH, Han N, Kim MG, Yun H-Y, Lee S, Bae E, Kim YS, Kim I-W, Oh JM.",
    "Increased Exposure of Tacrolimus by Co-administered Mycophenolate",
    "Mofetil: Population Pharmacokinetic Analysis in Healthy Volunteers.",
    "Sci Rep. 2018;8:1687. doi:10.1038/s41598-018-20071-3.",
    sep = " "
  )
  vignette <- "Kim_2018_tacrolimus"
  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "ng/mL (tacrolimus); ug/mL (MPA, MPAG, AcMPAG)"
  )

  covariateData <- list(
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser status: 1 = CYP3A5 expresser (CYP3A5*1/*1 or CYP3A5*1/*3 at rs776746), 0 = nonexpresser (CYP3A5*3/*3).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 nonexpresser, *3/*3; 13 of 17 subjects in the Kim 2018 cohort)",
      notes              = "Time-fixed germline genotype. In the Kim 2018 healthy-volunteer cohort 4 of 17 (23.5%) subjects were expressers (Table 1). Used as a multiplicative power-form effect on tacrolimus apparent oral clearance: cl_tac = cl_typ * e_cyp3a5_cl^CYP3A5_EXPR with e_cyp3a5_cl = 1.48 in the integrated model (equation (2); Table 3). The CYP3A5 genotype was the only clinical / genetic covariate retained in the final model.",
      source_name        = "CYP3A5"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 17L,
    n_studies      = 1L,
    n_observations = "1,082 concentrations across 4 analytes (TAC, MPA, MPAG, AcMPAG); average 63 concentrations per subject.",
    age_range      = "20-42 years",
    age_median     = "25 years (Table 1)",
    weight_range   = "57.4-88.3 kg",
    weight_median  = "69.7 kg (Table 1)",
    height_range   = "167.7-192.8 cm (median 173.4)",
    sex_female_pct = 0,
    race_ethnicity = "100% Korean (all subjects were of Korean ethnicity).",
    disease_state  = "Healthy male volunteers (no known history of previous disease; normal physical exam and laboratory tests).",
    dose_range     = "Single oral dose: 5 mg tacrolimus (Prograf) and 1,000 mg mycophenolate mofetil (CellCept), administered alone and in combination across a three-period fixed-sequence design.",
    regions        = "Republic of Korea (Seoul National University Hospital).",
    cyp3a5_expressers = "4 of 17 (23.5%) were CYP3A5 expressers (Table 1).",
    notes          = "Three-period, fixed-sequence, open-label, single-dose interaction study (NCT02743247, conducted 2015-2016). Period 1: 1000 mg MMF alone; period 2: 5 mg TAC alone; period 3: combination, with one-week washout between periods. Estimated GFR 74.1-122.3 mL/min/1.73m2 (median 104.8); serum creatinine 0.79-1.20 mg/dL; albumin 4.2-4.9 g/dL. Baseline demographics from Table 1. Variation in baseline demographics was small, so covariates other than CYP3A5 genotype were not significant."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final INTEGRATED-model parameter estimates from Kim 2018 Table 3
    # ("Integrated model" columns). The independent (TAC-alone / MMF-alone)
    # models in Table 3 are intermediate development steps; this file packages
    # the final integrated model that quantifies the TAC-MMF interaction.
    # IIV reported as %CV is converted to the log-normal variance scale via
    # omega^2 = log(1 + CV^2). Compartment / rate-constant numbering follows
    # the paper's Figure 1 (1 = TAC GI, 2 = TAC central, 3 = TAC peripheral,
    # 4 = MPA GI, 5 = MPA central, 6 = MPA peripheral, 7 = MPAG, 8 =
    # gallbladder, 9 = AcMPAG).
    # ----------------------------------------------------------------------

    # --- Tacrolimus (TAC) disposition: 2-compartment, first-order absorption
    #     with lag time. CL/F and V/F are apparent (function of bioavailability).
    lka      <- log(1.78);    label("Tacrolimus absorption rate constant Ka (1/h)")                       # Table 3 Ka TAC (integrated) = 1.78 1/h
    lcl      <- log(13.8);    label("Tacrolimus apparent oral clearance CL/F (L/h), typical value")       # Table 3 CL/F TAC (integrated) = 13.8 L/h
    lvc      <- log(93);      label("Tacrolimus apparent central volume V2/F (L)")                         # Table 3 V2/F (integrated) = 93 L
    lk23     <- log(0.313);   label("Tacrolimus central->peripheral rate constant k23 (1/h)")             # Table 3 k23 (integrated) = 0.313 1/h
    lk32     <- log(0.0719);  label("Tacrolimus peripheral->central rate constant k32 (1/h)")             # Table 3 k32 (integrated) = 0.0719 1/h
    ltlag    <- log(0.59);    label("Tacrolimus absorption lag time (h)")                                 # Table 3 Lag time (integrated) = 0.59 h

    # CYP3A5 expresser effect on TAC CL/F (multiplicative power-of-indicator).
    e_cyp3a5_cl <- 1.48;      label("CYP3A5 expresser multiplier on tacrolimus CL/F (unitless)")          # Table 3 CYP3A5 on CL/F (integrated) = 1.48; equation (2)
    # TAC-MMF interaction: MPA concentration suppresses TAC CL/F via an
    # inverse-exponential relationship cl_tac = cl_typ / exp(slope * Cmpa).
    # Cmpa is in ug/mL (slope therefore has units mL/ug); verified against the
    # Discussion: a 5 ug/mL rise in MPA lowers CL/F by 13.7% (exp(-0.0294*5) =
    # 0.863). This is a structural drug-interaction parameter driven by the
    # model-internal MPA concentration, not a data covariate column.
    e_cmpa_cl <- 0.0294;      label("Inverse-exponential interaction slope of MPA conc on TAC CL/F (mL/ug)") # Table 3 Interaction (integrated) = 0.0294; equation (2)

    # --- Mycophenolic acid (MPA) disposition: 2-compartment, first-order
    #     absorption. The dose is the administered MMF amount (mg); CL/F and
    #     V/F absorb both bioavailability and the MMF->MPA conversion.
    lka_mpa  <- log(2.29);    label("MPA absorption rate constant Ka (1/h)")                              # Table 3 Ka MPA (integrated) = 2.29 1/h
    lcl_mpa  <- log(16.3);    label("MPA apparent oral clearance CL/F (L/h)")                             # Table 3 CL/F MPA (integrated) = 16.3 L/h
    lvc_mpa  <- log(19.7);    label("MPA apparent central volume V5/F (L)")                               # Table 3 V5/F (integrated) = 19.7 L
    lk56     <- log(1.12);    label("MPA central->peripheral rate constant k56 (1/h)")                    # Table 3 k56 (integrated) = 1.12 1/h
    lk65     <- log(0.131);   label("MPA peripheral->central rate constant k65 (1/h)")                    # Table 3 k65 (integrated) = 0.131 1/h
    # Fraction of MPA metabolism that forms MPAG (vs AcMPAG); fixed at 0.85.
    f_mpag   <- fixed(0.85);  label("Fraction of MPA metabolised to MPAG (fMPA), fixed")                  # Table 3 fMPA = 0.85 (fixed)

    # --- MPAG (MPA 7-O-glucuronide): 1-compartment with enterohepatic
    #     recirculation through a gallbladder compartment.
    lvc_mpag <- log(5.83);    label("MPAG apparent volume V7/F (L)")                                      # Table 3 V7/F (integrated) = 5.83 L
    lk70     <- log(0.251);   label("MPAG elimination rate constant k70 (1/h)")                           # Table 3 k70 (integrated) = 0.251 1/h
    lehc     <- log(0.367);   label("Fraction of MPAG undergoing enterohepatic recirculation (EHC)")      # Table 3 EHC (integrated) = 0.367; equation (1)
    lk84     <- log(18.4);    label("Gallbladder emptying rate constant k84 (1/h)")                       # Table 3 k84 (integrated) = 18.4 1/h
    lmtime1  <- log(7.96);    label("Meal time / gallbladder emptying start MTIME1 (h post-dose)")        # Table 3 MTIME1 (integrated) = 7.96 h
    lmtime2  <- fixed(log(1.00)); label("Gallbladder emptying duration MTIME2 (h), fixed")                # Table 3 MTIME2 = 1 h (fixed)

    # --- AcMPAG (MPA acyl glucuronide): 1-compartment; structural parameters
    #     fixed after initial estimation (LAPL+I) because AcMPAG had a
    #     negligible effect on MPA and on the TAC-MMF interaction.
    lvc_acmpag <- fixed(log(23));   label("AcMPAG apparent volume V9/F (L), fixed")                       # Table 3 V9/F = 23 L (fixed)
    lk90       <- fixed(log(2.15)); label("AcMPAG elimination rate constant k90 (1/h), fixed")            # Table 3 k90 = 2.15 1/h (fixed)

    # --- Inter-individual variability (exponential / log-normal). Reported as
    #     %CV in Table 3 (integrated model); omega^2 = log(1 + CV^2).
    etalka     ~ 0.6233       # Table 3 IIV Ka TAC   = 93%   CV -> log(1 + 0.93^2)
    etalcl     ~ 0.0674       # Table 3 IIV CL/F TAC = 26.4% CV -> log(1 + 0.264^2)
    etalvc     ~ 0.0906       # Table 3 IIV V2/F TAC = 30.8% CV -> log(1 + 0.308^2)
    etalka_mpa ~ 0.2779       # Table 3 IIV Ka MPA   = 56.6% CV -> log(1 + 0.566^2)
    etalcl_mpa ~ 0.0344       # Table 3 IIV CL/F MPA = 18.7% CV -> log(1 + 0.187^2)
    etalvc_mpa ~ 0.0326       # Table 3 IIV V5/F MPA = 18.2% CV -> log(1 + 0.182^2)
    etalehc    ~ 0.1187       # Table 3 IIV EHC      = 35.5% CV -> log(1 + 0.355^2)

    # --- Residual error. TAC, MPA and AcMPAG: proportional only; MPAG:
    #     combined additive (ug/mL) + proportional.
    propSd        <- 0.131;   label("Proportional residual error, tacrolimus (fraction)")                # Table 3 sigma_propTAC (integrated) = 0.131
    propSd_mpa    <- 0.524;   label("Proportional residual error, MPA (fraction)")                       # Table 3 sigma_propMPA (integrated) = 0.524
    addSd_mpag    <- 0.104;   label("Additive residual error, MPAG (ug/mL)")                             # Table 3 sigma_addMPAG (integrated) = 0.104
    propSd_mpag   <- 0.237;   label("Proportional residual error, MPAG (fraction)")                      # Table 3 sigma_propMPAG (integrated) = 0.237
    propSd_acmpag <- 0.651;   label("Proportional residual error, AcMPAG (fraction)")                    # Table 3 sigma_propAcMPAG (integrated) = 0.651
  })

  model({
    # ----------------------------------------------------------------------
    # Individual parameters (log-normal IIV).
    # ----------------------------------------------------------------------
    ka      <- exp(lka      + etalka)
    vc      <- exp(lvc      + etalvc)
    k23     <- exp(lk23)
    k32     <- exp(lk32)
    tlag    <- exp(ltlag)

    ka_mpa  <- exp(lka_mpa  + etalka_mpa)
    cl_mpa  <- exp(lcl_mpa  + etalcl_mpa)
    vc_mpa  <- exp(lvc_mpa  + etalvc_mpa)
    k56     <- exp(lk56)
    k65     <- exp(lk65)

    vc_mpag <- exp(lvc_mpag)
    k70     <- exp(lk70)
    # EHC is a fraction in (0, 1). The source modelled its IIV as exponential
    # (log-normal), which can in principle exceed 1 in the upper stochastic
    # tail and would make the derived k78 = ehc*k70/(1-ehc) blow up. The
    # typical value (0.367) and the bulk of the log-normal distribution are
    # unaffected; the min() guard only clips the rare (~0.5%) tail above 0.99
    # to keep stochastic solves numerically stable (see vignette deviations).
    ehc     <- min(exp(lehc + etalehc), 0.99)
    k84     <- exp(lk84)
    mtime1  <- exp(lmtime1)
    mtime2  <- exp(lmtime2)

    vc_acmpag <- exp(lvc_acmpag)
    k90       <- exp(lk90)

    # ----------------------------------------------------------------------
    # MPA central concentration (ug/mL) drives the tacrolimus interaction.
    # MPA dose is in mg and vc_mpa in L, so central_mpa / vc_mpa is in mg/L =
    # ug/mL, the units in which the interaction slope was estimated.
    # ----------------------------------------------------------------------
    Cmpa    <- central_mpa / vc_mpa

    # Tacrolimus CL/F with the CYP3A5 effect and the inverse-exponential MPA
    # interaction (Kim 2018 equation (2)). kel from CL/F and V2/F.
    cl      <- exp(lcl + etalcl) / exp(e_cmpa_cl * Cmpa) * e_cyp3a5_cl^CYP3A5_EXPR
    kel     <- cl / vc

    # ----------------------------------------------------------------------
    # Derived MPA metabolite-formation and EHC rate constants. The paper
    # tabulates CL/F MPA, fMPA, EHC and k70 rather than k57 / k59 / k78
    # directly; these are recovered algebraically from the Figure 1 structure
    # (MPA central is eliminated solely by metabolism to MPAG (k57) and
    # AcMPAG (k59), so k57 + k59 = CL/F MPA / V5/F) and equation (1)
    # (EHC = k78 / (k70 + k78)):
    #   k57 = fMPA       * (CL/F MPA / V5/F)   -> MPA -> MPAG formation
    #   k59 = (1 - fMPA) * (CL/F MPA / V5/F)   -> MPA -> AcMPAG formation
    #   k78 = EHC * k70 / (1 - EHC)            -> MPAG -> gallbladder transport
    # ----------------------------------------------------------------------
    ke_mpa  <- cl_mpa / vc_mpa
    k57     <- f_mpag * ke_mpa
    k59     <- (1 - f_mpag) * ke_mpa
    k78     <- ehc * k70 / (1 - ehc)

    # Meal-triggered gallbladder emptying window (single-dose study): the
    # gallbladder empties into the MPA absorption compartment at rate k84
    # during [MTIME1, MTIME1 + MTIME2] post-dose, completing the EHC loop.
    empty_on <- (t >= mtime1) * (t <= (mtime1 + mtime2))
    k84_eff  <- k84 * empty_on

    # ----------------------------------------------------------------------
    # ODE system. Compartment declaration order matches the paper's Figure 1
    # numbering (1..9).
    # ----------------------------------------------------------------------
    # Tacrolimus (compartments 1-3)
    d/dt(depot)           <- -ka * depot
    d/dt(central)         <-  ka * depot - kel * central - k23 * central + k32 * peripheral1
    d/dt(peripheral1)     <-  k23 * central - k32 * peripheral1
    lag(depot)            <- tlag

    # Mycophenolic acid (compartments 4-6). The MPA absorption compartment
    # also receives the enterohepatic return from the gallbladder.
    d/dt(depot_mpa)       <- -ka_mpa * depot_mpa + k84_eff * gallbladder
    d/dt(central_mpa)     <-  ka_mpa * depot_mpa - k56 * central_mpa + k65 * peripheral1_mpa - k57 * central_mpa - k59 * central_mpa
    d/dt(peripheral1_mpa) <-  k56 * central_mpa - k65 * peripheral1_mpa

    # MPAG (compartment 7) and gallbladder (compartment 8): MPAG is eliminated
    # (k70) and transported to the gallbladder (k78); the gallbladder empties
    # back into the MPA absorption compartment during the meal window.
    d/dt(central_mpag)    <-  k57 * central_mpa - k70 * central_mpag - k78 * central_mpag
    d/dt(gallbladder)     <-  k78 * central_mpag - k84_eff * gallbladder

    # AcMPAG (compartment 9)
    d/dt(central_acmpag)  <-  k59 * central_mpa - k90 * central_acmpag

    # ----------------------------------------------------------------------
    # Observations. Tacrolimus is reported in ng/mL: central is in mg and vc
    # in L, so central / vc is in mg/L = 1000 ng/mL (hence the 1000 factor).
    # MPA, MPAG and AcMPAG are reported in ug/mL (mg/L), no scaling needed.
    # ----------------------------------------------------------------------
    Cc        <- 1000 * central / vc
    Cc_mpa    <- central_mpa / vc_mpa
    Cc_mpag   <- central_mpag / vc_mpag
    Cc_acmpag <- central_acmpag / vc_acmpag

    Cc        ~ prop(propSd)
    Cc_mpa    ~ prop(propSd_mpa)
    Cc_mpag   ~ add(addSd_mpag) + prop(propSd_mpag)
    Cc_acmpag ~ prop(propSd_acmpag)
  })
}

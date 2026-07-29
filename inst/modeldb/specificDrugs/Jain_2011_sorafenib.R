Jain_2011_sorafenib <- function() {
  description <- "One-compartment population PK model for orally administered sorafenib in patients with solid tumours (Jain 2011). Absorption is described by an Erlang-style chain of four catenary GI transit compartments downstream of an upstream absorption depot, all linked by a single first-order rate constant ka (mean absorption transit time MAT = 5 / ka). Enterohepatic recirculation is modelled by routing a fraction Fent of the drug leaving the central compartment into a gallbladder reservoir, with periodic release back to the most distal transit compartment gated by a smooth Hill switch Ehc = tad^40 / (tad^40 + t'^40), where tad is the time since the most recent dose; release becomes essentially full once tad exceeds the gallbladder-emptying onset time t'. The irreversible elimination rate constant ke equals the biliary excretion rate constant kb (= CL/V) per the published assumption kb = ke. Body weight is the only retained covariate (allometric exponent fixed to 1 on V/F, reference weight 80 kg)."
  reference <- paste(
    "Jain L, Woo S, Gardner ER, Dahut WL, Kohn EC, Kummar S, Mould DR,",
    "Giaccone G, Yarchoan R, Venitz J, Figg WD.",
    "Population pharmacokinetic analysis of sorafenib in patients with",
    "solid tumours.",
    "Br J Clin Pharmacol. 2011;72(2):294-305.",
    "doi:10.1111/j.1365-2125.2011.03963.x.",
    sep = " "
  )
  vignette <- "Jain_2011_sorafenib"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on V/F with the exponent FIXED to 1 (linear) and reference weight 80 kg (the cohort median). Body weight was the only covariate retained in the final PPK model (Table 3, p. 300). The paper notes that inclusion of WT on V/F explained only ~4% of IIV in V/F but reached the OFV-significance threshold during the stepwise forward-addition step.",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area (Dubois & Dubois formula)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened during covariate model building but not retained in the final model (Table 2 / Results, p. 300)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained (p. 300)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Mean body weight was higher for males (87 kg vs 71 kg, p < 0.0001). Tested as a covariate but not retained in the final model after stepwise backward elimination (p. 300)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Identified as significant on V/F during initial covariate evaluation but not retained in the final model because of statistically insignificant increases in OFV during stepwise backward elimination (p. 300)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Identified as significant on V/F during initial covariate evaluation but not retained in the final model (p. 300)."
    ),
    CLCR = list(
      description = "Creatinine clearance (Cockcroft-Gault formula)",
      units       = "mL/s",
      type        = "continuous",
      notes       = "Screened but not retained (p. 300)."
    ),
    CYP3A4_1B = list(
      description = "CYP3A4*1B genotype (heterozygous/homozygous variant indicator)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not retained; genetic polymorphisms in CYP3A4*1B, CYP3A5*3C, and UGT1A9*3 did not have any significant effect on PK model parameters (p. 300). UGT1A9*5 was non-polymorphic in the cohort."
    ),
    CYP3A5_3C = list(
      description = "CYP3A5*3C genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not retained (p. 300)."
    ),
    UGT1A9_3 = list(
      description = "UGT1A9*3 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not retained (p. 300). The UGT1A9*3 polymorphism did not follow Hardy-Weinberg equilibrium in the cohort due to small numbers of variant carriers."
    ),
    UGT1A9_5 = list(
      description = "UGT1A9*5 genotype",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "None of the patients carried the variant allele for UGT1A9*5; not evaluable (p. 300)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 111L,
    n_studies      = 5L,
    age_range      = "30.3-84.9 years",
    age_median     = "63.9 years",
    weight_range   = "50.1-132.5 kg",
    weight_median  = "81.4 kg",
    sex_female_pct = 31,
    race_ethnicity = c(Caucasian = 81, AfricanAmerican = 11, Other = 8),
    disease_state  = "Adults with metastatic castrate-resistant prostate cancer (PC, n=46), non-small cell lung cancer (NSCLC, n=18), colorectal cancer (CRC, n=18), other solid tumours (ST, n=28), or Kaposi's sarcoma (KS, n=2). Five phase I/II open-label, single-arm, single-centre trials at the US NCI.",
    dose_range     = "Sorafenib 200 or 400 mg twice daily (oral, continuous cycles) as a single agent or in combination with cetuximab, bevacizumab, or a protease inhibitor.",
    regions        = "USA (National Cancer Institute, Bethesda MD).",
    notes          = "Baseline demographics from Table 2 (n = 111 total). Median (range) bodyweight 81.4 (35.2-132.5) kg; the ST and KS cohorts in Table 2 also list small numbers within-cohort ranges (e.g. ST 92.3-94.0). The cohort was predominantly male Caucasian. All patients were genotyped for CYP3A4*1B, CYP3A5*3C, UGT1A9*3, and UGT1A9*5 SNPs. A total of 1249 sorafenib concentrations were used for PPK model building, of which 3.3% were below the LLOQ and reported as actual values (p. 298). One patient was excluded from analysis because of extremely low plasma concentrations (Results, p. 298)."
  )

  ini({
    # Structural parameters - reference values for an 80 kg patient (Table 3, p. 300).
    # The paper estimates apparent clearance CL/F and apparent central volume V/F
    # because sorafenib is dosed orally and absolute bioavailability is unknown.
    lmtt    <- log(1.98);  label("Mean absorption transit time MAT (h); ka = 5/MAT")            # Table 3 MAT = 1.98 h (CI 1.75-2.25 h)
    lcl     <- log(8.13);  label("Apparent clearance CL/F (L/h)")                                # Table 3 CL/F = 8.13 L/h (CI 6.98-10.7 L/h)
    lvc     <- log(213);   label("Apparent central volume V/F (L) at WT = 80 kg")                # Table 3 V/F = 213 L (CI 196-232 L)
    lkehc   <- log(0.857); label("EHC release rate constant kEhc (1/h)")                         # Table 3 kEhc = 0.857 1/h (CI 0.478-1.26)
    ltprime <- log(6.13);  label("Gallbladder-emptying onset time t' (h)")                       # Table 3 t' = 6.13 h (CI 5.78-7.34)
    fent    <- 0.498;      label("Fraction of dose undergoing enterohepatic recycling (unitless)") # Table 3 Fent = 0.498 (CI 0.464-0.498); derived in-paper as Fent = 1/(1 + Ferts) with the estimated logit-scale parameter Ferts (Table 3 footnote)

    # Allometric exponent on V/F, FIXED to 1 by the authors (Results, p. 300:
    # "Baseline bodyweight was found to be a statistically significant
    # covariate for V/F and was modelled as a power function with the exponent
    # fixed to 1"). Reference weight 80 kg (cohort median).
    e_wt_vc <- fixed(1.0); label("Allometric exponent on V/F (fixed)")                            # Results p. 300 ("exponent fixed to 1")

    # Inter-individual variability. Paper-reported %CV translates to log-scale
    # variance via omega^2 = log(1 + CV^2). IIV on CL/F and V/F are reported
    # together with their correlation coefficient r = 0.778 (Table 3).
    #
    #   CL/F IIV 18.0% -> log(1 + 0.180^2)  = 0.0319
    #   V/F  IIV 68.7% -> log(1 + 0.687^2)  = 0.3866
    #   ka   IIV 61.9% -> log(1 + 0.619^2)  = 0.3240
    #   cov(CL,V) = r * sqrt(omega^2_CL) * sqrt(omega^2_V)
    #             = 0.778 * sqrt(0.0319) * sqrt(0.3866) = 0.0864
    #
    # IIV on ka was reported on the linear ka scale (Table 3 row "IIV ka").
    # Because MAT = 5 / ka, var(log MAT) = var(-log ka + log 5) = var(log ka),
    # so the same variance applies whether the model parameterises MAT or ka.
    etalcl + etalvc ~ c(0.0319,
                        0.0864, 0.3866)                                                          # Table 3 IIV CL/F = 18.0%, IIV V/F = 68.7%, correlation 0.778
    etalmtt ~ 0.3240                                                                              # Table 3 IIV ka = 61.9% (variance equal on the log MAT scale)

    # Residual error - combined additive + proportional, as reported in
    # equation (1) on p. 297. The paper's log-scale form
    #   ln(Cij) = ln(Cij_hat) + sqrt(eps_p^2 + eps_a^2 / Cij_hat^2)
    # is a Taylor-equivalent of linear-scale combined error
    #   C = C_hat + N(0, sqrt(propSd^2 * C_hat^2 + addSd^2)),
    # so the values map directly onto the nlmixr2 add() + prop() pair.
    propSd <- 0.514;    label("Proportional residual error (fraction)")                          # Table 3 proportional residual = 51.4 %CV (CI 45.6-56%)
    addSd  <- 1.0003e-3; label("Additive residual error (mg/L)")                                  # Table 3 additive residual = 1.0003e-3 mg/L (essentially the LLOQ floor; CI 1.0001e-3 - 1.0003e-3)
  })
  model({
    # ---- Individual structural parameters ----
    # MAT (h) with IIV; ka (1/h) is derived as 5 / MAT because the transit
    # chain has 4 transit compartments + 1 absorption depot (Table 3 footnote:
    # "Mean absorption transit time = (number of transit compartments + 1) / ka
    # = 5 / 2.53 = 1.98").
    mtt    <- exp(lmtt + etalmtt)
    ka     <- 5 / mtt

    # CL/F, V/F and the EHC parameters.
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc
    kehc   <- exp(lkehc)
    tprime <- exp(ltprime)

    # Total first-order elimination rate constant from the central
    # compartment. The paper assumes the biliary excretion rate constant kb
    # equals the irreversible elimination rate constant ke (Results, p. 300:
    # "The model assuming kb = ke was stable and adequately described the EHC
    # for sorafenib."), and ke is parameterised as CL/V. The fraction Fent of
    # the central-compartment outflow is therefore routed to the gallbladder
    # via Fent * kel * central, while the irreversible loss (hepatic
    # metabolism plus non-recirculated biliary loss) is (1 - Fent) * kel *
    # central.
    kel <- cl / vc

    # ---- Enterohepatic gating function ----
    # Equation (10): Ehc = (t - DT)^40 / ((t - DT)^40 + t'^40), where
    # DT is the time of the most recent dose. tad() is the rxode2 helper
    # that returns the time since the most recent dose event. With the Hill
    # exponent fixed at 40 this acts as a smooth square-wave switch that
    # rises from ~0 to ~1 once tad exceeds t'. tad() resets at each dose so
    # the switch fires once per dose interval, exactly matching Figure 2.
    tad_h <- tad()
    Ehc <- tad_h^40 / (tad_h^40 + tprime^40)

    # ---- ODE system ----
    # Equations (4)-(9):
    #   dA0/dt   = -ka * A0
    #   dAk/dt   =  ka * (A(k-1) - Ak) for k = 1, 2, 3
    #   dA4/dt   =  ka * (A3 - A4) + Ehc * kEhc * Agb
    #   dAcc/dt  =  ka * A4 - Fent * kb * Acc - (1 - Fent) * (CL/V) * Acc
    #   dAgb/dt  =  Fent * kb * Acc - Ehc * kEhc * Agb
    # with kb = ke = CL/V (= kel here).
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot    - ka * transit1
    d/dt(transit2)    <-  ka * transit1 - ka * transit2
    d/dt(transit3)    <-  ka * transit2 - ka * transit3
    d/dt(transit4)    <-  ka * transit3 - ka * transit4 + Ehc * kehc * gallbladder
    d/dt(central)     <-  ka * transit4 - kel * central
    d/dt(gallbladder) <-  fent * kel * central - Ehc * kehc * gallbladder

    # Plasma concentration. Dose in mg, vc in L, central in mg, so
    # central / vc has units of mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

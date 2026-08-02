Brochot_2015_darunavir <- function() {
  description <- "Two-compartment population PK model for ritonavir-boosted darunavir with alpha-1 acid glycoprotein (AAG)-dependent apparent clearance and allometric weight scaling, in HIV-1-infected pediatric (>=3 to <18 years) and adult patients (Brochot 2015)."
  reference <- "Brochot A, Kakuda TN, Van De Casteele T, Opsomer M, Tomaka FL, Vermeulen A, Vis P. Model-Based Once-Daily Darunavir/Ritonavir Dosing Recommendations in Pediatric HIV-1-Infected Patients Aged >=3 to <12 Years. CPT Pharmacometrics Syst Pharmacol. 2015;4(7):406-414. doi:10.1002/psp4.44"
  vignette <- "Brochot_2015_darunavir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (per-visit / time-varying).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 70 kg (adult) in allometric power scaling on CL/F and V2/F. Time-varying per NONMEM column WTV in supplement control stream.",
      source_name        = "WTV"
    ),
    AAG = list(
      description        = "Serum alpha-1 acid glycoprotein (orosomucoid) concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Canonical AAG unit is g/L (100 mg/dL = 1 g/L). Paper Table 3 reports KAFF = 0.0304 dL/mg with AAG in mg/dL; converted to canonical L/g: e_aag_cl = 0.0304 * 100 = 3.04 L/g. Non-linear inverse-saturation effect on CL/F: CL/F is multiplied by 1 / (1 + e_aag_cl * AAG). Paper's pooled pediatric-cohort 5th/50th/95th percentiles were 0.564 / 0.912 / 1.57 g/L (56.4 / 91.2 / 157 mg/dL, Methods 'Simulations').",
      source_name        = "AAG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 102,
    n_studies      = 5,
    age_range      = "3-66 years",
    weight_range   = "12-96 kg",
    disease_state  = "HIV-1 infection (treatment-experienced and treatment-naive)",
    dose_range     = "300-800 mg darunavir with 50-100 mg ritonavir boosting; once- or twice-daily oral",
    notes          = paste(
      "Pooled population-PK dataset (Brochot 2015 Table 1) from five studies:",
      "DUET-1 + DUET-2 (adults 18-66 y, N=30, DRV/r 600/100 mg BID);",
      "DELPHI (children 6-17 y, N=41, DRV/r 300-600/50-100 mg BID);",
      "ARIEL (children 3-5 y, N=24, DRV/r 20/3 mg/kg BID);",
      "ARIEL once-daily substudy (children 3-5 y, N=10, DRV/r 40/7 mg/kg or 600/100 mg QD);",
      "DIONE (adolescents 12-17 y, N=12, DRV/r 800/100 mg QD).",
      "659 darunavir plasma concentrations total (423 pediatric).",
      "Ritonavir was co-administered as boosting agent but ritonavir concentrations were not modelled;",
      "darunavir concentrations without a time-matching quantifiable ritonavir concentration were excluded."
    )
  )

  ini({
    # Structural parameters -- Brochot 2015 Table 3 (final estimates after DIONE + ARIEL QD update)
    lcl     <- log(51.0);  label("Apparent intrinsic clearance CLint/F (L/h); before AAG and WT effects")   # Brochot 2015 Table 3, CLint/F = 51.0 L/h (SEE 4.7%)
    lvc     <- log(137);   label("Apparent central volume V2/F (L); at reference WT = 70 kg")               # Brochot 2015 Table 3, V2/F = 137 L (SEE 21%)
    lq      <- log(19.1);  label("Apparent intercompartmental clearance Q/F (L/h)")                          # Brochot 2015 Table 3, Q/F = 19.1 L/h (SEE 16%)
    lvp     <- log(254);   label("Apparent peripheral volume V3/F (L)")                                      # Brochot 2015 Table 3, V3/F = 254 L (SEE 41%)
    lka     <- log(0.528); label("First-order absorption rate KA (1/h)")                                     # Brochot 2015 Table 3, KA = 0.528 1/h (SEE 17%)
    lfdepot <- fixed(log(1.18)); label("Log relative bioavailability Frel (log-fraction, unitless); log(1.18)")  # Brochot 2015 Table 3, Frel = 1.18 FIXED (formulation correction from upstream Vis 2006 abstract; NONMEM supplement $THETA(8) FIXED, applied as F1 in $PK block on depot)

    # Covariate effects
    e_wt_cl  <- 0.504;         label("Allometric power exponent of body weight on CL/F (unitless)")          # Brochot 2015 Table 3, "Influence of WT on CL/F" = 0.504 (SEE 11%)
    e_wt_vc  <- 0.774;         label("Allometric power exponent of body weight on V2/F (unitless)")          # Brochot 2015 Table 3, "Influence of WT on V2/F" = 0.774 (SEE 18%)
    e_aag_cl <- fixed(3.04);   label("AAG effect on CL/F, e_aag_cl (L/g) (inverse-saturation form)")  # Brochot 2015 Table 3, KAFF = 0.0304 dL/mg FIXED = 3.04 L/g in canonical g/L units; NONMEM supplement $THETA(3) FIXED

    # IIV -- OMEGA reported as %CV in Table 3; interpreted as sqrt(OMEGA)*100 (matches NONMEM $OMEGA initial estimates 0.0858/0.422/0.5)
    etalcl ~ 0.0784   # Brochot 2015 Table 3, IIV CLint/F = 28% CV -> OMEGA = 0.28^2 = 0.0784 (SEE 20%)
    etalq  ~ 0.4096   # Brochot 2015 Table 3, IIV Q/F     = 64% CV -> OMEGA = 0.64^2 = 0.4096 (SEE 59%)
    etalka ~ 0.25     # Brochot 2015 Table 3, IIV KA      = 50% CV -> OMEGA = 0.50^2 = 0.25   (SEE 66%)

    # Residual error -- proportional; NONMEM supplement $SIGMA value converted to SD scale for nlmixr2
    propSd <- sqrt(0.0717); label("Proportional residual SD")   # Brochot 2015 Table 3, multiplicative residual error variance = 0.0717 (SEE 12%) -> SD = sqrt(0.0717) ~ 0.2678
  })
  model({
    # Individual PK parameters
    #   CL/F_i = CLint/F * 1 / (1 + KAFF * AAG_i) * (WT_i / 70)^e_wt_cl * exp(eta_cl)
    #   V2/F_i = V2/F * (WT_i / 70)^e_wt_vc
    # (Frel is applied as bioavailability of the depot compartment, matching NONMEM F1 = 1.18)
    cl <- exp(lcl + etalcl) / (1 + e_aag_cl * AAG) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq + etalq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Concentration: dose in mg, vc in L -> mg/L; multiply by 1000 -> ng/mL (matches LC-MS/MS reporting, LLOQ 5 ng/mL)
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}

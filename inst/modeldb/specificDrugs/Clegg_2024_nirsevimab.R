Clegg_2024_nirsevimab <- function() {
  description <- "Two-compartment population PK model for nirsevimab in preterm and term infants (Clegg 2024)"
  reference <- "Clegg L, Freshwater E, Leach A, Villafana T, Wahlby Hamren U. Population Pharmacokinetics of Nirsevimab in Preterm and Term Infants. J Clin Pharmacol. 2024;64(5):555-567. doi:10.1002/jcph.2401"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  covariates <- list(
    WT            = "Body weight (kg), time-varying",
    PAGE          = "Postmenstrual age (months; GA in weeks / 4.35 + postnatal age in months), time-varying",
    BLACK_OTH     = "1 if Black/African American or Other race, 0 otherwise (reference: White or Native Hawaiian/Pacific Islander)",
    ASIAN_AMIND_MULTI = "1 if Asian, American Indian/Alaskan Native, or Multiple races, 0 otherwise",
    SEASON2       = "1 if second RSV season, 0 if first RSV season",
    ADA_POS       = "1 if antidrug antibody positive, 0 if negative"
  )
  ini({
    # Structural parameters — reference values for a 70 kg adult
    # Units: CL and Q in L/day; V2 and V3 in L; Ka in 1/day
    # (Converted from published mL/day and mL by dividing by 1000)
    lka      <- log(0.401);   label("Absorption rate constant (Ka, 1/day)")
    lcl      <- log(0.0388);  label("Clearance for a 70 kg adult (CL, L/day)")
    lvc      <- log(1.98);    label("Central volume of distribution for a 70 kg adult (V2, L)")
    lvp      <- log(2.4);     label("Peripheral volume of distribution for a 70 kg adult (V3, L)")
    lq       <- log(0.709);   label("Intercompartmental clearance for a 70 kg adult (Q, L/day)")
    lfdepot  <- log(0.839);   label("Intramuscular bioavailability (F, fraction)")

    # Allometric exponents
    allo_cl <- 0.589; label("Allometric exponent on CL and Q (unitless)")
    allo_v  <- 0.84;  label("Allometric exponent on V2 and V3 (unitless)")

    # Maturation parameters for CL (asymptotic function centered at term birth: PAGE = 40/4.35 months)
    beta_cl <- 0.364; label("Fraction of mature CL at birth for a term infant (GA 40 weeks, unitless)")
    t50_cl  <- 14.8;  label("Maturation half-life for CL (months)")

    # Race effects on CL
    # Reference group: White or Native Hawaiian/Pacific Islander
    e_black_oth_cl       <-  0.132;  label("Race effect on CL: Black/African American or Other vs reference (fraction)")
    e_asian_amind_multi_cl <- -0.0894; label("Race effect on CL: Asian, American Indian/Alaskan Native, or Multiple vs reference (fraction)")

    # Race effect on V2
    # Reference group: White, Black/AA, Native Hawaiian/Pacific Islander, or Other
    e_asian_amind_multi_v2 <- -0.226; label("Race effect on V2: Asian, American Indian/Alaskan Native, or Multiple vs reference (fraction)")

    # Season and ADA effects on CL
    e_season2_cl <- -0.122; label("Season 2 effect on CL (fraction)")
    e_ada_cl     <-  0.124; label("Antidrug antibody (ADA) positive effect on CL (fraction)")

    # IIV: omega^2 = log(CV^2 + 1); CL and V2 are correlated (r = 0.785)
    etacl + etavc ~ c(0.06543,
                      0.08282, 0.17005)
    etaka ~ 0.17718  # 44% CV; note high eta shrinkage (83%) in the original analysis

    # Residual error (proportional; additive on log-scale in NONMEM)
    propSd <- 0.21; label("Proportional residual error (fraction)")
  })
  model({
    # Maturation function for CL, centered at term birth (PAGE = 40/4.35 months)
    page_centered <- PAGE - 40 / 4.35
    maturation_cl <- 1 - (1 - beta_cl) * exp(-page_centered * log(2) / t50_cl)

    # Race effects on CL and V2 (multiplicative; reference groups above)
    race_cl <- 1 + e_black_oth_cl * BLACK_OTH + e_asian_amind_multi_cl * ASIAN_AMIND_MULTI
    race_v2 <- 1 + e_asian_amind_multi_v2 * ASIAN_AMIND_MULTI

    # Season and ADA effects on CL
    season_cl <- 1 + e_season2_cl * SEASON2
    ada_cl    <- 1 + e_ada_cl * ADA_POS

    # PK parameters with allometric weight scaling (reference 70 kg)
    ka <- exp(lka + etaka)
    cl <- exp(lcl + etacl) * (WT / 70)^allo_cl * maturation_cl * race_cl * season_cl * ada_cl
    vc <- exp(lvc + etavc) * (WT / 70)^allo_v  * race_v2
    vp <- exp(lvp)         * (WT / 70)^allo_v
    q  <- exp(lq)          * (WT / 70)^allo_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

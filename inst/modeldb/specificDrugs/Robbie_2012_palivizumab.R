Robbie_2012_palivizumab <- function() {
  description <- "Two-compartment population PK model for palivizumab (anti-RSV humanized IgG1 kappa mAb) with first-order IM absorption in adults and children (Robbie 2012)"
  reference <- "Robbie GJ, Zhao L, Mondick J, Losonsky G, Roskos LK. Population Pharmacokinetics of Palivizumab, a Humanized Anti-Respiratory Syncytial Virus Monoclonal Antibody, in Adults and Children. Antimicrob Agents Chemother. 2012;56(9):4927-4936. doi:10.1128/AAC.06446-11. Erratum in: Antimicrob Agents Chemother. 2012;56(10):5431 (PMC3457364) adding Anderson BJ, Allegaert K, Holford NH. 2006. Eur J Pediatr. 165:819-829 as reference 1a for the maturation model."
  vignette <- "Robbie_2012_palivizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; used for allometric scaling on CL and Q (exponent 0.75) and on Vc and Vp (exponent 1.0) with reference weight 70 kg (Table 2).",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age = postnatal age (months) + gestational age (weeks) / 4.35, per Robbie 2012 equation (1).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the clearance maturation function (Robbie 2012 eq. 1a / eq. 4). Reference is term birth (PAGE = 40/4.35 months = 9.195 months).",
      source_name        = "PAGE"
    ),
    RACE_BLACK = list(
      description        = "Race indicator: 1 if Black / African American, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the paper's reference race).",
      notes              = "Multiplicative effect on CL; 95% CI included unity in the source analysis. Renamed from source column BLACK to canonical RACE_BLACK per covariate-columns.md.",
      source_name        = "BLACK"
    ),
    RACE_HISPANIC = list(
      description        = "Race indicator: 1 if Hispanic / Latino, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the paper's reference race).",
      notes              = "Multiplicative effect on both CL and Vc; 95% CI included unity in the source analysis. Renamed from source column HISPANIC to canonical RACE_HISPANIC per covariate-columns.md.",
      source_name        = "HISPANIC"
    ),
    RACE_ASIAN = list(
      description        = "Race indicator: 1 if Asian, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the paper's reference race).",
      notes              = "Multiplicative effect on CL; 95% CI included unity in the source analysis. Renamed from source column ASIAN to canonical RACE_ASIAN per covariate-columns.md.",
      source_name        = "ASIAN"
    ),
    RACE_OTHER = list(
      description        = "Race indicator: 1 if race = Other, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the paper's reference race).",
      notes              = "Multiplicative effect on CL; 95% CI included unity in the source analysis. Renamed from source column OTHER to canonical RACE_OTHER per covariate-columns.md.",
      source_name        = "OTHER"
    ),
    CLD_PREM = list(
      description        = "Chronic lung disease of prematurity (bronchopulmonary dysplasia).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no CLD of prematurity.",
      notes              = "Multiplicative +20% effect on CL (Robbie 2012 Table 2). Renamed from source column CLD to canonical CLD_PREM per covariate-columns.md.",
      source_name        = "CLD"
    ),
    ADA_TITER = list(
      description        = "Antidrug-antibody reciprocal-dilution titer at the matched PK sample time. American linear-titer convention: 0 = ADA negative.",
      units              = "(reciprocal dilution; 0 = negative)",
      type               = "continuous",
      reference_category = "0 (ADA negative).",
      notes              = "Robbie 2012 modelled ADA titer as a step-function multiplicative effect on CL with four non-reference categories (10, 20, 40, >=80); only the >=80 category CI excluded unity. In the implementation we derive the category indicators inside model() from the continuous ADA_TITER column using comparison expressions, so users supply a single ADA_TITER column rather than four separate indicator columns. Renamed from source column ADA titer to canonical ADA_TITER per covariate-columns.md.",
      source_name        = "ADA titer"
    )
  )

  population <- list(
    n_subjects     = 1883,
    n_studies      = 22,
    n_adults       = 116,
    n_pediatric    = 1767,
    age_range      = "Pediatric PAGE 6.89 to 34.7 months; adult 19.3 to 58.9 years.",
    age_median     = "Pediatric PAGE median 12.6 months; adult mean 32.8 years.",
    weight_range   = "Pediatric 0.92 to 16.3 kg; adult 48.9 to 104.3 kg.",
    weight_median  = "Pediatric median 5.64 kg; adult mean 72.5 kg.",
    ga_range       = "Pediatric gestational age 22 to 41 weeks (median 30 weeks); 62% preterm, 38% term.",
    sex_female_pct = "Pediatric 45.3% female (752 girls / 1660 reported); adult 65.5% female (76/116).",
    race_ethnicity = "Pediatric 54% White, 20% Black / African American, 18% Hispanic / Latino, 2% Asian, 6% Other.",
    disease_state  = "Pediatric cohort at high risk of RSV lower respiratory tract disease (preterm, bronchopulmonary dysplasia / chronic lung disease of prematurity, congenital heart disease); 36% had chronic lung disease of prematurity. Adult cohort: 7 healthy-volunteer studies and 2 hematopoietic / solid-organ transplant studies.",
    dose_range     = "3 to 30 mg/kg; label regimen 15 mg/kg IM monthly (up to 5 doses per RSV season).",
    regions        = "North America and Europe (22 pooled MedImmune studies).",
    notes          = "Baseline demographics taken from Robbie 2012 Table 1. Adult cohort dosed by a mix of IM and IV routes across 9 studies; pediatric cohort dosed exclusively IM. 1,661 adult PK samples and 4,095 pediatric PK samples."
  )

  ini({
    # Structural parameters for a 70 kg adult (white, 40-week PAGE, ADA titer = 0);
    # Robbie 2012 Table 2. Paper reports CL and Q in mL/day and Vc/Vp in mL;
    # converted to L/day and L here by dividing by 1000.
    lka      <- log(1.01);     label("Absorption rate constant (ka, 1/day)")     # Robbie 2012 Table 2
    lcl      <- log(0.198);    label("Clearance for a 70 kg adult (CL, L/day)")  # Robbie 2012 Table 2 (198 mL/day)
    lvc      <- log(4.09);     label("Central volume for a 70 kg adult (Vc, L)") # Robbie 2012 Table 2 (4,090 mL)
    lvp      <- log(2.23);     label("Peripheral volume for a 70 kg adult (Vp, L)")        # Robbie 2012 Table 2 (2,230 mL)
    lq       <- log(0.879);    label("Intercompartmental clearance for a 70 kg adult (Q, L/day)") # Robbie 2012 Table 2 (879 mL/day)
    lfdepot  <- log(0.694);    label("Intramuscular bioavailability (F, fraction)")        # Robbie 2012 Table 2

    # Allometric exponents (Robbie 2012 eq. 3 and Table 2)
    allo_cl <- 0.75; label("Allometric exponent on CL and Q (unitless)")
    allo_v  <- 1.0;  label("Allometric exponent on Vc and Vp (unitless)")

    # Maturation parameters for CL (Robbie 2012 eq. 1a; asymptotic-exponential form centered on
    # 40-week PAGE = 9.195 months). beta is CL at term as fraction of the allometrically scaled
    # adult CL; TCL is the maturation half-life. The equation as implemented below reproduces
    # Robbie 2012 Table 3 pediatric reference (CL = 11.0 mL/day at 4.5 kg, 12.3-mo PAGE) and
    # the abstract value 10.2 to 11.9 mL/day over PAGE 7 to 18 months. The paper's eq. 4/1a
    # text reads "(1 - beta*exp(...))" but numerical reconciliation requires "(1 - (1-beta)*exp(...))"
    # i.e. the Anderson/Allegaert/Holford 2006 parameterization added as reference 1a by
    # the 2012 erratum (PMC3457364).
    beta_cl <- 0.411; label("Fraction of mature CL at 40-week PAGE (term) (unitless)") # Robbie 2012 Table 2
    t50_cl  <- 62.3;  label("Maturation half-life for CL (TCL, months)")               # Robbie 2012 Table 2

    # Race effects on CL (multiplicative; white = reference). All 95% CIs include unity
    # in the source analysis (Robbie 2012 Table 2) but the estimates are reported as
    # part of the final model, so we carry them here verbatim.
    e_black_cl    <- 0.06; label("Race effect on CL: Black / African American vs white (fraction)") # Robbie 2012 Table 2 (1.06 -> 0.06)
    e_hispanic_cl <- 0.05; label("Race effect on CL: Hispanic / Latino vs white (fraction)")        # Robbie 2012 Table 2 (1.05 -> 0.05)
    e_asian_cl    <- 0.12; label("Race effect on CL: Asian vs white (fraction)")                    # Robbie 2012 Table 2 (1.12 -> 0.12)
    e_other_cl    <- 0.10; label("Race effect on CL: Other vs white (fraction)")                    # Robbie 2012 Table 2 (1.10 -> 0.10)

    # Race effect on Vc (Hispanic only reported; 95% CI includes unity).
    e_hispanic_vc <- 0.06; label("Race effect on Vc: Hispanic / Latino vs white (fraction)")        # Robbie 2012 Table 2 (1.06 -> 0.06)

    # Chronic lung disease of prematurity effect on CL (significant, +20%).
    e_cld_cl <- 0.20; label("Chronic lung disease of prematurity effect on CL (fraction)") # Robbie 2012 Table 2 (1.20 -> 0.20)

    # ADA titer category effects on CL (step function; 0 = reference).
    # Only the >=80 category CI excluded unity in the source analysis.
    e_ada10_cl   <- 0.15; label("ADA titer = 10 effect on CL (fraction)")  # Robbie 2012 Table 2 (1.15 -> 0.15)
    e_ada20_cl   <- 0.06; label("ADA titer = 20 effect on CL (fraction)")  # Robbie 2012 Table 2 (1.06 -> 0.06)
    e_ada40_cl   <- 0.08; label("ADA titer = 40 effect on CL (fraction)")  # Robbie 2012 Table 2 (1.08 -> 0.08)
    e_ada80_cl   <- 0.21; label("ADA titer >= 80 effect on CL (fraction)") # Robbie 2012 Table 2 (1.21 -> 0.21)

    # IIV: log-normal on CL and Vc (correlated, rho = 0.62). Omega variances derived from
    # reported CV% via omega^2 = log(CV^2 + 1): CL CV 48.7% -> omega^2 = 0.2128,
    # Vc CV 61.7% -> omega^2 = 0.3226, covariance = 0.62 * sqrt(0.2128*0.3226) = 0.1623.
    # Q, Vp, and ka have no reported IIV (sparse pediatric data; fixed to 0 in the final model).
    etalcl + etalvc ~ c(0.2128,
                        0.1623, 0.3226)                                    # Robbie 2012 Table 2 (IIV CL 48.7%, Vc 61.7%, rho 0.62)

    # Residual error: proportional. Paper reports sigma^2_prop = 0.0639 (Table 2);
    # propSd is the standard deviation: sqrt(0.0639) = 0.2528.
    propSd <- 0.2528; label("Proportional residual error (fraction)")      # Robbie 2012 Table 2 (sigma^2_prop = 0.0639)
  })
  model({
    # Clearance maturation centered at term birth (PAGE = 40 weeks = 9.195 months);
    # see ini() comment for the Anderson/Allegaert/Holford parameterization.
    page_centered <- PAGE - 40 / 4.35
    maturation_cl <- 1 - (1 - beta_cl) * exp(-page_centered * log(2) / t50_cl)

    # Race effects on CL (multiplicative; the reference patient is white so all
    # four indicators are zero for reference). Hispanic additionally has a Vc effect.
    race_cl <- 1 + e_black_cl * RACE_BLACK +
                    e_hispanic_cl * RACE_HISPANIC +
                    e_asian_cl * RACE_ASIAN +
                    e_other_cl * RACE_OTHER
    race_vc <- 1 + e_hispanic_vc * RACE_HISPANIC

    # CLD of prematurity effect on CL
    cld_cl <- 1 + e_cld_cl * CLD_PREM

    # ADA titer step-function effect on CL (four non-reference bins).
    # Bins: {10}, {20}, {40}, {>=80}; all other values (including 0) map to reference.
    ada_is10 <- (ADA_TITER >= 5) * (ADA_TITER < 15)
    ada_is20 <- (ADA_TITER >= 15) * (ADA_TITER < 30)
    ada_is40 <- (ADA_TITER >= 30) * (ADA_TITER < 60)
    ada_is80 <- (ADA_TITER >= 60)
    ada_cl   <- 1 + e_ada10_cl * ada_is10 +
                     e_ada20_cl * ada_is20 +
                     e_ada40_cl * ada_is40 +
                     e_ada80_cl * ada_is80

    # PK parameters with allometric weight scaling (reference 70 kg)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^allo_cl * maturation_cl * race_cl * cld_cl * ada_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^allo_v  * race_vc
    vp <- exp(lvp)          * (WT / 70)^allo_v
    q  <- exp(lq)           * (WT / 70)^allo_cl

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

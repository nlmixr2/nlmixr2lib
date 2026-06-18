Ravva_2009_varenicline <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag-time for varenicline in adult smokers (Ravva 2009): apparent clearance scales with creatinine clearance and race; central volume scales with body weight, age, and race; peripheral disposition uses fixed allometric exponents on weight."
  reference <- "Ravva P, Gastonguay MR, Tensfeldt TG, Faessel HM. Population pharmacokinetic analysis of varenicline in adult smokers. Br J Clin Pharmacol. 2009;68(5):669-681. doi:10.1111/j.1365-2125.2009.03520.x"
  vignette <- "Ravva_2009_varenicline"
  units <- list(time = "hr", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Estimated creatinine clearance by the Cockcroft-Gault formula (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F with reference 100 mL/min (Ravva 2009 Table 4 final-model q_CRCL; the reference equals the population-typical 'normal renal function' value cited in the Discussion). Estimated from serum creatinine, total body weight, age and sex via the Cockcroft-Gault formula (Ravva 2009 Methods); raw mL/min, NOT BSA-normalized. Per-paper documentation: same canonical CRCL column used in the precedent raw-Cockcroft-Gault encoding from Delattre_2010_amikacin.R. Cohort range: 15.6 to 268 mL/min, truncated at 150 to drop physiologically improbable Cockcroft-Gault upper values (Ravva 2009 Results and Table 2 footnote).",
      source_name        = "CLcr"
    ),
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on V2/F (estimated exponent), and on V3/F and Q/F (fixed allometric exponents 1 and 0.75) with reference 70 kg (Ravva 2009 Table 4 final-model q_WT and Methods/Results text on the V3/Q allometric scaling). Cohort range 41-129 kg.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on V2/F with reference 45 years (Ravva 2009 Table 4 final-model q_AGE; 45 years is the population reference cited in the Results paragraph defining the typical individual). Cohort range 18-76 years; %RSE 54% (95% CI includes 0), so the effect is poorly defined and should be viewed as exploratory.",
      source_name        = "AGE"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (reference race for Ravva 2009 typical individual)",
      notes              = "Power-of-categorical-indicator form: typical CL/F is multiplied by q_Black^RACE_BLACK and V2/F by q_Black,V2^RACE_BLACK (Ravva 2009 Table 4 q_Black; Methods 'Race entered the model as power functions with a separate dichotomous (0,1) covariate serving as an on-off switch for each category'). Cohort prevalence 12.6%.",
      source_name        = "Race (Black)"
    ),
    RACE_OTHER = list(
      description        = "Composite 'Other' race indicator pooling Hispanic, Asian, and Other (Ravva 2009 grouping)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (reference race for Ravva 2009 typical individual)",
      notes              = "Power-of-categorical-indicator form: typical CL/F is multiplied by q_Other^RACE_OTHER and V2/F by q_Other,V2^RACE_OTHER (Ravva 2009 Table 4 q_Other; Methods 'Because the number of subjects in the Hispanic, Asian and Other races was small (<=5% of total population studied), these categories were grouped together'). Composite cohort prevalence ~7% (Hispanic + Asian 1.22% + Other 5.22%).",
      source_name        = "Race (Hispanic + Asian + Other)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 1878L,
    n_studies       = 9L,
    age_range       = "18-76 years",
    age_median      = "44.2 years",
    age_mean        = "44.0 years",
    weight_range    = "41.0-129 kg",
    weight_median   = "77.0 kg",
    weight_mean     = "78.0 kg",
    height_range    = "135-202 cm",
    bmi_range       = "16.0-44.8 kg/m^2",
    crcl_range      = "15.6-268 mL/min (truncated at upper 150 mL/min for modelling)",
    crcl_median     = "107 mL/min",
    crcl_mean       = "112 mL/min",
    sex_female_pct  = 49.2,
    race_ethnicity  = c(White = 81.0, Black = 12.6, Asian = 1.22, Other = 5.22),
    disease_state   = "Adult smokers (target population for varenicline smoking-cessation therapy); about 15% with mild-to-severe renal impairment.",
    dose_range      = "0.3 to 3 mg/day oral immediate-release tablet, once-daily (q.d.) or twice-daily (b.i.d.); primary therapeutic regimen 1 mg b.i.d.",
    regions         = "Multinational; nine clinical trials pooled (4 phase-I, 2 phase-II, 3 phase-III).",
    samples         = "11935 plasma varenicline concentrations across 1878 subjects after removal of 664 BLQ records (4% of total).",
    notes           = "Baseline demographics from Ravva 2009 Table 2; study design summary from Table 1. Pooled population PK analysis across the varenicline clinical development programme."
  )

  ini({
    # Structural parameters - reference values for a White subject, 70 kg, 45 years,
    # CLcr = 100 mL/min (Ravva 2009 Results paragraph defining the 'typical individual').
    lka     <- log(1.69); label("First-order absorption rate constant (Ka, 1/h)")            # Ravva 2009 Table 4 final-model Ka
    lcl     <- log(10.4); label("Apparent clearance (CL/F, L/h) at reference covariates")   # Ravva 2009 Table 4 final-model q_CL
    lvc     <- log(337);  label("Apparent central volume of distribution (V2/F, L)")        # Ravva 2009 Table 4 final-model q_V2
    lvp     <- log(78.1); label("Apparent peripheral volume of distribution (V3/F, L)")     # Ravva 2009 Table 4 final-model q_V3
    lq      <- log(2.08); label("Apparent intercompartmental clearance (Q/F, L/h)")         # Ravva 2009 Table 4 final-model q_Q
    ltlag   <- log(0.43); label("Absorption lag time (Alag, h)")                              # Ravva 2009 Table 4 final-model q_Alag

    # Covariate effects on CL/F: power exponent on CRCL/100 plus categorical
    # race multipliers in power-of-indicator form (Ravva 2009 Table 4).
    e_crcl_cl   <- 0.54;  label("Power exponent of CRCL/100 on CL/F (unitless)")              # Ravva 2009 Table 4 final-model q_CRCL on CL/F
    e_black_cl  <- 1.16;  label("Multiplicative factor on CL/F for Black vs White (unitless)") # Ravva 2009 Table 4 final-model q_Black on CL/F
    e_other_cl  <- 1.11;  label("Multiplicative factor on CL/F for Other vs White (unitless)") # Ravva 2009 Table 4 final-model q_Other on CL/F

    # Covariate effects on V2/F: power exponents on WT/70 and AGE/45, plus
    # categorical race multipliers (Ravva 2009 Table 4).
    e_wt_vc     <- 0.77;  label("Power exponent of WT/70 on V2/F (unitless)")                  # Ravva 2009 Table 4 final-model q_WT on V2/F
    e_age_vc    <- 0.13;  label("Power exponent of AGE/45 on V2/F (unitless)")                 # Ravva 2009 Table 4 final-model q_AGE on V2/F
    e_black_vc  <- 0.92;  label("Multiplicative factor on V2/F for Black vs White (unitless)") # Ravva 2009 Table 4 final-model q_Black on V2/F
    e_other_vc  <- 0.71;  label("Multiplicative factor on V2/F for Other vs White (unitless)") # Ravva 2009 Table 4 final-model q_Other on V2/F

    # Allometric exponents on V3/F and Q/F (Ravva 2009 Table 4 reports
    # these as 'Fixed'; held at the canonical literature values).
    e_wt_vp     <- fixed(1);    label("Allometric exponent of WT/70 on V3/F (unitless)")        # Ravva 2009 Table 4 q_WT on V3/F = 1 (Fixed)
    e_wt_q      <- fixed(0.75); label("Allometric exponent of WT/70 on Q/F (unitless)")         # Ravva 2009 Table 4 q_WT on Q/F = 0.75 (Fixed)

    # IIV - log-normal on Ka, CL, V2 with a full BLOCK omega; IIV on V3 and Q
    # was fixed to 0 (Ravva 2009 Table 4 Interindividual variance block).
    # Row/column order: Ka, CL, V2; entries are the lower-triangular variances and
    # covariances reported in Table 4 (variances and absolute covariances; the r
    # values reported in the table are the standardised correlations).
    etalka + etalcl + etalvc ~ c(
       0.49,
      -0.009, 0.061,
       0.24,  0.006, 0.25
    )

    # Residual error (primary): Studies I (volunteer settings, intensive sampling).
    # Combined additive + proportional on plasma varenicline (ng/mL). Ravva 2009
    # Table 4 reports s^2_add,1 = 0.28 (SD = 0.5 ng/mL) and s^2_prop,1 = 0.030
    # (17.2% CV). The same final model also reports a distinct error set for the
    # large patient trials (Studies II + III): s^2_add,2 = 4.38 (SD = 2.1 ng/mL),
    # s^2_prop,2 = 0.046 (21.5% CV). The Studies-I values are encoded here as
    # primary; the Studies II + III values are documented in the vignette as an
    # alternative-deviation option. See vignette Assumptions and deviations.
    addSd  <- 0.5;   label("Additive residual error on Cc (ng/mL)")        # Ravva 2009 Table 4 s^2_add,1 = 0.28; reported SD = 0.5 ng/mL
    propSd <- 0.172; label("Proportional residual error on Cc (fraction)") # Ravva 2009 Table 4 s^2_prop,1 = 0.030; reported CV = 17.2%
  })
  model({
    # Individual PK parameters.
    # Apparent clearance: power on CRCL/100 plus power-of-indicator race effects.
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl) *
           (CRCL / 100)^e_crcl_cl *
           e_black_cl^RACE_BLACK *
           e_other_cl^RACE_OTHER
    vc  <- exp(lvc + etalvc) *
           (WT / 70)^e_wt_vc *
           (AGE / 45)^e_age_vc *
           e_black_vc^RACE_BLACK *
           e_other_vc^RACE_OTHER
    vp  <- exp(lvp) * (WT / 70)^e_wt_vp
    q   <- exp(lq)  * (WT / 70)^e_wt_q

    # Absorption lag-time (Ravva 2009 final model includes Alag on the depot
    # compartment).
    tlag <- exp(ltlag)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE system with first-order absorption from depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Plasma concentration (dose in mg, volumes in L -> mg/L; report in ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

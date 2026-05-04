Lin_2024_pozelimab <- function() {
  description <- "Two-compartment two-binding-site TMDD-QE population PK model of total pozelimab and total C5 in healthy volunteers, adults with paroxysmal nocturnal hemoglobinuria, and pediatric and adult patients with CHAPLE disease (Lin 2024)"
  reference <- "Lin K-J, Mendell J, Davis JD, Harnisch LO. Population pharmacokinetic analyses of pozelimab in patients with CD55-deficient protein-losing enteropathy (CHAPLE disease). J Pharmacokinet Pharmacodyn. 2024;51(6):905-917. doi:10.1007/s10928-024-09941-8"
  vignette <- "Lin_2024_pozelimab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; time-varying body weight evaluated as sensitivity analysis only and not retained in the final model)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, Vc, and Vp using (WT/70)^exponent. Reference weight 70 kg is the cohort-median adult body weight (Table 1). Exponents are estimated, not fixed at allometric values: 0.9989 on CL, 0.7560 on Vc and Vp.",
      source_name        = "WT"
    ),
    DIS_PNH = list(
      description        = "Paroxysmal nocturnal hemoglobinuria patient status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-PNH: healthy volunteer or CHAPLE disease patient)",
      notes              = "Additive-fractional effect on Vc: Vc = Vc_TV * (1 + e_pnh_vc * DIS_PNH); PNH patients have a 34.07% larger Vc than non-PNH subjects. The supplement's Supplemental Text 1 specifies only the structural ODEs and not the categorical-covariate equation form; the additive-fractional form is assumed because Table 2's 95% CI for the estimate (0.156-0.5254) is bracketed cleanly above zero and the alternative power form theta^DIS_PNH (= 0.3407) implies an implausible 66% reduction in Vc.",
      source_name        = "PNH (disease state, PNH versus healthy volunteer)"
    )
  )

  population <- list(
    n_subjects     = 116L,
    n_studies      = 4L,
    age_range      = "3-76 years (median 37 across all subjects; pediatric and adult)",
    age_median     = "37 years overall; 8.5 years in CHAPLE patients (range 3-19)",
    weight_range   = "11.0-108 kg overall; 11.0-53.8 kg in CHAPLE patients",
    weight_median  = "66.7 kg overall; 25.0 kg in CHAPLE patients",
    sex_female_pct = 60.3,
    race_ethnicity = "70.7% White, 22.4% Asian, 3.4% Black or African American, 0.9% American Indian or Alaska Native, 2.6% Other (Table 1).",
    disease_state  = "Healthy adult volunteers (n=82), adult patients with paroxysmal nocturnal hemoglobinuria (PNH; n=24), and pediatric and adult patients with CHAPLE disease (n=10; 9 children, 1 adult).",
    dose_range     = "1-30 mg/kg single IV; 300-600 mg single SC; 400 mg SC QW; 30 mg/kg IV loading + 800 mg SC QW (PNH study); 30 mg/kg IV loading + weight-tiered 125-800 mg SC QW (CHAPLE pediatric).",
    regions        = "Multi-regional phase 1-3 programme; specific regional breakdown not reported in the main text.",
    ada_status     = "All 116 subjects ADA-negative (0% positive in every cohort).",
    notes          = "Pooled phase 1 first-in-human (NCT03115996, n=42), phase 1 PK comparability (NCT04491838, n=40), phase 2 PNH (NCT03946748, n=24), and phase 2/3 CHAPLE (NCT04209634, n=10) trials. 2795 concentration samples (1640 total pozelimab, 1155 total C5). 100 (3.6%) post-dose BLQ samples excluded per Beal M1. Pozelimab MW assumed 150 kDa (typical IgG4); C5 MW 190 kDa (matches the paper's footnote conversion ksyn = 0.04922 uM/day = 9.352 mg/L/day and kD = 0.000189 uM = 0.03591 mg/L)."
  )

  ini({
    # ---------------------------------------------------------------------------
    # Structural PK parameters (Lin 2024 Table 2 final TMDD PopPK model)
    # ---------------------------------------------------------------------------
    lcl     <- log(0.1506); label("Linear clearance CL (L/day)")                   # Lin 2024 Table 2
    lvc     <- log(2.476);  label("Central volume of distribution Vc (L)")          # Lin 2024 Table 2
    lq      <- log(0.3931); label("Inter-compartmental clearance Q (L/day)")         # Lin 2024 Table 2
    lvp     <- log(9.901);  label("Peripheral volume of distribution Vp (L)")        # Lin 2024 Table 2
    lka     <- log(0.1726); label("First-order subcutaneous absorption rate ka (1/day)")  # Lin 2024 Table 2
    lfdepot <- log(0.6864); label("Subcutaneous bioavailability F1 (fraction)")     # Lin 2024 Table 2

    # ---------------------------------------------------------------------------
    # TMDD-QE target parameters (Lin 2024 Table 2; molar units uM as in
    # Supplemental Text 1; mass-equivalent values from Table 2 footnote (a))
    # ---------------------------------------------------------------------------
    lkD    <- fixed(log(0.000189));   label("Equilibrium dissociation constant kD (uM); FIXED per SPR-Biacore measurement [paper footnote: 0.03591 mg/L]")  # Lin 2024 Table 2
    lksyn  <- log(0.04922);           label("Synthesis rate of free C5 ksyn (uM/day) [paper footnote: 9.352 mg/L/day]")  # Lin 2024 Table 2
    lkdeg  <- log(0.1105);            label("Degradation rate constant of free C5 kdeg (1/day)")  # Lin 2024 Table 2
    lkint1 <- log(0.08086);           label("Internalization rate of pozelimab-C5 complex kint1 (1/day)")  # Lin 2024 Table 2
    lkint2 <- log(0.1149);            label("Internalization rate of pozelimab-C5-C5 complex kint2 (1/day)")  # Lin 2024 Table 2

    # Correction factor for baseline total C5 (theta_13 in the paper text);
    # R_tot(0) = ksyn/kdeg + theta_R0. Negative value reflects an early
    # (~1 hour post-dose) drop in observed total C5 not predicted by the
    # bare turnover steady state.
    theta_R0 <- -0.04824; label("Correction factor for baseline C5 (uM); R_tot(0) = ksyn/kdeg + theta_R0 [paper footnote: -9.166 mg/L]")  # Lin 2024 Table 2

    # ---------------------------------------------------------------------------
    # Covariate effects (Lin 2024 Table 2)
    # ---------------------------------------------------------------------------
    # Body weight on CL: power form CL = CL_TV * (WT/70)^e_wt_cl
    # Body weight on Vc and Vp (single shared exponent): Vc = Vc_TV * (WT/70)^e_wt_vc_vp
    # PNH on Vc: additive-fractional, Vc = Vc_TV * (1 + e_pnh_vc * DIS_PNH)
    e_wt_cl    <- 0.9989; label("Power exponent of body weight on CL (unitless; reference 70 kg)")  # Lin 2024 Table 2 ("Weight on CL")
    e_wt_vc_vp <- 0.7560; label("Shared power exponent of body weight on Vc and Vp (unitless; reference 70 kg)")  # Lin 2024 Table 2 ("Weight on Vc and Vp")
    e_pnh_vc   <- 0.3407; label("Additive-fractional effect of PNH disease state on Vc (unitless)")  # Lin 2024 Table 2 ("Patients with PNH on Vc")

    # ---------------------------------------------------------------------------
    # Inter-individual variability (Lin 2024 Table 2; %CV)
    # Convert to log-normal variance via omega^2 = log(CV^2 + 1):
    #   31.90 % -> log(0.319^2 + 1)  = 0.09732 (CL)
    #   16.10 % -> log(0.161^2 + 1)  = 0.02568 (Vc)
    #   80.37 % -> log(0.8037^2 + 1) = 0.48107 (Vp)
    #   16.48 % -> log(0.1648^2 + 1) = 0.02683 (ksyn)
    #   21.10 % -> log(0.211^2 + 1)  = 0.04362 (kint1)
    # ---------------------------------------------------------------------------
    etalcl    ~ 0.09732   # Lin 2024 Table 2: IIV on CL    = 31.90 %CV
    etalvc    ~ 0.02568   # Lin 2024 Table 2: IIV on Vc    = 16.10 %CV
    etalvp    ~ 0.48107   # Lin 2024 Table 2: IIV on Vp    = 80.37 %CV
    etalksyn  ~ 0.02683   # Lin 2024 Table 2: IIV on ksyn  = 16.48 %CV
    etalkint1 ~ 0.04362   # Lin 2024 Table 2: IIV on kint1 = 21.10 %CV

    # ---------------------------------------------------------------------------
    # Residual error (Lin 2024 Table 2). Pozelimab uses a log-transformed
    # additive error (ln(Y) = ln(C_tot) + eps), which is equivalent to
    # proportional in linear space; total C5 uses a proportional error
    # directly. The paper estimated separate %CV terms for adults vs.
    # children with CHAPLE disease; this model file applies the adult
    # residual variabilities. The pediatric values (pozelimab 23.17 %CV;
    # C5 14.83 %CV) are documented in the validation vignette under
    # Assumptions and deviations.
    # ---------------------------------------------------------------------------
    propSd     <- 0.3722; label("Proportional residual error for total pozelimab in adults (fraction)")  # Lin 2024 Table 2 ("RV for pozelimab in adults" = 37.22 %CV)
    propSd_Rtot  <- 0.1072; label("Proportional residual error for total C5 in adults (fraction)")          # Lin 2024 Table 2 ("RV for C5 in adults"        = 10.72 %CV)
  })

  model({
    # ---------------------------------------------------------------------------
    # Molecular weights (kDa = g/mol divided by 1000). Pozelimab is a fully
    # human IgG4 (typical 150 kDa); C5 MW (190 kDa) reproduces the paper's
    # own uM <-> mg/L conversions in the Table 2 footnote. The state space
    # runs in molar (uM / umol) units to match Supplemental Text 1; MW
    # converts mg dose -> umol amount and umol/L state -> mg/L observation.
    # ---------------------------------------------------------------------------
    MW_drug_kDa   <- 150
    MW_target_kDa <- 190
    WT_ref        <- 70

    # ---------------------------------------------------------------------------
    # Individual PK parameters with covariates
    # ---------------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT/WT_ref)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT/WT_ref)^e_wt_vc_vp * (1 + e_pnh_vc * DIS_PNH)
    vp <- exp(lvp + etalvp) * (WT/WT_ref)^e_wt_vc_vp
    q  <- exp(lq)
    ka <- exp(lka)

    # ---------------------------------------------------------------------------
    # Individual TMDD-QE parameters
    # ---------------------------------------------------------------------------
    kD    <- exp(lkD)                   # uM
    ksyn  <- exp(lksyn + etalksyn)      # uM/day
    kdeg  <- exp(lkdeg)                 # 1/day
    kint1 <- exp(lkint1 + etalkint1)    # 1/day
    kint2 <- exp(lkint2)                # 1/day

    # ---------------------------------------------------------------------------
    # QE algebraic relations (Supplemental Text 1)
    #   R = 0.5 * [-(2*Ctot + kD - Rtot) + sqrt((2*Ctot + kD - Rtot)^2 + 4*kD*Rtot)]
    #   Cu  = Ctot * kD^2 / (kD + R)^2
    #   RC  = Ctot * 2*kD*R / (kD + R)^2
    #   R2C = Ctot * R^2 / (kD + R)^2
    # All concentrations in uM; central state in umol; total_target state in uM.
    # ---------------------------------------------------------------------------
    Ctot   <- central / vc
    discR  <- 2*Ctot + kD - total_target
    R      <- 0.5 * (-discR + sqrt(discR*discR + 4*kD*total_target))
    denom2 <- (kD + R)^2
    Cu     <- Ctot * kD * kD / denom2

    # ---------------------------------------------------------------------------
    # Initial condition for total C5 (uM) per paper R0 = ksyn/kdeg + theta_R0.
    # No drug dosing prior to t = 0, so total_target(0) = R(0) = R0.
    # ---------------------------------------------------------------------------
    total_target(0) <- ksyn/kdeg + theta_R0

    # ---------------------------------------------------------------------------
    # Bioavailability and dose -> umol amount conversions (mg / MW_kDa = umol).
    # f(depot)   covers SC dosing (apparent F1 multiplied into the depot input)
    # f(central) covers IV dosing (no F adjustment; just mg -> umol scaling)
    # ---------------------------------------------------------------------------
    f(depot)   <- exp(lfdepot) / MW_drug_kDa
    f(central) <- 1 / MW_drug_kDa

    # ---------------------------------------------------------------------------
    # ODE system (Supplemental Text 1; central rewritten in umol amount space:
    # multiply Supplemental Text 1's dC_tot/dt equation by Vc).
    # ---------------------------------------------------------------------------
    d/dt(depot)       <- -ka*depot
    d/dt(central)     <-  ka*depot + (q/vp)*peripheral1 - (cl + q)*Cu -
                          vc * (Ctot*R/denom2) * (2*kint1*kD + kint2*R)
    d/dt(peripheral1) <-  q*Cu - (q/vp)*peripheral1
    d/dt(total_target) <- ksyn - kdeg*total_target -
                          (2*Ctot*R/denom2) * (kD*(kint1 - kdeg) + R*(kint2 - kdeg))

    # ---------------------------------------------------------------------------
    # Observation variables (mg/L). Cc and Rtot are the bioassay measurements
    # (total pozelimab and total C5). Free pozelimab and free C5 are exported
    # for simulation use (they are not observation-fit endpoints).
    # ---------------------------------------------------------------------------
    Cc           <- Ctot           * MW_drug_kDa
    Rtot         <- total_target   * MW_target_kDa
    free_drug    <- Cu             * MW_drug_kDa
    free_target  <- R              * MW_target_kDa

    # ---------------------------------------------------------------------------
    # Residual-error models
    # ---------------------------------------------------------------------------
    Cc   ~ prop(propSd)
    Rtot ~ prop(propSd_Rtot)
  })
}

Ng_2018_exendin939 <- function() {
  description <- "Two-compartment intravenous-infusion population PK model for exendin-(9-39) in patients with congenital hyperinsulinism (Ng 2018). Pooled paediatric (neonates and children) and adult cohort with allometric scaling fixed at 0.75 on CL and Q and 1.0 on Vc and Vp (reference WT 70 kg); inter-individual variability retained only on CL. Residual variability follows the NONMEM Poisson error model (Var(Y|F) = F * sigma^2), encoded as a power-error with fixed exponent 0.5."
  reference <- "Ng CM, Tang F, Seeholzer SH, Zou Y, De Leon DD. Population pharmacokinetics of exendin-(9-39) and clinical dose selection in patients with congenital hyperinsulinism. Br J Clin Pharmacol. 2018;84(3):520-528. doi:10.1111/bcp.13463"
  vignette <- "Ng_2018_exendin939"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject for this dataset (single-occasion analyses). Allometric scaling with reference weight 70 kg; exponents fixed at 0.75 on CL and Q and 1.0 on Vc and Vp. Ng 2018 Table 1: paediatric weight median 6.25 kg (neonates, range 3.92-6.60) and 21.0 kg (children, range 11.8-69.2); adult median 69.1 kg (range 58.2-130).",
      source_name        = "WT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 26L,
    n_studies        = 3L,
    age_range        = "0.06-47 years",
    age_median       = "Neonates 0.11 years; children 5 years; adults 19 years (Table 1)",
    weight_range     = "3.92-130 kg",
    weight_median    = "Neonates 6.25 kg; children 21.0 kg; adults 69.1 kg (Table 1)",
    sex_female_pct   = 53.8,
    race_ethnicity   = "Caucasian 92%; unknown 8% (Ng 2018 Table 1: 24/26 Caucasian, 2 unknown across the three cohorts)",
    disease_state    = "Congenital hyperinsulinism (KATP-HI subtype predominantly), spanning neonates (<12 months), children (6 months to 18 years), and older adolescents / adults (>=15 years).",
    dose_range       = "Intravenous infusion 0.02-0.1 mg/kg/h; cohorts received 2 h infusions per dose level for 6 h total (adults, children) or 6-12 h single-rate infusions (neonates).",
    regions          = "USA (Children's Hospital of Philadelphia)",
    n_concentrations = 182L,
    n_below_lloq     = 16L,
    studies          = "NCT00571324 (Adult Study); NCT00897676 (Children Study); NCT00835328 (Neonate Study)",
    notes            = "Baseline demographics from Ng 2018 Table 1. Bioanalytical: LC-MS/MS, linear range 10-1390 ng/mL, LLOQ 10 ng/mL (limit of detection 1.3 ng/mL). Below-LLOQ samples (16/182) were retained via the Beal M3 method during the original NONMEM 7.3 fit (FOCE-I). PopPK is human-only; preclinical TK in rats (Sprague-Dawley) and beagle dogs was used for NOAEL / dose-selection only and is not part of this model."
  )

  ini({
    # Structural parameters (Ng 2018 Table 2 final-model column; reference adult
    # body weight 70 kg). Log-transformed for positivity.
    lcl <- log(11.6); label("Typical clearance at 70 kg (CL_TV, L/h)")               # Ng 2018 Table 2: CL_TV = 11.6 L/h (SE 0.728)
    lvc <- log(9.59); label("Typical central volume at 70 kg (V1_TV, L)")            # Ng 2018 Table 2: V1_TV = 9.59 L (SE 0.449)
    lq  <- log(2.20); label("Typical inter-compartmental clearance at 70 kg (Q_TV, L/h)") # Ng 2018 Table 2: Q_TV  = 2.20 L/h (SE 0.576)
    lvp <- log(8.89); label("Typical peripheral volume at 70 kg (V2_TV, L)")         # Ng 2018 Table 2: V2_TV = 8.89 L (SE 0.996)

    # Allometric weight effects -- fixed by the authors at the physiology-based
    # canonical exponents (0.75 on clearances, 1.0 on volumes) per Ng 2018
    # Methods ("theta_a was ... fixed to 0.75 for clearance and 1 for volumes
    # based on physiologic consideration of the impact of size on metabolic
    # rate") and Table 2 ("fixed" annotation on both rows).
    e_wt_cl_q  <- fixed(0.75); label("Allometric WT exponent shared across CL and Q (unitless)") # Ng 2018 Table 2: fixed
    e_wt_vc_vp <- fixed(1.00); label("Allometric WT exponent shared across Vc and Vp (unitless)") # Ng 2018 Table 2: fixed

    # Inter-individual variability. Ng 2018 retained IIV only on CL: V1, Q, and
    # V2 IIVs failed convergence or yielded %CV > 50% and were dropped from the
    # final model (Results paragraph "Including interindividual variability
    # terms for V1, Q and V2 failed to achieve model convergence ...").
    etalcl ~ 0.0572 # Ng 2018 Table 2: omega^2_CL = 0.0572 (SE 0.0262); ~24% CV

    # Residual variability -- NONMEM Poisson error model. Standard NONMEM
    # encoding for continuous data: W = sqrt(F); Y = F + W * EPS(1), so
    # Var(Y|F) = F * sigma^2 and SD(Y|F) = sigma * sqrt(F). Equivalent to a
    # power-error with fixed exponent 0.5; ini parameter powSd is the
    # square-root-scale SD (sqrt(Sigma_poisson)). Ng 2018 Table 2 reports
    # Sigma_poisson = 3.66 (the EPS variance; SE 0.508; bootstrap 95% CI
    # 2.73-4.47), so powSd = sqrt(3.66) ~= 1.913. The exponent powExp is held
    # at 0.5 (fixed) to lock in the Poisson interpretation; rxode2 requires
    # the exponent as an ini parameter (it does not accept a numeric literal
    # in `pow(..., 0.5)` during forward simulation). See vignette Assumptions
    # and deviations section.
    powSd  <- sqrt(3.66); label("Poisson residual error scale (SD = powSd * sqrt(Cc); (ng/mL)^0.5)") # Ng 2018 Table 2: Sigma_poisson = 3.66 (SE 0.508)
    powExp <- fixed(0.5); label("Power exponent on Cc in the residual error (fixed for NONMEM Poisson interpretation)") # Ng 2018 Methods: Poisson error model implies W = sqrt(F)
  })
  model({
    # Two-compartment IV-infusion PK with allometric scaling on all four
    # structural parameters (Ng 2018 Equations 6-9: CL = CL_TV * (WT/70)^0.75,
    # V1 = V1_TV * (WT/70)^1, Q = Q_TV * (WT/70)^0.75, V2 = V2_TV * (WT/70)^1).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc)          * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV bolus / infusion goes directly to central; no depot for this dataset.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> central/vc has units mg/L = ug/mL = 1000 * ng/mL.
    # Multiply by 1000 to express concentration in ng/mL, matching the paper's
    # bioanalytical units (LC-MS/MS calibration 10-1390 ng/mL).
    Cc <- 1000 * central / vc

    # Poisson error model: SD(Y|Cc) = powSd * Cc^powExp with powExp = 0.5
    # (fixed). See ini() comment.
    Cc ~ pow(powSd, powExp)
  })
}

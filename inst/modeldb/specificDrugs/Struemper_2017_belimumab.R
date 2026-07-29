Struemper_2017_belimumab <- function() {
  description <- "Linear two-compartment subcutaneous population PK model for belimumab in healthy volunteers and adult patients with systemic lupus erythematosus, with first-order absorption + lag time, allometric body-weight scaling on CL/Vc/Q/Vp, and baseline BMI on Vc and baseline albumin and IgG on CL (Struemper 2017)"
  reference <- "Struemper H, Thapar M, Roth D. Population Pharmacokinetic and Pharmacodynamic Analysis of Belimumab Administered Subcutaneously in Healthy Volunteers and Patients with Systemic Lupus Erythematosus. Clin Pharmacokinet. 2018;57(6):717-728. doi:10.1007/s40262-017-0586-5"
  vignette <- "Struemper_2017_belimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight (BWT in the source NONMEM code); fixed allometric power effects on CL (0.75), Vc (1.00), Q (0.75), and Vp (0.8) with reference 67 kg (population median).",
      source_name        = "BWT"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline BMI (BBMI in the source NONMEM code); estimated power effect on Vc with reference 24.7 kg/m^2 (population median).",
      source_name        = "BBMI"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline albumin (BALB in the source NONMEM code); estimated power effect on CL with reference 41 g/L (population median). Renamed from source column BALB to canonical ALB per covariate-columns.md.",
      source_name        = "BALB"
    ),
    IGG = list(
      description        = "Baseline serum immunoglobulin G",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline IgG (BIGG in the source NONMEM code); estimated power effect on CL with reference 13.7 g/L (population median). Renamed from source column BIGG to canonical IGG per covariate-columns.md.",
      source_name        = "BIGG"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 688,
    n_studies       = 3,
    n_observations  = 4958,
    age_range       = "18-77 years",
    age_median      = "37 years",
    age_mean        = "37.6 (SD 11.8) years",
    weight_range    = "34.1-138 kg",
    weight_median   = "67.0 kg",
    weight_mean     = "70.0 (SD 17.6) kg",
    bmi_range       = "14.8-72.7 kg/m^2",
    bmi_median      = "24.7 kg/m^2",
    sex_female_pct  = 85,
    race_ethnicity  = c(White = 61, Asian = 23.5, Black = 11, AIAN = 6, Other = 1.7),
    disease_state   = "Pooled across 3 studies: 134 healthy volunteers (two phase I studies in Japan and the USA) and 554 adult patients with active autoantibody-positive systemic lupus erythematosus (one phase III global trial).",
    dose_range      = "Single 200 mg SC, four weekly 200 mg SC doses (phase I), or 200 mg SC weekly for 51 weeks (phase III BLISS-SC).",
    regions         = "Multinational pooled analysis (Japan, USA, global phase III).",
    baseline_labs   = "Mean (SD): IgG 14.9 (5.43) g/L; albumin 40.7 (4.47) g/L; creatinine clearance 115 (38.3) mL/min; haemoglobin 126 (15.8) g/L; WBC 6.09 (2.41) Gi/L.",
    studies         = "BEL114448 (NCT01583530, phase I), BEL116119 (NCT01516450, phase I), BEL112341 (NCT01484496, BLISS-SC phase III).",
    notes           = "Demographics from Struemper 2017 Table 2 (population-pharmacokinetic analyses: N = 688). Belimumab is an IgG1 monoclonal antibody targeting B-lymphocyte stimulator (BLyS); no evidence for substantial target-mediated disposition was detected at therapeutic exposures."
  )

  ini({
    # Structural parameters -- typical values from Struemper 2017 Table 3 (final SC popPK model)
    # Reported in mL / mL per day; converted to L / L per day by dividing by 1000.
    lka      <- log(0.235);    label("Absorption rate constant (Kabs, 1/day)")                                    # Table 3 THETA(5)
    lcl      <- log(0.204);    label("Clearance for the reference adult (CL, L/day)")                             # Table 3 THETA(1): 204 mL/day
    lvc      <- log(2.300);    label("Central volume of distribution for the reference adult (Vc, L)")            # Table 3 THETA(2): 2300 mL
    lq       <- log(0.698);    label("Intercompartmental clearance (Q, L/day)")                                   # Table 3 THETA(3): 698 mL/day
    lvp      <- log(2.650);    label("Peripheral volume of distribution (Vp, L)")                                 # Table 3 THETA(4): 2650 mL
    lfdepot  <- log(0.742);    label("Subcutaneous bioavailability (F, fraction)")                                # Table 3 THETA(6)
    lalag    <- log(0.179);    label("Absorption lag time (ALAG, day)")                                           # Table 3 THETA(7)

    # Fixed allometric body-weight exponents (reference 67 kg)
    # Effect column in Table 3 lists each as 9 (BWT/67)^e with no estimate ("-"), i.e. fixed.
    # The Vp exponent (0.8) is printed in the Table 3 row header "V_p [mL] 0.8";
    # see vignette Assumptions and deviations for the typesetting note.
    e_wt_cl  <- fixed(0.75);   label("Allometric exponent of BWT on CL (unitless)")                               # Table 3 fixed
    e_wt_vc  <- fixed(1.00);   label("Allometric exponent of BWT on Vc (unitless)")                               # Table 3 fixed; text Section 3.1.2
    e_wt_q   <- fixed(0.75);   label("Allometric exponent of BWT on Q (unitless)")                                # Table 3 fixed
    e_wt_vp  <- fixed(0.8);    label("Allometric exponent of BWT on Vp (unitless)")                               # Table 3 row header "V_p [mL] 0.8" (fixed)

    # Estimated covariate effects on CL and Vc (Struemper 2017 Table 3)
    e_alb_cl <- -0.736;  label("Power exponent of baseline albumin on CL (reference 41 g/L)")                     # Table 3 BALB effect on CL
    e_igg_cl <-  0.347;  label("Power exponent of baseline IgG on CL (reference 13.7 g/L)")                       # Table 3 BIGG effect on CL
    e_bmi_vc <- -0.610;  label("Power exponent of baseline BMI on Vc (reference 24.7 kg/m^2)")                    # Table 3 BBMI effect on Vc

    # Inter-individual variability (Struemper 2017 Table 3)
    # OMEGA(1,1) on CL, OMEGA(2,2) on V1 (Vc), OMEGA(2,1) = covariance(CL,V1); no Ka / F / ALAG IIV reported.
    etalcl + etalvc ~ c(0.0910,
                        0.0630, 0.497)   # Table 3 OMEGA(1,1), OMEGA(2,1), OMEGA(2,2)
    etalq  ~ 1.07                         # Table 3 OMEGA(3,3) on Q (diagonal)
    etalvp ~ 0.110                        # Table 3 OMEGA(4,4) on V2/Vp (diagonal)

    # Residual error -- combined proportional + additive (Struemper 2017 Table 3, SIGMA)
    # SIGMA(1) = 0.0327 is the proportional variance -> SD = sqrt(0.0327) = 0.1808 (fraction).
    # SIGMA(2) = 0.134 is the additive variance -> SD = sqrt(0.134) = 0.3661 ug/mL.
    propSd <- 0.1808;    label("Proportional residual error (fraction)")                                          # Table 3 SIGMA(1) = 0.0327 (variance)
    addSd  <- 0.3661;    label("Additive residual error (ug/mL)")                                                 # Table 3 SIGMA(2) = 0.134 (variance)
  })

  model({
    # Individual PK parameters with allometric BWT scaling and estimated covariate effects
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * (WT / 67)^e_wt_cl * (ALB / 41)^e_alb_cl * (IGG / 13.7)^e_igg_cl
    vc   <- exp(lvc + etalvc) * (WT / 67)^e_wt_vc * (BMI / 24.7)^e_bmi_vc
    q    <- exp(lq  + etalq)  * (WT / 67)^e_wt_q
    vp   <- exp(lvp + etalvp) * (WT / 67)^e_wt_vp
    fsc  <- exp(lfdepot)
    alag <- exp(lalag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Linear two-compartment ODE with first-order SC absorption from depot (NONMEM ADVAN3 TRANS4)
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability and absorption lag applied at the depot
    f(depot)    <- fsc
    alag(depot) <- alag

    # Concentration: dose in mg, vc in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}

Chen_2014_nab_paclitaxel <- function() {
  description <- "Three-compartment population PK with saturable (Michaelis-Menten) distribution between the central and first peripheral compartments and saturable elimination from the central compartment, coupled with a Friberg-style 5-compartment semi-mechanistic PD model for paclitaxel-induced neutropenia, fit to 150 adult patients with advanced solid tumors who received nab-paclitaxel (Abraxane) 80-375 mg/m^2 as a 30-minute IV infusion (Chen 2014). The first peripheral compartment exchanges with central via the saturable Vmtr / Kmtr process; the second peripheral compartment exchanges via linear intercompartmental clearance Q2. Baseline albumin lowers the maximal elimination rate VMEL via a power-form covariate; advanced age (>= 65 years) potentiates the linear Slope of paclitaxel-driven inhibition of the proliferating neutrophil precursor pool, and baseline albumin also modifies the baseline ANC via a power-form covariate. PK observation is plasma paclitaxel concentration (ug/L = ng/mL); PD observation is absolute neutrophil count (10^9 cells/L)."
  reference <- paste(
    "Chen N, Li Y, Ye Y, Palmisano M, Chopra R, Zhou S. (2014).",
    "Pharmacokinetics and Pharmacodynamics of nab-Paclitaxel in Patients With Solid Tumors:",
    "Disposition Kinetics and Pharmacology Distinct From Solvent-Based Paclitaxel.",
    "The Journal of Clinical Pharmacology 54(10):1097-1107.",
    "doi:10.1002/jcph.304.",
    "PD structure follows Friberg LE et al. (2002) J Clin Oncol 20(24):4713-4721",
    "(see modellib('Friberg_2002_paclitaxel') for the leukocyte arm of the original).",
    sep = " "
  )
  vignette <- "Chen_2014_nab_paclitaxel"
  units <- list(
    time          = "hour",
    dosing        = "ug",
    concentration = "ug/L",
    anc           = "10^9 cells/L"
  )

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin concentration (baseline).",
      units = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters two parameters as a power-form effect centred on the population median (Chen 2014 Methods, 'Continuous covariates were centered to their median values and included as power models'). On VMEL the reference value is the PK-population median 3.9 g/dL (Table 1, N = 150): VMEL_i = VMEL * (ALB/3.9)^0.554. On baseline ANC the reference value is the PD-population median 4.0 g/dL (Table 1, N = 125): BASE_i = BASE * (ALB/4.0)^(-0.998). The two power exponents and reference values come from the final PK and PD models in Table 2. The paper notes (Results, Population PD Model -> Covariate analysis) that the physiological and clinical relevance of the baseline-ANC / albumin correlation is unclear but the covariate improved the goodness-of-fit (P < 0.001) and was retained in the final model. Cohort range 2.1-4.7 g/dL.",
      source_name        = "ALB"
    ),
    AGE = list(
      description        = "Subject age in years (baseline).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only via the dichotomous indicator AGE >= 65, derived inside model() as `age_gte65 <- (AGE >= 65)`. Chen 2014 evaluated age on the PD Slope both as a continuous variable and as a binary cutoff at 65 years; the binary form was selected on the basis of greater statistical significance (dOFV = -12.8 vs -9.9) and consistency with clinical practice (Results, Population PD Model -> Covariate analysis). The effect is multiplicative on the linear drug-effect slope: Slope_i = Slope * (1 + 0.501 * (AGE >= 65)). At the typical-value level Slope = 0.00253 ng/mL^-1 for patients < 65 years and 0.00380 ng/mL^-1 for patients >= 65 years (Table 2 footnote and Results text). Cohort range 24-85 years; ~57 yr median.",
      source_name        = "AGE"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 150L,
    n_studies       = 8L,
    age_range       = "24-85 years",
    age_median      = "57 years",
    weight_range    = "40-143 kg",
    weight_median   = "74 kg",
    sex_female_pct  = 60,
    race_ethnicity  = c(White = 91, Asian = 9, Black = 1, Other = 0),
    disease_state   = "Advanced or metastatic solid tumors (breast cancer 16%, melanoma 29%, other solid tumors 55%). Five Phase I studies (one of which enrolled patients with hepatic impairment), one Phase II study, and two Phase III studies pooled per Chen 2014 Methods and Supplemental Table S1. Baseline laboratory characteristics include median total bilirubin 8.6 umol/L (range 3.4-81.0; 13.3% above ULN of 17 umol/L; 8% with moderate-to-severe hepatic impairment, total bilirubin > 1.5 to 5 x ULN), median CrCl 82.7 mL/min (range 29.6-150; 15.3% with moderate renal impairment, CrCl 30 to < 60 mL/min), median albumin 3.9 g/dL (range 2.1-4.7), median BSA 1.9 m^2 (range 1.3-2.4).",
    dose_range      = "nab-paclitaxel (nanoparticle albumin-bound paclitaxel; Abraxane) 80-375 mg/m^2 IV as a 30-minute infusion on day 1, 8, 15 of each 28-day cycle or day 1 of each 21-day cycle (q3w). One 135 mg/m^2 q3w cohort received a 3-hour infusion; all other doses were 30-minute infusions.",
    regions         = "Pooled across the eight contributing studies; specific country breakdown not detailed in the main paper text on disk.",
    notes           = "Pharmacokinetic analysis used 1418 paclitaxel concentration records from 150 patients (Results, Population PK Model). The pharmacodynamic analysis subset comprised 558 ANC records from 125 patients in the first treatment cycle (Results, Population PD Model). The PK analysis pooled whole-blood (WB) and plasma concentration data with an estimated WB / plasma ratio of 1.00 (95% CI 0.93-1.10); the model below produces plasma concentrations directly and omits the WB / plasma bridging multiplier (see vignette Assumptions and deviations). The original residual error model used a two-subpopulation NONMEM MIXEST mixture with fractions 0.341 / 0.659 and proportional CV 38.0% / 17.9%; this implementation collapses the mixture to a single variance-weighted lognormal residual (see vignette Assumptions and deviations)."
  )

  ini({
    # ---- Population PK structural parameters (Chen 2014 Table 2 'PK parameter' block) ----
    # Three-compartment IV PK: central exchanges with peripheral1 via saturable Michaelis-Menten
    # distribution (VMTR, KMTR); central exchanges with peripheral2 via linear Q2; central also
    # has saturable Michaelis-Menten elimination (VMEL, KMEL). Volumes V1, V2, V3 are linear.
    # See Methods 'Population PK Model Development' and Figure 1.
    lvc   <- log(15.8);   label("V1: central volume of distribution (L)")                                # Table 2: V1 = 15.8 L (95% CI 13.71-17.85; RSE 6.20%)
    lvp   <- log(1650);   label("V2: first peripheral volume of distribution (L)")                       # Table 2: V2 = 1650 L (95% CI 1396-1935; RSE 7.9%)
    lvp2  <- log(75.4);   label("V3: second peripheral volume of distribution (L)")                      # Table 2: V3 = 75.4 L (95% CI 59.8-99.1; RSE 11.7%)
    lq2   <- log(41.6);   label("Q2: linear intercompartmental clearance central<->peripheral2 (L/h)")   # Table 2: Q2 = 41.6 L/h (95% CI 35.1-50.0; RSE 8.5%)
    lvmtr <- log(325000); label("VMTR: maximum saturable distribution rate central<->peripheral1 (ug/h)") # Table 2: VMTR = 325000 ug/h (95% CI 190694-540445; RSE 25.8%)
    lkmtr <- log(4260);   label("KMTR: concentration at half-maximal distribution rate (ug/L)")          # Table 2: KMTR = 4260 ug/L (95% CI 2210-7910; RSE 33.1%)
    lvmax <- log(8070);   label("VMEL: maximum elimination rate from central compartment (ug/h)")        # Table 2: VMEL = 8070 ug/h (95% CI 6500-9836; RSE 10.1%)
    lkmel <- log(40.2);   label("KMEL: concentration at half-maximal elimination rate (ug/L)")           # Table 2: KMEL = 40.2 ug/L (95% CI 24.9-58.9; RSE 20.0%)

    # ---- Covariate effects on PK ----
    # Albumin on VMEL: power model centred at the PK-population median albumin 3.9 g/dL
    # (Methods 'Population PK Model Development and Covariate Analysis': continuous covariates
    # are centered at their median and entered as power models). Form: VMEL_i = VMEL * (ALB/3.9)^e
    e_alb_vmax <- 0.554; label("Power exponent of (ALB / 3.9 g/dL) on VMEL (unitless)") # Table 2: 'Albumin on VMEL' = 0.554 (95% CI 0.15-0.94; RSE 37.5%)

    # ---- PK inter-individual variability (Chen 2014 Table 2 'Interindividual variability' block) ----
    # Paper reports CV%; converted to log-normal eta variance via omega^2 = log(1 + CV^2).
    # No correlation structure is reported in Table 2 -- treat as diagonal.
    #   V1   46.7% -> log(1 + 0.467^2) = 0.1958
    #   V2   25.5% -> log(1 + 0.255^2) = 0.0628
    #   VMTR 26.7% -> log(1 + 0.267^2) = 0.0688
    #   VMEL 22.3% -> log(1 + 0.223^2) = 0.0486
    # The 16.7% IIV on the WB / plasma ratio is not encoded here -- see population$notes.
    etalvc   ~ 0.1958
    etalvp   ~ 0.0628
    etalvmtr ~ 0.0688
    etalvmax ~ 0.0486

    # ---- PK residual error (Chen 2014 Table 2 'Residual variability' block) ----
    # Source model used a two-subpopulation mixture (NONMEM MIXEST) with proportional CVs of
    # 38.0% in subpopulation 1 (fraction 0.341) and 17.9% in subpopulation 2 (fraction 0.659).
    # Paclitaxel concentrations were Ln-transformed (Methods); a proportional error on log-
    # transformed concentration is the LTBS / lognormal pattern in nlmixr2. We collapse the
    # mixture to a single lognormal residual using a variance-weighted SD:
    #   var = 0.341 * 0.380^2 + 0.659 * 0.179^2 = 0.0703 -> SD ~= 0.265
    # See vignette Assumptions and deviations.
    expSd <- 0.265; label("Lognormal residual SD on plasma paclitaxel concentration (mixture-collapsed)") # Table 2: subpop1 = 38.0% CV (95% CI 31.8-46.9; RSE 17.1%); subpop2 = 17.9% CV (95% CI 15.1-19.7; RSE 17.8%)

    # ---- Population PD structural parameters (Chen 2014 Table 2 'PD parameter' block) ----
    # Friberg 2002 myelosuppression structure: proliferating pool (Prol), three transit
    # compartments (M1, M2, M3), and circulating ANC (Circ). KTR = (n+1)/MTT with n = 3 transits;
    # KCIRC = KTR (parameter-reduction assumption). Drug effect is linear in Cc; feedback term
    # is (BASE / Circ)^gamma. See Methods 'Population PD Model Development' and Figure 1.
    lmtt    <- log(117);     label("MTT: mean transit time Prol->Circ (h)")                                # Table 2: MTT = 117 h (95% CI 109-125; RSE 3.4%)
    lslope  <- log(0.00253); label("Slope: linear drug-effect slope on proliferation (1/(ng/mL))")          # Table 2: Slope = 0.00253 ng/mL^-1 (95% CI 0.00216-0.00290; RSE 7.47%)
    lrbase   <- log(4.28);    label("BASE: baseline circulating ANC (10^9 cells/L)")                        # Table 2: Baseline ANC = 4.28 x 10^9/L (95% CI 3.94-4.62; RSE 4.09%)
    lgamma  <- log(0.187);   label("gamma: feedback exponent on (BASE / Circ) (unitless)")                 # Table 2: Feedback parameter = 0.187 (95% CI 0.171-0.203; RSE 4.44%)

    # ---- Covariate effects on PD ----
    # Age on Slope: dichotomous AGE >= 65 indicator, multiplicative form
    #   Slope_i = Slope * (1 + e_age_slope * (AGE >= 65))
    # Results text: 'The typical value for paclitaxel effect slope was 0.00253 for patients
    # aged < 65 years and 0.0038 for patients aged >= 65 years' (Slope * (1 + 0.501) = 0.00380).
    e_age_slope <- 0.501; label("Multiplicative effect of (AGE >= 65) on Slope (unitless)")             # Table 2: 'Age on drug slope' = 0.501 (95% CI 0.172-0.830; RSE 33.5%)

    # Albumin on baseline ANC: power model centred at the PD-population median 4.0 g/dL
    #   BASE_i = BASE * (ALB / 4.0)^e
    # Negative exponent: higher albumin -> lower baseline ANC.
    e_alb_rbase <- -0.998; label("Power exponent of (ALB / 4.0 g/dL) on baseline ANC (unitless)") # Table 2: 'Albumin on baseline ANC' = -0.998 (95% CI -1.5 to -0.494; RSE 25.8%)

    # ---- PD inter-individual variability (Chen 2014 Table 2 'Interindividual variability' PD block) ----
    # Paper reports CV%; converted to log-normal eta variance via omega^2 = log(1 + CV^2).
    #   MTT          19.0% -> log(1 + 0.190^2) = 0.0353
    #   Slope        42.9% -> log(1 + 0.429^2) = 0.1697
    #   Baseline ANC 35.1% -> log(1 + 0.351^2) = 0.1162
    etalmtt   ~ 0.0353
    etalslope ~ 0.1697
    etalrbase  ~ 0.1162

    # ---- PD residual error (Chen 2014 Table 2 'Residual variability' PD row) ----
    # ANC values were Ln-transformed; proportional error on log-transformed ANC = lognormal
    # residual in nlmixr2. CV% maps to log-scale SD as expSd ~= CV/100 in the small-CV regime.
    expSd_ANC <- 0.291; label("Lognormal residual SD on absolute neutrophil count")                # Table 2: PD residual variability = 29.1% CV (95% CI 24.0-33.3; RSE 16.1%)
  })

  model({
    # SI -> US-convention unit conversion (canonical ALB is in SI g/L per the
    # 2026-06-19 register standardization audit; the original calibration
    # used the g/dL reference value, so convert inline here).
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL (factor 0.1)

    # ---- Derived dichotomous covariate ----
    # Chen 2014 dichotomised AGE at 65 years for the Slope covariate; the canonical AGE column
    # is continuous and the binary indicator is derived inline (see Robbie_2012_palivizumab.R
    # for the same pattern with ADA_TITER).
    age_gte65 <- (AGE >= 65)

    # ---- Individual PK parameters ----
    vc   <- exp(lvc   + etalvc)
    vp   <- exp(lvp   + etalvp)
    vp2  <- exp(lvp2)
    q2   <- exp(lq2)
    vmtr <- exp(lvmtr + etalvmtr)
    kmtr <- exp(lkmtr)
    vmax <- exp(lvmax + etalvmax) * (alb_gdL / 3.9)^e_alb_vmax
    kmel <- exp(lkmel)

    # Concentrations driving the saturable processes (ug/L = ng/mL)
    Cc  <- central     / vc
    Cp1 <- peripheral1 / vp
    Cp2 <- peripheral2 / vp2

    # Saturable Michaelis-Menten flux from central to peripheral1 (mass conservation: the
    # symmetric reverse flux is VMTR * Cp1 / (KMTR + Cp1)).
    rate_mtr_out <- vmtr * Cc  / (kmtr + Cc)
    rate_mtr_in  <- vmtr * Cp1 / (kmtr + Cp1)

    # Saturable Michaelis-Menten elimination from central
    rate_elim <- vmax * Cc / (kmel + Cc)

    # Linear intercompartmental flux central <-> peripheral2 (Q2 in L/h)
    rate_q2_out <- q2 * Cc
    rate_q2_in  <- q2 * Cp2

    # ---- PK ODE system ----
    d/dt(central)     <- -rate_elim - rate_mtr_out + rate_mtr_in - rate_q2_out + rate_q2_in
    d/dt(peripheral1) <-  rate_mtr_out - rate_mtr_in
    d/dt(peripheral2) <-  rate_q2_out  - rate_q2_in

    # ---- Individual PD parameters ----
    mtt   <- exp(lmtt   + etalmtt)
    slope <- exp(lslope + etalslope) * (1 + e_age_slope * age_gte65)
    rbase  <- exp(lrbase  + etalrbase)  * (alb_gdL / 4.0)^e_alb_rbase
    gamma <- exp(lgamma)

    # Transit-rate constant for the Friberg chain: KTR = (n + 1) / MTT with n = 3
    ktr <- 4 / mtt

    # Linear drug effect on proliferation (Methods 'Population PD Model Development': drug
    # effect Slope is 'a linear proportionality constant relating paclitaxel concentration in
    # the central compartment to its effect on bone marrow stem or progenitor cells').
    edrug <- slope * Cc

    # Feedback from circulating cells (Friberg 2002 form)
    feedback <- (rbase / circ)^gamma

    # ---- PD ODE system (Friberg 2002 myelosuppression chain) ----
    # precursor1 = Prol (proliferating pool); precursor2..4 = M1..M3 (maturation chain).
    d/dt(precursor1) <- ktr * precursor1 * (1 - edrug) * feedback - ktr * precursor1
    d/dt(precursor2) <- ktr * precursor1 - ktr * precursor2
    d/dt(precursor3) <- ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr * precursor3 - ktr * precursor4
    d/dt(circ)       <- ktr * precursor4 - ktr * circ

    # Initial conditions: system at steady state at the per-subject baseline ANC.
    precursor1(0) <- rbase
    precursor2(0) <- rbase
    precursor3(0) <- rbase
    precursor4(0) <- rbase
    circ(0)       <- rbase

    # ---- Observations and error model ----
    Cc  ~ lnorm(expSd)
    ANC <- circ
    ANC ~ lnorm(expSd_ANC)
  })
}

Jonsson_2016_tanezumab <- function() {
  description <- "Two-compartment population PK model for intravenous tanezumab, an anti-nerve-growth-factor IgG2 monoclonal antibody, in 1608 adults with moderate to severe osteoarthritis of the knee or hip pooled across four phase 3 trials (Jonsson 2016). Elimination from the central compartment is the sum of a linear clearance and a parallel Michaelis-Menten pathway (Vmax 8.03 ug/day, Km 27.7 ng/mL) attributed to target-mediated disposition; the saturable route supplies only 18%, 10% and 5% of total clearance at the 2.5, 5 and 10 mg dose levels. Clearance and both volumes scale with body weight as power functions centred at 84.7 kg (exponents 0.77, 0.554 and 0.302). Clearance additionally carries a Cockcroft-Gault creatinine-clearance power effect centred at 93.5 mL/min, a +14.3% male effect, and a +6.69% effect for the 2.5 and 5 mg dose groups relative to 10 mg; central volume carries a +17.5% male effect. Inter-individual variability is log-normal on CL, Vc, Vp and Vmax with a correlated CL-Vc block. Residual error is a two-class per-subject mixture on the log scale: 76.3% of subjects take the 13% component and the remainder the 54% component, selected by the MIX_LARGE_PROPRUV indicator."
  reference   <- "Jonsson EN, Xie R, Marshall SF, Arends RH. Population pharmacokinetics of tanezumab in phase 3 clinical trials for osteoarthritis pain. Br J Clin Pharmacol. 2016;81(4):688-699. doi:10.1111/bcp.12850"
  vignette    <- "Jonsson_2016_tanezumab"

  units       <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight (time-fixed). Selected over body surface area, body mass index and baseline lean body weight as the body-size descriptor because it gave the largest OFV drop on CL and because the fits on CL, V1 and V2 differed little between descriptors (Results, body-size paragraph). Enters CL, V1 and V2 as separate estimated power functions centred at the 84.7 kg model reference (Equations 3-5), which is the rounded cohort median and is slightly below the 86.6 kg cohort mean of Table 1.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance, Cockcroft-Gault computed on TOTAL body weight and NOT BSA-normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Raw un-normalized Cockcroft-Gault creatinine clearance in mL/min, the same size normalisation used by Delattre_2010_amikacin.R and Chen_2023_nemonoxacin.R rather than the BSA-normalized mL/min/1.73 m^2 default of this canonical. The source Methods truncate the column at a MAXIMUM of 150 mL/min before it enters the model, because Cockcroft-Gault on total body weight returns unreasonably high values in heavy subjects (Methods, PK analysis paragraph, citing reference 21); the cohort range in Table 1 runs to 301 mL/min, so the cap is active in practice. That cap is reproduced inside model() as min(CRCL, 150) so a user cannot silently extrapolate past it. Enters CL as the power function (CRCL_capped / 93.5)^0.108 (Equation 3); 93.5 mL/min is the model reference and is slightly below the 97.7 mL/min cohort mean of Table 1.",
      source_name        = "CLcr"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) -- NOTE this is the SOURCE PAPER's reference, the inverse of the canonical",
      notes              = "The source paper's reference subject is FEMALE: Equations 3 and 4 both read 'CLGENDER = 1 if female or 1 + theta if male', and the Table 2 footnote states the CL and V1 estimates are for a female. To keep the canonical 1 = female orientation while preserving the verbatim published coefficients, the effects are applied as (1 + e_sexf_cl * (1 - SEXF)) and (1 + e_sexf_vc * (1 - SEXF)) -- the same construction used by Bajaj_2017_nivolumab.R and Wada_2023_sparsentan.R. Men have 14.3% higher CL and 17.5% higher V1 than women (Table 3, rows 'Gender on CL' and 'Gender on V1').",
      source_name        = "GENDER"
    ),
    DOSE_HIGH = list(
      description        = "Highest-dose-cohort indicator: 1 = the 10 mg tanezumab arm, 0 = the 2.5 mg or 5 mg arm",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (the 10 mg arm) -- NOTE this is the SOURCE PAPER's reference, the inverse of the canonical",
      notes              = "Time-fixed per subject: the four phase 3 trials are parallel-group designs in which each patient stayed on one dose level for the whole study. Equation 3 reads 'CLDOSE = 1 if dose = 10 mg or 1 + theta12 if dose = 2.5 or 5 mg', so the source reference cohort is the HIGHEST dose, the inverse of this canonical's reference category. The effect is therefore applied as (1 + e_dose_high_cl * (1 - DOSE_HIGH)) to preserve the verbatim +0.0669 coefficient. Clearance is 6.69% higher in the 2.5 and 5 mg arms than in the 10 mg arm (Table 3, row 'Dose on CL'), and the Discussion is explicit that this step is over and above the concentration-dependent saturable pathway already in the structural model. The 2.5 and 5 mg arms share one level; the paper reports no separate estimate for them.",
      source_name        = "DOSE"
    ),
    MIX_LARGE_PROPRUV = list(
      description        = "Latent mixture-model class indicator for the log-scale residual error magnitude: 1 = subject assigned to the minority large-residual subpopulation (54% CV), 0 = subject assigned to the majority small-residual subpopulation (13% CV)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (small-residual subpopulation; 76.3% of the source cohort)",
      notes              = "Not a measured patient covariate -- this is the per-subject latent class index of the NONMEM mixture the source fitted on the residual error alone (Equation 2: 'Y = Yhat + eps1 if subpopulation 1 or eps2 if subpopulation 2'). Structural PK, covariate effects and IIV are shared by the two classes; only the residual magnitude switches. Table 2 estimates the mixture probability for the low-residual class at 0.763 (95% CI 0.738, 0.789). For typical-value simulation set MIX_LARGE_PROPRUV = 0; for population simulation draw MIX_LARGE_PROPRUV ~ Bernoulli(1 - 0.763) per subject. The eta shrinkages the paper reports separately per class (CL 11%/10%, V1 15%/23%, VM 66%/71%, V2 57%/79% for low/high) confirm the mixture is assigned at the subject level, not per observation.",
      source_name        = "$MIX class assignment"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL, V1, V2, VM and KM in the stepwise covariate model but not retained in the final model (Results, covariate-model paragraph). Cohort mean 61.4 years (SD 10.4), range 21-93 (Table 1). No point estimate is reported, so the effect cannot be implemented."
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Race was screened on CL, V1, V2, VM and KM but not retained in the final model. The cohort is 86.4% White, 11.2% Black, 0.9% Asian and 1.5% Other (Table 1). No point estimate is reported."
    ),
    OA_HIP = list(
      description = "Index osteoarthritis joint indicator (1 = hip, 0 = knee)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Site of OA was screened on all structural model parameters and not retained (Results, covariate-model paragraph). 69.0% knee / 31.0% hip (Table 1). Missing values were imputed as knee because they occurred only in study A4091011, which enrolled knee OA only. No point estimate is reported. OA_HIP is a descriptive placeholder name for a screened-but-dropped covariate and is deliberately not registered as a canonical column."
    ),
    ADA_POS = list(
      description = "Anti-drug-antibody-positive indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not tested as a covariate at all. Only 8 of 1601 patients were ADA-positive and the paper states their PK, pain response and safety profile did not differ from ADA-negative patients, so ADA status was excluded from covariate model building (Methods, PK analysis paragraph)."
    )
  )

  compartmentData <- list(
    central     = list(analyte = "tanezumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tanezumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1608,
    n_studies      = 4,
    age_range      = "21-93 years",
    age_median     = "mean 61.4 years (SD 10.4)",
    weight_range   = "34-170 kg",
    weight_median  = "mean 86.6 kg (SD 17.8); model reference 84.7 kg",
    sex_female_pct = 60.5,
    race_ethnicity = c(White = 86.4, Black = 11.2, Asian = 0.9, Other = 1.5),
    disease_state  = "moderate to severe osteoarthritis of the knee (69.0%) or hip (31.0%)",
    dose_range     = "2.5, 5 or 10 mg intravenously over a 5 min infusion every 8 weeks, for a total of two doses (studies A4091015, A4091018) or three doses (studies A4091011, A4091014)",
    regions        = "not reported by region; four multicentre phase 3 trials",
    notes          = "Baseline demographics from Table 1. Four randomized, double-blind, placebo-controlled, multicentre, parallel-group phase 3 trials: NCT00733902 (A4091011), NCT00744471 (A4091014), NCT00830063 (A4091015) and NCT00863304 (A4091018). Arm sizes 289 / 655 / 664 at 2.5 / 5 / 10 mg. The final analysis data set held 7592 plasma concentrations from 1608 patients after three pre-specified data-cleaning rules (n = 4, 188 and 190 removals) and removal of observations with |CWRES| > 5; re-including the cleaning-rule-3 patients and the outliers changed most fixed-effect estimates by less than 10%, with a maximum 23% change in inter-compartmental clearance. Samples were drawn pre-dose, 1 h post-dose and at weeks 4, 8 (pre and post), 16 (pre and post) and 24, plus week 32 in A4091011 and A4091014. Validated ELISA with an LLOQ of 12.0 ng/mL. NONMEM 7.1, FOCE(I), ADVAN6; covariate selection by the PsN stepwise covariate model procedure."
  )

  ini({
    # ---- Structural PK ----
    # Typical values are for the model reference subject: a FEMALE weighing 84.7 kg
    # with a creatinine clearance of 93.5 mL/min receiving the 10 mg dose
    # (Table 2 footnote; Equation 3 definition of theta1).
    lcl   <- log(0.135);   label("Linear elimination clearance, CL (L/day)")                 # Table 2 row 'CL (l day-1)' = 0.135 (95% CI 0.129, 0.14)
    lvc   <- log(2.71);    label("Central volume, V1 (L)")                                   # Table 2 row 'V1 (l)' = 2.71 (95% CI 2.66, 2.76)
    lq    <- log(0.371);   label("Inter-compartmental clearance, Q (L/day)")                 # Table 2 row 'Q (l day-1)' = 0.371 (95% CI 0.198, 0.545); no covariates
    lvp   <- log(1.98);    label("Peripheral volume, V2 (L)")                                # Table 2 row 'V2 (l)' = 1.98 (95% CI 1.72, 2.24)
    # The saturable pathway is reported in ug/day and ng/mL; both are converted to the
    # mg / L / day system of this file. VM/KM = 0.00803/0.0277 = 0.29 L/day reproduces the
    # low-concentration non-linear clearance quoted in the Discussion.
    lvmax <- log(0.00803); label("Maximum non-linear elimination capacity, VM (mg/day)")     # Table 2 row 'VM (ug day-1)' = 8.03 (95% CI 5.72, 10.3); 8.03 ug/day = 0.00803 mg/day
    lkm   <- log(0.0277);  label("Concentration at half-maximal non-linear elimination, KM (mg/L)") # Table 2 row 'KM (ng ml-1)' = 27.7 (95% CI 7.8, 47.7); 27.7 ng/mL = 0.0277 mg/L

    # ---- Covariate effects ----
    # Body size: separate estimated power exponents on WT/84.7 for CL, V1 and V2
    # (Equations 3, 4 and 5; theta8, theta9 and theta10).
    e_wt_cl   <- 0.77;    label("Power exponent on CL for WT/84.7 (unitless)")                # Table 2 row 'WT on CL' = 0.77 (95% CI 0.682, 0.858); Table 3: 10% change in WT gives 8% change in CL
    e_wt_vc   <- 0.554;   label("Power exponent on V1 for WT/84.7 (unitless)")                # Table 2 row 'WT on V1' = 0.554 (95% CI 0.489, 0.62); Table 3: 10% change in WT gives 5% change in V1
    e_wt_vp   <- 0.302;   label("Power exponent on V2 for WT/84.7 (unitless)")                # Table 2 row 'WT on V2' = 0.302 (95% CI 0.15, 0.454); Table 3: 10% change in WT gives 3% change in V2
    # Renal function: power exponent on the 150 mL/min-capped Cockcroft-Gault CLcr,
    # centred at 93.5 mL/min (Equation 3, theta11).
    e_crcl_cl <- 0.108;   label("Power exponent on CL for capped CRCL/93.5 (unitless)")       # Table 2 row 'CL cr on CL' = 0.108 (95% CI 0.0738, 0.141); Table 3: 10% change in CLcr gives 1% change in CL
    # Categorical effects are printed as fractional differences from the reference
    # category: ParCov = 1 for the common category, 1 + theta otherwise (Equation 1).
    # Both reference categories here are the INVERSE of the canonical column's reference,
    # so each effect is applied against (1 - covariate) in model().
    e_dose_high_cl <- 0.0669; label("Fractional increase in CL for the 2.5 and 5 mg arms relative to the 10 mg arm (unitless)") # Table 2 row 'Dose on CL' = 0.0669 (95% CI 0.0346, 0.0992); Table 3: CL is 7% higher with 2.5 and 5 mg
    e_sexf_cl      <- 0.143;  label("Fractional increase in CL for males relative to females (unitless)")                       # Table 2 row 'Gender on CL' = 0.143 (95% CI 0.106, 0.181); Table 3: males have 14% higher CL
    e_sexf_vc      <- 0.175;  label("Fractional increase in V1 for males relative to females (unitless)")                       # Table 2 row 'Gender on V1' = 0.175 (95% CI 0.143, 0.208); Table 3: males have 18% higher V1

    # ---- Inter-individual variability ----
    # Methods, PK analysis paragraph, and the Table 2 footnote both define the reported
    # %CV as the SQUARE ROOT OF THE VARIANCE estimated by NONMEM, so omega = %CV / 100 and
    # the variance below is (%CV / 100)^2. The 'Cov CL-V1' row is reported directly on the
    # raw covariance scale, not as a %CV, which is what fixes the scale of the whole block:
    # 0.034 / (0.26 * 0.20) = 0.65, close to the 0.67 correlation quoted in the Results
    # (the residual gap is rounding of 26 and 20 to two significant figures).
    etalcl + etalvc ~ c(0.0676,
                        0.034, 0.04)  # Table 2 rows 'IIV CL, %CV' = 26 (95% CI 25, 27), 'Cov CL-V1' = 0.034 (95% CI 0.03, 0.038) and 'IIV V1, %CV' = 20 (95% CI 19, 21)
    etalvp   ~ 0.04    # Table 2 row 'IIV V2, %CV' = 20 (95% CI 15, 24)
    etalvmax ~ 0.1681  # Table 2 row 'IIV VM, %CV' = 41 (95% CI 26, 52)

    # ---- Residual error ----
    # The base model used log transformation of both sides, i.e. an additive residual on the
    # log scale, and the Methods state the residual %CV is reported 'under the log
    # transformation'. The reported percentages are therefore log-scale SDs and map directly
    # onto lnorm(). The two components are a per-subject NONMEM mixture (Equation 2), gated
    # in model() by MIX_LARGE_PROPRUV; the same pattern as Kappelhoff_2005_ritonavir.R.
    expSd_p1 <- 0.13; label("Log-normal residual SD, small-residual subpopulation P1 (unitless)")  # Table 2 row 'Low RSV, %CV' = 13 (95% CI 13, 13); mixture probability 0.763
    expSd_p2 <- 0.54; label("Log-normal residual SD, large-residual subpopulation P2 (unitless)")  # Table 2 row 'High RSV, %CV' = 54 (95% CI 52, 55); mixture probability 1 - 0.763
  })

  model({
    # ---- Derived covariate terms ----
    # The source truncates Cockcroft-Gault CLcr at 150 mL/min BEFORE it enters the model
    # (Methods, PK analysis paragraph); the cap is reproduced here so the covariate cannot
    # be extrapolated past the range the effect was estimated over.
    crclCapped <- min(CRCL, 150)

    # ---- Individual PK parameters (reference: female, 84.7 kg, CLcr 93.5 mL/min, 10 mg) ----
    # Equation 3: TVCL = theta1 * CL_WT * CL_CLCR * CL_DOSE * CL_GENDER.
    cl   <- exp(lcl + etalcl) *
            (WT / 84.7)^e_wt_cl *
            (crclCapped / 93.5)^e_crcl_cl *
            (1 + e_dose_high_cl * (1 - DOSE_HIGH)) *
            (1 + e_sexf_cl * (1 - SEXF))
    # Equation 4: TVV1 = theta2 * V1_WT * V1_GENDER.
    vc   <- exp(lvc + etalvc) *
            (WT / 84.7)^e_wt_vc *
            (1 + e_sexf_vc * (1 - SEXF))
    # Equation 5: TVV2 = theta4 * V2_WT. Q carries no covariate.
    vp   <- exp(lvp + etalvp) * (WT / 84.7)^e_wt_vp
    q    <- exp(lq)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)

    # Plasma tanezumab concentration (mg/L = ug/mL).
    Cc <- central / vc

    # ---- ODE system: two compartments, parallel linear and Michaelis-Menten elimination ----
    # Doses are given as 5 min intravenous infusions directly into the central compartment.
    # The saturable term vmax * Cc / (km + Cc) is mg/day because vmax is mg/day and km and Cc
    # share units; at the trough concentrations of an 8-week interval it is the arm that makes
    # the non-linear contribution to total clearance dose-dependent.
    d/dt(central)     <- -(cl / vc) * central -
                          (q  / vc) * central +
                          (q  / vp) * peripheral1 -
                          vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # ---- Observation and residual error ----
    # Mixture-gated log-normal residual (Equation 2). Only the magnitude switches between the
    # two subpopulations; the structural and IIV models are shared.
    expSdMix <- expSd_p1 * (1 - MIX_LARGE_PROPRUV) + expSd_p2 * MIX_LARGE_PROPRUV
    Cc ~ lnorm(expSdMix)
  })
}

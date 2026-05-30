Lodise_2018_iclaprim <- function() {
  description <- "Two-compartment IV-infusion population PK model for iclaprim, a bacterial dihydrofolate reductase inhibitor, in adult patients with complicated skin and skin-structure infections from the pooled ASSIST-1 and ASSIST-2 phase 3 trials (Lodise 2018). Structural typical-value equations are additive-linear (NONMEM theta-sum form rather than power form): central volume V1 carries a body-weight slope; clearance CL carries age + sex (male shift) + sampling-occasion (day 1-2 vs day 4 +/- 1) shifts; peripheral volume V2 has no covariates; inter-compartmental clearance Q carries a severe-cSSSI-infection shift. Block-correlated log-normal IIV on V1, CL, V2 was retained in the source paper but only diagonal CV% values are tabulated -- off-diagonal covariances are not reported and are implemented here as diagonal-only (documented in the vignette Assumptions and deviations section). Combined proportional + additive residual error."
  reference <- paste(
    "Lodise TP, Bosso J, Kelly C, Williams PJ, Lane JR, Huang DB.",
    "Pharmacokinetic and pharmacodynamic analyses to determine the optimal",
    "fixed dosing regimen of iclaprim for treatment of patients with serious",
    "infections caused by Gram-positive pathogens.",
    "Antimicrob Agents Chemother. 2018;62(4):e01184-17.",
    "doi:10.1128/AAC.01184-17.",
    sep = " "
  )
  vignette <- "Lodise_2018_iclaprim"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Additive-linear effect on CL (Lodise 2018 Results final-model equation",
        "CL = theta2 + theta5 * AGE + ...): theta5 = -0.210 L/h per year. Older",
        "patients have lower CL. Cohort range 18-87 years; median 47.5 years",
        "(Table 3 baseline demographics). No reference centering -- the slope is",
        "applied directly to AGE so the NONMEM equation reproduces verbatim."
      ),
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Total body weight at study entry",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Additive-linear effect on V1 (Lodise 2018 Results final-model equation",
        "V1 = theta1 + theta6 * WT): theta6 = 0.353 L per kg. Cohort range 42-143",
        "kg; median 76 kg (Table 3 baseline demographics). No effect on CL --",
        "the paper explicitly notes the lack of association between weight and",
        "clearance was the motivation for the subsequent fixed-dose regimen",
        "selection. No reference centering -- the slope is applied directly to",
        "WT so the NONMEM equation reproduces verbatim."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Lodise 2018 encodes sex as a male-indicator (1 = male, 0 = female) on",
        "CL: theta7 = +2.78 L/h additive shift in males. To store under the",
        "canonical SEXF (1 = female, 0 = male) while preserving Lodise's literal",
        "+2.78 male-indicator coefficient, the effect is applied in model() as",
        "e_sex_cl * (1 - SEXF) -- a male (SEXF = 0) adds 2.78 L/h to CL and a",
        "female (SEXF = 1) does not. Follows the Bajaj_2017_nivolumab convention."
      ),
      source_name        = "SEX (1 = male, 0 = female)"
    ),
    OCC = list(
      description        = "Sampling-occasion indicator (1 = day 1 or 2 of treatment; 2 = day 4 +/- 1)",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Lodise 2018 encodes a binary occasion indicator (0 = day 1 or 2 sample,",
        "1 = day 4 +/- 1 sample) as a fixed-effect categorical on CL: theta9 =",
        "+3.97 L/h additive shift for day-4 samples. Interoccasion variability",
        "(IOV) was tested and removed during model building (Lodise 2018",
        "Results, model-building paragraph), so the surviving 'occasion' effect",
        "is a typical-value time-varying covariate, not an IOV random effect.",
        "Stored under the canonical integer-valued OCC column with",
        "OCC = 1 for day 1-2 samples and OCC = 2 for day 4 +/- 1 samples;",
        "decomposed inside model() to a binary day-4 indicator",
        "occ_day4 <- (OCC == 2). Per-occasion application: hold OCC = 1 across",
        "doses 1-3 (first ~36 h) and OCC = 2 from the day-3-4 boundary onward."
      ),
      source_name        = "occasion (0 / 1)"
    ),
    DIS_INFECT_CSSSI_SEV = list(
      description        = "Complicated skin and skin-structure infection severity indicator (1 = severe; 0 = not severe)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not severe cSSSI)",
      notes              = paste(
        "Lodise 2018 encodes a clinical severity-of-infection indicator (SOI;",
        "0 = not severe, 1 = severe) on inter-compartmental clearance Q only:",
        "theta8 = +13.5 L/h additive shift in severe-cSSSI patients (i.e., Q",
        "rises from a typical 1.85 L/h to 15.35 L/h, an eight-fold increase).",
        "The paper does not detail the exact clinical criteria that classified a",
        "patient as severe vs not severe in the ASSIST trials; the criteria are",
        "protocol-defined within the trial dataset. Time-fixed per subject in",
        "the ASSIST analysis."
      ),
      source_name        = "SOI (0 / 1)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 470L,
    n_studies        = 2L,
    n_concentrations = 3061L,
    age_range        = "18-87 years",
    age_median       = "47.5 years",
    weight_range     = "42-143 kg",
    weight_median    = "76 kg",
    bmi_range        = "16.8-54.7 kg/m^2",
    bmi_median       = "26.0 kg/m^2",
    sex_female_pct   = NA_real_,
    race_ethnicity   = "Not tabulated by category in Lodise 2018; ethnicity tested as a covariate and found not significant",
    disease_state    = "Adult patients with complicated skin and skin-structure infections (cSSSI). Patients with BMI > 40 and estimated creatinine clearance < 30 mL/min were excluded from the ASSIST studies. cSSSI severity (severe vs not severe) was retained as a covariate on inter-compartmental clearance Q in the final model.",
    dose_range       = "Iclaprim 0.8 mg/kg infused intravenously over 0.5 h every 12 h for 8 to 14 days (per ASSIST-1 / ASSIST-2 protocol).",
    regions          = "ASSIST-1 and ASSIST-2 were multicenter phase 3 trials; geographic distribution not tabulated in Lodise 2018",
    renal_function   = "Cockcroft-Gault creatinine clearance 18-172 mL/min (mean 99 +/- 29, median 102; Table 3). Patients with CLcr < 30 mL/min were excluded.",
    hepatic_function = "ALT 0.04-245 SI units (mean 13.4 +/- 21.6, median 1.21); total bilirubin 1.7-40.9 SI units (mean 7.5 +/- 5.4, median 6.0). ALT and total bilirubin were tested as covariates and found not significant.",
    notes            = paste(
      "Pooled population PK analysis of two phase 3 cSSSI trials (ASSIST-1,",
      "ASSIST-2; iclaprim 0.8 mg/kg infused IV over 0.5 h every 12 h for 8-14",
      "days). Of 492 patients who received iclaprim, PK data were available",
      "for 476 patients (241 patients from ASSIST-1 and 235 patients from",
      "ASSIST-2; 1,528 + 1,533 = 3,061 concentrations). 470 patients had day-4",
      "samples and could contribute to the Monte Carlo simulation analysis;",
      "6 patients had no day-4 concentrations. Sampling per protocol on two",
      "occasions: the first occasion after the first dose (day 1) and the",
      "second occasion on the 4th day (+/- 1 day) of treatment. Baseline",
      "demographics per Lodise 2018 Table 3; final population PK parameter",
      "estimates per Lodise 2018 Table 1 and the Results section's structural",
      "equations. NONMEM VI, FOCE+I, $OMEGA BLOCK for V1 / CL / V2 (Q has no",
      "random effect)."
    )
  )

  ini({
    # Structural-parameter intercepts (Lodise 2018 Table 1 final-model estimates;
    # additive-linear NONMEM theta-sum equations reproduced from Results).
    # Log-transformed in ini() to keep the canonical lcl / lvc / lvp / lq names
    # and ensure exp(.) is positive; the additive covariate slopes are applied
    # in model() so the typical-value equations match the source verbatim.
    lcl <- log(36.1);  label("CL intercept at AGE = 0, female (canonical reference), day 1-2 sample, non-severe cSSSI (L/h)")  # Lodise 2018 Table 1: theta2 = 36.1 (SE 2.04; 95% CI 32.02-40.18)
    lvc <- log(48.2);  label("V1 intercept at WT = 0 (L)")                                                                     # Lodise 2018 Table 1: theta1 = 48.2 (SE 4.22; 95% CI 39.76-56.64)
    lvp <- log(36.8);  label("Peripheral volume of distribution V2 (L)")                                                       # Lodise 2018 Table 1: theta3 = 36.8 (SE 2.62; 95% CI 31.56-42.04)
    lq  <- log(1.85);  label("Q intercept at non-severe cSSSI (L/h)")                                                          # Lodise 2018 Table 1: theta4 = 1.85 (SE 0.851; 95% CI 0.148-3.552)

    # Additive-linear covariate slopes (Lodise 2018 Table 1 + Results equations).
    # Effects are applied in model() as (intercept + slope * covariate), not as
    # multiplicative (1 + slope * covariate) factors, matching the published
    # NONMEM theta-sum form verbatim.
    e_age_cl  <- -0.210; label("Additive slope of age on CL (L/h per year)")                                # Lodise 2018 Table 1: theta5 = -0.210 (SE 0.031; 95% CI -0.272 to -0.148)
    e_wt_vc   <-  0.353; label("Additive slope of body weight on V1 (L per kg)")                            # Lodise 2018 Table 1: theta6 =  0.353 (SE 0.046; 95% CI  0.261 to  0.445)
    e_sex_cl  <-  2.78;  label("Additive shift in CL for male sex (L/h; applied as e_sex_cl * (1 - SEXF))") # Lodise 2018 Table 1: theta7 =  2.78  (SE 1.06;  95% CI  0.66  to  4.90)
    e_occ_cl  <-  3.97;  label("Additive shift in CL for day-4 sampling occasion (L/h)")                    # Lodise 2018 Table 1: theta9 =  3.97  (SE 0.694; 95% CI  2.58  to  5.36)
    e_infect_csssi_sev_q <- 13.5; label("Additive shift in Q for severe cSSSI (L/h)")                       # Lodise 2018 Table 1: theta8 = 13.5  (SE 2.29;  95% CI  8.92  to 18.08)

    # Inter-individual variability (log-normal exponential IIV on CL, V1, V2;
    # no IIV on Q in the final model). Lodise 2018 Table 1 reports CV% with 95%
    # CI; converted to variance via omega^2 = log(CV^2 + 1) (exact log-normal
    # form). The source paper retained the $OMEGA BLOCK structure (correlated
    # etas across V1, CL, V2) but does NOT tabulate off-diagonal covariances;
    # the diagonal-only implementation here is documented in the vignette
    # Assumptions and deviations section.
    etalcl ~ 0.15543  # log(0.41^2 + 1) = 0.15543; CV% 41% on CL (Lodise 2018 Table 1; 95% CI 38-44)
    etalvc ~ 0.23107  # log(0.51^2 + 1) = 0.23107; CV% 51% on V1 (Lodise 2018 Table 1; 95% CI 45-56)
    etalvp ~ 0.55362  # log(0.86^2 + 1) = 0.55362; CV% 86% on V2 (Lodise 2018 Table 1; 95% CI 65-104)

    # Residual error (combined proportional + additive; Lodise 2018 Table 1).
    # Concentrations measured in ng/mL, so the additive component is in ng/mL.
    propSd <- 0.365; label("Proportional residual error (fraction)")    # Lodise 2018 Table 1: residual error CV% = 36.5 (SE 0.0935; 95% CI 36.31-36.69)
    addSd  <- 5.19;  label("Additive residual error (ng/mL)")            # Lodise 2018 Table 1: additive residual error = 5.19 (SE 0.87; 95% CI 3.45-6.93)
  })

  model({
    # ----- Source-paper covariate decompositions -----
    # Lodise 2018 encodes sex as a male-indicator (1 = male, 0 = female). The
    # canonical SEXF stores 1 = female / 0 = male, so the paper's male effect
    # is recovered by computing (1 - SEXF) and applying e_sex_cl to it. This
    # preserves the Lodise 2018 theta7 = +2.78 literal coefficient and the
    # female-as-reference baseline (CL intercept theta2 = 36.1 is the female
    # value, matching the published equation).
    sex_male <- 1 - SEXF

    # Lodise 2018 encodes a binary sampling-occasion indicator (0 = day 1 or 2
    # sample, 1 = day 4 +/- 1 sample). The canonical OCC takes integer values
    # 1..N; the data assembler sets OCC = 1 for day 1-2 records and OCC = 2 for
    # day-4 records. The binary day-4 indicator is decomposed here.
    occ_day4 <- (OCC == 2)

    # ----- Typical-value parameters (Lodise 2018 Results final-model equations,
    # additive-linear NONMEM theta-sum form). At canonical reference covariates
    # (AGE = 0, female / SEXF = 1, day-1-2 / OCC = 1, non-severe / DIS = 0,
    # WT = 0), cl_typ = exp(lcl) = 36.1 and vc_typ = exp(lvc) = 48.2 and
    # q_typ = exp(lq) = 1.85, matching Lodise 2018 Table 1's theta1, theta2,
    # theta3, theta4 verbatim.
    cl_typ <- exp(lcl) + e_age_cl * AGE + e_sex_cl * sex_male + e_occ_cl * occ_day4
    vc_typ <- exp(lvc) + e_wt_vc * WT
    vp_typ <- exp(lvp)
    q_typ  <- exp(lq)  + e_infect_csssi_sev_q * DIS_INFECT_CSSSI_SEV

    # ----- Individual PK parameters -----
    # Multiplicative log-normal IIV applied to typical values (CL, V1, V2);
    # Q has no IIV in the final model.
    cl <- cl_typ * exp(etalcl)
    vc <- vc_typ * exp(etalvc)
    vp <- vp_typ * exp(etalvp)
    q  <- q_typ

    # ----- Micro-constants (ADVAN3 TRANS4 equivalent) -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- Two-compartment IV PK -----
    # Iclaprim is administered as an IV infusion (0.5 h in ASSIST trials, 2 h in
    # the proposed REVIVE-1 / -2 regimen) into the central compartment via the
    # data set cmt = central with rate / dur encoded on the dose record.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- Observation -----
    # Dose units mg, vc units L -> central / vc has units mg/L (= ug/mL).
    # Lodise 2018 reports plasma iclaprim concentrations in ng/mL, so multiply
    # by 1000 to convert mg/L to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}

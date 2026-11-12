Zhu_2026_ciprofol <- function() {
  description <- "Population PK/PD model for ciprofol (HSK3486) and the bispectral index (BIS) in elderly (>= 65 years) Chinese surgical patients receiving a 0.4 mg/kg intravenous loading dose followed by a continuous infusion (Zhu 2026; 20 patients, 386 arterial plasma samples with time-matched BIS). PK is a three-compartment model with first-order elimination from the central compartment; clearance carries power covariate effects of body weight (exponent 0.74, centred on 66 kg) and age (exponent -0.21, centred on 72.5 years). PD is an effect-compartment (ke0) link to an inhibitory sigmoid Imax model on BIS with no covariates. Log-normal inter-individual variability on all six PK disposition parameters and on IC50 and the Hill coefficient; proportional residual error on concentration and additive residual error on BIS."
  reference <- paste(
    "Zhu J, He J, Zhong B, Cao Y, Zhang X, Xu B. From pharmacokinetics to",
    "precision dosing: optimizing",
    "continuous infusion regimens of ciprofol for elderly patients.",
    "Front Pharmacol. 2026;17:1764590. doi:10.3389/fphar.2026.1764590.",
    "PMID: 41704278. PMCID: PMC12907304.",
    "Chinese Clinical Trial Registry ChiCTR2400093796.",
    sep = " "
  )
  vignette <- "Zhu_2026_ciprofol"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds. Verified against Zhu 2026 Figure 1 (schematic
  # of the PK/PD model: central V1, shallow peripheral V2, deep peripheral V3,
  # and an effect compartment linked to V1 by Ke0) and the Methods ('Arterial
  # blood sample collection and analysis': ciprofol was measured in plasma
  # from arterial blood).
  compartmentData <- list(
    central = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "ciprofol", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on clearance only, CL = 1.01 * (WT/66)^0.74 * (AGE/72.5)^-0.21 (Zhu 2026 Results, 'Population PK model', final clearance equation). The reference 66 kg is the cohort median weight (Table 1: 65.63 +/- 9.34, median 66, range 49-80 kg).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age at surgery.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on clearance only (same equation as WT). The reference 72.5 years is the cohort median age (Table 1: 72.95 +/- 4.47, median 72.5, range 67-82 years). The cohort was restricted to patients aged >= 65 years.",
      source_name = "Age"
    )
  )

  # Covariates the paper screened in the forward-inclusion / backward-
  # elimination search but did not retain in either the PK or the PD model
  # (Methods, 'Construction of the PK model'; Results, 'Population PD model':
  # 'common physiological covariates such as age, body weight, height, and
  # gender did not significantly improve model fitting'). Documentation only
  # -- these are not referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on PK and PD, not retained. 11 male / 9 female (Table 1)."
    ),
    HT = list(
      description = "Body height.",
      units = "cm",
      type = "continuous",
      notes = "Screened on PD, not retained. Cohort 151-172 cm, median 164 (Table 1)."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on PK, not retained. Cohort 21.20-29.38 kg/m^2, median 24.91 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened on PK, not retained. Cohort 16-28 U/L, median 21.5 (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened on PK, not retained. Cohort 8-73 U/L, median 21.5 (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on PK, not retained. Cohort 5-28.3 umol/L, median 9.35 (Table 1)."
    ),
    TPRO = list(
      description = "Total serum protein.",
      units = "g/L",
      type = "continuous",
      notes = "Screened on PK, not retained. Cohort 61.6-81.6 g/L, median 69 (Table 1)."
    ),
    CRCL = list(
      description = "Renal function (estimated glomerular filtration rate).",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "The Discussion names 'creatinine clearance' among the PK covariates without a significant effect; Table 1 reports it as GFR, 56.4-98.3 mL/min/1.73 m^2 (median 80.15). Serum creatinine 37-107 umol/L (median 59)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "67-82 years",
    age_median = "72.5 years",
    weight_range = "49-80 kg",
    weight_median = "66 kg",
    sex_female_pct = 45,
    race_ethnicity = "Chinese (single-centre study in China; race not tabulated)",
    disease_state = "Elderly patients (>= 65 years) with ASA physical status II (75%) or III (25%) undergoing elective orthopaedic (40%), general (30%), urological (15%) or thoracic (15%) surgery under general anaesthesia. Exclusions included severe cardiac, pulmonary, hepatic or renal dysfunction, coagulation disorders, dementia or other psychiatric illness, chronic psychotropic or analgesic use, and expected blood loss > 800 mL.",
    dose_range = "Ciprofol 0.4 mg/kg intravenous bolus over 1 min for induction, then a continuous infusion of 0.8 mg/kg/h starting 5 min later until skin suturing began.",
    regions = "China (General Hospital of Southern Theater Command of PLA), December 2024 - May 2025.",
    co_medication = "Remifentanil by effect-site target-controlled infusion (Minto model) in all patients; cisatracurium for intubation; midazolam (3 patients) when BIS exceeded 60; phenylephrine or ephedrine as needed (Supplementary Table S3).",
    notes = "386 arterial plasma samples, all above the 5 ng/mL LLOQ except pre-dose, with 386 time-matched BIS readings. Sampling at 0, 1, 3 and 5 min after the 1-min bolus; 5, 10, 15, 30, 45 and 60 min after starting the infusion; and 3, 5, 10, 30 min and 1, 1.5, 2, 4 and 6 h after stopping it. Ciprofol quantified by UPLC-MS/MS (linear 5-5000 ng/mL). Fitted in NONMEM 7.2 by FOCE-I; the PD parameters are reported as a separate final PD model (Table 3)."
  )

  ini({
    # Final PK model estimates (Zhu 2026 Table 2, 'Final model' column).
    # Time in minutes, dose in mg, volumes in L. The paper's three-compartment
    # parameterisation maps to the nlmixr2lib canonical names as:
    #   V1 (central)            -> vc
    #   V2 (shallow peripheral) -> vp
    #   V3 (deep peripheral)    -> vp2
    #   Q2 (central <-> V2)     -> q
    #   Q3 (central <-> V3)     -> q2
    lcl  <- log(1.01);  label("Clearance at 66 kg and 72.5 years, CL (L/min)")          # Zhu 2026 Table 2: CL = 1.01 L/min (RSE 7.23%)
    lvc  <- log(2.95);  label("Central volume of distribution, V1 (L)")                 # Zhu 2026 Table 2: V1 = 2.95 L (RSE 8.54%)
    lvp  <- log(45.15); label("Shallow peripheral volume of distribution, V2 (L)")      # Zhu 2026 Table 2: V2 = 45.15 L (RSE 6.20%)
    lvp2 <- log(76.79); label("Deep peripheral volume of distribution, V3 (L)")         # Zhu 2026 Table 2: V3 = 76.79 L (RSE 5.15%)
    lq   <- log(0.76);  label("Intercompartmental clearance to V2, Q2 (L/min)")         # Zhu 2026 Table 2: Q2 = 0.76 L/min (RSE 6.83%)
    lq2  <- log(0.66);  label("Intercompartmental clearance to V3, Q3 (L/min)")         # Zhu 2026 Table 2: Q3 = 0.66 L/min (RSE 3.84%)

    # Covariate exponents on clearance (Zhu 2026 Results, 'Population PK
    # model': CL = 1.01 * (WT/66)^0.74 * (Age/72.5)^-0.21).
    e_wt_cl  <- 0.74;  label("Power exponent of (WT/66 kg) on CL (unitless)")       # Zhu 2026 Table 2: 'WT on CL' = 0.74 (RSE 12.42%)
    e_age_cl <- -0.21; label("Power exponent of (AGE/72.5 years) on CL (unitless)") # Zhu 2026 Table 2: 'Age on CL' = -0.21 (RSE 8.63%)

    # Final PD model estimates (Zhu 2026 Table 3).
    lke0   <- log(1.09);   label("Effect-compartment equilibration rate constant, ke0 (1/min)")           # Zhu 2026 Table 3: Ke0 = 1.09 1/min (RSE 25.20%)
    lrbase <- log(93.40);  label("Baseline BIS before ciprofol, E0 (BIS units)")                           # Zhu 2026 Table 3: E0 = 93.40 (RSE 0.71%)
    limax  <- log(45.77);  label("Maximum reduction in BIS attributable to ciprofol, Imax (BIS units)")    # Zhu 2026 Table 3: Imax = 45.77 (RSE 3.83%)
    lec50  <- log(233.91); label("Effect-site concentration producing half of Imax, IC50 (ng/mL)")         # Zhu 2026 Table 3: IC50 = 233.91 ng/mL (RSE 8.06%)
    lhill  <- log(3.00);   label("Hill coefficient of the sigmoid Imax model, gamma (unitless)")           # Zhu 2026 Table 3: gamma = 3.00 (RSE 12.92%)

    # Inter-individual variability. Zhu 2026 Methods ('Construction of the
    # PK model') define theta_i = theta_pop * exp(eta_i) with eta ~ N(0,
    # omega^2) and 'omega is the standard deviation of the inter-individual
    # variability model', so the printed omega values are log-scale SDs and
    # the variances below are omega^2. See the vignette 'Assumptions and
    # deviations' for the variance reading the Results prose also suggests.
    etalvc  ~ 0.0036 # Zhu 2026 Table 2: omega(V1) = 0.06 (RSE 22.18%); 0.06^2
    etalvp  ~ 0.0121 # Zhu 2026 Table 2: omega(V2) = 0.11 (RSE 9.12%); 0.11^2
    etalvp2 ~ 0.0121 # Zhu 2026 Table 2: omega(V3) = 0.11 (RSE 10.36%); 0.11^2
    etalcl  ~ 0.0009 # Zhu 2026 Table 2: omega(CL) = 0.03 (RSE 18.58%); 0.03^2
    etalq   ~ 0.0009 # Zhu 2026 Table 2: omega(Q2) = 0.03 (RSE 21.64%); 0.03^2
    etalq2  ~ 0.0004 # Zhu 2026 Table 2: omega(Q3) = 0.02 (RSE 8.75%); 0.02^2

    # Zhu 2026 Table 3 prints eta-shrinkage for Ke0 (29.72%), E0 (23.51%) and
    # Imax (26.85%) as well as for IC50 and gamma, implying random effects on
    # those three parameters, but its BSV section lists omega only for IC50
    # and gamma. The unreported variances are carried as fixed zero.
    etalke0   ~ fixed(0)  # Zhu 2026 Table 3: eta-shrinkage 29.72% printed, omega(Ke0) not reported
    etalrbase ~ fixed(0)  # Zhu 2026 Table 3: eta-shrinkage 23.51% printed, omega(E0) not reported
    etalimax  ~ fixed(0)  # Zhu 2026 Table 3: eta-shrinkage 26.85% printed, omega(Imax) not reported
    etalec50  ~ 0.0081    # Zhu 2026 Table 3: omega(IC50) = 0.09 (RSE 35.99%); 0.09^2
    etalhill  ~ 0.0004    # Zhu 2026 Table 3: omega(gamma) = 0.02 (RSE 88.18%); 0.02^2

    # Residual error. PK: proportional, C_obs = C_pred * (1 + eps) with
    # eps ~ N(0, sigma^2) (Methods; selected over additive and combined,
    # Supplementary Table S2). PD: Table 3 prints sigma = 7.70 without naming
    # the error model; a value in BIS units can only be an additive SD.
    propSd    <- 0.26; label("Proportional residual error on ciprofol concentration (fraction)") # Zhu 2026 Table 2: sigma = 0.26 (RSE 4.41%)
    addSd_BIS <- 7.70; label("Additive residual error on BIS (BIS units)")                       # Zhu 2026 Table 3: sigma = 7.70 (RSE 4.81%)
  })

  model({
    # Covariate reference values: cohort medians printed in the clearance
    # equation (Zhu 2026 Results, 'Population PK model').
    ref_wt  <- 66
    ref_age <- 72.5

    # Amounts in mg and volumes in L give mg/L; the paper reports
    # concentrations and IC50 in ng/mL (1 mg/L = 1000 ng/mL).
    mgL_to_ngmL <- 1000

    cl  <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl * (AGE / ref_age)^e_age_cl
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp + etalvp)
    vp2 <- exp(lvp2 + etalvp2)
    q   <- exp(lq + etalq)
    q2  <- exp(lq2 + etalq2)

    ke0   <- exp(lke0 + etalke0)
    rbase <- exp(lrbase + etalrbase)
    imax  <- exp(limax + etalimax)
    ec50  <- exp(lec50 + etalec50)
    hill  <- exp(lhill + etalhill)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment disposition with intravenous input into central
    # (Zhu 2026 Figure 1).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    Cc <- central / vc * mgL_to_ngmL

    # Effect compartment, dCe/dt = ke0 * C1 - ke0 * Ce (Zhu 2026 Methods,
    # 'Construction of the population PD model'); the state is the
    # effect-site concentration in ng/mL.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Inhibitory sigmoid Imax model, E = E0 - Emax * Ce^gamma /
    # (EC50^gamma + Ce^gamma) (Zhu 2026 Methods; Table 3 reports the same
    # parameters as Imax and IC50).
    BIS <- rbase - imax * effect^hill / (ec50^hill + effect^hill)

    Cc  ~ prop(propSd)
    BIS ~ add(addSd_BIS)
  })
}

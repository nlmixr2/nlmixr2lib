Ma_2020_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in 56 adult Chinese kidney transplant recipients receiving prophylactic vancomycin in the first postoperative weeks (Ma 2020). Clearance scales by power functions of body weight (reference 59.95 kg) and raw estimated glomerular filtration rate (mL/min, reference 36.67); volume of distribution scales by a power function of body weight. IIV on CL only; proportional residual error. Estimated from 195 trough concentrations collected by routine therapeutic drug monitoring."
  reference <- "Ma K-f, Liu Y-x, Jiao Z, Lv J-h, Yang P, Wu J-y, Yang S. Population Pharmacokinetics of Vancomycin in Kidney Transplant Recipients: Model Building and Parameter Optimization. Front Pharmacol. 2020;11:563967. doi:10.3389/fphar.2020.563967"
  vignette <- "Ma_2020_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Ma 2020 Results final CL and V equations normalize by 59.95 kg, the cohort median in Table 1 (mean 58.27, SD 8.47, range 37.7-79 kg). Figure 2B shows WT varying over the postoperative period, so the column may be time-varying.",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate from serum creatinine (Chinese-modified MDRD-type equation, raw mL/min, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Ma 2020 Methods, 'Patients and Data Collection': GFR was estimated from serum creatinine with the",
        "published formula of Chu et al. 2020, reported in ml min-1 and including a 0.742 female factor and a",
        "1.233 'correction for Chinese'. The typeset formula has lost its exponent formatting (it reads",
        "'GFR = 2.104*sCr(uM) - 1.154*age - 1.154*(0.742 for female)*1.233'), so the exact coefficients cannot",
        "be recovered from this paper; supply the GFR estimate directly. The final CL equation normalizes by",
        "36.67 mL/min, which is NOT the Table 1 per-patient median (39.91 mL/min; mean 41.95, SD 25.46, range",
        "3.38-108.61) -- the paper does not say what 36.67 is. Figure 2D shows GFR rising over the",
        "posttransplantation period, so the column is time-varying (renal-graft function recovery). Stored",
        "under canonical CRCL on the raw mL/min scale, following the register's documented raw-scale",
        "precedents (e.g. Georges_2009_ceftazidime.R, Zhou_2019_vancomycin.R)."
      ),
      source_name = "GFR"
    )
  )

  # Screened in the Ma 2020 stepwise covariate search but not retained. Results,
  # 'Assessment of Covariates and Evaluation of Models': "categorical variables
  # such as combined medications including imipenem, and continuous variables
  # such as endogenous creatinine concentration and age, were excluded." No
  # point estimates are published for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Ma 2020 Table 1: mean 43.72 (SD 9.92), median 43.5, range 24-70 years. Screened; not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Ma 2020 Table 1: median 164 (IQR 141-283, range 60-1490) umol/L. Screened; not retained (its information enters through the eGFR)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "categorical",
      notes = "Ma 2020 Table 1: 35 male / 21 female. Recorded as a basal characteristic; not retained (enters only through the eGFR equation's female factor)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Ma 2020 Table 1: median 16 (IQR 13-19) IU/L. Recorded covariate; not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Ma 2020 Table 1: median 14.5 (IQR 10-24.25) IU/L. Recorded covariate; not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Ma 2020 Table 1 'Proteinaemia (PROT)': median 59.15 (IQR 55.42-63.92) g/L. Recorded covariate; not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 56L,
    n_studies = 1L,
    age_range = "24-70 years",
    age_median = "43.5 years (mean 43.72, SD 9.92)",
    weight_range = "37.7-79 kg",
    weight_median = "59.95 kg (mean 58.27, SD 8.47)",
    sex_female_pct = 37.5,
    race_ethnicity = "Chinese (single center, The First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou)",
    disease_state = "Adult kidney transplant recipients (all grafts from brain-dead donors) receiving IV vancomycin as postoperative prophylaxis; pretransplant dialysis median 48 months. Concomitant prednisone (55/56), mycophenolate mofetil (55/56), tacrolimus (52/56) or cyclosporine A (4/56).",
    dose_range = "500 mg per administration; a single 500 mg dose on the first postoperative day, then 1-4 administrations per day (500-2000 mg/day) adjusted by trough TDM. Infusion duration not reported.",
    regions = "China (Hangzhou)",
    renal_function = "eGFR (raw mL/min) median 39.91 (IQR 32.53-59.44, range 3.38-108.61); serum creatinine median 164 umol/L",
    n_concentrations = 195L,
    notes = "Retrospective single-center TDM study of transplant surgeries March-June 2017. All 195 observations are troughs drawn about 30 min before the morning dose; first monitoring median postoperative day 4. HPLC-UV assay (calibration range 1.5625-100 ug/mL). NONMEM 7.4, FOCE-I. Proportional residual error was selected over additive and combined. Evaluated with goodness-of-fit plots and NPDE (Figure 5; mean 0.093, variance 0.990)."
  )

  ini({
    # Structural parameters (Ma 2020 Table 2 'Final model'). The reference
    # subject has WT = 59.95 kg and GFR = 36.67 mL/min (Results final CL and V
    # equations).
    lcl <- log(2.08)
    label("Clearance at WT = 59.95 kg and GFR = 36.67 mL/min (L/h)") # Ma 2020 Table 2 final model: CL = 2.08 L/h (CV 3.4%)
    lvc <- log(63.2)
    label("Volume of distribution at WT = 59.95 kg (L)") # Ma 2020 Table 2 final model: V = 63.2 L (CV 6.7%)

    # Covariate effects (Ma 2020 Results):
    #   CL = 2.08 * (WT/59.95)^1.07 * (GFR/36.67)^0.698
    #   V  = 63.2 * (WT/59.95)^0.934
    e_crcl_cl <- 0.698
    label("Power exponent on (CRCL/36.67) for CL") # Ma 2020 Table 2: theta_1 'Influential factor for GFR on CL' = 0.698 (CV 7.7%)
    e_wt_cl <- 1.07
    label("Power exponent on (WT/59.95) for CL") # Ma 2020 Table 2: theta_2 'Influential factor for WT on CL' = 1.07 (CV 20.3%)
    e_wt_vc <- 0.934
    label("Power exponent on (WT/59.95) for V") # Ma 2020 Table 2: theta_3 'Influential factor for WT on V' = 0.934 (CV 43.1%)

    # Inter-individual variability, exponential (Methods: P_ij = TV(P_j) *
    # exp(eta_ij)), on CL only -- Table 2 carries no omega for V. Table 2
    # prints omega_1 = 21.5% under the heading 'Intersubject variance of CL';
    # it is read as a CV% and converted to log(0.215^2 + 1) = 0.04517. The
    # variance reading (omega^2 = 0.215) is falsified by Table 3: under it
    # every one of the 15 recommended weight/GFR regimens gives only 80-88%
    # attainment of AUC0-24 >= 400, whereas the paper reports >90% for all.
    # The CV reading gives 97-99.5%. See the vignette.
    etalcl ~ 0.04517 # Ma 2020 Table 2 final model, omega_1 = 21.5% CV (RSE 33.7%) -> log(0.215^2 + 1)

    # Proportional residual error (Methods residual equation with the additive
    # term dropped; Results: the proportional model had the lowest OFV).
    propSd <- 0.242
    label("Proportional residual error (fraction)") # Ma 2020 Table 2 final model: sigma_1 = 24.2% (CV 20.2%)
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 59.95)^e_wt_cl * (CRCL / 36.67)^e_crcl_cl
    vc <- exp(lvc) * (WT / 59.95)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

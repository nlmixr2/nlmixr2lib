Brooks_2026_thiotepa <- function() {
  description <- "Parent-metabolite population PK model for intravenous thiotepa and its active metabolite TEPA in children and young adults undergoing allogeneic or autologous hematopoietic cell transplantation (Brooks 2026): two-compartment thiotepa with linear clearance plus saturable (Michaelis-Menten) conversion to a one-compartment TEPA metabolite (volume fixed to 1 L); fat-free mass (Al-Sallami) allometric scaling referenced to 15 kg and two sigmoid age-maturation functions, one on the clearance terms and one on the thiotepa central volume."
  reference <- paste(
    "Brooks J, Ismaeil S, Hughes MS, Apsel-Winger B, Lu Y, Huang L, Higham C,",
    "Nelson K, Long S, Gupta A, Lalefar NR, Anderson E, Pulsipher MA, Pai SY,",
    "Keizer RJ, Dvorak CC, Long-Boyle J. Prospective Evaluation of Thiotepa",
    "Pharmacokinetics in Children and Young Adults Undergoing Allogeneic or",
    "Autologous Hematopoietic Cell Transplantation. J Clin Pharmacol.",
    "2026;66(9):e70296. doi:10.1002/jcph.70296.",
    sep = " "
  )
  vignette <- "Brooks_2026_thiotepa"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Verified against Brooks 2026 Methods 'Bioanalysis of TT and TEPA' (both
  # analytes quantified in plasma by LC-MS/MS) and the Supporting Information
  # NONMEM model file ($MODEL COMP=(CENTR), (PERIPH), (METAB); IV infusion
  # into CENTR). TEPA amount is in mg of TEPA formed per mg of thiotepa
  # converted: the model file applies no molecular-weight conversion.
  compartmentData <- list(
    central = list(analyte = "thiotepa", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "thiotepa", units = "mg", specimen = "plasma", verified = TRUE),
    central_tepa = list(
      analyte = "TEPA (N,N',N''-triethylenethiophosphoramide)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the model only through fat-free mass (FFM), computed inside model() with the Al-Sallami paediatric equations exactly as in the Supporting Information NONMEM model file. FFM scales CL, Vmax, Q and CL_TEPA as (FFM/15)^0.75 and Vc, Vp as (FFM/15)^1; the 15 kg reference and both exponents are fixed (Results 'Final Population PK Model'; Table 2 units '/15 kg').",
      source_name = "WT"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to compute BMI = WT / (HT/100)^2 inside the Al-Sallami FFM equation (Supporting Information model file 'BMI = WT/((HT*HT)/10000)').",
      source_name = "HT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives two sigmoid maturation functions AGE^hill / (AGE^hill + TM50^hill): one shared by thiotepa CL, Vmax and TEPA CL (TM50 0.370 years, Hill 1.056) and one on the thiotepa central volume (TM50 0.3067 years, Hill 2.712) (Table 2; Supporting Information model file FMAT and FMAT2). Also enters the Al-Sallami FFM age correction. Must be > 0: at AGE = 0 the maturation fractions are 0. Observed range 0.23-22.5 years.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 (male)",
      notes = "Selects the sex-specific Al-Sallami FFM equation and has no other effect. The Supporting Information model file branches on a SEX column with 'IF(SEX.EQ.0)' selecting the female equation (1.11, 7.1, -1.1; 8780 + 244*BMI) and ELSE the male one (0.88, 13.4, -12.7; 6680 + 216*BMI); the paper never states its SEX coding, so the female/male assignment here follows the published Al-Sallami 2015 equations themselves.",
      source_name = "SEX"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance (bedside Schwartz under 17 years; Cockcroft-Gault with actual body weight from 17 years)",
      units = "mL/min/1.73 m^2 (Schwartz) or mL/min (Cockcroft-Gault)",
      type = "continuous",
      notes = "Screened in the stepwise covariate analysis (Methods 'Covariate Modeling') but not retained in the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Methods 'Covariate Modeling') but not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Methods 'Covariate Modeling') but not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Methods 'Covariate Modeling') but not retained; the Discussion notes it was a covariate on thiotepa clearance in the earlier Huitema adult analysis."
    ),
    BILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods 'Covariate Modeling') but not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Tested on TEPA clearance following the Huitema analysis but not retained (Discussion)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened as a size descriptor (Methods 'Covariate Modeling'); FFM was retained instead."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 91L,
    n_studies = 2L,
    n_observations = 445L,
    age_range = "0.23-22.5 years",
    age_median = "4.5 years",
    weight_range = "4.6-110 kg",
    weight_median = "16.7 kg",
    sex_female_pct = 44,
    race_ethnicity = c(White = 47, Multiple = 20, Asian = 12, Unknown = 12, Other = 5, Black = 2, NativeAmerican = 1),
    disease_state = "Children and young adults undergoing allogeneic (n = 48; 33 malignant, 15 non-malignant) or autologous (n = 43; all malignant) hematopoietic cell transplantation with thiotepa-containing conditioning",
    dose_range = "Allogeneic: 4, 5 or 7.5 mg/kg as a 1-h IV infusion (5 mg/kg every 12 h x 2 doses is the usual regimen); autologous: 10 mg/kg or 300 mg/m^2 as a 2-h IV infusion (300 mg/m^2 every 24 h x 3 doses)",
    regions = "United States (UCSF Benioff Children's Hospitals San Francisco and Oakland, University of Minnesota Masonic Children's Hospital, Rady Children's Hospital San Diego, and the multi-centre PTCTC CSIDE infant SCID study)",
    notes = "Brooks 2026 Table 1. 51 male / 40 female; Hispanic or Latino 28, not Hispanic or Latino 56, unknown 7. PK sampled after a single dose: intensive sampling (4-5 samples to 24 h) in the first 17 subjects, two samples (0-2 h and 2-24 h after the end of infusion) in the remaining 74. 210 thiotepa and 235 TEPA quantifiable concentrations; values below the 2.5 ng/mL LLOQ were used as reported. NONMEM 7.6 FOCE-I."
  )

  ini({
    # Thiotepa structural parameters, typical values for a 15 kg FFM subject
    # with mature clearance and volume. Values from the Supporting Information
    # NONMEM model file $THETA (the final run; they equal Table 2's 'Final
    # Model Value' column to the printed precision).
    lcl <- log(6.665); label("Thiotepa linear clearance CL_TT (L/h per 15 kg FFM)") # SI $THETA(1) 'TVCL' = 6.665; Table 2 'Typical value for thiotepa CL (L/h/15 kg)' 6.665 (RSE 9.9%)
    lvc <- log(15.837); label("Thiotepa central volume Vc (L per 15 kg FFM)") # SI $THETA(2) 'TVV' = 15.837; Table 2 'Volume of the central compartment (L/15 kg)' 15.837 (RSE 7.0%)
    lq <- log(1.876); label("Thiotepa intercompartmental clearance Q (L/h per 15 kg FFM)") # SI $THETA(3) 'TVQ' = 1.876; Table 2 'Intercompartmental CL (L/h/15 kg)' 1.876 (RSE 29.6%)
    lvp <- log(7.139); label("Thiotepa peripheral volume Vp (L per 15 kg FFM)") # SI $THETA(4) 'TVVP' = 7.139; Table 2 'Volume of the peripheral compartment' 7.139 (RSE 18.2%)
    lvmax <- log(3.3164); label("Maximum rate of saturable thiotepa-to-TEPA conversion Vmax (mg/h per 15 kg FFM)") # SI $THETA(6) 'TVVMAX' = 3.3164; Table 2 'Vmax, mg/h/15 kg' 3.316 (RSE 12.9%)
    lkm <- log(0.5056); label("Thiotepa concentration at half-maximal conversion Km (mg/L)") # SI $THETA(7) 'TVKM' = 0.5056; Table 2 'Km, mg/L' 0.5056 (RSE 20.4%)

    # TEPA metabolite
    lcl_tepa <- log(2.034); label("TEPA clearance CL_TEPA (L/h per 15 kg FFM)") # SI $THETA(5) 'TVCLM' = 2.034; Table 2 'Typical value of TEPA CL (L/h/15 kg)' 2.035 (RSE 6.3%)
    lvc_tepa <- fixed(log(1)); label("TEPA volume Vm (L)") # SI $PK 'VM = 1'; Table 2 'Volume of the central compartment (L)' 1 Fixed

    # Fat-free-mass allometry, fixed (Results 'Final Population PK Model';
    # SI $PK 'wtcl = (FFM/15)**0.75', 'wtv = FFM/15')
    e_ffm_cl <- fixed(0.75); label("Allometric FFM exponent on CL_TT, Vmax, Q and CL_TEPA (unitless)") # SI $PK wtcl = (FFM/15)**0.75
    e_ffm_vc <- fixed(1); label("Allometric FFM exponent on Vc and Vp (unitless)") # SI $PK wtv = FFM/15

    # Age maturation
    lhill_cl <- log(1.056); label("Hill coefficient of the clearance maturation function (unitless)") # SI $THETA(10) 'HILL' = 1.056; Table 2 'Fmat, CL' 1.056 (RSE 23.3%)
    ltm50_cl <- log(0.370); label("Age at 50% clearance maturation TM50 (years)") # SI $THETA(11) 'TM50' = 0.370; Table 2 'TM50, CL, years' 0.370 (RSE 29.9%)
    lhill_vc <- log(2.712); label("Hill coefficient of the central-volume maturation function (unitless)") # SI $THETA(12) 'HILL2' = 2.712; Table 2 'Fmat,Vc' 2.712 (RSE 30.5%)
    ltm50_vc <- log(0.3067); label("Age at 50% central-volume maturation TM50 (years)") # SI $THETA(13) 'TM502' = 0.3067; Table 2 'TM50, Vc, years' 0.3067 (RSE 17.7%)

    # Inter-individual variability: the SI model file's $OMEGA (exponential
    # ETA(1) on CL, ETA(2) on VMAX; no Km eta). Table 2 prints different
    # values and a Km row; see the vignette errata.
    etalcl ~ 0.08612 # SI $OMEGA '1 CL' = 0.08612 (29.3% CV)
    etalvmax ~ 0.1755 # SI $OMEGA '2 VMAX' = 0.1755 (41.9% CV)

    # Residual error: proportional, $SIGMA 1 FIX with the SD carried as a THETA
    propSd <- 0.2842; label("Thiotepa proportional residual error (fraction)") # SI $THETA(8) 'PROP ERROR TT' = 0.2842; Table 2 prints 29.4%
    propSd_tepa <- 0.2248; label("TEPA proportional residual error (fraction)") # SI $THETA(9) 'PROP ERROR TEPA' = 0.2248; Table 2 22.5%
  })

  model({
    # Al-Sallami 2015 paediatric fat-free mass (SI $PK)
    bmi <- WT / (HT / 100)^2
    ffm_female <- (1.11 + (1 - 1.11) / (1 + (AGE / 7.1)^(-1.1))) * 9270 * WT / (8780 + 244 * bmi)
    ffm_male <- (0.88 + (1 - 0.88) / (1 + (AGE / 13.4)^(-12.7))) * 9270 * WT / (6680 + 216 * bmi)
    ffm <- SEXF * ffm_female + (1 - SEXF) * ffm_male

    size_cl <- (ffm / 15)^e_ffm_cl
    size_v <- (ffm / 15)^e_ffm_vc

    # Sigmoid age maturation (SI $PK FMAT, FMAT2)
    hill_cl <- exp(lhill_cl)
    tm50_cl <- exp(ltm50_cl)
    fmat_cl <- AGE^hill_cl / (AGE^hill_cl + tm50_cl^hill_cl)
    hill_vc <- exp(lhill_vc)
    tm50_vc <- exp(ltm50_vc)
    fmat_vc <- AGE^hill_vc / (AGE^hill_vc + tm50_vc^hill_vc)

    # Individual parameters
    cl <- exp(lcl + etalcl) * size_cl * fmat_cl
    vmax <- exp(lvmax + etalvmax) * size_cl * fmat_cl
    km <- exp(lkm)
    vc <- exp(lvc) * size_v * fmat_vc
    q <- exp(lq) * size_cl
    vp <- exp(lvp) * size_v
    cl_tepa <- exp(lcl_tepa) * size_cl * fmat_cl
    vc_tepa <- exp(lvc_tepa)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Thiotepa (TT) plasma concentration drives the saturable conversion
    ctt <- central / vc
    conversion <- vmax * ctt / (km + ctt)

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - conversion
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_tepa) <- conversion - cl_tepa / vc_tepa * central_tepa

    Cc <- central / vc
    Cc_tepa <- central_tepa / vc_tepa

    Cc ~ prop(propSd)
    Cc_tepa ~ prop(propSd_tepa)
  })
}

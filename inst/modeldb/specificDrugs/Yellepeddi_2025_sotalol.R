Yellepeddi_2025_sotalol <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order oral absorption",
    "for sotalol in 22 adults with atrial fibrillation or atrial flutter",
    "(AFIB/AFL) enrolled in the PK/PD substudy of the PEAKS Registry, who",
    "received an expedited 1 h intravenous loading dose followed by two",
    "oral maintenance doses 12 h apart. Creatinine clearance is a power",
    "covariate on clearance (CL * [CRCL/92.4]^0.65) and body weight is a",
    "power covariate on the central volume (Vc * [WT/104]^1.2). Relative",
    "oral bioavailability was estimated above unity (Foral = 1.6), which",
    "the authors flag as an unexplained finding requiring further study.",
    "Interindividual variability is carried on Ka, CL and Vc; residual",
    "error is proportional. Plasma concentrations are returned in ng/mL",
    "(the source control stream scales the central compartment as",
    "S2 = Vc/1000 with doses in mg and volumes in L).",
    "The file also packages the paper's three concentration-QTc linear",
    "regressions, which were fitted in R (not in NONMEM) against the",
    "individual model-predicted sotalol concentrations at the times of the",
    "Bazett-corrected QTc observations, and are returned as the derived",
    "outputs QTc (ms), dQTc (change from baseline, ms) and pctQTc",
    "(percent change from baseline). They carry no residual error term",
    "because the source reports only slope, intercept, R^2 and p-value;",
    "the paper applies them deterministically to simulated concentrations",
    "in its Monte Carlo dosing evaluation (Methods 2.7).",
    sep = " "
  )

  reference <- paste(
    "Yellepeddi VK, Ismail M, Bunch TJ, Deering TF, Holubkov R, Kennedy R,",
    "Mittal S, Perez M, Piccini JP, Pokharel P, Savona S, Verma N,",
    "Steinberg B, Watt K. (2025). Population Pharmacokinetics and",
    "Pharmacodynamics of Sotalol Following Expedited Intravenous Loading",
    "in Patients With Atrial Arrhythmias.",
    "CPT: Pharmacometrics & Systems Pharmacology 14(4):658-666.",
    "doi:10.1002/psp4.13302. PMCID PMC12001255.",
    "Final NONMEM control stream: Supplementary Datafile S1",
    "(file PSP4-14-658-s004.docx of the publisher supplementary bundle).",
    sep = " "
  )

  vignette <- "Yellepeddi_2025_sotalol"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Yellepeddi 2025 Supplementary
  # Datafile S1 ($SUBROUTINES ADVAN4 TRANS4; depot = oral dosing site with
  # F1, central = plasma-sampled compartment scaled S2 = V2/1000, peripheral
  # = the distribution compartment scaled S3 = V3/1000).
  compartmentData <- list(
    depot       = list(analyte = "sotalol", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "sotalol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sotalol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (a single admission weight was recorded).",
        "Enters the central volume as the power term (WT / 104)^e_wt_vc.",
        "The reference 104 kg is the cohort median body weight",
        "(Yellepeddi 2025 Table 1: mean 107 +/- 29.9 kg, median 104 kg,",
        "range 68.7-185 kg) and is written literally into the source",
        "control stream as WTCOVVOL = (WTKG/104)**THETA(7).",
        "Retained in the final model after stepwise forward addition",
        "(p < 0.001; delta OFV 6.8 for weight on Vc alone, 14.8 for the",
        "joint CrCl-on-CL + WT-on-Vc model; Table S2 models 10 and 11).",
        "Note the cohort is heavy: the median 104 kg is well above a",
        "typical 70 kg reference, so supplying weights from a leaner",
        "population extrapolates below the observed range."
      ),
      source_name        = "WTKG"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault formula",
        "and reported WITHOUT body-surface-area normalization, i.e. in raw",
        "mL/min, not mL/min/1.73 m^2."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject; serum creatinine was required within 24 h",
        "prior to dosing (Yellepeddi 2025 Methods 2.1).",
        "SIZE NORMALIZATION - this model uses the raw un-normalized",
        "Cockcroft-Gault clearance in mL/min, the same normalization used",
        "by Delattre_2010_amikacin.R and Georges_2009_ceftazidime.R, NOT",
        "the mL/min/1.73 m^2 default of the CRCL register entry. Supplying",
        "a BSA-normalized value would silently rescale the renal-function",
        "term. The source paper computed the BSA-normalized clearance as a",
        "separate covariate using DuBois and DuBois BSA and tested it too:",
        "'Adding BSA-normalized CrCl on CL did not show any influence on",
        "the model' (Results 3.2), so the BSA-normalized column is listed",
        "in covariatesDataExcluded and must not be substituted here.",
        "Enters clearance as the power term (CRCL / 92.4)^e_crcl_cl. The",
        "reference 92.4 mL/min is the cohort median (Table 1: mean",
        "116 +/- 57.4, median 92.4, range 64.3-306 mL/min) and is written",
        "literally into the source control stream as",
        "CRCLCOV = (CRCL/92.4)**THETA(8).",
        "Enrollment targeted normal renal function (CrCl > 60 mL/min;",
        "Methods 2.1), so the model is not informed below ~64 mL/min",
        "even though sotalol labelling reduces the dose in renal",
        "impairment."
      ),
      source_name        = "CRCL"
    )
  )

  # Covariates screened during stepwise covariate model development but NOT
  # retained in the final model (Yellepeddi 2025 Methods 2.5.2, Results 3.2,
  # Table S2, and the $INPUT record of Supplementary Datafile S1). Documented
  # here for provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at enrollment",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a power model on CL (Table S2 model 7, TVCL = theta1 * [Age/69]^theta2); delta OFV 0, not retained."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional model on CL (Table S2 model 2); delta OFV 0, not retained. Cohort was 19 male / 3 female."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Carried in the analysis dataset ($INPUT HTCM of Datafile S1) to compute body surface area; never tested as a covariate in its own right."
    ),
    BSA = list(
      description = "Body surface area (DuBois and DuBois)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Carried in the analysis dataset ($INPUT BSA) solely to normalize creatinine clearance; never tested as a covariate in its own right."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Carried in the analysis dataset ($INPUT SCR) as the input to the Cockcroft-Gault calculation; the renal covariate that was tested is the derived CRCL, not SCR itself."
    ),
    CRCL_BSA = list(
      description = "Body-surface-area-normalized creatinine clearance (1.73 x Cockcroft-Gault CrCl / BSA)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened as an alternative renal covariate on CL ($INPUT BSACRCL). 'Adding BSA-normalized CrCl on CL did not show any influence on the model' (Results 3.2); the raw mL/min CRCL was retained instead. Table 1 reports 87.9 +/- 32, median 75.3, range 59.6-188 mL/min/1.73 m^2."
    ),
    POT = list(
      description = "Serum potassium",
      units       = "mg/dL as reported by the source (Table 1); values 3.7-5 are on the mmol/L scale",
      type        = "continuous",
      notes       = "Carried in the analysis dataset ($INPUT Potassium); not reported as a tested covariate in Table S2. Table 1 reports 4.2 +/- 0.3, median 4.1, range 3.7-5."
    ),
    MAGNESIUM = list(
      description = "Serum magnesium",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Carried in the analysis dataset ($INPUT Magnesium); not reported as a tested covariate in Table S2. Table 1 reports 2 +/- 0.2, median 2, range 1.7-2.5."
    ),
    CONMED_NSAID = list(
      description = "Concomitant non-steroidal anti-inflammatory drug use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional model on CL (Table S2 model 3); delta OFV 0, not retained."
    ),
    CONMED_BETABLOCKER = list(
      description = "Concomitant beta-blocker use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional model on CL (Table S2 model 5); delta OFV 0, not retained. 19 of 22 patients (86.4%) were on a beta-blocker."
    ),
    CONMED_CCB = list(
      description = "Concomitant calcium channel blocker use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional model on CL (Table S2 model 6); delta OFV 0, not retained. 3 of 22 patients (13.6%)."
    ),
    CONMED_ANTIPLATELET = list(
      description = "Concomitant antiplatelet agent use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Listed among the screened concomitant medications (Methods 2.5.2) and carried in the dataset ($INPUT Antiplatelet); not retained. 2 of 22 patients (9%)."
    ),
    CONMED_ANTICOAGULANT = list(
      description = "Concomitant anticoagulant use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Listed among the screened concomitant medications (Methods 2.5.2) and carried in the dataset ($INPUT Anticoagulant); not retained."
    ),
    DIS_CHF = list(
      description = "History of heart failure indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional model on CL (Table S2 model 4); delta OFV 0, not retained. 7 of 22 patients (31.8%)."
    ),
    DIS_CAD = list(
      description = "History of coronary artery disease indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a proportional model on CL (Table S2 model 8); delta OFV 0, not retained. 12 of 22 patients (54.5%)."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Named in Methods 2.5.2 as a planned categorical covariate but not testable: 'Race was not included in the covariate analysis as all patients in the dataset were white' (20 of 22 white, 2 unknown)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 22L,
    n_studies        = 1L,
    n_observations   = 99L,
    age_range        = "48-79 years",
    age_median       = "69 years (mean 67.8 +/- 8.13)",
    weight_range     = "68.7-185 kg",
    weight_median    = "104 kg (mean 107 +/- 29.9)",
    sex_female_pct   = 13.6,
    race_ethnicity   = c(White = 91, Unknown = 9),
    disease_state    = paste(
      "Adults (18 years or older) admitted primarily for intravenous",
      "sotalol loading -- initiation or dose titration -- to treat atrial",
      "fibrillation or atrial flutter, with no other planned procedure",
      "other than cardioversion. Excluded: active ventricular-arrhythmia",
      "treatment, prior class III antiarrhythmic intolerance, bradycardia",
      "below 40 bpm, 2nd/3rd degree heart block without a pacemaker,",
      "baseline QTc at or above 450 ms (500 ms with bundle branch block),",
      "or severe left ventricular hypertrophy above 1.5 cm",
      "(Yellepeddi 2025 Methods 2.1).",
      "Baseline cardiac status: LVEF 53.9 +/- 8 % (median 55, range",
      "25-65), heart rate 71.8 +/- 16 bpm, 14 of 22 in sinus rhythm at",
      "dosing, 4 of 22 with baseline QRS above 120 ms, 7 with heart",
      "failure history, 12 with coronary artery disease history."
    ),
    renal_function   = paste(
      "Cockcroft-Gault creatinine clearance 116 +/- 57.4 mL/min (median",
      "92.4, range 64.3-306); BSA-normalized 87.9 +/- 32 mL/min/1.73 m^2",
      "(median 75.3, range 59.6-188); serum creatinine 0.96 +/- 0.16",
      "mg/dL. PK-substudy enrollment deliberately targeted normal renal",
      "function (CrCl above 60 mL/min), so no patient with clinically",
      "significant renal impairment contributed data."
    ),
    dose_range       = paste(
      "One intravenous loading dose infused over 1 h, then at least a 4 h",
      "delay to the first oral dose and a 12 h delay to the second oral",
      "dose. IV doses 45.6-150 mg (median 92.3); oral doses 80 mg (50% of",
      "oral administrations) or 120 mg (50%), median 100 mg (Table 1 and",
      "Table S1). All 22 patients received all three doses (one IV, two",
      "oral). The parent PEAKS registry evaluated IV doses of 60-125 mg",
      "followed by two oral doses of 40-160 mg."
    ),
    regions          = "United States; 10 academic, private and hybrid health systems.",
    notes            = paste(
      "PK/PD substudy of the prospective PEAKS Registry (Prospective",
      "Evaluation Analysis and Kinetics of IV Sotalol), enrolled",
      "4 February 2022 to 13 June 2023; 210 screened, 167 in the parent",
      "registry, 22 in the PK/PD substudy.",
      "99 plasma sotalol concentrations, all above the 5 ng/mL LLOQ",
      "(UPLC-MS/MS, NMS Labs; between-run precision 4.89%, recovery",
      "79.1%), sampled 0-30 min before the IV dose, 0-5 min after end of",
      "infusion, 3 h after end of infusion, 0-30 min before the first",
      "oral dose, 2-4 h after the first oral dose and 2-4 h after the",
      "second oral dose (Methods 2.2). Observed concentrations ranged",
      "220-1900 ng/mL.",
      "104 Bazett-corrected QTc values from the same 22 patients support",
      "the concentration-QTc regressions; QT correction was applied only",
      "when heart rate exceeded 60 bpm (about 80% of the ECGs). Baseline",
      "QTc median 435 ms (range 386-482); post-dose median 459 ms (range",
      "388-548).",
      "Estimation: NONMEM 7.5 FOCE with interaction, via PsN 5.3.0 and",
      "Finch Studio; pcVPC (1000 replicates) and a 1000-sample bootstrap",
      "(96.7% success) were used for evaluation.",
      "Sponsor: AltaThera Inc. sponsored the PEAKS Registry."
    )
  )

  ini({
    # ==================================================================
    # PARAMETER SOURCING NOTE (read before changing any value below).
    #
    # Three source artifacts report the final fixed effects and they do
    # not fully agree:
    #   (T2)  Table 2 "Final population PK model parameter estimates",
    #         printed to one decimal place.
    #   (T3)  Table 3, whose first column repeats the final estimate
    #         alongside the bootstrap summary.
    #   (S1)  Supplementary Datafile S1, "The final NONMEM model code",
    #         whose $THETA / $OMEGA / $SIGMA records carry values to two
    #         significant figures.
    #
    # S1's records are the final estimates rounded to two significant
    # figures, not hand-chosen starting values. The falsifier is the
    # variance block, which S1 could not have matched by chance:
    #   $SIGMA 0.023          vs T2 RUV sigma^2 = 0.0231
    #   $OMEGA 0.068 (CL)     vs log(1 + 0.266^2) = 0.06837
    #   $OMEGA 0.21  (Vc)     vs log(1 + 0.481^2) = 0.20811
    #   $OMEGA 0.23  (Ka)     vs log(1 + 0.505^2) = 0.22715
    # All four agree to two significant figures under the log-normal CV
    # convention; none agrees under the alternative reading in which
    # T2's IIV percentages are log-scale standard deviations (which
    # would give 0.0708 / 0.2314 / 0.2550). That simultaneously fixes
    # the IIV scale used below AND establishes that S1's $THETA record
    # holds the final estimates.
    #
    # Value selected per parameter = the most precise rendering that the
    # majority of artifacts support:
    #   CL   9.5   T2 = T3 = S1
    #   Q    71.3  T2 = T3 (S1's 71 is the 2-s.f. rounding)
    #   Vc   38.2  T2 = T3 (S1's 38)
    #   Vp   89.1  T2 = T3 (S1's 89)
    #   Ka   0.16  S1 only at 2 s.f.; T2 and T3 print 0.2 at 1 d.p.
    #   Foral 1.6  T2 = S1;  T3 and the Discussion print 1.5  <-- CONFLICT
    #   CrCl exponent 0.65  S1; T2 and T3 print 0.7 at 1 d.p.
    #   WT   exponent 1.2   T3 = S1;  T2 prints 1.1          <-- CONFLICT
    # Each of T2 and T3 carries exactly one value that the other two
    # artifacts contradict. Both conflicts are recorded in the vignette
    # Errata section.
    # ==================================================================

    lka <- log(0.16)
    label("Sotalol first-order oral absorption rate constant (Ka, 1/h)")
    # Yellepeddi 2025 Datafile S1 $THETA (0, 0.16) ; KA. Table 2 prints
    # 0.2 (%RSE 34.2) and Table 3 prints 0.2 with bootstrap mean 0.2 and
    # 95% CI 0.08-0.4; both are the one-decimal rendering of 0.16.

    lcl <- log(9.5)
    label("Sotalol clearance at the reference creatinine clearance of 92.4 mL/min (CL, L/h)")
    # Yellepeddi 2025 Table 2 Clearance (CL) = 9.5 L/h (%RSE 33.1);
    # Table 3 final estimate 9.5 (bootstrap mean 8.2, 95% CI 3-13.2);
    # Datafile S1 $THETA (0, 9.5) ; CL. All three agree.

    lvc <- log(38.2)
    label("Sotalol central volume of distribution at the reference body weight of 104 kg (Vc, L)")
    # Yellepeddi 2025 Table 2 Central Volume of distribution (Vc) = 38.2 L
    # (%RSE 31.2); Table 3 final estimate 38.2 (bootstrap mean 48.2,
    # 95% CI 25.6-70.6); Datafile S1 $THETA (0, 38) ; V2.

    lq <- log(71.3)
    label("Sotalol intercompartmental clearance (Q, L/h)")
    # Yellepeddi 2025 Table 2 Intercompartmental clearance (Q) = 71.3 L/h
    # (%RSE 20.3); Table 3 final estimate 71.3 (bootstrap mean 53.8,
    # 95% CI 22.6-83.7); Datafile S1 $THETA (0, 71) ; Q.

    lvp <- log(89.1)
    label("Sotalol peripheral volume of distribution (Vp, L)")
    # Yellepeddi 2025 Table 2 Peripheral volume of distribution (Vp) =
    # 89.1 L (%RSE 25.4); Table 3 final estimate 89.1 (bootstrap mean
    # 102.4, 95% CI 54.7-130.2); Datafile S1 $THETA (0, 89) ; V3.

    lfdepot <- log(1.6)
    label("Oral bioavailability of sotalol relative to intravenous (Foral, unitless)")
    # Yellepeddi 2025 Table 2 'Oral bioavailability relative to IV
    # (Foral)' = 1.6 (%RSE 28.4) and Datafile S1 $THETA (0, 1.6) ; F1.
    # Table 3 and the Discussion instead report 1.5 (bootstrap mean 1.5,
    # 95% CI 0.9-2.2). Estimated, not fixed: F1 = TVF1*EXP(ETA(6)) in S1
    # with $OMEGA 0.0 FIX on ETA(6), i.e. the typical value is estimated
    # while its interindividual variability is fixed to zero. A value
    # above 1 is not physically interpretable as a fraction absorbed; the
    # authors state 'The reason for this higher oral bioavailability must
    # be evaluated further in future studies' (Discussion). It is kept as
    # published because it is load-bearing for the oral-dose predictions.

    e_crcl_cl <- 0.65
    label("Power exponent of creatinine clearance on CL, referenced to 92.4 mL/min (unitless)")
    # Yellepeddi 2025 Datafile S1 $THETA (0, 0.65) ; CRCL on CL, entering
    # as CRCLCOV = (CRCL/92.4)**THETA(8). Table 2 prints 0.7 (%RSE 31.8)
    # for 'theta CrCl on CL (CL*[CrCl/92.4]^theta)' and Table 3 prints
    # 0.7 (bootstrap mean 0.7, 95% CI 0.06-1.2); both are the
    # one-decimal rendering of 0.65. Retained in the covariate model with
    # delta OFV 9.7 alone (Table S2 model 9) and p < 0.001.

    e_wt_vc <- 1.2
    label("Power exponent of body weight on Vc, referenced to 104 kg (unitless)")
    # Yellepeddi 2025 Datafile S1 $THETA (0, 1.2) ; WT on V, entering as
    # WTCOVVOL = (WTKG/104)**THETA(7), and Table 3 'Power coefficient WT
    # on Vc' final estimate 1.2 (bootstrap mean 1, 95% CI 0.4-2).
    # Table 2 instead prints 1.1 (%RSE 49.6). Retained with delta OFV 6.8
    # alone (Table S2 model 10) and p < 0.001. Note this is close to the
    # linear (exponent 1) scaling the Discussion describes, not the
    # allometric 0.75/1 convention.

    # ==================================================================
    # Interindividual variability. Table 2 reports IIV as a percentage;
    # the sourcing note above establishes that these are log-normal
    # coefficients of variation, so omega^2 = log(1 + CV^2). NONMEM
    # placement (Datafile S1 $PK): ETA(1) on KA, ETA(2) on CL, ETA(3) on
    # V2; ETA(4) on Q, ETA(5) on V3 and ETA(6) on F1 were all fixed to
    # zero variance and so are absent here. The paper reports no
    # off-diagonal covariance, so the etas are independent.
    # ==================================================================

    etalka ~ 0.227149
    # Yellepeddi 2025 Table 2 'IIV - Ka (%)' = 50.5 (%RSE 64.6,
    # shrinkage 24.9%); omega^2 = log(1 + 0.505^2) = 0.227149, matching
    # Datafile S1 $OMEGA 0.23 ; KA.

    etalcl ~ 0.068366
    # Yellepeddi 2025 Table 2 'IIV - CL(%)' = 26.6 (%RSE 74.3,
    # shrinkage 24.4%); omega^2 = log(1 + 0.266^2) = 0.068366, matching
    # Datafile S1 $OMEGA 0.068 ; CL. Table S2 shows this fell from 57.6%
    # in the base model once CrCl was added to CL.

    etalvc ~ 0.208106
    # Yellepeddi 2025 Table 2 'IIV - Vc (%)' = 48.1 (%RSE 90.9,
    # shrinkage 24.9%); omega^2 = log(1 + 0.481^2) = 0.208106, matching
    # Datafile S1 $OMEGA 0.21 ; V2.

    propSd <- 0.151987
    label("Proportional residual error on plasma sotalol concentration (fraction)")
    # Yellepeddi 2025 Table 2 'RUV (sigma^2)' = 0.0231 (%RSE 19.8,
    # shrinkage 18.5%); propSd = sqrt(0.0231) = 0.151987. Datafile S1
    # $ERROR uses Y = F + F*EPS(1) with $SIGMA 0.023, i.e. a pure
    # proportional model on the linear concentration scale.

    # ==================================================================
    # Concentration-QTc relationships (Yellepeddi 2025 Methods 2.6 and
    # Results 3.3). These were fitted in R by ordinary linear regression
    # of the observed Bazett-corrected QTc endpoints on the individual
    # (empirical Bayes) model-predicted sotalol concentration, NOT
    # jointly in NONMEM, and no residual standard deviation is reported
    # for any of the three -- only the slope, the intercept, R^2 and the
    # p-value. They are therefore packaged as deterministic derived
    # outputs with no error model, exactly as the paper uses them in its
    # Monte Carlo dosing evaluation (Methods 2.7).
    # ==================================================================

    e0_qtc <- 440
    label("Intercept of the linear concentration-QTc regression, i.e. QTc at zero sotalol concentration (ms)")
    # Yellepeddi 2025 Results 3.3: 'For the QTc vs. sotalol plasma levels
    # model, the R^2 = 0.1, p = 0.0011 with slope (0.015) and intercept
    # (440 ms).' Consistent with the observed baseline QTc median of
    # 435 ms (Results 3.1).

    e_sotalol_qtc <- 0.015
    label("Slope of the linear concentration-QTc regression (ms per ng/mL)")
    # Yellepeddi 2025 Results 3.3, same sentence as e0_qtc. Equivalent to
    # 15 ms per ug/mL; for comparison the healthy-volunteer
    # concentration-QTcF slope in Darpo_2014_racSotalol_QTcF.R is 24 ms
    # per ug/mL.

    e0_dqtc <- 11
    label("Intercept of the linear concentration-versus-change-from-baseline-QTc regression (ms)")
    # Yellepeddi 2025 Results 3.3: 'For change in QTc from baseline vs.
    # sotalol plasma levels, model the R^2 = 0.29, p = 1.1e-07 with slope
    # (0.019) and intercept (11 ms).'

    e_sotalol_dqtc <- 0.019
    label("Slope of the linear concentration-versus-change-from-baseline-QTc regression (ms per ng/mL)")
    # Yellepeddi 2025 Results 3.3, same sentence as e0_dqtc. Equivalent
    # to 19 ms per ug/mL.

    e0_pctqtc <- 1.7
    label("Intercept of the linear concentration-versus-percent-change-from-baseline-QTc regression (percent)")
    # Yellepeddi 2025 Results 3.3: 'For % change in QTc vs. sotalol
    # plasma levels, model has R^2 = 0.27, p = 1.6e-08 with slope
    # (0.0045) and intercept (1.7%; Figure 4).'

    e_sotalol_pctqtc <- 0.0045
    label("Slope of the linear concentration-versus-percent-change-from-baseline-QTc regression (percent per ng/mL)")
    # Yellepeddi 2025 Results 3.3, same sentence as e0_pctqtc. Internally
    # consistent with the absolute change model: 0.019 ms per ng/mL over
    # the 435 ms baseline median is 0.00437 percent per ng/mL.
  })

  model({
    # ==================================================================
    # 1. Individual parameters. Covariate forms are transcribed from
    #    Yellepeddi 2025 Datafile S1 $PK:
    #      WTCOVVOL = (WTKG/104)**THETA(7);  V2 = TVV2*WTCOVVOL*EXP(ETA(3))
    #      CRCLCOV  = (CRCL/92.4)**THETA(8); CL = TVCL*CRCLCOV*EXP(ETA(2))
    #    Q, Vp and Foral carry no covariate.
    # ==================================================================
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (CRCL / 92.4)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 104)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # ==================================================================
    # 2. Micro-constants for the ADVAN4/TRANS4 two-compartment system.
    # ==================================================================
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ==================================================================
    # 3. ODE system. Dose intravenous doses into `central` (as a 1 h
    #    infusion) and oral doses into `depot`.
    # ==================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ==================================================================
    # 4. Relative oral bioavailability, applied to the depot only so the
    #    intravenous dose keeps F = 1 (Datafile S1 applies F1 to the
    #    ADVAN4 depot compartment).
    # ==================================================================
    f(depot) <- exp(lfdepot)

    # ==================================================================
    # 5. Observation. Datafile S1 scales the central compartment as
    #    S2 = V2/1000, so with doses in mg and volumes in L the predicted
    #    concentration is in ng/mL rather than mg/L.
    # ==================================================================
    Cc <- central / vc * 1000

    # ==================================================================
    # 6. Concentration-QTc outputs (Results 3.3). Deterministic
    #    functions of the predicted concentration; see the ini() note on
    #    why they carry no residual error term. QTc is the absolute
    #    Bazett-corrected interval, dQTc its change from the patient's
    #    pre-IV-dose baseline, and pctQTc that change as a percentage of
    #    the baseline. Each comes from its own regression, so the three
    #    are not algebraically consistent with one another.
    # ==================================================================
    QTc    <- e0_qtc    + e_sotalol_qtc    * Cc
    dQTc   <- e0_dqtc   + e_sotalol_dqtc   * Cc
    pctQTc <- e0_pctqtc + e_sotalol_pctqtc * Cc

    Cc ~ prop(propSd)
  })
}

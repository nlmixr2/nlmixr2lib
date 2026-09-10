Wright_2025_colchicine <- function() {
  description <- "Two-compartment population pharmacokinetic model for oral colchicine in people with gout (Wright 2025). Absorption is zero-order into the central compartment over a duration D1 = 0.99 h, preceded by an absorption lag that is fixed to 1.3 h in the subpopulation showing delayed absorption (Cmax,ss at or below Cmin,ss) and to 0 h otherwise; oral bioavailability F1 is fixed at 0.469 and is multiplied by (1 + theta_FORM) for the tablet formulations used in the literature-extracted studies. Clearance and inter-compartmental clearance are allometrically scaled on total body weight with a fixed exponent of 0.75; the central and peripheral volumes are scaled with a fixed exponent of 1. Concomitant statin use multiplies clearance by 0.66 (a ~30% reduction) and female sex multiplies both volumes by 0.53 (a ~50% reduction). Between-subject variability is estimated on clearance and fixed from Karatza 2021 on both volumes; residual error is combined proportional plus a small fixed additive term."
  reference   <- "Wright DFB, Hishe HZ, Dalbeth N, Horne A, Drake J, Haslett J, Stamp LK. The Influence of Patient Factors on the Population Pharmacokinetics of Colchicine: Implications for Safe and Effective Dosing. Clin Pharmacokinet. 2025;64(10):1519-1531. doi:10.1007/s40262-025-01551-y. Structure and fixed values taken from the final NONMEM control stream ($PROBLEM 66_FINAL_colchicine) reproduced in the Electronic Supplementary Information; final estimates from Table 2. Between-subject variances on V1 and V2 and the structural starting values are inherited from Karatza E, Ismailos G, Karalis V. Xenobiotica. 2021;51:643-656, doi:10.1080/00498254.2021.1909782 (ESI Table S1)."
  vignette    <- "Wright_2025_colchicine"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Standardised to 70 kg and allometrically scaled with a fixed exponent of 0.75 on CL and Q and a fixed exponent of 1 on V1 and V2 (Wright 2025 Table 2 footnotes a and b; ESI control stream SIZECL = (WTKG/70)**0.75, SIZEV = WTKG/70). Fat-free mass and normal fat mass were also tested as body-size descriptors and an estimated power exponent on clearance was tested; total body weight with the fixed exponents was retained. Observed range 57-150 kg, median 97 kg in the NZ Gout Study (Table 1).",
      source_name        = "WTKG"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed. Multiplies both V1 and V2 by e_sexf_vc_vp = 0.53, i.e. a ~50% lower apparent volume in women. The source column SEX is coded 1 = male, 0 = female (ESI control stream: IF(SEX.EQ.0) FSEX=THETA(13)); the canonical SEXF column reverses that orientation, so SEXF = 1 - SEX. Only 6 of 78 NZ Gout Study participants (7.7%) were women, so this effect rests on a small subgroup (RSE 25.8%, SIR 95% CI 0.29-0.86).",
      source_name        = "SEX (1 = male, 0 = female)"
    ),
    CONMED_STATIN = list(
      description        = "Concomitant statin (HMG-CoA reductase inhibitor) therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant statin)",
      notes              = "Time-fixed. Multiplies CL by e_conmed_statin_cl = 0.66, i.e. a ~30% lower clearance in statin users. Wright 2025 codes CONMED_STATIN = 1 for any concomitant statin at any dose; the retained cohort comprised atorvastatin (n = 12) and simvastatin (n = 4) (Table 2 footnote c), and rosuvastatin was named as an eligible statin in the Table 1 footnote. The effect is on clearance only -- Q carries the allometric term but no statin term (ESI control stream: TVCL = THETA(1)*STATINCL*SIZECL, TVQ = THETA(3)*SIZECL). 17 of 78 NZ Gout Study participants (22%) were taking a statin.",
      source_name        = "STATIN"
    ),
    MIX_LAGGED_ABS = list(
      description        = "Delayed-absorption subpopulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no absorption lag; TLAG2 fixed to 0 h)",
      notes              = "Time-fixed per subject. 1 = subject assigned to the delayed-absorption class, whose lag time is fixed at 1.3 h (TLAG1); 0 = no lag. Unlike the Bonate 2004 founding example this is NOT a NONMEM $MIXTURE class index -- Wright 2025 states that a $MIX mixture model for the bimodal lag was not supported, so the class was assigned deterministically in the data as a LAGFLAG column: 17 of the 78 NZ Gout Study participants had a steady-state Cmax at or below their Cmin, which the authors took as evidence of delayed absorption (Wright 2025 Results paragraph 2). 1.3 h is the longest post-dose sampling time recorded in the NZ Gout Study, not an estimate. For population simulation of an NZ-Gout-like cohort draw MIX_LAGGED_ABS ~ Bernoulli(17/78 = 0.218); the lag shifts the profile in time and does not change AUC or clearance.",
      source_name        = "LAGFLAG"
    ),
    FORM_COL_LIT = list(
      description        = "Literature-extracted colchicine tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the colchicine tablets used in the NZ Gout Study)",
      notes              = "Per-dose-occasion in principle, per-study in the source data. 1 = the dose was one of the colchicine tablet products used in the three published studies from which concentration-time data were digitised (Thomas 1989, Ferron 1996, Rochdi 1994); 0 = the NZ Gout Study tablets. Multiplies bioavailability as F1 = 0.469 * (1 + e_form_col_lit_fcentral), i.e. F1 = 0.727 for the literature formulations versus 0.469 for the NZ Gout Study. Derived from the source STUDY column (ESI control stream: IF(STUDY.EQ.1) FORM=1; IF(STUDY.GT.1) FORM=1+THETA(11)), so FORM_COL_LIT = as.integer(STUDY > 1). Imprecisely estimated (RSE 54.8%, SIR 95% CI 0.14-1.17). Set to 0 to simulate the gout-patient population the model was built to describe.",
      source_name        = "STUDY (1 = NZ Gout Study, > 1 = extracted literature studies)"
    )
  )

  # Covariates Wright 2025 screened but did not retain in the final model. These
  # are documentation only -- none is referenced in model(). The names are
  # descriptive labels for the paper's screen, not a claim that each is a
  # registered canonical: WT, SEXF, AGE, FFM and CRCL are registered, while
  # NFM, CONMED_ACEI, ADHERENCE, ETHNICITY, CONMED_CYP3A4_PGP_INHIB,
  # CONMED_CYP3A4_PGP_IND and CONMED_OTHER_CARDIO are not, because no model in
  # the library uses them.
  covariatesDataExcluded <- list(
    FFM = list(
      description = "Fat-free mass (Janmahasatian 2005 formula)",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as an alternative body-size descriptor against total body weight; not retained in the final model (Wright 2025 Methods 'Model Development')."
    ),
    NFM = list(
      description = "Normal fat mass (Anderson and Holford)",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as an alternative body-size descriptor against total body weight; not retained in the final model (Wright 2025 Methods 'Model Development')."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested on CL/F standardised to 6 L/h/70 kg using linear and power models and as a split renal / non-renal clearance model; not retained. The authors note the analysis was underpowered for a kidney-function effect -- the NZ Gout Study excluded eGFR < 30 mL/min and only nine participants had CLcr < 50 mL/min (Wright 2025 Discussion 'limitations')."
    ),
    CONMED_ACEI = list(
      description = "Concomitant angiotensin-converting enzyme inhibitor indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "A ~25% fractional reduction in colchicine clearance was observed with concomitant ACEI use but was deliberately deleted from the final model: ACEIs are not P-glycoprotein or CYP3A4 substrates, inhibitors or inducers, the ACEI cohort was older and had lower kidney function, and the authors judged the effect to have no mechanistic basis (Wright 2025 Discussion paragraph 4)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested with a linear model and as a fractional effect on binned categories; not retained (Wright 2025 Methods 'Model Development')."
    ),
    ADHERENCE = list(
      description = "Treatment adherence from pill counts at month 3",
      units       = "%",
      type        = "continuous",
      notes       = "Tested with a linear model and as a fractional effect on binned categories; not retained. Median 97% (range 30-113%) (Wright 2025 Table 1)."
    ),
    ETHNICITY = list(
      description = "Self-reported ethnicity (Maori, Pacific Peoples, NZ European, Other)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Tested as a fractional effect; not retained (Wright 2025 Methods 'Model Development' and Results)."
    ),
    CONMED_CYP3A4_PGP_INHIB = list(
      description = "Concomitant strong or moderate CYP3A4 and/or P-glycoprotein inhibitor indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Grouped per the Flockhart Table and UpToDate and tested as a fractional effect on colchicine clearance; not retained. n = 8 (10%): amlodipine (4), diltiazem (2), verapamil (1), carvedilol (1) (Wright 2025 Table 1 footnote)."
    ),
    CONMED_CYP3A4_PGP_IND = list(
      description = "Concomitant strong or moderate CYP3A4 and/or P-glycoprotein inducer indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a fractional effect on colchicine clearance; not retained. n = 14 (18%): 13 who had taken prednisone for a gout flare since the previous visit and one taking carbamazepine (Wright 2025 Table 1 footnote)."
    ),
    CONMED_OTHER_CARDIO = list(
      description = "Concomitant diuretic, angiotensin receptor blocker, alpha-blocker, beta-blocker or calcium channel blocker indicators",
      units       = "(binary)",
      type        = "binary",
      notes       = "Each drug class was tested individually as a fractional effect on colchicine clearance; none was retained in the final model (Wright 2025 Methods 'Model Development', Table 1)."
    )
  )

  compartmentData <- list(
    central     = list(analyte = "colchicine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "colchicine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 91L,
    n_studies        = 4L,
    n_observations   = "364 colchicine plasma concentrations: 153 from 78 NZ Gout Study participants plus 211 digitised from 13 individuals in three published studies (Wright 2025 Results paragraph 1).",
    age_range        = "27-90 years (median 59) in the NZ Gout Study; not reported for the extracted individuals",
    weight_range     = "57-150 kg (median 97) in the NZ Gout Study; not reported for the extracted individuals",
    sex_female_pct   = 7.7,
    race_ethnicity   = c(`NZ European` = 56, Other = 17, `Pacific Peoples` = 14, Maori = 13),
    disease_state    = "Adults with gout commencing allopurinol urate-lowering therapy (NZ Gout Study, ACTRN 12618001179224), plus ten healthy volunteers and one participant each with liver disease, gout and kidney disease whose concentration-time data were digitised from three published studies.",
    dose_range       = "Colchicine 0.5 mg orally once daily for 6 months in the NZ Gout Study (sampled at month 3); single oral doses of 0.5, 1 and 2 mg and single intravenous doses of 0.5 and 2 mg in the extracted literature data.",
    regions          = "New Zealand (Auckland and Christchurch) for the NZ Gout Study; the extracted studies were conducted elsewhere.",
    renal_function   = "Creatinine clearance (Cockcroft-Gault) 81 mL/min (range 32-132) in the NZ Gout Study; participants with eGFR < 30 mL/min were excluded, and only nine had CLcr < 50 mL/min.",
    co_medication    = "NZ Gout Study: statins 22%, ACEIs 21%, beta-blockers 13%, ARBs 12%, diuretics 10%, calcium channel blockers 10%, CYP3A4 and/or P-glycoprotein inhibitors 10%, CYP3A4 and/or P-glycoprotein inducers 18% (Wright 2025 Table 1).",
    notes            = "Sampling in the NZ Gout Study was sparse -- two samples per participant at month 3, one just before the daily dose and one 30-60 min after it. The literature-extracted intensively sampled profiles were added specifically to stabilise the structural model. Baseline demographics are in Wright 2025 Table 1; individual-level demographics were not available for the extracted individuals, so the covariate model is identified by the NZ Gout Study cohort. The virtual gout population used for the paper's dosing simulations was built from 309 people with gout across five New Zealand studies (ESI Table S2: weight 101 kg (51-172), 45% on statins, 14% female)."
  )

  ini({
    # ======================================================================
    # Structural parameters -- Wright 2025 Table 2, "Final model" column.
    # Values are for a 70 kg male not taking a statin. The paper's Eqs. in
    # Results write these as CL/F, V1/F, Q/F and V2/F, but the final NONMEM
    # control stream in the ESI estimates CL, V1, Q and V2 with the oral
    # bioavailability handled separately as F1 = 0.469 * FORM (ESI $PK:
    # TVCL = THETA(1)*STATINCL*SIZECL ... TVF1 = THETA(9)*FORM). This file
    # follows the control stream, so f(central) below carries F1 explicitly
    # and the apparent oral clearance is cl / f1 = 19.1 / 0.469 = 40.7 L/h.
    # ======================================================================
    lcl <- log(19.1);   label("Clearance CL for a 70 kg male not taking a statin (L/h)")                     # Wright 2025 Table 2: theta_CL = 19.1 L/h/70 kg (RSE 9.0%; SIR 19.5 [16.6-22.5])
    lvc <- log(245.5);  label("Central volume of distribution V1 for a 70 kg male (L)")                      # Wright 2025 Table 2: theta_V1 = 245.5 L/70 kg (RSE 16.2%; SIR 243.1 [185.5-300.7])
    lq  <- log(29.9);   label("Inter-compartmental clearance Q for a 70 kg subject (L/h)")                   # Wright 2025 Table 2: theta_Q = 29.9 L/h/70 kg (RSE 21.9%; SIR 32.8 [20.7-45.0])
    lvp <- log(821.8);  label("Peripheral volume of distribution V2 for a 70 kg male (L)")                   # Wright 2025 Table 2: theta_V2 = 821.8 L/70 kg (RSE 28.4%; SIR 893.6 [472.2-1314])
    ld1 <- log(0.99);   label("Duration of the zero-order input into the central compartment D1 (h)")        # Wright 2025 Table 2: D1 = 0.99 h (RSE 15.2%; SIR 1.0 [0.8-1.2])

    # Lag time and bioavailability are both fixed, not estimated.
    ltlag     <- fixed(log(1.3));    label("Absorption lag time in the delayed-absorption subpopulation (h)") # Wright 2025 Table 2: TLAG1 = 1.3 h (fixed); TLAG2 = 0 h (fixed) applies when MIX_LAGGED_ABS = 0. 1.3 h is the longest post-dose sampling time in the NZ Gout Study, not an estimate (Results paragraph 2)
    lfcentral <- fixed(log(0.469));  label("Oral bioavailability F1 for the NZ Gout Study tablets (fraction)") # ESI final control stream $THETA 9: "0.469 FIX ; 9. F1_oral"; Wright 2025 Table 2 rounds this to F1 = 0.47 (fixed) and the Results state "Bioavailability was fixed to 0.47 in the final model"

    # ======================================================================
    # Allometric size model -- fixed exponents, standardised to 70 kg.
    # ======================================================================
    e_wt_cl_q  <- fixed(0.75); label("Body-weight allometric exponent shared by CL and Q (unitless)")        # Wright 2025 Table 2 footnote a: "allometrically scaled to an exponent of 3/4, modelled as CL or Q = THETA*(weight/70)^3/4"; ESI $PK SIZECL = (WTKG/70)**0.75
    e_wt_vc_vp <- fixed(1);    label("Body-weight allometric exponent shared by V1 and V2 (unitless)")       # Wright 2025 Table 2 footnote b: "Volume expressed per 70 kg total body weight, modelled as V1 or V2 = THETA*(weight/70)"; ESI $PK SIZEV = WTKG/70

    # ======================================================================
    # Covariate effects -- all coded as multiplicative fractional effects,
    # i.e. the parameter is multiplied by theta when the indicator is 1 and
    # left unchanged when it is 0 (ESI $PK IF-blocks).
    # ======================================================================
    e_conmed_statin_cl      <- 0.66; label("Multiplicative fractional effect of concomitant statin use on CL (unitless)")               # Wright 2025 Table 2: theta_statin = 0.66 (RSE 12.4%; SIR 0.66 [0.50-0.82]) -- a ~30% reduction in clearance
    e_sexf_vc_vp            <- 0.53; label("Multiplicative fractional effect of female sex on V1 and V2 (unitless)")                    # Wright 2025 Table 2: theta_SEX = 0.53 (RSE 25.8%; SIR 0.57 [0.29-0.86]) -- a ~50% reduction in both volumes
    e_form_col_lit_fcentral <- 0.55; label("Fractional increment in F1 for the literature-extracted tablet formulations (unitless)")    # Wright 2025 Table 2: theta_FORM = 0.55 (RSE 54.8%; SIR 0.65 [0.14-1.17]); Table 2 footnote d: "FORM = 1 for the NZ Gout Study and FORM = 1 + THETA for the extracted data"

    # ======================================================================
    # Between-subject variability. Wright 2025 reports omega as a CV% that
    # equals sqrt(variance) * 100, not the log-normal sqrt(exp(var) - 1):
    # the two fixed volume omegas are 35.5% and 53.3% in Table 2 and 0.126
    # and 0.284 in the ESI $OMEGA block, and sqrt(0.126) = 0.355 while
    # sqrt(0.284) = 0.533. The clearance omega therefore reads back as
    # variance = 0.338^2 = 0.1142. Q, D1, TLAG and F1 carry no eta.
    # ======================================================================
    etalcl ~ 0.1142;          label("BSV on clearance (variance on the log scale)")            # Wright 2025 Table 2: omega_CL = 33.8 CV% (RSE 21.6%; SIR 35.5 [28.2-41.6]); 0.338^2 = 0.1142
    etalvc ~ fixed(0.126);    label("BSV on central volume (variance on the log scale)")       # Wright 2025 Table 2: omega_V1 = 35.5 CV% (fixed); ESI $OMEGA "0.126 FIX ; PPV_V1", inherited from Karatza 2021 (ESI Table S1: omega_V1 variance 0.126)
    etalvp ~ fixed(0.284);    label("BSV on peripheral volume (variance on the log scale)")    # Wright 2025 Table 2: omega_V2 = 53.3 CV% (fixed); ESI $OMEGA "0.284 FIX ; PPV_V2", inherited from Karatza 2021

    # ======================================================================
    # Residual unexplained variability -- combined proportional plus a small
    # fixed additive term. ESI $ERROR: W = SQRT(THETA(6)**2 + THETA(7)**2 *
    # F*F); Y = IPRED + W*EPS(1) with $SIGMA 1 FIX, so THETA(6) and THETA(7)
    # are standard deviations on the additive and proportional arms.
    # ======================================================================
    propSd <- 0.377;         label("Proportional residual SD (fraction)")     # Wright 2025 Table 2: sigma_prop = 37.7 CV% (RSE 8.3%; SIR 38.1 [33.8-42.2])
    addSd  <- fixed(0.006);  label("Additive residual SD (ng/mL)")            # Wright 2025 Table 2: sigma_add = 0.006 (fixed); ESI $THETA 6 "0.006 FIX ; 6. RUV_ADD". Table 2 tags the units as umol/L, which cannot be right -- 0.006 umol/L is 2.4 ng/mL, comparable to the whole therapeutic range; the assay and every reported concentration are in ng/mL, so the tag is read as ug/L = ng/mL. See vignette Errata.
  })

  model({
    # ---- Individual parameters --------------------------------------------
    # Wright 2025 Results, final covariate model:
    #   CL/F = theta_CL * (WT/70)^0.75 * theta_statin
    #   V1/F = theta_V1 * (WT/70)      * theta_SEX
    #   Q/F  = theta_Q  * (WT/70)^0.75
    #   V2/F = theta_V2 * (WT/70)      * theta_SEX
    # The binary effects are written here as (1 - (1 - theta) * IND) so the
    # multiplier is exactly 1 in the reference category (male, no statin)
    # and exactly theta in the covariate-positive category, matching the
    # ESI $PK IF-blocks.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q *
      (1 - (1 - e_conmed_statin_cl) * CONMED_STATIN)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp *
      (1 - (1 - e_sexf_vc_vp) * SEXF)
    q  <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp *
      (1 - (1 - e_sexf_vc_vp) * SEXF)

    d1 <- exp(ld1)

    # Deterministic delayed-absorption class: lag = 1.3 h when
    # MIX_LAGGED_ABS = 1 and 0 h otherwise (ESI $PK: IF(LAGFLAG.EQ.1)
    # TVLAG = THETA(8) ELSE TVLAG = THETA(10), with THETA(10) = 0 FIX).
    tlag <- exp(ltlag) * MIX_LAGGED_ABS

    # F1 = 0.469 for the NZ Gout Study tablets and 0.469 * 1.55 = 0.727 for
    # the literature-extracted formulations (ESI $PK: TVF1 = THETA(9)*FORM).
    f1 <- exp(lfcentral) * (1 + e_form_col_lit_fcentral * FORM_COL_LIT)

    # ---- Micro-constants ---------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Disposition -------------------------------------------------------
    # NONMEM ADVAN3 TRANS4 with a zero-order input into compartment 1: there
    # is no depot state, the dose enters the central compartment at a
    # constant rate over d1 hours after the lag. Dose records must carry
    # rate = -2 so that rxode2 uses the modelled duration dur(central) = d1;
    # without it the dose collapses to an instantaneous bolus.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    f(central)    <- f1
    dur(central)  <- d1
    alag(central) <- tlag

    # Amounts are in mg and volumes in L, so central / vc is mg/L; multiply
    # by 1000 to report ng/mL, the unit of every concentration in the paper.
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd) + add(addSd)
  })
}

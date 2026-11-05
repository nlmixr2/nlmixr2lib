Hishe_2023_oxypurinol <- function() {
  description <- "One-compartment population PK model for oxypurinol (the active metabolite of allopurinol) in adults with gout (Hishe 2023 Clin Transl Sci). Pooled analysis of five New Zealand studies (NONMEM 7.3, FOCEI): first-order absorption with Ka fixed at 0.993 1/h, one-compartment distribution and first-order elimination. Apparent oral clearance CL/F is power-scaled on Cockcroft-Gault creatinine clearance standardised to 70 kg (CRCL in L/h, reference 6 L/h), multiplicatively reduced by concomitant diuretic use (CONMED_DIURETIC) and increased in Pacific peoples (RACE_PACIFIC), and power-scaled on total body weight; apparent volume V/F is linearly weight-scaled. The dose entered into the model is allopurinol oral mg; the implicit 1:1 molar conversion to oxypurinol uses a 0.9 conversion factor and the allopurinol molecular weight 136.11 g/mol, so Cc is oxypurinol concentration in umol/L. Seven urate-transporter gene variants (ABCG2, SLC2A9, SLC17A1, SLC22A12, SLC22A11, ABCC4) were screened but none were retained in the final model. A separate combined additive + proportional residual error was estimated for each of the five studies."
  reference <- "Hishe HZ, Stocker SL, Stamp LK, Dalbeth N, Merriman TR, Phipps-Green A, Wright DFB. The impact of genetic variability in urate transporters on oxypurinol pharmacokinetics. Clin Transl Sci. 2023 Mar;16(3):422-428. doi:10.1111/cts.13460"
  vignette <- "Hishe_2023_oxypurinol"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "allopurinol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "oxypurinol", units = "mg (allopurinol-equivalent)", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault formula and standardised to 70 kg (NOT BSA-normalised), expressed as L/h.",
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = "Hishe 2023 Methods: 'Creatinine clearance (CrCl) was calculated according the Cockcroft-Gault formula. CrCl was normalized to a standard of 6 L/h/70 kg (100 ml/min/70 kg) and was assessed on oxypurinol clearance as linear and power models.' The final model uses the power form. Values are stored in L/h (so the normalisation standard 6 L/h = 100 mL/min/70 kg; the cohort median 63 mL/min ~= 3.78 L/h). Enters CL/F via power scaling (CRCL/6)^0.58. Same Cockcroft-Gault-standardised-to-70 kg form as the sibling Wright_2016_allopurinol.R; the NONMEM column is CLCRCG.",
      source_name = "CLCRCG"
    ),
    WT = list(
      description = "Total body weight (TBW).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Hishe 2023 Methods / Results. Enters CL/F via power scaling (WT/70)^0.34 and V/F via linear scaling (WT/70)^1 (the volume exponent is fixed at the theoretical value 1.0; Table 2 reports no V weight exponent). NONMEM column WTKG.",
      source_name = "WTKG"
    ),
    CONMED_DIURETIC = list(
      description = "Concomitant diuretic indicator; class composition not individually itemised by Hishe 2023 (the NONMEM dataset carries a DIURETICTYPE column but the final covariate uses the pooled any-diuretic indicator).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant diuretic or no data)",
      notes = "Hishe 2023 Appendix S2 $PK: 'IF(DIURETIC.EQ.1) THEN FRUCL=THETA(14) ELSE FRUCL=1 ; NO DIURETIC OR NO DATA'. Multiplicative effect on CL/F: cl *= 0.79^CONMED_DIURETIC (Table 2 final theta_d = 0.79), i.e. 21% lower apparent oxypurinol clearance on a diuretic, in line with Hishe 2023 Results ('Diuretic use was found to decrease oxypurinol clearance by 21%'). 37% of the analysis cohort (112 of 300) were on a diuretic (Table S4). Narrower-class diuretic definitions are used by the sibling Wright_2016_allopurinol.R (thiazide + loop only) and a broader one by Stocker_2012_oxypurinol.R (thiazide + loop + spironolactone); Hishe 2023 does not state the class composition, so users must populate the column per the diuretic classes they treat as relevant.",
      source_name = "DIURETIC"
    ),
    RACE_PACIFIC = list(
      description = "Pacific peoples ethnicity indicator (1 = Pacific peoples, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (European / NZ European, Maori, Chinese, Indian, or Other)",
      notes = "Hishe 2023 Appendix S2 $PK: 'IF(ETHNICITY.EQ.3) THEN FRACE=THETA(15) ; PACIFIC ISLAND ETHNICITY ELSE FRACE=1 ; EVERY BODY ELSE; REFERENCE, MAORI, OTHER'. Multiplicative effect on CL/F: cl *= 1.25^RACE_PACIFIC (Table 2 final theta_e = 1.25), i.e. ~25% higher apparent oxypurinol clearance in Pacific peoples (Hishe 2023 Results: 'Pacific Island ethnicity increased oxypurinol clearance by around 25%'; the mechanism is stated to be unclear). 21% of the analysis cohort (62 of 300) were Pacific peoples (Table S4). The reference group explicitly includes Maori, so RACE_PACIFIC is distinct from an indigenous-New-Zealand indicator. NONMEM column ETHNICITY (coded 3 = Pacific).",
      source_name = "ETHNICITY"
    ),
    STUDY_STAMP2011 = list(
      description = "Study-1 indicator in the Hishe 2023 five-study pooled analysis (Stamp 2011 Arthritis Rheum; ANZCTR 012606000276550).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (another of the five pooled studies; Study 5 / Stamp 2017 is the all-zeros reference)",
      notes = "Selects the Study-1 combined residual error (additive fixed at 0.0001 umol/L, proportional 43.5% CV; Hishe 2023 Table 2). Acts on no structural parameter. NONMEM STUDY.EQ.1.",
      source_name = "STUDY"
    ),
    STUDY_STAMP2012 = list(
      description = "Study-2 indicator in the Hishe 2023 five-study pooled analysis (Stamp 2012 Rheumatology furosemide-interaction study; ANZCTR 12609000529246).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (another of the five pooled studies; Study 5 / Stamp 2017 is the all-zeros reference)",
      notes = "Selects the Study-2 combined residual error (additive 10.0 umol/L, proportional fixed at 0.0001; Hishe 2023 Table 2). Acts on no structural parameter. NONMEM STUDY.EQ.2.",
      source_name = "STUDY"
    ),
    STUDY_STAMP2013 = list(
      description = "Study-3 indicator in the Hishe 2023 five-study pooled analysis (Stamp 2013 Arthritis Rheum vitamin-C study; ANZCTR 12610000545066).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (another of the five pooled studies; Study 5 / Stamp 2017 is the all-zeros reference)",
      notes = "Selects the Study-3 combined residual error (additive fixed at 0.0001 umol/L, proportional 42.2% CV; Hishe 2023 Table 2). Acts on no structural parameter. NONMEM STUDY.EQ.3.",
      source_name = "STUDY"
    ),
    STUDY_WRIGHT2016 = list(
      description = "Study-4 indicator in the Hishe 2023 five-study pooled analysis (the low-dose data reported in Wright 2016 / Wright 2018 Br J Clin Pharmacol).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (another of the five pooled studies; Study 5 / Stamp 2017 is the all-zeros reference)",
      notes = "Selects the Study-4 combined residual error (additive 0.93 umol/L, proportional 34.8% CV; Hishe 2023 Table 2). Acts on no structural parameter. NONMEM STUDY.EQ.4.",
      source_name = "STUDY"
    )
  )

  covariatesDataExcluded <- list(
    SNP_ABCG2_RS2231142 = list(
      description = "ABCG2 (BCRP) rs2231142 (Q141K) risk-allele (T) carrier indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a fractional effect on oxypurinol clearance. Univariate effect 1.24 [1.08, 1.40] (+24% CL for T-allele carriers) but NOT retained in the multivariate final model (Hishe 2023 Results / Discussion: the univariate signal is confounded by correlated CrCl, weight, age, diuretic use and ethnicity). Participants with >=1 variant allele were pooled; missing genotypes assigned to the reference allele."
    ),
    SNP_ABCG2_RS10011796 = list(
      description = "ABCG2 (BCRP) rs10011796 risk-allele (T) carrier indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; univariate effect 1.09 [0.91, 1.27]; not retained (Hishe 2023 Results)."
    ),
    SNP_SLC2A9_RS11942223 = list(
      description = "SLC2A9/GLUT9 rs11942223 risk-allele (T) carrier indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; univariate effect 1.63 [0.61, 2.65]; not retained (Hishe 2023 Results)."
    ),
    SNP_SLC17A1_RS1183201 = list(
      description = "SLC17A1/NPT1 rs1183201 risk-allele (T) carrier indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; univariate effect 1.22 [1.02, 1.42] (+22% CL for T-allele carriers) but NOT retained in the multivariate final model (Hishe 2023 Results / Discussion; confounded by correlated participant factors)."
    ),
    SNP_SLC22A12_RS3825018 = list(
      description = "SLC22A12/URAT1 rs3825018 risk-allele (G) carrier indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; univariate effect 0.97 [0.85, 1.09]; not retained (Hishe 2023 Results)."
    ),
    SNP_SLC22A11_RS17300741 = list(
      description = "SLC22A11/OAT4 rs17300741 risk-allele (A) carrier indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; univariate effect 0.90 [0.76, 1.04]; not retained (Hishe 2023 Results)."
    ),
    SNP_ABCC4_RS4148500 = list(
      description = "ABCC4/MRP4 rs4148500 risk-allele (A) carrier indicator (population-specific variant absent in Caucasians).",
      units = "(binary)",
      type = "binary",
      notes = "Screened; univariate effect 1.02 [0.75, 1.29]; not retained (Hishe 2023 Results). 60% missing genotype in this cohort (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 300L,
    n_studies = 5L,
    age_range = "23 to 86 years",
    age_median = "61 years",
    weight_range = "51 to 186 kg total body weight",
    weight_median = "102 kg total body weight",
    sex_female_pct = 12.0,
    race_ethnicity = "European / NZ European 59% (177/300), Pacific peoples 21% (62/300), Maori 16% (48/300), Chinese 2% (5/300), Indian 2% (6/300), Other 1% (2/300) (Hishe 2023 Table S4 total column).",
    disease_state = "Adults with gout receiving allopurinol; participants on dialysis were excluded. Plasma oxypurinol concentrations <20 umol/L were excluded as presumed suboptimal adherence.",
    dose_range = "Allopurinol 50-900 mg/day oral (median 300 mg/day; Hishe 2023 Table S4).",
    regions = "New Zealand (Christchurch and Auckland).",
    crcl_range = "8 to 167 mL/min Cockcroft-Gault (median 63 mL/min ~= 3.78 L/h; Hishe 2023 Table S4).",
    co_medication = "Diuretic use (37%, retained on CL/F). Also recorded but not retained: ACE inhibitors (46%), beta-blockers (43%), ARBs (11%).",
    genotypes_tested = "Seven urate-transporter SNPs were tested (ABCG2 rs2231142 and rs10011796, SLC2A9 rs11942223, SLC17A1 rs1183201, SLC22A12 rs3825018, SLC22A11 rs17300741, ABCC4 rs4148500). None were retained in the final model; see covariatesDataExcluded. ABCG2 rs2231142 (+24%) and SLC17A1 rs1183201 (+22%) were significant univariately but confounded by correlated participant factors.",
    samples = "2841 plasma oxypurinol concentrations from 300 people with gout in the analysis set (3666 from 325 before exclusions; Hishe 2023 Table S1 / S4).",
    studies_notes = "Five-study pooled analysis. Study 1 (n=73; Stamp 2011 Arthritis Rheum, ANZCTR 012606000276550), Study 2 (n=10; Stamp 2012 Rheumatology furosemide interaction, ANZCTR 12609000529246), Study 3 (n=30; Stamp 2013 Arthritis Rheum vitamin C, ANZCTR 12610000545066), Study 4 (n=19; low-dose data reported in Wright 2016/2018 BJCP), Study 5 (n=176; Stamp 2017 Ann Rheum Dis dose-escalation, ANZCTR 12611000845932). A separate combined additive + proportional residual error was estimated per study; Study 5 is the all-zeros reference for the STUDY_* indicator columns because it is the largest sub-cohort (176 of 300).",
    notes = "Base PK model (one-compartment, first-order absorption and elimination, combined residual) carried from Stamp 2017 (Hishe 2023 Appendix S1). Allopurinol plasma concentrations were not measured in all studies, so a parent-metabolite model was not fitted; a 0.9 conversion factor relates the allopurinol dose to oxypurinol per prior work (Day 2007). 1000 nonparametric bootstraps supported the final estimates (Hishe 2023 Table 2)."
  )

  ini({
    # ----- Structural PK parameters (Hishe 2023 Table 2 'Final model'). -----
    # All structural PK parameters are APPARENT (CL/F, V/F): the dose is
    # allopurinol mg and only oxypurinol concentrations were measured, so the
    # bioavailability, the 0.9 allopurinol->oxypurinol conversion factor and
    # the molar conversion are absorbed into the apparent parameters. The
    # model() block converts mg of allopurinol to umol/L of oxypurinol using
    # the allopurinol molecular weight 136.11 g/mol and the 0.9 factor.
    lka <- fixed(log(0.993)); label("Apparent first-order absorption rate constant Ka (1/h)") # Hishe 2023 Table 2 / Appendix S2 $THETA: Ka = 0.993 fixed
    lcl <- log(0.97); label("Apparent oral clearance CL/F at reference CRCL = 6 L/h, no diuretic, non-Pacific, WT = 70 kg (L/h)") # Hishe 2023 Table 2 final: theta_CL = 0.97 L/h (RSE 4.7%); bootstrap 0.97 [0.88, 1.06]
    lvc <- log(49); label("Apparent oral volume V/F at reference WT = 70 kg (L)") # Hishe 2023 Table 2 final: theta_V = 49 L (RSE 5.4%); bootstrap 49 [43, 54]

    # Covariate effects on CL/F.
    e_crcl_cl <- 0.58; label("Power exponent of (CRCL/6) on CL/F (unitless)") # Hishe 2023 Table 2 final: theta_RFexp = 0.58 (RSE 7.7%); bootstrap 0.58 [0.49, 0.66]; power on (CRCL/6)
    e_conmed_diuretic_cl <- 0.79; label("Multiplicative factor on CL/F when CONMED_DIURETIC = 1 (unitless)") # Hishe 2023 Table 2 final: theta_d = 0.79 (RSE 3.5%); bootstrap 0.79 [0.74, 0.85]; CL reduced 21% on a diuretic
    e_race_pacific_cl <- 1.25; label("Multiplicative factor on CL/F when RACE_PACIFIC = 1 (unitless)") # Hishe 2023 Table 2 final: theta_e = 1.25 (RSE 4.9%); bootstrap 1.25 [1.12, 1.37]; CL raised 25% in Pacific peoples
    e_wt_cl <- 0.34; label("Power exponent of (WT/70) on CL/F (unitless)") # Hishe 2023 Table 2 final: theta_WTexp = 0.34 (RSE 27.2%); bootstrap 0.35 [0.17, 0.53]; power on (WT/70)

    # V/F weight exponent held fixed at the theoretical value 1.0: the Results
    # equation is V/F = theta_V * (WT/70) with no exponent reported in Table 2.
    e_wt_vc <- fixed(1.0); label("Power exponent of (WT/70) on V/F (unitless)") # Hishe 2023 Results equation: V/F = theta_V * (WT/70)^1 (linear; no Table 2 exponent row)

    # ----- IIV (Hishe 2023 Table 2 'Final model'). -----
    # The paper reports between-subject variability as 'omega CV%', which is
    # the log-scale standard deviation expressed as a percentage (the NONMEM
    # tradition CV ~= omega). This reading is confirmed by the Appendix S2
    # control stream, which FIXes the Ka variance at omega^2 = 0.388, and
    # sqrt(0.388) = 0.623 = the reported 62.3% CV. The encoding used here is
    # therefore omega^2 = (CV% / 100)^2:
    #   var CL = 0.298^2 = 0.088804
    #   var V  = 0.413^2 = 0.170569
    #   var Ka = 0.623^2 = 0.388 (control-stream fixed value)
    #   cov(CL, V) = 0.41 * sqrt(0.088804 * 0.170569) = 0.050464
    etalcl + etalvc ~ c(
      0.088804,
      0.050464, 0.170569
    ) # Hishe 2023 Table 2: omega CL 29.8% CV, omega V 41.3% CV, corr(CL,V) 0.41
    etalka ~ fixed(0.388) # Hishe 2023 Table 2: omega Ka 62.3% CV, not estimated; Appendix S2 $OMEGA 0.388

    # ----- Residual error (Hishe 2023 Table 2 'Final model'). -----
    # A separate combined additive + proportional error was estimated per
    # study. The NONMEM form is W = sqrt(add^2 + prop^2 * F^2), so each theta
    # is a standard deviation: 'add' in umol/L, 'prop' as a fraction. Within
    # each study one component was effectively switched off by fixing it at
    # 0.0001. Study 5 (Stamp 2017, the largest sub-cohort) is the all-zeros
    # reference for the STUDY_* indicators.
    addSd_s1 <- fixed(0.0001); label("Study-1 additive residual error (umol/L)") # Hishe 2023 Table 2 Study 1: eps add = 0.0001 (fixed) umol/L
    propSd_s1 <- 0.435; label("Study-1 proportional residual error (fraction)") # Hishe 2023 Table 2 Study 1: eps prop = 43.5% CV (RSE 6.4%)
    addSd_s2 <- 10.0; label("Study-2 additive residual error (umol/L)") # Hishe 2023 Table 2 Study 2: eps add = 10.0 umol/L (RSE 14.3%)
    propSd_s2 <- fixed(0.0001); label("Study-2 proportional residual error (fraction)") # Hishe 2023 Table 2 Study 2: eps prop = 0.0001 (fixed)
    addSd_s3 <- fixed(0.0001); label("Study-3 additive residual error (umol/L)") # Hishe 2023 Table 2 Study 3: eps add = 0.0001 (fixed) umol/L
    propSd_s3 <- 0.422; label("Study-3 proportional residual error (fraction)") # Hishe 2023 Table 2 Study 3: eps prop = 42.2% CV (RSE 8.4%)
    addSd_s4 <- 0.93; label("Study-4 additive residual error (umol/L)") # Hishe 2023 Table 2 Study 4: eps add = 0.93 umol/L (RSE 16.7%)
    propSd_s4 <- 0.348; label("Study-4 proportional residual error (fraction)") # Hishe 2023 Table 2 Study 4: eps prop = 34.8% CV (RSE 9.4%)
    addSd_s5 <- fixed(0.0001); label("Study-5 additive residual error (umol/L)") # Hishe 2023 Table 2 Study 5: eps add = 0.0001 (fixed) umol/L
    propSd_s5 <- 0.485; label("Study-5 proportional residual error (fraction)") # Hishe 2023 Table 2 Study 5: eps prop = 48.5% CV (RSE 3.6%)
  })

  model({
    # Reference values for the final covariate model (Hishe 2023 Results
    # equations). CRCL is in L/h standardised to 70 kg; the normalisation
    # standard is 6 L/h (= 100 mL/min/70 kg).
    ref_crcl <- 6
    ref_wt <- 70

    # Allopurinol molecular weight (g/mol = mg/mmol) and the 0.9
    # allopurinol->oxypurinol conversion factor (Hishe 2023 Appendix S1,
    # citing Day 2007). Used to convert mg of allopurinol (the dose unit) in
    # the central compartment to umol/L of oxypurinol (the observed unit). The
    # molar reading (0.9 of the allopurinol molecules become oxypurinol, 1:1
    # molar) matches the same research group's Wright 2013 convention; see the
    # validation vignette for the alternative mass reading and why the molar
    # form was chosen.
    mwAllop <- 136.11
    convFactor <- 0.9

    # ----- Individual PK parameters -----
    # Hishe 2023 final-model equations:
    #   CL/F = theta_CL * (CRCL/6)^0.58 * 0.79^DIURETIC * 1.25^PACIFIC * (WT/70)^0.34
    #   V/F  = theta_V * (WT/70)^1
    #   Ka   = 0.993 (fixed)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (CRCL / ref_crcl)^e_crcl_cl *
      e_conmed_diuretic_cl^CONMED_DIURETIC *
      e_race_pacific_cl^RACE_PACIFIC *
      (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    kel <- cl / vc

    # ----- ODE system -----
    # depot and central carry mg of allopurinol (the dose mass); the 0.9
    # conversion factor and molar conversion to oxypurinol are applied at the
    # observation. Linear one-compartment kinetics, so the constant scaling at
    # observation does not affect the disposition.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ----- Observation -----
    # Oxypurinol concentration in umol/L: mg allopurinol in central / vc gives
    # mg/L; * 1000 / 136.11 converts to umol/L of allopurinol-equivalent; * 0.9
    # applies the molar allopurinol->oxypurinol conversion factor.
    Cc <- convFactor * central / vc * 1000 / mwAllop

    # ----- Residual error: per-study combined additive + proportional -----
    # Study 5 (Stamp 2017) is the reference: all STUDY_* indicators 0 selects
    # the Study-5 error. Exactly one study is active per observation.
    s5frac <- 1 - STUDY_STAMP2011 - STUDY_STAMP2012 - STUDY_STAMP2013 - STUDY_WRIGHT2016
    addSdCc <- addSd_s1 * STUDY_STAMP2011 +
      addSd_s2 * STUDY_STAMP2012 +
      addSd_s3 * STUDY_STAMP2013 +
      addSd_s4 * STUDY_WRIGHT2016 +
      addSd_s5 * s5frac
    propSdCc <- propSd_s1 * STUDY_STAMP2011 +
      propSd_s2 * STUDY_STAMP2012 +
      propSd_s3 * STUDY_STAMP2013 +
      propSd_s4 * STUDY_WRIGHT2016 +
      propSd_s5 * s5frac
    Cc ~ add(addSdCc) + prop(propSdCc)
  })
}

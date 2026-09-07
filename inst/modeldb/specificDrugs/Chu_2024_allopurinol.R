Chu_2024_allopurinol <- function() {
  description <- paste(
    "Joint parent-metabolite population pharmacokinetic model for intravenous",
    "allopurinol and its active metabolite oxypurinol in 14 neonates with",
    "critical congenital heart disease (CCHD) undergoing cardiac surgery with",
    "cardiopulmonary bypass (CPB) in the CRUCIAL trial (Chu 2024). Structural",
    "model: two-compartment allopurinol disposition feeding a one-compartment",
    "oxypurinol disposition, with full (formation fraction 1) conversion of",
    "allopurinol to oxypurinol and auto-inhibition of that conversion by",
    "oxypurinol (Imax fixed to 1, IC50 fixed to 1.1 mg/L). The allopurinol",
    "central compartment (V1 = 0.1 L per 3.5 kg) and intercompartmental",
    "clearance (Q1 = 6.97 L/h per 3.5 kg) were fixed to capture the 10-minute",
    "post-dose peak sampled during CPB and carry no physiologic",
    "interpretation. Clearances and volumes are allometrically scaled to a",
    "3.5 kg reference body weight with fixed exponents 0.75 and 1. Three",
    "mutually exclusive perioperative periods carry different disposition:",
    "during the postnatal-preoperative period both clearances rise with",
    "postnatal age along a fixed sigmoidal recovery curve (TM50 4.2 days,",
    "Hill 2.98; maximum fold increase 3 for allopurinol and 1.35 for",
    "oxypurinol); during CPB and after CPB the clearances and volumes are",
    "instead fixed multiples of their at-birth values. Between-subject",
    "variability is carried on allopurinol clearance and oxypurinol volume,",
    "and between-occasion variability on both clearances across the three",
    "periods. Residual error is proportional plus an additive component fixed",
    "to LLOQ/2 for each analyte.")
  reference <- paste(
    "Chu WY, Nijman M, Stegeman R, Breur JMPJ, Jansen NJG, Nijman J, van Loon K,",
    "Koomen E, Allegaert K, Benders MJNL, Dorlo TPC, Huitema ADR;",
    "the CRUCIAL trial consortium.",
    "Population Pharmacokinetics and Target Attainment of Allopurinol and",
    "Oxypurinol Before, During, and After Cardiac Surgery with Cardiopulmonary",
    "Bypass in Neonates with Critical Congenital Heart Disease.",
    "Clin Pharmacokinet. 2024;63:1205-1220.",
    "doi:10.1007/s40262-024-01401-3",
    sep = " "
  )
  vignette <- "Chu_2024_allopurinol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the ESM S3 NONMEM control stream
  # ($MODEL COMP=(ALL1)/(OXY)/(ALL2) with S1 = V1, S2 = VO*(136.11/152.11),
  # S3 = VA), which fixes both the analyte and the mass units of each state.
  compartmentData <- list(
    central     = list(analyte = "allopurinol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "allopurinol", units = "mg", specimen = "plasma", verified = TRUE),
    central_oxy = list(analyte = "oxypurinol",  units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor on every clearance and every volume, normalised to 3.5 kg with fixed exponents 0.75 (clearances) and 1 (volumes) per Chu 2024 Sect. 2.3 and the ESM S3 control stream (`(WT/3.5)**0.75` / `(WT/3.5)**1`). The 3.5 kg reference was chosen to permit comparison with the authors' earlier models in neonates with hypoxic-ischemic encephalopathy, not because it is the cohort median: Chu 2024 Table 1 reports a median BIRTH weight of 3.16 kg (IQR 2.75-3.73). Time-varying in the source dataset (ESM S3 $INPUT `WT ; Weight (kg)`).",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Chu 2024 reports postnatal age in DAYS (TM50 = 4.2 days); the canonical PNA column carries months, so `model()` recovers days with `PNA * 30.4375` before forming the sigmoidal recovery term (the same reparameterisation used by Zhao_2018_omeprazole.R and Bardhi_2026_ampicillin_foal.R). Drives the recovery of BOTH allopurinol and oxypurinol clearance during the postnatal-preoperative period ONLY -- during CPB and after CPB the clearances are fixed multiples of the at-birth value and the recovery term drops out entirely (ESM S3: `TVCLA = (CLA_PNA**FLAG1) * (CLA_CPB**FLAG2) * (CLA_POST**FLAG3)` with `CLA_CPB = CLA_PRE * FA_CPB1`, i.e. the CPB and post-CPB fractions multiply the PNA-free baseline). Cohort postnatal age at start of surgery: median 5.60 days, IQR 4.78-7.81 (Chu 2024 Table 1).",
      source_name        = "PNA"
    ),
    CPB_ON = list(
      description        = "Cardiopulmonary bypass phase indicator (on bypass, before rewarming begins)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (postnatal-preoperative period, when CPB_POST is also 0)",
      notes              = "Time-varying. Chu 2024 did NOT separate a rewarming phase: its 'intraoperative period' is the whole bypass run, delimited by the start and end times of CPB (Sect. 2.2). This model therefore consumes the on-bypass window as the sum `CPB_ON + CPB_REWARM`, which is the idiom the covariate register prescribes for an effect spanning both sub-windows; widening CPB_ON itself is explicitly forbidden there. Because only the sum enters, it does not matter how a user's data splits the run -- setting CPB_ON = 1 for the entire bypass run with CPB_REWARM = 0 gives the same predictions as an accurate split. Corresponds to ESM S3 `FLAG2` (`IF(OOC.EQ.3) FLAG2=1`).",
      source_name        = "OCC == 3"
    ),
    CPB_REWARM = list(
      description        = "Cardiopulmonary bypass rewarming phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not in the rewarming phase)",
      notes              = "Time-varying and mutually exclusive with CPB_ON. Carried only so that the on-bypass window can be written as `CPB_ON + CPB_REWARM` without widening CPB_ON. Chu 2024 fitted no rewarming-specific effect, so this column has no independent effect in this model and may be left at 0 throughout when the rewarming boundary is unknown (see the CPB_ON notes). The cohort was cooled to a median lowest rectal temperature of 27.7 C during CPB (Chu 2024 Table 1), so rewarming did occur; it simply was not modelled separately.",
      source_name        = "OCC == 3"
    ),
    CPB_POST = list(
      description        = "Post-cardiopulmonary-bypass (postoperative) period indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (postnatal-preoperative period, when CPB_ON and CPB_REWARM are also 0)",
      notes              = "Time-varying; 1 from separation from the bypass circuit onward. Chu 2024 estimated a distinct postoperative clearance and volume for both analytes rather than letting the postoperative phase collapse onto the pre-CPB reference, which is the condition under which the covariate register directs a sibling CPB_POST to be registered. Corresponds to ESM S3 `FLAG3` (`IF(POC.EQ.4) FLAG3=1`). The postoperative oxypurinol clearance (0.05 L/h per 3.5 kg) is LOWER than the at-birth value, which the authors attribute to CPB-induced acute kidney injury (Chu 2024 Sect. 4).",
      source_name        = "OCC == 4"
    )
  )

  # Columns present in the ESM S3 $INPUT record but carrying no effect in the
  # published final model. Chu 2024 Sect. 2.3 reports a covariate analysis on
  # postnatal age only (body size entering as fixed-exponent allometry), and
  # the ESM S2 model-development table lists no screening step for these two.
  covariatesDataExcluded <- list(
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Present in the ESM S3 $INPUT record as `GAW ; Gestational age (week)` and summarised in Chu 2024 Table 1 (median 38.0 weeks, IQR 38.0-38.8), but no covariate effect of gestational age is reported and none is present in the published $PK block. The cohort is near-uniformly term, so the data would carry little information about prematurity in any case."
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Present in the ESM S3 $INPUT record as `SEX ; Sex` and summarised in Chu 2024 Table 1 (10 of 14, 71.4%, male; hence 28.6% female), but no covariate effect of sex is reported and none is present in the published $PK block."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 14L,
    n_studies      = 1L,
    age_range      = "Neonates. Gestational age at birth median 38.0 weeks (IQR 38.0-38.8). Postnatal age at the start of cardiac surgery median 5.60 days (IQR 4.78-7.81); dosing and sampling span birth to roughly 1.5 weeks of life.",
    weight_range   = "Birth weight median 3.16 kg (IQR 2.75-3.73). Model parameters are reported at a 3.5 kg reference weight.",
    sex_female_pct = 28.6,
    race_ethnicity = NA_character_,
    disease_state  = "Critical congenital heart disease (CCHD) requiring cardiac surgery with cardiopulmonary bypass within the first month of life. Cardiac pathology: transposition of the great arteries 6 (42.9%), single ventricle physiology 4 (28.6%), aortic arch anomaly 2 (14.3%), other 2 (14.3%). Median total duration of cardiac surgery with CPB 320 min (IQR 280-368); median lowest rectal temperature during CPB 27.7 C (IQR 23.7-28); deep hypothermic cardiac arrest in 1 of 13 (7.7%) and antegrade cerebral perfusion in 4 of 13 (30.8%).",
    dose_range     = "Five intravenous allopurinol doses of 20 mg/kg each, delivered over 10 min by syringe pump: within 45-60 min after birth (DOSE 1), 12 h later (DOSE 2), 12 h before cardiac surgery (DOSE 3), at the start of CPB (DOSE 4), and 24 h after surgery (DOSE 5). The 10-min infusion duration is from the CRUCIAL study protocol (Stegeman 2022 Trials, doi:10.1186/s13063-022-06098-y), which Chu 2024 cites as reference 6 and as the source of its Fig. 1 dosing schedule; the PK paper itself does not restate it.",
    regions        = "The Netherlands. The CRUCIAL trial runs in four Dutch academic centres; this PK substudy was performed exclusively at Wilhelmina Children's Hospital, University Medical Center Utrecht.",
    n_observations = "140 allopurinol and oxypurinol plasma observations: 60 postnatal-preoperative, 36 intraoperative, 44 postoperative; median 10 samples per patient (IQR 9-12). 5.8% of allopurinol concentrations were below the LLOQ (0.05 mg/L allopurinol, 0.0467 mg/L oxypurinol) and were imputed at LLOQ/2 with an additive residual component fixed to LLOQ/2.",
    notes          = "PK substudy of the CRUCIAL trial (ClinicalTrials.gov NCT04217421; EudraCT 2017-004596-31), a phase III randomised quadruple-blinded placebo-controlled multicentre trial of postnatal and perioperative allopurinol for postoperative brain injury in neonates with CCHD. One of the 14 neonates had cardiac surgery at a non-participating centre and contributed postnatal-period samples only, so the 13 surgical patients underpin the intraoperative and postoperative parameters. Estimation by NONMEM 7.5 FOCE-I with ADVAN13; parameter precision by sampling importance resampling."
  )

  ini({
    # ---------------------------------------------------------------------
    # Allopurinol (parent). Chu 2024 Table 2 reports these at the 3.5 kg
    # reference weight and, for clearance, at birth (PNA = 0, where the
    # sigmoidal recovery term of Eq. 6 equals 1).
    #
    # NOTE on which volume is which: the ESM S3 control stream sets
    # S3 = VA and K31 = Q1/VA, so the estimated "All Vd_Postnatal" (2.22 L)
    # is the PERIPHERAL volume and the fixed "All V1" (0.1 L) is the
    # observed CENTRAL volume. Total Vss = 0.1 + 2.22 = 2.32 L.
    # ---------------------------------------------------------------------
    lcl <- log(0.95);          label("Allopurinol clearance at birth (CL_Postnatal, L/h per 3.5 kg)")            # Chu 2024 Table 2: All CL_Postnatal = 0.95 L/h (95% CI 0.75-1.2)
    lvp <- log(2.22);          label("Allopurinol peripheral volume of distribution (Vd_Postnatal, L per 3.5 kg)")  # Chu 2024 Table 2: All Vd_Postnatal = 2.22 L (95% CI 2-2.49); ESM S3 S3 = VA
    lvc <- fixed(log(0.1));    label("Allopurinol central volume of distribution (V1, L per 3.5 kg)")             # Chu 2024 Table 2: All V1 = 0.1 L (fix); ESM S3 S1 = V1
    lq  <- fixed(log(6.97));   label("Allopurinol intercompartmental clearance (Q1, L/h per 3.5 kg)")             # Chu 2024 Table 2: All Q1 = 6.97 L/h (fix). Sect. 3.2 text prints 6.79 L/h; see vignette Errata.

    # ---------------------------------------------------------------------
    # Oxypurinol (metabolite). Chu 2024 Table 2 footnote b: "Oxypurinol
    # clearances and volume of distributions are relative to formation
    # fraction (assumed to be 1)".
    # ---------------------------------------------------------------------
    lcl_oxy <- log(0.21);      label("Oxypurinol clearance at birth (CL_Postnatal, L/h per 3.5 kg)")              # Chu 2024 Table 2: Oxy CL_Postnatal = 0.21 L/h (95% CI 0.17-0.27)
    lvc_oxy <- log(12);        label("Oxypurinol volume of distribution (Vd_Postnatal, L per 3.5 kg)")            # Chu 2024 Table 2: Oxy Vd_Postnatal = 12 L (95% CI 9.87-15.33)

    # ---------------------------------------------------------------------
    # Auto-inhibition of allopurinol metabolism by oxypurinol. ESM S3 $DES:
    #   AUTOI = EMAX*C2/(IC50+C2);  DADT(1) = -K12*A(1)*(1-AUTOI) + ...
    # Both parameters were fixed (the IC50 was estimated on the PNA <= 3 day
    # data and then held) -- Chu 2024 Sect. 3.2 and Table 2.
    # ---------------------------------------------------------------------
    limax <- fixed(log(1));    label("Maximum achievable auto-inhibition of allopurinol conversion (fraction of conversion clearance)")  # Chu 2024 Table 2: Maximum achievable autoinhibition effect = 1 (fix)
    lic50 <- fixed(log(1.1));  label("Oxypurinol concentration at 50% of maximum auto-inhibition (IC50, mg/L)")                          # Chu 2024 Table 2: IC50,auto-inhibition = 1.1 mg/L (fix)
    fm    <- fixed(1);         label("Molar fraction of allopurinol converted to oxypurinol (unitless)")                                 # Chu 2024 Table 2 footnote b and ESM S1 legend: "fm, formation fraction was assumed 1"

    # ---------------------------------------------------------------------
    # Allometric exponents. Fixed, not estimated; a single exponent is shared
    # across every clearance (CL, Q1, CL_oxy) and every volume (V1, Vd,
    # Vd_oxy) -- ESM S3 applies `(WT/3.5)**0.75` and `(WT/3.5)**1` to each in
    # turn.
    # ---------------------------------------------------------------------
    e_wt_cl <- fixed(0.75);    label("Allometric (WT) exponent shared by CL, Q1 and CL_oxy versus WT / 3.5 kg")   # Chu 2024 Sect. 2.3: "fixed power exponents of 0.75 and 1 for CL and Vd"
    e_wt_vc <- fixed(1);       label("Allometric (WT) exponent shared by V1, Vd and Vd_oxy versus WT / 3.5 kg")   # Chu 2024 Sect. 2.3: "fixed power exponents of 0.75 and 1 for CL and Vd"

    # ---------------------------------------------------------------------
    # Recovery of clearance with postnatal age during the postnatal-
    # preoperative period. Chu 2024 Eq. 6:
    #   Recovery = 1 + Emax * PNA^Hill / (TM50^Hill + PNA^Hill)
    # All four were fixed after the covariate step (Chu 2024 Sect. 3.2: "The
    # effect of recovery was fixed to the identified values for further model
    # development"). TM50 and Hill are shared by the two analytes (ESM S3
    # THETA(9), THETA(10) feed both EFFPNA1 and EFFPNA2).
    # ---------------------------------------------------------------------
    e_pna_cl     <- fixed(3);     label("Maximum fold increase in allopurinol CL_Postnatal with postnatal age (unitless)")   # Chu 2024 Table 2: Maximum fold of increase in All CL_Postnatal = 3 (fixed); ESM S3 MMAX1 = THETA(7)
    e_pna_cl_oxy <- fixed(1.35);  label("Maximum fold increase in oxypurinol CL_Postnatal with postnatal age (unitless)")    # Chu 2024 Table 2: Maximum fold of increase in Oxy CL_Postnatal = 1.35 (fixed); ESM S3 MMAX2 = THETA(8)
    tm50_cl      <- fixed(4.2);   label("Postnatal age at 50% of the maximum recovery effect on clearance (days)")           # Chu 2024 Table 2: Postnatal age at 50% of maximum recovery effect = 4.2 days (fixed)
    hill_cl      <- fixed(2.98);  label("Hill coefficient of the postnatal-age recovery effect on clearance (unitless)")     # Chu 2024 Table 2: Hill coefficient = 2.98 (fixed)

    # ---------------------------------------------------------------------
    # Fractional change during cardiopulmonary bypass (Chu 2024 Eqs. 1-2).
    # These multiply the AT-BIRTH baseline, not the postnatal-age-adjusted
    # value: 0.95 * 1.46 = 1.39 L/h reproduces the 1.38 L/h quoted in
    # Chu 2024 Sect. 3.2 and Fig. 3, whereas 2.97 * 1.46 would not.
    # ---------------------------------------------------------------------
    e_cpb_on_cl     <- 1.46;  label("Fractional change in allopurinol CL during CPB versus CL_Postnatal at birth (unitless)")       # Chu 2024 Table 2: All E_CL,CPB = 1.46 (95% CI 0.96-2.13)
    e_cpb_on_vp     <- 1.47;  label("Fractional change in allopurinol Vd during CPB versus Vd_Postnatal (unitless)")                # Chu 2024 Table 2: All E_Vd,CPB = 1.47 (95% CI 1.27-1.72)
    e_cpb_on_cl_oxy <- 0.54;  label("Fractional change in oxypurinol CL during CPB versus CL_Postnatal at birth (unitless)")        # Chu 2024 Table 2: Oxy E_CL,CPB = 0.54 (95% CI 0.21-1.02)
    e_cpb_on_vc_oxy <- 1.3;   label("Fractional change in oxypurinol Vd during CPB versus Vd_Postnatal (unitless)")                 # Chu 2024 Table 2: Oxy E_Vd,CPB = 1.3 (95% CI 1.2-1.43)

    # ---------------------------------------------------------------------
    # Fractional change after cardiopulmonary bypass (Chu 2024 Eqs. 3-4).
    # Same baseline convention: 0.95 * 2.33 = 2.21 L/h and
    # 0.21 * 0.23 = 0.048 L/h reproduce the 2.21 and 0.05 L/h of Fig. 3.
    # ---------------------------------------------------------------------
    e_cpb_post_cl     <- 2.33;  label("Fractional change in allopurinol CL after CPB versus CL_Postnatal at birth (unitless)")      # Chu 2024 Table 2: All E_CL,Postop = 2.33 (95% CI 1.95-2.84)
    e_cpb_post_vp     <- 1.42;  label("Fractional change in allopurinol Vd after CPB versus Vd_Postnatal (unitless)")               # Chu 2024 Table 2: All E_Vd,Postop = 1.42 (95% CI 1.19-1.69)
    e_cpb_post_cl_oxy <- 0.23;  label("Fractional change in oxypurinol CL after CPB versus CL_Postnatal at birth (unitless)")       # Chu 2024 Table 2: Oxy E_CL,Postop = 0.23 (95% CI 0.06-0.42)
    e_cpb_post_vc_oxy <- 1.48;  label("Fractional change in oxypurinol Vd after CPB versus Vd_Postnatal (unitless)")                # Chu 2024 Table 2: Oxy E_Vd,Postop = 1.48 (95% CI 1.34-1.63)

    # ---------------------------------------------------------------------
    # Between-subject variability. Chu 2024 Table 2 footnote states the
    # conversion explicitly: "CV coefficient of variation, approximated using
    # CV% = sqrt(omega^2) * 100" -- so the tabulated percentages ARE the
    # omega SDs and omega^2 = (CV/100)^2. This is NOT the log-normal
    # moment-match omega^2 = log(CV^2 + 1).
    # The ESM S3 control stream also declares ETA(2) on the allopurinol
    # volume and ETA(3) on the oxypurinol clearance, but Chu 2024 Table 2
    # reports no estimate for either, so neither is carried here.
    # ---------------------------------------------------------------------
    etalcl     ~ 0.1296  # Chu 2024 Table 2: BSV All CL = 36% CV (95% CI 24-50); omega^2 = 0.36^2
    etalvc_oxy ~ 0.1764  # Chu 2024 Table 2: BSV Oxy Vd = 42% CV (95% CI 30-54); omega^2 = 0.42^2

    # ---------------------------------------------------------------------
    # Between-occasion variability across the three perioperative periods.
    # ESM S3: BOV_CLA = FLAG1*ETA(5) + FLAG2*ETA(6) + FLAG3*ETA(7) and
    # BOV_CLO = FLAG1*ETA(8) + FLAG2*ETA(9) + FLAG3*ETA(10), i.e. one
    # occasion per period with a SHARED variance (NONMEM $OMEGA BLOCK SAME).
    # rxode2 has no `| occ` level, so the occasions are expanded into three
    # indicator-multiplexed etas each; occasion 1 carries the estimated
    # variance and occasions 2-3 are fixed to the same value.
    # ---------------------------------------------------------------------
    etaiov_cl_1 ~ 0.0324         # Chu 2024 Table 2: BOV All CL = 18% CV (95% CI 12-26); omega^2 = 0.18^2
    etaiov_cl_2 ~ fixed(0.0324)  # Chu 2024 Table 2: BOV All CL = 18% CV; shared variance
    etaiov_cl_3 ~ fixed(0.0324)  # Chu 2024 Table 2: BOV All CL = 18% CV; shared variance

    etaiov_cl_oxy_1 ~ 0.1849         # Chu 2024 Table 2: BOV Oxy CL = 43% CV (95% CI 30-59); omega^2 = 0.43^2
    etaiov_cl_oxy_2 ~ fixed(0.1849)  # Chu 2024 Table 2: BOV Oxy CL = 43% CV; shared variance
    etaiov_cl_oxy_3 ~ fixed(0.1849)  # Chu 2024 Table 2: BOV Oxy CL = 43% CV; shared variance

    # ---------------------------------------------------------------------
    # Residual error. ESM S3 $ERROR is combined proportional + additive per
    # analyte: Y = F*(1+EPS(1)) + EPS(2) for allopurinol and
    # Y = F*(1+EPS(3)) + EPS(4) for oxypurinol. Chu 2024 Sect. 3.1 states the
    # additive component was fixed to LLOQ/2 for the below-LLOQ imputation;
    # Table 2 tabulates only the proportional components.
    # ---------------------------------------------------------------------
    propSd     <- 0.25;            label("Allopurinol proportional residual SD (fraction)")   # Chu 2024 Table 2: Residual proportional error for allopurinol = 25% CV (95% CI 21-29)
    addSd      <- fixed(0.025);    label("Allopurinol additive residual SD (mg/L)")           # Chu 2024 Sect. 3.1: "an additive residual error component of LLOQ/2 fixed to the error model"; LLOQ_allopurinol = 0.05 mg/L (Sect. 2.1)
    propSd_oxy <- 0.16;            label("Oxypurinol proportional residual SD (fraction)")    # Chu 2024 Table 2: Residual proportional error for oxypurinol = 16% CV (95% CI 14-18)
    addSd_oxy  <- fixed(0.02335);  label("Oxypurinol additive residual SD (mg/L; value assumed, see vignette Errata)")  # ESM S3 declares EPS(4) but no value is published; set to LLOQ_oxypurinol/2 = 0.0467/2 following the stated allopurinol convention
  })

  model({
    # Molecular weights (g/mol) used to convert the mass flux of allopurinol
    # leaving its central compartment (mg/h) into the mass flux of oxypurinol
    # entering the metabolite compartment (mg/h).
    mwAllo <- 136.11
    mwOxy  <- 152.11

    # The ESM S3 control stream keeps A(2) in ALLOPURINOL-equivalent mass and
    # pushes the molar conversion into the observation scaling instead
    # (S2 = VO * 136.11/152.11). Applying (mwOxy/mwAllo) at the transfer, as
    # done here, is algebraically identical for both Cc_oxy and the
    # oxypurinol elimination rate, and has the advantage that `central_oxy`
    # genuinely holds mg of oxypurinol (as compartmentData asserts).

    # Postnatal age recovered in DAYS from the canonical PNA column (months).
    pna_d <- PNA * 30.4375

    # Perioperative period indicators. Chu 2024 splits the record into three
    # mutually exclusive periods at the start and end of CPB. The on-bypass
    # window is the SUM of the two bypass sub-phase canonicals because the
    # paper fitted no rewarming-specific effect (see covariateData notes).
    on_cpb   <- CPB_ON + CPB_REWARM
    post_cpb <- CPB_POST
    pre_cpb  <- 1 - on_cpb - post_cpb

    # Sigmoidal recovery of clearance with postnatal age (Chu 2024 Eq. 6).
    # Applies during the postnatal-preoperative period only.
    pna_frac        <- pna_d^hill_cl / (tm50_cl^hill_cl + pna_d^hill_cl)
    recovery_cl     <- 1 + e_pna_cl     * pna_frac
    recovery_cl_oxy <- 1 + e_pna_cl_oxy * pna_frac

    # Allometric size scaling to the 3.5 kg reference weight.
    wt_cl <- (WT / 3.5)^e_wt_cl
    wt_v  <- (WT / 3.5)^e_wt_vc

    # Between-occasion variability, one occasion per perioperative period.
    iov_cl     <- pre_cpb * etaiov_cl_1     + on_cpb * etaiov_cl_2     + post_cpb * etaiov_cl_3
    iov_cl_oxy <- pre_cpb * etaiov_cl_oxy_1 + on_cpb * etaiov_cl_oxy_2 + post_cpb * etaiov_cl_oxy_3

    # At-birth (postnatal-preoperative) baselines, before the period effects.
    cl_base     <- exp(lcl     + etalcl + iov_cl) * wt_cl
    cl_oxy_base <- exp(lcl_oxy + iov_cl_oxy)      * wt_cl
    vp_base     <- exp(lvp)                       * wt_v
    vc_oxy_base <- exp(lvc_oxy + etalvc_oxy)      * wt_v

    # Period-specific disposition (Chu 2024 Eqs. 1-4). Exactly one indicator
    # is 1, so this sum selects a single branch; it is the additive form of
    # the control stream's `(X_PRE**FLAG1) * (X_CPB**FLAG2) * (X_POST**FLAG3)`.
    cl <- pre_cpb  * cl_base * recovery_cl +
          on_cpb   * cl_base * e_cpb_on_cl +
          post_cpb * cl_base * e_cpb_post_cl
    vp <- pre_cpb  * vp_base +
          on_cpb   * vp_base * e_cpb_on_vp +
          post_cpb * vp_base * e_cpb_post_vp
    cl_oxy <- pre_cpb  * cl_oxy_base * recovery_cl_oxy +
              on_cpb   * cl_oxy_base * e_cpb_on_cl_oxy +
              post_cpb * cl_oxy_base * e_cpb_post_cl_oxy
    vc_oxy <- pre_cpb  * vc_oxy_base +
              on_cpb   * vc_oxy_base * e_cpb_on_vc_oxy +
              post_cpb * vc_oxy_base * e_cpb_post_vc_oxy

    # V1 and Q1 are period-invariant; only the allometric term applies.
    vc <- exp(lvc) * wt_v
    q  <- exp(lq)  * wt_cl

    imax <- exp(limax)
    ic50 <- exp(lic50)

    kel     <- cl     / vc
    k12     <- q      / vc
    k21     <- q      / vp
    kel_oxy <- cl_oxy / vc_oxy

    # Auto-inhibition of the allopurinol -> oxypurinol conversion, driven by
    # the oxypurinol concentration (ESM S3 $DES). With imax fixed to 1 the
    # remaining fraction (1 - autoi) collapses to ic50 / (ic50 + Cc_oxy).
    Cc_oxy <- central_oxy / vc_oxy
    autoi  <- imax * Cc_oxy / (ic50 + Cc_oxy)

    # ODE system. Allopurinol elimination IS the conversion to oxypurinol
    # (formation fraction fixed to 1), so the auto-inhibition term gates both
    # the loss of allopurinol and the appearance of oxypurinol.
    d/dt(central)     <- -kel * central * (1 - autoi) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_oxy) <-  fm * kel * central * (1 - autoi) * (mwOxy / mwAllo) -
                          kel_oxy * central_oxy

    # Plasma concentrations in mg/L (amounts in mg, volumes in L).
    Cc <- central / vc

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_oxy ~ add(addSd_oxy) + prop(propSd_oxy)
  })
}

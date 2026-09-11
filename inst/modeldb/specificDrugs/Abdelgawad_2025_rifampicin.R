Abdelgawad_2025_rifampicin <- function() {
  description <- paste(
    "Semi-mechanistic two-compartment population PK model for rifampicin in",
    "plasma and lumbar cerebrospinal fluid (CSF) in adults with",
    "HIV-associated tuberculous meningitis given standard-dose (10 mg/kg",
    "oral), high-dose (35 mg/kg oral) or intravenous (20 mg/kg) rifampicin",
    "in the LASER-TBM trial (Abdelgawad 2025). Oral absorption is a Savic",
    "analytical transit chain (19 transit compartments fixed, mean transit",
    "time 0.634 h) feeding a first-order absorption compartment (ka 0.486",
    "1/h) that empties into a liver compartment, so first-pass extraction is",
    "structural; prehepatic oral bioavailability is 0.934 and intravenous",
    "bioavailability is fixed to 1 with a modelled 1 h infusion duration.",
    "Elimination is a well-stirred liver with saturable intrinsic clearance,",
    "CLint = CLint,max * Km / (CH + Km) and EH = CLint * fu / (CLint * fu +",
    "QH), with liver volume 1 L, hepatic blood flow 90 L/h and fraction",
    "unbound 0.2 all fixed. The maximal intrinsic clearance is estimated as",
    "four separate typical values rather than by an autoinduction model,",
    "one per dose group and PK visit: the reported CLint,max * fu products",
    "are 33.1 L/h (standard dose, day 3), 41.4 L/h (standard dose, day 28),",
    "46.1 L/h (high dose, day 3) and 70.2 L/h (high dose, day 28).",
    "Disposition parameters are allometrically scaled on fat-free mass with",
    "fixed 0.75 / 1 exponents, referenced to 46 kg for CLint,max, Q, V and",
    "Vp and to 56.1 kg for the fixed hepatic physiology. CSF is a",
    "Sheiner-style effect compartment holding a concentration, equilibrating",
    "with plasma at a 3.20 h half-life toward a pseudo-partition coefficient",
    "of 0.0593. Random effects are between-subject variability on CLint,max",
    "(25.3%), central volume (17.2%) and infusion duration (17.0%), and",
    "five-occasion between-occasion variability on prehepatic",
    "bioavailability (18.2%), ka (78.1%) and mean transit time (111%); the",
    "reported percentages are omega standard deviations on the log scale.",
    "Residual error is combined proportional plus additive, separately for",
    "plasma (25.2%, 0.0234 mg/L) and CSF (98.4%, 0.0231 mg/L)."
  )
  reference <- paste(
    "Abdelgawad N, Wasserman S, Gausi K, Davis A, Stek C, Wiesner L,",
    "Meintjes G, Wilkinson RJ, Denti P (2025).",
    "Population Pharmacokinetics of Rifampicin in Plasma and Cerebrospinal",
    "Fluid in Adults With Tuberculosis Meningitis.",
    "J Infect Dis 232(4):e234-e241. doi:10.1093/infdis/jiaf178.",
    "Parameter estimates from Table 2; model equations from the Figure 1",
    "caption and from the NONMEM control stream reproduced verbatim in the",
    "supplementary material.",
    "The saturable-hepatic-extraction structure was adapted from",
    "Chirehwa et al. (2016) Antimicrob Agents Chemother 60(1):487-494",
    "doi:10.1128/AAC.01084-15, which also supplied the informative prior on",
    "the Michaelis-Menten constant.",
    "The CSF effect compartment follows Sheiner et al. (1979)",
    "Clin Pharmacol Ther 25(3):358-371 and Savic et al. (2015)",
    "Clin Pharmacol Ther 98(6):622-629 doi:10.1002/cpt.202.",
    "Fat-free mass follows Janmahasatian et al. (2005)",
    "Clin Pharmacokinet 44(10):1051-1065 doi:10.2165/00003088-200544100-00004.",
    sep = " "
  )
  vignette <- "Abdelgawad_2025_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The `csf` state is an exception to the usual "states hold an amount" rule:
  # the source control stream's $DES integrates a concentration directly,
  # DADT(5) = KE0*(PPC*CP_DES - A(5)) with CP_DES = A(2)/V, and $ERROR reads
  # the CSF prediction as CE = A(5) with no volume division. Its units are
  # therefore mg/L, not mg. This matches the sibling LASER-TBM linezolid model
  # `Abdelgawad_2024_linezolid.R` (issue #482).
  compartmentData <- list(
    depot = list(
      analyte = "rifampicin", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    liver = list(
      analyte = "rifampicin", units = "mg",
      specimen = "tissue", verified = TRUE
    ),
    central = list(
      analyte = "rifampicin", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "rifampicin", units = "mg",
      specimen = "tissue", verified = TRUE
    ),
    csf = list(
      analyte = "rifampicin", units = "mg/L",
      specimen = "CSF", verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass, computed from sex, total body weight and height by the Janmahasatian (2005) formula",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor for every disposition parameter, with",
        "exponents fixed a priori at 0.75 for the clearance terms and 1 for",
        "the volume terms (Methods, PK Modeling: 'Allometric scaling was",
        "applied for all disposition parameters via the fixed power exponents",
        "of 0.75 for clearance parameters and 1 for volume parameters').",
        "Two reference values are used in the same model, exactly as in the",
        "control stream $PK block: CLint,max, Q, V and Vp are normalised to",
        "46 kg (control stream 'TVFFM = 46'; ALLMCL_FFM = (FFM/TVFFM)**0.75,",
        "ALLMV_FFM = (FFM/TVFFM)) while the fixed hepatic physiology QH and",
        "VH are normalised to 56.1 kg, the fat-free mass of the reference",
        "adult from which the 90 L/h and 1 L values are taken",
        "(ALLMCL_H_FFM = (FFM/56.1)**0.75, ALLMV_H_FFM = (FFM/56.1)).",
        "Table 2 footnote c describes the typical participant as having a",
        "median fat-free mass of 45 kg, matching the Table 1 cohort medians",
        "of 45.2 kg (day 3) and 45.5 kg (day 28); the control stream that",
        "produced the estimates normalises to 46 kg, so 46 kg is used here",
        "and the discrepancy is recorded in the vignette Errata.",
        "Allometry on fat-free mass fitted better than allometry on total",
        "body weight (dOFV 33.1 points, P < .0001 for FFM versus 13.7 points,",
        "P < .001 for total body weight; Results, PK Modeling).",
        "Height was missing for 29 of 49 participants at the day-3 visit and",
        "19 of 34 at the day-28 visit (Table 1 footnote a), so the authors",
        "imputed it inside NONMEM by the Johansson and Karlsson (2013)",
        "multiple-imputation approach before computing FFM: supplementary",
        "'Imputation of missing covariates' and control stream $PK give",
        "HT = (0.00133 * WT + 1.51) * exp(eta) for females and",
        "HT = (0.00281 * WT + 1.53) * exp(eta) for males with fixed eta",
        "variances 0.00215 and 0.00170, then Janmahasatian",
        "FFM = 37.99 * HT^2 * WT / (35.98 * HT^2 + WT) for females and",
        "FFM = 42.92 * HT^2 * WT / (30.93 * HT^2 + WT) for males, HT in m and",
        "WT in kg. That imputation is a missing-data device for the original",
        "fit, not part of the structural model; users should supply a measured",
        "or Janmahasatian-derived FFM column directly.",
        "Cohort range 30.3-59.4 kg at the day-3 visit (Table 1)."
      ),
      source_name        = "FFM"
    ),
    DOSE_HIGH = list(
      description        = "1 = participant randomised to a high-dose rifampicin experimental arm (35 mg/kg orally, or 20 mg/kg intravenously for the first 3 days); 0 = participant in the control arm receiving the standard 10 mg/kg oral dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (standard-dose control arm, 10 mg/kg orally by World Health Organization weight bands)",
      notes              = paste(
        "Together with DAY28 this selects one of the four typical values of",
        "CLint,max * fu that the authors estimated in place of an",
        "autoinduction model. Control stream $PK:",
        "TVCL = THETA(14)*ALLMCL_FFM;",
        "IF (RIFHIGH.EQ.0.AND.PK_VISIT_2.EQ.3)  TVCL = THETA(15)*ALLMCL_FFM;",
        "IF (RIFHIGH.EQ.1.AND.PK_VISIT_2.EQ.28) TVCL = THETA(16)*ALLMCL_FFM;",
        "IF (RIFHIGH.EQ.0.AND.PK_VISIT_2.EQ.28) TVCL = THETA(17)*ALLMCL_FFM.",
        "The four thetas are not constrained to a multiplicative 2 x 2",
        "structure -- the high-dose / standard-dose ratio is 1.39 at day 3 and",
        "1.70 at day 28 -- so all four are carried as separate stratum-suffixed",
        "typical values rather than as a reference plus covariate offsets.",
        "The mg/kg threshold that maps to DOSE_HIGH = 1 in this study is",
        "35 mg/kg orally (20 mg/kg intravenously on days 1-3); after day 3 all",
        "experimental-arm participants continued oral 35 mg/kg to the end of",
        "the study (Methods, Parent Study and Interventions), so the indicator",
        "is time-fixed per participant. High-dose oral rifampicin was given as",
        "fixed-dose-combination tablets topped up with individual rifampicin",
        "tablets according to bespoke weight bands; the authors found no",
        "significant bioavailability difference between the two tablet types",
        "(Results, PK Modeling).",
        "Clearance was higher in participants receiving larger doses",
        "(Results, PK Modeling); the effect is confounded with the arm-level",
        "difference in prior rifampicin exposure and with the enzyme-inducing",
        "co-treatment, which the authors could not separate because",
        "observations from the uninduced state were unavailable (Discussion,",
        "limitations). Source column RIFHIGH."
      ),
      source_name        = "RIFHIGH"
    ),
    DAY28 = list(
      description        = "1 = the record belongs to the day-28 pharmacokinetic visit; 0 = the record belongs to the day-3 pharmacokinetic visit",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (day-3 PK visit, the first sampling visit after study enrolment)",
      notes              = paste(
        "Within-subject landmark indicator gating the step change in",
        "CLint,max * fu between the two PK visits, which stands in for",
        "rifampicin autoinduction. The authors first tried the exponential",
        "time-on-treatment autoinduction model of Chirehwa et al (2016) and",
        "an enzyme-turnover model in the style of Svensson et al (2018);",
        "neither converged, because participants had already been taking",
        "rifampicin for about a week before the first PK visit, so no",
        "uninduced-state observations were available and the actual dose and",
        "duration of pre-enrolment treatment were uncertain (Results, PK",
        "Modeling; Discussion, limitations). Separate typical values per",
        "visit were used instead.",
        "Sampling visits were day 3 (visit 1) and day 28 (visit 2) plus or",
        "minus 2 days after study enrolment; the median (range) time since the",
        "start of rifampicin-based treatment was 4 days (0-7) at the day-3",
        "visit and 30 days (26-38) at the day-28 visit (Table 1).",
        "Day-3 sampling was intensive (predose and 0.5, 1, 2, 3, 6, 8-10 and",
        "24 h postdose) and day-28 sampling sparse (predose, 2 and 4 h",
        "postdose), which the authors note limits the precision of the",
        "day-28 CLint,max * fu estimates (Discussion, limitations).",
        "Data assemblers derive DAY28 = as.integer(pk_visit_day >= 28).",
        "Source column PK_VISIT_2 (values 3 and 28)."
      ),
      source_name        = "PK_VISIT_2"
    ),
    OCC = list(
      description        = "Integer sampling-occasion index (1-5) used for the between-occasion random effects on prehepatic bioavailability, the absorption rate constant and the mean transit time",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Five occasions, matching the control stream $PK block's",
        "IF (OCC==1) ... IF (OCC==5) multiplexers over ETA(4)-ETA(8) for",
        "prehepatic bioavailability, ETA(9)-ETA(13) for ka and",
        "ETA(14)-ETA(18) for mean transit time, each backed by a single",
        "$OMEGA BLOCK(1) followed by four SAME repeats.",
        "The supplementary material defines an occasion as 'each dose and its",
        "following samples', so the dose given before a sampling visit",
        "together with the predose concentration is a different occasion from",
        "the dose administered during the visit and the concentrations that",
        "follow it; the control stream comment marks occasion 1 as the",
        "predose occasion of the day-3 visit. For simulation, set OCC to the",
        "occasion index of each dosing interval; a single-occasion simulation",
        "may use OCC = 1 throughout. Source column OCC."
      ),
      source_name        = "OCC"
    )
  )

  # Covariates the source paper screened but did not retain in the final model.
  # Documented here so the provenance of the covariate search survives without
  # declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as the allometric body-size descriptor and rejected in favour of fat-free mass (dOFV 13.7 points, P < .001 for total body weight versus 33.1 points, P < .0001 for FFM; Results, PK Modeling). Still required upstream of the model as an input to the Janmahasatian FFM formula and to the height-imputation regression. Cohort median 59.5 kg, range 30-107.2 kg at the day-3 visit (Table 1); the control stream records the dataset median as TVWT = 60."
    ),
    HT = list(
      description        = "Body height at the PK visit.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not a model covariate; an input to the Janmahasatian FFM formula. Missing for 29 of 49 participants at day 3 and 19 of 34 at day 28 and imputed inside NONMEM from sex and weight (Table 1 footnote a; supplementary 'Imputation of missing covariates'). Reported in metres in the source (median 1.60 m, range 1.48-1.80); the canonical column is in cm."
    ),
    SEXF = list(
      description        = "1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Not a model covariate in its own right, but required upstream as the switch between the sex-specific Janmahasatian FFM formulas and between the sex-specific height-imputation regressions (control stream $PK: IF (SEXF.EQ.0) selects the male coefficients). 27 of 49 participants (55.1%) at the day-3 visit were male (Table 1)."
    ),
    AGE = list(
      description        = "Age at enrolment.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. Cohort median 39 years, range 25-78 at the day-3 visit (Table 1)."
    ),
    CSF_TPRO = list(
      description        = "Total protein concentration in lumbar cerebrospinal fluid.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a marker of meningeal inflammation on the pseudo-partition coefficient and on the equilibration half-life; not retained. 'None of the covariates tested resulted in a statistically significant effect on the PPC or the equilibration half-life' (Results, PK Modeling), a null result the authors attribute to the small sample size and narrow range of CSF protein values (Discussion). The full screened list on those two parameters was CSF total protein, CSF albumin, CSF glucose, polymorphonuclear cells, lymphocytes and the Glasgow Coma Scale (Methods, PK Modeling); only CSF total protein is registered as a canonical column, so the other five are named here in prose rather than minted as canonical names for a screen that produced no retained effect. Cohort median 1.16 g/L, range 0.2-55 at the day-3 visit; missing for 17 participants at day 3 and 8 at day 28 (Table 1). The sibling LASER-TBM linezolid model Abdelgawad_2024_linezolid.R DID retain a CSF-protein effect on its PPC, so the null result here is drug-specific rather than a property of the cohort."
    ),
    CREAT = list(
      description        = "Serum creatinine.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on the plasma PK parameters and not retained: 'We did not find a statistically significant difference in bioavailability for the FDC and the individual top-up tablets or for biomarkers such as creatinine, aspartate aminotransferase, and alanine aminotransferase' (Results, PK Modeling). Values are not tabulated in the paper."
    ),
    AST = list(
      description        = "Aspartate aminotransferase.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on the plasma PK parameters and not retained (Results, PK Modeling). Values are not tabulated in the paper."
    ),
    ALT = list(
      description        = "Alanine aminotransferase.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on the plasma PK parameters and not retained (Results, PK Modeling). Values are not tabulated in the paper."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 48L,
    n_studies      = 1L,
    age_range      = "25-78 years (median 39 at the day-3 visit; 25-57, median 39 at the day-28 visit)",
    age_median     = "39 years",
    weight_range   = "30-107.2 kg (median 59.5 at the day-3 visit; 37.4-105.1, median 61.7 at the day-28 visit)",
    weight_median  = "59.5 kg",
    ffm_range      = "30.3-59.4 kg (median 45.2 at the day-3 visit); the control stream normalises allometry to 46 kg",
    sex_female_pct = 44.9,
    race_ethnicity = "Not reported; the cohort was enrolled at four public hospitals in South Africa",
    disease_state  = paste(
      "HIV-associated tuberculous meningitis (TBM). All participants were",
      "living with HIV: 14 of 49 (28.6%) had previously taken antiretroviral",
      "therapy, 20 (40.8%) were antiretroviral-naive and 15 (30.6%) were on",
      "treatment at the day-3 visit. Median CSF total protein 1.16 g/L,",
      "albumin 387 mg/L and glucose 3.05 mmol/L at the day-3 visit. All",
      "participants received adjunctive corticosteroids."
    ),
    dose_range     = paste(
      "Control arm: standard-of-care oral rifampicin 10 mg/kg once daily by",
      "World Health Organization weight bands, as fixed-dose-combination",
      "tablets with isoniazid 5 mg/kg, pyrazinamide 25 mg/kg and ethambutol",
      "15 mg/kg. Experimental arms: high-dose rifampicin plus oral linezolid",
      "1200 mg daily, with or without aspirin, randomised for the first 3",
      "days to either oral 35 mg/kg (fixed-dose-combination tablets topped up",
      "with individual rifampicin tablets by bespoke weight bands) or",
      "intravenous 20 mg/kg given as a 1 h infusion; from day 3 onward all",
      "experimental-arm participants took oral 35 mg/kg once daily."
    ),
    regions        = "South Africa (four hospitals)",
    notes          = paste(
      "Pharmacokinetic substudy of LASER-TBM, a phase 2A trial of intensified",
      "antibiotic therapy in adults with HIV-associated TBM",
      "(ClinicalTrials.gov NCT03927313). Forty-nine participants underwent PK",
      "sampling on day 3 and 34 on day 28, providing 411 plasma samples (56",
      "below the limit of quantification) and 46 CSF samples (13 below the",
      "limit of quantification); rifampicin concentrations from one",
      "participant were excluded after intravenous catheter dislocation and",
      "tissue extravasation (control stream IGNORE(ID.EQ.4019)), leaving the",
      "400 plasma and 44 CSF concentrations from 48 participants quoted in",
      "the abstract. Plasma sampling was predose and 0.5, 1, 2, 3, 6, 8-10",
      "and 24 h postdose on day 3 and predose and 2 and 4 h postdose on day",
      "28. One lumbar CSF sample was taken per visit, with the sampling time",
      "randomised across the 1-3, 3-6, 6-10 and 24 h postdose windows.",
      "Concentrations below the limit of quantification were handled by",
      "Beal's M6 method. Free plasma rifampicin was measured in a subset of",
      "participants; Deming regression through the origin gave a fraction",
      "unbound of 0.172, i.e. 82.8% plasma protein binding, with no trend of",
      "fraction unbound against total concentration (Results; supplementary",
      "Figure S2). Baseline characteristics are Table 1."
    )
  )

  ini({
    # =======================================================================
    # Structural parameters. Point estimates are Abdelgawad 2025 Table 2
    # (95% CIs by sampling importance resampling), cross-checked against the
    # $THETA block of the NONMEM control stream reproduced verbatim in the
    # supplementary material. The control stream's block is labelled
    # "initialization-of-theta(S)-from the previous run" but every one of its
    # values reproduces the corresponding Table 2 row to three significant
    # figures, so it is the converged final vector carried forward; where it
    # holds extra significant figures they are used here.
    # Time unit h, dose mg, concentration mg/L.
    # =======================================================================

    # Km was estimated on the log scale, so THETA(1) IS log(Km) and needs no
    # log() wrapper: exp(1.09075) = 2.977 mg/L, the 2.97 mg/L of Table 2.
    # Estimated under an informative prior taken from Chirehwa et al (2016)
    # ($THETAP 1.21 FIX with $THETAPV BLOCK(1) FIX 0.1, i.e. exp(1.21) =
    # 3.35 mg/L with sqrt(0.1) = 32% uncertainty on the log scale, matching
    # Table 2 footnote d "a prior value from Chirehwa et al of 3.35 mg/L with
    # 30% uncertainty").
    lkm <- 1.09075
    label("Michaelis-Menten constant Km, the hepatic rifampicin concentration at half of Vmax (mg/L, on the log scale)")  # Table 2 Km = 2.97 mg/L (2.01-4.56); control stream $THETA (0,1.09075,100) ; 1 LOGKM

    lvc <- log(27.3473)
    label("Central volume of distribution V at the reference fat-free mass of 46 kg (L)")                                 # Table 2 V = 27.3 L (22.8-34.0); control stream $THETA (0,27.3473,125) ; 2 V
    lvp <- log(31.5279)
    label("Peripheral volume of distribution Vp at the reference fat-free mass of 46 kg (L)")                             # Table 2 Vp = 31.5 L (25.3-37.1); control stream $THETA (0,31.5279,1000) ; 12 VP
    lq <- log(10.9971)
    label("Intercompartmental clearance Q at the reference fat-free mass of 46 kg (L/h)")                                 # Table 2 Q = 11.0 L/h (7.45-14.4); control stream $THETA (0,10.9971,1000) ; 11 Q
    lka <- log(0.485942)
    label("First-order absorption rate constant ka from the absorption compartment into the liver (1/h)")                 # Table 2 ka = 0.486 1/h (0.365-0.638); control stream $THETA (0,0.485942,3) ; 4 KA
    lmtt <- log(0.634202)
    label("Mean transit time MTT through the absorption transit chain (h)")                                               # Table 2 mean transit time = 0.634 h (0.467-0.791); control stream $THETA (0,0.634202,3) ; 5 MTT
    lntr <- fixed(log(19))
    label("Number of absorption transit compartments (unitless)")                                                         # Table 2 "No. of absorption transit compartments: 19 fixed" with footnote f; control stream $THETA (0,19,100) FIX ; 13 NN
    lfdepot <- log(0.934124)
    label("Prehepatic oral bioavailability, the fraction absorbed from the gastrointestinal tract before hepatic extraction (unitless)")  # Table 2 "Prehepatic oral 0.934 (0.852-0.991)" with footnote e; control stream $THETA (0,0.934124,1) ; 3 BIO_ORAL
    lfcentral <- fixed(log(1))
    label("Absolute intravenous bioavailability F for a dose given into the central compartment (unitless)")              # Table 2 "Intravenous, F: 1 fixed"; control stream $THETA 1 FIX ; 10 BIO_IV
    ldur <- fixed(log(1))
    label("Duration of the intravenous infusion into the central compartment (h)")                                        # Table 2 footnote g "The infusion duration is 1 hour according to the protocol"; control stream D2 = DUR*EXP(BSVD2) with the protocol DUR column

    # ----- Fixed hepatic physiology (Table 2 footnote a) -------------------
    lvh <- fixed(log(1))
    label("Liver volume VH at the reference fat-free mass of 56.1 kg (L)")                                                # Table 2 footnote a "Hepatic volume of distribution ... fixed to 1 L"; control stream $THETA 1 FIX ; 8 VH
    lqh <- fixed(log(90))
    label("Hepatic blood flow QH at the reference fat-free mass of 56.1 kg (L/h)")                                        # Table 2 footnote a "hepatic intercompartmental clearance ... fixed to ... 90 L/h"; control stream $THETA 90 FIX ; 9 QH
    fub <- fixed(0.2)
    label("Fraction of rifampicin unbound in blood, fu (unitless)")                                                       # Table 2 footnote a "fraction unbound ... fixed to ... 0.2"

    # ----- Maximal intrinsic clearance, one typical value per dose group ---
    # and PK visit. Table 2 reports the PRODUCT CLint,max * fu (its row label
    # is "CL int,max . f u , L/h") and footnote a fixes fu to 0.2, so
    # CLint,max = (printed product) / 0.2. The control stream carries the
    # product directly as its CL variable and never writes fu out
    # ($DES: SAT_CL = VMAX/(CH+EXP(LOGKM)); EH = SAT_CL/(SAT_CL+QH)); the
    # equivalent well-stirred form used in model() below multiplies by fub
    # explicitly, matching the sibling rifampicin model
    # Gafar_2026_rifampicin.R. The two forms are algebraically identical:
    # (CLint,max / fu) * fu = CLint,max * fu.
    #
    # Four stratum-suffixed typical values rather than a reference plus
    # covariate offsets, because the four thetas are unconstrained: the
    # high-dose / standard-dose ratio is 46.1/33.1 = 1.39 at day 3 but
    # 70.2/41.4 = 1.70 at day 28, so no 2 x 2 multiplicative form reproduces
    # all four.
    lclint_max_std_d3 <- log(33.1219 / 0.2)
    label("Maximal intrinsic hepatic clearance CLint,max, standard dose at the day-3 visit, at the reference fat-free mass of 46 kg (L/h)")   # Table 2 CLint,max*fu standard dose visit 1 = 33.1 L/h (25.2-42.9), / fu 0.2 = 165.6 L/h; control stream $THETA (0,33.1219,1000) ; 15 CL_Day3_Std
    lclint_max_std_d28 <- log(41.411 / 0.2)
    label("Maximal intrinsic hepatic clearance CLint,max, standard dose at the day-28 visit, at the reference fat-free mass of 46 kg (L/h)")  # Table 2 CLint,max*fu standard dose visit 2 = 41.4 L/h (28.0-58.3), / fu 0.2 = 207.1 L/h; control stream $THETA (0,41.411,1000) ; 17 CL_Day28_Std
    lclint_max_high_d3 <- log(46.1064 / 0.2)
    label("Maximal intrinsic hepatic clearance CLint,max, high dose at the day-3 visit, at the reference fat-free mass of 46 kg (L/h)")       # Table 2 CLint,max*fu high dose visit 1 = 46.1 L/h (34.2-62.3), / fu 0.2 = 230.5 L/h; control stream $THETA (0,46.1064,1000) ; 14 CL_Day3_High
    lclint_max_high_d28 <- log(70.2333 / 0.2)
    label("Maximal intrinsic hepatic clearance CLint,max, high dose at the day-28 visit, at the reference fat-free mass of 46 kg (L/h)")      # Table 2 CLint,max*fu high dose visit 2 = 70.2 L/h (51.0-95.5), / fu 0.2 = 351.2 L/h; control stream $THETA (0,70.2333,1000) ; 16 CL_Day28_High

    # ----- CSF effect compartment -----------------------------------------
    # The control stream parameterises the equilibration HALF-LIFE and
    # derives the rate constant from it (KE0 = LOG(2)/EQHR), which is why
    # Table 2 reports a half-life rather than a rate constant. The canonical
    # library name is the rate constant, so the half-life is converted here:
    # log(2) / 3.19641 h = 0.216845 1/h.
    lke0 <- log(log(2) / 3.19641)
    label("Plasma-to-CSF equilibration rate constant ke0 (1/h); equivalent to the reported equilibration half-life of 3.20 h")  # Table 2 "Equilibration half-life to CSF, HL plasma-CSF, h: 3.20 (2.06-4.93)"; control stream $THETA (0,3.19641,10) ; 18 EQHR --> KE0
    lppc <- log(0.059288)
    label("Pseudo-partition coefficient PPC, the steady-state ratio of total CSF to total plasma rifampicin (fraction)")        # Table 2 "Pseudo-partition coefficient to CSF, PPC plasma-CSF: 0.0593 (0.0544-0.0672)"; control stream $THETA (0,0.059288,1) ; 19 PPC

    # ----- Allometric exponents, fixed a priori ---------------------------
    # Methods, PK Modeling: "Allometric scaling was applied for all
    # disposition parameters via the fixed power exponents of 0.75 for
    # clearance parameters and 1 for volume parameters". Hardcoded in the
    # control stream $PK block rather than estimated as THETAs.
    e_ffm_clint_max <- fixed(0.75)
    label("Allometric exponent of fat-free mass on CLint,max (unitless)")            # control stream ALLMCL_FFM = (FFM/TVFFM)**0.75 with TVFFM = 46
    e_ffm_q <- fixed(0.75)
    label("Allometric exponent of fat-free mass on Q (unitless)")                    # control stream TVQ = THETA(11)*ALLMCL_FFM
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on V (unitless)")                    # control stream ALLMV_FFM = (FFM/TVFFM) with TVFFM = 46
    e_ffm_vp <- fixed(1)
    label("Allometric exponent of fat-free mass on Vp (unitless)")                   # control stream TVV2 = THETA(12)*ALLMV_FFM
    e_ffm_qh <- fixed(0.75)
    label("Allometric exponent of fat-free mass on QH (unitless)")                   # control stream ALLMCL_H_FFM = (FFM/56.1)**0.75
    e_ffm_vh <- fixed(1)
    label("Allometric exponent of fat-free mass on VH (unitless)")                   # control stream ALLMV_H_FFM = (FFM/56.1)

    # =======================================================================
    # Random effects. Table 2 reports each variability as a percentage that
    # is the omega STANDARD DEVIATION on the log scale, not a CV%. Verified
    # against the control stream $OMEGA block: sqrt(0.609918) = 0.781
    # reproduces the reported 78.1% for ka and sqrt(1.23318) = 1.111 the
    # reported 111% for mean transit time, neither of which is recoverable
    # through the log-normal CV formula. The variances below are the control
    # stream's values.
    #
    # Four random effects that the source model carries but fixes to zero are
    # omitted here rather than written as degenerate zero-variance etas:
    # $OMEGA BLOCK(1) FIX 0 for between-subject variability in prehepatic
    # bioavailability (ETA(3)), between-visit variability in CLint,max
    # (ETA(19)/ETA(20)), between-visit variability in PPC (ETA(24)/ETA(25))
    # and between-subject variability in PPC (ETA(26)). None appears in
    # Table 2. The two height-imputation etas ($OMEGA 0.00215 FIX and 0.00170
    # FIX) belong to the missing-covariate imputation, not the structural
    # model, and are likewise omitted.
    # =======================================================================
    etalclint_max ~ 0.0642298
    label("Between-subject variability in CLint,max (log-scale variance)")           # Table 2 BSV in CLint,max*fu = 25.3% (24.1-26.6); control stream $OMEGA BLOCK(1) 0.0642298 ; 1 BSVCL; sqrt(0.0642298) = 0.2534
    etalvc ~ 0.0295196
    label("Between-subject variability in central volume (log-scale variance)")      # Table 2 BSV in central volume = 17.2% (14.8-18.5); control stream $OMEGA BLOCK(1) 0.0295196 ; 2 BSVV; sqrt(0.0295196) = 0.1718
    etaldur ~ 0.0288515
    label("Between-subject variability in the intravenous infusion duration (log-scale variance)")  # Table 2 BSV in infusion duration = 17.0% (14.9-18.6) with footnote g; control stream $OMEGA BLOCK(1) 0.0288515 ; 21 BSV D2; sqrt(0.0288515) = 0.1699

    # Between-occasion variability over five occasions, each a single shared
    # variance ($OMEGA BLOCK(1) <value> followed by four SAME repeats).
    # nlmixr2 has no SAME keyword, so occasions 2-5 repeat the value with
    # fix() (the Abdelgawad_2024_linezolid.R / Gafar_2026_rifampicin.R
    # pattern).
    etaiov_fdepot_1 ~ 0.0329715
    label("Between-occasion variability in prehepatic oral bioavailability, occasion 1 (log-scale variance)")  # Table 2 BOV in F oral,prehepatic = 18.2% (16.7-20.1); control stream $OMEGA BLOCK(1) 0.0329715 ; 4-8 BOVBIO; sqrt(0.0329715) = 0.1816
    etaiov_fdepot_2 ~ fix(0.0329715)                                                                           # control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fix(0.0329715)                                                                           # control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fix(0.0329715)                                                                           # control stream $OMEGA BLOCK(1) SAME
    etaiov_fdepot_5 ~ fix(0.0329715)                                                                           # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_1 ~ 0.609918
    label("Between-occasion variability in ka, occasion 1 (log-scale variance)")                                # Table 2 BOV in ka = 78.1% (55.8-95.6); control stream $OMEGA BLOCK(1) 0.609918 ; 9-13 BOVKA; sqrt(0.609918) = 0.7810
    etaiov_ka_2 ~ fix(0.609918)                                                                                 # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fix(0.609918)                                                                                 # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fix(0.609918)                                                                                 # control stream $OMEGA BLOCK(1) SAME
    etaiov_ka_5 ~ fix(0.609918)                                                                                 # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 1.23318
    label("Between-occasion variability in mean transit time, occasion 1 (log-scale variance)")                 # Table 2 BOV in mean transit time = 111% (87.7-137); control stream $OMEGA BLOCK(1) 1.23318 ; 14-18 BOVMTT; sqrt(1.23318) = 1.1105
    etaiov_mtt_2 ~ fix(1.23318)                                                                                 # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fix(1.23318)                                                                                 # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_4 ~ fix(1.23318)                                                                                 # control stream $OMEGA BLOCK(1) SAME
    etaiov_mtt_5 ~ fix(1.23318)                                                                                 # control stream $OMEGA BLOCK(1) SAME

    # =======================================================================
    # Residual error, one combined proportional-plus-additive model per
    # matrix. The control stream $ERROR block builds each additive term from
    # a theta plus a fraction of the assay's lower limit of quantification:
    #   plasma: ADD_P = THETA(7) + 0.2*LLOQ_P, LLOQ_P = 0.117, THETA(7) = 0 FIX
    #   CSF:    ADD_E = THETA(21) + 0.5*LLOQ_E, LLOQ_E = 0.005
    # and combines them as W = SQRT(ADD**2 + PROP**2), which is exactly
    # nlmixr2's root-sum-square add() + prop().
    # =======================================================================
    propSd <- 0.251954
    label("Proportional residual error for plasma rifampicin (fraction)")            # Table 2 plasma proportional error = 25.2% (22.3-29.7); control stream $THETA (0,0.251954,1) ; 6 PROP
    addSd <- fixed(0.0234)
    label("Additive residual error for plasma rifampicin (mg/L)")                    # Table 2 plasma additive error = 0.0234 mg/L with footnote h "fixed to ... 20% of the lower limit of quantification (0.117 mg/L)"; control stream ADD_P = THETA(7) + 0.117*0.2 with THETA(7) = 0 FIX
    propSd_Ccsf <- 0.984216
    label("Proportional residual error for CSF rifampicin (fraction)")               # Table 2 CSF proportional error = 98.4% (91.8-99.9); control stream $THETA (0,0.984216,1) ; 20 PROP_CSF
    # Table 2 prints the CSF additive error as "2.31 (1.77-2.93)" under the
    # unit heading "ug/mL", which would be 100 times larger than the highest
    # CSF concentration the study ever measured. The control stream resolves
    # it: ADD_E = THETA(21) + 0.5*LLOQ_E = 0.0206258 + 0.5*0.005 = 0.0231258
    # mg/L, i.e. the printed 2.31 is in units of 1e-2 ug/mL. The arithmetic
    # reproduces all three printed digits, and Table 2 footnote i's statement
    # that "the lower boundary of the additive error was fixed to 50% of the
    # lower limit of quantification (0.005 mg/L)" is exactly the 0.5*LLOQ_E
    # term. See the vignette Errata.
    addSd_Ccsf <- 0.0231258
    label("Additive residual error for CSF rifampicin (mg/L)")                       # Table 2 CSF additive error printed as 2.31 (1.77-2.93) x 1e-2 ug/mL; control stream ADD_E = THETA(21) + 0.005*0.5 with $THETA (0,0.0206258,10) ; 21 ADD_CSF
  })

  model({
    # --- 1. Occasion indicators for the between-occasion variability. ------
    #     Control stream $PK: IF (OCC==1) ... IF (OCC==5) multiplexers.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5

    # --- 2. Allometric scaling on fat-free mass. Two reference values are --
    #     in play: 46 kg (the dataset median carried by the control stream as
    #     TVFFM) for the estimated disposition parameters and 56.1 kg for the
    #     fixed hepatic physiology.
    allm_cl <- (FFM / 46)^e_ffm_clint_max
    allm_q  <- (FFM / 46)^e_ffm_q
    allm_v  <- (FFM / 46)^e_ffm_vc
    allm_vp <- (FFM / 46)^e_ffm_vp
    allm_qh <- (FFM / 56.1)^e_ffm_qh
    allm_vh <- (FFM / 56.1)^e_ffm_vh

    # --- 3. Dose-group and visit selection of the typical CLint,max. -------
    #     The four control-stream branches, written as mutually exclusive
    #     indicator products so the solver sees no discontinuity in code path.
    lclint_max_tv <-
      (1 - DOSE_HIGH) * (1 - DAY28) * lclint_max_std_d3 +
      (1 - DOSE_HIGH) * DAY28       * lclint_max_std_d28 +
      DOSE_HIGH       * (1 - DAY28) * lclint_max_high_d3 +
      DOSE_HIGH       * DAY28       * lclint_max_high_d28

    # --- 4. Individual parameters. -----------------------------------------
    clint_max <- exp(lclint_max_tv + etalclint_max) * allm_cl
    vc     <- exp(lvc + etalvc) * allm_v
    vp     <- exp(lvp) * allm_vp
    q      <- exp(lq) * allm_q
    ka     <- exp(lka + iov_ka)
    mtt    <- exp(lmtt + iov_mtt)
    ntr    <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)
    qh     <- exp(lqh) * allm_qh
    vh     <- exp(lvh) * allm_vh
    km     <- exp(lkm)
    ke0    <- exp(lke0)
    ppc    <- exp(lppc)

    # Maximal elimination rate, mg/h. Control stream $PK:
    # VMAX = CL*EXP(LOGKM), where its CL variable is the CLint,max * fu
    # product; here clint_max is the unbound-corrected CLint,max, so fub is
    # reinstated in the extraction ratio below instead.
    vmax <- clint_max * km

    # --- 5. Well-stirred liver with saturable intrinsic clearance. ---------
    #     Figure 1 and control stream $DES:
    #       CH     = A(3)/VH
    #       SAT_CL = VMAX/(CH + KM)      -- i.e. CLint,max * Km / (CH + Km)
    #       EH     = SAT_CL/(SAT_CL + QH)
    #       FH     = 1 - EH
    #     At CH = 0 the intrinsic clearance equals CLint,max and the
    #     extraction ratio is maximal; as the liver concentration rises CLint
    #     falls, EH falls, and exposure grows faster than proportionally with
    #     dose. The control stream guards the division with IF (CH>0), a
    #     no-op because every term it feeds is multiplied by a liver amount;
    #     the smooth form is used here so the solver sees no discontinuity.
    c_liver <- liver / vh
    clint   <- vmax / (c_liver + km)
    eh      <- (clint * fub) / (clint * fub + qh)
    fh      <- 1 - eh

    k30 <- qh * eh / vh   # liver -> eliminated
    k32 <- qh * fh / vh   # liver -> central (the fraction escaping extraction)
    k23 <- qh / vc        # central -> liver (recirculation with hepatic blood flow)
    k24 <- q / vc         # central -> peripheral
    k42 <- q / vp         # peripheral -> central

    # --- 6. ODE system, matching $DES DADT(1)-DADT(5) term for term. -------
    #     transit() is rxode2's implementation of the Savic 2007 analytical
    #     transit-chain input rate and reproduces the control stream's
    #     PIZZA / TRANSIT construction exactly:
    #       TRANSIT = BIO_ORAL*PD*KTR * (KTR*TAD)^NN * exp(-KTR*TAD)/gamma(NN+1)
    #     with KTR = (NN+1)/MTT. Prehepatic bioavailability enters as the bio
    #     argument, exactly as BIO_ORAL enters PIZZA. Absorbed drug is
    #     delivered to the liver, so first-pass extraction is structural
    #     rather than a separate F term.
    d/dt(depot)       <- transit(ntr, mtt, fdepot) - ka * depot
    d/dt(liver)       <- ka * depot - k32 * liver + k23 * central - k30 * liver
    d/dt(central)     <- k32 * liver - k23 * central - k24 * central + k42 * peripheral1
    d/dt(peripheral1) <- k24 * central - k42 * peripheral1

    Cc <- central / vc

    # The CSF state holds a CONCENTRATION, not an amount: the effect
    # compartment is assumed to have negligible volume and negligible mass
    # transfer with plasma (supplementary "Effect compartment modelling for
    # CSF"), and $ERROR reads the CSF prediction as CE = A(5) directly.
    d/dt(csf) <- ke0 * (ppc * Cc - csf)

    # F1 = 0 in the control stream: the oral dose event must not deposit a
    # bolus into the absorption compartment because transit() already feeds
    # the whole dose in through the analytical chain. transit() reads the raw
    # amount from podo(depot) regardless of f(depot).
    f(depot) <- 0

    # Intravenous doses go straight into central with F fixed to 1, over a
    # modelled 1 h infusion carrying between-subject variability (control
    # stream: F2 = 0; IF(RIFIV.EQ.1) F2 = BIO_IV, and D2 = DUR*EXP(BSVD2)).
    # Dose records targeting `central` must therefore set rate = -2 so
    # rxode2 uses the modelled duration.
    f(central)   <- exp(lfcentral)
    dur(central) <- exp(ldur + etaldur)

    # --- 7. Observations. --------------------------------------------------
    Ccsf <- csf
    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}

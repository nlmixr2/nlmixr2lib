Riggs_2012_albinterferon <- function() {
  description <- "One-compartment population PK model with first-order absorption for subcutaneous albinterferon alfa-2b (albIFN, a genetic fusion of recombinant human albumin and recombinant interferon alfa-2b) in 1984 patients with chronic hepatitis C virus infection, from Riggs 2012. Apparent clearance, apparent central volume, and the absorption rate constant each carry a prespecified FULL covariate model (weight, age, sex, race, HCV genotype, baseline HCV RNA, and - on CL/F only - serum albumin, ALT, creatinine clearance, ribavirin dose, immunogenicity status, and CYP2D6 / CYP3A4 inhibitor comedication), plus a relative-bioavailability and an absorption-rate effect of thigh-versus-abdomen injection site. Interindividual variability is exponential on all three disposition parameters with a block CL/F-V/F OMEGA; residual error is combined proportional plus additive. Body weight was the only covariate the authors judged clinically meaningful. Sister model files from the same paper (post-hoc logistic sustained-virologic-response exposure-response): modellib('Riggs_2012_albinterferon_svr_gt1') and modellib('Riggs_2012_albinterferon_svr_gt23')."
  reference <- paste(
    "Riggs MM, Bergsma TT, Rogers JA, Gastonguay MR, Subramanian GM, Chen C,",
    "Devalaraja M, Corey AE, Sun H, Yu J, Stein DS.",
    "Population pharmacokinetics and exposure-response of albinterferon alfa-2b.",
    "J Clin Pharmacol. 2012;52(4):475-486.",
    "doi:10.1177/0091270011399576.",
    "Structural and covariate equations are the four displayed equations on p. 479;",
    "final parameter estimates are in Table I (p. 476-477);",
    "random-effect and residual-error estimates are in the Results,",
    "'Population PK Analysis' (p. 480).",
    "Sister model files from the same paper:",
    "modellib('Riggs_2012_albinterferon_svr_gt1'),",
    "modellib('Riggs_2012_albinterferon_svr_gt23').",
    sep = " "
  )
  vignette <- "Riggs_2012_albinterferon"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline, time-fixed per subject. Enters CL/F, V/F, and ka as a power function of TWT = WT / 75 kg (p. 479 'where' block). Reference 75 kg is the hypothetical reference individual of Table I footnote b. Observed range 38-166 kg (Results, 'Population PK Analysis'). The only covariate the authors judged clinically meaningful: relative to the 75 kg reference, exposure rises about 30% at 50 kg and falls about 20% and 30% at 100 and 125 kg.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline, time-fixed per subject. Enters CL/F, V/F, and ka as a power function of TAGE = AGE / 45 y (p. 479 'where' block). Observed range 18-79 years.",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (Table I footnote c, 'value at baseline'). The source calibrated the coefficients against US-convention g/dL with reference 4.3 g/dL (Table I footnote b and Figure 1 caption), so model() applies the register-mandated inline conversion alb_gdL <- ALB * 0.1 before forming TALB = alb_gdL / 4.3. NOTE: the p. 479 'where' block prints TALB with 'g/L' in both numerator and denominator while giving the denominator as 4.3, which is the g/dL value; the ratio is unit-free provided both sides use the same unit, and the reference individual definition in Table I footnote b resolves the intended unit as g/dL. Enters CL/F and V/F only.",
      source_name        = "ALB"
    ),
    ALT = list(
      description        = "Baseline serum alanine aminotransferase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (Table I footnote c). Enters CL/F only, as a power function of TALT = ALT / 34 IU/L (p. 479 'where' block; reference from Table I footnote b, which prints '34 IU/mL' - a unit typo for IU/L, corrected in the Figure 1 caption's 'alanine transaminase (ALT), 34 IU/mL' and resolved by the equation block's explicit 'ALT (IU/L)_i / 34 IU/L'). IU/L and U/L are used interchangeably.",
      source_name        = "ALT"
    ),
    CRCL = list(
      description        = "Baseline estimated creatinine clearance.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "RAW, NOT BSA-normalized (the source writes 'CRCN (mL/min)_i = CRCL (mL/min)_i', p. 479), so supply mL/min rather than the register's default mL/min/1.73 m^2 - same convention as Delattre_2010_amikacin.R and Chen_2023_nemonoxacin.R. Time-fixed per subject (Table I footnote c). The source CAPS the value at 150 mL/min before normalizing ('if CRCN >= 150 mL/min then CRCN = 150 mL/min'), and model() reproduces the cap. Enters CL/F only, as a power function of TCRL = CRCN / 120 mL/min.",
      source_name        = "CRCL"
    ),
    SEXF = list(
      description        = "Female sex indicator: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the source's reference individual is a woman)",
      notes              = "Time-fixed per subject. REFERENCE-CATEGORY INVERSION relative to the usual SEXF usage: the source's reference individual is a 45-year-old white WOMAN (Table I footnote b), and the estimated coefficients are MALE effects, so model() forms male <- 1 - SEXF and raises the coefficients to that indicator. Cohort 1189 men and 795 women (Results, 'Population PK Analysis'); note the source reports these counts in the order 'men and women' while its own total is 1984.",
      source_name        = "Male (1 = male)"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator: 1 = Asian, 0 = other race category.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; the source's reference individual)",
      notes              = "Time-fixed per subject. One of three mutually exclusive non-White indicators (RACE_ASIAN, RACE_BLACK, RACE_OTHER) whose common reference is White. Cohort 13% Asian. Enters CL/F, V/F, and ka.",
      source_name        = "Asian"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator: 1 = Black, 0 = other race category.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; the source's reference individual)",
      notes              = "Time-fixed per subject. Cohort 5% Black. Enters CL/F, V/F, and ka.",
      source_name        = "Black"
    ),
    RACE_OTHER = list(
      description        = "Race-category 'Other' indicator: 1 = a race other than White, Asian, or Black; 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White; the source's reference individual)",
      notes              = "Time-fixed per subject. The source writes this indicator as 'Race = other' in the p. 479 equations. Cohort composition given as 81% White, 13% Asian, 5% Black, leaving about 1% 'Other'. Enters CL/F, V/F, and ka.",
      source_name        = "Race = other"
    ),
    HCV_GT2 = list(
      description        = "Hepatitis C virus genotype-2 indicator: 1 = HCV genotype 2, 0 = any other genotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HCV genotype 1; the source's reference individual)",
      notes              = "Time-fixed per subject (genotype is determined at infection). One of three mutually exclusive non-genotype-1 indicators (HCV_GT2, HCV_GT3, HCV_GT4) whose common reference is genotype 1. 66% of the cohort had genotype 1; the remainder had genotype 2 or 3 apart from 5 patients with genotype 4. Enters CL/F, V/F, and ka.",
      source_name        = "HCV genotype 2"
    ),
    HCV_GT3 = list(
      description        = "Hepatitis C virus genotype-3 indicator: 1 = HCV genotype 3, 0 = any other genotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HCV genotype 1; the source's reference individual)",
      notes              = "Time-fixed per subject. Enters CL/F, V/F, and ka.",
      source_name        = "HCV genotype 3"
    ),
    HCV_GT4 = list(
      description        = "Hepatitis C virus genotype-4 indicator: 1 = HCV genotype 4, 0 = any other genotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HCV genotype 1; the source's reference individual)",
      notes              = "Time-fixed per subject. Only 5 of 1984 patients carried genotype 4, and the authors retained the indicator specifically to quantify that lack of information rather than to declare the effect insignificant (Discussion, p. 484). The corresponding coefficients are correspondingly imprecise - the ka effect 0.202 has a bootstrap 95% CI of 0.102 to 1.4 - so simulations at HCV_GT4 = 1 are extrapolation, not estimation. Enters CL/F, V/F, and ka.",
      source_name        = "HCV genotype 4"
    ),
    HCV_VLOAD = list(
      description        = "Baseline hepatitis C virus RNA concentration in serum.",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (Table I footnote c, 'value at baseline'). The popPK model uses it only DICHOTOMIZED at 800,000 IU/mL: model() forms hcvHigh <- (HCV_VLOAD >= 800000), matching the source's '(theta)^{Baseline HCV RNA >= 800,000 IU/mL}' term. The reference individual has baseline HCV RNA < 800,000 IU/mL (Table I footnote b). Supplied as the continuous concentration rather than as a pre-binarized flag so the same column also drives the 400,000 / 800,000 three-level categorization used by the sister SVR exposure-response models. Enters CL/F, V/F, and ka.",
      source_name        = "Baseline HCV RNA"
    ),
    DOSE_RBV_MGD = list(
      description        = "Total daily oral ribavirin dose administered concomitantly.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The source parameterizes ribavirin as TRIB = (number of 200 mg tablets per day) / 6, i.e. the daily dose normalized to 1200 mg/day (p. 479 'where' block); model() forms trib <- (DOSE_RBV_MGD / 200) / 6, so DOSE_RBV_MGD = 1200 gives TRIB = 1 and reproduces the reference individual. The effect is a SPLIT term: when TRIB > 0 the CL/F multiplier is TRIB^e_rbv_cl (Table I 'RBV present', 0.0436); when TRIB = 0 the multiplier is the separate constant e_norbv_cl (Table I 'RBV not present', 1.1). Ribavirin was given daily per the standard of care for each HCV genotype. Enters CL/F only.",
      source_name        = "RIBA (# of 200 mg tabs/day)"
    ),
    ADA_POS = list(
      description        = "Anti-albinterferon immunogenicity status: 1 = positive for anti-interferon antibodies, 0 = negative.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (negative for immunogenicity; the source's reference individual)",
      notes              = "The source's term is '(theta18)^{Positive for IFN immunogenicity}'. Treated as time-fixed per subject here because the source's equation carries no time index beyond the subject subscript and Table I reports a single coefficient. Enters CL/F only.",
      source_name        = "Positive for IFN immunogenicity"
    ),
    CONMED_CYP2D6_INH = list(
      description        = "Concomitant medication classified as a CYP2D6 inhibitor: 1 = the patient was using one, 0 = not.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP2D6-inhibitor comedication)",
      notes              = "The source's term is '(theta42)^{Patient using con med classified as CYP2D6 inhibitor}'. The paper does not enumerate which agents were classified into the class. Note that albinterferon is an albumin-interferon fusion protein cleared by proteolysis rather than by cytochrome-P450 metabolism, so this is best read as a marker of comedication burden or of the underlying condition prompting it; the estimated effect (1.05, bootstrap 95% CI 0.979-1.12) sits well inside the authors' +/-25% no-clinical-relevance band. Enters CL/F only.",
      source_name        = "Patient using con med classified as CYP2D6 inhibitor"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant medication classified as a CYP3A inhibitor: 1 = the patient was using one, 0 = not.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A-inhibitor comedication)",
      notes              = "The source's term is '(theta43)^{Patient using con med classified as CYP3A inhibitor}'; Table I labels the row 'CYP 3A4 inhibitor present'. Same interpretive caveat as CONMED_CYP2D6_INH - albinterferon is not a CYP substrate - and the estimated effect (0.961, bootstrap 95% CI 0.9-1.04) is likewise inside the no-clinical-relevance band. Enters CL/F only.",
      source_name        = "Patient using con med classified as CYP3A inhibitor"
    ),
    INJSITE_THIGH = list(
      description        = "Subcutaneous injection-site indicator: 1 = thigh, 0 = abdomen.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen; the source's reference individual)",
      notes              = "albIFN was self-administered subcutaneously in the thigh or the abdomen. The site acts in TWO places: a multiplicative effect on ka (Table I 'Thigh as injection site', 1.09) and a relative-bioavailability effect on the depot (Table I 'Frel for thigh as injection site', 1.07, footnote d 'bioavailability relative to abdomen as injection site'). The source encodes the latter as the switch 'F1 = 1; if injection site = thigh, then F1 = theta44' (p. 479). Because CL/F and V/F are APPARENT parameters already carrying the abdomen bioavailability, e_injsite_thigh_fdepot is a RELATIVE bioavailability and lfdepot is anchored at log(1) for the abdomen reference.",
      source_name        = "Injection site = thigh"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "albinterferon alfa-2b", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "albinterferon alfa-2b", units = "ug", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1984L,
    n_studies      = 5L,
    age_range      = "18-79 years",
    weight_range   = "38-166 kg",
    sex_female_pct = 40.1,
    race_ethnicity = c(White = 81, Asian = 13, Black = 5, Other = 1),
    disease_state  = "Chronic hepatitis C virus infection with detectable serum HCV RNA and compensated liver disease. 66% HCV genotype 1, the remainder genotype 2 or 3 apart from 5 patients with genotype 4. Interferon-treatment-naive except for one phase 2 trial (NCT00097435), which enrolled patients who had failed previous interferon alfa treatment.",
    dose_range     = "900-1800 ug albinterferon alfa-2b, self-administered subcutaneously in the thigh or abdomen once every 2 weeks or once every 4 weeks, for 48 weeks (genotype 1) or 24 weeks (genotype 2/3); up to 72 weeks for late responders in NCT00097435. All patients also received daily oral ribavirin per the standard of care for their genotype. The 1200 ug every-2-weeks arms of the phase 3 trials were reduced to 900 ug during the studies because of serious pulmonary adverse events.",
    regions        = "multinational",
    notes          = "Pooled from three phase 2 and two phase 3 randomized trials (NCT00656006, NCT00097435, NCT00115908, NCT00411385, NCT00402428), contributing 12,042 serum albIFN concentrations. Full baseline demographics are in Supplementary Table I, which is not on disk; the figures given here are those stated in the main-text Results, 'Population PK Analysis'. sex_female_pct is computed from the reported counts (795 women of 1984). Serum albIFN was measured by ELISA with a lower limit of quantitation of 0.53 ng/mL for all phase 3 specimens and most phase 2 specimens, and 0.26 ng/mL for the remainder. Hypothetical reference individual for the Table I estimates: 45-year-old white woman, 75 kg, HCV genotype 1, negative for immunogenicity, baseline HCV RNA < 800,000 IU/mL, ribavirin 1200 mg/day, abdomen as injection site, albumin 4.3 g/dL, ALT 34 IU/L, estimated creatinine clearance 120 mL/min."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters, for the Table I footnote b reference individual.
    # The source reports CL/F in mL/h and V/F in L; both are converted to
    # L here so that central (ug) / vc (L) = ug/L = ng/mL, the reported
    # concentration unit, with no further scaling.
    # ------------------------------------------------------------------
    lka <- log(0.0148)   ; label("Absorption rate constant ka (1/h)")           # Table I, row "None", ka column: theta3 = 0.0148 1/h (%RSE 8.3; bootstrap median 0.0149, 95% CI 0.0107-0.0173)
    lcl <- log(0.0389)   ; label("Apparent clearance CL/F (L/h)")               # Table I, row "None", CL/F column: theta1 = 38.9 mL/h = 0.0389 L/h (%RSE 1.74; bootstrap median 38.8, 95% CI 37.5-40.2 mL/h)
    lvc <- log(11.6)     ; label("Apparent central volume of distribution V/F (L)") # Table I, row "None", V/F column: theta2 = 11.6 L (%RSE 4.6; bootstrap median 11.5, 95% CI 9.69-12.7 L)

    # Bioavailability anchor. CL/F and V/F are APPARENT parameters, so the
    # absolute bioavailability of the reference (abdomen) site is already
    # folded into them and F1 is anchored at 1 by the source itself:
    # "F1 = 1; if injection site = thigh, then F1 = theta44" (p. 479).
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability of the depot at the abdomen reference injection site (unitless)") # p. 479 structural equations: F1 = 1 for the abdomen reference

    # ------------------------------------------------------------------
    # Covariate effects on CL/F (p. 479 equation 1; Table I CL/F columns).
    # Continuous covariates enter as power exponents on a normalized ratio;
    # categorical covariates enter as multiplicative factors raised to a
    # 0/1 indicator.
    # ------------------------------------------------------------------
    e_wt_cl                <- 0.686   ; label("Exponent on (WT / 75 kg) for CL/F (unitless)")                        # Table I, row "Weight, kg", CL/F column: theta4 = 0.686 (%RSE 6.12; 95% CI 0.607-0.761)
    e_age_cl               <- 0.0221  ; label("Exponent on (AGE / 45 y) for CL/F (unitless)")                        # Table I, row "Age, y", CL/F column: theta5 = 0.0221 (%RSE 128; 95% CI -0.0395-0.0727)
    e_alb_cl               <- 0.507   ; label("Exponent on (ALB / 4.3 g/dL) for CL/F (unitless)")                    # Table I, row "Albumin, g/dL", CL/F column: theta10 = 0.507 (%RSE 13.2; 95% CI 0.367-0.65)
    e_alt_cl               <- 0.0224  ; label("Exponent on (ALT / 34 IU/L) for CL/F (unitless)")                     # Table I, row "ALT, IU/L", CL/F column: theta11 = 0.0224 (%RSE 29.1; 95% CI 0.00912-0.0338)
    e_crcl_cl              <- 0.116   ; label("Exponent on (min(CRCL, 150) / 120 mL/min) for CL/F (unitless)")       # Table I, row "CrCl, mL/min", CL/F column: theta12 = 0.116 (%RSE 23; 95% CI 0.0588-0.17)
    e_male_cl              <- 1.05    ; label("Multiplicative effect of male sex on CL/F (unitless)")                # Table I, row "Male sex", CL/F column: theta6 = 1.05 (%RSE 1.68; 95% CI 1.02-1.08)
    e_race_asian_cl        <- 1.03    ; label("Multiplicative effect of Asian race on CL/F (unitless)")              # Table I, row "Asian race", CL/F column: theta7 = 1.03 (%RSE 2.47; 95% CI 0.975-1.08)
    e_race_black_cl        <- 0.926   ; label("Multiplicative effect of Black race on CL/F (unitless)")              # Table I, row "Black race", CL/F column: theta8 = 0.926 (%RSE 5.79; 95% CI 0.827-1.03)
    e_race_other_cl        <- 0.877   ; label("Multiplicative effect of race category Other on CL/F (unitless)")     # Table I, row "Other race", CL/F column: theta9 = 0.877 (%RSE 5.73; 95% CI 0.788-0.968)
    e_hcv_gt2_cl           <- 0.971   ; label("Multiplicative effect of HCV genotype 2 on CL/F (unitless)")          # Table I, row "HCV genotype 2", CL/F column: theta13 = 0.971 (%RSE 2.21; 95% CI 0.934-1.02)
    e_hcv_gt3_cl           <- 0.994   ; label("Multiplicative effect of HCV genotype 3 on CL/F (unitless)")          # Table I, row "HCV genotype 3", CL/F column: theta14 = 0.994 (%RSE 1.89; 95% CI 0.957-1.03)
    e_hcv_gt4_cl           <- 0.886   ; label("Multiplicative effect of HCV genotype 4 on CL/F (unitless)")          # Table I, row "HCV genotype 4", CL/F column: theta15 = 0.886 (%RSE 24.7; 95% CI 0.462-1.25)
    e_hcv_vload_cl         <- 0.988   ; label("Multiplicative effect of baseline HCV RNA >= 800,000 IU/mL on CL/F (unitless)") # Table I, row "HCV RNA >= 800,000 IU/mL", CL/F column: theta16 = 0.988 (%RSE 1.48; 95% CI 0.965-1.02)
    e_rbv_cl               <- 0.0436  ; label("Exponent on (ribavirin dose / 1200 mg/day) for CL/F when ribavirin is present (unitless)") # Table I, row "RBV present", CL/F column: theta17 = 0.0436 (%RSE 58.2; 95% CI -0.00115-0.0943)
    e_norbv_cl             <- 1.1     ; label("Multiplicative effect on CL/F when no ribavirin is present (unitless)") # Table I, row "RBV not present", CL/F column: theta41 = 1.1 (%RSE 3.67; 95% CI 1.01-1.18)
    e_ada_cl               <- 0.942   ; label("Multiplicative effect of positive interferon immunogenicity on CL/F (unitless)") # Table I, row "Immunogenicity positive", CL/F column: theta18 = 0.942 (%RSE 1.83; 95% CI 0.911-0.973)
    e_conmed_cyp2d6_inh_cl <- 1.05    ; label("Multiplicative effect of CYP2D6-inhibitor comedication on CL/F (unitless)") # Table I, row "CYP 2D6 inhibitor present", CL/F column: theta42 = 1.05 (%RSE 3.16; 95% CI 0.979-1.12)
    e_conmed_cyp3a4_inh_cl <- 0.961   ; label("Multiplicative effect of CYP3A-inhibitor comedication on CL/F (unitless)") # Table I, row "CYP 3A4 inhibitor present", CL/F column: theta43 = 0.961 (%RSE 3.95; 95% CI 0.9-1.04)

    # ------------------------------------------------------------------
    # Covariate effects on V/F (p. 479 equation 2; Table I V/F columns).
    # No ALT, creatinine-clearance, ribavirin, immunogenicity, or CYP
    # comedication terms appear on V/F - those Table I cells are dashes.
    # ------------------------------------------------------------------
    e_wt_vc         <- 0.888   ; label("Exponent on (WT / 75 kg) for V/F (unitless)")                        # Table I, row "Weight, kg", V/F column: theta19 = 0.888 (%RSE 10.8; 95% CI 0.66-1.09)
    e_age_vc        <- -0.0713 ; label("Exponent on (AGE / 45 y) for V/F (unitless)")                        # Table I, row "Age, y", V/F column: theta20 = -0.0713 (%RSE 88.1; 95% CI -0.19-0.0687)
    e_alb_vc        <- 0.358   ; label("Exponent on (ALB / 4.3 g/dL) for V/F (unitless)")                    # Table I, row "Albumin, g/dL", V/F column: theta25 = 0.358 (%RSE 42.3; 95% CI 0.0419-0.678)
    e_male_vc       <- 0.891   ; label("Multiplicative effect of male sex on V/F (unitless)")                # Table I, row "Male sex", V/F column: theta21 = 0.891 (%RSE 4.29; 95% CI 0.81-0.981)
    e_race_asian_vc <- 1.15    ; label("Multiplicative effect of Asian race on V/F (unitless)")              # Table I, row "Asian race", V/F column: theta22 = 1.15 (%RSE 5.32; 95% CI 0.3-1.33)
    e_race_black_vc <- 0.932   ; label("Multiplicative effect of Black race on V/F (unitless)")              # Table I, row "Black race", V/F column: theta23 = 0.932 (%RSE 12.6; 95% CI 0.291-1.14)
    e_race_other_vc <- 0.996   ; label("Multiplicative effect of race category Other on V/F (unitless)")     # Table I, row "Other race", V/F column: theta24 = 0.996 (%RSE 12.1; 95% CI 0.455-1.32)
    e_hcv_gt2_vc    <- 0.952   ; label("Multiplicative effect of HCV genotype 2 on V/F (unitless)")          # Table I, row "HCV genotype 2", V/F column: theta26 = 0.952 (%RSE 4.42; 95% CI 0.829-1.06)
    e_hcv_gt3_vc    <- 0.891   ; label("Multiplicative effect of HCV genotype 3 on V/F (unitless)")          # Table I, row "HCV genotype 3", V/F column: theta27 = 0.891 (%RSE 4.54; 95% CI 0.211-0.982)
    e_hcv_gt4_vc    <- 0.595   ; label("Multiplicative effect of HCV genotype 4 on V/F (unitless)")          # Table I, row "HCV genotype 4", V/F column: theta28 = 0.595 (%RSE 29.1; 95% CI 0.321-1.35)
    e_hcv_vload_vc  <- 0.989   ; label("Multiplicative effect of baseline HCV RNA >= 800,000 IU/mL on V/F (unitless)") # Table I, row "HCV RNA >= 800,000 IU/mL", V/F column: theta29 = 0.989 (%RSE 3.66; 95% CI 0.91-1.08)

    # ------------------------------------------------------------------
    # Covariate effects on ka (p. 479 equation 3; Table I ka columns).
    # Note the source's theta numbering is not monotonic here: the weight
    # exponent is theta39 and the age exponent is theta30.
    # ------------------------------------------------------------------
    e_wt_ka             <- 0.0548 ; label("Exponent on (WT / 75 kg) for ka (unitless)")                        # Table I, row "Weight, kg", ka column: theta39 = 0.0548 (%RSE 366; 95% CI -0.432-0.424)
    e_age_ka            <- -0.541 ; label("Exponent on (AGE / 45 y) for ka (unitless)")                        # Table I, row "Age, y", ka column: theta30 = -0.541 (%RSE 28.3; 95% CI -0.849--0.14)
    e_male_ka           <- 1.16   ; label("Multiplicative effect of male sex on ka (unitless)")                # Table I, row "Male sex", ka column: theta31 = 1.16 (%RSE 8.06; 95% CI 0.981-1.35)
    e_race_asian_ka     <- 0.867  ; label("Multiplicative effect of Asian race on ka (unitless)")              # Table I, row "Asian race", ka column: theta32 = 0.867 (%RSE 14.4; 95% CI 0.215-1.15)
    e_race_black_ka     <- 0.918  ; label("Multiplicative effect of Black race on ka (unitless)")              # Table I, row "Black race", ka column: theta33 = 0.918 (%RSE 44.3; 95% CI 0.22-1.91)
    e_race_other_ka     <- 1.22   ; label("Multiplicative effect of race category Other on ka (unitless)")     # Table I, row "Other race", ka column: theta34 = 1.22 (%RSE 42.7; 95% CI 0.322-2.98)
    e_hcv_gt2_ka        <- 1.14   ; label("Multiplicative effect of HCV genotype 2 on ka (unitless)")          # Table I, row "HCV genotype 2", ka column: theta35 = 1.14 (%RSE 9.48; 95% CI 0.905-1.49)
    e_hcv_gt3_ka        <- 0.854  ; label("Multiplicative effect of HCV genotype 3 on ka (unitless)")          # Table I, row "HCV genotype 3", ka column: theta36 = 0.854 (%RSE 9.84; 95% CI 0.208-1.15)
    e_hcv_gt4_ka        <- 0.202  ; label("Multiplicative effect of HCV genotype 4 on ka (unitless)")          # Table I, row "HCV genotype 4", ka column: theta37 = 0.202 (%RSE 32.8; 95% CI 0.102-1.4)
    e_hcv_vload_ka      <- 1.09   ; label("Multiplicative effect of baseline HCV RNA >= 800,000 IU/mL on ka (unitless)") # Table I, row "HCV RNA >= 800,000 IU/mL", ka column: theta38 = 1.09 (%RSE 7.92; 95% CI 0.917-1.28)
    e_injsite_thigh_ka  <- 1.09   ; label("Multiplicative effect of thigh (vs abdomen) injection site on ka (unitless)") # Table I, row "Thigh as injection site", ka column: theta40 = 1.09 (%RSE 5.94; 95% CI 0.964-1.22)

    # Relative bioavailability effect of injection site (p. 479: F1 = theta44 for thigh).
    e_injsite_thigh_fdepot <- 1.07 ; label("Bioavailability of the thigh injection site relative to the abdomen (unitless)") # Table I, row "Frel for thigh as injection site" (footnote d), ka block: theta44 = 1.07 (%RSE 1.2; 95% CI 1.04-1.09)

    # ------------------------------------------------------------------
    # Interindividual variability. The source reports EXPONENTIAL random
    # effects summarized as coefficients of variation (Results, p. 480):
    # CL/F 21% (95% CI 20-22%), V/F 34% (23-39%), ka 24% (1-60%). Those
    # are converted here with the exact log-normal identity
    #   omega^2 = log(CV^2 + 1)
    # giving 0.0431534, 0.1093834, and 0.0559616 respectively.
    #
    # The source states that a COVARIANCE term between the CL/F (ETA1)
    # and V/F (ETA2) random effects was estimated, but reports it only
    # qualitatively - "There was minimal covariance between interindividual
    # random effects for CL/F (ETA1) and V/F (ETA2)" (p. 480) - with no
    # numeric value anywhere in the paper. The block structure is preserved
    # here with the off-diagonal entered as 0. See the vignette's
    # Assumptions and deviations section.
    # ------------------------------------------------------------------
    # Source trace for the block below. Kept ABOVE the c(...) rather than as
    # a trailing comment: a trailing comment attached to a multi-line omega
    # block breaks the comment-to-label parsing used by the convention lint.
    #   etalcl variance 0.0431534 = log(0.21^2 + 1); Results p. 480, CL/F
    #     interindividual CV 21% (95% CI 20%, 22%)
    #   etalvc variance 0.1093834 = log(0.34^2 + 1); Results p. 480, V/F
    #     interindividual CV 34% (95% CI 23%, 39%)
    #   off-diagonal 0: estimated by the source but reported only
    #     qualitatively as minimal, with no numeric value published anywhere
    etalcl + etalvc ~ c(0.0431534,
                        0.0000000, 0.1093834)
    # etalka variance 0.0559616 = log(0.24^2 + 1); Results p. 480, ka
    #   interindividual CV 24% (95% CI 1%, 60%)
    etalka ~ 0.0559616

    # ------------------------------------------------------------------
    # Residual error: combined proportional plus additive (Results p. 479
    # "residual random effects were described with a combined additive and
    # proportional model").
    # ------------------------------------------------------------------
    propSd <- 0.27 ; label("Proportional residual error (fraction)")  # Results p. 480: proportional residual CV 27% (95% CI 27%, 28%)
    addSd  <- 1.51 ; label("Additive residual error (ng/mL)")         # Results p. 480: additive residual SD 1.51 ng/mL (95% CI 0.91, 2.19 ng/mL)
  })

  model({
    # ----------------------------------------------------------------
    # 1. Derived covariate terms (p. 479 "where" block)
    # ----------------------------------------------------------------
    twt  <- WT / 75                       # WT (kg) / 75 kg
    tage <- AGE / 45                      # AGE (y) / 45 y

    # Serum albumin: the register holds ALB in SI g/L, the source's
    # coefficients are calibrated against US-convention g/dL.
    alb_gdL <- ALB * 0.1                  # SI g/L -> US-convention g/dL
    talb    <- alb_gdL / 4.3              # ALB (g/dL) / 4.3 g/dL

    talt <- ALT / 34                      # ALT (IU/L) / 34 IU/L

    # Creatinine clearance is CAPPED at 150 mL/min before normalizing:
    # "if CRCN >= 150 mL/min then CRCN = 150 mL/min" (p. 479).
    crcn <- min(CRCL, 150)
    tcrl <- crcn / 120                    # CRCN (mL/min) / 120 mL/min

    # Ribavirin: TRIB is the daily dose in 200 mg tablets normalized to 6
    # tablets/day (= 1200 mg/day). The CL/F effect is a SPLIT term:
    #   if TRIB > 0 then CLRB = TRIB^theta17
    #   if TRIB = 0 then CLRB = theta41
    trib  <- (DOSE_RBV_MGD / 200) / 6
    rbvOn <- (trib > 0)
    clrb  <- rbvOn * trib^e_rbv_cl + (1 - rbvOn) * e_norbv_cl

    # The source's reference individual is a WOMAN and the estimated
    # coefficients are MALE effects, so invert the canonical SEXF column.
    male <- 1 - SEXF

    # The popPK model uses baseline HCV RNA only dichotomized at 800,000 IU/mL.
    hcvHigh <- (HCV_VLOAD >= 800000)

    # ----------------------------------------------------------------
    # 2. Individual parameters (p. 479 equations 1-3)
    # ----------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      twt^e_wt_cl * tage^e_age_cl * talb^e_alb_cl * talt^e_alt_cl * tcrl^e_crcl_cl *
      e_male_cl^male *
      e_race_asian_cl^RACE_ASIAN * e_race_black_cl^RACE_BLACK * e_race_other_cl^RACE_OTHER *
      e_hcv_gt2_cl^HCV_GT2 * e_hcv_gt3_cl^HCV_GT3 * e_hcv_gt4_cl^HCV_GT4 *
      e_hcv_vload_cl^hcvHigh *
      clrb *
      e_ada_cl^ADA_POS *
      e_conmed_cyp2d6_inh_cl^CONMED_CYP2D6_INH *
      e_conmed_cyp3a4_inh_cl^CONMED_CYP3A4_INH

    vc <- exp(lvc + etalvc) *
      twt^e_wt_vc * tage^e_age_vc * talb^e_alb_vc *
      e_male_vc^male *
      e_race_asian_vc^RACE_ASIAN * e_race_black_vc^RACE_BLACK * e_race_other_vc^RACE_OTHER *
      e_hcv_gt2_vc^HCV_GT2 * e_hcv_gt3_vc^HCV_GT3 * e_hcv_gt4_vc^HCV_GT4 *
      e_hcv_vload_vc^hcvHigh

    ka <- exp(lka + etalka) *
      twt^e_wt_ka * tage^e_age_ka *
      e_male_ka^male *
      e_race_asian_ka^RACE_ASIAN * e_race_black_ka^RACE_BLACK * e_race_other_ka^RACE_OTHER *
      e_hcv_gt2_ka^HCV_GT2 * e_hcv_gt3_ka^HCV_GT3 * e_hcv_gt4_ka^HCV_GT4 *
      e_hcv_vload_ka^hcvHigh *
      e_injsite_thigh_ka^INJSITE_THIGH

    # ----------------------------------------------------------------
    # 3. Micro-constants
    # ----------------------------------------------------------------
    kel <- cl / vc

    # ----------------------------------------------------------------
    # 4. ODE system: one compartment with first-order absorption
    # ----------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----------------------------------------------------------------
    # 5. Bioavailability: F1 = 1 (abdomen) or theta44 (thigh)
    # ----------------------------------------------------------------
    f(depot) <- exp(lfdepot) * e_injsite_thigh_fdepot^INJSITE_THIGH

    # ----------------------------------------------------------------
    # 6. Observation. central is in ug and vc in L, so central / vc is
    #    ug/L, which is numerically identical to the reported ng/mL.
    # ----------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

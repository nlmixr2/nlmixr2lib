Wurthwein_2025_pegasparaginase_r2ea <- function() {
  description <- "Fourteen-compartment de-PEGylation transit population PK model for intravenous PEGylated asparaginase (PEG-ASNase) in non-high-risk children with acute lymphoblastic leukemia randomised to the experimental arm of the R2 randomisation (R2-EA) of the AIEOP-BFM ALL 2009 trial, German/Czech group (Wurthwein 2025, ESM Table S9 model 493158). Covers twelve administrations: induction protocol IA days 12 and 26, then the ten biweekly re-induction and maintenance doses from protocol II day 8 through maintenance dose M10. Asparaginase activity is carried by a chain of 14 serial compartments sharing one serum volume; drug moves down the chain with intercompartmental clearance Qtr, which mimics stepwise de-PEGylation, and every compartment is eliminated with the initial clearance CLinitial while the terminal compartment is eliminated with CLinitial + Qtr. Body surface area enters volume and the two clearance terms linearly, centred on 0.79 m^2. Initial clearance rises linearly with age above 8 years, is lower in females, and rises with a pre-existing anti-polyethylene-glycol IgM antibody level above a hockey-stick cut point, the antibody effect acting only on the first induction dose. Repeated dosing drives initial clearance down in three steps to 60.1 percent below the first induction dose, the largest accumulation effect reported in the trial. Inter-individual variability on initial clearance, inter-occasion variability on both initial clearance and volume, and combined proportional plus additive residual error."
  reference   <- "Wurthwein G, Siebel C, Lanvers-Kaminsky C, Smisek P, Nath CE, Matteo C, Rizzari C, Schrappe M, Boos J. PEGylated Asparaginase in Children with Acute Lymphoblastic Leukemia Treated within the AIEOP-BFM ALL 2009 Trial: Population Pharmacokinetics and Drug Exposure. Eur J Drug Metab Pharmacokinet. 2025;50(6):683-696. doi:10.1007/s13318-025-00962-3"
  vignette    <- "Wurthwein_2025_pegasparaginase"
  units       <- list(time = "day", dosing = "U", concentration = "U/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed by the Mosteller formula (Wurthwein 2025 Methods, Pharmacokinetic Model, citing reference 14). German/Czech non-high-risk median 0.74 m^2, range 0.39-2.44 (Wurthwein 2025 Table 1). BSA enters volume and the two clearance terms LINEARLY and CENTRED, not allometrically -- ESM Table S9 footnote (a), 'linear increase in V / CLinitial / Qtr with BSA (centred on the median)'. The centring constant 0.79 m^2 is hard-coded in the Wurthwein 2021 ESM Section 4 control stream this model descends from (TVV1 = THETA(1)*(1+SCV1*(BSA-0.79))). IMPORTANT: the control stream contains no multiplication by BSA, so V, CLinitial and Qtr are ABSOLUTE quantities in L and L/day; the 'L/m^2' unit tags in Table S9 record the authors' convention of quoting the absolute value for a child with BSA = 1 m^2, not a per-square-metre normalisation. model() renormalises the linear factor to 1 at BSA = 1 m^2 so the ini() values are exactly the printed table entries.",
      source_name        = "BSA"
    ),
    AGE = list(
      description        = "Age at the time of the PEG-asparaginase administration",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "German/Czech non-high-risk median 4.7 years, range 1.06-18.0 (Wurthwein 2025 Table 1). Enters initial clearance as a hockey stick with the break point FIXED at 8 years: flat up to 8 years, rising linearly above it (Wurthwein 2021 ESM Section 4 control stream, two AGEGTP branches with the sub-8-year slope fixed to 0; ESM Section 2, 'only a minor increase in OFV for break points between 7 and 11 years; thus, the break point was fixed to 8 years'). Age at administration, not enrolment, which matters here because the R2-EA schedule runs from induction into maintenance dose M10 and so spans many months per subject. Wurthwein 2025 Section 3.2.3 confirms retention: 'The impact of age, sex, and anti-PEG IgM prior PIAd12 as covariates on CLinitial were confirmed.'",
      source_name        = "AGEGTP"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "German/Czech non-high-risk cohort 1119 male / 873 female (Wurthwein 2025 Table 1). Females have lower initial clearance -- ESM Table S9 footnote (d), 'fractional change in CLinitial for females' -- which the paper connects to the observation that 'In all trial groups, female patients showed higher DIPs compared with male patients' (Section 3.4.2). The source column is a 1/2 coded SEX (Wurthwein 2021 ESM Section 4: IF(SEX.EQ.1) CLSEX = 1 for males, IF(SEX.EQ.2) CLSEX = (1 + THETA(15)) for females), so SEXF = SEX - 1 preserves both the coefficient sign and the male reference category.",
      source_name        = "SEX (1 = male, 2 = female)"
    ),
    ABPEG_IGM = list(
      description        = "Pre-existing anti-polyethylene-glycol IgM antibody level measured in the sample drawn before the first PEG-asparaginase dose of induction (protocol IA day 12)",
      units              = "(assay-specific mean fluorescence intensity, linear scale)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "BASELINE, NOT TIME-VARYING: a single per-subject value from the sample taken at or within 3 days before the first induction dose, carried unchanged for the whole profile even though this model spans induction through maintenance (ESM Section 2.2 defines 'prior' as within 3 days before the first dose of the phase). Assayed by the plate-reader method of Khalil et al. 2022, which binds antibodies to immobilised methoxy-PEG chains and reports the duplicate mean fluorescence intensity; an arbitrary assay scale, not a mass concentration, and every subject has a measurable positive value. Anti-PEG antibodies are PRE-EXISTING (from environmental polyethylene-glycol exposure) rather than treatment-emergent, which is why they can act on the very first dose. Enters initial clearance as a HOCKEY STICK on the LOG scale with no effect at or below the cut point -- ESM Section 2.2 equation (7): 'For ABprior > cut point: CLinitial = theta_pop x (1 + theta_AB x (ABprior - cut point))'. The cut point is 1.3 on the log scale, equivalently 3.67 on the linear scale, and ESM Table S9 footnote (e) records it as held fixed ('> cut point (1.3 FIX)'). This model stores the LINEAR level and takes the logarithm inside model(), so the stored cut point is 3.67. The effect applies ONLY to the first administration in induction: footnote (e) reads 'following the 1st dose in induction', and Section 3.2.3 explains that no further antibody testing was done here because 'anti-PEG antibody levels prior to or after the first dose in re-induction were not identified as predictive covariates'. NAME PROVISIONAL: no anti-PEG antibody canonical existed in inst/references/covariate-columns.md when this model was written; the name is proposed there and awaits operator ratification (see the vignette's Assumptions and deviations section).",
      source_name        = "anti-PEG IgM prior PIAd12"
    ),
    OCC = list(
      description        = "Integer administration-occasion index; drives both the treatment-phase effects on initial clearance and volume and the inter-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (first induction dose, protocol IA day 12) is the model reference for every fractional change",
      notes              = "One occasion per PEG-asparaginase administration, matching the authors' own control-stream idiom of a single OCC column that simultaneously selects the phase covariate factors and the IOV etas (Wurthwein 2021 ESM Section 4), and their statement that 'definition of occasion always matched the grouping of administrations' (Wurthwein 2025 ESM Section 2.1). Mapping for the R2-EA schedule described in Wurthwein 2025 Section 2.1 and Figure 1: OCC = 1 protocol IA day 12 (induction dose 1, the reference); 2 protocol IA day 26 (induction dose 2); 3 protocol II day 8 (re-induction); 4, 5, 6 the three biweekly protocol II experimental-arm doses on days 22, 36 and 50; 7, 8, 9, 10, 11, 12 the six biweekly maintenance doses M5 through M10 (protocol days 64, 78, 92, 106, 120 and 134). Volume pools occasions 3-12 into a single factor while initial clearance groups them into three steps -- see the ini() comments. Records with OCC outside 1-12 receive the reference phase factors and no IOV contribution.",
      source_name        = "OCC / administration (PIAd12, PIAd26, PIId8, PII-ASP+ d22/d36/d50, M5-M10)"
    )
  )

  covariatesDataExcluded <- list(
    ABPEG_IGG = list(
      description        = "Pre-existing anti-polyethylene-glycol IgG antibody level before the first dose",
      units              = "(assay-specific mean fluorescence intensity, linear scale)",
      type               = "continuous",
      notes              = "Not tested in this arm. Wurthwein 2025 Section 3.2.3: 'As anti-PEG antibody levels prior to or after the first dose in re-induction were not identified as predictive covariates, the impact of these antibodies in R2-EA was not tested further.' In the parallel high-risk model the IgG isotype was screened and rejected as 'less pronounced' than IgM (ESM Section 2.2, Table S6)."
    ),
    ADA_POS = list(
      description        = "Anti-E. coli-asparaginase and anti-PEG-asparaginase antibody positivity (indirect ELISA, sample optical density over floating cut-off > 1.1)",
      units              = "(binary)",
      type               = "binary",
      notes              = "Screened across the trial and abandoned for lack of positive samples: after the filter criteria only 5 of 694 samples after protocol II day 8 were anti-E. coli-asparaginase positive and none were anti-PEG-asparaginase positive, so 'evaluation of these antibody levels as potential covariates in the popPK models was considered not to be possible' (Wurthwein 2025 Section 3.2.4, ESM Table S10). Distinct from ABPEG_IGM, which targets the polyethylene-glycol moiety rather than the asparaginase protein."
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      notes              = "Screened as an alternative body-size descriptor during development of the parent model and rejected in favour of BSA. Wurthwein 2021 ESM Section 2: 'Allometric weight scaling with fixed exponents (0.75 on clearance-terms, 1 on V), estimated exponents or linear scaling with body weight did not improve the model compared to BSA scaling.'"
    )
  )

  compartmentData <- list(
    central   = list(analyte = "PEG-asparaginase, fully PEGylated", units = "U", specimen = "serum", verified = TRUE),
    transit1  = list(analyte = "PEG-asparaginase, de-PEGylation step 1",  units = "U", specimen = "serum", verified = TRUE),
    transit2  = list(analyte = "PEG-asparaginase, de-PEGylation step 2",  units = "U", specimen = "serum", verified = TRUE),
    transit3  = list(analyte = "PEG-asparaginase, de-PEGylation step 3",  units = "U", specimen = "serum", verified = TRUE),
    transit4  = list(analyte = "PEG-asparaginase, de-PEGylation step 4",  units = "U", specimen = "serum", verified = TRUE),
    transit5  = list(analyte = "PEG-asparaginase, de-PEGylation step 5",  units = "U", specimen = "serum", verified = TRUE),
    transit6  = list(analyte = "PEG-asparaginase, de-PEGylation step 6",  units = "U", specimen = "serum", verified = TRUE),
    transit7  = list(analyte = "PEG-asparaginase, de-PEGylation step 7",  units = "U", specimen = "serum", verified = TRUE),
    transit8  = list(analyte = "PEG-asparaginase, de-PEGylation step 8",  units = "U", specimen = "serum", verified = TRUE),
    transit9  = list(analyte = "PEG-asparaginase, de-PEGylation step 9",  units = "U", specimen = "serum", verified = TRUE),
    transit10 = list(analyte = "PEG-asparaginase, de-PEGylation step 10", units = "U", specimen = "serum", verified = TRUE),
    transit11 = list(analyte = "PEG-asparaginase, de-PEGylation step 11", units = "U", specimen = "serum", verified = TRUE),
    transit12 = list(analyte = "PEG-asparaginase, de-PEGylation step 12", units = "U", specimen = "serum", verified = TRUE),
    transit13 = list(analyte = "PEG-asparaginase, de-PEGylation step 13", units = "U", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species             = "human",
    n_subjects          = 267L,
    n_studies           = 1L,
    n_observations      = 1794L,
    age_range           = "1.06-18.0 years across the German/Czech non-high-risk cohort; median 4.7",
    bsa_range           = "0.39-2.44 m^2 across the German/Czech non-high-risk cohort; median 0.74",
    disease_state       = "Newly diagnosed pediatric acute lymphoblastic leukemia enrolled in the AIEOP-BFM ALL 2009 trial (EudraCT 2007-004270-43, NCT01117441), non-high-risk, randomised in the R2 randomisation to the experimental arm PII-ASP+ receiving nine PEG-asparaginase doses in addition to the protocol II day 8 dose.",
    dose_range          = "2500 U/m^2 per dose as a 2-hour intravenous infusion, capped at an absolute maximum of 3750 U/dose. Twelve administrations per subject: protocol IA days 12 and 26, protocol II day 8, protocol II experimental arm days 22, 36 and 50, and maintenance M5-M10 on protocol days 64, 78, 92, 106, 120 and 134. Administered dose per m^2, median [range]: induction 1850 [620, 4850]; protocol II experimental arm through M10 1830 [625, 4080] (Table 1).",
    regions             = "Germany and Czech Republic.",
    bioanalytic_methods = "Asparaginase activity in serum by the aspartic acid beta-hydroxamate (AHA) test; lower limit of quantification 5 U/L, calibration ranges 5-100 and 100-1000 U/L. Anti-PEG IgM by the plate-reader method of Khalil et al. 2022 (immobilised methoxy-PEG, duplicate mean fluorescence intensity).",
    notes               = "267 R2-EA patients contributing 1794 administrations (Wurthwein 2025 Table 1, PIIASP+ d22-M10 row); as with the high-risk model, data from all induction patients were carried alongside for model stabilisation ('For model stabilization, data from all patients treated within induction together with data in the respective treatment phases were considered for model development'). Dose-normalized activity climbs from a median 1440 U/L on day 7 after the first re-induction dose to a plateau of about 2070 U/L after the last, which the authors read as approaching steady state (Section 3.2.3, ESM Table S3). Records indicating silent inactivation and every record at or after a hypersensitivity reaction were EXCLUDED before fitting, so this model describes standard elimination only. NONMEM 7.4.4, FOCE with interaction. The bootstrap used 500 rather than 1000 replicates, 89.6 percent successful, because 'Due to the large data set and the rather complex transit-model the run-time of this bootstrap analysis was extremely long: on the high-performance computer of the University Muenster this bootstrap analysis took 9 days' (ESM Table S9 footnote). IIV shrinkage 18.2 percent, condition number 198.8 -- the highest of the paper's models."
  )

  ini({
    # ---- Structural PK -----------------------------------------------------
    # ESM Table S9, 'Final popPK model' column (model 493158). Values are
    # quoted for a child with BSA = 1 m^2; the BSA factors in model() are
    # renormalised to 1 at BSA = 1 m^2 so these are used exactly as printed.
    lvc <- log(1.62);  label("Volume of the serum compartment shared by all 14 species, for a child with BSA 1 m^2 (L)")                    # ESM Table S9 (V 1.62, RSE 1.1%; bootstrap 1.62, 95% CI 1.59-1.66)
    lcl <- log(0.127); label("Initial clearance of fully PEGylated asparaginase, for a child with BSA 1 m^2 (L/day)")               # ESM Table S9 (CLinitial 0.127, RSE 1.0%; bootstrap 0.127, 95% CI 0.125-0.13)
    lq  <- log(0.822); label("Intercompartmental clearance driving the de-PEGylation transit chain, for a child with BSA 1 m^2 (L/day)")   # ESM Table S9 (Qrt 0.822, RSE 2.5%; bootstrap 0.820, 95% CI 0.78-0.857)

    # ---- Body surface area -------------------------------------------------
    e_bsa_vc   <- 1.57; label("Linear slope of BSA on V, centred on 0.79 m^2 (1/m^2)")                  # ESM Table S9 (FBSA on V 1.57, RSE 1.6%; bootstrap 1.56, 95% CI 1.51-1.61)
    e_bsa_cl_q <- 1.44; label("Linear slope of BSA on CLinitial and Qtr, centred on 0.79 m^2 (1/m^2)")  # ESM Table S9 (FBSA on CLinitial + Qtr 1.44, RSE 1.7%; bootstrap 1.44, 95% CI 1.4-1.49)

    # ---- Treatment-phase fractional changes on volume ----------------------
    # V resolves only three groups: the reference first induction dose, the
    # second induction dose, and everything from protocol II day 8 onward
    # pooled into one factor. ESM Table S8 records that 'Additional changes in
    # V did not improve the fit of the data'.
    e_piad26_vc   <- -0.175; label("Fractional change in V at the second induction dose, protocol IA day 26 (unitless)")                # ESM Table S9 (F(V PIAd26) -0.175, RSE 5.4%; bootstrap -0.176, 95% CI -0.195 to -0.156)
    e_piid8m10_vc <- -0.302; label("Fractional change in V from protocol II day 8 through maintenance dose M10 (unitless)")             # ESM Table S9 (F(V PIId8 - M10) -0.302, RSE 3.3%; bootstrap -0.302, 95% CI -0.323 to -0.282)

    # ---- Treatment-phase fractional changes on initial clearance -----------
    # Clearance falls in three steps across the ten biweekly doses. ESM
    # Table S8 shows how the grouping was chosen: model 493111 estimated a
    # separate factor for every administration, the confidence intervals
    # overlapped except between protocol II days 22 and 36, and the retained
    # model 493151 keeps exactly the two breaks whose intervals separated --
    # after the 2nd re-induction dose and after M6. The terminal -0.601 is the
    # largest accumulation effect in the paper and corresponds to the roughly
    # steady-state activity levels seen from M7 onward.
    e_piad26_cl        <- -0.0758; label("Fractional change in CLinitial at the second induction dose, protocol IA day 26 (unitless)")                      # ESM Table S9 (F(CLinitial PIAd26) -0.0758, RSE 15.6%; bootstrap -0.0762, 95% CI -0.0994 to -0.0504)
    e_piid8piiasp22_cl <- -0.377;  label("Fractional change in CLinitial from protocol II day 8 through the experimental-arm day 22 dose (unitless)")        # ESM Table S9 (F(CLinitial PIId8-PIIASP+d22) -0.377, RSE 2.5%; bootstrap -0.377, 95% CI -0.395 to -0.359)
    e_piiasp36m6_cl    <- -0.515;  label("Fractional change in CLinitial from the experimental-arm day 36 dose through maintenance dose M6 (unitless)")      # ESM Table S9 (F(CLinitial PIIASP+d36-M6) -0.515, RSE 2.6%; bootstrap -0.514, 95% CI -0.537 to -0.485)
    e_m7m10_cl         <- -0.601;  label("Fractional change in CLinitial from maintenance dose M7 through M10 (unitless)")                                  # ESM Table S9 (F(CLinitial M7-M10) -0.601, RSE 2.6%; bootstrap -0.599, 95% CI -0.629 to -0.568)

    # ---- Age, sex and pre-existing anti-PEG IgM on initial clearance -------
    e_age_cl  <- 0.0178;  label("Linear slope of age above 8 years on CLinitial (1/year)")  # ESM Table S9 (F(Age > 8 years) 0.0178, RSE 18.3%; bootstrap 0.0178, 95% CI 0.011-0.0238)
    e_sexf_cl <- -0.0651; label("Fractional change in CLinitial for females (unitless)")    # ESM Table S9 (F(Sex) -0.0651, RSE 16.6%; bootstrap -0.0651, 95% CI -0.0866 to -0.0423)

    # Hockey stick on the LOG antibody level, active only above the cut point
    # and only on the first induction dose -- ESM Section 2.2 equation (7) and
    # Table S9 footnote (e).
    e_abpeg_igm_cl <- 0.449; label("Linear slope of log anti-PEG IgM above the cut point on CLinitial at the first induction dose (per log unit)")  # ESM Table S9 (F(anti-PEG-IgMprior PIAd12) 0.449, RSE 14.5%; bootstrap 0.4427, 95% CI 0.3296-0.5789)

    # ---- Inter-individual variability --------------------------------------
    # Exponential IIV; the printed percentage is the log-normal CV, so
    # omega^2 = log(1 + CV^2). No IIV on V or Qtr.
    etalcl ~ 0.045955  # ESM Table S9 (IIV CLinitial 21.6%, RSE 3.4%, shrinkage 18.2%; bootstrap 21.5%, 95% CI 20.1-23.0) converted as log(1 + 0.216^2)

    # ---- Inter-occasion variability ----------------------------------------
    # On BOTH V and CLinitial, one occasion per administration (twelve here).
    # Occasions 2-12 carry their own etas with the variance fixed equal to the
    # occasion-1 estimate, encoding NONMEM's `$OMEGA BLOCK(1) SAME`.
    etaiov_vc_1  ~ 0.014085        # ESM Table S9 (IOV V 11.9%, RSE 7.9%; bootstrap 11.7%, 95% CI 9.5-13.7) converted as log(1 + 0.119^2)
    etaiov_vc_2  ~ fixed(0.014085) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_3  ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_4  ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_5  ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_6  ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_7  ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_8  ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_9  ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_10 ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_11 ~ fixed(0.014085) # SAME-equivalent
    etaiov_vc_12 ~ fixed(0.014085) # SAME-equivalent
    etaiov_cl_1  ~ 0.040558        # ESM Table S9 (IOV CLinitial 20.3%, RSE 3.2%; bootstrap 20.3%, 95% CI 19.0-21.8) converted as log(1 + 0.203^2)
    etaiov_cl_2  ~ fixed(0.040558) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3  ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_4  ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_5  ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_6  ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_7  ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_8  ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_9  ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_10 ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_11 ~ fixed(0.040558) # SAME-equivalent
    etaiov_cl_12 ~ fixed(0.040558) # SAME-equivalent

    # ---- Residual error ----------------------------------------------------
    # Combined proportional + additive SD, per the Wurthwein 2021 ESM
    # Section 4 $ERROR block (W = SQRT(prop^2*IPRED^2 + add^2), $SIGMA 1 FIX).
    propSd <- 0.209; label("Proportional residual error (fraction)")  # ESM Table S9 (Proportional error 20.9%, RSE 1.6%; bootstrap 20.8%, 95% CI 20.0-21.8)
    addSd  <- 15.9;  label("Additive residual error (U/L)")           # ESM Table S9 (Additive error 15.9 U/L, RSE 42.9%; bootstrap 15.7, 95% CI 6.7-29.5)
  })

  model({
    # ---- Administration occasion ------------------------------------------
    # OCC = 1 protocol IA day 12 (reference), 2 protocol IA day 26,
    # 3 protocol II day 8, 4-6 protocol II experimental arm days 22/36/50,
    # 7-12 maintenance M5-M10. Records outside 1-12 fall back to the reference
    # phase and contribute no IOV.
    occ1  <- (OCC == 1)
    occ2  <- (OCC == 2)
    occ3  <- (OCC == 3)
    occ4  <- (OCC == 4)
    occ5  <- (OCC == 5)
    occ6  <- (OCC == 6)
    occ7  <- (OCC == 7)
    occ8  <- (OCC == 8)
    occ9  <- (OCC == 9)
    occ10 <- (OCC == 10)
    occ11 <- (OCC == 11)
    occ12 <- (OCC == 12)

    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 + occ3 * etaiov_vc_3 +
      occ4 * etaiov_vc_4 + occ5 * etaiov_vc_5 + occ6 * etaiov_vc_6 +
      occ7 * etaiov_vc_7 + occ8 * etaiov_vc_8 + occ9 * etaiov_vc_9 +
      occ10 * etaiov_vc_10 + occ11 * etaiov_vc_11 + occ12 * etaiov_vc_12
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3 +
      occ4 * etaiov_cl_4 + occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6 +
      occ7 * etaiov_cl_7 + occ8 * etaiov_cl_8 + occ9 * etaiov_cl_9 +
      occ10 * etaiov_cl_10 + occ11 * etaiov_cl_11 + occ12 * etaiov_cl_12

    # ---- Body surface area -------------------------------------------------
    # Linear and centred on 0.79 m^2, renormalised to 1 at BSA = 1 m^2 so the
    # ini() values are the printed table entries unaltered.
    fbsa_vc <- (1 + e_bsa_vc   * (BSA - 0.79)) / (1 + e_bsa_vc   * 0.21)
    fbsa_cl <- (1 + e_bsa_cl_q * (BSA - 0.79)) / (1 + e_bsa_cl_q * 0.21)

    # ---- Age, sex and anti-PEG IgM on initial clearance --------------------
    fage_cl <- 1 + e_age_cl * (AGE - 8) * (AGE > 8)
    fsex_cl <- 1 + e_sexf_cl * SEXF

    # Hockey stick: no effect at or below the cut point. The stored covariate
    # is the LINEAR antibody level, so the cut point is 3.67 = exp(1.3); the
    # authors work on the log scale where it reads 1.3 ('1.3 FIX', ESM
    # Table S9 footnote (e)). Gated to the first induction dose only.
    ab_igm  <- log(ABPEG_IGM)
    fabm_cl <- 1 + e_abpeg_igm_cl * (ab_igm - 1.3) * (ab_igm > 1.3) * occ1

    # ---- Treatment-phase factors -------------------------------------------
    # V resolves three groups; CLinitial resolves four.
    fphase_vc <- 1 + occ2 * e_piad26_vc +
      (occ3 + occ4 + occ5 + occ6 + occ7 + occ8 + occ9 + occ10 + occ11 +
         occ12) * e_piid8m10_vc
    fphase_cl <- 1 + occ2 * e_piad26_cl +
      (occ3 + occ4) * e_piid8piiasp22_cl +
      (occ5 + occ6 + occ7 + occ8) * e_piiasp36m6_cl +
      (occ9 + occ10 + occ11 + occ12) * e_m7m10_cl

    # ---- Individual parameters ---------------------------------------------
    vc <- exp(lvc + iov_vc) * fbsa_vc * fphase_vc
    cl <- exp(lcl + etalcl + iov_cl) * fbsa_cl * fage_cl * fsex_cl *
      fabm_cl * fphase_cl
    q  <- exp(lq) * fbsa_cl

    # ---- De-PEGylation transit chain ---------------------------------------
    # Wurthwein 2021 Figure 1 and ESM Section 4 $DES. Fourteen compartments in
    # series sharing the single serum volume vc; drug steps down the chain at
    # Qtr/V (progressive hydrolysis of the PEG moiety) and every compartment
    # is eliminated at CLinitial/V. The terminal compartment additionally
    # loses drug at Qtr/V because the published simplification sets the
    # induced clearance CLinduced equal to Qtr, which is what makes apparent
    # clearance climb from CLinitial toward CLinitial + Qtr within a dosing
    # interval as mass migrates down the chain.
    kel <- cl / vc
    ktr <- q / vc

    d/dt(central)   <-              -(kel + ktr) * central
    d/dt(transit1)  <- ktr * central   - (kel + ktr) * transit1
    d/dt(transit2)  <- ktr * transit1  - (kel + ktr) * transit2
    d/dt(transit3)  <- ktr * transit2  - (kel + ktr) * transit3
    d/dt(transit4)  <- ktr * transit3  - (kel + ktr) * transit4
    d/dt(transit5)  <- ktr * transit4  - (kel + ktr) * transit5
    d/dt(transit6)  <- ktr * transit5  - (kel + ktr) * transit6
    d/dt(transit7)  <- ktr * transit6  - (kel + ktr) * transit7
    d/dt(transit8)  <- ktr * transit7  - (kel + ktr) * transit8
    d/dt(transit9)  <- ktr * transit8  - (kel + ktr) * transit9
    d/dt(transit10) <- ktr * transit9  - (kel + ktr) * transit10
    d/dt(transit11) <- ktr * transit10 - (kel + ktr) * transit11
    d/dt(transit12) <- ktr * transit11 - (kel + ktr) * transit12
    d/dt(transit13) <- ktr * transit12 - (kel + ktr) * transit13

    # ---- Observation --------------------------------------------------------
    # The AHA assay measures total catalytic asparaginase activity, so every
    # species in the chain contributes regardless of how far de-PEGylation has
    # progressed. Control stream $ERROR: TOT = A(1)+...+A(14); IPRED = TOT/V1.
    Cc <- (central + transit1 + transit2 + transit3 + transit4 + transit5 +
             transit6 + transit7 + transit8 + transit9 + transit10 +
             transit11 + transit12 + transit13) / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}

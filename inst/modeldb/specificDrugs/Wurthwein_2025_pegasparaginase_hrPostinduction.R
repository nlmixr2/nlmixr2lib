Wurthwein_2025_pegasparaginase_hrPostinduction <- function() {
  description <- "Fourteen-compartment de-PEGylation transit population PK model for intravenous PEGylated asparaginase (PEG-ASNase) in high-risk children with acute lymphoblastic leukemia across induction and the whole of high-risk post-induction therapy in the German/Czech group of the AIEOP-BFM ALL 2009 trial (Wurthwein 2025, ESM Table S7 model 354103). Extends the induction / re-induction model to seven administration groups: induction protocol IA days 12 and 26, the four weekly experimental protocol IB doses, the three high-risk blocks HR-1 to HR-3, and the three protocol III re-inductions. Asparaginase activity is carried by a chain of 14 serial compartments sharing one serum volume; drug moves down the chain with intercompartmental clearance Qtr, which mimics stepwise de-PEGylation, and every compartment is eliminated with the initial clearance CLinitial while the terminal compartment is eliminated with CLinitial + Qtr. Body surface area enters volume and the two clearance terms linearly, centred on 0.79 m^2. Initial clearance rises linearly with age above 8 years, is lower in females, and rises with a pre-existing anti-polyethylene-glycol IgM antibody level above a hockey-stick cut point, the antibody effect acting only on the first induction dose. Initial clearance and volume fall stepwise across the protocol, by up to 47.3 and 27.2 percent respectively relative to the first induction dose. Inter-individual variability on initial clearance, inter-occasion variability on both initial clearance and volume, and combined proportional plus additive residual error."
  reference   <- "Wurthwein G, Siebel C, Lanvers-Kaminsky C, Smisek P, Nath CE, Matteo C, Rizzari C, Schrappe M, Boos J. PEGylated Asparaginase in Children with Acute Lymphoblastic Leukemia Treated within the AIEOP-BFM ALL 2009 Trial: Population Pharmacokinetics and Drug Exposure. Eur J Drug Metab Pharmacokinet. 2025;50(6):683-696. doi:10.1007/s13318-025-00962-3"
  vignette    <- "Wurthwein_2025_pegasparaginase"
  units       <- list(time = "day", dosing = "U", concentration = "U/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed by the Mosteller formula (Wurthwein 2025 Methods, Pharmacokinetic Model, citing reference 14). German/Czech high-risk median 0.99 m^2, range 0.41-2.42; whole German/Czech cohort median 0.78 m^2 (Wurthwein 2025 Table 1). BSA enters volume and the two clearance terms LINEARLY and CENTRED, not allometrically -- ESM Table S7 footnote (a), 'linear increase in V / CLinitial / Qtr with BSA (centred on the median)'. The centring constant 0.79 m^2 is hard-coded in the Wurthwein 2021 ESM Section 4 control stream this model descends from (TVV1 = THETA(1)*(1+SCV1*(BSA-0.79))). IMPORTANT: the control stream contains no multiplication by BSA, so V, CLinitial and Qtr are ABSOLUTE quantities in L and L/day; the 'L/m^2' unit tags in Table S7 record the authors' convention of quoting the absolute value for a child with BSA = 1 m^2 ('Typical values for V, CLinitial and Qtr are reported for a child with BSA = 1 m^2 for better comparison'), not a per-square-metre normalisation. model() renormalises the linear factor to 1 at BSA = 1 m^2 so the ini() values are exactly the printed table entries.",
      source_name        = "BSA"
    ),
    AGE = list(
      description        = "Age at the time of the PEG-asparaginase administration",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "German/Czech high-risk median 8.0 years, range 1.07-18.0 (Wurthwein 2025 Table 1); high-risk patients 'tended to be older than non-HR patients', so this cohort straddles the 8-year break point far more than the whole-trial cohort does. Enters initial clearance as a hockey stick with the break point FIXED at 8 years: flat up to 8 years, rising linearly above it (Wurthwein 2021 ESM Section 4 control stream, two AGEGTP branches with the sub-8-year slope fixed to 0; ESM Section 2, 'only a minor increase in OFV for break points between 7 and 11 years; thus, the break point was fixed to 8 years'). Age at administration, so time-varying across a protocol that spans years. Wurthwein 2025 Section 3.2.1 confirms retention here: 'forward inclusion and backward elimination confirmed the inclusion of covariates age and sex'.",
      source_name        = "AGEGTP"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "German/Czech high-risk cohort 349 male / 194 female (Wurthwein 2025 Table 1); 'HR patients were more often male'. Females have lower initial clearance -- ESM Table S7 footnote (d), 'fractional change in CLinitial for females'. The source column is a 1/2 coded SEX (Wurthwein 2021 ESM Section 4: IF(SEX.EQ.1) CLSEX = 1 for males, IF(SEX.EQ.2) CLSEX = (1 + THETA(15)) for females), so SEXF = SEX - 1 preserves both the coefficient sign and the male reference category.",
      source_name        = "SEX (1 = male, 2 = female)"
    ),
    ABPEG_IGM = list(
      description        = "Pre-existing anti-polyethylene-glycol IgM antibody level measured in the sample drawn before the first PEG-asparaginase dose of induction (protocol IA day 12)",
      units              = "(assay-specific mean fluorescence intensity, linear scale)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "BASELINE, NOT TIME-VARYING: a single per-subject value from the sample taken at or within 3 days before the first induction dose, carried unchanged for the whole profile (ESM Section 2.2 defines 'prior' as 'analyzed antibody levels prior (<=3 days) first dose of each exposure phase'). Assayed by the plate-reader method of Khalil et al. 2022, which binds antibodies to immobilised methoxy-PEG chains and reports the duplicate mean fluorescence intensity; the readout is an arbitrary assay scale, not a mass concentration, and every subject has a measurable positive value. Anti-PEG antibodies are PRE-EXISTING (from environmental polyethylene-glycol exposure) rather than treatment-emergent, which is why they can act on the very first dose. Enters initial clearance as a HOCKEY STICK on the LOG scale, with no effect at or below the cut point -- ESM Section 2.2 equation (7): 'For ABprior > cut point: CLinitial = theta_pop x (1 + theta_AB x (ABprior - cut point))', with ABprior the log-transformed antibody level. The cut point is 1.3 on the log scale, equivalently 3.67 on the linear scale; ESM Table S5 prints both, and ESM Table S9 footnote (e) records that it was held fixed ('> cut point (1.3 FIX)'). This model stores the LINEAR level and takes the logarithm inside model(), so the stored cut point is 3.67. The effect applies ONLY to the first administration in induction: model 354103 is described in ESM Section 2.2 as 'IgMprior effect only on 1st admin. in induction', and Table S6 shows the corresponding HR-1 coefficient fixed to zero after it estimated with 54-68 percent RSE. 30.1 percent of the 1444 induction patients were above the Khalil cut point of 2, but only 10.1 percent above the modelled 3.67 (ESM Table S5). Founding example of the ABPEG_<isotype> canonical family in inst/references/covariate-columns.md, which is separate from the ADA_ family because the antigen is the polyethylene-glycol carrier rather than the drug.",
      source_name        = "anti-PEG IgM prior PIAd12"
    ),
    OCC = list(
      description        = "Integer administration-group index; drives both the treatment-phase effects on initial clearance and volume and the inter-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (first induction dose, protocol IA day 12) is the model reference for every fractional change",
      notes              = "One occasion per administration GROUP, not per dose: Wurthwein 2025 ESM Section 2.1 states 'definition of occasion always matched the grouping of administrations', and the authors' own control-stream idiom is a single OCC column that simultaneously selects the phase covariate factors and the IOV etas (Wurthwein 2021 ESM Section 4). Mapping for this model, following the AIEOP-BFM ALL 2009 outline in Figure 1: OCC = 1 protocol IA day 12 (induction dose 1, the reference); 2 protocol IA day 26 (induction dose 2); 3 protocol IB experimental arm, the four weekly doses on days 40, 47, 54 and 61 (all four share one factor); 4 high-risk block HR-1 (day 6); 5 high-risk block HR-2; 6 high-risk block HR-3; 7 protocol III, all three re-inductions 1.PIII, 2.PIII and 3.PIII on day 1 of each (all three share one factor). Note that CLinitial resolves all seven groups separately while V pools them into only four -- see the ini() comments. Records with OCC outside 1-7 receive the reference phase factors and no IOV contribution.",
      source_name        = "OCC / treatment phase (PIAd12, PIAd26, PIB-ASP+, HR-1, HR-2, HR-3, 1.-3.PIII)"
    )
  )

  covariatesDataExcluded <- list(
    ABPEG_IGG = list(
      description        = "Pre-existing anti-polyethylene-glycol IgG antibody level before the first dose of the exposure phase",
      units              = "(assay-specific mean fluorescence intensity, linear scale)",
      type               = "continuous",
      notes              = "Screened alongside the IgM isotype with the same hockey-stick form and the Khalil cut point of 2.079 on the log scale, and NOT retained: ESM Table S6 step 1 shows the best IgG model (353003) reaching BIC 73259 against 73239 for the retained IgM model 354103, and ESM Section 2.2 concludes 'As already shown in Siebel et al. the influence of anti-PEG IgGprior was less pronounced.'"
    ),
    ADA_POS = list(
      description        = "Anti-E. coli-asparaginase and anti-PEG-asparaginase antibody positivity (indirect ELISA, sample optical density over floating cut-off > 1.1)",
      units              = "(binary)",
      type               = "binary",
      notes              = "Screened as potential covariates and abandoned for lack of positive samples. Wurthwein 2025 Section 3.2.4 counts them: in induction 10 of 1011 samples before the first dose were anti-E. coli-asparaginase positive and 31 of 1011 anti-PEG-asparaginase positive, with even fewer in later phases, so 'evaluation of these antibody levels as potential covariates in the popPK models was considered not to be possible' (ESM Table S10). Distinct from ABPEG_IGM, which targets the polyethylene-glycol moiety rather than the asparaginase protein."
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      notes              = "Screened as an alternative body-size descriptor during development of the parent model and rejected in favour of BSA. Wurthwein 2021 ESM Section 2: 'Allometric weight scaling with fixed exponents (0.75 on clearance-terms, 1 on V), estimated exponents or linear scaling with body weight did not improve the model compared to BSA scaling', and for the highly correlated pair (r = 0.988) 'the covariate with the best improvement in OFV was retained'."
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
    n_subjects          = 1444L,
    n_studies           = 1L,
    age_range           = "1.06-18.0 years across the German/Czech cohort; high-risk subgroup median 8.0 years, range 1.07-18.0",
    bsa_range           = "0.39-2.44 m^2 across the German/Czech cohort; high-risk subgroup median 0.99 m^2, range 0.41-2.42",
    disease_state       = "Newly diagnosed pediatric acute lymphoblastic leukemia enrolled in the AIEOP-BFM ALL 2009 trial (EudraCT 2007-004270-43, NCT01117441). High-risk patients receive PEG-asparaginase in the high-risk blocks and protocol III re-inductions in addition to induction; those randomised to the protocol IB experimental arm also receive four extra weekly doses.",
    dose_range          = "2500 U/m^2 per dose as a 2-hour intravenous infusion, capped at an absolute maximum of 3750 U/dose; proportionally more high-risk patients hit the absolute cap. Administered dose per m^2, median [range] by phase: induction 2500 [720, 5300]; protocol IB experimental arm days 40-61 2130 [1150, 3750]; high-risk blocks HR-1 to HR-3 2240 [275, 4800]; protocol III re-inductions 2300 [275, 4580] (Table 1).",
    regions             = "Germany and Czech Republic.",
    bioanalytic_methods = "Asparaginase activity in serum by the aspartic acid beta-hydroxamate (AHA) test; lower limit of quantification 5 U/L, calibration ranges 5-100 and 100-1000 U/L. Anti-PEG IgG and IgM by the plate-reader method of Khalil et al. 2022 (immobilised methoxy-PEG, duplicate mean fluorescence intensity).",
    notes               = "The 1444 subjects are the induction / protocol IB exposure-phase analysis set of ESM Table S5; the later exposure phases contribute far fewer patients (HR-1 to 1.PIII 115, 2.PIII 85, 3.PIII 71), which is the stated reason no antibody effect could be estimated beyond induction (Section 3.2.2, Section 4.2). Model built stepwise: ESM Table S4 shows induction plus protocol IB first (retaining model 333105), then the high-risk blocks (334004), then protocol III (335000), with covariates removed at the start and re-tested at the end because 'the influences of the treatment phases or repetitive doses within a phase were expected to be much higher than the impact of these covariates' (Methods, Pharmacokinetic Model). Records indicating silent inactivation and every record at or after a hypersensitivity reaction were EXCLUDED before fitting, so this model describes standard elimination only; the authors note this is a material exclusion here because 'the proportion of samples positive for anti-PEG antibodies was particularly large among patients developing an inactivation reaction during the HR courses' (Section 4.2). A further 11 samples with absolute conditional weighted residuals above 4 were excluded before the antibody analysis (ESM Section 2.1). NONMEM 7.4.4, FOCE with interaction; 1000-replicate bootstrap, 99.9 percent successful; IIV shrinkage 22.9 percent, condition number 101.18 (ESM Table S7)."
  )

  ini({
    # ---- Structural PK -----------------------------------------------------
    # ESM Table S7, 'Final covariate model' column (model 354103). Values are
    # quoted for a child with BSA = 1 m^2; the BSA factors in model() are
    # renormalised to 1 at BSA = 1 m^2 so these are used exactly as printed.
    lvc <- log(1.73);  label("Volume of the serum compartment shared by all 14 species, for a child with BSA 1 m^2 (L)")                    # ESM Table S7 (V 1.73, RSE 1.4%; bootstrap 1.73, 95% CI 1.66-1.80)
    lcl <- log(0.127); label("Initial clearance of fully PEGylated asparaginase, for a child with BSA 1 m^2 (L/day)")               # ESM Table S7 (CLinitial 0.127, RSE 1.4%; bootstrap 0.127, 95% CI 0.122-0.132)
    lq  <- log(0.960); label("Intercompartmental clearance driving the de-PEGylation transit chain, for a child with BSA 1 m^2 (L/day)")   # ESM Table S7 (Qrt 0.960, RSE 2.2%; bootstrap 0.958, 95% CI 0.908-1.01)

    # ---- Body surface area -------------------------------------------------
    e_bsa_vc   <- 1.61; label("Linear slope of BSA on V, centred on 0.79 m^2 (1/m^2)")                  # ESM Table S7 (FBSA on V 1.61, RSE 2.2%; bootstrap 1.61, 95% CI 1.54-1.68)
    e_bsa_cl_q <- 1.48; label("Linear slope of BSA on CLinitial and Qtr, centred on 0.79 m^2 (1/m^2)")  # ESM Table S7 (FBSA on CLinitial + Qtr 1.48, RSE 2.4%; bootstrap 1.48, 95% CI 1.40-1.55)

    # ---- Treatment-phase fractional changes on volume ----------------------
    # V resolves only FOUR groups: the reference first induction dose, then
    # protocol IA day 26 pooled with the whole protocol IB experimental arm,
    # then the three high-risk blocks pooled, then the three protocol III
    # re-inductions pooled. ESM Table S4 records why -- modelling further
    # per-dose changes in V 'did not improve the fit of the data'.
    e_piad26pibasp_vc <- -0.159; label("Fractional change in V from protocol IA day 26 through the protocol IB experimental arm (unitless)")  # ESM Table S7 (F(V PIAd26 - PIB-ASP+) -0.159, RSE 6.6%; bootstrap -0.159, 95% CI -0.179 to -0.139)
    e_hr_vc           <- -0.136; label("Fractional change in V across the high-risk blocks HR-1 to HR-3 (unitless)")                          # ESM Table S7 (F(V HR-1 - HR-3) -0.136, RSE 12.4%; bootstrap -0.139, 95% CI -0.174 to -0.105)
    e_piii_vc         <- -0.272; label("Fractional change in V across the protocol III re-inductions 1.PIII to 3.PIII (unitless)")            # ESM Table S7 (F(V 1. - 3.PIII) -0.272, RSE 5.3%; bootstrap -0.273, 95% CI -0.302 to -0.244)

    # ---- Treatment-phase fractional changes on initial clearance -----------
    # CLinitial resolves all SEVEN groups. The high-risk-block effects are the
    # interesting ones: HR-1 came back at -0.031 with 70.4 percent RSE
    # (ESM Table S4 model 334003), so the authors FIXED it to zero, i.e. HR-1
    # clearance equals the first-induction-dose reference. HR-2 then drops and
    # HR-3 partially rebounds -- a non-monotone pattern the authors keep
    # because it fits, while cautioning that these phase terms 'have to be
    # regarded as a kind of dummy variables' standing in for unknown
    # physiology (Section 4.1).
    e_piad26_cl <- -0.110;    label("Fractional change in CLinitial at the second induction dose, protocol IA day 26 (unitless)")   # ESM Table S7 (F(CLinitial PIAd26) -0.110, RSE 12.5%; bootstrap -0.109, 95% CI -0.135 to -0.079)
    e_pibasp_cl <- -0.421;    label("Fractional change in CLinitial across the protocol IB experimental arm, days 40-61 (unitless)")  # ESM Table S7 (F(CLinitial PIB-ASP+) -0.421, RSE 6.3%; bootstrap -0.420, 95% CI -0.471 to -0.370)
    e_hr1_cl    <- fixed(0);  label("Fractional change in CLinitial in high-risk block HR-1 (unitless)")                            # ESM Table S7 (F(CLinitial HR-1) 0, footnote (f) 'the fractional change for HR-1 was fixed to 0'; freely estimated it was -0.031 with 70.4% RSE, ESM Table S4 model 334003)
    e_hr2_cl    <- -0.287;    label("Fractional change in CLinitial in high-risk block HR-2 (unitless)")                            # ESM Table S7 (F(CLinitial HR-2) -0.287, RSE 11.5%; bootstrap -0.289, 95% CI -0.354 to -0.222)
    e_hr3_cl    <- -0.185;    label("Fractional change in CLinitial in high-risk block HR-3 (unitless)")                            # ESM Table S7 (F(CLinitial HR-3) -0.185, RSE 22.9%; bootstrap -0.184, 95% CI -0.270 to -0.098)
    e_piii_cl   <- -0.473;    label("Fractional change in CLinitial across the protocol III re-inductions 1.PIII to 3.PIII (unitless)")  # ESM Table S7 (F(CLinitial 1. - 3.PIII) -0.473, RSE 4.3%; bootstrap -0.474, 95% CI -0.512 to -0.435)

    # ---- Age, sex and pre-existing anti-PEG IgM on initial clearance -------
    e_age_cl  <- 0.014;  label("Linear slope of age above 8 years on CLinitial (1/year)")  # ESM Table S7 (F(Age > 8 years) 0.014, RSE 33.5%; bootstrap 0.014, 95% CI 0.005-0.024)
    e_sexf_cl <- -0.060; label("Fractional change in CLinitial for females (unitless)")    # ESM Table S7 (F(Sex) -0.060, RSE 28.9%; bootstrap -0.061, 95% CI -0.093 to -0.028)

    # Hockey stick on the LOG antibody level, active only above the cut point
    # and only on the first induction dose -- ESM Section 2.2 equation (7)
    # and Table S6 model 354103.
    e_abpeg_igm_cl <- 0.396; label("Linear slope of log anti-PEG IgM above the cut point on CLinitial at the first induction dose (per log unit)")  # ESM Table S7 (F(anti-PEG-IgMprior PIAd12) 0.396, RSE 16.4%; bootstrap 0.401, 95% CI 0.271-0.534)

    # ---- Inter-individual variability --------------------------------------
    # Exponential IIV; the printed percentage is the log-normal CV, so
    # omega^2 = log(1 + CV^2). No IIV on V or Qtr.
    etalcl ~ 0.067649  # ESM Table S7 (IIV CLinitial 26.4%, RSE 4.4%, shrinkage 22.9%; bootstrap 26.3%, 95% CI 24.0-28.7) converted as log(1 + 0.264^2)

    # ---- Inter-occasion variability ----------------------------------------
    # On BOTH V and CLinitial, one occasion per administration group (seven
    # here). Occasions 2-7 carry their own etas with the variance fixed equal
    # to the occasion-1 estimate, encoding NONMEM's `$OMEGA BLOCK(1) SAME`.
    etaiov_vc_1 ~ 0.008251        # ESM Table S7 (IOV V 9.1%, RSE 9.7%; bootstrap 9.0%, 95% CI 7.2-10.9) converted as log(1 + 0.091^2)
    etaiov_vc_2 ~ fixed(0.008251) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_3 ~ fixed(0.008251) # SAME-equivalent
    etaiov_vc_4 ~ fixed(0.008251) # SAME-equivalent
    etaiov_vc_5 ~ fixed(0.008251) # SAME-equivalent
    etaiov_vc_6 ~ fixed(0.008251) # SAME-equivalent
    etaiov_vc_7 ~ fixed(0.008251) # SAME-equivalent
    etaiov_cl_1 ~ 0.049805        # ESM Table S7 (IOV CLinitial 22.6%, RSE 4.6%; bootstrap 22.5%, 95% CI 20.5-24.5) converted as log(1 + 0.226^2)
    etaiov_cl_2 ~ fixed(0.049805) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fixed(0.049805) # SAME-equivalent
    etaiov_cl_4 ~ fixed(0.049805) # SAME-equivalent
    etaiov_cl_5 ~ fixed(0.049805) # SAME-equivalent
    etaiov_cl_6 ~ fixed(0.049805) # SAME-equivalent
    etaiov_cl_7 ~ fixed(0.049805) # SAME-equivalent

    # ---- Residual error ----------------------------------------------------
    # Combined proportional + additive SD, per the Wurthwein 2021 ESM
    # Section 4 $ERROR block (W = SQRT(prop^2*IPRED^2 + add^2), $SIGMA 1 FIX).
    propSd <- 0.208; label("Proportional residual error (fraction)")  # ESM Table S7 (Proportional error 20.8%, RSE 2.0%; bootstrap 20.7%, 95% CI 19.8-21.7)
    addSd  <- 7.99;  label("Additive residual error (U/L)")           # ESM Table S7 (Additive error 7.99 U/L, RSE 27.5%; bootstrap 7.91, 95% CI 4.36-17.7)
  })

  model({
    # ---- Administration group ----------------------------------------------
    # OCC = 1 protocol IA day 12 (reference), 2 protocol IA day 26,
    # 3 protocol IB experimental arm (days 40/47/54/61), 4 HR-1, 5 HR-2,
    # 6 HR-3, 7 protocol III (1./2./3.PIII). Records outside 1-7 fall back to
    # the reference phase and contribute no IOV.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)
    occ7 <- (OCC == 7)

    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 + occ3 * etaiov_vc_3 +
      occ4 * etaiov_vc_4 + occ5 * etaiov_vc_5 + occ6 * etaiov_vc_6 +
      occ7 * etaiov_vc_7
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3 +
      occ4 * etaiov_cl_4 + occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6 +
      occ7 * etaiov_cl_7

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
    # authors work on the log scale where it reads 1.3 (ESM Table S5, Table S9
    # footnote (e) '1.3 FIX'). Gated to the first induction dose only.
    ab_igm  <- log(ABPEG_IGM)
    fabm_cl <- 1 + e_abpeg_igm_cl * (ab_igm - 1.3) * (ab_igm > 1.3) * occ1

    # ---- Treatment-phase factors -------------------------------------------
    # V pools the seven groups into four; CLinitial resolves all seven.
    fphase_vc <- 1 + (occ2 + occ3) * e_piad26pibasp_vc +
      (occ4 + occ5 + occ6) * e_hr_vc + occ7 * e_piii_vc
    fphase_cl <- 1 + occ2 * e_piad26_cl + occ3 * e_pibasp_cl +
      occ4 * e_hr1_cl + occ5 * e_hr2_cl + occ6 * e_hr3_cl + occ7 * e_piii_cl

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

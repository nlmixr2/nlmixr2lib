Wurthwein_2025_pegasparaginase_germanCzech <- function() {
  description <- "Fourteen-compartment de-PEGylation transit population PK model for intravenous PEGylated asparaginase (PEG-ASNase) in children with acute lymphoblastic leukemia treated in induction (protocol IA days 12 and 26) and re-induction (protocol II day 8) in the German/Czech group of the AIEOP-BFM ALL 2009 trial (Wurthwein 2025, ESM Table S13 model 223015). Asparaginase activity is carried by a chain of 14 serial compartments sharing one serum volume; drug moves down the chain with intercompartmental clearance Qtr, which mimics stepwise de-PEGylation, and every compartment is eliminated with the initial clearance CLinitial while the terminal compartment is eliminated with CLinitial + Qtr. Apparent clearance therefore rises within a dosing interval from CLinitial toward CLinitial + Qtr as the enzyme de-PEGylates. Body surface area enters volume and the two clearance terms linearly, centred on the model-building median of 0.79 m^2. Initial clearance rises linearly with age above 8 years and is lower in females. Initial clearance and volume drop stepwise between administrations, encoded as fractional changes against the first induction dose. Inter-individual variability on initial clearance, inter-occasion variability on both initial clearance and volume, and combined proportional plus additive residual error."
  reference   <- "Wurthwein G, Siebel C, Lanvers-Kaminsky C, Smisek P, Nath CE, Matteo C, Rizzari C, Schrappe M, Boos J. PEGylated Asparaginase in Children with Acute Lymphoblastic Leukemia Treated within the AIEOP-BFM ALL 2009 Trial: Population Pharmacokinetics and Drug Exposure. Eur J Drug Metab Pharmacokinet. 2025;50(6):683-696. doi:10.1007/s13318-025-00962-3"
  vignette    <- "Wurthwein_2025_pegasparaginase"
  units       <- list(time = "day", dosing = "U", concentration = "U/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed by the Mosteller formula (Wurthwein 2025 Methods, Pharmacokinetic Model, citing reference 14). German/Czech median 0.78 m^2, range 0.39-2.44 (Wurthwein 2025 Table 1). BSA enters volume and the two clearance terms LINEARLY and CENTRED, not allometrically: the Wurthwein 2021 ESM Section 4 control stream (the direct predecessor of this model, reference 7 of the 2025 supplement) writes TVV1 = THETA(1)*(1+SCV1*(BSA-0.79)) and TVCLP = THETA(2)*(1+SCCL*(BSA-0.79)), with the same SCCL slope applied to Qtr. The centring constant 0.79 m^2 is the model-building dataset median and is hard-coded in that control stream. Wurthwein 2025 ESM Table S13 footnote (a) confirms the form for the present fit: 'linear increase in V / CLinitial / Qtr with BSA (centred on the median)'. IMPORTANT: the control stream contains no multiplication by BSA anywhere, so V, CLinitial and Qtr are ABSOLUTE quantities in L and L/day; the 'L/m^2' and 'L/day/m^2' unit tags printed in Table S13 record the authors' convention of quoting the absolute value for a child with BSA = 1 m^2 ('Typical values for V, CLinitial and Qtr are reported for a child with BSA = 1 m^2 for better comparison'), NOT a per-square-metre normalisation. The model() code therefore renormalises the linear factor to 1 at BSA = 1 m^2 so the ini() values are exactly the printed table entries.",
      source_name        = "BSA"
    ),
    AGE = list(
      description        = "Age at the time of the PEG-asparaginase administration",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "German/Czech median 5.1 years, range 1.06-18.0 (Wurthwein 2025 Table 1). Enters initial clearance as a hockey stick with the break point FIXED at 8 years: clearance is flat up to 8 years and rises linearly above it. The Wurthwein 2021 ESM Section 4 control stream encodes this as two branches, IF(AGEGTP.LE.8) CLAGEGTP = (1 + THETA(13)*(AGEGTP-8)) with THETA(13) fixed to 0, and IF(AGEGTP.GT.8) CLAGEGTP = (1 + THETA(14)*(AGEGTP-8)). The 2021 ESM Section 2 records why the break point is fixed rather than estimated: 'Further analyses indicated only a minor increase in OFV for break points between 7 and 11 years; thus, the break point was fixed to 8 years.' The control stream's column name AGEGTP is age at the time of administration, not at enrolment, so this covariate is time-varying across a multi-year protocol.",
      source_name        = "AGEGTP"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "German/Czech cohort 1468 male / 1067 female (Wurthwein 2025 Table 1). Females have lower initial clearance: Table S13 footnote (d) reads 'fractional change in CLinitial for females' and the estimate is negative. The source column is NOT SEXF: the Wurthwein 2021 ESM Section 4 control stream uses a 1/2 coded SEX column with IF(SEX.EQ.1) CLSEX = 1 (male, the reference) and IF(SEX.EQ.2) CLSEX = (1 + THETA(15)) (female). Mapping to the canonical column is SEXF = SEX - 1, which preserves both the sign of the coefficient and the male reference category, so the printed estimate is used unchanged.",
      source_name        = "SEX (1 = male, 2 = female)"
    ),
    OCC = list(
      description        = "Integer administration-occasion index; drives both the treatment-phase effects on initial clearance and volume and the inter-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (first induction dose, protocol IA day 12) is the model reference for every fractional change",
      notes              = "One occasion per PEG-asparaginase administration, exactly as the authors encode it. The Wurthwein 2021 ESM Section 4 control stream carries a single OCC column that simultaneously selects the treatment-phase covariate factors (V1KOV, CLKOV) and the IOV etas (IOVV1, IOVCL), with the header comment ';; OCC: 2=PIA d12, 3=PIA d26, 5=PII d8'. This model renumbers those three occasions to the 1..N form the register prescribes for OCC: OCC = 1 is protocol IA day 12 (induction dose 1, the reference), OCC = 2 is protocol IA day 26 (induction dose 2), OCC = 3 is protocol II day 8 (re-induction). Wurthwein 2025 ESM Section 2.1 states the occasion definition for the extended models in the same terms: 'Here, definition of occasion always matched the grouping of administrations.' Records with OCC outside 1-3 receive the reference phase factors and no IOV contribution.",
      source_name        = "OCC (2 = PIA d12, 3 = PIA d26, 5 = PII d8)"
    )
  )

  # Covariates screened by Wurthwein 2025 but NOT retained in THIS model.
  # Documented for provenance only; deliberately absent from model().
  covariatesDataExcluded <- list(
    ABPEG_IGM = list(
      description        = "Pre-existing anti-polyethylene-glycol IgM antibody level in the sample drawn before the first induction dose",
      units              = "(assay-specific mean fluorescence intensity, linear scale)",
      type               = "continuous",
      notes              = "Retained in the two extended German/Czech models of this paper (Wurthwein_2025_pegasparaginase_hrPostinduction.R from ESM Table S7 and Wurthwein_2025_pegasparaginase_r2ea.R from ESM Table S9) but NOT part of the induction / re-induction model in ESM Table S13, whose covariate list is BSA, age and sex only. Wurthwein 2025 Section 3.4.2 gives the authors' reason for being content without it: 'there was almost no difference in DIPs derived from the popPK models without versus with anti-PEG IgM prior PIAd12 as covariate'."
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      notes              = "Screened as an alternative body-size descriptor and rejected in favour of BSA. Wurthwein 2021 ESM Section 2 records both the screening and the rule used to choose between them: 'Allometric weight scaling with fixed exponents (0.75 on clearance-terms, 1 on V), estimated exponents or linear scaling with body weight did not improve the model compared to BSA scaling', and 'For the highly correlated covariates body weight and BSA (r=0.988), the covariate with the best improvement in OFV was retained in the model.'"
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
    n_subjects          = 2535L,
    n_studies           = 1L,
    n_observations      = 17221L,
    age_range           = "1.06-18.0 years; median 5.1",
    bsa_range           = "0.39-2.44 m^2; median 0.78",
    sex_female_pct      = 42.1,
    race_ethnicity      = "Not reported (German and Czech trial sites of the AIEOP-BFM ALL 2009 trial).",
    disease_state       = "Newly diagnosed pediatric acute lymphoblastic leukemia enrolled in the AIEOP-BFM ALL 2009 trial (EudraCT 2007-004270-43, NCT01117441) between 1 June 2010 and 28 February 2017.",
    dose_range          = "2500 U/m^2 per dose as a 2-hour intravenous infusion, capped at an absolute maximum of 3750 U/dose (applied for about 12 percent of administrations). This model covers the two induction doses on protocol IA days 12 and 26 and the single re-induction dose on protocol II day 8. Administered dose per m^2, median [range]: induction 1940 [620, 5300]; re-induction 1880 [713, 4500] (Table 1).",
    regions             = "Germany and Czech Republic.",
    bioanalytic_methods = "Asparaginase activity in serum by the aspartic acid beta-hydroxamate (AHA) test; lower limit of quantification 5 U/L, calibration ranges 5-100 and 100-1000 U/L. Samples below the limit of quantification (1.9 percent) were excluded from the fit (ESM Section 1.2).",
    notes               = "Demographics from Wurthwein 2025 Table 1 (German/Czech, All column); the female percentage is 1067/2535. Samples were scheduled before each treatment phase and 7 and 14 days after each dose, so 85.5-93.3 percent of the analysed samples fall in the 7 +/- 1 or 14 +/- 1 day windows (ESM Table S2) and the design is sparse rather than a classical PK sampling scheme. Records indicating silent inactivation (asparaginase below 100 U/L within 8 days and/or undetectable within 15 days, per van der Sluis 2016) and every record at or after a hypersensitivity reaction were EXCLUDED before fitting, so this model describes standard elimination only and must not be used to predict the pharmacokinetics of an inactivating patient (Wurthwein 2025 Section 4.3). NONMEM 7.4.4, FOCE with interaction; evaluated by goodness-of-fit plots, prediction-corrected VPC (n = 1000) and a 1000-replicate bootstrap of which 90.6 percent were successful; IIV shrinkage 23.3 percent, condition number 111.81 (ESM Table S13 footnote (e))."
  )

  ini({
    # ---- Structural PK -----------------------------------------------------
    # ESM Table S13, 'German/Czech data' Estimate column (model 223015).
    # The table quotes V, CLinitial and Qtr for a child with BSA = 1 m^2; the
    # BSA factors in model() are renormalised to 1 at BSA = 1 m^2 so these are
    # used exactly as printed. See covariateData$BSA on why the printed
    # 'L/m^2' tag is a quoting convention, not a normalisation.
    lvc <- log(1.68);  label("Volume of the serum compartment shared by all 14 species, for a child with BSA 1 m^2 (V, L)")                      # ESM Table S13 German/Czech (V 1.68, RSE 1.1%; bootstrap 1.68, 95% CI 1.62-1.72)
    lcl <- log(0.126); label("Initial clearance of fully PEGylated asparaginase, for a child with BSA 1 m^2 (CLinitial, L/day)")                 # ESM Table S13 German/Czech (CLinitial 0.126, RSE 1.1%; bootstrap 0.126, 95% CI 0.123-0.130)
    lq  <- log(0.926); label("Intercompartmental clearance driving the de-PEGylation transit chain, for a child with BSA 1 m^2 (Qtr, L/day)")    # ESM Table S13 German/Czech (Qtr 0.926, RSE 1.9%; bootstrap 0.925, 95% CI 0.883-0.968)

    # ---- Body surface area -------------------------------------------------
    # Linear and centred on the model-building median 0.79 m^2, with a single
    # shared slope on CLinitial and Qtr (Table S13 rows 'FBSA on V' and
    # 'FBSA on CLinitial + Qtr'). Shared-exponent naming per
    # parameter-names.md: e_<cov>_<param1>_<param2>.
    e_bsa_vc   <- 1.57; label("Linear slope of BSA on V, centred on 0.79 m^2 (1/m^2)")                     # ESM Table S13 German/Czech (FBSA on V 1.57, RSE 1.6%; bootstrap 1.57, 95% CI 1.52-1.62)
    e_bsa_cl_q <- 1.45; label("Linear slope of BSA on CLinitial and Qtr, centred on 0.79 m^2 (1/m^2)")     # ESM Table S13 German/Czech (FBSA on CLinitial + Qtr 1.45, RSE 1.6%; bootstrap 1.45, 95% CI 1.40-1.50)

    # ---- Treatment-phase fractional changes --------------------------------
    # Multiplicative (1 + F) factors relative to the FIRST induction dose
    # (protocol IA day 12), which is the reference and carries no parameter --
    # ESM Table S13 footnote (b), 'fractional change in V / CLinitial compared
    # to the 1st administration in induction'.
    e_piad26_vc <- -0.15;  label("Fractional change in V at the second induction dose, protocol IA day 26 (unitless)")   # ESM Table S13 German/Czech (F(V PIAd26) -0.15, RSE 5.6%; bootstrap -0.150, 95% CI -0.169 to -0.134)
    e_piid8_vc  <- -0.272; label("Fractional change in V at the re-induction dose, protocol II day 8 (unitless)")        # ESM Table S13 German/Czech (F(V PIId8) -0.272, RSE 3.1%; bootstrap -0.274, 95% CI -0.294 to -0.255)
    e_piad26_cl <- -0.111; label("Fractional change in CLinitial at the second induction dose, protocol IA day 26 (unitless)")  # ESM Table S13 German/Czech (F(CLinitial PIAd26) -0.111, RSE 9.6%; bootstrap -0.111, 95% CI -0.131 to -0.089)
    e_piid8_cl  <- -0.411; label("Fractional change in CLinitial at the re-induction dose, protocol II day 8 (unitless)")       # ESM Table S13 German/Czech (F(CLinitial PIId8) -0.411, RSE 2.4%; bootstrap -0.411, 95% CI -0.432 to -0.390)

    # ---- Age and sex on initial clearance ----------------------------------
    e_age_cl  <- 0.017;  label("Linear slope of age above 8 years on CLinitial (1/year)")  # ESM Table S13 German/Czech (F(Age > 8 years) 0.017, RSE 20.1%; bootstrap 0.017, 95% CI 0.011-0.024)
    e_sexf_cl <- -0.073; label("Fractional change in CLinitial for females (unitless)")    # ESM Table S13 German/Czech (F(Sex) -0.073, RSE 16.5%; bootstrap -0.073, 95% CI -0.097 to -0.048)

    # ---- Inter-individual variability --------------------------------------
    # Table S13 reports IIV on CLinitial as 24.2 percent. The variability is
    # exponential (Wurthwein 2021 control stream: CLP = TVCLP*EXP(ETA(2)+IOVCL)),
    # and the percentage is the log-normal CV, so omega^2 = log(1 + CV^2).
    # Cross-check: the control stream's $OMEGA for IIV CLP is 0.0567, and
    # sqrt(exp(0.0567) - 1) = 24.15 percent, matching the printed 24.2.
    # There is NO IIV on V or Qtr -- the control stream fixes both to 0 and the
    # 2021 ESM Section 2 states 'inter-individual variability (IIV) on Qtr or V
    # could be neglected'.
    etalcl ~ 0.056918  # ESM Table S13 German/Czech (IIV CLinitial 24.2%, RSE 3.4%; bootstrap 24.2%, 95% CI 22.4-25.8) converted as log(1 + 0.242^2)

    # ---- Inter-occasion variability ----------------------------------------
    # On BOTH V and CLinitial, one occasion per administration. nlmixr2 has no
    # equivalent of NONMEM's `$OMEGA BLOCK(1) SAME`, so occasions 2 and 3 get
    # their own etas with the variance fixed equal to the occasion-1 estimate
    # (Kloos_2021_pegasparaginase.R / Chen_2023_nemonoxacin.R precedent).
    etaiov_vc_1 ~ 0.020524        # ESM Table S13 German/Czech (IOV V 14.4%, RSE 7.8%; bootstrap 14.2%, 95% CI 11.2-16.8) converted as log(1 + 0.144^2)
    etaiov_vc_2 ~ fixed(0.020524) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_3 ~ fixed(0.020524) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_1 ~ 0.049805        # ESM Table S13 German/Czech (IOV CLinitial 22.6%, RSE 3.5%; bootstrap 22.6%, 95% CI 20.8-24.3) converted as log(1 + 0.226^2)
    etaiov_cl_2 ~ fixed(0.049805) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fixed(0.049805) # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---- Residual error ----------------------------------------------------
    # Wurthwein 2021 ESM Section 4 $ERROR block:
    #   W = SQRT(THETA(7)**2 * IPRED**2 + THETA(8)**2); Y = IPRED + W*EPS(1)
    # with $SIGMA 1 FIX, i.e. a combined proportional + additive standard
    # deviation, which is exactly nlmixr2's `prop() + add()`.
    propSd <- 0.187; label("Proportional residual error (fraction)")  # ESM Table S13 German/Czech (Proportional error 18.7%, RSE 2.6%; bootstrap 18.8%, 95% CI 17.5-20.0)
    addSd  <- 9.03;  label("Additive residual error (U/L)")           # ESM Table S13 German/Czech (Additive error 9.03 U/L, RSE 35.2%; bootstrap 9.19, 95% CI 3.93-18.1)
  })

  model({
    # ---- Administration occasion ------------------------------------------
    # OCC = 1 protocol IA day 12 (reference), 2 protocol IA day 26,
    # 3 protocol II day 8. Records outside 1-3 fall back to the reference
    # phase and contribute no IOV.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)

    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 + occ3 * etaiov_vc_3
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3

    # ---- Body surface area -------------------------------------------------
    # TVV1 = THETA(1) * (1 + SCV1 * (BSA - 0.79)) in the authors' control
    # stream. Because ini() holds the value printed for BSA = 1 m^2 rather
    # than THETA(1) itself, the factor is divided by its own value at
    # BSA = 1 m^2 (i.e. by 1 + slope * 0.21). At BSA = 1 both factors are
    # exactly 1, so the ini() values are the table values unaltered.
    fbsa_vc <- (1 + e_bsa_vc   * (BSA - 0.79)) / (1 + e_bsa_vc   * 0.21)
    fbsa_cl <- (1 + e_bsa_cl_q * (BSA - 0.79)) / (1 + e_bsa_cl_q * 0.21)

    # ---- Age, sex and treatment phase on initial clearance -----------------
    # Age is a hockey stick with the break point fixed at 8 years: the
    # sub-8-year branch has its slope fixed to 0 in the control stream, so the
    # factor is 1 there.
    fage_cl <- 1 + e_age_cl * (AGE - 8) * (AGE > 8)
    fsex_cl <- 1 + e_sexf_cl * SEXF

    fphase_vc <- 1 + occ2 * e_piad26_vc + occ3 * e_piid8_vc
    fphase_cl <- 1 + occ2 * e_piad26_cl + occ3 * e_piid8_cl

    # ---- Individual parameters ---------------------------------------------
    # V and CLinitial carry IOV; CLinitial additionally carries IIV. Qtr has
    # neither (both were fixed to zero by the authors) but does take the same
    # BSA factor as CLinitial.
    vc <- exp(lvc + iov_vc) * fbsa_vc * fphase_vc
    cl <- exp(lcl + etalcl + iov_cl) * fbsa_cl * fage_cl * fsex_cl * fphase_cl
    q  <- exp(lq) * fbsa_cl

    # ---- De-PEGylation transit chain ---------------------------------------
    # Wurthwein 2021 Figure 1 and the ESM Section 4 $DES block. Fourteen
    # compartments in series, all sharing the single serum volume vc. Drug
    # steps down the chain at Qtr/V, standing in for progressive hydrolysis of
    # the PEG moiety, and every compartment is eliminated at CLinitial/V. The
    # terminal compartment additionally loses drug at Qtr/V: the authors'
    # original model gave it a separate induced clearance CLinduced, and the
    # published simplification sets CLinduced = Qtr (2021 ESM Section 2,
    # dOFV = 0.6), which is why its total loss rate matches every other
    # compartment's. That is what makes apparent clearance climb from
    # CLinitial toward CLinitial + Qtr over a dosing interval as mass migrates
    # down the chain.
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

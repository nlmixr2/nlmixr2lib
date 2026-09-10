Wurthwein_2025_pegasparaginase_italian <- function() {
  description <- "Fourteen-compartment de-PEGylation transit population PK model for intravenous PEGylated asparaginase (PEG-ASNase) in children with acute lymphoblastic leukemia treated in induction (protocol IA days 12 and 26) and re-induction (protocol II day 8) in the Italian group of the AIEOP-BFM ALL 2009 trial (Wurthwein 2025, ESM Table S13 model 223005). Structurally identical to the German/Czech model, refit independently on Italian data measured with the medac asparaginase activity test (MAAT) after conversion to the aspartic acid beta-hydroxamate (AHA) scale with the published two-band factors 1.23 for levels at or below 600 U/L and 1.42 above. Asparaginase activity is carried by a chain of 14 serial compartments sharing one serum volume; drug moves down the chain with intercompartmental clearance Qtr, which mimics stepwise de-PEGylation, and every compartment is eliminated with the initial clearance CLinitial while the terminal compartment is eliminated with CLinitial + Qtr. Body surface area enters volume and the two clearance terms linearly, centred on 0.79 m^2. Initial clearance rises linearly with age above 8 years and is lower in females. Initial clearance and volume drop stepwise between administrations, encoded as fractional changes against the first induction dose. Inter-individual variability on initial clearance, inter-occasion variability on both initial clearance and volume, and combined proportional plus additive residual error."
  reference   <- "Wurthwein G, Siebel C, Lanvers-Kaminsky C, Smisek P, Nath CE, Matteo C, Rizzari C, Schrappe M, Boos J. PEGylated Asparaginase in Children with Acute Lymphoblastic Leukemia Treated within the AIEOP-BFM ALL 2009 Trial: Population Pharmacokinetics and Drug Exposure. Eur J Drug Metab Pharmacokinet. 2025;50(6):683-696. doi:10.1007/s13318-025-00962-3"
  vignette    <- "Wurthwein_2025_pegasparaginase"
  units       <- list(time = "day", dosing = "U", concentration = "U/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed by the Mosteller formula (Wurthwein 2025 Methods, Pharmacokinetic Model, citing reference 14). Italian median 0.77 m^2, range 0.28-2.30 (Wurthwein 2025 Table 1). BSA enters volume and the two clearance terms LINEARLY and CENTRED, not allometrically -- ESM Table S13 footnote (a), 'linear increase in V / CLinitial / Qtr with BSA (centred on the median)'. The centring constant is taken as 0.79 m^2, the value hard-coded in the Wurthwein 2021 ESM Section 4 control stream (TVV1 = THETA(1)*(1+SCV1*(BSA-0.79))) that this model is an adaptation of; the Italian cohort median is 0.77 m^2 and the paper does not state whether the constant was re-centred for this fit. The choice is not load-bearing: re-centring from 0.79 to 0.77 moves typical-value predictions by less than 1 percent across the whole observed BSA range, because the factor is renormalised to 1 at BSA = 1 m^2 in model(). Flagged in the vignette's Assumptions and deviations section. IMPORTANT: the control stream contains no multiplication by BSA, so V, CLinitial and Qtr are ABSOLUTE quantities in L and L/day; the 'L/m^2' unit tags in Table S13 record the authors' convention of quoting the absolute value for a child with BSA = 1 m^2, not a per-square-metre normalisation.",
      source_name        = "BSA"
    ),
    AGE = list(
      description        = "Age at the time of the PEG-asparaginase administration",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Italian median 5.1 years, range 1.04-18.0 (Wurthwein 2025 Table 1). Enters initial clearance as a hockey stick with the break point FIXED at 8 years: clearance is flat up to 8 years and rises linearly above it (Wurthwein 2021 ESM Section 4 control stream, two AGEGTP branches with the sub-8-year slope fixed to 0; ESM Section 2 explains the fixing, 'only a minor increase in OFV for break points between 7 and 11 years'). Age at administration, not at enrolment, so time-varying across the protocol. Wurthwein 2025 Section 3.3.2 confirms the covariate carried over to the Italian fit: 'the covariates BSA, age, and sex were confirmed; parameter estimates were comparable in both trial groups.'",
      source_name        = "AGEGTP"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Italian cohort 927 male / 676 female (Wurthwein 2025 Table 1). Females have lower initial clearance -- ESM Table S13 footnote (d), 'fractional change in CLinitial for females'. The source column is a 1/2 coded SEX (Wurthwein 2021 ESM Section 4 control stream: IF(SEX.EQ.1) CLSEX = 1 for males, IF(SEX.EQ.2) CLSEX = (1 + THETA(15)) for females), so SEXF = SEX - 1 preserves both the coefficient sign and the male reference category and the printed estimate is used unchanged.",
      source_name        = "SEX (1 = male, 2 = female)"
    ),
    OCC = list(
      description        = "Integer administration-occasion index; drives both the treatment-phase effects on initial clearance and volume and the inter-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (first induction dose, protocol IA day 12) is the model reference for every fractional change",
      notes              = "One occasion per PEG-asparaginase administration. The Wurthwein 2021 ESM Section 4 control stream carries a single OCC column that simultaneously selects the treatment-phase covariate factors (V1KOV, CLKOV) and the IOV etas (IOVV1, IOVCL), with the header comment ';; OCC: 2=PIA d12, 3=PIA d26, 5=PII d8'. Renumbered here to the 1..N form the register prescribes: OCC = 1 protocol IA day 12 (induction dose 1, the reference), OCC = 2 protocol IA day 26 (induction dose 2), OCC = 3 protocol II day 8 (re-induction). Records with OCC outside 1-3 receive the reference phase factors and no IOV contribution.",
      source_name        = "OCC (2 = PIA d12, 3 = PIA d26, 5 = PII d8)"
    )
  )

  covariatesDataExcluded <- list(
    ABPEG_IGM = list(
      description        = "Pre-existing anti-polyethylene-glycol IgM antibody level in the sample drawn before the first induction dose",
      units              = "(assay-specific mean fluorescence intensity, linear scale)",
      type               = "continuous",
      notes              = "Not estimable in the Italian cohort and therefore absent from this model. Wurthwein 2025 Section 3.3.2: 'The subgroup of 356 patients with analyzed anti-PEG IgM antibody levels prior to the first dose in induction (324 levels below and only 32 levels above the cut-point ...) was too small to evaluate the impact of this covariate in the Italian group.' Section 3.4.2 records that omitting it is adequate because 'anti-PEG antibody levels explained only a small part of variability in CLinitial'. Retained in the two extended German/Czech models of this paper (Wurthwein_2025_pegasparaginase_hrPostinduction.R and Wurthwein_2025_pegasparaginase_r2ea.R)."
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      notes              = "Screened as an alternative body-size descriptor during development of the parent German/Czech model and rejected in favour of BSA; carried over unchanged. Wurthwein 2021 ESM Section 2: 'Allometric weight scaling with fixed exponents (0.75 on clearance-terms, 1 on V), estimated exponents or linear scaling with body weight did not improve the model compared to BSA scaling.'"
    ),
    COUNTRY = list(
      description        = "Trial-group / country indicator (German/Czech versus Italian)",
      units              = "(binary)",
      type               = "binary",
      notes              = "Deliberately NOT modelled. The authors tried and abandoned a pooled model carrying country as a covariate -- Wurthwein 2025 Section 3.3.2: 'PopPK models for the combined data set that included the covariate <country> became unstable: some runs showed extremely high relative standard error (RSE) values in one or the other parameter estimate. ... we decided not to combine the two data sets.' That decision is precisely why this model and Wurthwein_2025_pegasparaginase_germanCzech.R are separate independent fits rather than one model with a country term."
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
    n_subjects          = 1603L,
    n_studies           = 1L,
    n_observations      = 6894L,
    age_range           = "1.04-18.0 years; median 5.1",
    bsa_range           = "0.28-2.30 m^2; median 0.77",
    sex_female_pct      = 42.2,
    race_ethnicity      = "Not reported (Italian trial sites of the AIEOP-BFM ALL 2009 trial).",
    disease_state       = "Newly diagnosed pediatric acute lymphoblastic leukemia enrolled in the AIEOP-BFM ALL 2009 trial (EudraCT 2007-004270-43, NCT01117441) between 1 June 2010 and 28 February 2017.",
    dose_range          = "2500 U/m^2 per dose as a 2-hour intravenous infusion, capped at an absolute maximum of 3750 U/dose. This model covers the two induction doses on protocol IA days 12 and 26 and the single re-induction dose on protocol II day 8. Administered dose per m^2, median [range]: induction 1900 [705, 4500]; re-induction 1930 [590, 9000] (Table 1).",
    regions             = "Italy.",
    bioanalytic_methods = "Asparaginase activity in serum by the medac asparaginase activity test (MAAT); lower limit of quantification 30 U/L, calibration range 30-600 U/L. The MAAT calibrates against native E. coli asparaginase and therefore OVERESTIMATES activity relative to the AHA test used for the German/Czech and Australian groups. All activities were divided to the AHA scale before modelling using the two-band conversion factors of Lanvers-Kaminsky et al.: 1.23 for MAAT results at or below 600 U/L and 1.42 above 600 U/L. Wurthwein 2025 Section 3.3.1 selected these over the single factors 1.24 (MAAT manual) and 1.37 by external evaluation of the published German/Czech model, and ESM Table S12 shows the two-band pair minimising bias and imprecision on individual predictions (-3.7 percent bias, 7.1 percent precision) and on CLinitial (5.4 percent bias) versus 20.6 and 8.4 percent bias for 1.24 and 1.37. Samples below the limit of quantification (1.6 percent) were excluded from the fit.",
    notes               = "Demographics from Wurthwein 2025 Table 1 (Italian column); the female percentage is 676/1603. Sampling was scheduled before each treatment phase and 7 and 14 days after each dose (ESM Table S2). Records indicating silent inactivation and every record at or after a hypersensitivity reaction were EXCLUDED before fitting, so this model describes standard elimination only. Fitted independently of the German/Czech data: the authors explicitly rejected a pooled model with a country covariate as unstable and over-parameterised (Section 3.3.2). NONMEM 7.4.4, FOCE with interaction; 1000-replicate bootstrap, 100 percent successful; IIV shrinkage 26.0 percent, condition number 62.45 (ESM Table S13 footnote (f)). Both the residual error and the inter-individual variability are smaller than in the German/Czech fit, which the authors offer as the likely reason for this cohort's slightly better predictivity when monitoring samples are missing (Section 3.3.3, ESM Table S14)."
  )

  ini({
    # ---- Structural PK -----------------------------------------------------
    # ESM Table S13, 'Italian data' Estimate column (model 223005). Values are
    # quoted for a child with BSA = 1 m^2; the BSA factors in model() are
    # renormalised to 1 at BSA = 1 m^2 so these are used exactly as printed.
    lvc <- log(1.53);  label("Volume of the serum compartment shared by all 14 species, for a child with BSA 1 m^2 (L)")                    # ESM Table S13 Italian (V 1.53, RSE 1.0%; bootstrap 1.53, 95% CI 1.49-1.57)
    lcl <- log(0.120); label("Initial clearance of fully PEGylated asparaginase, for a child with BSA 1 m^2 (L/day)")               # ESM Table S13 Italian (CLinitial 0.120, RSE 1.1%; bootstrap 0.120, 95% CI 0.116-0.124)
    lq  <- log(0.787); label("Intercompartmental clearance driving the de-PEGylation transit chain, for a child with BSA 1 m^2 (L/day)")  # ESM Table S13 Italian (Qtr 0.787, RSE 2.3%; bootstrap 0.787, 95% CI 0.744-0.828)

    # ---- Body surface area -------------------------------------------------
    e_bsa_vc   <- 1.48; label("Linear slope of BSA on V, centred on 0.79 m^2 (1/m^2)")                  # ESM Table S13 Italian (FBSA on V 1.48, RSE 1.7%; bootstrap 1.48, 95% CI 1.43-1.53)
    e_bsa_cl_q <- 1.39; label("Linear slope of BSA on CLinitial and Qtr, centred on 0.79 m^2 (1/m^2)")  # ESM Table S13 Italian (FBSA on CLinitial + Qtr 1.39, RSE 2.0%; bootstrap 1.39, 95% CI 1.34-1.45)

    # ---- Treatment-phase fractional changes --------------------------------
    # Relative to the first induction dose (protocol IA day 12), which is the
    # reference and carries no parameter -- ESM Table S13 footnote (b).
    e_piad26_vc <- -0.155; label("Fractional change in V at the second induction dose, protocol IA day 26 (unitless)")          # ESM Table S13 Italian (F(V PIAd26) -0.155, RSE 6.0%; bootstrap -0.155, 95% CI -0.174 to -0.137)
    e_piid8_vc  <- -0.256; label("Fractional change in V at the re-induction dose, protocol II day 8 (unitless)")               # ESM Table S13 Italian (F(V PIId8) -0.256, RSE 3.8%; bootstrap -0.256, 95% CI -0.277 to -0.239)
    e_piad26_cl <- -0.095; label("Fractional change in CLinitial at the second induction dose, protocol IA day 26 (unitless)")  # ESM Table S13 Italian (F(CLinitial PIAd26) -0.095, RSE 12.7%; bootstrap -0.094, 95% CI -0.117 to -0.071)
    e_piid8_cl  <- -0.293; label("Fractional change in CLinitial at the re-induction dose, protocol II day 8 (unitless)")       # ESM Table S13 Italian (F(CLinitial PIId8) -0.293, RSE 3.4%; bootstrap -0.292, 95% CI -0.313 to -0.272)

    # ---- Age and sex on initial clearance ----------------------------------
    e_age_cl  <- 0.017;  label("Linear slope of age above 8 years on CLinitial (1/year)")  # ESM Table S13 Italian (F(Age > 8 years) 0.017, RSE 23.9%; bootstrap 0.017, 95% CI 0.010-0.025)
    e_sexf_cl <- -0.070; label("Fractional change in CLinitial for females (unitless)")    # ESM Table S13 Italian (F(Sex) -0.070, RSE 18.2%; bootstrap -0.070, 95% CI -0.095 to -0.044)

    # ---- Inter-individual variability --------------------------------------
    # Exponential IIV; the printed percentage is the log-normal CV, so
    # omega^2 = log(1 + CV^2). No IIV on V or Qtr (Wurthwein 2021 ESM
    # Section 2: 'inter-individual variability (IIV) on Qtr or V could be
    # neglected').
    etalcl ~ 0.039631  # ESM Table S13 Italian (IIV CLinitial 20.1%, RSE 4.0%; bootstrap 20.0%, 95% CI 18.5-21.5) converted as log(1 + 0.201^2)

    # ---- Inter-occasion variability ----------------------------------------
    # On BOTH V and CLinitial, one occasion per administration. Occasions 2
    # and 3 carry their own etas with the variance fixed equal to the
    # occasion-1 estimate, encoding NONMEM's `$OMEGA BLOCK(1) SAME`.
    etaiov_vc_1 ~ 0.018341        # ESM Table S13 Italian (IOV V 13.6%, RSE 7.5%; bootstrap 13.5%, 95% CI 11.4-15.3) converted as log(1 + 0.136^2)
    etaiov_vc_2 ~ fixed(0.018341) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_3 ~ fixed(0.018341) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_1 ~ 0.045193        # ESM Table S13 Italian (IOV CLinitial 21.5%, RSE 3.7%; bootstrap 21.4%, 95% CI 19.7-23.1) converted as log(1 + 0.215^2)
    etaiov_cl_2 ~ fixed(0.045193) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fixed(0.045193) # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---- Residual error ----------------------------------------------------
    # Combined proportional + additive SD, per the Wurthwein 2021 ESM
    # Section 4 $ERROR block. On the AHA-converted scale: the Italian
    # additive term is four times the German/Czech one, consistent with the
    # MAAT's higher lower limit of quantification (30 versus 5 U/L).
    propSd <- 0.134; label("Proportional residual error (fraction)")  # ESM Table S13 Italian (Proportional error 13.4%, RSE 3.7%; bootstrap 13.4%, 95% CI 12.2-14.5)
    addSd  <- 37.1;  label("Additive residual error (U/L)")           # ESM Table S13 Italian (Additive error 37.1 U/L, RSE 15.2%; bootstrap 36.9, 95% CI 25.7-51.3)
  })

  model({
    # ---- Administration occasion ------------------------------------------
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)

    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 + occ3 * etaiov_vc_3
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3

    # ---- Body surface area -------------------------------------------------
    # Linear and centred on 0.79 m^2, renormalised to 1 at BSA = 1 m^2 so the
    # ini() values are the printed table entries unaltered.
    fbsa_vc <- (1 + e_bsa_vc   * (BSA - 0.79)) / (1 + e_bsa_vc   * 0.21)
    fbsa_cl <- (1 + e_bsa_cl_q * (BSA - 0.79)) / (1 + e_bsa_cl_q * 0.21)

    # ---- Age, sex and treatment phase on initial clearance -----------------
    fage_cl <- 1 + e_age_cl * (AGE - 8) * (AGE > 8)
    fsex_cl <- 1 + e_sexf_cl * SEXF

    fphase_vc <- 1 + occ2 * e_piad26_vc + occ3 * e_piid8_vc
    fphase_cl <- 1 + occ2 * e_piad26_cl + occ3 * e_piid8_cl

    # ---- Individual parameters ---------------------------------------------
    vc <- exp(lvc + iov_vc) * fbsa_vc * fphase_vc
    cl <- exp(lcl + etalcl + iov_cl) * fbsa_cl * fage_cl * fsex_cl * fphase_cl
    q  <- exp(lq) * fbsa_cl

    # ---- De-PEGylation transit chain ---------------------------------------
    # Wurthwein 2021 Figure 1 and ESM Section 4 $DES. Fourteen compartments in
    # series sharing the single serum volume vc; drug steps down the chain at
    # Qtr/V (progressive hydrolysis of the PEG moiety) and every compartment
    # is eliminated at CLinitial/V. The terminal compartment additionally
    # loses drug at Qtr/V because the published simplification sets the
    # induced clearance CLinduced equal to Qtr.
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
    # Total catalytic activity: every species in the chain contributes.
    # Predictions are on the AHA scale, so a MAAT measurement must be divided
    # by 1.23 (at or below 600 U/L) or 1.42 (above 600 U/L) before comparison.
    Cc <- (central + transit1 + transit2 + transit3 + transit4 + transit5 +
             transit6 + transit7 + transit8 + transit9 + transit10 +
             transit11 + transit12 + transit13) / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}

Wurthwein_2021_pegasparaginase <- function() {
  description <- "Fourteen-compartment de-PEGylation transit population PK model for intravenous PEGylated asparaginase (PEG-ASNase) in children with acute lymphoblastic leukemia treated in induction (protocol IA days 12 and 26) and re-induction (protocol II day 8) in the German and Czech part of the AIEOP-BFM ALL 2009 trial (Wurthwein 2021 Final Pharmacokinetic Model). Asparaginase activity is carried by a chain of 14 serial compartments sharing one serum volume; drug moves down the chain with intercompartmental clearance Qtr, which mimics stepwise de-PEGylation, and every compartment is eliminated with the initial clearance CLinitial while the terminal compartment is eliminated with CLinitial + Qtr. Apparent clearance therefore rises within a dosing interval from CLinitial toward CLinitial + Qtr. Body surface area enters volume and the two clearance terms linearly, centred on the model-building median of 0.79 m^2. Initial clearance rises linearly with age above 8 years and is lower in females. Initial clearance and volume drop stepwise between administrations, encoded as fractional changes against the first induction dose. Inter-individual variability on initial clearance, inter-occasion variability on both initial clearance and volume, and combined proportional plus additive residual error."
  reference <- "Wurthwein G, Lanvers-Kaminsky C, Siebel C, Gerss J, Moricke A, Zimmermann M, Stary J, Smisek P, Schrappe M, Rizzari C, Zucchetti M, Hempel G, Wicha SG, Boos J, on behalf of the AIEOP-BFM ALL 2009 Asparaginase Working Party. Population Pharmacokinetics of PEGylated Asparaginase in Children with Acute Lymphoblastic Leukemia: Treatment Phase Dependency and Predictivity in Case of Missing Data. Eur J Drug Metab Pharmacokinet. 2021;46(2):289-300. doi:10.1007/s13318-021-00670-8"
  vignette <- "Wurthwein_2021_pegasparaginase"
  units <- list(time = "day", dosing = "U", concentration = "U/L")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Computed by the Mosteller formula (Wurthwein 2021 Section 2.4) from weight and height recorded before each dose (Section 2.1). Model-building median 0.78 m^2, range 0.41-2.42 (Table 1). BSA enters volume and the two clearance terms LINEARLY and CENTRED on 0.79 m^2, not allometrically: the ESM Section 4 control stream writes TVV1 = THETA(1)*(1+SCV1*(BSA-0.79)), TVCLP = THETA(2)*(1+SCCL*(BSA-0.79)) and TVQ = THETA(4)*(1+SCCL*(BSA-0.79)); Table 3 footnote (b) prints the same form, F_i = F*(1 + F_BSA*(BSA_i - 0.79)). The control stream multiplies nothing by BSA, so THETA(1), THETA(2) and THETA(4) are absolute values (L, L/day) for a child at the 0.79 m^2 centring point. Table 3 quotes them per m^2: 'During NONMEM estimation, values for CLinitial, CLinduced, Qtr and V are reported for BSA=0.79 m2; values were converted from L/day/0.79m2 or L/0.79m2 to L/day/m2 or L/m2 for better comparison'. The model() code therefore multiplies the printed per-m^2 value by 0.79 m^2 to recover the NONMEM THETA.",
      source_name = "BSA"
    ),
    AGE = list(
      description = "Age at the time of the PEG-asparaginase administration",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Model-building median 5.13 years, range 1.03-17.9 (Table 1). Enters initial clearance as a hockey stick with the break point FIXED at 8 years: flat up to 8 years and rising linearly above it (Table 3 footnote (d), F_i = F*(1 + F_age>8years*(age_i - 8))). The ESM Section 4 control stream encodes the two branches as IF(AGEGTP.LE.8) CLAGEGTP = (1 + THETA(13)*(AGEGTP - 8)) with THETA(13) fixed to 0, and IF(AGEGTP.GT.8) CLAGEGTP = (1 + THETA(14)*(AGEGTP - 8)). ESM Section 2 gives the reason the break point is fixed: 'Further analyses indicated only a minor increase in OFV for break points between 7 and 11 years; thus, the break point was fixed to 8 years.' The control stream comment defines AGEGTP as 'Age at time of PEG-ASNase administration', so the covariate is time-varying across the protocol.",
      source_name = "AGEGTP"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Model-building dataset 779 male / 595 female (Table 1). Females have 7.1 percent lower initial clearance (Results Section 3.2.1; Table 3 footnote (e), 'fractional change in CLinitial for females'). The source column is a 1/2-coded SEX: the ESM Section 4 control stream has IF(SEX.EQ.1) CLSEX = 1 (male, 'Most common', the reference) and IF(SEX.EQ.2) CLSEX = (1 + THETA(15)). SEXF = SEX - 1 keeps both the sign of the estimate and the male reference, so the printed value is used unchanged.",
      source_name = "SEX (1 = male, 2 = female)"
    ),
    OCC = list(
      description = "Integer administration-occasion index; drives both the treatment-phase effects on initial clearance and volume and the inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = "1 (first induction dose, protocol IA day 12) is the model reference for every fractional change",
      notes = "One occasion per PEG-asparaginase administration, exactly as the authors encode it. The ESM Section 4 control stream carries a single OCC column that selects both the treatment-phase factors (V1KOV, CLKOV) and the inter-occasion etas (IOVV1, IOVCL), with the header comment ';; OCC: 2=PIA d12, 3=PIA d26, 5=PII d8'. This model renumbers those occasions to the 1..N form of the covariate register: OCC = 1 is protocol IA day 12 (induction dose 1, the reference), OCC = 2 is protocol IA day 26 (induction dose 2) and OCC = 3 is protocol II day 8 (re-induction, non-high-risk patients only). Records with OCC outside 1-3 receive the reference phase factors and no IOV contribution.",
      source_name = "OCC (2 = PIA d12, 3 = PIA d26, 5 = PII d8)"
    )
  )

  # Covariates screened by Wurthwein 2021 but NOT retained in the final model.
  # Documented for provenance only; deliberately absent from model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as an alternative body-size descriptor and rejected in favour of BSA. ESM Section 2: 'Allometric weight scaling with fixed exponents (0.75 on clearance-terms, 1 on V), estimated exponents or linear scaling with body weight did not improve the model compared to BSA scaling'; Section 2.5: 'For the highly correlated covariates, i.e., body weight and BSA (r = 0.988), the covariate with the best improvement in objective function value (OFV) was retained in the model.' Model-building median 19.55 kg, range 8.5-128 (Table 1)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "PEG-asparaginase, fully PEGylated", units = "U", specimen = "serum", verified = TRUE),
    transit1 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 1",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 2",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 3",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 4",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 5",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit6 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 6",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit7 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 7",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit8 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 8",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit9 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 9",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit10 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 10",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit11 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 11",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit12 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 12",
      units = "U",
      specimen = "serum",
      verified = TRUE
    ),
    transit13 = list(
      analyte = "PEG-asparaginase, de-PEGylation step 13",
      units = "U",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2545L,
    n_studies = 1L,
    n_observations = 11486L,
    age_range = "1.03-18.0 years; median 5.13 (model-building) and 5.06 (testing)",
    weight_range = "7.7-128 kg; median 19.55 (model-building) and 19 (testing)",
    bsa_range = "0.39-2.44 m^2; median 0.78 (model-building) and 0.77 (testing)",
    sex_female_pct = 42.1,
    race_ethnicity = "Not reported (German and Czech trial sites of the AIEOP-BFM ALL 2009 trial).",
    disease_state = "Newly diagnosed pediatric acute lymphoblastic leukemia enrolled in the AIEOP-BFM ALL 2009 trial (EudraCT 2007-004270-43, NCT01117441).",
    dose_range = "2500 U/m^2 per dose as a 2-hour intravenous infusion, maximal absolute dose 3750 U. Two induction doses (protocol IA days 12 and 26) and, for non-high-risk patients, one re-induction dose (protocol II day 8, scheduled 18 weeks after the last induction dose). Administered dose median 1940 U (range 710-4850) in the model-building dataset (Table 1).",
    regions = "Germany and Czech Republic.",
    bioanalytic_methods = "Asparaginase activity in serum, centrally analysed (ESM Section 1); lower limit of quantification 5 U/L. Samples below the limit (1.9 percent of the Final Dataset) were omitted from the fit (Section 2.4).",
    notes = "The Final Pharmacokinetic Model was fitted to the Final Dataset of 2545 patients and 11,486 samples (Discussion Section 4.1), after 26 observations with |CWRES| > 4 were excluded for stability (ESM Section 2). Covariate selection was first done on the Model Building Dataset (1374 patients, 6069 samples; diagnosis on or before 31 December 2015) and externally validated on the Testing Dataset (1253 patients, 5523 samples), then repeated and confirmed on the Final Dataset. Demographics are from Table 1 (model-building and testing columns); the female percentage pools both, (595 + 512)/(1374 + 1253). Monitoring was scheduled 7 and 14 days after each dose, so 45 percent of samples fall in the day 7 +/- 1 and 41.9 percent in the day 14 +/- 1 window and only 0.9 percent in days 0-5 (Sections 3.1 and 3.2.3). Samples after a hypersensitivity reaction, samples indicating silent inactivation (below 100 U/L within 8 days and/or undetectable within 15 days) and pharmacologically implausible samples were EXCLUDED (Section 2.3), so the model describes standard elimination only. NONMEM 7.3.0 / 7.4.4, FOCE with interaction; 1000-replicate bootstrap, 92.9 percent successful (Table 3)."
  )

  ini({
    # ---- Structural PK -----------------------------------------------------
    # Table 3, 'Final PK model' column. Values are printed per m^2: the NONMEM
    # THETAs are absolute values for a child at the 0.79 m^2 centring point and
    # were divided by 0.79 for the table (Table 3 footnote). model() multiplies
    # them back by 0.79. See covariateData$BSA.
    lvc <- log(1.69);  label("Serum volume shared by all 14 species, per m^2 at the 0.79 m^2 centring point (L/m^2)")                      # Table 3 Final PK model, V 1.69 L/m2 (RSE 1.1%; bootstrap 1.68, 95% CI 1.63-1.73)
    lcl <- log(0.126); label("Initial clearance of fully PEGylated asparaginase, per m^2 at the 0.79 m^2 centring point (L/day/m^2)")     # Table 3 Final PK model, CLinitial 0.126 L/day/m2 (RSE 1.0%; bootstrap 0.126, 95% CI 0.123-0.130)
    lq  <- log(0.918); label("Intercompartmental clearance driving the de-PEGylation chain, per m^2 at the 0.79 m^2 centring point (L/day/m^2)") # Table 3 Final PK model, Qtr 0.918 L/day/m2 (RSE 1.9%; bootstrap 0.918, 95% CI 0.870-0.960)

    # ---- Body surface area -------------------------------------------------
    # Linear and centred on 0.79 m^2; one slope on V and one shared by
    # CLinitial and Qtr (Table 3 footnote (b); ESM Section 4 SCV1 and SCCL).
    e_bsa_vc   <- 1.56; label("Linear slope of BSA on V, centred on 0.79 m^2 (1/m^2)")                  # Table 3 Final PK model, F_BSA on V 1.56 (RSE 1.5%; bootstrap 1.56, 95% CI 1.51-1.61)
    e_bsa_cl_q <- 1.44; label("Linear slope of BSA on CLinitial and Qtr, centred on 0.79 m^2 (1/m^2)")  # Table 3 Final PK model, F_BSA on CLinitial + Qtr 1.44 (RSE 1.6%; bootstrap 1.44, 95% CI 1.40-1.49)

    # ---- Treatment-phase fractional changes --------------------------------
    # (1 + F) factors against the FIRST induction dose (Table 3 footnote (c)).
    e_piad26_vc <- -0.159; label("Fractional change in V at the second induction dose, protocol IA day 26 (unitless)")            # Table 3 Final PK model, F(V Induction 2. admin) -0.159 (RSE 5.1%; bootstrap -0.159, 95% CI -0.174 to -0.142)
    e_piid8_vc  <- -0.284; label("Fractional change in V at the re-induction dose, protocol II day 8 (unitless)")                 # Table 3 Final PK model, F(V Reinduction) -0.284 (RSE 2.9%; bootstrap -0.283, 95% CI -0.300 to -0.266)
    e_piad26_cl <- -0.110; label("Fractional change in CLinitial at the second induction dose, protocol IA day 26 (unitless)")     # Table 3 Final PK model, F(CLinitial Induction 2. admin) -0.110 (RSE 9.5%; bootstrap -0.110, 95% CI -0.131 to -0.087)
    e_piid8_cl  <- -0.412; label("Fractional change in CLinitial at the re-induction dose, protocol II day 8 (unitless)")          # Table 3 Final PK model, F(CLinitial Reinduction) -0.412 (RSE 2.4%; bootstrap -0.411, 95% CI -0.430 to -0.391)

    # ---- Age and sex on initial clearance ----------------------------------
    e_age_cl  <- 0.018;  label("Linear slope of age above 8 years on CLinitial (1/year)")  # Table 3 Final PK model, F(Age > 8 years) 0.018 (RSE 19.3%; bootstrap 0.018, 95% CI 0.011-0.025)
    e_sexf_cl <- -0.071; label("Fractional change in CLinitial for females (unitless)")    # Table 3 Final PK model, F Sex -0.071 (RSE 16.7%; bootstrap -0.071, 95% CI -0.095 to -0.046)

    # ---- Inter-individual variability --------------------------------------
    # Exponential (ESM Section 4: CLP = TVCLP*EXP(ETA(2)+IOVCL)); the printed
    # percentage is taken as the log-normal CV, omega^2 = log(1 + CV^2). The
    # control stream fixes IIV on V and Qtr to 0 ('$OMEGA 0 FIX ; IIV V1',
    # '0 FIX ; IIV Q'), so there is no eta on either.
    etalcl ~ 0.056002  # Table 3 Final PK model, IIV CLinitial 24.0% (RSE 3.4%) converted as log(1 + 0.240^2)

    # ---- Inter-occasion variability ----------------------------------------
    # On V and CLinitial, one occasion per administration. The control stream
    # uses '$OMEGA BLOCK(1) SAME' for occasions 2 and 3; nlmixr2 has no SAME,
    # so those etas are fixed to the occasion-1 variance.
    etaiov_vc_1 ~ 0.017534        # Table 3 Final PK model, IOV V 13.3% (RSE 6.7%) converted as log(1 + 0.133^2)
    etaiov_vc_2 ~ fixed(0.017534) # SAME-equivalent: equal to the occasion-1 IOV variance on V
    etaiov_vc_3 ~ fixed(0.017534) # SAME-equivalent: equal to the occasion-1 IOV variance on V
    etaiov_cl_1 ~ 0.048958        # Table 3 Final PK model, IOV CLinitial 22.4% (RSE 3.4%) converted as log(1 + 0.224^2)
    etaiov_cl_2 ~ fixed(0.048958) # SAME-equivalent: equal to the occasion-1 IOV variance on CLinitial
    etaiov_cl_3 ~ fixed(0.048958) # SAME-equivalent: equal to the occasion-1 IOV variance on CLinitial

    # ---- Residual error ----------------------------------------------------
    # ESM Section 4 $ERROR: W = SQRT(THETA(7)**2*IPRED**2 + THETA(8)**2),
    # Y = IPRED + W*EPS(1), $SIGMA 1 FIX -- combined proportional + additive SD.
    propSd <- 0.189; label("Proportional residual error (fraction)")  # Table 3 Final PK model, proportional error 18.9% (RSE 2.2%; bootstrap 18.9, 95% CI 18.0-19.7)
    addSd  <- 8.74;  label("Additive residual error (U/L)")           # Table 3 Final PK model, additive error 8.74 U/L (RSE 32.2%; bootstrap 8.64, 95% CI 3.96-14.6)
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
    # ESM Section 4: TVV1 = THETA(1)*(1 + SCV1*(BSA - 0.79)), and the same
    # form with SCCL for CLP and Q. THETA = printed per-m^2 value * 0.79 m^2
    # (Table 3 footnote: values 'converted from L/day/0.79m2 or L/0.79m2 to
    # L/day/m2 or L/m2').
    fbsa_vc <- 0.79 * (1 + e_bsa_vc * (BSA - 0.79))
    fbsa_cl <- 0.79 * (1 + e_bsa_cl_q * (BSA - 0.79))

    # ---- Age, sex and treatment phase on initial clearance -----------------
    # Hockey stick at 8 years; the sub-8-year slope is fixed to 0 in the
    # control stream, so the factor is 1 there.
    fage_cl <- 1 + e_age_cl * (AGE - 8) * (AGE > 8)
    fsex_cl <- 1 + e_sexf_cl * SEXF

    fphase_vc <- 1 + occ2 * e_piad26_vc + occ3 * e_piid8_vc
    fphase_cl <- 1 + occ2 * e_piad26_cl + occ3 * e_piid8_cl

    # ---- Individual parameters ---------------------------------------------
    vc <- exp(lvc + iov_vc) * fbsa_vc * fphase_vc
    cl <- exp(lcl + etalcl + iov_cl) * fbsa_cl * fage_cl * fsex_cl * fphase_cl
    q  <- exp(lq) * fbsa_cl

    # ---- De-PEGylation transit chain ---------------------------------------
    # Figure 1 and ESM Section 4 $DES: 14 compartments in series sharing V1.
    # Drug steps down the chain at Qtr/V and every compartment is eliminated at
    # CLinitial/V. The terminal compartment also loses Qtr/V, because the
    # published simplification sets CLinduced = Qtr (Table 3 footnote (a);
    # ESM Section 2, dOFV = 0.6), so every compartment has the same total loss
    # rate CLinitial/V + Qtr/V.
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
    # Total asparaginase activity over all 14 species (ESM Section 4 $ERROR:
    # TOT = A(1)+...+A(14); IPRED = TOT/V1).
    Cc <- (central + transit1 + transit2 + transit3 + transit4 + transit5 +
             transit6 + transit7 + transit8 + transit9 + transit10 +
             transit11 + transit12 + transit13) / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}

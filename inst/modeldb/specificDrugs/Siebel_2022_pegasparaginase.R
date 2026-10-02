Siebel_2022_pegasparaginase <- function() {
  description <- "Fourteen-compartment de-PEGylation transit population PK model for intravenous PEGylated asparaginase (PEG-ASNase) in children with acute lymphoblastic leukemia treated in induction (protocol IA days 12 and 26) and re-induction (protocol II day 8) in the German/Czech part of the AIEOP-BFM ALL 2009 trial, extended with the effect of pre-existing anti-polyethylene-glycol (anti-PEG) IgM antibodies (Siebel 2022, final anti-PEG IgMprior covariate model). Asparaginase activity is carried by a chain of 14 serial compartments sharing one serum volume; drug moves down the chain with intercompartmental clearance Qtr, which mimics stepwise de-PEGylation, and every compartment is eliminated with the initial clearance CLinitial while the terminal compartment is eliminated with CLinitial + Qtr. Body surface area enters volume and the two clearance terms linearly, centred on 0.79 m^2. Initial clearance rises linearly with age above 8 years, is lower in females, and at the first induction dose only rises linearly with the natural-log anti-PEG IgM level above an estimated hockey-stick cut point (1.30 on the log scale, 3.67 mean fluorescence intensity). Initial clearance and volume drop stepwise between administrations, encoded as fractional changes against the first induction dose. Inter-individual variability on initial clearance, inter-occasion variability on both initial clearance and volume, and combined proportional plus additive residual error."
  reference <- "Siebel C, Lanvers-Kaminsky C, Alten J, Smisek P, Nath CE, Rizzari C, Boos J, Wurthwein G. Impact of Antibodies Against Polyethylene Glycol on the Pharmacokinetics of PEGylated Asparaginase in Children with Acute Lymphoblastic Leukaemia: A Population Pharmacokinetic Approach. Eur J Drug Metab Pharmacokinet. 2022;47(2):187-198. doi:10.1007/s13318-021-00741-w"
  vignette <- "Siebel_2022_pegasparaginase"
  units <- list(time = "day", dosing = "U", concentration = "U/L")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 0.78 m^2, range 0.41-2.58 (Siebel 2022 Table 1). BSA enters volume and the two clearance terms LINEARLY and CENTRED on 0.79 m^2, not allometrically: the Siebel 2022 ESM Section 3 control stream writes TVV1 = THETA(1)*(1+SCV1*(BSA-0.79)), TVCLP = THETA(2)*(1+SCCL*(BSA-0.79)) and TVQ = THETA(4)*(1+SCCL*(BSA-0.79)), and Table 3 footnote (a) reads 'Linear increase in V/CLinitial/Qtr with BSA (centred on the median)'. The control stream multiplies nothing else by BSA, so V, CLinitial and Qtr are ABSOLUTE quantities (L, L/day). The 'L/m^2' and 'L/day/m^2' tags in Table 3 record the authors' quoting convention, stated under the table: 'Typical values for V, CLinitial and Qtr are reported for a child with BSA = 1 m^2 for better comparison'. model() therefore renormalises each linear factor to 1 at BSA = 1 m^2 so the ini() values are exactly the printed table entries.",
      source_name = "BSA"
    ),
    AGE = list(
      description = "Age at the time of the PEG-asparaginase administration",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 5.13 years, range 1.06-18.32 (Siebel 2022 Table 1). Enters initial clearance as a hockey stick with the break point fixed at 8 years: the ESM Section 3 control stream writes IF(AGEGTP.LE.8) CLAGEGTP = (1 + THETA(13)*(AGEGTP-8)) with THETA(13) '(0) FIX' and IF(AGEGTP.GT.8) CLAGEGTP = (1 + THETA(14)*(AGEGTP-8)). The column name AGEGTP is age at the time of administration, so the covariate is time-varying across the protocol.",
      source_name = "AGEGTP"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "842 male / 602 female (Siebel 2022 Table 1). The source column is SEX coded 1 = male, 2 = female: the ESM Section 3 control stream writes IF(SEX.EQ.1) CLSEX = 1 ('Most common') and IF(SEX.EQ.2) CLSEX = (1 + THETA(15)). Mapping to the canonical column is SEXF = SEX - 1, which preserves the sign of the coefficient and the male reference, so the printed estimate is used unchanged. Table 3 footnote (d): 'Fractional change in CLinitial for females'.",
      source_name = "SEX (1 = male, 2 = female)"
    ),
    OCC = list(
      description = "Integer administration-occasion index; drives the treatment-phase effects on initial clearance and volume, the inter-occasion variability, and the restriction of the anti-PEG IgM effect to the first induction dose",
      units = "(count)",
      type = "categorical",
      reference_category = "1 (first induction dose, protocol IA day 12) is the model reference for every fractional change",
      notes = "One occasion per PEG-asparaginase administration, exactly as the authors encode it: the ESM Section 3 control stream carries a single OCC column that selects the treatment-phase factors (V1KOV, CLKOV), the IOV etas (IOVV1, IOVCL) and the occasion-specific anti-PEG IgM coefficients (THETA(16)-THETA(18)). The source codes the three administrations 2, 3 and 5; this model renumbers them to the register's 1..N form: OCC = 1 is protocol IA day 12 (induction dose 1, the reference), OCC = 2 is protocol IA day 26 (induction dose 2), OCC = 3 is protocol II day 8 (re-induction, non-high-risk patients only). Records with OCC outside 1-3 receive the reference phase factors, no IOV contribution and no antibody effect.",
      source_name = "OCC (2 = PIA d12, 3 = PIA d26, 5 = PII d8)"
    ),
    ABPEG_IGM = list(
      description = "Pre-existing anti-polyethylene-glycol IgM antibody level in the sample drawn before the first PEG-asparaginase dose of induction",
      units = "(assay-specific mean fluorescence intensity, linear scale)",
      type = "continuous",
      reference_category = NULL,
      notes = "Stored on the LINEAR mean-fluorescence-intensity (MFI) scale; model() takes the natural logarithm, matching the control-stream column IGMPLOG (the paper's cut point 1.30 on the log scale is 3.67 on the linear scale, and exp(1.30) = 3.67; one log unit above it, 2.30, is 9.97). Assayed by the plate-reader method of Khalil et al., which binds antibodies to immobilised methoxy-PEG (5000 Da) on TentaGel particles and reports the duplicate mean fluorescence intensity (ESM Section 1); the readout is an arbitrary assay scale, not a mass concentration. Induction median 1.37 MFI, range 0.19-18.1 (Table 2). Anti-PEG antibodies are PRE-EXISTING (from environmental PEG exposure), which is why they act on the very first dose. Enters initial clearance as a HOCKEY STICK: no effect at or below the estimated cut point, and above it CLIGM = 1 + THETA(16)*(IGMPLOG - CP) (ESM Section 3, CL-IGMPLOG per OCC block). The effect is estimated only at the first induction dose; the control stream fixes the protocol IA day 26 and re-induction coefficients THETA(17) and THETA(18) to 0 because 'the covariate effect could only be estimated with sufficient precision following the first drug administration in induction' (Section 3.2.1). The control stream's IGMPLOG is the level prior to the first dose of the current treatment phase, but because the coefficient is non-zero only at OCC = 1, only the pre-induction value affects the prediction; supply that value on every record. 146 of 1444 induction patients were above the estimated cut point (Discussion).",
      source_name = "IGMPLOG (log of IGMP, anti-PEG IgM prior)"
    )
  )

  # Covariates screened by Siebel 2022 but NOT retained in the final model.
  # Documented for provenance only; deliberately absent from model().
  covariatesDataExcluded <- list(
    ABPEG_IGG = list(
      description = "Pre-existing anti-polyethylene-glycol IgG antibody level in the sample drawn before the first PEG-asparaginase dose",
      units = "(assay-specific mean fluorescence intensity, linear scale)",
      type = "continuous",
      notes = "Screened with the same forms as IgM (categorical, linear on the log scale, per-administration, hockey stick at the Khalil cut point 8 MFI = 2.079 log units). The best IgG-only model (ESM model 1.10) lowered the OFV by 22.899 from the reference model against 49.184 for the retained IgM model, and adding IgG to the IgM hockey stick (ESM model 3.5) gave only dOFV = -4.509 with df = 1 (Section 3.2.1, 'Combination of both antibody isotypes'), so the final model carries IgM alone."
    ),
    ABPEG_IGM_AFTER = list(
      description = "Anti-PEG IgM (and IgG) antibody levels measured after PEG-asparaginase administration, as time-varying covariates or as the after/prior ratio",
      units = "(assay-specific mean fluorescence intensity, linear scale)",
      type = "continuous",
      notes = "Screened as time-varying covariates on CLinitial (absolute log level, log ratio after versus prior, and the Wahlby baseline-plus-change split) and not retained: 'Antibody levels post administration were not associated with an effect on clearance' (Abstract; Section 3.2.1, ESM Section 2). Named here for provenance only; this is not a register canonical."
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
    n_subjects = 1444L,
    n_studies = 1L,
    n_observations = 6261L,
    age_range = "1.06-18.32 years; median 5.13",
    bsa_range = "0.41-2.58 m^2; median 0.78",
    sex_female_pct = 41.7,
    race_ethnicity = "Not reported (German and Czech trial sites of the AIEOP-BFM ALL 2009 trial).",
    disease_state = "Newly diagnosed pediatric acute lymphoblastic leukemia enrolled in the AIEOP-BFM ALL 2009 trial (EudraCT 2007-004270-43), German and Czech part.",
    dose_range = "2500 U/m^2 per dose as a 2-hour intravenous infusion (maximum 3750 U/day), on protocol IA days 12 and 26 (all patients) and protocol II day 8 (non-high-risk patients only). Administered dose median 1950 U (range 720-5300), 2500 U/m^2 (985-4464) (Table 1). 3403 administrations in total.",
    regions = "Germany and Czech Republic.",
    bioanalytic_methods = "Asparaginase activity in serum by the aspartic acid beta-hydroxamate (AHA) assay, lower limit of quantification 5 U/L; samples below it (3.1 percent) were excluded. Anti-PEG IgG and IgM by a plate-reader assay on methoxy-PEG-coated TentaGel particles, reported as mean fluorescence intensity (ESM Section 1).",
    notes = "Demographics from Siebel 2022 Table 1; the female percentage is 602/1444. Samples were scheduled before the first dose of each treatment phase and 7 and 14 days after each dose (ESM Table S1). Anti-PEG antibody levels were available for 2082 samples before and 6412 samples after administration (Table 2). Samples at or after a hypersensitivity reaction or silent inactivation (activity below 100 U/L within 8 days and/or below the LLOQ within 15 days) were EXCLUDED, so the model describes standard elimination only (Section 2.2, Fig. 1). NONMEM 7.4.4, FOCE with interaction; 1000-replicate bootstrap, 97.8 percent successful (Table 3)."
  )

  ini({
    # ---- Structural PK -----------------------------------------------------
    # Siebel 2022 Table 3, 'Final covariate model' column. The table quotes V,
    # CLinitial and Qtr for a child with BSA = 1 m^2; the BSA factors in
    # model() are renormalised to 1 at BSA = 1 m^2 so these are used exactly
    # as printed. See covariateData$BSA.
    lvc <- log(1.71);  label("Volume of the serum compartment shared by all 14 species, for a child with BSA 1 m^2 (L)")                # Table 3 final covariate model (V 1.71, RSE 1.4%; bootstrap 1.71, 95% CI 1.65-1.78)
    lcl <- log(0.128); label("Initial clearance of fully PEGylated asparaginase, for a child with BSA 1 m^2 (L/day)")                  # Table 3 final covariate model (CLinitial 0.128, RSE 1.4%; bootstrap 0.128, 95% CI 0.123-0.133)
    lq  <- log(0.954); label("Intercompartmental clearance driving the de-PEGylation transit chain, for a child with BSA 1 m^2 (L/day)") # Table 3 final covariate model (Qtr 0.954, RSE 2.4%; bootstrap 0.957, 95% CI 0.901-1.01)

    # ---- Body surface area -------------------------------------------------
    # Linear and centred on 0.79 m^2, one slope on V and one shared by
    # CLinitial and Qtr (ESM Section 3: SCV1 = THETA(5), SCCL = THETA(6)).
    e_bsa_vc   <- 1.60; label("Linear slope of BSA on V, centred on 0.79 m^2 (1/m^2)")                 # Table 3 final covariate model (FBSA on V 1.60, RSE 2.1%; bootstrap 1.60, 95% CI 1.54-1.67)
    e_bsa_cl_q <- 1.49; label("Linear slope of BSA on CLinitial and Qtr, centred on 0.79 m^2 (1/m^2)") # Table 3 final covariate model (FBSA on CLinitial + Qtr 1.49, RSE 2.1%; bootstrap 1.49, 95% CI 1.43-1.55)

    # ---- Treatment-phase fractional changes --------------------------------
    # Multiplicative (1 + F) factors relative to the first induction dose,
    # Table 3 footnote (b).
    e_piad26_vc <- -0.158; label("Fractional change in V at the second induction dose, protocol IA day 26 (unitless)")          # Table 3 final covariate model (F(V Induction 2.admin) -0.158, RSE 6.7%; bootstrap -0.157, 95% CI -0.179 to -0.136)
    e_piid8_vc  <- -0.288; label("Fractional change in V at the re-induction dose, protocol II day 8 (unitless)")               # Table 3 final covariate model (F(V Reinduction) -0.288, RSE 3.7%; bootstrap -0.287, 95% CI -0.306 to -0.265)
    e_piad26_cl <- -0.107; label("Fractional change in CLinitial at the second induction dose, protocol IA day 26 (unitless)")  # Table 3 final covariate model (F(CLinitial induction 2.admin) -0.107, RSE 12.8%; bootstrap -0.108, 95% CI -0.135 to -0.078)
    e_piid8_cl  <- -0.438; label("Fractional change in CLinitial at the re-induction dose, protocol II day 8 (unitless)")       # Table 3 final covariate model (F(CLinitial reinduction) -0.438, RSE 3.4%; bootstrap -0.439, 95% CI -0.469 to -0.407)

    # ---- Age and sex on initial clearance ----------------------------------
    e_age_cl  <- 0.011;  label("Linear slope of age above 8 years on CLinitial (1/year)") # Table 3 final covariate model (F(age > 8 years) 0.011, RSE 38.8%; bootstrap 0.011, 95% CI 0.003-0.021)
    e_sexf_cl <- -0.070; label("Fractional change in CLinitial for females (unitless)")   # Table 3 final covariate model (F(sex) -0.070, RSE 23.4%; bootstrap -0.070, 95% CI -0.101 to -0.036)

    # ---- Pre-existing anti-PEG IgM on initial clearance (hockey stick) -----
    # ESM Section 3: CP = THETA(19); IF(OCC.EQ.2.AND.IGMPLOG.GT.CP)
    # CLIGM = (1 + THETA(16)*(IGMPLOG - CP)), with the IA day 26 and
    # re-induction coefficients THETA(17), THETA(18) fixed to 0. The printed
    # cut point is on the natural-log scale, which is exactly the log of the
    # hinge on the linear MFI scale (exp(1.30) = 3.67), so it is stored
    # directly as the log-form hinge parameter.
    e_abpeg_igm_cl   <- 0.414; label("Linear slope of log anti-PEG IgM above the cut point on CLinitial at the first induction dose (per log unit)") # Table 3 final covariate model (F(IgMPrior induction 1.admin > CP) 0.414, RSE 18.0%; bootstrap 0.416, 95% CI 0.237-0.831)
    labpeg_igm_hinge <- 1.30;  label("Cut point of the anti-PEG IgM hockey stick, natural log of mean fluorescence intensity (log MFI)")              # Table 3 final covariate model (CP 1.30 on log scale, RSE 2.5%; bootstrap 1.30, 95% CI 0.984-1.68); Section 3.2.1 'MFI = 3.67 instead of 2; linear scale'

    # ---- Inter-individual variability --------------------------------------
    # Exponential (ESM Section 3: CLP = TVCLP*EXP(ETA(2)+IOVCL)); the printed
    # percentage is taken as the log-normal CV, omega^2 = log(1 + CV^2).
    # No IIV on V or Qtr: $OMEGA fixes both to 0.
    etalcl ~ 0.062520  # Table 3 final covariate model (IIV CLinitial 25.4%, RSE 4.5%; bootstrap 25.4%, 95% CI 22.9-27.7) converted as log(1 + 0.254^2)

    # ---- Inter-occasion variability ----------------------------------------
    # On V and CLinitial, one occasion per administration. NONMEM
    # `$OMEGA BLOCK(1) SAME` is encoded as separate etas with the variance
    # fixed equal to the occasion-1 estimate.
    etaiov_vc_1 ~ 0.014062        # Table 3 final covariate model (IOV V 11.9%, RSE 9.6%; bootstrap 11.8%, 95% CI 9.6-14.5) converted as log(1 + 0.119^2)
    etaiov_vc_2 ~ fixed(0.014062) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_3 ~ fixed(0.014062) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_1 ~ 0.050678        # Table 3 final covariate model (IOV CLinitial 22.8%, RSE 4.9%; bootstrap 22.9%, 95% CI 20.7-25.2) converted as log(1 + 0.228^2)
    etaiov_cl_2 ~ fixed(0.050678) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fixed(0.050678) # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---- Residual error ----------------------------------------------------
    # ESM Section 3 $ERROR: W = SQRT(THETA(7)**2*IPRED**2 + THETA(8)**2),
    # Y = IPRED + W*EPS(1), $SIGMA 1 FIX -> nlmixr2 prop() + add().
    propSd <- 0.191; label("Proportional residual error (fraction)") # Table 3 final covariate model (Proportional error 0.191, RSE 2.8%; bootstrap 0.191, 95% CI 0.180-0.202)
    addSd  <- 5.36;  label("Additive residual error (U/L)")          # Table 3 final covariate model (Additive error 5.36 U/L, RSE 45.1%; bootstrap 5.16, 95% CI 2.23-11.5)
  })

  model({
    # ---- Administration occasion ------------------------------------------
    # OCC = 1 protocol IA day 12 (reference), 2 protocol IA day 26,
    # 3 protocol II day 8.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)

    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 + occ3 * etaiov_vc_3
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3

    # ---- Body surface area -------------------------------------------------
    # Control stream: TVV1 = THETA(1)*(1 + SCV1*(BSA - 0.79)). ini() holds the
    # value printed for BSA = 1 m^2, so each factor is divided by its own value
    # at BSA = 1 m^2 (1 + slope*0.21).
    fbsa_vc <- (1 + e_bsa_vc * (BSA - 0.79)) / (1 + e_bsa_vc * 0.21)
    fbsa_cl <- (1 + e_bsa_cl_q * (BSA - 0.79)) / (1 + e_bsa_cl_q * 0.21)

    # ---- Age, sex, treatment phase and anti-PEG IgM on CLinitial -----------
    fage_cl <- 1 + e_age_cl * (AGE - 8) * (AGE > 8)
    fsex_cl <- 1 + e_sexf_cl * SEXF

    fphase_vc <- 1 + occ2 * e_piad26_vc + occ3 * e_piid8_vc
    fphase_cl <- 1 + occ2 * e_piad26_cl + occ3 * e_piid8_cl

    # Hockey stick on the natural-log IgM level, no effect at or below the
    # cut point, first induction dose only.
    abpeg_igm_log <- log(ABPEG_IGM)
    fabm_cl <- 1 + e_abpeg_igm_cl * max(abpeg_igm_log - labpeg_igm_hinge, 0) * occ1

    # ---- Individual parameters ---------------------------------------------
    vc <- exp(lvc + iov_vc) * fbsa_vc * fphase_vc
    cl <- exp(lcl + etalcl + iov_cl) * fbsa_cl * fage_cl * fsex_cl * fphase_cl * fabm_cl
    q <- exp(lq) * fbsa_cl

    # ---- De-PEGylation transit chain ---------------------------------------
    # ESM Section 3 $DES and Fig. S1: 14 compartments in series sharing V;
    # every compartment is eliminated at CLinitial/V and passes drug on at
    # Qtr/V; the terminal compartment loses its Qtr/V out of the system.
    kel <- cl / vc
    ktr <- q / vc

    d/dt(central)   <-                   -(kel + ktr) * central
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
    # ESM Section 3 $ERROR: TOT = A(1)+...+A(14); IPRED = TOT/V1. The AHA
    # assay measures total catalytic activity of every species in the chain.
    Cc <- (central + transit1 + transit2 + transit3 + transit4 + transit5 +
             transit6 + transit7 + transit8 + transit9 + transit10 +
             transit11 + transit12 + transit13) / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}

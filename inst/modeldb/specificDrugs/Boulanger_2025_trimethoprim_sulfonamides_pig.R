Boulanger_2025_trimethoprim_sulfonamides_pig <- function() {
  description <- paste(
    "Veterinary (pig). Joint fourteen-compartment population PK model for",
    "trimethoprim (TMP) co-administered with each of three sulfonamides -",
    "sulfadiazine (SDZ), sulfadimethoxine (SDMX) and sulfamethoxazole (SMX) -",
    "in growing pigs, fitted simultaneously in Monolix to three cross-over",
    "studies of the licensed 1:5 TMP:S combination products plus re-analysed",
    "raw data from De Smet 2017 (Boulanger 2025). Each of the four drugs has",
    "its own two-compartment disposition with first-order absorption; TMP and",
    "SDZ additionally carry a separate intramuscular depot with its own",
    "absorption rate and bioavailability. Doses are expressed per kilogram, so",
    "clearances are in L/h/kg and volumes in L/kg and concentrations come out",
    "directly in ug/mL. Body weight acts as a power function on CL and V1 of",
    "TMP and SDZ, normalised to the 31.1 kg median pig. Random effects are",
    "inter-occasion (IOV) on CL, V1 and V2 and inter-individual (IIV) on ka;",
    "the IOV terms for CL_SDZ, CL_TMP, V1_SDMX, V1_TMP, V2_SDMX, V2_SDZ and",
    "V2_TMP form a single correlated cross-drug 7x7 block, which is what lets",
    "the model reproduce the paper's central result: the distribution of the",
    "unbound TMP:sulfonamide concentration ratio against the 1:19 target.",
    "TMP uses the unsuffixed canonical compartment / parameter set; the three",
    "sulfonamides carry the sibling-drug suffixes _sdz, _sdmx and _smx."
  )
  reference <- paste(
    "Boulanger M, Taillandier JF, Henri J, Devreese M, De Baere S, Lacroix M,",
    "Ferran AA, Viel A. Population pharmacokinetic modeling of",
    "sulfadimethoxine, sulfadiazine and sulfamethoxazole combined to",
    "trimethoprim in pigs. Vet Q. 2025;45(1):1-17.",
    "doi:10.1080/01652176.2025.2565351.",
    "Fixed effects, random-effect standard deviations and residual error from",
    "Table 2; random-effect correlations from Supplementary Table S6;",
    "protein binding from the Results 'Protein binding experiment' section.",
    sep = " "
  )
  vignette <- "Boulanger_2025_trimethoprim_sulfonamides_pig"
  units    <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")
  # Declared explicitly: buildModelDb()'s dosing heuristic only recognises the
  # literal names `depot` and `central`, so without this field the registry
  # would record just "depot,central" and hide the eight sibling-drug and
  # intramuscular dosing targets. Table 1 doses every one of these: IV into the
  # four central compartments, oral into the four `depot`s, IM into `depot2`
  # (TMP) and `depot2_sdz` (SDZ).
  dosing   <- c(
    "depot", "depot2", "central",
    "depot_sdz", "depot2_sdz", "central_sdz",
    "depot_sdmx", "central_sdmx",
    "depot_smx", "central_smx"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Doses are administered per kilogram of body weight
  # (Table 1 reports every regimen in mg/kg), so the states carry dose-
  # normalised amounts in mg/kg and dividing by the L/kg volumes yields
  # ug/mL directly.
  compartmentData <- list(
    depot            = list(analyte = "trimethoprim",    units = "mg/kg", specimen = "administration site", verified = TRUE),
    depot2           = list(analyte = "trimethoprim",    units = "mg/kg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "trimethoprim",    units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1      = list(analyte = "trimethoprim",    units = "mg/kg", specimen = "plasma", verified = TRUE),
    depot_sdz        = list(analyte = "sulfadiazine",    units = "mg/kg", specimen = "administration site", verified = TRUE),
    depot2_sdz       = list(analyte = "sulfadiazine",    units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_sdz      = list(analyte = "sulfadiazine",    units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_sdz  = list(analyte = "sulfadiazine",    units = "mg/kg", specimen = "plasma", verified = TRUE),
    depot_sdmx       = list(analyte = "sulfadimethoxine", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_sdmx     = list(analyte = "sulfadimethoxine", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_sdmx = list(analyte = "sulfadimethoxine", units = "mg/kg", specimen = "plasma", verified = TRUE),
    depot_smx        = list(analyte = "sulfamethoxazole", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_smx      = list(analyte = "sulfamethoxazole", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_smx  = list(analyte = "sulfamethoxazole", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight of the pig",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Occasion-varying: pigs were re-weighed before each period of the",
        "cross-over (Methods, 'Experimental PK study and blood collection').",
        "Enters as a power function normalised to the 31.1 kg median of all",
        "pigs (Methods Equation 3) on CL and V1 of TMP and SDZ only; no",
        "weight effect was retained for SDMX or SMX. Weights spanned 23.8-46.8",
        "kg across the three combination studies."
      ),
      source_name        = "BW"
    ),
    OCC = list(
      description        = "Occasion index for inter-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Decomposed inside model() into binary indicators oc1/oc2/oc3. An",
        "occasion is one administration period of the cross-over design",
        "(Table 1): OCC = 1 and OCC = 2 are the two cross-over periods, and",
        "OCC = 3 is the additional intramuscular period given to 8 pigs on",
        "the TMP/SDZ combination. The paper reports a single IOV standard",
        "deviation per parameter, so the three occasion slots share one",
        "value (occasions 2 and 3 are fixed() to the occasion-1 estimate,",
        "the equivalent of NONMEM $OMEGA BLOCK SAME). Records taken before",
        "any dose may carry OCC = 0, which zeroes every indicator."
      ),
      source_name        = "occasion"
    )
  )

  population <- list(
    species        = "pig (Large White x Landrace)",
    n_subjects     = 34,
    n_studies      = 4,
    age_range      = "7-8 weeks at enrolment",
    weight_range   = "23.8-46.8 kg",
    weight_median  = "31.1 kg",
    sex_female_pct = 100,
    disease_state  = "healthy growing pigs",
    dose_range     = paste(
      "single dose of licensed 1:5 TMP:S products - TMP/SDZ 2.5+12.5 mg/kg IV",
      "and IM, 5+25 mg/kg oral; TMP/SDMX 4+18.6 mg/kg IV, 8+37.36 mg/kg oral;",
      "TMP/SMX 6+30 mg/kg IV and oral (Table 1). Pooled De Smet 2017 data add",
      "multiple-dose oral (25+5 or 12.5+2.5 mg/kg q12h x 5 d) and IM",
      "(12.5+2.5 or 25+5 mg/kg q24h x 5 d) SDZ/TMP records."
    ),
    regions        = "France and Belgium",
    notes          = paste(
      "Three independent two-period cross-over studies with 10 (TMP/SDZ), 10",
      "(TMP/SMX) and 14 (TMP/SDMX) pigs; all animals female, housed",
      "individually after jugular catheterisation (Methods, 'Animals' and",
      "Table 1). Raw individual data from De Smet et al. 2017 were pooled in",
      "to support the SDZ and TMP models. Of the 1858 concentrations included",
      "in the modelling, 155 were considered outliers, as were 3 animals (1",
      "on TMP/SDMX, 2 on TMP/SMX), mainly for catheter problems; SMX",
      "concentrations beyond 12 h were discarded as analytically implausible.",
      "Mean in vivo plasma protein binding measured in a subset after IV",
      "dosing was 29.2% (SDZ), 57.3% (SMX), 94.1% (SDMX) and 51.2% (TMP),",
      "i.e. unbound fractions of 0.708, 0.427, 0.059 and 0.489 respectively",
      "(Results, 'Protein binding experiment')."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Fixed effects - Table 2, 'Fixed effects', median of the bootstrap
    # (n = 100 runs). Clearances are L/h/kg and volumes L/kg because the
    # doses in Table 1 are expressed per kilogram of body weight.
    # ------------------------------------------------------------------

    # --- Trimethoprim (unsuffixed parent) ---
    lka          <- log(0.66)  ; label("Oral absorption rate constant, trimethoprim (ka, 1/h)")                 # Table 2, ka_TMP
    lka2         <- log(1.14)  ; label("Intramuscular absorption rate constant, trimethoprim (ka, 1/h)")        # Table 2, ka_TMP_IM
    lfdepot      <- log(0.58)  ; label("Oral bioavailability, trimethoprim (F, fraction)")                      # Table 2, F_TMP_oral = 58%
    lfdepot2     <- log(0.99)  ; label("Intramuscular bioavailability, trimethoprim (F, fraction)")             # Table 2, F_TMP_IM = 99%
    lcl          <- log(0.48)  ; label("Clearance, trimethoprim (CL, L/h/kg)")                                  # Table 2, Cl_TMP
    lvc          <- log(0.92)  ; label("Central volume, trimethoprim (V1, L/kg)")                               # Table 2, V1_TMP
    lq           <- log(1.26)  ; label("Intercompartmental clearance, trimethoprim (Q, L/h/kg)")                # Table 2, Q_TMP
    lvp          <- log(0.86)  ; label("Peripheral volume, trimethoprim (V2, L/kg)")                            # Table 2, V2_TMP

    # --- Sulfadiazine ---
    lka_sdz      <- log(0.50)  ; label("Oral absorption rate constant, sulfadiazine (ka, 1/h)")                 # Table 2, ka_SDZ
    lka2_sdz     <- log(1.3)   ; label("Intramuscular absorption rate constant, sulfadiazine (ka, 1/h)")        # Table 2, ka_SDZ_IM
    lfdepot_sdz  <- log(0.93)  ; label("Oral bioavailability, sulfadiazine (F, fraction)")                      # Table 2, F_SDZ_oral = 93%
    lfdepot2_sdz <- log(0.99)  ; label("Intramuscular bioavailability, sulfadiazine (F, fraction)")             # Table 2, F_SDZ_IM = 99%
    lcl_sdz      <- log(0.12)  ; label("Clearance, sulfadiazine (CL, L/h/kg)")                                  # Table 2, Cl_SDZ
    lvc_sdz      <- log(0.3)   ; label("Central volume, sulfadiazine (V1, L/kg)")                               # Table 2, V1_SDZ
    lq_sdz       <- log(0.32)  ; label("Intercompartmental clearance, sulfadiazine (Q, L/h/kg)")                # Table 2, Q_SDZ
    lvp_sdz      <- log(0.29)  ; label("Peripheral volume, sulfadiazine (V2, L/kg)")                            # Table 2, V2_SDZ

    # --- Sulfadimethoxine ---
    lka_sdmx     <- log(0.60)  ; label("Oral absorption rate constant, sulfadimethoxine (ka, 1/h)")             # Table 2, ka_SDMX
    lfdepot_sdmx <- log(0.68)  ; label("Oral bioavailability, sulfadimethoxine (F, fraction)")                  # Table 2, F_SDMX_oral = 68%
    lcl_sdmx     <- log(0.015) ; label("Clearance, sulfadimethoxine (CL, L/h/kg)")                              # Table 2, Cl_SDMX
    lvc_sdmx     <- log(0.13)  ; label("Central volume, sulfadimethoxine (V1, L/kg)")                           # Table 2, V1_SDMX
    lq_sdmx      <- log(0.20)  ; label("Intercompartmental clearance, sulfadimethoxine (Q, L/h/kg)")            # Table 2, Q_SDMX
    lvp_sdmx     <- log(0.16)  ; label("Peripheral volume, sulfadimethoxine (V2, L/kg)")                        # Table 2, V2_SDMX

    # --- Sulfamethoxazole ---
    lka_smx      <- log(1.82)  ; label("Oral absorption rate constant, sulfamethoxazole (ka, 1/h)")             # Table 2, ka_SMX
    lfdepot_smx  <- log(0.64)  ; label("Oral bioavailability, sulfamethoxazole (F, fraction)")                  # Table 2, F_SMX_oral = 64%
    lcl_smx      <- log(0.21)  ; label("Clearance, sulfamethoxazole (CL, L/h/kg)")                              # Table 2, Cl_SMX
    lvc_smx      <- log(0.48)  ; label("Central volume, sulfamethoxazole (V1, L/kg)")                           # Table 2, V1_SMX
    lq_smx       <- log(0.83)  ; label("Intercompartmental clearance, sulfamethoxazole (Q, L/h/kg)")            # Table 2, Q_SMX
    lvp_smx      <- log(0.17)  ; label("Peripheral volume, sulfamethoxazole (V2, L/kg)")                        # Table 2, V2_SMX

    # ------------------------------------------------------------------
    # Covariate effects - Methods Equation 3, power model on body weight
    # normalised to the 31.1 kg median pig. Only TMP and SDZ retained a
    # weight effect, and only on CL and V1.
    # ------------------------------------------------------------------
    e_wt_cl      <- 0.78       ; label("Power exponent on (WT/31.1) for CL, trimethoprim (unitless)")           # Table 2, beta_Cl_TMP_BW
    e_wt_vc      <- 1.33       ; label("Power exponent on (WT/31.1) for V1, trimethoprim (unitless)")           # Table 2, beta_V1_TMP_BW
    e_wt_cl_sdz  <- 0.51       ; label("Power exponent on (WT/31.1) for CL, sulfadiazine (unitless)")           # Table 2, beta_Cl_SDZ_BW
    e_wt_vc_sdz  <- 1.1        ; label("Power exponent on (WT/31.1) for V1, sulfadiazine (unitless)")           # Table 2, beta_V1_SDZ_BW

    # ------------------------------------------------------------------
    # Inter-individual variability on the absorption rate constants.
    # Table 2 reports STANDARD DEVIATIONS of the random effects (column
    # header) and Monolix reports omega on the SD scale, so each entry
    # below is the square of the tabulated value.
    # IOV was not implemented for ka because each pig received the oral
    # solution only once (Methods, 'Population PK modeling').
    # ------------------------------------------------------------------
    etalka       ~ 0.3969   # omega_ka_TMP  = 0.63 -> 0.63^2
    etalka_sdz   ~ 0.2704   # omega_ka_SDZ  = 0.52 -> 0.52^2
    etalka_sdmx  ~ 0.1369   # omega_ka_SDMX = 0.37 -> 0.37^2
    # No eta on ka_SMX: IIV for ka of SMX was fixed to 0 (Results,
    # 'Population PK analysis'). Likewise no eta on any Q: the IIV and IOV
    # for Q were fixed to 0 for all four drugs.

    # ------------------------------------------------------------------
    # Inter-occasion variability, correlated cross-drug block.
    # Supplementary Table S6 reports 21 pairwise correlations, which is
    # exactly the number of off-diagonal pairs of a 7x7 matrix over
    # CL_SDZ, CL_TMP, V1_SDMX, V1_TMP, V2_SDMX, V2_SDZ and V2_TMP. The
    # covariances below are corr(i,j) * sd_i * sd_j using the Table 2
    # standard deviations 0.33, 0.34, 0.40, 0.44, 0.34, 0.62 and 0.42.
    # The resulting matrix is positive definite (min eigenvalue 0.0124).
    #
    # Occasion order within the block:
    #   1 CL_SDZ   2 CL_TMP   3 V1_SDMX   4 V1_TMP
    #   5 V2_SDMX  6 V2_SDZ   7 V2_TMP
    # Occasions 2 and 3 repeat the occasion-1 block (one IOV magnitude per
    # parameter in Table 2), so they are fixed() to the same values.
    # ------------------------------------------------------------------
    etaiov_cl_sdz_1 + etaiov_cl_1 + etaiov_vc_sdmx_1 + etaiov_vc_1 +
      etaiov_vp_sdmx_1 + etaiov_vp_sdz_1 + etaiov_vp_1 ~ c(
        0.1089,
        0.030294, 0.1156,
        0.0792, 0.05304, 0.16,
        0.0726, 0.113696, 0.1144, 0.1936,
        -0.060588, -0.071672, -0.09384, -0.101728, 0.1156,
        0.178002, 0.033728, 0.14632, 0.130944, -0.096968, 0.3844,
        0.026334, 0.032844, 0.04872, 0.053592, -0.022848, 0.109368, 0.1764
      )
    etaiov_cl_sdz_2 + etaiov_cl_2 + etaiov_vc_sdmx_2 + etaiov_vc_2 +
      etaiov_vp_sdmx_2 + etaiov_vp_sdz_2 + etaiov_vp_2 ~ fixed(c(
        0.1089,
        0.030294, 0.1156,
        0.0792, 0.05304, 0.16,
        0.0726, 0.113696, 0.1144, 0.1936,
        -0.060588, -0.071672, -0.09384, -0.101728, 0.1156,
        0.178002, 0.033728, 0.14632, 0.130944, -0.096968, 0.3844,
        0.026334, 0.032844, 0.04872, 0.053592, -0.022848, 0.109368, 0.1764
      ))
    etaiov_cl_sdz_3 + etaiov_cl_3 + etaiov_vc_sdmx_3 + etaiov_vc_3 +
      etaiov_vp_sdmx_3 + etaiov_vp_sdz_3 + etaiov_vp_3 ~ fixed(c(
        0.1089,
        0.030294, 0.1156,
        0.0792, 0.05304, 0.16,
        0.0726, 0.113696, 0.1144, 0.1936,
        -0.060588, -0.071672, -0.09384, -0.101728, 0.1156,
        0.178002, 0.033728, 0.14632, 0.130944, -0.096968, 0.3844,
        0.026334, 0.032844, 0.04872, 0.053592, -0.022848, 0.109368, 0.1764
      ))

    # ------------------------------------------------------------------
    # Inter-occasion variability, uncorrelated terms. These four random
    # effects appear in Table 2 but not in Supplementary Table S6, so they
    # are independent of the correlated block and of each other. Values
    # are the squares of the tabulated standard deviations.
    # ------------------------------------------------------------------
    etaiov_cl_sdmx_1 ~ 0.0324          # gamma_Cl_SDMX = 0.18 -> 0.18^2
    etaiov_cl_sdmx_2 ~ fixed(0.0324)
    etaiov_cl_sdmx_3 ~ fixed(0.0324)
    etaiov_vc_sdz_1  ~ 0.0144          # gamma_V1_SDZ  = 0.12 -> 0.12^2
    etaiov_vc_sdz_2  ~ fixed(0.0144)
    etaiov_vc_sdz_3  ~ fixed(0.0144)
    etaiov_vc_smx_1  ~ 0.0081          # gamma_V1_SMX  = 0.09 -> 0.09^2
    etaiov_vc_smx_2  ~ fixed(0.0081)
    etaiov_vc_smx_3  ~ fixed(0.0081)
    etaiov_vp_smx_1  ~ 0.16            # gamma_V2_SMX  = 0.40 -> 0.40^2
    etaiov_vp_smx_2  ~ fixed(0.16)
    etaiov_vp_smx_3  ~ fixed(0.16)

    # ------------------------------------------------------------------
    # Residual error. Concentrations were log transformed and the residual
    # error was best described by a constant model for all drugs (Methods,
    # 'Population PK modeling'), i.e. a constant SD on the log scale, which
    # is the lnorm / expSd structure here.
    # ------------------------------------------------------------------
    expSd      <- 0.40 ; label("Log-scale residual SD, trimethoprim (unitless)")      # Table 2, a_TMP
    expSd_sdz  <- 0.32 ; label("Log-scale residual SD, sulfadiazine (unitless)")      # Table 2, a_SDZ
    expSd_sdmx <- 0.22 ; label("Log-scale residual SD, sulfadimethoxine (unitless)")  # Table 2, a_SDMX
    expSd_smx  <- 0.18 ; label("Log-scale residual SD, sulfamethoxazole (unitless)")  # Table 2, a_SMX
  })

  model({
    # 1. Occasion indicators (Methods Equation 2). An occasion is one
    #    administration period of the cross-over design.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    # 2. Inter-occasion random-effect terms. Kept on their own lines so the
    #    structural parameters below stay in mu-referenced form.
    iov_cl        <- oc1 * etaiov_cl_1        + oc2 * etaiov_cl_2        + oc3 * etaiov_cl_3
    iov_vc        <- oc1 * etaiov_vc_1        + oc2 * etaiov_vc_2        + oc3 * etaiov_vc_3
    iov_vp        <- oc1 * etaiov_vp_1        + oc2 * etaiov_vp_2        + oc3 * etaiov_vp_3
    iov_cl_sdz    <- oc1 * etaiov_cl_sdz_1    + oc2 * etaiov_cl_sdz_2    + oc3 * etaiov_cl_sdz_3
    iov_vc_sdz    <- oc1 * etaiov_vc_sdz_1    + oc2 * etaiov_vc_sdz_2    + oc3 * etaiov_vc_sdz_3
    iov_vp_sdz    <- oc1 * etaiov_vp_sdz_1    + oc2 * etaiov_vp_sdz_2    + oc3 * etaiov_vp_sdz_3
    iov_cl_sdmx   <- oc1 * etaiov_cl_sdmx_1   + oc2 * etaiov_cl_sdmx_2   + oc3 * etaiov_cl_sdmx_3
    iov_vc_sdmx   <- oc1 * etaiov_vc_sdmx_1   + oc2 * etaiov_vc_sdmx_2   + oc3 * etaiov_vc_sdmx_3
    iov_vp_sdmx   <- oc1 * etaiov_vp_sdmx_1   + oc2 * etaiov_vp_sdmx_2   + oc3 * etaiov_vp_sdmx_3
    iov_vc_smx    <- oc1 * etaiov_vc_smx_1    + oc2 * etaiov_vc_smx_2    + oc3 * etaiov_vc_smx_3
    iov_vp_smx    <- oc1 * etaiov_vp_smx_1    + oc2 * etaiov_vp_smx_2    + oc3 * etaiov_vp_smx_3

    # 3. Individual parameters. Body weight enters as a power function
    #    normalised to the 31.1 kg median pig (Methods Equation 3).
    ka        <- exp(lka + etalka)
    ka2       <- exp(lka2)
    cl        <- exp(lcl + iov_cl) * (WT / 31.1)^e_wt_cl
    vc        <- exp(lvc + iov_vc) * (WT / 31.1)^e_wt_vc
    q         <- exp(lq)
    vp        <- exp(lvp + iov_vp)

    ka_sdz    <- exp(lka_sdz + etalka_sdz)
    ka2_sdz   <- exp(lka2_sdz)
    cl_sdz    <- exp(lcl_sdz + iov_cl_sdz) * (WT / 31.1)^e_wt_cl_sdz
    vc_sdz    <- exp(lvc_sdz + iov_vc_sdz) * (WT / 31.1)^e_wt_vc_sdz
    q_sdz     <- exp(lq_sdz)
    vp_sdz    <- exp(lvp_sdz + iov_vp_sdz)

    ka_sdmx   <- exp(lka_sdmx + etalka_sdmx)
    cl_sdmx   <- exp(lcl_sdmx + iov_cl_sdmx)
    vc_sdmx   <- exp(lvc_sdmx + iov_vc_sdmx)
    q_sdmx    <- exp(lq_sdmx)
    vp_sdmx   <- exp(lvp_sdmx + iov_vp_sdmx)

    ka_smx    <- exp(lka_smx)
    cl_smx    <- exp(lcl_smx)
    vc_smx    <- exp(lvc_smx + iov_vc_smx)
    q_smx     <- exp(lq_smx)
    vp_smx    <- exp(lvp_smx + iov_vp_smx)

    # 4. ODE system. Two compartments per drug; TMP and SDZ each carry a
    #    separate intramuscular depot (depot2) with its own ka and F.
    d/dt(depot)            <- -ka * depot
    d/dt(depot2)           <- -ka2 * depot2
    d/dt(central)          <-  ka * depot + ka2 * depot2 -
      (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1)      <-  (q / vc) * central - (q / vp) * peripheral1

    d/dt(depot_sdz)        <- -ka_sdz * depot_sdz
    d/dt(depot2_sdz)       <- -ka2_sdz * depot2_sdz
    d/dt(central_sdz)      <-  ka_sdz * depot_sdz + ka2_sdz * depot2_sdz -
      (cl_sdz / vc_sdz) * central_sdz - (q_sdz / vc_sdz) * central_sdz +
      (q_sdz / vp_sdz) * peripheral1_sdz
    d/dt(peripheral1_sdz)  <-  (q_sdz / vc_sdz) * central_sdz -
      (q_sdz / vp_sdz) * peripheral1_sdz

    d/dt(depot_sdmx)       <- -ka_sdmx * depot_sdmx
    d/dt(central_sdmx)     <-  ka_sdmx * depot_sdmx -
      (cl_sdmx / vc_sdmx) * central_sdmx - (q_sdmx / vc_sdmx) * central_sdmx +
      (q_sdmx / vp_sdmx) * peripheral1_sdmx
    d/dt(peripheral1_sdmx) <-  (q_sdmx / vc_sdmx) * central_sdmx -
      (q_sdmx / vp_sdmx) * peripheral1_sdmx

    d/dt(depot_smx)        <- -ka_smx * depot_smx
    d/dt(central_smx)      <-  ka_smx * depot_smx -
      (cl_smx / vc_smx) * central_smx - (q_smx / vc_smx) * central_smx +
      (q_smx / vp_smx) * peripheral1_smx
    d/dt(peripheral1_smx)  <-  (q_smx / vc_smx) * central_smx -
      (q_smx / vp_smx) * peripheral1_smx

    # 5. Bioavailability. IV doses go straight to the central compartments
    #    and are complete by definition.
    f(depot)      <- exp(lfdepot)
    f(depot2)     <- exp(lfdepot2)
    f(depot_sdz)  <- exp(lfdepot_sdz)
    f(depot2_sdz) <- exp(lfdepot2_sdz)
    f(depot_sdmx) <- exp(lfdepot_sdmx)
    f(depot_smx)  <- exp(lfdepot_smx)

    # 6. Observations. Volumes are L/kg and the dose amounts are mg/kg, so
    #    each ratio is already mg/L = ug/mL.
    Cc      <- central / vc
    Cc_sdz  <- central_sdz / vc_sdz
    Cc_sdmx <- central_sdmx / vc_sdmx
    Cc_smx  <- central_smx / vc_smx

    Cc      ~ lnorm(expSd)
    Cc_sdz  ~ lnorm(expSd_sdz)
    Cc_sdmx ~ lnorm(expSd_sdmx)
    Cc_smx  ~ lnorm(expSd_smx)
  })
}

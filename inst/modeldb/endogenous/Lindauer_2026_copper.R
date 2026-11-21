Lindauer_2026_copper <- function() {
  description <- paste(
    "QSP (mechanistic 14-compartment linear ODE). Whole-body copper kinetics",
    "in healthy adults, heterozygous ATP7B carriers and Wilson disease (WD)",
    "patients, fitted to radiocopper (64Cu / 67Cu / 65Cu) plasma, ceruloplasmin-",
    "bound plasma, liver, faecal and urinary data. Distinguishes",
    "ceruloplasmin-bound copper from non-ceruloplasmin-bound copper (NCC) and",
    "carries the two ATP7B-dependent processes (biliary excretion and",
    "ceruloplasmin loading) as activity factors: 1 in healthy subjects,",
    "estimated in heterozygotes, 0 in untreated WD. The ATP7B gene therapy",
    "VTX-801 enters as a restored activity factor in WD patients, dosed into a",
    "constant state as in the authors' simulation model. Chelator and",
    "zinc background therapy and a COMET-study faecal-transit shift are",
    "covariates. The fitted model has no IIV; the fixed 15% CV IIV on the rate",
    "constants is the authors' simulation assumption.",
    sep = " "
  )
  reference <- paste(
    "Lindauer A, Benichou B, Gonzalez Aseguinolaza G, Combal JP. From",
    "Radiocopper to Cold Copper: Mechanistic Modeling and Simulation to Define",
    "Clinical Response Criteria and Biomarkers for VTX-801 in Wilson Disease.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70153.",
    "doi:10.1002/psp4.70153. NONMEM control stream in Supporting Information",
    "Data S1 section 3; campsis simulation model in Data S2.",
    sep = " "
  )
  vignette <- "Lindauer_2026_copper"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Copper doses go into gut_abs (oral / dietary) or plasma_ncc (intravenous
  # radiocopper). trt_af takes the unitless VTX-801 activity factor (0-1) as a
  # bolus, not a copper amount.
  dosing <- c("gut_abs", "plasma_ncc", "trt_af")

  # Every state is a copper amount. Radiocopper studies dose a unit (or 100%)
  # tracer dose, so amounts read directly as a fraction (percent) of the
  # injected dose; cold-copper simulations dose mg of dietary copper.
  paper_specific_compartments <- c(
    "gut_abs",
    "gut_distal",
    "liver_storage",
    "liver_precursor",
    "plasma_ncc",
    "plasma_cp",
    "tissue_ncc",
    "tissue_deep",
    "tissue_cp",
    "gut_transit",
    "feces",
    "sweat",
    "trt_af"
  )

  compartmentData <- list(
    gut_abs = list(
      analyte = "copper",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    gut_distal = list(
      analyte = "copper",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    liver = list(
      analyte = "copper (hepatic extracellular / vascular space)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    liver_storage = list(
      analyte = "copper (hepatocellular storage pool)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    liver_precursor = list(
      analyte = "copper (hepatocellular precursor pool)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    plasma_ncc = list(
      analyte = "non-ceruloplasmin-bound copper",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    plasma_cp = list(
      analyte = "ceruloplasmin-bound copper",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    tissue_ncc = list(
      analyte = "non-ceruloplasmin-bound copper",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    tissue_deep = list(
      analyte = "non-ceruloplasmin-bound copper",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    tissue_cp = list(
      analyte = "ceruloplasmin-bound copper",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    gut_transit = list(
      analyte = "copper",
      units = "mg",
      specimen = "faeces",
      verified = TRUE
    ),
    feces = list(
      analyte = "copper (cumulative faecal excretion)",
      units = "mg",
      specimen = "faeces",
      verified = TRUE
    ),
    urine = list(
      analyte = "copper (cumulative urinary excretion)",
      units = "mg",
      specimen = "urine",
      verified = TRUE
    ),
    sweat = list(
      analyte = "copper (cumulative loss in sweat)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    trt_af = list(
      analyte = "ATP7B activity factor restored by VTX-801 gene therapy",
      units = "(fraction of normal ATP7B activity)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    DIS_WILSON = list(
      description = "Wilson disease patient indicator (loss of both functional ATP7B alleles)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = Wilson disease patient (NONMEM isWD = 1). Sets both activity",
        "factors from the gene-therapy state trt_af (0 if untreated),",
        "switches on k_livPre_liv, and applies the WD shifts on k_liv_livPre,",
        "k_tis_pl and k_tis2_pl. 0 = healthy subject or heterozygous carrier.",
        "Mutually exclusive with ATP7B_HET."
      ),
      source_name = "isWD"
    ),
    ATP7B_HET = list(
      description = "Heterozygous ATP7B carrier indicator (one functional allele; Wilson disease relative)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = heterozygous carrier (NONMEM isWD = 0.5): activity factors",
        "af_het_cp and af_het_bil replace 1, and no WD distribution shifts",
        "apply. 0 = healthy non-carrier or WD patient. Mutually exclusive",
        "with DIS_WILSON."
      ),
      source_name = "isWD == 0.5 (campsis isHet)"
    ),
    CONMED_COPPER_CHELATOR = list(
      description = "Copper-chelator background therapy (penicillamine, trientine or dimercaprol)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = on penicillamine, trientine or dimercaprol (BAL), pooled as one",
        "class effect: k_pl_urine * (1 + e_chelator_k_pl_urine). Time-varying",
        "allowed."
      ),
      source_name = "onChelator"
    ),
    CONMED_ZINC = list(
      description = "Zinc background therapy",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = on zinc: k_abs_liv * (1 + e_zinc_k_abs_liv) and",
        "k_pl_tis * (1 + e_zinc_k_pl_tis). Time-varying allowed."
      ),
      source_name = "onZinc"
    ),
    STUDY_COMET = list(
      description = "COMET study (EudraCT 2019-001157-13) participant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = COMET participant: k_gutD_feces * (1 + e_comet_k_gutd_feces), a",
        "60.7% slower faecal transit than in the literature data, and the",
        "COMET plasma residual error propSd_comet on Cc. Set 0 for",
        "simulation of new patients, as the paper does."
      ),
      source_name = "isCOMET"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 221L,
    n_studies = 19L,
    age_range = "adults (not reported)",
    disease_state = paste(
      "Healthy volunteers, heterozygous ATP7B carriers and Wilson disease",
      "patients (untreated, on zinc, or on chelators)."
    ),
    dose_range = paste(
      "Single tracer doses of 64Cu, 67Cu or 65Cu, intravenous or oral",
      "(literature); single intravenous [64Cu]CuCl2 in COMET."
    ),
    regions = "Europe and North America (literature since 1954); Austria (COMET)",
    notes = paste(
      "Digitised mean / median / individual time-activity data from 18",
      "published studies (Data S1 Table S1 lists 215 digitised subjects: 62",
      "untreated WD, 7 WD on zinc, 16 WD on chelators, 36 heterozygotes, 94",
      "healthy; four further PET / scintillation sources are listed",
      "without counts) plus individual data from the COMET study (3 healthy",
      "adults and 3 heterozygous carriers, 68 plasma samples above the limit",
      "of detection up to 96 h). Only intravenous data were used in the final",
      "fit; oral absorption parameters were fixed from an earlier model",
      "version. Fitted in NONMEM 7.4.3 with first-order estimation and no IIV."
    )
  )

  ini({
    # All rate constants are on the log scale (Table 1 footnote: 'All rate
    # constants were estimated on the log-scale'). Fixed values use the exact
    # log-scale value from the NONMEM $THETA block (Data S1 section 3), which
    # Table 1 prints rounded.
    lk_abs_liv <- fixed(1.79)
    label("Log absorption rate, absorption compartment to liver (1/day)") # NONMEM THETA(1) 1.79 FIX; Table 1 k_abs_liv = 6 (fixed)
    lk_abs_gutd <- fixed(2.64)
    label("Log transfer rate, absorption compartment to distal gut (1/day)") # NONMEM THETA(2) 2.64 FIX; Table 1 k_abs_gutD = 14.0 (fixed)
    lk_livpre_out <- log(1.13)
    label("Log biliary elimination rate from the hepatic precursor pool (1/day)") # Table 1 k_livPre_out = 1.13
    f_bile_out <- fixed(0.75)
    label("Fraction of biliary copper NOT reabsorbed (unitless)") # Table 1 E_bill = 0.75 (fixed, Buckley); NONMEM THETA(23)
    lk_liv_livpre <- log(169)
    label("Log rate, extracellular liver to precursor pool (1/day)") # Table 1 k_liv_livPre = 169
    lk_livpre_liv <- log(4.9)
    label("Log rate, precursor pool to extracellular liver; WD patients only (1/day)") # Table 1 k_livPre_liv = 4.9 (fixed to 0 in non-WD)
    lk_livpre_plcp <- log(0.516)
    label("Log rate, precursor pool to ceruloplasmin-bound plasma copper (1/day)") # Table 1 k_livPre_plCP = 0.516
    lk_livpre_livstore <- log(5.26)
    label("Log rate, precursor to storage pool (1/day)") # Table 1 k_livPre_livStore = 5.26
    lk_livstore_livpre <- log(0.711)
    label("Log rate, storage to precursor pool (1/day)") # Table 1 k_livStore_livPre = 0.711
    lk_gutd_feces <- log(1.00)
    label("Log elimination rate, distal gut to faeces (1/day)") # Table 1 k_gutD_feces = 1.00
    lk_pl_liv <- fixed(6.13)
    label("Log plasma to liver distribution rate of NCC (1/day)") # NONMEM THETA(10) 6.13 FIX; Table 1 k_pl_liv = 459 (fixed to physiological value)
    lk_liv_pl <- log(1339)
    label("Log liver to plasma distribution rate of NCC (1/day)") # Table 1 k_liv_pl = 1339
    lk_pl_tis <- log(91.8)
    label("Log plasma to tissue distribution rate of NCC (1/day)") # Table 1 k_pl_tis = 91.8
    lk_tis_pl <- log(29.7)
    label("Log tissue to plasma distribution rate of NCC (1/day)") # Table 1 k_tis_pl = 29.7
    lk_pl_tis2 <- log(42.9)
    label("Log plasma to deep tissue distribution rate of NCC (1/day)") # Table 1 k_pl_tis2 = 42.9
    lk_tis2_pl <- log(0.0863)
    label("Log deep tissue to plasma distribution rate of NCC (1/day)") # Table 1 k_tis2_pl = 0.0863
    lk_plcp_tiscp <- fixed(-0.968)
    label("Log plasma to tissue distribution rate of ceruloplasmin-bound copper (1/day)") # NONMEM THETA(20) -0.968 FIX; Table 1 k_plCP_tisCP = 0.380 (Waldmann)
    lk_tiscp_plcp <- fixed(-0.223)
    label("Log tissue to plasma distribution rate of ceruloplasmin-bound copper (1/day)") # NONMEM THETA(21) -0.223 FIX; Table 1 k_tisCP_plCP = 0.800 (Waldmann)
    lk_plcp_liv <- fixed(-2.07)
    label("Log conversion rate of ceruloplasmin-bound copper to NCC in the liver (1/day)") # NONMEM THETA(16) -2.07 FIX; Table 1 k_plCP_liv = 0.126 (CP half-life ~5.5 days)
    lk_tis_sweat <- fixed(-5.42)
    label("Log tissue to sweat excretion rate of NCC (1/day)") # NONMEM THETA(17) -5.42 FIX; Table 1 k_tis_sweat = 0.00443 (Buckley)
    lk_pl_urine <- log(0.411)
    label("Log plasma to urine excretion rate of NCC (1/day)") # Table 1 k_pl_urine = 0.411
    lk_tis_abs <- log(0.983)
    label("Log elimination rate, tissue to the absorption compartment (1/day)") # Table 1 k_tis_abs = 0.983
    lvpl <- fixed(log(2.8))
    label("Log plasma volume (L)") # Methods 'average plasma volume of 2800 mL'; Data S2 theta.csv V_pl = 2800 mL

    # ATP7B activity factors
    af_het_cp <- 0.369
    label("Ceruloplasmin-loading activity factor in heterozygous carriers (unitless)") # Table 1 AFcp = 0.369
    af_het_bil <- 0.101
    label("Biliary-excretion activity factor in heterozygous carriers (unitless)") # Table 1 AFbil = 0.101
    e_af_bil <- fixed(2.28)
    label("Exponent linking biliary to ceruloplasmin-loading activity under gene therapy (unitless)") # Data S2 model.campsis: AFbill = AF_WD_BILL + pow(AFcp, 2.28)

    # Covariate effects (fractional changes)
    e_chelator_k_pl_urine <- 15.8
    label("Fractional increase of k_pl_urine on chelators (unitless)") # Table 1 Chel_eff = 15.8 (16.8-fold)
    e_wilson_k_liv_livpre <- -0.153
    label("Fractional change of k_liv_livPre in WD patients (unitless)") # Table 1 WD_liv_livPre = -0.153
    e_wilson_k_tis_pl <- 1.43
    label("Fractional change of k_tis_pl in WD patients (unitless)") # Table 1 WD_tp = 1.43 (2.43-fold)
    e_wilson_k_tis2_pl <- 26.6
    label("Fractional change of k_tis2_pl in WD patients (unitless)") # Table 1 WD_t2p = 26.6 (27.6-fold)
    e_comet_k_gutd_feces <- -0.607
    label("Fractional change of k_gutD_feces in COMET subjects (unitless)") # Table 1 COMET_eff_gut = -0.607
    # Table 1 prints Zinc_eff_abs = -0.401 and Zinc_eff_dist = 1.2, but the
    # NONMEM control stream that produced Table 1 fixes them at -0.8 and 1, and
    # the campsis simulation model uses -0.8 and 1; the code values are used.
    e_zinc_k_abs_liv <- fixed(-0.8)
    label("Fractional change of k_abs_liv on zinc (unitless)") # NONMEM THETA(26) -0.8 FIX; Data S2 theta.csv Zinc_eff_abs = -0.8
    e_zinc_k_pl_tis <- fixed(1)
    label("Fractional change of k_pl_tis on zinc (unitless)") # NONMEM THETA(27) 1 FIX; Data S2 theta.csv Zinc_eff_dist = 1

    # IIV was not estimated ($OMEGA 0 FIX). The paper's simulations set 15% CV
    # on the rate constants (Methods, Simulations); the eta set follows the
    # campsis model (Data S2 omega.csv: every rate constant except k_plCP_liv,
    # plus V_pl), with variance 0.15^2 on the campsis SD scale.
    etalk_abs_liv ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_abs_gutd ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_livpre_out ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_liv_livpre ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_livpre_liv ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_livpre_plcp ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_livpre_livstore ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_livstore_livpre ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_gutd_feces ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_pl_liv ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_liv_pl ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_pl_tis ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_tis_pl ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_pl_tis2 ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_tis2_pl ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_plcp_tiscp ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_tiscp_plcp ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_tis_sweat ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_pl_urine ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalk_tis_abs ~ fixed(0.0225) # simulation assumption, Methods 'IIV was set to 15% CV'
    etalvpl ~ fixed(0.0225) # simulation assumption, Data S2 omega.csv carries V_pl

    # Residual error: proportional per matrix (Table 1). The fit weighted each
    # SD by 1/sqrt(N) for mean data; these are the per-individual values.
    propSd <- 0.698
    label("Proportional residual error, total plasma copper, literature data (fraction)") # Table 1 propRUVpl = 0.698
    propSd_comet <- 0.176
    label("Proportional residual error, total plasma copper, COMET study (fraction)") # Table 1 propRUVplCo = 0.176
    propSd_Ccp <- 0.441
    label("Proportional residual error, ceruloplasmin-bound plasma copper (fraction)") # Table 1 propRUVcp = 0.441
    propSd_cum_feces <- 0.711
    label("Proportional residual error, cumulative faecal excretion (fraction)") # Table 1 propRUVfe = 0.711
    propSd_cum_urine <- 0.611
    label("Proportional residual error, cumulative urinary excretion (fraction)") # Table 1 propRUVur = 0.611
    propSd_liver_total <- 0.203
    label("Proportional residual error, total liver copper (fraction)") # Table 1 propRUVliv = 0.203
  })

  model({
    # Genotype / treatment status. Healthy non-carriers have both activity
    # factors equal to 1; heterozygotes carry the estimated factors; WD patients
    # carry the gene-therapy-restored factor (0 if untreated), with biliary
    # activity tied to ceruloplasmin-loading activity by a power function.
    # The restored factor is the amount in trt_af, a constant state that starts
    # at 0 and is set by a bolus of the activity factor at vector injection
    # (Data S2 model.campsis: TRT state, A_TRT; Bolus(amount = 0.4,
    # compartment = 'TRT')).
    # trt_af is clamped at 0 because the stiff solver's trial values can dip
    # just below 0, where the non-integer power is NaN.
    healthy <- 1 - DIS_WILSON - ATP7B_HET
    af_wd_cp <- max(trt_af, 0)
    af_wd_bil <- af_wd_cp^e_af_bil
    af_cp <- healthy + ATP7B_HET * af_het_cp + DIS_WILSON * af_wd_cp
    af_bil <- healthy + ATP7B_HET * af_het_bil + DIS_WILSON * af_wd_bil

    # Individual rate constants (1/day)
    k_abs_liv <- exp(lk_abs_liv + etalk_abs_liv) * (1 + e_zinc_k_abs_liv * CONMED_ZINC)
    k_abs_gutd <- exp(lk_abs_gutd + etalk_abs_gutd)
    k_livpre_out <- exp(lk_livpre_out + etalk_livpre_out) * af_bil
    k_liv_livpre <- exp(lk_liv_livpre + etalk_liv_livpre) * (1 + e_wilson_k_liv_livpre * DIS_WILSON)
    k_livpre_liv <- exp(lk_livpre_liv + etalk_livpre_liv) * DIS_WILSON
    k_livpre_plcp <- exp(lk_livpre_plcp + etalk_livpre_plcp) * af_cp
    k_livpre_livstore <- exp(lk_livpre_livstore + etalk_livpre_livstore)
    k_livstore_livpre <- exp(lk_livstore_livpre + etalk_livstore_livpre)
    k_gutd_feces <- exp(lk_gutd_feces + etalk_gutd_feces) * (1 + e_comet_k_gutd_feces * STUDY_COMET)
    k_pl_liv <- exp(lk_pl_liv + etalk_pl_liv)
    k_liv_pl <- exp(lk_liv_pl + etalk_liv_pl)
    k_pl_tis <- exp(lk_pl_tis + etalk_pl_tis) * (1 + e_zinc_k_pl_tis * CONMED_ZINC)
    k_tis_pl <- exp(lk_tis_pl + etalk_tis_pl) * (1 + e_wilson_k_tis_pl * DIS_WILSON)
    k_pl_tis2 <- exp(lk_pl_tis2 + etalk_pl_tis2)
    k_tis2_pl <- exp(lk_tis2_pl + etalk_tis2_pl) * (1 + e_wilson_k_tis2_pl * DIS_WILSON)
    k_plcp_liv <- exp(lk_plcp_liv)
    k_tis_sweat <- exp(lk_tis_sweat + etalk_tis_sweat)
    k_pl_urine <- exp(lk_pl_urine + etalk_pl_urine) * (1 + e_chelator_k_pl_urine * CONMED_COPPER_CHELATOR)
    k_tis_abs <- exp(lk_tis_abs + etalk_tis_abs)
    k_plcp_tiscp <- exp(lk_plcp_tiscp + etalk_plcp_tiscp)
    k_tiscp_plcp <- exp(lk_tiscp_plcp + etalk_tiscp_plcp)
    vpl <- exp(lvpl + etalvpl)

    # Biliary output split: 25% reabsorbed in the proximal gut, 75% to distal gut
    k_liv_gutd <- k_livpre_out * f_bile_out
    k_liv_abs <- k_livpre_out * (1 - f_bile_out)

    d/dt(gut_abs) <- -k_abs_liv * gut_abs + k_liv_abs * liver_precursor - k_abs_gutd * gut_abs + k_tis_abs * tissue_ncc
    d/dt(gut_distal) <- k_abs_gutd * gut_abs + k_liv_gutd * liver_precursor - k_gutd_feces * gut_distal
    d/dt(liver) <- k_abs_liv * gut_abs + k_pl_liv * plasma_ncc - k_liv_pl * liver - k_liv_livpre * liver + k_plcp_liv * plasma_cp + k_livpre_liv * liver_precursor
    d/dt(liver_storage) <- k_livpre_livstore * liver_precursor - k_livstore_livpre * liver_storage
    d/dt(liver_precursor) <- k_liv_livpre * liver + k_livstore_livpre * liver_storage - k_livpre_livstore * liver_precursor - k_livpre_plcp * liver_precursor - k_liv_abs * liver_precursor - k_liv_gutd * liver_precursor - k_livpre_liv * liver_precursor
    d/dt(plasma_ncc) <- k_liv_pl * liver - k_pl_liv * plasma_ncc - k_pl_tis * plasma_ncc + k_tis_pl * tissue_ncc - k_pl_urine * plasma_ncc - k_pl_tis2 * plasma_ncc + k_tis2_pl * tissue_deep
    d/dt(plasma_cp) <- k_livpre_plcp * liver_precursor - k_plcp_tiscp * plasma_cp + k_tiscp_plcp * tissue_cp - k_plcp_liv * plasma_cp
    d/dt(tissue_ncc) <- k_pl_tis * plasma_ncc - k_tis_pl * tissue_ncc - k_tis_sweat * tissue_ncc - k_tis_abs * tissue_ncc
    d/dt(tissue_deep) <- k_pl_tis2 * plasma_ncc - k_tis2_pl * tissue_deep
    d/dt(tissue_cp) <- k_plcp_tiscp * plasma_cp - k_tiscp_plcp * tissue_cp
    d/dt(gut_transit) <- k_gutd_feces * gut_distal - k_gutd_feces * gut_transit
    d/dt(feces) <- k_gutd_feces * gut_transit
    d/dt(urine) <- k_pl_urine * plasma_ncc
    d/dt(sweat) <- k_tis_sweat * tissue_ncc
    d/dt(trt_af) <- 0

    # Fitted outputs (amounts in dose units; concentrations in mg/L)
    Cc <- (plasma_ncc + plasma_cp) / vpl
    Ccp <- plasma_cp / vpl
    cum_feces <- feces
    cum_urine <- urine
    liver_total <- liver + liver_storage + liver_precursor

    # Cold-copper biomarkers (Data S2 model.campsis; molar mass 63.5 g/mol,
    # liver 1500 g wet weight with 75% water)
    Cncc <- plasma_ncc / vpl
    cu_total_umol_l <- 1000 * Cc / 63.5
    ncc_umol_l <- 1000 * Cncc / 63.5
    rec <- 100 * Cncc / Cc
    liver_ug_g <- 1000 * liver_total / (1500 * (1 - 0.75))
    urine_rate <- k_pl_urine * plasma_ncc

    # COMET plasma data carry their own residual error (NONMEM DVID 1 split by isCOMET)
    sd_cc <- propSd * (1 - STUDY_COMET) + propSd_comet * STUDY_COMET
    Cc ~ prop(sd_cc)
    Ccp ~ prop(propSd_Ccp)
    cum_feces ~ prop(propSd_cum_feces)
    cum_urine ~ prop(propSd_cum_urine)
    liver_total ~ prop(propSd_liver_total)
  })
}

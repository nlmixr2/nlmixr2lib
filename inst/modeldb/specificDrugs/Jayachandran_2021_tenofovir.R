Jayachandran_2021_tenofovir <- function() {
  description <- paste(
    "Seven-compartment multicompartment population PK model for tenofovir (TFV)",
    "and its intracellular anabolite tenofovir-diphosphate (TFVdp) after oral",
    "tenofovir disoproxil fumarate (TDF) or rectally applied 1% TFV gel, from the",
    "RMP-02/MTN-006 phase 1 pre-exposure-prophylaxis trial. Plasma TFV is a",
    "two-compartment model with first-order oral or rectal absorption; the tissue",
    "and cellular matrices are added as biophase effect compartments driven by the",
    "plasma TFV concentration, each with its own first-order equilibration rate",
    "and a plasma-to-matrix concentration ratio: TFV in rectal-tissue homogenate,",
    "and TFVdp in rectal-tissue homogenate, in rectal mucosal mononuclear cells",
    "(MMCs) and in peripheral blood mononuclear cells (PBMCs). Both the",
    "equilibration rate and the ratio are stratified by administration route",
    "(oral vs rectal). The companion ex vivo viral-dynamics PK/PD model that this",
    "PK model feeds is Jayachandran_2021_tenofovir_p24.",
    sep = " "
  )
  reference <- paste(
    "Jayachandran P, Garcia-Cremades M, Vucicevic K, Bumpus NN, Anton P,",
    "Hendrix C, Savic R. A Mechanistic In Vivo/Ex Vivo",
    "Pharmacokinetic-Pharmacodynamic Model of Tenofovir for HIV Prevention.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(3):179-187.",
    "doi:10.1002/psp4.12583. Structural parameters from Table 1; plasma",
    "two-compartment and biophase effect-compartment equations from Methods",
    "(Eqs 1-4 and the plasma disposition); route stratification and fixed uptake",
    "rates from the Supplementary Model Code 1 NONMEM control stream.",
    sep = " "
  )
  vignette <- "Jayachandran_2021_tenofovir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    ROUTE_ORAL = list(
      description = "Administration-route indicator selecting the absorption and biophase parameter set",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (rectal 1% TFV gel)",
      notes = paste(
        "1 = dose administered orally as tenofovir disoproxil fumarate (TDF);",
        "0 = dose administered rectally as 1% TFV gel. Carried per dose record.",
        "In the source Supplementary Model Code 1 the same information is the",
        "integer column ROUTE, branched as IF(ROUTE.EQ.2) for the rectal arm, so",
        "ROUTE_ORAL = as.integer(ROUTE == 1). The indicator selects, structurally",
        "rather than as a covariate coefficient: the oral vs rectal first-order",
        "absorption rate constant, oral (F = 1) vs rectal (F = 0.102)",
        "bioavailability, the residual-error structure (proportional for oral,",
        "combined additive-plus-proportional for rectal), and for every biophase",
        "matrix the route-specific equilibration rate and plasma-to-matrix ratio.",
        sep = " "
      ),
      source_name = "ROUTE (1 = oral, 2 = rectal)"
    )
  )

  # Biophase (effect-compartment) states are carried directly as CONCENTRATIONS,
  # not amounts (Methods Eqs 1-4 are written in dC/dt form), so no volume divides
  # them at the observation step -- the same encoding as the intracellular pools
  # in Yu_2026_tenofovir.R. `pbmc_tfvdp` uses the registered pbmc + _tfvdp
  # canonical; the rectal-tissue and MMC biophase states have no canonical and
  # are declared paper-specific.
  paper_specific_compartments <- c("rt_tfv", "rt_tfvdp", "mmc_tfvdp")

  compartmentData <- list(
    depot = list(analyte = "tenofovir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tenofovir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tenofovir", units = "mg", specimen = "plasma", verified = TRUE),
    rt_tfv = list(analyte = "tenofovir", units = "ng/mg", specimen = "tissue", verified = TRUE),
    rt_tfvdp = list(analyte = "tenofovir diphosphate", units = "fmol/mg", specimen = "tissue", verified = TRUE),
    mmc_tfvdp = list(
      analyte = "tenofovir diphosphate",
      units = "fmol/million cells",
      specimen = "blood cell",
      verified = TRUE
    ),
    pbmc_tfvdp = list(
      analyte = "tenofovir diphosphate",
      units = "fmol/million cells",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "22-66 years",
    sex_female_pct = 22.2,
    disease_state = "HIV-1 seronegative healthy adults (HIV pre-exposure-prophylaxis study)",
    dose_range = paste(
      "Single oral 300 mg TDF (= 136 mg TFV); single and multiple (7 daily) rectal",
      "1% TFV gel (44 mg TFV per dose)",
      sep = " "
    ),
    regions = "USA (Los Angeles, CA and Pittsburgh, PA)",
    notes = paste(
      "RMP-02/MTN-006 (NCT00984971): a two-site phase 1 partially blinded,",
      "placebo-controlled safety, acceptability and PK trial. Eighteen subjects",
      "(14 men, 4 women), aged 22-66 years. Arm 1 single oral TDF; arm 2 single",
      "rectal TFV or placebo gel (2:1); arm 3 seven consecutive daily rectal",
      "doses. PK sampled in five matrices: TFV (plasma, rectal-tissue homogenate)",
      "and TFVdp (PBMCs, rectal mucosal mononuclear cells, rectal-tissue",
      "homogenate). Tissue and cellular matrices were sparse and heavily censored",
      "(70-94% below the limit of quantification), so IIV was retained only on",
      "oral Ka and on central volume (Table 1); rectal PBMC uptake was not",
      "estimable because all rectal-arm PBMC concentrations were below the limit",
      "of quantification.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma tenofovir -- two-compartment disposition with first-order
    # oral or rectal absorption. Jayachandran 2021 Table 1 (Plasma, TFV).
    # Concentrations are in ng/mL (= mcg/L); doses and amounts in mg, so
    # Cc = 1000 * central / vc converts mg/L to ng/mL (control-stream
    # S2 = V2/1000).
    # ------------------------------------------------------------------
    lka_oral <- log(1.78); label("First-order oral absorption rate constant (1/h)")   # Table 1 'k a oral' = 1.78 /h (RSE 30%)
    lka_rectal <- log(22.2); label("First-order rectal absorption rate constant (1/h)") # Table 1 'k a rectal' = 22.2 /h (RSE 45%)
    lcl <- log(52.1); label("Apparent clearance CL/F of plasma TFV (L/h)")            # Table 1 'CL/F' = 52.1 L/h (RSE 8%)
    lvc <- log(679); label("Apparent central volume of distribution V2/F (L)")        # Table 1 'V 2 /F' = 679 L (RSE 11%)
    lq <- log(7.90); label("Apparent intercompartmental clearance Q/F (L/h)")         # Table 1 'Q/F' = 7.90 L/h (RSE 14%)
    lvp <- log(399); label("Apparent peripheral volume of distribution V3/F (L)")     # Table 1 'V 3 /F' = 399 L (RSE 12%)
    lfrectal <- log(0.102); label("Relative bioavailability of the rectal route (oral = 1)") # Table 1 'F rectal' = 0.102 (RSE 17%)

    # ------------------------------------------------------------------
    # Biophase effect compartments (Methods Eqs 1-4):
    #   dC_matrix/dt = ke0_matrix * (ppc_matrix * Cp - C_matrix)
    # ke0_matrix is the route-specific equilibration rate (source k_P-matrix)
    # and ppc_matrix the route-specific plasma-to-matrix concentration ratio
    # (source R_P-matrix). Cp is the plasma TFV concentration in ng/mL.
    # Both terms are stratified oral vs rectal (Table 1). Rectal-tissue
    # uptake rates were FIXED in the source (Table 1 footnote 'a').
    # ------------------------------------------------------------------

    # Rectal tissue, TFV
    lke0_rt_oral <- fixed(log(0.01)); label("Rectal-tissue TFV equilibration rate, oral (1/h)")     # Table 1 'k o_P-RT' = 0.01 (fixed)
    lke0_rt_rectal <- fixed(log(0.0001)); label("Rectal-tissue TFV equilibration rate, rectal (1/h)") # Table 1 'k r_P-RT' = 0.0001 (fixed)
    lppc_rt_oral <- log(2.05); label("Rectal-tissue:plasma TFV ratio, oral ((ng/mg)/(mcg/L))")      # Table 1 'R o_P-RT' = 2.05 (RSE 19%)
    lppc_rt_rectal <- log(171); label("Rectal-tissue:plasma TFV ratio, rectal ((ng/mg)/(mcg/L))")   # Table 1 'R r_P-RT' = 171 (RSE 95%)

    # Rectal tissue, TFVdp
    lke0_rtm_oral <- fixed(log(0.0001)); label("Rectal-tissue TFVdp equilibration rate, oral (1/h)")   # Table 1 'k o_P-RTm' = 0.0001 (fixed)
    lke0_rtm_rectal <- fixed(log(0.01)); label("Rectal-tissue TFVdp equilibration rate, rectal (1/h)") # Table 1 'k r_P-RTm' = 0.01 (fixed)
    lppc_rtm_oral <- log(530); label("Rectal-tissue:plasma TFVdp ratio, oral ((fmol/mg)/(mcg/L))")     # Table 1 'R o_P-RTm' = 530 (RSE 51%)
    lppc_rtm_rectal <- log(540); label("Rectal-tissue:plasma TFVdp ratio, rectal ((fmol/mg)/(mcg/L))") # Table 1 'R r_P-RTm' = 540 (RSE 2%)

    # MMC (rectal mucosal mononuclear cells), TFVdp
    lke0_mmc_oral <- log(0.0102); label("MMC TFVdp equilibration rate, oral (1/h)")    # Table 1 'k o_P-MMC' = 0.0102 (RSE 13%)
    lke0_mmc_rectal <- log(0.0606); label("MMC TFVdp equilibration rate, rectal (1/h)") # Table 1 'k r_P-MMC' = 0.0606 (RSE 71%)
    lppc_mmc_oral <- log(11.4); label("MMC:plasma TFVdp ratio, oral ((fmol/million cells)/(mcg/L))")    # Table 1 'R o_P-MMC' = 11.4 (RSE 16%)
    lppc_mmc_rectal <- log(1530); label("MMC:plasma TFVdp ratio, rectal ((fmol/million cells)/(mcg/L))") # Table 1 'R r_P-MMC' = 1530 (RSE 24%)

    # PBMC (peripheral blood mononuclear cells), TFVdp -- oral only; the
    # rectal rate and ratio were not estimable (all rectal PBMC data BLQ).
    lke0_pbmc_oral <- log(0.0320); label("PBMC TFVdp equilibration rate, oral (1/h)")  # Table 1 'k o_P-PBMC' = 0.0320 (RSE 49%)
    lppc_pbmc_oral <- log(0.169); label("PBMC:plasma TFVdp ratio, oral ((fmol/million cells)/(mcg/L))") # Table 1 'R o_P-PBMC' = 0.169 (RSE 61%)

    # ------------------------------------------------------------------
    # Between-subject variability. Table 1 reports IIV only on oral Ka and
    # on V2/F; the source $OMEGA gives the variances directly and the
    # printed %CV is sqrt(variance) (0.262 -> 51.2%; 0.045 -> 21.2%).
    # ------------------------------------------------------------------
    etalka_oral ~ 0.262  # variance 0.262; sqrt(0.262) = 0.512 = Table 1 'IIV k a oral 51.2% (RSE 33%)'
    etalvc ~ 0.045       # variance 0.045; sqrt(0.045) = 0.212 = Table 1 'IIV V 2 /F 21.2% (RSE 31%)'

    # ------------------------------------------------------------------
    # Residual error. Plasma TFV is proportional after oral dosing and
    # combined additive + proportional after rectal dosing; each tissue /
    # cellular matrix is proportional (Table 1). %CV entries are the
    # proportional SD as a fraction.
    # ------------------------------------------------------------------
    propSd_oral <- 0.515; label("Proportional residual error, plasma TFV, oral (fraction)")   # Table 1 'Proportional error (oral)' 51.5% CV (RSE 8%)
    propSd_rectal <- 0.520; label("Proportional residual error, plasma TFV, rectal (fraction)") # Table 1 'Proportional error (rectal)' 52.0% CV (RSE 17%)
    addSd_rectal <- 3.71; label("Additive residual error, plasma TFV, rectal (ng/mL)")         # Table 1 'Additive error (rectal)' 3.71 mcg/L (RSE 18%)
    propSd_Crt_tfv <- 1.00; label("Proportional residual error, rectal-tissue TFV (fraction)")   # Table 1 'Proportional error' RT TFV 100% CV (RSE 14%)
    propSd_Crt_tfvdp <- 2.38; label("Proportional residual error, rectal-tissue TFVdp (fraction)") # Table 1 'Proportional error' RT TFVdp 238% CV (RSE 7%)
    propSd_Cmmc_tfvdp <- 1.33; label("Proportional residual error, MMC TFVdp (fraction)")        # Table 1 'Proportional error' MMC 133% CV (RSE 12%)
    propSd_Cpbmc_tfvdp <- 4.61; label("Proportional residual error, PBMC TFVdp (fraction)")      # Table 1 'Proportional error' PBMC 461% CV (RSE 22%)
  })

  model({
    # 1. Route-selected absorption and bioavailability.
    ka <- exp(lka_oral * ROUTE_ORAL + lka_rectal * (1 - ROUTE_ORAL) + etalka_oral * ROUTE_ORAL)
    cl <- exp(lcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    # 2. Micro-constants for the plasma two-compartment model.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Route-selected biophase equilibration rates (ke0) and
    #    plasma-to-matrix ratios (ppc). For PBMC only the oral parameters
    #    were estimated; a rectal dose drives no PBMC uptake (rectal PBMC
    #    data were all BLQ), so the ratio is gated to zero for rectal.
    ke0_rt <- exp(lke0_rt_oral * ROUTE_ORAL + lke0_rt_rectal * (1 - ROUTE_ORAL))
    ppc_rt <- exp(lppc_rt_oral * ROUTE_ORAL + lppc_rt_rectal * (1 - ROUTE_ORAL))
    ke0_rtm <- exp(lke0_rtm_oral * ROUTE_ORAL + lke0_rtm_rectal * (1 - ROUTE_ORAL))
    ppc_rtm <- exp(lppc_rtm_oral * ROUTE_ORAL + lppc_rtm_rectal * (1 - ROUTE_ORAL))
    ke0_mmc <- exp(lke0_mmc_oral * ROUTE_ORAL + lke0_mmc_rectal * (1 - ROUTE_ORAL))
    ppc_mmc <- exp(lppc_mmc_oral * ROUTE_ORAL + lppc_mmc_rectal * (1 - ROUTE_ORAL))
    ke0_pbmc <- exp(lke0_pbmc_oral)
    ppc_pbmc <- exp(lppc_pbmc_oral) * ROUTE_ORAL

    # 4. Plasma TFV concentration (ng/mL) drives every biophase compartment.
    cp <- 1000 * central / vc

    # 5. ODE system. depot / central / peripheral1 are amounts (mg); the
    #    four biophase states are concentrations in their matrix units.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(rt_tfv) <- ke0_rt * (ppc_rt * cp - rt_tfv)
    d/dt(rt_tfvdp) <- ke0_rtm * (ppc_rtm * cp - rt_tfvdp)
    d/dt(mmc_tfvdp) <- ke0_mmc * (ppc_mmc * cp - mmc_tfvdp)
    d/dt(pbmc_tfvdp) <- ke0_pbmc * (ppc_pbmc * cp - pbmc_tfvdp)

    # 6. Rectal bioavailability (oral F = 1).
    f(depot) <- 1 * ROUTE_ORAL + exp(lfrectal) * (1 - ROUTE_ORAL)

    # 7. Observations. Plasma TFV in ng/mL; biophase states are already
    #    concentrations in their matrix units.
    Cc <- cp
    Crt_tfv <- rt_tfv
    Crt_tfvdp <- rt_tfvdp
    Cmmc_tfvdp <- mmc_tfvdp
    Cpbmc_tfvdp <- pbmc_tfvdp

    # Plasma residual error switches by route: proportional for oral,
    # combined additive + proportional for rectal.
    ruvProp <- propSd_oral * ROUTE_ORAL + propSd_rectal * (1 - ROUTE_ORAL)
    ruvAdd <- addSd_rectal * (1 - ROUTE_ORAL)
    Cc ~ add(ruvAdd) + prop(ruvProp)
    Crt_tfv ~ prop(propSd_Crt_tfv)
    Crt_tfvdp ~ prop(propSd_Crt_tfvdp)
    Cmmc_tfvdp ~ prop(propSd_Cmmc_tfvdp)
    Cpbmc_tfvdp ~ prop(propSd_Cpbmc_tfvdp)
  })
}

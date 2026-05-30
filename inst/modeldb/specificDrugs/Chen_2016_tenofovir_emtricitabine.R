Chen_2016_tenofovir_emtricitabine <- function() {
  description <- "Linked population PKPD model for daily oral co-administered tenofovir (TFV, given as the prodrug TDF 300 mg = TFV 136 mg) and emtricitabine (FTC 200 mg) in HIV-positive and HIV-negative adults (Chen 2016 Cell-PrEP study). Each parent drug is described by a two-compartment first-order-absorption plasma popPK model. Each parent feeds a hybrid first-order-formation + saturation link into its intracellular triphosphate anabolite in peripheral blood mononuclear cells (TFV-DP, FTC-TP), modelled with a two-compartment 'recycle' elimination structure where a fraction R of the eliminated drug re-enters the central intracellular compartment. Each anabolite inhibits the zero-order production rate of two endogenous deoxynucleoside triphosphates via an Emax indirect-response model with Kout fixed to 1/day and Emax fixed to 1: TFV-DP inhibits dATP and dGTP (deoxypurines); FTC-TP inhibits dCTP and TTP (deoxypyrimidines). The dGTP effect waned over time and is described by an additional 1/(1+t^gamma) time factor. Sex is a covariate on FTC plasma Vc/F and HIV-infection status is a covariate on FTC-TP Kf. Intended for simulating analog:dNTP molar ratios (TFV-DP:dATP, FTC-TP:dCTP) for various dosing strategies, e.g., the IPERGAY on-demand PrEP regimen."
  reference <- "Chen X, Seifert SM, Castillo-Mancilla JR, Bushman LR, Zheng J-H, Kiser JJ, MaWhinney S, Anderson PL. Model Linking Plasma and Intracellular Tenofovir/Emtricitabine with Deoxynucleoside Triphosphates. PLoS ONE. 2016;11(11):e0165505. doi:10.1371/journal.pone.0165505"
  vignette <- "Chen_2016_tenofovir_emtricitabine"
  paper_specific_compartments <- c("tfvdp", "ftctp", "datp", "dgtp", "dctp", "ttp")

  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) per the canonical SEXF register; matches Chen 2016 'sex = 0 (female)' once values are inverted",
      notes              = "Chen 2016 reports the FTC plasma Vc/F covariate effect using a male indicator: tvVc/F = 99.4 + 24.3 * SEX with male = 1, female = 0 (Table 1B). The canonical column in inst/references/covariate-columns.md is SEXF (1 = female), so values invert: SEXF = 1 - source_SEX. The model applies the +24.3 L additive shift via (1 - SEXF) so that males (SEXF = 0) get Vc/F = 123.7 L and females (SEXF = 1) get Vc/F = 99.4 L, preserving the paper's parameterisation. No effect on TFV plasma or any intracellular compartment.",
      source_name        = "SEX (male = 1, female = 0)"
    ),
    HIV_POS = list(
      description        = "HIV-1 serostatus at study entry (1 = HIV-positive, 0 = HIV-negative)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative)",
      notes              = "Chen 2016 reports the FTC-TP intracellular Kf covariate effect as tvKf = 41.6 + 31.3 * HIV with HIV-positive = 1 and HIV-negative = 0 (Table 1B, FTC-TP IC section). Source column orientation matches the canonical HIV_POS (1 = positive). The +31.3 1/day additive shift gives HIV-positive subjects a 75.2% higher FTC-TP formation rate (the paper's reported effect). HIV-positive subjects in this cohort received daily co-administered efavirenz 600 mg in addition to TDF/FTC; the model does not separate the HIV-status indicator from possible efavirenz-mediated effects (see vignette Assumptions and deviations). No effect on TFV plasma, FTC plasma, TFV-DP IC, or any of the four dNTP responses.",
      source_name        = "HIV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "20-52 years (median 31)",
    age_median     = "31 years",
    weight_range   = "56.5-127 kg (median 81.1)",
    weight_median  = "81.1 kg",
    bmi_range      = "19.9-37.7 kg/m^2 (median 26.6)",
    bmi_median     = "26.6 kg/m^2",
    egfr_range     = "66.0-131 mL/min/1.73m^2 (median 93.3)",
    egfr_median    = "93.3 mL/min/1.73m^2",
    sex_female_pct = 32.5,
    race_ethnicity = c(White = 47.5, Black_or_African_American = 40.0, Hispanic = 12.5),
    disease_state  = "21 HIV-negative healthy adults and 19 HIV-positive adults (the Cell-PrEP study; protocol NCT01040091). HIV-negative subjects received daily oral co-formulated TDF 300 mg / FTC 200 mg for 30 days followed by a washout period (days 35, 45, 60). HIV-positive subjects received daily TDF/FTC + efavirenz 600 mg for 60 days. Baseline dNTP samples taken prior to first dose.",
    dose_range     = "TDF 300 mg (= TFV 136 mg) + FTC 200 mg by mouth once daily for 30 days (HIV-negative) or 60 days (HIV-positive). HIV-positive group co-administered efavirenz 600 mg PO QD.",
    regions        = "University of Colorado Anschutz Medical Campus, USA (single-center; enrollment 2010-2013)",
    notes          = "Demographics from Chen 2016 Results 'Study demographics' (page 6 of the PLOS ONE PDF). Sampling on day 1 and day 30 at 1, 2, 4, 8, 24 hours post-dose; on days 3, 7, 20 pre-dose and 2, 8 hours post-dose; single samples on days 35, 45, 60. dNTP pool components (dATP, dGTP, dCTP, TTP) measured in PBMC at baseline, 1, 2, 4, 8, 24 hours post-dose on day 1, and 8 hours post-dose on subsequent visits. 34 of 40 subjects completed all visits. Below-limit-of-quantitation (BLQ) samples were treated as missing during NONMEM estimation."
  )

  ini({
    # ============================================================
    # Tenofovir (TFV) plasma 2-compartment popPK. Chen 2016 Table 1A
    # (TFV plasma section). Apparent parameters CL/F, Vc/F, Vp/F,
    # Q/F: bioavailability is implicit in F-adjusted CL and V. The
    # daily oral dose is the TFV-equivalent of TDF: 300 mg TDF
    # delivers 136 mg TFV (Methods 'Study design').
    # ============================================================
    lka  <- log(80.1);  label("TFV first-order oral absorption rate (1/day)")        # Chen 2016 Table 1A (Ka_TFV = 80.1 /day)
    lcl  <- log(1410);  label("TFV apparent oral clearance CL/F (L/day)")             # Chen 2016 Table 1A (CL/F_TFV = 1410 L/day)
    lvc  <- log(390);   label("TFV apparent central volume Vc/F (L)")                 # Chen 2016 Table 1A (Vc/F_TFV = 390 L)
    lq   <- log(5390);  label("TFV apparent intercompartmental clearance Q/F (L/day)") # Chen 2016 Table 1A (Q/F_TFV = 5390 L/day)
    lvp  <- log(877);   label("TFV apparent peripheral volume Vp/F (L)")              # Chen 2016 Table 1A (Vp/F_TFV = 877 L)

    # ============================================================
    # Emtricitabine (FTC) plasma 2-compartment popPK. Chen 2016
    # Table 1B (FTC plasma section). Daily oral dose: 200 mg FTC.
    # ============================================================
    lka_ftc  <- log(55.7); label("FTC first-order oral absorption rate (1/day)")       # Chen 2016 Table 1B (Ka_FTC = 55.7 /day)
    lcl_ftc  <- log(482);  label("FTC apparent oral clearance CL/F (L/day)")           # Chen 2016 Table 1B (CL/F_FTC = 482 L/day)
    lvc_ftc  <- log(99.4); label("FTC apparent central volume Vc/F at female reference (L)")  # Chen 2016 Table 1B (Vc/F_FTC = 99.4 L; female reference)
    lq_ftc   <- log(141);  label("FTC apparent intercompartmental clearance Q/F (L/day)") # Chen 2016 Table 1B (Q/F_FTC = 141 L/day)
    lvp_ftc  <- log(166);  label("FTC apparent peripheral volume Vp/F (L)")             # Chen 2016 Table 1B (Vp/F_FTC = 166 L)

    # Sex effect on FTC Vc/F. Paper: tvVc/F = 99.4 + 24.3 * SEX with
    # male = 1, female = 0 (Table 1B "Sex on Vc/F (linear)" row).
    # Canonical SEXF = 1 - source_SEX, so the model applies the
    # +24.3 L additive shift via (1 - SEXF).
    e_sexm_vc_ftc <- 24.3; label("Additive shift on FTC Vc/F for male sex (L); applied as (1 - SEXF)") # Chen 2016 Table 1B (Sex on Vc/F = 24.3 L)

    # ============================================================
    # Intracellular TFV-DP. Chen 2016 Table 1A (TFV-DP IC section)
    # and Methods 'Population pharmacokinetics modeling of TFV-DP/
    # FTC-TP'. Hybrid first-order formation + saturation link from
    # plasma TFV: formation rate = Kf * Cp_TFV * SC50/(SC50 + sat),
    # where the virtual saturation compartment (effect_tfvdp) is
    # driven by Kf * Cp_TFV and eliminates at CL_TFV/Vc_TFV (matching
    # the plasma central compartment). Elimination of TFV-DP follows
    # a two-compartment "recycle" structure: K46 = K64 = R * Kel
    # routes between IC central and recycle (peripheral1_tfvdp); the
    # IC central elimination is K40 = (1 - R) * Kel.
    # Concentrations: fmol per 10^6 PBMC. Saturation compartment
    # SC50 has plasma-concentration units (ng/mL).
    # ============================================================
    lkf_tfvdp   <- log(1.4);   label("TFV-DP first-order formation rate constant Kf from plasma TFV (1/day)")   # Chen 2016 Table 1A (Kf_TFVDP = 1.4 /day)
    lsc50_tfvdp <- log(6.55);  label("TFV-DP saturation half-effect concentration SC50 (ng/mL of plasma TFV)") # Chen 2016 Table 1A (SC50_TFV = 6.55)
    lkel_tfvdp  <- log(0.228); label("TFV-DP total elimination rate constant Kel from IC central (1/day)")    # Chen 2016 Table 1A (Kel_TFVDP = 0.228 /day)
    r_tfvdp     <- 0.0582;     label("TFV-DP recycle ratio R (unitless, range 0-1); K46 = K64 = R * Kel; K40 = (1 - R) * Kel") # Chen 2016 Table 1A (R_TFVDP = 5.82%)

    # ============================================================
    # Intracellular FTC-TP. Chen 2016 Table 1B (FTC-TP IC section).
    # Same hybrid first-order + saturation structure as TFV-DP.
    # ============================================================
    lkf_ftctp   <- log(41.6);  label("FTC-TP first-order formation rate constant Kf from plasma FTC at HIV-negative reference (1/day)") # Chen 2016 Table 1B (Kf_FTCTP = 41.6 /day; HIV-negative reference)
    lsc50_ftctp <- log(3320);  label("FTC-TP saturation half-effect concentration SC50 (ng/mL of plasma FTC)") # Chen 2016 Table 1B (SC50_FTC = 3320)
    lkel_ftctp  <- log(1.6);   label("FTC-TP total elimination rate constant Kel from IC central (1/day)")    # Chen 2016 Table 1B (Kel_FTCTP = 1.6 /day)
    r_ftctp     <- 0.160;      label("FTC-TP recycle ratio R (unitless, range 0-1)")                          # Chen 2016 Table 1B (R_FTCTP = 16.0%)

    # HIV serostatus effect on FTC-TP Kf. Paper: tvKf = 41.6 + 31.3
    # * HIV with HIV-positive = 1 (Table 1B "HIV on Kf (Linear)").
    e_hiv_pos_kf_ftctp <- 31.3; label("Additive shift on FTC-TP Kf for HIV-positive (1/day); applied as HIV_POS") # Chen 2016 Table 1B (HIV on Kf = 31.3 /day)

    # ============================================================
    # Endogenous deoxypurine indirect-response models driven by
    # TFV-DP. Chen 2016 Table 1A (dATP, dGTP). Indirect-response
    # Emax structure with Emax fixed to 1 and Kout fixed to 1/day
    # (Methods 'Population pharmacodynamics modeling'). K0in is
    # interconvertible with the baseline level R0 via K0in = R0 *
    # Kout. The paper estimated tvEC50 = exp(theta) to ensure
    # positive estimates; the values below are the back-transformed
    # population-mean EC50 from Table 1A and are stored on the log
    # scale for consistency with the other paper-named parameters.
    # ============================================================
    lr0_datp   <- log(155);  label("Baseline dATP level R0 (fmol per 10^6 PBMC)")                        # Chen 2016 Table 1A (R0_dATP = 155)
    lec50_datp <- log(1020); label("EC50 of TFV-DP for 50% inhibition of dATP production (fmol/10^6 PBMC)") # Chen 2016 Table 1A (EC50_TFVDP_dATP = 1020)

    lr0_dgtp   <- log(245);  label("Baseline dGTP level R0 (fmol per 10^6 PBMC)")                        # Chen 2016 Table 1A (R0_dGTP = 245)
    lec50_dgtp <- log(54.6); label("EC50 of TFV-DP for 50% inhibition of dGTP production (fmol/10^6 PBMC)") # Chen 2016 Table 1A (EC50_TFVDP_dGTP = 54.6)
    gamma_dgtp <- 0.928;     label("Power exponent of time (day^gamma) in the waning factor 1/(1 + t^gamma) applied to the dGTP inhibition term (unitless)") # Chen 2016 Table 1A (gamma_dGTP = 0.928)

    # ============================================================
    # Endogenous deoxypyrimidine indirect-response models driven by
    # FTC-TP. Chen 2016 Table 1B (dCTP, TTP). Same indirect-response
    # Emax structure with Emax fixed to 1 and Kout fixed to 1/day.
    # ============================================================
    lr0_dctp   <- log(771);   label("Baseline dCTP level R0 (fmol per 10^6 PBMC)")                          # Chen 2016 Table 1B (R0_dCTP = 771)
    lec50_dctp <- log(44400); label("EC50 of FTC-TP for 50% inhibition of dCTP production (fmol/10^6 PBMC)") # Chen 2016 Table 1B (EC50_FTCTP_dCTP = 44400)

    lr0_ttp    <- log(335);   label("Baseline TTP level R0 (fmol per 10^6 PBMC)")                           # Chen 2016 Table 1B (R0_TTP = 335)
    lec50_ttp  <- log(18800); label("EC50 of FTC-TP for 50% inhibition of TTP production (fmol/10^6 PBMC)") # Chen 2016 Table 1B (EC50_FTCTP_TTP = 18800)

    # ============================================================
    # Inter-individual variability (omega^2 on the log / variance
    # scale; exponential random-effect model on the typical value).
    # Chen 2016 Tables 1A and 1B report estimates directly as
    # variances (the "Estimate" column under "Interindividual
    # Variability"; the %CV column is sqrt(omega^2) * 100). Only
    # parameters for which the paper reports an IIV estimate are
    # included; the rest were not estimated.
    # ============================================================
    etalvc ~ 0.288   # Chen 2016 Table 1A (omega^2 Vc/F_TFV = 0.288; %CV 53.7)
    etalcl ~ 0.117   # Chen 2016 Table 1A (omega^2 CL/F_TFV = 0.117; %CV 34.2)
    etalq  ~ 0.0693  # Chen 2016 Table 1A (omega^2 Q/F_TFV  = 0.0693; %CV 26.3)

    etalvc_ftc ~ 0.0319  # Chen 2016 Table 1B (omega^2 Vc/F_FTC = 0.0319; %CV 17.9)
    etalcl_ftc ~ 0.0942  # Chen 2016 Table 1B (omega^2 CL/F_FTC = 0.0942; %CV 30.7)
    etalvp_ftc ~ 0.0335  # Chen 2016 Table 1B (omega^2 Vp/F_FTC = 0.0335; %CV 18.3)

    etalkf_tfvdp  ~ 0.238   # Chen 2016 Table 1A (omega^2 Kf_TFVDP = 0.238; %CV 48.8)
    etalkel_tfvdp ~ 0.316   # Chen 2016 Table 1A (omega^2 Kel_TFVDP = 0.316; %CV 56.2)

    etalkf_ftctp  ~ 0.0358  # Chen 2016 Table 1B (omega^2 Kf_FTCTP = 0.0358; %CV 18.9)
    etalkel_ftctp ~ 0.0561  # Chen 2016 Table 1B (omega^2 Kel_FTCTP = 0.0561; %CV 23.7)

    etalr0_datp   ~ 0.220   # Chen 2016 Table 1A (omega^2 R0_dATP = 0.220; %CV 46.9)
    etalec50_datp ~ 1.70    # Chen 2016 Table 1A (omega^2 EC50_dATP = 1.70; %CV 130; bootstrapping unstable, see Errata)

    etalr0_dgtp   ~ 0.203   # Chen 2016 Table 1A (omega^2 R0_dGTP = 0.203; %CV 45.1) -- IIV on EC50_dGTP not estimable per the paper

    etalr0_dctp   ~ 0.195   # Chen 2016 Table 1B (omega^2 R0_dCTP = 0.195; %CV 44.2)
    etalec50_dctp ~ 0.68    # Chen 2016 Table 1B (omega^2 EC50_dCTP = 0.68; %CV 82.5; bootstrapping unstable, see Errata)

    etalr0_ttp    ~ 0.383   # Chen 2016 Table 1B (omega^2 R0_TTP = 0.383; %CV 61.9)
    etalec50_ttp  ~ 1.02    # Chen 2016 Table 1B (omega^2 EC50_TTP = 1.02; %CV 101)

    # ============================================================
    # Residual error. Chen 2016 reports the variance (the "Estimate"
    # column for sigma) and indicates the error model in the row
    # label (exponential = log-normal lnorm; proportional = linear
    # proportional). SDs below are sqrt(variance).
    # ============================================================
    expSd        <- 0.273; label("TFV plasma log-scale residual SD (unitless; sqrt of variance 0.0745)") # Chen 2016 Table 1A (sigma_TFVplasma = 0.0745)
    expSd_ftc    <- 0.349; label("FTC plasma log-scale residual SD (unitless; sqrt of variance 0.122)")  # Chen 2016 Table 1B (sigma_FTCplasma = 0.122)
    propSd_tfvdp <- 0.339; label("TFV-DP IC proportional residual SD (fraction; sqrt of variance 0.115)") # Chen 2016 Table 1A (sigma_TFVDP = 0.115)
    propSd_ftctp <- 0.307; label("FTC-TP IC proportional residual SD (fraction; sqrt of variance 0.0942)") # Chen 2016 Table 1B (sigma_FTCTP = 0.0942)
    expSd_datp   <- 0.508; label("dATP log-scale residual SD (unitless; sqrt of variance 0.258)")        # Chen 2016 Table 1A (sigma_dATP = 0.258)
    expSd_dgtp   <- 0.517; label("dGTP log-scale residual SD (unitless; sqrt of variance ~0.267 derived from reported %CV 51.7; see Errata)") # Chen 2016 Table 1A (sigma_dGTP %CV = 51.7; estimate value missing from the extracted table)
    expSd_dctp   <- 0.479; label("dCTP log-scale residual SD (unitless; sqrt of variance 0.229)")        # Chen 2016 Table 1B (sigma_dCTP = 0.229)
    expSd_ttp    <- 0.597; label("TTP log-scale residual SD (unitless; sqrt of variance 0.356)")         # Chen 2016 Table 1B (sigma_TTP = 0.356)
  })

  model({
    # ------------------------------------------------------------
    # 1. Derived covariate terms
    # ------------------------------------------------------------
    # Reconstruct the paper's male indicator from canonical SEXF
    # so the additive +24.3 L Vc/F shift is applied to males
    # exactly as in Chen 2016 Table 1B (tvVc/F_FTC = 99.4 + 24.3 *
    # source_SEX; source_SEX = 1 - SEXF).
    sexm <- 1 - SEXF

    # ------------------------------------------------------------
    # 2. Individual PK parameters (typical-value * exp(eta) where
    # an IIV is estimated; otherwise typical-value).
    # ------------------------------------------------------------
    # TFV plasma
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    # FTC plasma. The Vc/F covariate effect is additive linear on
    # the typical value: tvVc/F = exp(lvc_ftc) + e_sexm_vc_ftc *
    # SEXM. IIV is then applied multiplicatively on the typical
    # value (the paper's exponential random-effect model).
    ka_ftc <- exp(lka_ftc)
    cl_ftc <- exp(lcl_ftc + etalcl_ftc)
    tv_vc_ftc <- exp(lvc_ftc) + e_sexm_vc_ftc * sexm
    vc_ftc <- tv_vc_ftc * exp(etalvc_ftc)
    q_ftc  <- exp(lq_ftc)
    vp_ftc <- exp(lvp_ftc + etalvp_ftc)

    # ------------------------------------------------------------
    # 3. Micro-constants for the two plasma 2-cmt models.
    # ------------------------------------------------------------
    kel     <- cl / vc
    k12     <- q  / vc
    k21     <- q  / vp
    kel_ftc <- cl_ftc / vc_ftc
    k12_ftc <- q_ftc  / vc_ftc
    k21_ftc <- q_ftc  / vp_ftc

    # ------------------------------------------------------------
    # 4. Intracellular TFV-DP / FTC-TP parameters.
    # ------------------------------------------------------------
    kf_tfvdp   <- exp(lkf_tfvdp + etalkf_tfvdp)
    sc50_tfvdp <- exp(lsc50_tfvdp)
    kel_tfvdp  <- exp(lkel_tfvdp + etalkel_tfvdp)

    # FTC-TP Kf: additive HIV effect on the typical value
    # (tvKf = 41.6 + 31.3 * HIV_POS), then multiplicative IIV on
    # the resulting typical value.
    tv_kf_ftctp <- exp(lkf_ftctp) + e_hiv_pos_kf_ftctp * HIV_POS
    kf_ftctp    <- tv_kf_ftctp * exp(etalkf_ftctp)
    sc50_ftctp  <- exp(lsc50_ftctp)
    kel_ftctp   <- exp(lkel_ftctp + etalkel_ftctp)

    # ------------------------------------------------------------
    # 5. Endogenous dNTP individual baseline (R0) and EC50.
    # The indirect-response Kout is fixed to 1/day and Emax fixed
    # to 1; baseline equals R0 (since K0in = R0 * Kout = R0).
    # ------------------------------------------------------------
    r0_datp   <- exp(lr0_datp   + etalr0_datp)
    ec50_datp <- exp(lec50_datp + etalec50_datp)

    r0_dgtp   <- exp(lr0_dgtp   + etalr0_dgtp)
    ec50_dgtp <- exp(lec50_dgtp)

    r0_dctp   <- exp(lr0_dctp   + etalr0_dctp)
    ec50_dctp <- exp(lec50_dctp + etalec50_dctp)

    r0_ttp    <- exp(lr0_ttp    + etalr0_ttp)
    ec50_ttp  <- exp(lec50_ttp  + etalec50_ttp)

    # ------------------------------------------------------------
    # 6. Plasma 2-compartment ODEs for TFV and FTC. Doses enter
    # the respective depot; users dose 'depot' with the TFV-
    # equivalent of TDF (= 136 mg per 300 mg TDF tablet, per
    # Chen 2016 Methods 'Study design') and 'depot_ftc' with the
    # FTC dose (= 200 mg per tablet).
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    d/dt(depot_ftc)       <- -ka_ftc * depot_ftc
    d/dt(central_ftc)     <-  ka_ftc * depot_ftc - kel_ftc * central_ftc -
                              k12_ftc * central_ftc + k21_ftc * peripheral1_ftc
    d/dt(peripheral1_ftc) <-  k12_ftc * central_ftc - k21_ftc * peripheral1_ftc

    # ------------------------------------------------------------
    # 7. Virtual saturation compartments for the hybrid TFV-DP and
    # FTC-TP formation link. Per Chen 2016 Methods 'Population
    # pharmacokinetics modeling of TFV-DP/FTC-TP', the saturation
    # compartment receives Kf * Cp_parent (Cp = central / vc) and
    # eliminates at the plasma elimination rate constant
    # (kel_parent = CL/Vc), so its trajectory mirrors plasma with
    # a delay. The compartment is virtual (apparent volume = 1
    # unit), so the state value is numerically the SC50-scaled
    # concentration used in the formation-rate inhibition term.
    #
    # Unit alignment: dose is in mg and vc is in L, so central/vc is
    # in mg/L = ug/mL. The paper reports plasma TFV/FTC in ng/mL and
    # calibrated Kf and SC50 in those native units; multiply by 1000
    # to convert mg/L -> ng/mL for the saturation and intracellular
    # formation drivers.
    # ------------------------------------------------------------
    Cp_tfv <- (central     / vc)     * 1000   # mg/L -> ng/mL
    Cp_ftc <- (central_ftc / vc_ftc) * 1000   # mg/L -> ng/mL

    d/dt(effect_tfvdp) <- kf_tfvdp * Cp_tfv - kel     * effect_tfvdp
    d/dt(effect_ftctp) <- kf_ftctp * Cp_ftc - kel_ftc * effect_ftctp

    # Formation-rate saturation factor for the IC central compart-
    # ments: SC50 / (SC50 + sat_concentration).
    sat_factor_tfvdp <- sc50_tfvdp / (sc50_tfvdp + effect_tfvdp)
    sat_factor_ftctp <- sc50_ftctp / (sc50_ftctp + effect_ftctp)

    # ------------------------------------------------------------
    # 8. Intracellular TFV-DP and FTC-TP ODEs with 2-compartment
    # "recycle" elimination. Per Chen 2016: K46 = K64 = R * Kel
    # (rate between IC central and recycle); K40 = (1 - R) * Kel
    # (true elimination from IC central). Consolidated, the net
    # rate on IC central is -kel * tfvdp + R * kel * recycle, and
    # the recycle compartment has +R*kel*tfvdp - R*kel*recycle.
    # Concentrations: fmol per 10^6 PBMC; the compartment volume
    # is implicit (= 1 unit of PBMC pool) so state values are
    # numerically equal to per-cell concentrations.
    # ------------------------------------------------------------
    d/dt(tfvdp)             <-  kf_tfvdp * Cp_tfv * sat_factor_tfvdp -
                                 kel_tfvdp * tfvdp +
                                 r_tfvdp * kel_tfvdp * peripheral1_tfvdp
    d/dt(peripheral1_tfvdp) <-  r_tfvdp * kel_tfvdp * tfvdp -
                                 r_tfvdp * kel_tfvdp * peripheral1_tfvdp

    d/dt(ftctp)             <-  kf_ftctp * Cp_ftc * sat_factor_ftctp -
                                 kel_ftctp * ftctp +
                                 r_ftctp * kel_ftctp * peripheral1_ftctp
    d/dt(peripheral1_ftctp) <-  r_ftctp * kel_ftctp * ftctp -
                                 r_ftctp * kel_ftctp * peripheral1_ftctp

    # ------------------------------------------------------------
    # 9. Endogenous dNTP indirect-response ODEs. Indirect-response
    # Emax with Emax = 1 and Kout = 1/day fixed (Chen 2016 Methods
    # 'Population pharmacodynamics modeling'). K0in = R0 * Kout =
    # R0, so the baseline is R0 and the deviation is driven by
    # the analog-mediated inhibition term.
    #
    # The dGTP inhibition is modulated by a time-waning factor
    # 1/(1 + t^gamma) (Eq D in S1 File; simplified sigmoidal
    # E'max-vs-time form): the effective inhibition fraction
    # approaches the unmodulated value 1/(1 + sc50/(sc50 + tfvdp))
    # at t = 0 and decays toward zero as t increases.
    # ------------------------------------------------------------
    kout_dntp <- 1.0   # Kout fixed to 1/day per Chen 2016 Methods
    emax_dntp <- 1.0   # Emax fixed to 1   per Chen 2016 Methods

    inhib_datp <- emax_dntp * tfvdp / (ec50_datp + tfvdp)
    inhib_dgtp <- emax_dntp * tfvdp / (ec50_dgtp + tfvdp) / (1 + t ^ gamma_dgtp)
    inhib_dctp <- emax_dntp * ftctp / (ec50_dctp + ftctp)
    inhib_ttp  <- emax_dntp * ftctp / (ec50_ttp  + ftctp)

    d/dt(datp) <- r0_datp * kout_dntp * (1 - inhib_datp) - kout_dntp * datp
    d/dt(dgtp) <- r0_dgtp * kout_dntp * (1 - inhib_dgtp) - kout_dntp * dgtp
    d/dt(dctp) <- r0_dctp * kout_dntp * (1 - inhib_dctp) - kout_dntp * dctp
    d/dt(ttp)  <- r0_ttp  * kout_dntp * (1 - inhib_ttp)  - kout_dntp * ttp

    # Initial conditions: each dNTP starts at its individual
    # baseline R0 (pre-dose steady state in the absence of drug).
    datp(0) <- r0_datp
    dgtp(0) <- r0_dgtp
    dctp(0) <- r0_dctp
    ttp(0)  <- r0_ttp

    # ------------------------------------------------------------
    # 10. Observations.
    # Cc        -- TFV plasma concentration (ng/mL).
    # Cc_ftc    -- FTC plasma concentration (ng/mL).
    # Cc_tfvdp  -- intracellular TFV-DP (fmol per 10^6 PBMC).
    # Cc_ftctp  -- intracellular FTC-TP (fmol per 10^6 PBMC).
    # datp,...  -- dNTP concentrations (fmol per 10^6 PBMC).
    # ------------------------------------------------------------
    Cc       <- Cp_tfv
    Cc_ftc   <- Cp_ftc
    Cc_tfvdp <- tfvdp
    Cc_ftctp <- ftctp

    Cc       ~ lnorm(expSd)
    Cc_ftc   ~ lnorm(expSd_ftc)
    Cc_tfvdp ~ prop(propSd_tfvdp)
    Cc_ftctp ~ prop(propSd_ftctp)
    datp     ~ lnorm(expSd_datp)
    dgtp     ~ lnorm(expSd_dgtp)
    dctp     ~ lnorm(expSd_dctp)
    ttp      ~ lnorm(expSd_ttp)
  })
}

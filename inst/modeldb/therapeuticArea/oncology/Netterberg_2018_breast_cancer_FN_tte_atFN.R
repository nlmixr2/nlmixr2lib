Netterberg_2018_breast_cancer_FN_tte_atFN <- function() {
  description <- "Parametric time-to-event (TTE) submodel for febrile neutropenia (FN) in early breast cancer patients (n=49) receiving adjuvant FEC and docetaxel chemotherapy: the 'when-FN-occurs' variant whose hazard depends on the instantaneous log of normalized CRP. The hazard form is h(t) = h0 * exp(beta4 * LN_CRP(t)) where LN_CRP(t) = log(CRP(t) / CRP0) (Netterberg 2018 Eq. 3c). The underlying joint IL-6 / CRP biomarker turnover model (see Netterberg_2018_breast_cancer_FN_biomarkers) is embedded inline so this submodel can be simulated standalone. The hazard is fixed to 0 for t < 84 h (3.5 days) per the paper."
  reference <- paste(
    "Netterberg I, Karlsson MO, Nielsen EI, Quartino AL, Lindman H, Friberg LE.",
    "The risk of febrile neutropenia in breast cancer patients following adjuvant",
    "chemotherapy is predicted by the time course of interleukin-6 and C-reactive",
    "protein by modelling.",
    "Br J Clin Pharmacol. 2018;84(3):490-500.",
    "doi:10.1111/bcp.13477.",
    "Companion biomarker model (same paper):",
    "modellib('Netterberg_2018_breast_cancer_FN_biomarkers').",
    sep = " "
  )
  vignette <- "Netterberg_2018_breast_cancer_FN"
  units <- list(
    time = "hour",
    dosing = "n/a (no drug doses; the cycle anchors t = 0 at the start of chemotherapy)",
    concentration = "IL-6 in pg/mL; CRP in mg/L; `sur` is a survival probability"
  )

  covariateData <- list(
    MIX_ELEV_IL6 = list(
      description = "Binary cycle-level indicator: 1 = elevated IL-6 production surge active in this cycle; 0 = no elevated IL-6 production. Time-fixed within the cycle.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Carried here because CRP production is regulated by IL-6 RCFB through the Slope coefficient. For simulation, draw as a Bernoulli random variable with probability Pelevation_IL6 = 0.634 (Netterberg 2018 Table 2).",
      source_name = "mixture component (NONMEM MIXTURE/PMIX)"
    ),
    MIX_ELEV_CRP = list(
      description = "Binary cycle-level indicator: 1 = elevated CRP production surge active in this cycle; 0 = no elevated CRP production.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "For simulation, draw as a Bernoulli random variable with probability Pelevation_CRP = 0.443 (Netterberg 2018 Table 2).",
      source_name = "mixture component (NONMEM MIXTURE/PMIX)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "31-73 years; median 54 years (Netterberg 2018 Table 1)",
    weight_range = "54-111 kg; median 70 kg (Netterberg 2018 Table 1)",
    sex_female_pct = 100,
    disease_state = "Early breast cancer; receiving adjuvant chemotherapy.",
    dose_range = "FEC (epirubicin 75 mg/m^2 + 5-FU 600 mg/m^2 + cyclophosphamide 600 mg/m^2) and docetaxel 80 mg/m^2 Q3W (typical regimen).",
    regions = "Sweden (Uppsala University Hospital).",
    notes = "Hazard delayed: h(t) = 0 for t < 84 h (3.5 days). Per Netterberg 2018 Discussion, the hazard ratio is ~5.0 for a CRP of 10 mg/L vs 5 mg/L (relative to CRP0 = 1.88 mg/L: exp(2.33 * (log(10/1.88) - log(5/1.88))) = exp(2.33 * log(2)) = 5.03)."
  )

  ini({
    # ----- Biomarker turnover parameters (inherited from the biomarker model; Table 2) -----
    lbl_il6   <- log(2.50);   label("Typical baseline IL-6 IL-60 (pg/mL)")                # Netterberg 2018 Table 2 IL-60 = 2.50
    lbl_crp   <- log(1.88);   label("Typical baseline CRP CRP0 (mg/L)")                   # Netterberg 2018 Table 2 CRP0 = 1.88
    lkout_il6 <- log(0.0141); label("IL-6 turnover rate constant kout_IL-6 (1/h)")          # Netterberg 2018 Table 2 kout,IL-6 = 0.0141
    lkout_crp <- log(0.0224); label("CRP turnover rate constant kout_CRP (1/h)")            # Netterberg 2018 Table 2 kout,CRP = 0.0224
    lsa_il6   <- log(7.99);   label("IL-6 surge amplitude SA_IL-6 (unitless)")                          # Netterberg 2018 Table 2 SA_IL-6 = 7.99
    lsa_crp   <- log(4.40);   label("CRP surge amplitude SA_CRP (unitless)")                            # Netterberg 2018 Table 2 SA_CRP = 4.40
    lsw_il6   <- log(32.4);   label("IL-6 surge width SW_IL-6 (h)")                          # Netterberg 2018 Table 2 SW_IL-6 = 32.4
    lsw_crp   <- log(53.8);   label("CRP surge width SW_CRP (h)")                            # Netterberg 2018 Table 2 SW_CRP = 53.8
    lpt_il6   <- log(137.0);  label("IL-6 surge peak time PT_IL-6 (h)")                      # Netterberg 2018 Table 2 PT_IL-6 = 137
    lpt_crp_plus <- log(50.3); label("PT_CRP+ added to PT_IL-6 to get PT_CRP (h)")           # Netterberg 2018 Table 2 PT_CRP+ = 50.3
    slope_il6_crp <- 1.05;    label("Slope of RCFB_IL-6(t) into CRP production")             # Netterberg 2018 Table 2 Slope = 1.05

    # ----- TTE when-FN-occurs parameters (Eq. 3c; Table 2 when-FN-occurs column) -----
    lh0   <- log(7.61e-5); label("Log baseline FN hazard h0 (1/h)")                       # Netterberg 2018 Table 2 when-FN-occurs h0 = 7.61e-5 (RSE 120%)
    beta4 <- 2.33;          label("Hazard log-linear coefficient on LN_CRP(t)")             # Netterberg 2018 Table 2 when-FN-occurs beta4 = 2.33 (RSE 12%); paper Table units LN_CRP(t)^-1

    # ----- Biomarker IIV / IOV (mirror the biomarker model) -----
    etalbl_il6     ~ 0.3858  # IIV(IL-60) = 68.0% CV
    etalbl_crp     ~ 0.4862  # IIV(CRP0) = 80.5% CV
    etalkout_il6   ~ 0.9883  # IIV(kout,IL-6) = 130% CV
    etalsa_crp     ~ 0.3198  # IOV(SA_CRP) = 61.4% CV
    etalsw_crp     ~ 0.5219  # IOV(SW_CRP) = 83.8% CV
    etalpt_il6     ~ 0.3036  # IOV(PT_IL-6) = 59.7% CV
    etalpt_crp_plus ~ 0.4920 # IOV(PT_CRP+) = 81.3% CV

    # ----- Biomarker residual errors -----
    propSd_Cc_il6 <- 0.547; label("Proportional residual error on IL-6 (fraction)")  # Netterberg 2018 Table 2
    propSd_Cc_crp <- 0.530; label("Proportional residual error on CRP (fraction)")   # Netterberg 2018 Table 2
  })

  model({
    # ----- Biomarker layer (same as Netterberg_2018_breast_cancer_FN_biomarkers) -----
    bl_il6     <- exp(lbl_il6   + etalbl_il6)
    bl_crp     <- exp(lbl_crp   + etalbl_crp)
    kout_il6_i <- exp(lkout_il6 + etalkout_il6)
    kout_crp_i <- exp(lkout_crp)

    sa_il6_i      <- exp(lsa_il6)
    sa_crp_i      <- exp(lsa_crp      + etalsa_crp)
    sw_il6_i      <- exp(lsw_il6)
    sw_crp_i      <- exp(lsw_crp      + etalsw_crp)
    pt_il6_i      <- exp(lpt_il6      + etalpt_il6)
    pt_crp_plus_i <- exp(lpt_crp_plus + etalpt_crp_plus)
    pt_crp_i      <- pt_il6_i + pt_crp_plus_i

    g_il6 <- MIX_ELEV_IL6 * sa_il6_i / (((t - pt_il6_i) / sw_il6_i)^4 + 1)
    g_crp <- MIX_ELEV_CRP * sa_crp_i / (((t - pt_crp_i) / sw_crp_i)^4 + 1)

    rin_il6 <- kout_il6_i * bl_il6
    rin_crp <- kout_crp_i * bl_crp

    d/dt(il6) <- rin_il6 * (1 + g_il6) - kout_il6_i * il6
    il6(0)    <- bl_il6

    rcfb_il6 <- (il6 - bl_il6) / bl_il6
    d/dt(crp) <- rin_crp * (1 + g_crp + slope_il6_crp * rcfb_il6) - kout_crp_i * crp
    crp(0)    <- bl_crp

    # ----- Hazard (Eq. 3c) with 3.5-day delay -----
    ln_crp  <- log(crp / bl_crp)
    h0      <- exp(lh0)
    haz_raw <- h0 * exp(beta4 * ln_crp)
    hazard  <- ifelse(t < 84, 0, haz_raw)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    # ----- Biomarker observations -----
    Cc_il6 <- il6
    Cc_crp <- crp
    Cc_il6 ~ prop(propSd_Cc_il6)
    Cc_crp ~ prop(propSd_Cc_crp)
  })
}

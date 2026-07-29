Luu_2017_nusinersen <- function() {
  description <- "Four-compartment population PK model for nusinersen (antisense oligonucleotide) following intrathecal administration in pediatric patients with spinal muscular atrophy (Luu 2017): a CSF + CNS-tissue subsystem (intrathecal bolus enters the CSF) coupled by a unidirectional CSF-to-plasma transport to a plasma + systemic-tissue subsystem, with baseline body weight as a power covariate on CL_p and V_CSF and a linear covariate on V_p."
  reference <- "Luu KT, Norris DA, Gunawan R, Henry S, Geary R, Wang Y. Population Pharmacokinetics of Nusinersen in the Cerebral Spinal Fluid and Plasma of Pediatric Patients With Spinal Muscular Atrophy Following Intrathecal Administrations. J Clin Pharmacol. 2017;57(8):1031-1041. doi:10.1002/jcph.884"
  vignette <- "Luu_2017_nusinersen"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # cns_tissue is the lumped CNS-tissue compartment that pairs with csf in
  # the CSF subsystem. The paper describes it as "a lumped compartment
  # consisting possibly of the spinal cord tissue, subarachnoid space, and
  # brain tissue" (Luu 2017 Discussion), so it does not map onto the
  # canonical brain_<region> namespace or the standard peripheral1.
  paper_specific_compartments <- c("cns_tissue")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight per Luu 2017 Table 2 / Equations 5-7: power scaling on CL_p (exponent 0.689) and V_CSF (exponent 0.596), and linear proportional-deviation scaling on V_p (coefficient 0.047/kg). Reference body weight MBWT = 15.2 kg (median baseline BWT in the pooled dataset, Table 1).",
      source_name        = "BWT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 72L,
    n_studies      = 5L,
    studies        = "Pooled: ISIS 396443-CS1, -CS2, -CS3A, -CS10, -CS12",
    age_range      = "0.10 to 17 years (median 5)",
    weight_range   = "5.1 to 83 kg (median 15.2)",
    sex_female_pct = 48.6,
    race_ethnicity = c(White = 86.1, Black = 5.6, Asian = 4.2, Other = 4.2),
    disease_state  = "Pediatric patients with spinal muscular atrophy (SMA) spanning presymptomatic / infantile-onset (likely type I) and later-onset (likely type II/III) phenotypes",
    dose_range     = "Single or repeated intrathecal doses of 1, 3, 6, 9, or 12 mg",
    administration_routes = "Intrathecal bolus injection",
    notes          = "Pooled CSF and plasma data: 279 CSF and 1181 plasma concentration data points across 5 trials (Luu 2017 Table 1). Estimation with NONMEM 7.2 first-order conditional estimation with interaction (FOCE-I)."
  )

  ini({
    # ---- Structural parameters: final-model typical values (Luu 2017 Table 2) ----
    lcl          <- log(2.36);   label("Plasma clearance CL_p (L/h)")                                                  # Luu 2017 Table 2 final: 2.36 L/h (%RSE 5.04)
    lcl_csf      <- log(0.136);  label("CSF clearance CL_CSF (one-way CSF to plasma transport; L/h)")                  # Luu 2017 Table 2 final: 0.136 L/h (%RSE 9.12)
    lq           <- log(0.568);  label("Plasma inter-compartmental clearance Q_p (central <-> peripheral1; L/h)")      # Luu 2017 Table 2 final: 0.568 L/h (%RSE 11.7)
    lqcsf        <- log(0.0712); label("CSF inter-compartmental clearance Q_CSF (csf <-> cns_tissue; L/h)")            # Luu 2017 Table 2 final: 0.0712 L/h (%RSE 19.5)
    lvc          <- log(29.0);   label("Plasma volume V_p (L)")                                                        # Luu 2017 Table 2 final: 29.0 L (%RSE 11.4)
    lvp          <- log(418);    label("Systemic-tissue volume V_systemic_tissue (L)")                                 # Luu 2017 Table 2 final: 418 L (%RSE 20.9)
    lvcsf        <- log(0.433);  label("CSF volume V_CSF (L)")                                                         # Luu 2017 Table 2 final: 0.433 L (%RSE 21.7)
    lvcns        <- log(263);    label("CNS-tissue volume V_CNS_tissue (L)")                                           # Luu 2017 Table 2 final: 263 L (%RSE 18.9)

    # ---- Baseline body-weight covariate effects (Luu 2017 Table 2; Equations 5-7) ----
    # MBWT (median baseline body weight) = 15.2 kg (Luu 2017 Table 1).
    e_wt_cl     <- 0.689;  label("Power exponent on (BWT/MBWT) for CL_p")                          # Luu 2017 Table 2 final: 0.689 (%RSE 12.2)
    e_wt_vcsf   <- 0.596;  label("Power exponent on (BWT/MBWT) for V_CSF")                         # Luu 2017 Table 2 final: 0.596 (%RSE 52.3)
    e_wt_vc     <- 0.047;  label("Linear coefficient on (BWT - MBWT) for V_p (1/kg)")              # Luu 2017 Table 2 final: 0.047 (%RSE 21.7)

    # ---- IIV (Luu 2017 Table 2 final-model %IIV; omega^2 = ln(1 + CV^2)) ----
    etalcl     ~ 0.0834   # CL_p   IIV 29.5%  -> ln(1 + 0.295^2)
    etalcl_csf ~ 0.0578   # CL_CSF IIV 24.4%  -> ln(1 + 0.244^2)
    etalvcsf   ~ 0.5746   # V_CSF  IIV 88.1%  -> ln(1 + 0.881^2)
    etalvc     ~ 0.1416   # V_p    IIV 39.0%  -> ln(1 + 0.390^2)
    # Inter-occasion variability of 38.1% CV was estimated on CL_CSF in the
    # original NONMEM model (Luu 2017 Table 2 IOV row) but is not encoded
    # here -- nlmixr2 would require an explicit per-occasion OCC variable
    # and a parallel family of occasion etas. See vignette "Assumptions and
    # deviations".

    # ---- Residual unexplained variability (Luu 2017 Table 2) ----
    # The paper reports a single proportional residual term: epsilon
    # (Proportional) = 0.494 (%RSE 2.45). Per the nlmixr2lib convention
    # for NONMEM 7-source residuals, the reported value is interpreted as
    # the $SIGMA variance, so the proportional SD is sqrt(0.494) = 0.703
    # (~70% CV). It is applied identically to plasma Cc and CSF Ccsf
    # because the paper does not separate the two endpoints. See vignette
    # "Assumptions and deviations".
    propSd      <- sqrt(0.494); label("Plasma proportional residual SD (fraction; sqrt(NONMEM $SIGMA = 0.494))")
    propSd_Ccsf <- sqrt(0.494); label("CSF proportional residual SD (fraction; sqrt(NONMEM $SIGMA = 0.494))")
  })

  model({
    # Reference body weight MBWT = median baseline BWT in the pooled
    # dataset (Luu 2017 Table 1: All cohort weight median 15.2 kg).
    MBWT <- 15.2

    # Individual structural parameters (Luu 2017 Equations 5, 6, 7).
    cl     <- exp(lcl     + etalcl)     * (WT / MBWT)^e_wt_cl        # Equation 7
    cl_csf <- exp(lcl_csf + etalcl_csf)
    q      <- exp(lq)
    qcsf   <- exp(lqcsf)
    vc     <- exp(lvc     + etalvc)     * (1 + e_wt_vc * (WT - MBWT)) # Equation 6
    vp     <- exp(lvp)
    vcsf   <- exp(lvcsf   + etalvcsf)   * (WT / MBWT)^e_wt_vcsf       # Equation 5
    vcns   <- exp(lvcns)

    # Concentrations. Volume in L, mass in ug after f(csf) unit conversion:
    # mass [ug] / volume [L] = ng/mL.
    Cc   <- central     / vc     # plasma concentration (ng/mL)
    Ccsf <- csf         / vcsf   # CSF concentration (ng/mL)
    Ccns <- cns_tissue  / vcns   # CNS-tissue concentration (derived; not observed)
    Cp   <- peripheral1 / vp     # systemic-tissue concentration (derived; not observed)

    # Structural ODE system (Luu 2017 Figure 1).
    # CSF subsystem: csf <-> cns_tissue via Q_CSF (bidirectional). Intrathecal
    # bolus dose enters csf via the dosing event in the user data
    # (cmt = "csf"). Inter-system coupling: csf -> central via CL_CSF
    # (one-way; Luu 2017 Discussion: "a single transport rate constant
    # from CSF to plasma but no transport rate constant from the plasma
    # to the CSF"). Plasma subsystem: central <-> peripheral1 via Q_p;
    # central -> elimination via CL_p.
    d/dt(csf)         <- -qcsf * Ccsf + qcsf * Ccns - cl_csf * Ccsf
    d/dt(cns_tissue)  <-  qcsf * Ccsf - qcsf * Ccns
    d/dt(central)     <-  cl_csf * Ccsf - q * Cc + q * Cp - cl * Cc
    d/dt(peripheral1) <-  q * Cc - q * Cp

    # Intrathecal dose unit conversion: input dose in mg, csf state in ug
    # so that csf [ug] / vcsf [L] yields concentrations in ng/mL.
    f(csf) <- 1000

    Cc   ~ prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}

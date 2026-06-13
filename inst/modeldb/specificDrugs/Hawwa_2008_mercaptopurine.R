Hawwa_2008_mercaptopurine <- function() {
  description <- paste(
    "Population PK / pharmacogenetic model for oral 6-mercaptopurine (6-MP)",
    "and its two active intracellular metabolites 6-thioguanine nucleotides",
    "(6-TGNs) and 6-methylmercaptopurine nucleotides (6-mMPNs) measured in",
    "erythrocytes (RBCs) of 19 paediatric patients (n = 75 samples; 150",
    "concentrations) with acute lymphoblastic leukaemia receiving",
    "maintenance chemotherapy at a target oral dose of 75 mg/m^2/day. The",
    "structural model is a one-compartment first-order absorption + first-",
    "order elimination model for 6-MP whose plasma concentration is NOT",
    "observed; transformation of bioavailable 6-MP into the two RBC",
    "metabolites occurs at a common metabolic rate constant kme = 0.78 *",
    "k20 (78% of total 6-MP elimination), with the fractional split",
    "between metabolites governed by FM3 (the fractional metabolic",
    "transformation of 6-MP into 6-TGNs; the complementary fraction 1 -",
    "FM3 goes to 6-mMPNs). The fixed structural anchors ka = 1.3 1/h, F",
    "= 0.22, k20 = 0.53 1/h, and the kme : k20 ratio of 0.78 are taken",
    "from the prior 6-MP literature (Hawwa 2008 cites Zimm 1983 and",
    "Lennard 1990; values listed in Methods page 4 of the British",
    "Journal of Clinical Pharmacology article). The apparent",
    "distribution volume of 6-MP central and of each metabolite",
    "compartment is not identifiable from the RBC sampling design and is",
    "fixed to 1 L by the ADVAN6 implementation convention used by",
    "Hawwa 2008 (analogous to Urien 2005 capecitabine; see vignette",
    "Errata for the dimensional-analysis discussion). The only",
    "estimated structural parameters are FM3, CL_6TGNs, and CL_6mMPNs;",
    "the only retained covariates are TPMT genotype (any TPMT*3*",
    "mutation, pooled binary) on FM3 (theta = 2.56) and body surface",
    "area (m^2) on the apparent clearance of 6-TGNs as a power-law",
    "(theta = 1.16, anchored to BSA = 1 m^2). IIV is estimated only on",
    "the two metabolite clearances; no IIV on FM3 was estimable."
  )
  reference <- "Hawwa AF, Collier PS, Millership JS, McCarthy A, Dempsey S, Cairns C, McElnay JC. Population pharmacokinetic and pharmacogenetic analysis of 6-mercaptopurine in paediatric patients with acute lymphoblastic leukaemia. Br J Clin Pharmacol. 2008;66(6):826-837. doi:10.1111/j.1365-2125.2008.03281.x"
  vignette  <- "Hawwa_2008_mercaptopurine"
  units     <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline, computed from height and weight by the formula the source paper applied (Hawwa 2008 Table 1 lists BSA as a per-patient summary; the formula was not stated explicitly).",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on the apparent clearance of 6-TGNs anchored to BSA = 1 m^2: CL_6TGNs = theta_CL_6TGNs * BSA^theta_BSA with theta_BSA = 1.16 (Hawwa 2008 Table 4). The cohort BSA range was 0.59 - 2.00 m^2 (median 1.14 m^2); the model can be evaluated outside that range but extrapolation to adult-size BSA (>= 2 m^2) is not supported by the source data. Time-fixed at baseline within the source dataset.",
      source_name        = "BSA"
    ),
    TPMT_MUT = list(
      description        = "Pooled binary indicator of any reduced-function TPMT variant allele (TPMT*3A, TPMT*3B, or TPMT*3C; heterozygous or homozygous) -- i.e., TPMT genotype-derived thiopurine S-methyltransferase deficiency status.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (TPMT*1/*1 wild-type across the TPMT*3A, TPMT*3B, and TPMT*3C panel; full canonical TPMT activity).",
      notes              = "Hawwa 2008 Table 1 reports the cohort genotyping panel as TPMT*3A (1 heterozygote / 0 homozygotes), TPMT*3B (none observed), and TPMT*3C (2 heterozygotes / 0 homozygotes); the TPMT*2 allele was not assayed. The covariate enters the model as a power-of-binary multiplicative factor on the fractional metabolic transformation of 6-MP into 6-TGNs: FM3 = TVFM3 * theta_TPMT^TPMT_MUT with theta_TPMT = 2.56 (Hawwa 2008 Table 4). Carriers therefore have FM3 = 0.0489 (vs 0.0191 in wild-type), a 156% increase consistent with shunting of 6-MP away from TPMT-mediated mMPN methylation towards 6-TGN production.",
      source_name        = "TPMT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 19L,
    n_studies      = 1L,
    age_range      = "3 - 17 years",
    age_median     = "10 years",
    weight_range   = "13.2 - 77.5 kg",
    weight_median  = "33.4 kg",
    bsa_range      = "0.59 - 2.00 m^2",
    bsa_median     = "1.14 m^2",
    sex_female_pct = 32.0,
    disease_state  = "Paediatric patients with acute lymphoblastic leukaemia (ALL) on continuous / maintenance 6-mercaptopurine + weekly methotrexate maintenance chemotherapy. Children had to have been on a constant 6-MP daily dose for at least one week and had received no intensification therapy or red blood cell transfusion within the preceding two months.",
    dose_range     = "Target 6-MP oral dose 75 mg/m^2/day; per-patient daily-dose range 10 - 100 mg (median 50 mg) titrated to leucocyte count and clinically relevant infection status. Co-medications during 6-MP chemotherapy: methotrexate 5 - 25 mg/week (median 15), cotrimoxazole 120 - 480 mg b.d. twice weekly (median 360); monthly IV vincristine.",
    regions        = "Northern Ireland (Royal Belfast Hospital for Sick Children, Belfast Health and Social Care Trust).",
    n_observations = "150 erythrocyte metabolite concentrations across 75 sampling occasions (one sample per occasion, up to five occasions per patient, sampled at least 12 h after the preceding 6-MP dose and prior to monthly vincristine administration). Patients were assigned randomly to an index group (n = 15) for model development and a validation group (n = 4) for internal predictive-performance assessment; the final-model parameters in this file come from refitting the FINAL covariate model to the full dataset (n = 19).",
    baseline_chem  = "Haemoglobin 10.9 - 16.5 g/dL (median 12.9), WBC 1.2 - 9.1 x 10^9/L (median 3.3), platelets 66 - 648 x 10^9/L (median 282), absolute neutrophil count 0.3 - 8.3 x 10^9/L (median 1.68). No clinically significant renal impairment in any patient.",
    notes          = "Demographics from Hawwa 2008 Table 1. Modelling software NONMEM VI level 1.1 (double precision), DIGITAL Visual Fortran compiler v5.0.A, PREDPP subroutine ADVAN6 (general nonlinear ODE solver), first-order conditional estimation (FOCE) with INTERACTION. The genotyped polymorphisms outside TPMT (XO A1936G and A2107G; ITPA C94A and IVS2+21A>C) were screened in univariate analysis but did not reach the OFV-retention threshold and were not included in the FINAL model -- they are not represented in this file."
  )

  ini({
    # Structural anchors (fixed per Hawwa 2008 Methods page 4 of the British
    # Journal of Clinical Pharmacology article, citing Zimm 1983 and
    # Lennard 1990 as the literature sources). The 6-MP plasma concentration
    # is NOT observed in this study; only the two RBC metabolites are
    # measured. The 6-MP central compartment carries only the absorption
    # and first-order elimination dynamics that feed the two metabolite
    # compartments.
    lka      <- fixed(log(1.3))
    label("6-MP first-order absorption rate constant ka (1/h) -- FIXED per Hawwa 2008 Methods (literature)") # Hawwa 2008 Methods page 4: "ka was fixed at 1.3 1/h according to the literature [1, 9]"
    lfdepot  <- fixed(log(0.22))
    label("6-MP oral bioavailability F (unitless) -- FIXED per Hawwa 2008 Methods (literature)") # Hawwa 2008 Methods page 4: "the bioavailability factor (F) of the model was fixed at 22% according to the literature [1, 9]"
    lvc      <- fixed(log(1))
    label("6-MP central apparent volume V_central (L) -- FIXED to 1 L (not identifiable; 6-MP plasma not observed)") # Hawwa 2008 ADVAN6 implementation convention; same as Urien 2005 capecitabine metabolites (see vignette Errata)
    lcl      <- fixed(log(0.53))
    label("6-MP total apparent clearance CL_central (L/h) -- FIXED per Hawwa 2008 Methods (k20 = 0.53 1/h)") # Hawwa 2008 Methods page 5 equations: k20 = 0.53 1/h fixed per literature [1, 9]; with V_central = 1 L, CL_central = k20 * V_central = 0.53 L/h
    e_fmet   <- fixed(0.78)
    label("Fraction of 6-MP elimination that is metabolic (kme / k20) -- FIXED per Hawwa 2008 Methods") # Hawwa 2008 Methods page 5: "k_other = 0.22 * k20 = 0.1166 1/h" so kme = k20 - k_other = (1 - 0.22) * k20 = 0.78 * k20

    # 6-TGN compartment. Apparent V is not identifiable and is fixed to
    # 1 L per the ADVAN6 implementation convention (see vignette Errata).
    # The apparent clearance CL_6TGNs is the only estimated parameter for
    # the elimination side of this compartment.
    lvc_tgn  <- fixed(log(1))
    label("6-TGN apparent volume V_6TGNs (L) -- FIXED to 1 L (not identifiable)") # Hawwa 2008 ADVAN6 implementation convention (see vignette Errata)
    lcl_tgn  <- log(0.00914)
    label("6-TGN apparent clearance at BSA = 1 m^2 (L/h)") # Hawwa 2008 Table 4 FINAL model: theta_CL_6TGNs = 0.00914 L/h (SE 56.8%); BSA power-law equation TVCL_6TGNs = theta_CL_6TGNs * BSA^theta_BSA

    # 6-mMPN compartment. Apparent V not identifiable and is fixed to
    # 1 L per the same ADVAN6 convention; CL_6mMPNs is estimated with
    # no retained covariate effect.
    lvc_mmpn <- fixed(log(1))
    label("6-mMPN apparent volume V_6mMPNs (L) -- FIXED to 1 L (not identifiable)") # Hawwa 2008 ADVAN6 implementation convention (see vignette Errata)
    lcl_mmpn <- log(0.0228)
    label("6-mMPN apparent clearance (L/h)") # Hawwa 2008 Table 4 FINAL model: theta_CL_6mMPNs = 0.0228 L/h (SE 14.2%); no covariate retained

    # Fractional transformation of 6-MP to 6-TGNs (the complementary
    # fraction FM4 = 1 - FM3 goes to 6-mMPNs).
    lfm3     <- log(0.0191)
    label("Fractional metabolic transformation of 6-MP to 6-TGNs in TPMT wild-type (unitless)") # Hawwa 2008 Table 4 FINAL model: theta_FM3 = 0.0191 (SE 64.1%); covariate equation TVFM3 = theta_FM3 * theta_TPMT^TPMT_MUT

    # Covariate effects (BSA on CL_6TGNs as power-law; TPMT on FM3
    # as a power-of-binary multiplicative factor).
    e_bsa_cl_tgn     <- 1.16
    label("Power-law BSA exponent on 6-TGN apparent clearance (unitless)") # Hawwa 2008 Table 4 FINAL model: theta_BSA = 1.16 (SE 49.0%); equation CL_6TGNs = theta_CL_6TGNs * BSA^1.16
    e_tpmt_mut_fm3   <- log(2.56)
    label("Log of multiplicative TPMT-mutation effect on FM3 (= log(theta_TPMT))") # Hawwa 2008 Table 4 FINAL model: theta_TPMT = 2.56 (SE 35.9%); equation TVFM3 = theta_FM3 * theta_TPMT^TPMT_MUT; e_tpmt_mut_fm3 = log(2.56) so that exp(e_tpmt_mut_fm3 * TPMT_MUT) reproduces 2.56^TPMT_MUT

    # IIV (the variance column in Hawwa 2008 Table 4 is reported as
    # omega^2 on the internal log-scale; the paired CV% column equals
    # sqrt(omega^2) per the Methods page 5 definition "the square root
    # of the variance was interpreted as the CV"). No IIV on FM3 was
    # reported in the FINAL model -- only on the two metabolite
    # clearances.
    etalcl_tgn  ~ 0.113
    label("IIV variance on log 6-TGN clearance (omega^2)") # Hawwa 2008 Table 4 FINAL model: w_CL_6TGNs = 0.113 (sqrt = 33.6% reported as the paired CV%)
    etalcl_mmpn ~ 0.11
    label("IIV variance on log 6-mMPN clearance (omega^2)") # Hawwa 2008 Table 4 FINAL model: w_CL_6mMPNs = 0.11 (sqrt = 33.2% reported as the paired CV%)

    # Residual error. Hawwa 2008 Methods page 5 selected an additive
    # error model on each metabolite RBC concentration after testing
    # additive, proportional, and combined structures. Residual SDs
    # reported in Table 4 are in the linear-concentration space of each
    # metabolite (mg/L of packed RBC).
    addSd_tgn  <- 0.177
    label("6-TGN additive residual SD (mg/L of packed RBC)") # Hawwa 2008 Table 4 FINAL model: s_6TGNs = 0.177 mg/L (SE 17.5%)
    addSd_mmpn <- 8.42
    label("6-mMPN additive residual SD (mg/L of packed RBC)") # Hawwa 2008 Table 4 FINAL model: s_6mMPNs = 8.42 mg/L (SE 19.3%)
  })

  model({
    # 1. Structural anchors (all fixed at their log-transformed ini()
    # values; no IIV).
    ka       <- exp(lka)
    vc       <- exp(lvc)
    cl       <- exp(lcl)
    fmet     <- e_fmet                # kme / k20 ratio (fixed at 0.78)

    # 2. Individual parameters with covariate effects and IIV.
    fm3      <- exp(lfm3 + e_tpmt_mut_fm3 * TPMT_MUT)
    vc_tgn   <- exp(lvc_tgn)
    cl_tgn   <- exp(lcl_tgn + etalcl_tgn) * BSA^e_bsa_cl_tgn
    vc_mmpn  <- exp(lvc_mmpn)
    cl_mmpn  <- exp(lcl_mmpn + etalcl_mmpn)

    # 3. Micro-constants derived from the fixed structural anchors. kel
    # is the total 6-MP elimination rate (= k20 in the paper's notation)
    # and kme is the fraction of that flux destined for either
    # metabolite (the complementary k_other arm is consumed by the
    # central -> outside pathway and is not tracked as a separate
    # compartment).
    kel      <- cl / vc                       # = 0.53 1/h (paper's k20)
    kme      <- fmet * kel                    # = 0.78 * k20 = 0.4134 1/h (paper's kme)

    # 4. ODE system (Hawwa 2008 Methods page 5 differential equations
    # for the ADVAN6 implementation). The depot equation reproduces the
    # analytical ka * dose input; central tracks 6-MP amount in plasma
    # (not observed); the two metabolite compartments receive their
    # input from the metabolic-flux fraction kme and have first-order
    # elimination at the apparent rate cl/vc.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot       - kel * central
    d/dt(central_tgn)  <-  fm3       * kme * central - (cl_tgn  / vc_tgn ) * central_tgn
    d/dt(central_mmpn) <- (1 - fm3)  * kme * central - (cl_mmpn / vc_mmpn) * central_mmpn

    # 5. Bioavailability (oral) on the depot compartment.
    f(depot) <- exp(lfdepot)

    # 6. Observation variables. The 6-MP plasma concentration is NOT
    # measured in the source study; only the two RBC metabolites are
    # observed. With the apparent metabolite volumes fixed to 1 L the
    # RBC concentration is numerically equal to the compartment amount
    # in mg.
    Cc_tgn   <- central_tgn  / vc_tgn
    Cc_mmpn  <- central_mmpn / vc_mmpn

    Cc_tgn   ~ add(addSd_tgn)
    Cc_mmpn  ~ add(addSd_mmpn)
  })
}

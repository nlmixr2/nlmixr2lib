Diep_2026_donidalorsen <- function() {
  description <- "Two-compartment population PK and indirect-response PD model for the GalNAc3-conjugated antisense oligonucleotide donidalorsen targeting prekallikrein (PKK) mRNA, fit to pooled data from phase 1 to phase 3 studies in healthy volunteers and patients with hereditary angioedema (Diep 2026). First-order SC absorption with categorical covariates on ka (arm vs abdomen/thigh injection site; autoinjector vs vial drug presentation), allometric scaling of CL/F, Vc/F, Q/F, and Vp/F on total body weight with paper-estimated exponents, multiplicative disease-status effects on Vc/F and Q/F, full 5x5 omega block on PK random effects, and an indirect-response model with donidalorsen-driven inhibition of PKK production carrying multiplicative disease-status effects on baseline PKK and IC50."
  reference   <- "Diep JK, Liu M, Singh P, Dorow S, Cohn DM, Bordone L, Newman KB, Gao X. Population pharmacokinetic/pharmacodynamic modeling of donidalorsen, an antisense oligonucleotide in development for prophylaxis of hereditary angioedema. CPT Pharmacometrics Syst Pharmacol. 2026;15(2):e70206. doi:10.1002/psp4.70206"
  vignette    <- "Diep_2026_donidalorsen"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power exponents on CL/F (1.52), Vc/F (2.34), Q/F (1.79), and Vp/F (1.60); reference 70 kg = paper's stated reference subject (Diep 2026 Table 1 footnotes a-d and the simulation-reference description on page 4). Cohort median body weight was 78 kg with range 37-151.9 kg (Section 3 overview / Table S2).",
      source_name        = "WTKG"
    ),
    DIS_HAE = list(
      description        = "Hereditary angioedema patient indicator (1 = patient with HAE, 0 = healthy volunteer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = "Time-fixed per subject. Pooled across HAE-C1INH-Type1, HAE-C1INH-Type2, and HAE-nC1INH; reference category is the pooled healthy-volunteer cohort from NCT03263507 and ISIS 721744-CS9. Diep 2026 Table 1 footnotes b-c carry linear effects on Vc/F (+42.6%) and Q/F (-26.1%); Table 2 footnotes a-b carry linear effects on baseline PKK (-13.2%) and IC50 (+77.0%).",
      source_name        = "Disease status (patient with HAE indicator)"
    ),
    INJSITE_ARM = list(
      description        = "SC injection-site indicator: 1 = arm, 0 = abdomen or thigh",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen or thigh)",
      notes              = "Per-dose-record covariate. Diep 2026 Table 1 footnote e carries a linear effect on ka with multiplier 0.662 (= 1 - 0.338) when INJSITE_ARM = 1. The paper's reference category is 'abdomen or thigh' (paper text: 'injection into the arm vs. abdomen ... injection into the arm was predicted to have a 12% lower Cmax,ss compared with injection into the abdomen or thigh'); this is consistent with the canonical INJSITE_ARM definition whose reference includes abdomen.",
      source_name        = "Site of administration (arm vs abdomen/thigh)"
    ),
    DEVICE_AI = list(
      description        = "SC drug-presentation indicator: 1 = autoinjector (AI), 0 = vial and syringe",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (vial and syringe)",
      notes              = "Per-dose-record covariate. Diep 2026 Table 1 footnote e carries a linear effect on ka with multiplier 1.262 (= 1 + 0.262) when DEVICE_AI = 1. The autoinjector / vial-and-syringe contrast was characterized in the ISIS 721744-CS9 bioequivalence cohort.",
      source_name        = "Drug presentation (autoinjector vs vial and syringe)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 177L,
    n_studies      = 4L,
    n_healthy      = 101L,
    n_patients_hae = 76L,
    age_range      = "12 - 68 years",
    age_median     = "43 years",
    weight_range   = "37 - 151.9 kg",
    weight_median  = "78 kg",
    sex_female_pct = 55.4,
    race_ethnicity = c(White = 68.4),
    disease_state  = "Pooled adult and adolescent patients with hereditary angioedema (HAE-C1INH-Type1, HAE-C1INH-Type2, or HAE-nC1INH) together with healthy volunteers across the phase 1 to phase 3 program for donidalorsen.",
    dose_range     = "Subcutaneous donidalorsen at 20, 40, 60, or 80 mg every 4 weeks (Q4W) or 80 mg every 8 weeks (Q8W); also a single 80 mg dose in the ISIS 721744-CS9 vial-vs-autoinjector bioequivalence study.",
    regions        = "Global; phase 3 OASIS-HAE (NCT05139810) enrolled adult (>= 18 y) and adolescent (>= 12 to < 18 y) HAE patients.",
    studies        = "NCT03263507 (phase 1 multiple-ascending-dose 20/40/60/80 mg Q4W in healthy volunteers; n = 32, 3:1 active:placebo), NCT04030598 (phase 2 in HAE patients; 80 mg Q4W in Part A n = 20 [2:1] and Part B n = 3 open-label HAE-nC1INH), NCT05139810 / OASIS-HAE (phase 3 HAE; 80 mg Q4W cohort A n = 60 [3:1] and 80 mg Q8W cohort B n = 30 [3:1]), and ISIS 721744-CS9 (phase 1 single-dose vial-vs-autoinjector bioequivalence in healthy volunteers; n = 78).",
    notes          = "PK/PD analysis pooled n = 177 active-arm subjects (101 healthy volunteers + 76 HAE patients) after excluding placebo data, missing data, and PK samples collected at or beyond the onset of antidrug antibodies (3.1% of PK samples). Final dataset: 4242 plasma donidalorsen concentrations and 1159 plasma PKK concentrations (Diep 2026 Section 3 overview). PK quantification used a validated hybridization electrochemiluminescence assay with LLOQ = 0.1735 ng/mL; PKK quantification used a validated ELISA with LLOQ = 3.76 mg/L at minimum required dilution (lowest reportable value 0.47 mg/L). Postdose data below LLOQ (3.4% of PK, 1.5% of PD) were excluded per M1 methodology."
  )

  ini({
    # ---- Structural PK parameters (Diep 2026 Table 1 final-model column) ----
    # Reference subject: total body weight 70 kg (paper-stated reference per
    # page 4 / Figure 3 caption: "reference individual, who was defined as a
    # healthy volunteer with a body weight of 70 kg receiving donidalorsen
    # 80 mg Q4W subcutaneously into the abdomen or thigh, using a vial and
    # syringe as the drug presentation"), DIS_HAE = INJSITE_ARM = DEVICE_AI = 0.
    lka <- log(0.952); label("First-order SC absorption rate constant at reference covariates (ka, 1/h; INJSITE_ARM = DEVICE_AI = 0 reference)") # Diep 2026 Table 1 final-model ka = 0.952 1/h
    lcl <- log(12.8);  label("Apparent clearance at reference WT (CL/F, L/h)")              # Diep 2026 Table 1 final-model CL/F = 12.8 L/h
    lvc <- log(69.8);  label("Apparent central volume at reference WT (Vc/F, L)")           # Diep 2026 Table 1 final-model Vc/F = 69.8 L
    lq  <- log(2.58);  label("Apparent intercompartmental clearance at reference WT (Q/F, L/h)") # Diep 2026 Table 1 final-model Q/F = 2.58 L/h
    lvp <- log(1840);  label("Apparent peripheral volume at reference WT (Vp/F, L)")        # Diep 2026 Table 1 final-model Vp/F = 1840 L

    # ---- Covariate effects on PK (Diep 2026 Table 1 footnotes a-e) ----
    # Functional forms (paper footnotes, with WTKG -> WT and explicit indicators):
    #   CL/F = TVCL/F * (WT/70)^e_wt_cl
    #   Vc/F = TVVc/F * (WT/70)^e_wt_vc * (1 + e_dis_hae_vc * DIS_HAE)
    #   Q/F  = TVQ/F  * (WT/70)^e_wt_q  * (1 + e_dis_hae_q  * DIS_HAE)
    #   Vp/F = TVVp/F * (WT/70)^e_wt_vp
    #   ka   = TVka  * (1 + e_injsite_arm_ka * INJSITE_ARM) * (1 + e_device_ai_ka * DEVICE_AI)
    # The (1 + theta * COV) shape for categorical covariates is the Phoenix
    # NLME linear parameterisation; the Table 1 "Estimate" column reports
    # the theta value, and the page 4 prose / Table 1 footnotes spell out
    # the corresponding multipliers (e.g., 1.426 = 1 + 0.426 for HAE on Vc/F).
    e_wt_cl          <- 1.52;   label("Allometric exponent of WT on CL/F (unitless)")            # Diep 2026 Table 1 footnote a: CL/F = 12.8 x (WTKG/70)^1.52
    e_wt_vc          <- 2.34;   label("Allometric exponent of WT on Vc/F (unitless)")            # Diep 2026 Table 1 footnote b: Vc/F = 69.8 x (WTKG/70)^2.34 x 1.426 [HAE]
    e_wt_q           <- 1.79;   label("Allometric exponent of WT on Q/F (unitless)")             # Diep 2026 Table 1 footnote c: Q/F  = 2.58 x (WTKG/70)^1.79 x 0.739 [HAE]
    e_wt_vp          <- 1.60;   label("Allometric exponent of WT on Vp/F (unitless)")            # Diep 2026 Table 1 footnote d: Vp/F = 1840 x (WTKG/70)^1.60
    e_dis_hae_vc     <- 0.426;  label("Linear effect of HAE-patient indicator on Vc/F (fraction)") # Diep 2026 Table 1 Estimate; multiplier 1 + 0.426 = 1.426 for HAE
    e_dis_hae_q      <- -0.261; label("Linear effect of HAE-patient indicator on Q/F (fraction)")  # Diep 2026 Table 1 Estimate; multiplier 1 - 0.261 = 0.739 for HAE
    e_injsite_arm_ka <- -0.338; label("Linear effect of arm-injection-site indicator on ka (fraction)") # Diep 2026 Table 1 Estimate; multiplier 1 - 0.338 = 0.662 for arm
    e_device_ai_ka   <-  0.262; label("Linear effect of autoinjector-device indicator on ka (fraction)") # Diep 2026 Table 1 Estimate; multiplier 1 + 0.262 = 1.262 for AI

    # ---- Structural PD parameters (Diep 2026 Table 2 final-model column) ----
    # Indirect-response model with inhibition of PKK production:
    #   d/dt(effect) = kin * (1 - Cp * Imax / (IC50 + Cp)) - kout * effect
    #   kin = BL * kout (so effect = BL at the no-drug steady state; Diep 2026
    #   Section 3.2: "The model was parameterized with baseline PKK (BL)
    #   estimated as a parameter: BL = kin/kout"). The PD compartment is the
    #   canonical `effect`; the observation variable is aliased to `pkk`.
    # Functional forms (paper footnotes):
    #   BL   = TVBL   * (1 + e_dis_hae_bl   * DIS_HAE)
    #   IC50 = TVIC50 * (1 + e_dis_hae_ic50 * DIS_HAE)
    lbl   <- log(139);     label("Baseline plasma prekallikrein at reference DIS_HAE = 0 (BL, mg/L)")               # Diep 2026 Table 2 final-model BL = 139 mg/L
    lkout <- log(0.00266); label("First-order PKK loss rate constant (kout, 1/h)")                                  # Diep 2026 Table 2 final-model kout = 0.00266 1/h
    imax  <- 0.992;        label("Maximum fractional inhibition of PKK production by donidalorsen (Imax, unitless)") # Diep 2026 Table 2 final-model Imax = 0.992
    lic50 <- log(0.158);   label("Plasma donidalorsen concentration yielding half-maximum inhibition (IC50, ng/mL)") # Diep 2026 Table 2 final-model IC50 = 0.158 ng/mL

    # ---- Covariate effects on PD (Diep 2026 Table 2 footnotes a-b) ----
    e_dis_hae_bl   <- -0.132; label("Linear effect of HAE-patient indicator on BL (fraction)")   # Diep 2026 Table 2 Estimate; multiplier 1 - 0.132 = 0.868 for HAE
    e_dis_hae_ic50 <-  0.770; label("Linear effect of HAE-patient indicator on IC50 (fraction)") # Diep 2026 Table 2 Estimate; multiplier 1 + 0.770 = 1.770 for HAE

    # ---- IIV (Diep 2026 Tables 1 and 2; log-normal eta, omega^2 = log(CV^2 + 1)) ----
    # PK: full 5x5 omega block on (etalcl, etalvc, etalq, etalvp, etalka)
    # using the BSV%CV and pairwise correlations from Table 1.
    #   sigma_i = sqrt(log(CV_i^2 + 1))
    #   cov(i,j) = rho_ij * sigma_i * sigma_j
    # BSV%CV values from Table 1: CL = 19.6%, Vc = 59.0%, Q = 10.2%, Vp = 42.9%, ka = 42.8%.
    # Variances (omega^2):
    #   etalcl: log(1 + 0.196^2) = 0.037708
    #   etalvc: log(1 + 0.590^2) = 0.298499
    #   etalq:  log(1 + 0.102^2) = 0.010350
    #   etalvp: log(1 + 0.429^2) = 0.168994
    #   etalka: log(1 + 0.428^2) = 0.168258
    # Standard deviations (sigma):
    #   sd_cl = 0.194186, sd_vc = 0.546351, sd_q = 0.101737, sd_vp = 0.411089, sd_ka = 0.410193
    # Correlations from Table 1 and resulting covariances:
    #   rho(CL,Vc) = 0.753  -> cov = 0.753 * 0.194186 * 0.546351 = 0.079891
    #   rho(CL,Q)  = 0.682  -> cov = 0.682 * 0.194186 * 0.101737 = 0.013474
    #   rho(CL,Vp) = 0.438  -> cov = 0.438 * 0.194186 * 0.411089 = 0.034966
    #   rho(CL,ka) = 0.380  -> cov = 0.380 * 0.194186 * 0.410193 = 0.030266
    #   rho(Vc,Q)  = 0.741  -> cov = 0.741 * 0.546351 * 0.101737 = 0.041189
    #   rho(Vc,Vp) = 0.564  -> cov = 0.564 * 0.546351 * 0.411089 = 0.126676
    #   rho(Vc,ka) = 0.721  -> cov = 0.721 * 0.546351 * 0.410193 = 0.161589
    #   rho(Q,Vp)  = 0.921  -> cov = 0.921 * 0.101737 * 0.411089 = 0.038514
    #   rho(Q,ka)  = 0.439  -> cov = 0.439 * 0.101737 * 0.410193 = 0.018320
    #   rho(Vp,ka) = 0.303  -> cov = 0.303 * 0.411089 * 0.410193 = 0.051098
    # Block c(...) is the lower-triangle, row-major, of the 5x5 omega matrix
    # in the eta order (etalcl, etalvc, etalq, etalvp, etalka).
    etalcl + etalvc + etalq + etalvp + etalka ~ c(
      0.037708,
      0.079891, 0.298499,
      0.013474, 0.041189, 0.010350,
      0.034966, 0.126676, 0.038514, 0.168994,
      0.030266, 0.161589, 0.018320, 0.051098, 0.168258
    )

    # PD: independent etas (Table 2; no off-diagonal omega block reported).
    etalbl   ~ 0.064920  # Diep 2026 Table 2 BSV%CV BL  = 25.9% -> log(0.259^2 + 1)
    etalkout ~ 0.125749  # Diep 2026 Table 2 BSV%CV kout = 36.6% -> log(0.366^2 + 1)
    etalic50 ~ 0.525122  # Diep 2026 Table 2 BSV%CV IC50 = 83.1% -> log(0.831^2 + 1)

    # ---- Residual error ----
    # PK: Diep 2026 Table 1 reports sigma_logadd = 0.296 as the additive
    # component of the residual error on log-transformed plasma donidalorsen
    # concentrations. NONMEM/Phoenix "additive on log scale" maps to a
    # proportional residual in nlmixr2's linear space (see
    # references/naming-conventions.md NONMEM error-block table:
    # Y = LOG(F) + EPS(1) -> Cc ~ prop(propSd)).
    # PD: Diep 2026 Table 2 reports sigma_prop = 0.159 as the proportional
    # component on linear PKK concentrations.
    propSd  <- 0.296; label("Proportional residual error on plasma donidalorsen Cc (fraction)") # Diep 2026 Table 1 sigma_logadd = 0.296
    propSd_pkk <- 0.159; label("Proportional residual error on plasma prekallikrein pkk (fraction)") # Diep 2026 Table 2 sigma_prop  = 0.159
  })

  model({
    # ---- 1. Individual PK parameters ----
    # ka carries categorical covariates only (no allometry); Diep 2026 Table 1
    # footnote e: ka = 0.952 * 0.662^[arm] * 1.262^[autoinjector], encoded as
    # the (1 + theta * COV) Phoenix linear-effect form.
    ka <- exp(lka + etalka) *
      (1 + e_injsite_arm_ka * INJSITE_ARM) *
      (1 + e_device_ai_ka   * DEVICE_AI)

    # CL/F, Vc/F, Q/F, Vp/F carry allometric WT effects and (for Vc/F, Q/F)
    # a multiplicative HAE-patient effect (Diep 2026 Table 1 footnotes a-d).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (1 + e_dis_hae_vc * DIS_HAE)
    q  <- exp(lq  + etalq ) * (WT / 70)^e_wt_q  * (1 + e_dis_hae_q  * DIS_HAE)
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- 2. Individual PD parameters ----
    # kin is derived so that effect = bl is the no-drug steady state
    # (Diep 2026 Section 3.2: BL = kin / kout).
    bl   <- exp(lbl   + etalbl)   * (1 + e_dis_hae_bl   * DIS_HAE)
    kout <- exp(lkout + etalkout)
    ic50 <- exp(lic50 + etalic50) * (1 + e_dis_hae_ic50 * DIS_HAE)
    kin  <- bl * kout

    # ---- 3. ODE system ----
    # depot / central / peripheral1: amounts in mg (subcutaneous dose enters
    # depot in mg). Vc and Vp in L; central / vc has units mg/L = ug/mL.
    # Multiply by 1000 to express plasma donidalorsen in ng/mL so that the
    # paper's IC50 (0.158 ng/mL) shares units with Cc.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    Cc <- central / vc * 1000

    # ---- 4. PD: indirect response with inhibition of PKK production ----
    # Inhibition factor 1 - Cc * Imax / (IC50 + Cc) with Cc and IC50 in ng/mL
    # (Diep 2026 Section 3.2 inset equation). At Cc = 0 the factor is 1
    # (uninhibited basal production); at Cc >> IC50 the factor approaches
    # 1 - Imax = 0.008 (>99% inhibition).
    inh          <- 1 - Cc * imax / (ic50 + Cc)
    d/dt(effect) <- kin * inh - kout * effect
    effect(0)    <- bl

    # pkk is the paper-named alias of the canonical `effect` PD compartment
    # (plasma prekallikrein concentration in mg/L).
    pkk <- effect

    # ---- 5. Observations ----
    Cc  ~ prop(propSd)
    pkk ~ prop(propSd_pkk)
  })
}

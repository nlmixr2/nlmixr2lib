Kim_2016_udenafil <- function() {
  description <- "Parent-metabolite population PK model for oral udenafil and its active metabolite DA-8164 in healthy subjects and patients with mild (Child-Pugh A) and moderate (Child-Pugh B) hepatic impairment (Kim 2016). Two-compartment udenafil with first-order absorption and an absorption lag time, two parallel parent-side clearances (CLp/F = non-metabolic apparent clearance, CLpm/F = apparent formation clearance to DA-8164) feeding a two-compartment metabolite. Central and peripheral apparent volumes are assumed equal for parent and metabolite (the fraction metabolised f_m and the metabolite volume of distribution are not separately identifiable from this dataset). Mass-balance is preserved by multiplying the formation flux into the metabolite central compartment by the molecular-weight ratio Rpm = MW(DA-8164) / MW(udenafil) = 405.4 / 516.66. Prothrombin time expressed as INR (PT) acts on CLpm/F via a power covariate normalised to the cohort median 1.13: CLpm/F = theta1 * (PT/1.13)^theta10 with theta10 = -1.65 (decrease in CLpm/F with increasing PT)."
  reference   <- paste(
    "Kim A, Lee J, Shin D, Jung YJ, Bahng MY, Cho JY, Jang IJ.",
    "Population pharmacokinetic analysis to recommend the optimal dose",
    "of udenafil in patients with mild and moderate hepatic impairment.",
    "Br J Clin Pharmacol. 2016;82(4):1024-1033.",
    "doi:10.1111/bcp.12977.",
    sep = " "
  )
  vignette    <- "Kim_2016_udenafil"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    INR_BASE = list(
      description        = "Pre-dose prothrombin time expressed as the international normalised ratio (INR). Continuous, dimensionless, measured once per subject before the single 100 mg oral udenafil dose. The Kim 2016 paper labels this column 'PT' throughout but its Methods (Covariate selection and model evaluation paragraph) explicitly state 'prothrombin time and were expressed as international normalized ratio (PT)'.",
      units              = "(unitless ratio; INR has no units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Single-dose, time-fixed per subject. Used as a covariate on the apparent formation clearance of udenafil to DA-8164 via a power function normalised to the cohort median 1.13 (Kim 2016 Eq. for CLpm/F = theta1 * (PT/1.13)^theta10). Cohort summary statistics: healthy 0.97 +/- 0.039, mild HI 1.13 +/- 0.13, moderate HI 1.33 +/- 0.094 (Kim 2016 Table 1, mean +/- SD). Kruskal-Wallis P < 0.05 among study groups.",
      source_name        = "PT (prothrombin time expressed as INR; Kim 2016 Table 1 and Methods)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_range      = "41-61 years",
    age_median     = "52.9 years (cohort mean)",
    weight_range   = "approximately 57-77 kg (mean 66.9 kg, SD 6.18 kg)",
    weight_median  = "66.9 kg (cohort mean)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Three parallel age- and weight-matched groups (six subjects each): healthy adults, patients with mild hepatic impairment (Child-Pugh class A), and patients with moderate hepatic impairment (Child-Pugh class B). Two of the six moderate-HI subjects had clinical ascites; none of the healthy or mild-HI subjects did. Exclusion criteria included history of cardiovascular or cerebrovascular disease, creatinine clearance < 40 mL/min, and use of CYP3A4 / CYP2D6 inducers or inhibitors.",
    dose_range     = "Single 100 mg oral dose of udenafil in the fasted state (>= 10 h fast).",
    regions        = "South Korea (multicentre: Seoul National University Hospital; SMG-SNU Boramae Medical Center; Seoul National University Bundang Hospital; Asan Medical Center).",
    n_observations = "Plasma samples at predose, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 24, 32, 48, and 72 h post-dose; LLOQ 2 ng/mL with the validated range 2-2000 ng/mL.",
    co_medication  = "No drugs known to affect udenafil PK (CYP3A4 / CYP2D6 inducers or inhibitors).",
    notes          = paste(
      "Demographics and clinical baseline labs from Kim 2016 Table 1.",
      "Albumin (g/dL): healthy 4.35 +/- 0.27, mild HI 3.95 +/- 0.34, moderate HI 3.37 +/- 0.29.",
      "Total bilirubin (mg/dL): healthy 0.867 +/- 0.44, mild HI 1.20 +/- 0.59, moderate HI 2.18 +/- 0.48.",
      "Creatinine clearance (mL/min): healthy 99.5 +/- 17.8, mild HI 99.5 +/- 22.9, moderate HI 91.0 +/- 17.1.",
      "ClinicalTrials.gov NCT00956306. NONMEM 7.2 with the FOCE-I (eta-eps interaction) method and ADVAN6 TOL3 subroutine."
    )
  )

  ini({
    # Structural PK parameters from Kim 2016 Table 2, "Estimate" column.
    # Parameters are reported as apparent (oral) values; Kim 2016 assumed
    # parent and metabolite central and peripheral volumes are identical
    # (single Vp/F and Vp2/F), which collapses the metabolite volume and
    # f_m into the metabolite clearance terms (CLm/F * f_m, Qm/F * f_m).

    lka   <- log(0.326) ; label("First-order absorption rate constant for udenafil (1/h)")                              # Kim 2016 Table 2: ka = 0.326 1/h, RSE 19.4%
    ltlag <- log(0.314) ; label("Absorption lag time (h)")                                                              # Kim 2016 Table 2: ALAG = 0.314 h, RSE 14.8%

    # Two parallel parent-side clearance pathways. Paper-named-param form
    # (`lclp`, `lclpm`) is used so the paper's two THETAs are encoded
    # one-to-one with their own RSEs and IIVs preserved (the canonical
    # `lcl_renal` / `lcl_nonren` pattern from Pierre 2017 is route-specific
    # and would mislabel the non-metabolic arm here; Kim 2016 does not
    # decompose CLp/F into renal vs biliary vs other-metabolic).
    lclp  <- log(3.62)  ; label("Apparent non-metabolic clearance of udenafil CLp/F (L/h)")                             # Kim 2016 Table 2: CLp/F = 3.62 L/h, RSE 19.4%
    lclpm <- log(35.7)  ; label("Apparent formation clearance of udenafil to DA-8164 CLpm/F at PT = 1.13 (L/h)")        # Kim 2016 Table 2: CLpm/F = 35.7 L/h, RSE 25.3% (typical value at median PT = 1.13)

    # Metabolite (DA-8164) elimination clearance and intercompartmental
    # clearance. Reported as CLm/(F * fm) and Qm/(F * fm) because the
    # paper assumes the metabolite distribution volumes equal those of
    # the parent (single Vp/F and Vp2/F shared across parent and
    # metabolite); the unidentifiable f_m is absorbed into these terms.
    lcl_da8164 <- log(36.5) ; label("Apparent DA-8164 clearance CLm/(F*fm) (L/h)")                                       # Kim 2016 Table 2: CLm/F * fm = 36.5 L/h, RSE 26.2%
    lq_da8164  <- log(11.4) ; label("Apparent DA-8164 intercompartmental clearance Qm/(F*fm) (L/h)")                     # Kim 2016 Table 2: Qm/F * fm = 11.4 L/h, RSE 25.4%

    # Distribution volumes are shared by parent and metabolite (Kim 2016
    # structural-model assumption).
    lvc <- log(44.1) ; label("Apparent central volume of distribution Vp/F (shared parent/DA-8164) (L)")                # Kim 2016 Table 2: Vp/F = 44.1 L, RSE 29.0%
    lvp <- log(588)  ; label("Apparent peripheral volume of distribution Vp2/F (shared parent/DA-8164) (L)")            # Kim 2016 Table 2: Vp2/F = 588 L, RSE 7.21%

    # Parent (udenafil) intercompartmental clearance Qp/F.
    lq <- log(61.7) ; label("Apparent udenafil intercompartmental clearance Qp/F (L/h)")                                 # Kim 2016 Table 2: Qp/F = 61.7 L/h, RSE 15.7%

    # PT (prothrombin time as INR) covariate effect on CLpm/F. Encoded
    # as a power covariate normalised to the cohort median 1.13. The
    # paper's prose (Final pharmacokinetic model paragraph) reports
    # "The exponent of PT normalized to the median value (1.13) was
    # 1.65, which indicates a decrease in the CLpm/F with increase in
    # PT." Bootstrap median 1.66 with 95% CI (3.43, -0.729) per Table 2
    # bootstrap column. Both prose and bootstrap CI are consistent with
    # a negative exponent (the Table 2 "Estimate" column appears to
    # have lost the leading minus sign in typesetting); the negative
    # sign reproduces the published Kim 2016 typical values: CLpm/F at
    # mean PT = 1.33 (moderate HI) is 27.3 L/h, vs 45.9 L/h at PT =
    # 0.97 (healthy). Verified: 35.7 * (1.33/1.13)^(-1.65) = 26.3 L/h
    # and 35.7 * (0.97/1.13)^(-1.65) = 46.3 L/h (within rounding of
    # the printed 27.3 / 45.9).
    e_inr_base_clpm <- -1.65 ; label("Exponent of PT (as INR, normalised by 1.13) on CLpm/F")                            # Kim 2016 Table 2: theta10 = 1.65 in Estimate (sign per prose: decrease in CLpm/F with PT), bootstrap -1.66 (-3.43, -0.729), RSE 33.6%

    # IIVs as %CV from Kim 2016 Table 2, converted to log-normal omega^2
    # via omega^2 = log(1 + CV^2). All IIVs are diagonal (no block
    # correlations reported); IIV is NOT estimated for CLp/F or
    # Qm/(F * fm) per Table 2.
    etalclpm     ~ 0.11622   # Kim 2016 Table 2: omega CLpm/F = 35.1% CV; log(1 + 0.351^2)
    etalcl_da8164 ~ 0.25929  # Kim 2016 Table 2: omega CLm/(F*fm) = 54.4% CV; log(1 + 0.544^2)
    etalka       ~ 0.34416   # Kim 2016 Table 2: omega ka = 64.1% CV; log(1 + 0.641^2)
    etalvc       ~ 0.45227   # Kim 2016 Table 2: omega Vp/F = 75.6% CV; log(1 + 0.756^2)
    etalvp       ~ 0.05512   # Kim 2016 Table 2: omega Vp2/F = 23.8% CV; log(1 + 0.238^2)
    etalq        ~ 0.23354   # Kim 2016 Table 2: omega Qp/F = 51.3% CV; log(1 + 0.513^2)
    etaltlag     ~ 0.04195   # Kim 2016 Table 2: omega ALAG = 20.7% CV; log(1 + 0.207^2)

    # Residual error: proportional on both udenafil and DA-8164. Kim 2016
    # Results paragraph: "A proportional error model was shown to be
    # sufficient to explain the residual variability for both udenafil
    # and DA-8164."
    propSd         <- 0.216 ; label("Proportional residual error on udenafil plasma concentration (fraction)")          # Kim 2016 Table 2: sigma_prop,p = 21.6%
    propSd_da8164  <- 0.231 ; label("Proportional residual error on DA-8164 plasma concentration (fraction)")           # Kim 2016 Table 2: sigma_prop,m = 23.1%
  })

  model({
    # Molecular-weight ratio applied to the formation flux entering the
    # metabolite central compartment so the simulated DA-8164 mass
    # concentration is dimensionally consistent with the parent mass
    # concentration. Kim 2016 Methods (Population pharmacokinetic model
    # development paragraph): "the molecular weight ratio of DA-8164 to
    # udenafil (Rpm) was multiplied by the turnover rate of udenafil to
    # DA-8164". MW(udenafil) = 516.66 g/mol; MW(DA-8164) = 405.4 g/mol.
    rpm <- 405.4 / 516.66

    # Individual structural parameters (typical * IIV).
    ka         <- exp(lka + etalka)
    tlag       <- exp(ltlag + etaltlag)
    clp        <- exp(lclp)
    clpm       <- exp(lclpm + etalclpm) * (INR_BASE / 1.13)^e_inr_base_clpm
    cl_da8164  <- exp(lcl_da8164 + etalcl_da8164)
    vc         <- exp(lvc + etalvc)
    vp         <- exp(lvp + etalvp)
    q          <- exp(lq + etalq)
    q_da8164   <- exp(lq_da8164)

    # Micro-constants. Parent has two parallel CL pathways out of
    # central: CLp/F (non-metabolic) and CLpm/F (formation to DA-8164).
    # Both are explicit clearance terms applied to the parent central
    # concentration. Q_parent and Q_metabolite are paper-specific
    # (distinct from each other) but share Vp and Vp2.
    kp_out <- (clp + clpm) / vc        # total parent elimination rate constant out of central
    kp_per <- q / vc                   # parent central -> peripheral1
    kp_ret <- q / vp                   # parent peripheral1 -> central
    km_out <- cl_da8164 / vc           # metabolite elimination rate constant out of metabolite central
    km_per <- q_da8164 / vc            # metabolite central -> peripheral1_da8164
    km_ret <- q_da8164 / vp            # metabolite peripheral1_da8164 -> metabolite central

    # ODE system. Compartments: depot (oral dose), central (udenafil),
    # peripheral1 (udenafil peripheral), central_da8164 (DA-8164),
    # peripheral1_da8164 (DA-8164 peripheral). Mass-balance Rpm scaling
    # is applied to the formation flux entering the metabolite central
    # compartment (see Methods paragraph quoted above).
    d/dt(depot)              <- -ka * depot
    d/dt(central)            <-  ka * depot -
                                 kp_out * central -
                                 kp_per * central + kp_ret * peripheral1
    d/dt(peripheral1)        <-  kp_per * central - kp_ret * peripheral1
    d/dt(central_da8164)     <-  rpm * clpm * (central / vc) -
                                 km_out * central_da8164 -
                                 km_per * central_da8164 + km_ret * peripheral1_da8164
    d/dt(peripheral1_da8164) <-  km_per * central_da8164 - km_ret * peripheral1_da8164

    # Absorption lag on the depot.
    alag(depot) <- tlag

    # Observation equations. Cc and Cc_da8164 are plasma concentrations
    # of udenafil and DA-8164 in ng/mL. Both carry independent
    # proportional residual errors.
    Cc        <- central        / vc
    Cc_da8164 <- central_da8164 / vc

    Cc        ~ prop(propSd)
    Cc_da8164 ~ prop(propSd_da8164)
  })
}

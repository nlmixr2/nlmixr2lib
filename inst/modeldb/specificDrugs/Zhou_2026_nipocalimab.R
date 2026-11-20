Zhou_2026_nipocalimab <- function() {
  description <- "Two-compartment quasi-steady-state TMDD model for nipocalimab with FcRn occupancy and an indirect-response total IgG model, fit to pooled healthy adults and pregnant participants at risk of early-onset severe hemolytic disease of the fetus and newborn (HDFN), with time-varying body weight on CL and Vc (Zhou 2026)"
  reference <- "Zhou J, Xiong Y, Lam E, Ling LE, Saeed-Khawaja S, Komatsu Y, Vermeulen A, Samtani MN, Leu JH. Nipocalimab Phase 3 Dose Selection for Severe Hemolytic Disease of the Fetus and Newborn. J Clin Pharmacol. 2026;66:e70259. doi:10.1002/jcph.70259"
  vignette <- "Zhou_2026_nipocalimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (gestational weight gain during pregnancy). Allometric scaling (WT/75)^0.75 on CL and (WT/75)^1 on Vc (Supplement Equation 3); Q and Vp are not weight-scaled. Body weight is also how the mg/kg dose is converted to the mg amount in the event table.",
      source_name = "WT"
    ),
    ELISA = list(
      description = "Bioanalytical assay indicator: 1 = ELISA ('old assay'), 0 = Meso Scale Discovery electrochemiluminescence immunoassay ('new assay')",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (electrochemiluminescence)",
      notes = "Switches the additive PK residual SD only (0.424 nmol/L ELISA vs 0.0306 nmol/L ECLIA, Supplement PK parameter table). ELISA was used in three healthy-volunteer studies; the electrochemiluminescence assay was used in healthy-volunteer study NCT04848558 and in the UNITY HDFN study, so set ELISA = 0 to simulate pregnant participants.",
      source_name = "ADD1 / ADD2 RUV PK"
    )
  )

  compartmentData <- list(
    central = list(analyte = "nipocalimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "nipocalimab", units = "mg", specimen = "tissue", verified = TRUE),
    total_target = list(analyte = "FcRn", units = "nmol/L", specimen = "not applicable", verified = TRUE),
    total_igg = list(analyte = "IgG", units = "umol/L", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 144L,
    n_studies = 5L,
    age_range = "18-64 years",
    age_median = "35 years",
    weight_range = "50.3-110 kg",
    weight_median = "71.0 kg (PK/FcRn data set); 70.5 kg (IgG data set)",
    sex_female_pct = 66.0,
    race_ethnicity = c(White = 73.6, Asian = 18.8, Black = 4.2, Other = 3.5),
    disease_state = "Healthy adult volunteers pooled with pregnant participants at high risk of recurrent early-onset severe hemolytic disease of the fetus and newborn (EOS-HDFN; prior pregnancy with severe fetal anemia, hydrops or stillbirth, anti-D titer >= 32 or anti-Kell titer >= 4)",
    dose_range = "Healthy volunteers: single IV 0.3-60 mg/kg or 300 / 1200 mg, or 15 / 30 mg/kg IV weekly x 4. Pregnancy (UNITY): 30 or 45 mg/kg IV weekly from gestational week 14 to 35, on baseline or current body weight",
    regions = "Multinational (UNITY Phase 2, NCT03842189) plus four Phase 1 healthy-volunteer studies including NCT04848558",
    n_pregnancies = 14L,
    n_healthy_pk = 114L,
    n_healthy_igg = 130L,
    ada_positive_pct = 36.8,
    notes = "Zhou 2026 Table 1. The PK/FcRn-occupancy data set had 128 participants (114 healthy, 14 pregnancies from 13 UNITY participants; 1513 PK and 1209 FcRn observations); the IgG data set had 144 (130 healthy, 14 pregnancies; 2110 observations). Demographics above are the IgG data set's; the PK/FcRn data set had 64.8% female and race 75.8% White, 16.4% Asian, 4.6% Black, 3.1% Other."
  )

  ini({
    # ---- PK and FcRn binding (Supplement 'Parameter estimates for the final nipocalimab HDFN PK model'; reference weight 75 kg) ----
    lcl <- log(0.687); label("Linear clearance for a 75 kg adult (L/day)") # Supplement PK table: CL = 0.687 L/day
    lvc <- log(3.36); label("Central volume of distribution for a 75 kg adult (L)") # Supplement PK table: Vc = 3.36 L
    lq <- log(0.238); label("Intercompartmental clearance (L/day)") # Supplement PK table: Q = 0.238 L/day
    lvp <- log(0.548); label("Peripheral volume of distribution (L)") # Supplement PK table: Vp = 0.548 L
    lFcRn0 <- log(135); label("Total FcRn concentration at baseline, R_total,0 (nmol/L)") # Supplement PK table: R_total,0 = 135 nmol/L
    lkss <- log(43.9); label("Quasi-steady-state constant of the nipocalimab-FcRn complex (nmol/L)") # Supplement PK table: k_ss = 43.9 nmol/L (footnote b: 6.23 ug/mL)
    lkint <- log(65.6); label("Internalization rate constant of the nipocalimab-FcRn complex (1/day)") # Supplement PK table: k_int = 65.6 /day
    lkdeg <- fixed(log(1.3)); label("Degradation rate constant of free FcRn (1/day)") # Supplement PK table: k_deg = 1.3 /day, RSE column 'FIX'
    FRmax <- 0.948; label("Maximal fraction of FcRn that nipocalimab can bind, FR_R (fraction)") # Supplement PK table: FR_R = 0.948
    e_wt_cl <- fixed(0.75); label("Allometric exponent of time-varying body weight on CL (unitless)") # Supplement PK table and Equation 3: 0.75, 'FIX'
    e_wt_vc <- fixed(1); label("Allometric exponent of time-varying body weight on Vc (unitless)") # Supplement PK table and Equation 3: 1, 'FIX'

    # ---- Total IgG indirect response (Supplement 'Parameter estimates for the final nipocalimab HDFN PK-IgG model') ----
    lIgG0 <- log(75.2); label("Baseline total serum IgG (umol/L)") # Supplement IgG table: IgG baseline = 75.2 umol/L (footnote b: 11.3 g/L)
    lkdeg_IgG <- log(0.218); label("IgG degradation rate constant without FcRn recycling (1/day)") # Supplement IgG table: k_deg_IgG = 0.218 /day
    lIgK <- log(4.90); label("FcRn-mediated IgG recycling parameter IgK (unitless)") # Supplement IgG table: IgK = 4.90

    # ---- IIV (log-normal; omega^2 = log(CV^2 + 1)) ----
    etalcl ~ 0.041560 # Supplement PK table: IIV CL 20.6 %CV
    etalvc ~ 0.021958 # Supplement PK table: IIV Vc 14.9 %CV
    etalFcRn0 ~ 0.044370 # Supplement PK table: IIV R_total,0 21.3 %CV
    etalIgG0 ~ 0.050245 # Supplement IgG table: IIV IgG baseline 22.7 %CV
    etalkdeg_IgG ~ 0.086178 # Supplement IgG table: IIV k_deg_IgG 30.0 %CV

    # ---- Residual error ----
    addSdELISA <- 0.424; label("Additive PK residual SD, ELISA assay (nmol/L; converted to ug/mL in model)") # Supplement PK table: ADD1 RUV PK, old assay = 0.424 nmol/L
    addSdECLIA <- 0.0306; label("Additive PK residual SD, electrochemiluminescence assay (nmol/L; converted to ug/mL in model)") # Supplement PK table: ADD2 RUV PK, new assay = 0.0306 nmol/L
    propSd <- 0.141; label("Proportional PK residual SD (fraction)") # Supplement PK table: PROP RUV IV PK = 14.1 %
    addSd_pctUnoccupiedFcRn <- 4.95; label("Additive residual SD for percent unoccupied FcRn (percentage points)") # Supplement PK table: ADD RUV RO = 4.95 %
    propSd_pctUnoccupiedFcRn <- 0.241; label("Proportional residual SD for percent unoccupied FcRn (fraction)") # Supplement PK table: PROP RUV RO = 24.1 %
    propSd_IgG_obs <- 0.079; label("Proportional residual SD for total serum IgG (fraction)") # Supplement IgG table: Prop error = 7.90 %
    addSd_IgG_obs <- 4.13; label("Additive residual SD for total serum IgG (umol/L)") # Supplement IgG table: Additive error = 4.13 umol/L
  })

  model({
    # Nipocalimab molecular weight 142,000 g/mol (Supplement PK table footnote b).
    # 1 ug/mL = 1e6 / 142000 nmol/L, so nmol/L x 0.142 = ug/mL.
    mw <- 142000
    nM_per_ugmL <- 1e6 / mw

    # ---- Individual PK parameters (Equation 3; WT is time-varying) ----
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

    # ---- FcRn parameters, kept in nmol/L as estimated ----
    FcRn0 <- exp(lFcRn0 + etalFcRn0)
    kss <- exp(lkss)
    kint <- exp(lkint)
    kdeg <- exp(lkdeg)
    # FcRn synthesis from the pre-dose steady state of the bindable pool,
    # R_total(time 0) = R_total,0 x FR_R (Equation 2): dR_total/dt = 0 at t = 0.
    ksyn <- kdeg * FcRn0 * FRmax

    # ---- IgG parameters (Equations 4 and 5) ----
    IgG_BL <- exp(lIgG0 + etalIgG0)
    kdeg_IgG <- exp(lkdeg_IgG + etalkdeg_IgG)
    IgK <- exp(lIgK)
    krec_IgG <- kdeg_IgG * IgK / (1 + IgK)
    ksyn_IgG <- (kdeg_IgG - krec_IgG) * IgG_BL

    total_target(0) <- FcRn0 * FRmax
    total_igg(0) <- IgG_BL

    # ---- Quasi-steady-state binding in serum (Equation 2), in nmol/L ----
    # central holds the total (free + FcRn-bound) nipocalimab amount in mg.
    ctot_nM <- central / vc * nM_per_ugmL
    disc <- ctot_nM - kss - total_target
    cfree_nM <- 0.5 * (disc + sqrt(disc * disc + 4 * ctot_nM * kss))
    complex <- total_target * cfree_nM / (kss + cfree_nM)
    cfree <- cfree_nM / nM_per_ugmL

    # ---- ODEs (Equations 1 and 4) ----
    # The internalised complex amount kint x Vc x RC is converted from nmol back to mg.
    d/dt(central) <- -cl * cfree - kint * vc * complex / nM_per_ugmL - q * cfree + q / vp * peripheral1
    d/dt(peripheral1) <- q * cfree - q / vp * peripheral1
    d/dt(total_target) <- ksyn - kdeg * total_target - (kint - kdeg) * complex
    # F_free: free bindable FcRn relative to the bindable pool at time 0, so it
    # reaches 0 at full blockade, as Equation 5's minimal-IgG derivation requires.
    f_free <- (total_target - complex) / (FcRn0 * FRmax)
    d/dt(total_igg) <- ksyn_IgG - (kdeg_IgG - krec_IgG * f_free) * total_igg

    # ---- Observations ----
    # Total nipocalimab concentration (ug/mL).
    Cc <- central / vc
    # Percent unoccupied FcRn (Equation 2): free FcRn, including the
    # (1 - FR_R) fraction nipocalimab cannot bind, relative to R_total,0.
    Rfree <- total_target - complex + FcRn0 * (1 - FRmax)
    pctUnoccupiedFcRn <- Rfree / FcRn0 * 100
    # Total serum IgG (umol/L).
    IgG_obs <- total_igg

    addSd <- (addSdELISA * ELISA + addSdECLIA * (1 - ELISA)) / nM_per_ugmL
    Cc ~ add(addSd) + prop(propSd)
    pctUnoccupiedFcRn ~ add(addSd_pctUnoccupiedFcRn) + prop(propSd_pctUnoccupiedFcRn)
    IgG_obs ~ add(addSd_IgG_obs) + prop(propSd_IgG_obs)
  })
}

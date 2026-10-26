Stachanow_2023_hydrocortisone <- function() {
  description <- "Pediatric PK/PD model of oral hydrocortisone (Alkindi granules) in young children and infants with congenital adrenal hyperplasia: total plasma cortisol from a two-compartment model with saturable (Michaelis-Menten) absorption, saturable corticosteroid-binding-globulin and linear albumin plasma binding acting on unbound cortisol, and allometric weight scaling (Stachanow 2022) inhibits the synthesis of 17-alpha-hydroxyprogesterone (17-OHP) measured in venous dried blood spots through an indirect-response model whose baseline is the subject's observed pre-dose DBS 17-OHP (Stachanow 2023)."
  reference <- paste(
    "Stachanow V, Neumann U, Blankenstein O, Alder-Baerens N, Bindellini D, Hindmarsh P, Ross RJ, Whitaker MJ, Melin J, Huisinga W, Michelet R, Kloft C.",
    "Model-Informed Target Morning 17alpha-Hydroxyprogesterone Concentrations in Dried Blood Spots for Pediatric Congenital Adrenal Hyperplasia Patients.",
    "Pharmaceuticals (Basel). 2023;16(3):464. doi:10.3390/ph16030464.",
    "Cortisol PK layer from Stachanow V, Neumann U, Blankenstein O, Bindellini D, Melin J, Ross R, Whitaker MJ, Huisinga W, Michelet R, Kloft C.",
    "Exploring Dried Blood Spot Cortisol Concentrations as an Alternative for Monitoring Pediatric Adrenal Insufficiency Patients: A Model-Based Analysis.",
    "Front Pharmacol. 2022;13:819590. doi:10.3389/fphar.2022.819590.",
    sep = " "
  )
  vignette <- "Stachanow_2023_hydrocortisone"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of the cortisol PK with reference 70 kg: exponent 0.75 on CL and Q, 1 on Vc and Vp (Stachanow 2022 Methods and Supplementary control stream, 'TVCL = THETA(1)*((BW/70)**0.75)', 'TVV1 = THETA(2)*(BW/70)'). PK/PD cohort median 15 kg, range 7-21 kg (Stachanow 2023 Table S1).",
      source_name = "BW"
    ),
    TRANSCORTIN = list(
      description = "Plasma corticosteroid-binding globulin (CBG) concentration",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Sets the capacity of the saturable CBG binding of cortisol with one binding site per CBG molecule: Bmax (nmol/L) = CBG (ug/mL = mg/L) * 1000 / 52 (MW 52 kDa) (Stachanow 2022 Supplementary control stream, 'CBGmol = CBG2*1000/52', 'BMAX = CBGmol*BS', BS = 1 FIX). Measured CBG was used for the pediatric data; the authors substituted 22.4 ug/mL when CBG was not measured (control stream THETA(15) = 22.4 FIX).",
      source_name = "CBG"
    ),
    OHP17_DBS_BL = list(
      description = "Observed pre-dose (baseline) 17-alpha-hydroxyprogesterone concentration in a dried blood spot",
      units = "nmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject baseline that initialises the 17-OHP turnover state and sets ksyn = baseline * kdeg (Stachanow 2023 Supplementary NONMEM code, 'BASEDV = ODV', 'OHPBASE = BASEDV*EXP(ETA(3)*THETA(3))', 'KSYN = OHPBASE*KDEG'). Baselines below the LLOQ of 1.3 nmol/L were entered as 0.65 nmol/L (LLOQ / 2; 2 of 18 subjects). For the Step-C simulations the authors used the cohort median of 14.4 nmol/L.",
      source_name = "ODV (baseline record)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "hydrocortisone", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cortisol", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cortisol", units = "nmol", specimen = "tissue", verified = TRUE),
    effect = list(analyte = "17-alpha-hydroxyprogesterone", units = "nmol/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "4 months - 5 years (cohorts: young children 2-6 years, infants 28 days - 2 years)",
    age_median = "2.5 years",
    weight_range = "7-21 kg",
    weight_median = "15 kg",
    disease_state = "Pediatric adrenal insufficiency, predominantly congenital adrenal hyperplasia (the phase 3 trial enrolled 23 CAH patients and 1 patient with hypopituitarism; the PK/PD analysis used its 12 young children and 6 infants).",
    dose_range = "Single individualised morning oral dose of hydrocortisone granules (Alkindi), 1-4 mg (Stachanow 2022 Methods).",
    regions = "Germany (Charite Berlin)",
    n_observations = "88 venous DBS 17-OHP concentrations (22 below the LLOQ of 1.3 nmol/L, handled with the M3 method) and 88 plasma cortisol concentrations",
    notes = "Stachanow 2023 Section 2.1.1 and Table S1; phase 3 Alkindi trial NCT02720952. Sex and race are not reported. The PD was fitted sequentially: individual cortisol PK parameters from the Stachanow 2022 model were data items. That PK model was fitted to 30 dexamethasone-suppressed healthy adults plus 24 pediatric patients (including 6 neonates who are not part of the PD cohort)."
  )

  ini({
    # ---- Cortisol PK (Stachanow 2022 Table 1, SIR medians; allometric reference 70 kg) ----
    lcl <- log(400); label("Cortisol clearance of unbound drug at 70 kg (L/h)") # Stachanow 2022 Table 1 'CL [L/h]' = 400
    lvc <- log(10.6); label("Cortisol central volume at 70 kg (L)") # Stachanow 2022 Table 1 'Vc [L]' = 10.6
    lq <- log(160); label("Cortisol intercompartmental clearance at 70 kg (L/h)") # Stachanow 2022 Table 1 'Q [L/h]' = 160
    lvp <- log(124); label("Cortisol peripheral volume at 70 kg (L)") # Stachanow 2022 Table 1 'Vp [L]' = 124
    ltmax_abs <- log(21388); label("Maximum saturable absorption rate (nmol/h)") # Stachanow 2022 Table 1 'Vmax [nmol/h]' = 21,388
    lkt_abs <- fixed(log(4810)); label("Depot amount giving half-maximal absorption rate (nmol)") # Stachanow 2022 Table 1 'Km [nmol]' = 4,810 (held at the Michelet 2020 estimate)
    lfdepot <- fixed(log(1)); label("Oral bioavailability (fraction)") # Stachanow 2022 Table 1 'F [-]' = 1
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of weight on CL and Q (unitless)") # Stachanow 2022 Methods and control stream '(BW/70)**0.75'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of weight on Vc and Vp (unitless)") # Stachanow 2022 Methods and control stream '(BW/70)'

    # Plasma protein binding (Melin 2017 binding model, carried into Stachanow 2022)
    lkcbg <- fixed(log(9.71)); label("Equilibrium dissociation constant of cortisol for CBG (nmol/L)") # Stachanow 2022 Table 1 'Kd [nM]' = 9.71
    nsalb <- fixed(4.15); label("Ratio of albumin-bound to unbound cortisol (unitless)") # Stachanow 2022 Table 1 'NSAlb [-]' = 4.15
    lbs_cbg <- fixed(log(1)); label("Cortisol binding sites per CBG molecule (unitless)") # Stachanow 2022 Supplementary control stream THETA(6) 'BS [-]' = 1 (from Melin 2017)

    # Pre-dose plasma cortisol in children (used as the central initial condition)
    lrbase_cortisol <- log(9.41); label("Pediatric pre-dose plasma cortisol (nmol/L)") # Stachanow 2022 Table 1 'BASEchild, pla [nM]' = 9.41

    # ---- 17-OHP indirect response (Stachanow 2023 Table 1) ----
    lkdeg <- log(1.22); label("First-order degradation rate constant of DBS 17-OHP (1/h)") # Stachanow 2023 Table 1 'kdeg [1/h]' = 1.22
    lic50 <- log(21.0); label("Total plasma cortisol giving half-maximal inhibition of 17-OHP synthesis (nmol/L)") # Stachanow 2023 Table 1 'IC50 [nmol/L]' = 21.0
    imax <- fixed(1); label("Maximum fractional inhibition of 17-OHP synthesis (fraction)") # Stachanow 2023 Table 1 'Imax [-]' = 1 (pre-defined)
    hill <- fixed(1); label("Hill coefficient of the inhibition (unitless)") # Stachanow 2023 Table 1 'Hill [-]' = 1 (estimates were close to 1)

    # ---- IIV (omega^2 = log(CV^2 + 1); exact conversion confirmed by 131.1 %CV for omega^2 = 1) ----
    etalcl ~ 0.064442 # Stachanow 2022 Table 1 'omega CL, %CV' = 25.8
    etaltmax_abs ~ 0.195752 # Stachanow 2022 Table 1 'omega Vmax, %CV' = 46.5
    etalkt_abs ~ 0.270217 # Stachanow 2022 Table 1 'omega Km, %CV' = 55.7
    etalfdepot ~ 0.122502 # Stachanow 2022 Table 1 'omega F, %CV' = 36.1
    etalbs_cbg ~ fixed(0.049) # Stachanow 2022 Supplementary control stream $OMEGA 6 'BS' = 0.049 (from Melin 2017)
    etalrbase_cortisol ~ fixed(0.020736) # Stachanow 2022 control stream 'IBASE = OBASE*EXP(ETA(13)*THETA(14))' with omega 1 and THETA(14) = Table 1 sigma exp 0.144; 0.144^2
    etalkdeg ~ fixed(0.0025) # Stachanow 2023 Supplementary $OMEGA 'KDEG' = 0.0025; Table 1 'IIV on kdeg, %CV' = 5.0
    etalic50 ~ 0.733137 # Stachanow 2023 Table 1 'IIV on IC50, %CV' = 104
    etalrbase ~ fixed(0.145161) # Stachanow 2023 Supplementary 'OHPBASE = BASEDV*EXP(ETA(3)*THETA(3))' with omega 1 and THETA(3) = RUV SD 0.381; 0.381^2

    # ---- Residual error (additive on the log scale) ----
    expSd <- 0.144; label("Log-scale residual SD of total plasma cortisol (SD)") # Stachanow 2022 Table 1 'sigma exp [CV%]' = 14.4
    expSd_ohp17_dbs <- 0.381; label("Log-scale residual SD of DBS 17-OHP (SD)") # Stachanow 2023 Table 1 'RUV [%CV]' = 38.1
  })

  model({
    # ---- Cortisol PK (Stachanow 2022 Supplementary equations and control stream) ----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp
    tmax_abs <- exp(ltmax_abs + etaltmax_abs)
    kt_abs <- exp(lkt_abs + etalkt_abs)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # CBG capacity: CBG (mg/L) / 52 kDa -> umol/L, x 1000 -> nmol/L; amounts
    # in nmol throughout, so the binding constants are scaled by vc.
    bmax <- TRANSCORTIN * 1000 / 52 * exp(lbs_cbg + etalbs_cbg)
    kcbg <- exp(lkcbg)
    kcbg_amt <- kcbg * vc
    bmax_amt <- bmax * vc

    # Unbound amount in central: positive root of the quadratic for saturable
    # CBG binding plus linear albumin binding (bound_alb = nsalb * unbound).
    bterm <- central - kcbg_amt * (1 + nsalb) - bmax_amt
    central_u <- (bterm + sqrt(bterm^2 + 4 * kcbg_amt * central * (1 + nsalb))) / (2 * (1 + nsalb))

    rbase_cortisol <- exp(lrbase_cortisol + etalrbase_cortisol)

    # ---- 17-OHP indirect response (Stachanow 2023 Supplementary) ----
    kdeg <- exp(lkdeg + etalkdeg)
    ic50 <- exp(lic50 + etalic50)
    rbase <- OHP17_DBS_BL * exp(etalrbase)
    ksyn <- rbase * kdeg

    d/dt(depot) <- -tmax_abs * depot / (kt_abs + depot)
    d/dt(central) <- tmax_abs * depot / (kt_abs + depot) - kel * central_u - k12 * central_u + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central_u - k21 * peripheral1
    central(0) <- rbase_cortisol * vc
    f(depot) <- fdepot

    # Total plasma cortisol (nmol/L) drives the inhibition.
    Cc <- central / vc
    inh <- 1 - imax * Cc^hill / (ic50^hill + Cc^hill)
    d/dt(effect) <- ksyn * inh - kdeg * effect
    effect(0) <- rbase

    # 17-OHP in dried blood spots (nmol/L), the paper-named alias of `effect`.
    ohp17_dbs <- effect

    Cc ~ lnorm(expSd)
    ohp17_dbs ~ lnorm(expSd_ohp17_dbs)
  })
}

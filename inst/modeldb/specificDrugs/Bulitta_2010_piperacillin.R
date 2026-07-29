Bulitta_2010_piperacillin <- function() {
  description <- "Three-compartment population PK model for piperacillin in healthy adult volunteers after a single intravenous infusion, with first-order non-renal clearance and parallel first-order plus mixed-order (Michaelis-Menten) renal elimination, allometrically scaled to 70 kg; a urine compartment accumulates the renally excreted amount (Bulitta 2010 Model 3, final model, NONMEM estimates)"
  reference <- paste(
    "Bulitta JB, Kinzig M, Jakob V, Holzgrabe U, Sorgel F, Holford NHG.",
    "Nonlinear pharmacokinetics of piperacillin in healthy volunteers --",
    "implications for optimal dosage regimens.",
    "Br J Clin Pharmacol. 2010;70(5):682-693.",
    "doi:10.1111/j.1365-2125.2010.03750.x"
  )
  vignette <- "Bulitta_2010_piperacillin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight; used for allometric scaling with reference weight 70 kg. Exponent 0.75 on clearances (CL_NR, CL_R, CL_ic_shallow, CL_ic_deep) and on V_max (Table 3 footnote *); exponent 1.0 on volumes V_1, V_2, V_3 (Table 3 footnote +). Source paper cohort: median 77.5 kg, range 67-85 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 4L,
    n_studies       = 1L,
    n_occasions     = 5L,
    age_range       = "22 to 24 years (Table 1)",
    age_median      = "23 years (Table 1)",
    weight_range    = "67 to 85 kg (Table 1)",
    weight_median   = "77.5 kg (Table 1)",
    height_range    = "164 to 178 cm (Table 1)",
    sex_female_pct  = 50,
    race_ethnicity  = c(White = 100),
    disease_state   = "Healthy Caucasian adult volunteers",
    dose_range      = "Single 4 g intravenous 5-min infusion of piperacillin in each of five study periods (occasions) at days 1, 3, 10, 24, 52",
    regions         = "Single-centre Germany (University Hospital of Essen)",
    notes           = "4 subjects (2 female, 2 male) each studied on five replicate occasions; the analysis treated the dataset as 20 effective profiles. Baseline demographics from Table 1. Sampling: 21 plasma and 9 urine windows per occasion after each 5-min infusion. Piperacillin given alone (no tazobactam). Renal function inferred normal (healthy volunteers); creatinine clearance not tabulated."
  )

  ini({
    # Structural parameters: Bulitta 2010 Table 3 (NONMEM V FOCE+I, used for
    # simulation in the paper). Model 3 has first-order non-renal elimination
    # (cl_nonren) and parallel first-order (cl_renal) plus mixed-order
    # Michaelis-Menten (vmax / km) renal elimination; three-compartment
    # disposition with central, shallow peripheral, and deep peripheral.
    # All disposition parameters are scaled allometrically to a 70 kg
    # reference subject (Methods, "Size"; Table 3 footnotes).
    # Concentration Cc = central / vc gives mg/L; vmax * Cc / (km + Cc),
    # cl_renal * Cc, and cl_nonren * Cc are all mass-rate fluxes in mg/h.

    lcl_nonren <- log(4.40);  label("First-order non-renal clearance CL_NR for a 70 kg subject (L/h)")              # Table 3 column "NONMEM"
    lcl_renal  <- log(5.70);  label("First-order renal clearance CL_R for a 70 kg subject (L/h)")                   # Table 3 column "NONMEM"
    lvmax      <- log(170);   label("Maximum rate of mixed-order renal elimination V_max for a 70 kg subject (mg/h)") # Table 3 column "NONMEM"
    lkm        <- log(49.7);  label("Michaelis-Menten constant K_m for mixed-order renal elimination (mg/L)")       # Table 3 column "NONMEM"
    lvc        <- log(7.00);  label("Volume of central compartment V_1 for a 70 kg subject (L)")                    # Table 3 column "NONMEM"
    lvp        <- log(2.95);  label("Volume of shallow peripheral compartment V_2 for a 70 kg subject (L)")         # Table 3 column "NONMEM"
    lvp2       <- log(2.71);  label("Volume of deep peripheral compartment V_3 for a 70 kg subject (L)")            # Table 3 column "NONMEM"
    lq         <- log(12.7);  label("Inter-compartmental clearance central <-> shallow peripheral CLic_shallow for a 70 kg subject (L/h)") # Table 3 column "NONMEM"
    lq2        <- log(1.28);  label("Inter-compartmental clearance central <-> deep peripheral CLic_deep for a 70 kg subject (L/h)")       # Table 3 column "NONMEM"

    # Allometric exponents are reported as fixed structural assumptions
    # (Table 3 footnotes * and +). Exponent 0.75 on clearance terms and on
    # V_max; exponent 1.0 on volumes.
    e_wt_cl <- fixed(0.75); label("Allometric exponent applied to all clearance-like terms (CL_NR, CL_R, V_max, CL_ic_shallow, CL_ic_deep) (unitless)")  # Table 3 footnote *
    e_wt_vc <- fixed(1.0);  label("Allometric exponent applied to all volume terms (V_1, V_2, V_3) (unitless)")                                          # Table 3 footnote +

    # Inter-individual variability (Table 3, NONMEM column). The paper
    # reports PPV (= BSV + BOV) as a single apparent CV% for each of four
    # quantities: total CL combined for CL_NR + CL_R (Table 3 footnote ||),
    # V_ss, V_max, and K_m. Mechanically the NONMEM model carried a single
    # eta on CL (applied multiplicatively to both CL_NR and CL_R) and a
    # single eta on V_ss (applied multiplicatively to V_1, V_2, V_3), so the
    # shared-eta names below pair each random effect with the
    # parameters it multiplies. Off-diagonal correlations: r(CL_tot, V_ss) =
    # 0.84 (Table 3 footnote ++); r(V_max, K_m) = 0.99 (Table 3 footnote
    # SS). The conversion omega^2 = log(CV^2 + 1) takes the reported CV%
    # to the internal log-scale variance. Inter-compartmental clearances
    # CL_ic_shallow and CL_ic_deep had no PPV reported in the NONMEM
    # column (Table 3) and are typical-value-only.

    # Shared CL eta (multiplies both CL_NR and CL_R) and shared V_ss eta
    # (multiplies V_1, V_2, V_3) with their cross-correlation. Variances:
    # log(0.0962^2 + 1) = 0.009211; log(0.135^2 + 1) = 0.018061; covariance
    # = 0.84 * sqrt(0.009211 * 0.018061) = 0.010835.
    etalcl_renal_nonren + etalvc_vp_vp2 ~ c(0.009211, 0.010835, 0.018061) # Table 3 (9.62 pct PPV for CL combined, 13.5 pct PPV for V_ss; r = 0.84 footnote ++)

    # Correlated V_max / K_m etas. Variances: log(0.504^2 + 1) = 0.226374;
    # log(1.50^2 + 1) = 1.178655; covariance = 0.99 * sqrt(0.226374 *
    # 1.178655) = 0.511422.
    etalvmax + etalkm ~ c(0.226374, 0.511422, 1.178655)  # Table 3 (50.4 pct PPV for V_max, 150 pct PPV for K_m; r = 0.99 footnote SS)

    # Combined additive plus proportional residual error on plasma drug
    # concentrations and on cumulative amounts excreted unchanged in urine
    # (Methods, "Observation model"; Table 3 column "NONMEM").
    propSd        <- 0.125;  label("Proportional residual error on plasma concentrations (fraction)") # Table 3, CV_C
    addSd         <- 0.447;  label("Additive residual error on plasma concentrations (mg/L)")        # Table 3, SD_C
    propSd_Aurine <- 0.246;  label("Proportional residual error on cumulative urine amounts (fraction)") # Table 3, CV_AU
    addSd_Aurine  <- 3.89;   label("Additive residual error on cumulative urine amounts (mg)")        # Table 3, SD_AU
  })

  model({
    # Individual structural parameters. The shared etalcl_renal_nonren is
    # applied to both CL_NR and CL_R per the Bulitta 2010 NONMEM
    # parameterisation (Table 3 footnote ||); the shared etalvc_vp_vp2 is
    # applied to all three volumes per the V_ss IIV structure in the same
    # column. Allometric scaling is fixed at 0.75 on clearance terms and
    # V_max and at 1.0 on volumes (Table 3 footnotes * and +).
    cl_nonren <- exp(lcl_nonren + etalcl_renal_nonren) * (WT / 70)^e_wt_cl
    cl_renal  <- exp(lcl_renal  + etalcl_renal_nonren) * (WT / 70)^e_wt_cl
    vmax      <- exp(lvmax      + etalvmax)            * (WT / 70)^e_wt_cl
    km        <- exp(lkm        + etalkm)
    vc        <- exp(lvc        + etalvc_vp_vp2)       * (WT / 70)^e_wt_vc
    vp        <- exp(lvp        + etalvc_vp_vp2)       * (WT / 70)^e_wt_vc
    vp2       <- exp(lvp2       + etalvc_vp_vp2)       * (WT / 70)^e_wt_vc
    q         <- exp(lq)                                * (WT / 70)^e_wt_cl
    q2        <- exp(lq2)                               * (WT / 70)^e_wt_cl

    # Inter-compartmental rate constants
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Plasma drug concentration in the central compartment (mg/L)
    Cc <- central / vc

    # Elimination fluxes (mg/h). The renal pathway is parallel first-order
    # (cl_renal) plus mixed-order Michaelis-Menten (vmax, km); the non-renal
    # pathway is first-order (cl_nonren).
    flux_renal_lin <- cl_renal  * Cc
    flux_renal_sat <- vmax * Cc / (km + Cc)
    flux_nonrenal  <- cl_nonren * Cc

    # Three-compartment disposition with parallel renal and non-renal
    # elimination from the central compartment. The renal pathways feed the
    # urine compartment so the cumulative amount excreted unchanged can be
    # observed; the non-renal pathway leaves the system.
    d/dt(central)     <- -flux_renal_lin - flux_renal_sat - flux_nonrenal -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(urine)       <-  flux_renal_lin + flux_renal_sat

    # Observable for the cumulative amount excreted in urine (mg).
    Aurine <- urine

    Cc     ~ add(addSd)        + prop(propSd)
    Aurine ~ add(addSd_Aurine) + prop(propSd_Aurine)
  })
}

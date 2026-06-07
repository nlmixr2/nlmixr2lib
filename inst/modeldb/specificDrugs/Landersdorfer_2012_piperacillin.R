Landersdorfer_2012_piperacillin <- function() {
  description <- "Three-compartment population PK model for piperacillin in healthy adult volunteers after intravenous infusion, with parallel first-order plus mixed-order (Michaelis-Menten) renal clearance and first-order non-renal clearance; a urine compartment accumulates the renally excreted amount (Landersdorfer 2012 Model 3, the final model)"
  reference <- paste(
    "Landersdorfer CB, Bulitta JB, Kirkpatrick CMJ, Kinzig M, Holzgrabe U,",
    "Drusano GL, Stephan U, Sorgel F.",
    "Population pharmacokinetics of piperacillin at two dose levels:",
    "influence of nonlinear pharmacokinetics on the pharmacodynamic profile.",
    "Antimicrob Agents Chemother. 2012;56(11):5715-5723.",
    "doi:10.1128/AAC.00937-12"
  )
  vignette <- "Landersdorfer_2012_piperacillin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 10L,
    n_studies       = 1L,
    age_range       = "23 to 31 years; mean 25.7 (SD 3.1) years",
    age_median      = "mean 25.7 years (Results, baseline demographics paragraph)",
    weight_range    = "mean 69.6 (SD 9.7) kg",
    weight_median   = "mean 69.6 kg (Results)",
    height_range    = "mean 177.5 (SD 8.0) cm",
    sex_female_pct  = 50,
    race_ethnicity  = c(White = 100),
    crcl_range      = "76 to 125 mL/min (Cockcroft-Gault); all subjects with normal renal function",
    disease_state   = "Healthy Caucasian adult volunteers",
    dose_range      = "Single intravenous 5-min infusion of 1500 mg or 3000 mg piperacillin; randomized two-way crossover with a washout of at least 4 days",
    regions         = "Single-centre Germany (University Hospital of Essen)",
    notes           = "5 male and 5 female. Piperacillin was administered alone in the study; the paper notes piperacillin PK is unaffected by tazobactam at clinical 4:1 or 8:1 ratios. Baseline demographics from Results section. Sampling schedule and assay details in Methods."
  )

  ini({
    # Structural parameters: Landersdorfer 2012 Table 2 (NONMEM, FOCE+I), Model 3
    # final model. The renal-elimination pathway is parallel first-order
    # (cl_renal) plus mixed-order Michaelis-Menten (vmax / km); the non-renal
    # pathway is first-order (cl_nonren). Concentration in the central
    # compartment is Cc = central / vc with dose in mg and vc in L, giving Cc
    # in mg/L. With vmax in mg/h and km in mg/L, the Michaelis-Menten flux
    # vmax * Cc / (km + Cc) is in mg/h; the linear arms cl_renal * Cc and
    # cl_nonren * Cc are also in mg/h.
    lcl_renal  <- log(4.42);  label("First-order renal clearance CL_R (L/h)")             # Table 2 Model 3
    lvmax      <- log(219);   label("Maximum rate of mixed-order renal elimination V_maxR (mg/h)") # Table 2 Model 3
    lkm        <- log(36.1);  label("Michaelis-Menten constant for mixed-order renal elimination Km_R (mg/L)") # Table 2 Model 3
    lcl_nonren <- log(5.44);  label("First-order non-renal clearance CL_NR (L/h)")        # Table 2 Model 3
    lvc        <- log(6.32);  label("Central volume of distribution V_1 (L)")             # Table 2 Model 3
    lvp        <- log(3.59);  label("Shallow peripheral volume of distribution V_2 (L)")  # Table 2 Model 3
    lvp2       <- log(2.69);  label("Deep peripheral volume of distribution V_3 (L)")     # Table 2 Model 3
    lq         <- log(15.2);  label("Inter-compartmental clearance central <-> shallow peripheral CLic_shallow (L/h)") # Table 2 Model 3
    lq2        <- log(1.65);  label("Inter-compartmental clearance central <-> deep peripheral CLic_deep (L/h)")     # Table 2 Model 3

    # Inter-individual variability. Table 2 reports BSV %CV from a full
    # covariance matrix on all structural parameters except CLic_shallow and
    # CLic_deep (footnote a; NONMEM, FOCE+I, Model 3). The off-diagonal
    # covariances of that full matrix are not tabulated, so only the diagonal
    # CV%s are encoded here; see vignette Assumptions and deviations.
    # omega^2 = log(CV^2 + 1) converts each reported CV% to the internal
    # log-scale variance.
    etalcl_renal  ~ 0.199533  # CL_R   47% CV -- Table 2 Model 3; log(0.47^2 + 1)
    etalvmax      ~ 0.534007  # V_maxR 84% CV -- Table 2 Model 3; log(0.84^2 + 1)
    etalkm        ~ 0.812973  # Km_R  112% CV -- Table 2 Model 3; log(1.12^2 + 1)
    etalcl_nonren ~ 0.031881  # CL_NR  18% CV -- Table 2 Model 3; log(0.18^2 + 1)
    etalvc        ~ 0.031881  # V_1    18% CV -- Table 2 Model 3; log(0.18^2 + 1)
    etalvp        ~ 0.207302  # V_2    48% CV -- Table 2 Model 3; log(0.48^2 + 1)
    etalvp2       ~ 0.022254  # V_3    15% CV -- Table 2 Model 3; log(0.15^2 + 1)
    # No BSV reported (NONMEM, Table 2) on lq (CLic_shallow) and lq2 (CLic_deep).

    # Residual error: combined additive + proportional on plasma drug
    # concentrations and on amounts excreted unchanged in urine (Methods,
    # 'Observation model and computation'). Table 2 Model 3: CV_C = 12.8%,
    # SD_C = 0.26 mg/L; CV_AU = 24.7%, SD_AU = 4.17 mg.
    propSd        <- 0.128;  label("Proportional residual error on plasma concentrations (fraction)") # Table 2 Model 3 (CV_C)
    addSd         <- 0.26;   label("Additive residual error on plasma concentrations (mg/L)")        # Table 2 Model 3 (SD_C)
    propSd_Aurine <- 0.247;  label("Proportional residual error on cumulative urine amounts (fraction)") # Table 2 Model 3 (CV_AU)
    addSd_Aurine  <- 4.17;   label("Additive residual error on cumulative urine amounts (mg)")        # Table 2 Model 3 (SD_AU)
  })

  model({
    # Individual structural parameters
    cl_renal  <- exp(lcl_renal  + etalcl_renal)
    vmax      <- exp(lvmax      + etalvmax)
    km        <- exp(lkm        + etalkm)
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    vc        <- exp(lvc        + etalvc)
    vp        <- exp(lvp        + etalvp)
    vp2       <- exp(lvp2       + etalvp2)
    q         <- exp(lq)
    q2        <- exp(lq2)

    # Inter-compartmental rate constants
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Plasma drug concentration in the central compartment (mg/L)
    Cc <- central / vc

    # Elimination fluxes (mg/h). The renal pathway is parallel first-order
    # plus mixed-order Michaelis-Menten; the non-renal pathway is first-order.
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

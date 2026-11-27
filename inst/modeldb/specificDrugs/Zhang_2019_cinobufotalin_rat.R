Zhang_2019_cinobufotalin_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous-bolus pharmacokinetic model",
    "for cinobufotalin, a bufadienolide constituent of Cinobufacini injection,",
    "in male Sprague-Dawley rats with and without diethylnitrosamine",
    "(DEN)-induced hepatocellular carcinoma (Zhang 2019). Each rat received a",
    "single 2.5 mg/kg jugular-vein bolus and was fitted individually with the",
    "WinNonlin 3.3 two-compartment IV-bolus model; the central volume,",
    "peripheral volume, elimination clearance and intercompartmental clearance",
    "are the per-group means of Table 4, selected by TUMTP_HCC (0 = saline",
    "control rats, 1 = DEN-injured rats) because the two groups were fitted",
    "separately rather than through an estimated covariate function. Parameters",
    "are expressed per kg body weight (L/kg, L/h/kg). No between-subject",
    "variability or residual-error model was reported (two-stage individual",
    "fits), so parameters are fixed at the published means and the residual SD",
    "is fixed at zero."
  )
  reference <- paste(
    "Zhang X, Liu T, Zhang Y, Liu F, Li H, Fang D, Wang C, Sun H, Xie S (2019).",
    "Elucidation of the Differences in Cinobufotalin's Pharmacokinetics Between",
    "Normal and Diethylnitrosamine-Injured Rats: The Role of P-Glycoprotein.",
    "Front Pharmacol 10:521. doi:10.3389/fphar.2019.00521"
  )
  vignette <- "Zhang_2019_cinobufotalin"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list(
    TUMTP_HCC = list(
      description = "Diethylnitrosamine-induced hepatocellular carcinoma indicator; selects the group-specific two-compartment disposition parameter set",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = saline-treated normal rats",
      notes = paste(
        "Time-fixed per rat. Zhang 2019 fitted the normal and DEN-injured groups",
        "as two independent two-compartment models (Table 4), so this column is",
        "a group selector rather than a continuous covariate with an estimated",
        "effect. 1 = DEN-injured rat (50 mg/kg intraperitoneal DEN twice weekly",
        "for 14 weeks, an HCC model); 0 = normal saline control."
      ),
      source_name = "group (normal vs DEN-injured)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "cinobufotalin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cinobufotalin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 12,
    n_studies = 1,
    sex = "male",
    weight_range = "200-220 g",
    disease_state = "diethylnitrosamine-induced hepatocellular carcinoma (DEN-injured) vs normal control",
    dose_range = "2.5 mg/kg single intravenous bolus (jugular vein)",
    regions = "China",
    notes = paste(
      "Zhang 2019 Pharmacokinetic Studies: male Sprague-Dawley rats, 200-220 g.",
      "The DEN-injured HCC model was induced with 50 mg/kg intraperitoneal",
      "diethylnitrosamine twice weekly for 14 weeks; controls received saline on",
      "the same schedule. Cinobufotalin was given as a single 2.5 mg/kg bolus",
      "through a jugular-vein cannula and plasma sampled at 5, 10, 15, 30,",
      "45 min and 1, 1.5, 2, 4, 6, 8, 12, 24 h. Table 4 PK parameters are",
      "reported as mean +/- SD with n = 6 per group (the two groups contribute",
      "the two parameter sets of this model)."
    )
  )

  ini({
    # Two-compartment IV-bolus disposition, fitted independently per group
    # (Zhang 2019 Table 4). The macro set V (central volume), V2 (peripheral
    # volume), CL (elimination clearance) and CL2 (intercompartmental
    # clearance) reproduces the printed K10/K12/K21 micro-constants and both
    # half-lives, and dose/CL reproduces the printed AUC, so it is the set
    # used here. Values are per kg body weight.
    # --- Normal control rats (TUMTP_HCC = 0) ---
    lvc_normal <- fixed(log(2.1)); label("Central volume V, normal rats (L/kg)")  # Table 4, V, normal = 2.1 +/- 0.6
    lvp_normal <- fixed(log(13.1)); label("Peripheral volume V2, normal rats (L/kg)")  # Table 4, V2, normal = 13.1 +/- 1.3
    lcl_normal <- fixed(log(1.0)); label("Elimination clearance CL, normal rats (L/h/kg)")  # Table 4, CL, normal = 1.0 +/- 0.2
    lq_normal <- fixed(log(3.3)); label("Intercompartmental clearance CL2, normal rats (L/h/kg)")  # Table 4, CL2, normal = 3.3 +/- 0.4
    # --- DEN-injured HCC rats (TUMTP_HCC = 1) ---
    lvc_den <- fixed(log(2.4)); label("Central volume V, DEN-injured rats (L/kg)")  # Table 4, V, DEN-injured = 2.4 +/- 0.6
    lvp_den <- fixed(log(31.5)); label("Peripheral volume V2, DEN-injured rats (L/kg)")  # Table 4, V2, DEN-injured = 31.5 +/- 10.4
    lcl_den <- fixed(log(2.4)); label("Elimination clearance CL, DEN-injured rats (L/h/kg)")  # Table 4, CL, DEN-injured = 2.4 +/- 0.6
    lq_den <- fixed(log(8.1)); label("Intercompartmental clearance CL2, DEN-injured rats (L/h/kg)")  # Table 4, CL2, DEN-injured = 8.1 +/- 1.7

    propSd <- fixed(0); label("Proportional residual SD on cinobufotalin plasma concentration (not reported)")  # Zhang 2019 reports no residual-error model (two-stage individual fits)
  })

  model({
    # Group selector: DEN-injured vs normal control.
    isDen <- (TUMTP_HCC > 0.5)

    vc <- isDen * exp(lvc_den) + (1 - isDen) * exp(lvc_normal)
    vp <- isDen * exp(lvp_den) + (1 - isDen) * exp(lvp_normal)
    cl <- isDen * exp(lcl_den) + (1 - isDen) * exp(lcl_normal)
    q <- isDen * exp(lq_den) + (1 - isDen) * exp(lq_normal)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}

Huang_2019_tilmicosin_invitro_pk <- function() {
  description <- "In vitro (dialysis-type dynamic PK/PD apparatus simulating chicken lung tissue). One-compartment pharmacokinetic model with first-order absorption for tilmicosin in the in vitro dynamic model of Huang 2019, an apparatus of an absorption chamber (61.54 mL, the administration site) draining into a central chamber (500 mL external compartment holding a 10 mL dialysis-tube internal compartment that carries the Mycoplasma gallisepticum culture), with fresh medium pumped in and spent medium pumped out at 0.13 mL/min so that the medium concentration mimics the lung-tissue time course of tilmicosin in chickens (absorption half-life 5.37 h, elimination half-life 42.83 h). Huang 2019 fitted the one-compartment open model with first-order absorption separately to each of seven dose groups in WinNonlin (Table 2, R^2 0.993-0.999); the packaged typical values are the across-group means printed in Table 2 (absorption half-life 5.25 h, elimination half-life 41.51 h, central volume 530.19 mL). Doses are the in vitro amounts injected into the absorption chamber (Table 2 row D: 0.089 to 1.786 mg, i.e. 0.0893 x the nominal 1 to 20 mg clinical dose group); there is no bioavailability term, since the whole injected amount reaches the central chamber. No between-replicate variability or residual error magnitude was reported, so no eta parameters are present and propSd is FIXED at 0. Sibling models from the same paper: Huang_2019_tilmicosin_killrate (static time-kill sigmoid Emax kill-rate model), Huang_2019_tilmicosin_aucmic and Huang_2019_tilmicosin_tmic (inhibitory sigmoid Emax PK/PD-index models fitted to the dynamic-model kill data)."
  reference <- paste(
    "Huang Z, Wu Y, Zhou Z, Xia X, Gu X, Cai Q, Shen X, Yang H, Ding H. (2019).",
    "Pharmacokinetic and Pharmacodynamic Integration and Resistance Analysis",
    "of Tilmicosin Against Mycoplasma gallisepticum in an In Vitro Dynamic Model.",
    "Frontiers in Pharmacology 10:670. doi:10.3389/fphar.2019.00670.",
    "PMCID: PMC6598723.",
    "Apparatus: Materials and Methods, 'Description of the Model' and Figure 1.",
    "Target chicken-lung PK and dose groups: 'In Vitro Pharmacokinetic and",
    "Pharmacodynamic Model and Dosing Regimens'. Fitted PK parameters: Table 2.",
    "Concentration-time curves: Figure 5.",
    sep = " "
  )
  vignette <- "Huang_2019_tilmicosin"
  units <- list(
    time = "h",
    dosing = "mg (tilmicosin injected into the absorption chamber of the in vitro apparatus)",
    concentration = "ug/mL (= mg/L, tilmicosin in the central-chamber medium)"
  )

  # In vitro apparatus: the central chamber holds culture medium, not a
  # biological matrix.
  compartmentData <- list(
    depot = list(
      analyte = "tilmicosin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tilmicosin",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (dynamic PK/PD apparatus; Mycoplasma gallisepticum strain S6 culture)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Three-chamber in vitro dynamic model driven by a multi-channel peristaltic pump at 0.13 mL/min: an absorption chamber (61.54 mL drug-containing medium, the administration site), a central chamber (500 mL sterile medium external compartment containing a 10 mL dialysis-tube internal compartment) and a reserve chamber of fresh medium, with a waste bottle; contents kept at 37 C and mixed by magnetic stirrers",
    target_pk = "Tilmicosin lung-tissue PK in chickens: absorption half-life 5.37 h, elimination half-life 42.83 h, apparent volume of distribution in lung 5.6 L/kg (laboratory data, unpublished)",
    dose_range = "Seven nominal clinical dose groups of 1, 2.5, 5, 7.5, 10, 15 and 20 mg, realised in vitro as 0.089, 0.223, 0.446, 0.670, 0.893, 1.339 and 1.786 mg injected once into the absorption chamber at time 0",
    sampling = "1.5 mL samples from the external compartment at 1, 3, 6, 9, 12, 24, 36, 48, 72, 96, 120, 144 and 168 h; tilmicosin by HPLC-MS/MS",
    disease_state = "not applicable (in vitro)",
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = "Huang 2019 Table 2 reports per-dose-group WinNonlin estimates of the absorption and elimination half-lives, Tmax, Cmax and central volume; the per-group values are tabulated in the vignette and reproduce the printed Tmax to within 0.01 h when substituted in the closed-form solution. The packaged model uses the across-group means printed in the Table 2 mean column."
  )

  ini({
    # Huang 2019 Table 2, mean column (X-bar) across the seven dose groups.
    # The paper reports half-lives; the rate constants are log(2)/t1/2.
    lka <- log(log(2) / 5.25)
    label("Log first-order absorption rate constant from the absorption chamber (1/h)") # Huang 2019 Table 2, t1/2Ka mean = 5.25 h
    lkel <- log(log(2) / 41.51)
    label("Log first-order elimination rate constant from the central chamber (1/h)") # Huang 2019 Table 2, t1/2Kel mean = 41.51 h
    lvc <- log(530.19 / 1000)
    label("Log central-chamber volume (L)") # Huang 2019 Table 2, VCentral mean = 530.19 mL

    # Huang 2019 reports only R^2 per dose group (0.993-0.999, Table 2),
    # no residual standard deviation.
    propSd <- fixed(0)
    label("Proportional residual SD on tilmicosin concentration (not reported in Huang 2019)") # Huang 2019 Table 2 reports R^2 only
  })

  model({
    ka <- exp(lka)
    kel <- exp(lkel)
    vc <- exp(lvc)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

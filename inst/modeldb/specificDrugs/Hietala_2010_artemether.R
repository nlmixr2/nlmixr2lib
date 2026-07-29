# Joint parent-metabolite population PK model of oral artemether (ARM) and its
# active metabolite dihydroartemisinin (DHA) in Tanzanian children with
# uncomplicated falciparum malaria treated with the artemether-lumefantrine
# combination Coartem (Hietala 2010, Antimicrob Agents Chemother 54:4780-4788;
# doi:10.1128/AAC.00252-10).

Hietala_2010_artemether <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral artemether (ARM)",
    "and its active metabolite dihydroartemisinin (DHA) in 50 Tanzanian",
    "children (ages 1-10 years, weights 8-30 kg) with uncomplicated",
    "Plasmodium falciparum malaria treated with the standard six-dose",
    "weight-based Coartem (artemether 20 mg + lumefantrine 120 mg per",
    "tablet) regimen at 0, 8, 24, 36, 48, and 60 hours (Hietala 2010).",
    "Absorption is first-order with ka fixed at 1/h. Disposition is",
    "two-compartment for ARM with complete in-vivo conversion to DHA",
    "(bioavailability of DHA fixed at 1 to render the metabolite model",
    "identifiable); DHA disposition is one-compartment. The apparent",
    "oral clearance of ARM is time-dependent with a linear increase per",
    "dose number occasion (CL/F_ARM = theta1 * (1 + theta2 * (OCC - 1)),",
    "OCC = 1..6), reproducing a ~3.4-fold rise over the six-dose",
    "regimen attributed to enzyme induction. All PK parameters are",
    "reported per kilogram body weight (linear weight normalisation",
    "applied inside model()).",
    sep = " "
  )
  reference <- paste(
    "Friberg Hietala S, Martensson A, Ngasala B, Dahlstrom S, Lindegardh N,",
    "Annerberg A, Premji Z, Farnert A, Gil P, Bjorkman A, Ashton M (2010).",
    "Population pharmacokinetics and pharmacodynamics of artemether and",
    "lumefantrine during combination treatment in children with uncomplicated",
    "falciparum malaria in Tanzania. Antimicrob Agents Chemother",
    "54(11):4780-4788. doi:10.1128/AAC.00252-10.",
    sep = " "
  )
  vignette <- "Hietala_2010_artemether_lumefantrine_malaria"
  units <- list(time = "hour", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (per-kg) weight normalisation of CL / Q / V parameters.",
        "Hietala 2010 reports all population PK parameters on a per-kg",
        "basis (Table 1: 'CL/F_ARM (liters/h/kg)', 'V_c/F_ARM (liters/kg)',",
        "etc.), so individual CL = lcl * WT, Vc = lvc * WT, etc. inside",
        "model(). Population mean WT 14 kg (range 8-30 kg, Results)."
      ),
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer-valued dose-occasion (dose number) indicator",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Values 1..6 identify the dose-occasion of the six-dose Coartem",
        "regimen (doses at 0, 8, 24, 36, 48, and 60 hours). Used to",
        "encode the linear time-dependent increase in apparent oral",
        "artemether clearance described by the Methods / Results:",
        "CL/F_ARM = theta1 * (1 + theta2 * (OCC - 1)), with",
        "theta1 = 2.6 L/h/kg (first dose) and theta2 = 0.57 (fractional",
        "increase per occasion); reaches ~10 L/h/kg at OCC = 6, a",
        "~3.4-fold rise attributed to enzyme induction (Discussion).",
        "Time-varying within subject: the per-dose-record OCC propagates",
        "through to the immediately-following sampling occasion until",
        "the next dose increments it."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 50L,
    n_studies       = 1L,
    n_observations  = 797L,
    age_range       = "1-10 years (mean 4; Results)",
    weight_range    = "8-30 kg (mean 14; Results)",
    sex_female_pct  = 62,
    disease_state   = paste(
      "Acute uncomplicated Plasmodium falciparum malaria with asexual",
      "parasite density 2,000-200,000 / microL at admission and either",
      "axillary temperature >= 37.5 degC or history of fever within 24 h",
      "(Methods)."
    ),
    dose_range      = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per",
      "tablet. Weight-based dosing: 5-14 kg -> 1 tablet/dose,",
      "15-24 kg -> 2 tablets/dose, 25-34 kg -> 3 tablets/dose. Six",
      "doses (oral) at 0, 8, 24, 36, 48, and 60 hours. Half of patients",
      "received each dose with 200 mL full-fat (3.4%) cow's milk; the",
      "other half with water (no milk effect on PK was retained for",
      "lumefantrine; ARM PK was not covariate-tested for milk)."
    ),
    regions            = "Tanzania (Fukayosi Primary Health Care Centre, Bagamoyo District)",
    trial_registration = "ClinicalTrials.gov NCT00336375",
    notes              = paste(
      "Demographics from Hietala 2010 Results. 397 ARM and 400 DHA plasma",
      "concentrations from 50 patients were included in the PK model",
      "(plus 9 ARM, 21 DHA below limit-of-detection samples excluded; 12",
      "ARM and 47 DHA below-LLOQ samples set to LLOQ/2 = 2 and 3 nM",
      "respectively). Companion lumefantrine PK model from the same",
      "cohort: modellib('Hietala_2010_lumefantrine'). Companion",
      "parasitemia PD model driven by the same ARM / DHA PK:",
      "modellib('Hietala_2010_artemether_parasitemia')."
    )
  )

  ini({
    # Structural PK parameters from Hietala 2010 Table 1 ("Estimate (95% CI)" column).
    # All clearances and volumes are reported per kilogram body weight; linear
    # (per-kg) weight normalisation is applied inside model() (cl <- lcl * WT, etc.).
    # The Table 1 caption for "Q ARM (liters/kg)" omits "/h" -- it is treated here as
    # L/h/kg consistent with its description "Intercompartment clearance" (a flux)
    # and consistent with the other CL / Q clearance entries in the same table.
    lka      <- fixed(log(1))    ; label("Absorption rate constant ka (1/h); fixed")                  # Hietala 2010 Table 1: ka = 1 (fixed); Results "the absorption rate constant was fixed to 1/h"
    lcl      <- log(2.6)         ; label("Apparent oral artemether clearance at OCC = 1, theta1 (L/h/kg)") # Hietala 2010 Table 1: theta1 = 2.6 (95% CI 1.5-2.6)
    e_occ_cl <- 0.57             ; label("Fractional increase in CL/F_ARM per dose occasion, theta2 (per OCC step)") # Hietala 2010 Table 1: theta2 = 0.57 (95% CI 0.39-0.75)
    lvc      <- log(5.2)         ; label("Apparent artemether central volume Vc/F_ARM (L/kg)")           # Hietala 2010 Table 1: Vc/F_ARM = 5.2 (95% CI 3.5-7.1)
    lq       <- log(1.4)         ; label("Apparent artemether intercompartmental clearance Q/F_ARM (L/h/kg)") # Hietala 2010 Table 1: Q_ARM = 1.4 (95% CI 1.1-1.8); units L/h/kg per Discussion / Methods (Table caption omits '/h')
    lvp      <- log(41.4)        ; label("Apparent artemether peripheral volume Vp/F_ARM (L/kg)")        # Hietala 2010 Table 1: Vp/F_ARM = 41.4 (95% CI 29.0-58.1)
    lcl_dihydroart  <- log(6.8)         ; label("Apparent DHA metabolite clearance CL/F_DHA (L/h/kg); F_DHA fixed to 1") # Hietala 2010 Table 1: CL/F_DHA = 6.8 (95% CI 5.8-8.0)
    lvc_dihydroart  <- log(3.7)         ; label("Apparent DHA metabolite volume V/F_DHA (L/kg); F_DHA fixed to 1")      # Hietala 2010 Table 1: V/F_DHA = 3.7 (95% CI 2.3-8.7)

    # Inter-individual variability (IIV). Table 1 reports IIV as %CV per the
    # NONMEM log-normal convention; the corresponding internal variance is
    #   omega^2 = log((CV/100)^2 + 1)
    #     CL/F_ARM CV 41% -> log(0.41^2 + 1) = 0.155649
    #     CL/F_DHA CV 47% -> log(0.47^2 + 1) = 0.198101
    # IIV on the other parameters (ka, Vc, Q, Vp, V_DHA) was not estimated
    # (Table 1 "NE" entries) and is therefore omitted here.
    etalcl     ~ 0.155649    # Hietala 2010 Table 1: IIV on CL/F_ARM = 41% CV (95% CI 37-50)
    etalcl_dihydroart ~ 0.198101    # Hietala 2010 Table 1: IIV on CL/F_DHA = 47% CV (95% CI 35-57)

    # Residual error. The paper reports separate proportional and additive
    # residual errors for ARM and DHA (Table 1). nlmixr2's combined error
    # form 'Cc ~ add(addSd) + prop(propSd)' takes the additive SD on the
    # concentration scale (nM here, matching the LC-MS/MS quantification
    # used) and the proportional SD as a fraction. The additive components
    # were fixed by the paper (Table 1 caption marks them "fixed").
    propSd     <- 0.61          ; label("Proportional residual SD for artemether plasma concentration (fraction)") # Hietala 2010 Table 1: sigma_prop_ARM = 61% (95% CI 54-67)
    addSd      <- fixed(2)      ; label("Additive residual SD for artemether plasma concentration (nM); fixed")    # Hietala 2010 Table 1: sigma_add_ARM = 2 nM (fixed)
    propSd_dihydroart <- 0.82          ; label("Proportional residual SD for DHA plasma concentration (fraction)")        # Hietala 2010 Table 1: sigma_prop_DHA = 82% (95% CI 73-90)
    addSd_dihydroart  <- fixed(3)      ; label("Additive residual SD for DHA plasma concentration (nM); fixed")           # Hietala 2010 Table 1: sigma_add_DHA = 3 nM (fixed)
  })

  model({
    # Molecular weights (g/mol). Artemether 298.4, dihydroartemisinin 284.3.
    # Used at the ARM -> DHA formation step to convert mass-rate of artemether
    # elimination into mass-rate of DHA formation under the source paper's
    # assumption of complete in-vivo 1:1 molar conversion of ARM to DHA
    # (Results: "it was assumed that all ARM is eliminated via this pathway").
    # With doses in mg, central and central_dihydroart are mass amounts (mg).
    mw_arm <- 298.4
    mw_dihydroart <- 284.3

    # Individual PK parameters with linear (per-kg) weight normalisation.
    # The time-dependent CL/F_ARM is modelled as
    #   CL/F_ARM(OCC, individual) = theta1 * (1 + theta2 * (OCC - 1)) * exp(eta)
    # (Hietala 2010 Results / Table 1). OCC = 1 at the first dose-occasion
    # and steps to 6 at the sixth dose; in the absence of any dose record
    # (OCC = 0 in the dataset) the (OCC - 1) factor produces -1, so use of
    # the model requires the dataset to set OCC >= 1 on every record.
    ka     <-  exp(lka)
    cl     <- (exp(lcl     + etalcl)     * (1 + e_occ_cl * (OCC - 1))) * WT
    vc     <-  exp(lvc)                                                * WT
    q      <-  exp(lq)                                                 * WT
    vp     <-  exp(lvp)                                                * WT
    cl_dihydroart <-  exp(lcl_dihydroart + etalcl_dihydroart)                               * WT
    vc_dihydroart <-  exp(lvc_dihydroart)                                            * WT

    # Two-compartment ARM disposition with first-order absorption from depot,
    # complete conversion to DHA at rate cl/vc * central (with stoichiometric
    # MW factor), and one-compartment DHA disposition.
    kel     <- cl     / vc
    k12     <- q      / vc
    k21     <- q      / vp
    kel_dihydroart <- cl_dihydroart / vc_dihydroart

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                              k12 * central - k21 * peripheral1
    d/dt(central_dihydroart) <-  kel * central * (mw_dihydroart / mw_arm) - kel_dihydroart * central_dihydroart

    # Plasma concentrations in nM. Dose units mg; vc units L; central / vc has
    # units mg/L. Multiplication by 1e6 / MW converts mg/L to nmol/L = nM
    # (1 mg/L / (g/mol) = mmol/L; * 1e6 -> nmol/L).
    Cc     <- 1e6 * central     / vc     / mw_arm
    Cc_dihydroart <- 1e6 * central_dihydroart / vc_dihydroart / mw_dihydroart

    # Combined proportional + additive residual error on the linear-nM scale.
    Cc     ~ add(addSd)     + prop(propSd)
    Cc_dihydroart ~ add(addSd_dihydroart) + prop(propSd_dihydroart)
  })
}

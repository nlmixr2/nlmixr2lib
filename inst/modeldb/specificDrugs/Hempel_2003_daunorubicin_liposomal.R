Hempel_2003_daunorubicin_liposomal <- function() {
  description <- paste(
    "One-compartment IV-infusion population PK model for total daunorubicin",
    "(free plus liposome-encapsulated) following liposomal daunorubicin",
    "(Daunoxome) in paediatric and adolescent oncology patients",
    "(Hempel 2003). Clearance and volume of distribution scale linearly",
    "with total body weight (CL = theta_CL * WT; V = theta_V * WT, i.e.",
    "the source paper's per-kg parameterisation with allometric exponent",
    "fixed to 1 and no reference-weight normalisation). The final model",
    "(Table 2 model 15) retains inter-individual variability on CL",
    "(51% CV) and V (27% CV), inter-occasion variability on CL (16.7%",
    "CV) -- documented but NOT encoded structurally here, per the",
    "Andrews 2017 / Brooks 2021 nlmixr2lib precedent for IOV without an",
    "operational occasion column -- and a proportional residual error",
    "(22%). Distinct from Varatharajan_2016_daunorubicin (free",
    "daunorubicin + daunorubicinol metabolite in adult AML)."
  )
  reference <- "Hempel G, Reinhardt D, Creutzig U, Boos J. Population pharmacokinetics of liposomal daunorubicin in children. Br J Clin Pharmacol. 2003;56(4):370-377. doi:10.1046/j.1365-2125.2003.01886.x"
  vignette <- "Hempel_2003_daunorubicin_liposomal"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at the time of the dose. Enters as a linear (proportional) scalar on both CL and V: CL_ind = theta_CL * WT and V_ind = theta_V * WT (i.e., allometric exponent fixed to 1, no reference-weight normalisation -- the published typical-value parameters are reported per kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = "n/a -- linear scaling with no reference weight; the paper's typical-value CL = 6.41 mL/h/kg and V = 65.4 mL/kg multiply by the individual's weight directly.",
      notes              = paste(
        "Hempel 2003 Methods: 'total body weight was used throughout';",
        "cohort range 14-76.5 kg, median 48.8 kg (Table 1). The source",
        "paper does not state whether weight was time-varying across",
        "cycles; the published parameters do not depend on the choice.",
        "Inter-occasion variability (IOV) of 16.7% CV on CL was retained",
        "in the final model (Table 2 model 15) but is NOT encoded",
        "structurally here -- nlmixr2lib has no canonical occasion-column",
        "convention and the Andrews 2017 / Brooks 2021 precedent is to",
        "omit IOV when no operational OCC mapping exists. Downstream",
        "users who want IOV can add an OCC indicator and a per-occasion",
        "eta in rxode2; see the validation vignette Assumptions and",
        "deviations section."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species          = "human (paediatric and adolescent)",
    n_subjects       = 24L,
    n_studies        = 1L,
    age_range        = "2.84-23.2 years",
    age_median       = "15.4 years",
    weight_range     = "14-76.5 kg",
    weight_median    = "48.8 kg",
    height_range     = "0.89-1.98 m",
    height_median    = "1.61 m",
    bsa_range        = "0.58-1.98 m^2",
    bsa_median       = "1.48 m^2",
    bmi_range        = "12.6-23.5 kg/m^2",
    bmi_median       = "17.6 kg/m^2",
    sex_female_pct   = NA_real_,
    disease_state    = "Predominantly relapsed acute myeloic leukaemia (19/24); five other malignancies (two relapsed acute lymphoblastic leukaemia, one relapsed osteosarcoma, one relapsed Ewing sarcoma, one multiple endocrine neoplasia). All patients had previously received conventional daunorubicin and/or doxorubicin.",
    dose_range       = "Liposomal daunorubicin (Daunoxome) 30-60 mg/m^2 as a 1- to 2.5-hour IV infusion on days 1 and 5 (induction); four patients in the 60 mg/m^2 group intensified to days 1, 3, and 5 in subsequent cycles.",
    regions          = "Germany (16 paediatric haematology/oncology centres enrolling in the AML-REZ 97 protocol of the German Society for Paediatric Oncology and Haematology, GPOH).",
    sampling_window  = "Plasma sampling suggested at end of infusion, 2-4 h, 4-8 h, 12-17 h, 20-28 h after first dose, predose and ~24/48 h after second dose. Overall 214 plasma concentrations from 72 cycles (mean 3 per cycle, 9 per patient).",
    assay            = "Total daunorubicin (free plus liposome-encapsulated; organic-solvent sample preparation destroys liposomes prior to analysis) quantified by capillary electrophoresis with laser-induced fluorescence detection. Within-day accuracy / precision 3.9-12.9% and 3.8-10.7% (n = 8); LOQ 2 ug/L (= 2 ng/mL).",
    notes            = paste(
      "Patient demographics from Table 1. Sex split is not tabulated as a",
      "count; the paper notes 'boys are on average heavier than girls",
      "(median 58.0 vs 37.5 kg)' but no significant sex difference in PK",
      "(differences explained by weight). Cohort skews adolescent (median",
      "15.4 years) and includes one young-adult patient (23.2 years);",
      "the analysis pools paediatric and adolescent patients into a",
      "single typical-value model. NONMEM Version V, FOCE without",
      "interaction (ADVAN1 TRAN2 for the final one-compartment model)."
    )
  )

  ini({
    # Final model (Table 2 model 15): one-compartment IV with linear
    # weight scaling on CL and V, inter-occasion variability on CL only,
    # proportional residual error. The published typical-value
    # parameters are reported per kg body weight so the model below
    # multiplies the typical-value parameter by the individual's WT
    # directly (no reference-weight normalisation).
    lcl <- log(0.00641)  ; label("Clearance per kg body weight (L/h/kg)")               # Table 2 model 15: CL = 6.41 mL/h/kg = 0.00641 L/h/kg; abstract value identical
    lvc <- log(0.0654)   ; label("Volume of distribution per kg body weight (L/kg)")    # Table 2 model 15: V  = 65.4 mL/kg  = 0.0654  L/kg; abstract value identical

    # Inter-individual variability (lognormal; omega^2 = log(1 + CV^2)).
    # IIV on CL = 51% CV -> log(1 + 0.51^2) = 0.23133
    # IIV on V  = 27% CV -> log(1 + 0.27^2) = 0.07033
    etalcl ~ 0.23133   # Table 2 model 15: omega_CL 51 percent CV; abstract +/- 0.5 51 percent
    etalvc ~ 0.07033   # Table 2 model 15: omega_V  27 percent CV; abstract +/- 0.5 27 percent

    # Inter-occasion variability (IOV) on CL was estimated at 16.7% CV in
    # the final model (Table 2 model 15; Results paragraph "The
    # combination of the model with weight as a covariate and IOV on CL
    # (16.7%) resulted in the best model"). NOT encoded structurally
    # here -- see covariateData[["WT"]]$notes and vignette.

    # Residual error: proportional only. A combined proportional plus
    # additive error model gave no improvement (additive term close to
    # zero per Results paragraph "For the residual error, a proportional
    # error model was chosen").
    propSd <- 0.22  ; label("Proportional residual error (fraction)")                   # Table 2 model 15: residual error 22%
  })

  model({
    # Individual PK parameters: linear weight scaling per Hempel 2003
    # Table 2 footnote ("CL = q1 x Weight and V = q2 x Weight").
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT

    # First-order elimination rate constant.
    kel <- cl / vc

    # One-compartment IV. Dose record carries the infusion amount (mg)
    # and rate (mg/h) for the 1- to 2.5-hour Daunoxome infusion; no
    # depot compartment is required.
    d/dt(central) <- -kel * central

    # Plasma total daunorubicin concentration. central is in mg and vc
    # in L, so central / vc has units mg/L. Convert mg/L to ng/mL by
    # multiplying by 1000 (1 mg/L = 1000 ng/mL = 1000 ug/L); the paper
    # reports concentrations in ug/L (numerically identical to ng/mL).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}

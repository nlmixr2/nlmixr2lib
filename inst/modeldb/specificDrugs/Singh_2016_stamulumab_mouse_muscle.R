Singh_2016_stamulumab_mouse_muscle <- function() {
  description <- "Preclinical (mouse, SCID). Landmark direct-effect Emax exposure-response model for the increase in skeletal muscle weight after 12 weeks of weekly intraperitoneal stamulumab (MYO-029), a human IgG1 anti-myostatin (GDF-8) monoclonal antibody, in female SCID mice (Singh 2016). Three independent fits, one per muscle (extensor digitorum longus, gastrocnemius, quadriceps): muscle weight = E0 * (1 + Emax/100 * AUC / (AUC50 + AUC)), where AUC is the steady-state weekly AUC(0-168 h) supplied as the AUC_STAMULUMAB covariate. Typical-value model: no IIV or residual error is reported."
  reference <- "Singh P, Rong H, Gordi T, Bosley J, Bhattacharya I. Translational Pharmacokinetic/Pharmacodynamic Analysis of MYO-029 Antibody for Muscular Dystrophy. Clin Transl Sci. 2016;9(6):302-310. doi:10.1111/cts.12420. PMID 27700008. PMCID PMC5351001. MYO-029 is the development code of the antibody assigned the INN stamulumab."
  vignette <- "Singh_2016_stamulumab"
  units <- list(
    time = "week",
    dosing = "mg",
    concentration = "muscle weight (muscle_edl in mg; muscle_gastro and muscle_quad in g)"
  )

  covariateData <- list(
    AUC_STAMULUMAB = list(
      description = "Steady-state serum AUC of stamulumab over the weekly (0-168 h) dosing interval",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "0 for vehicle-treated (PBS) mice, which gives muscle weight = E0.",
        "Singh 2016 does not say how the mouse AUC was derived for the",
        "intraperitoneal efficacy doses; Figure 3a-c place the 1, 2.5, 5 and",
        "10 mg/kg/week groups at about 2200, 6500, 13000 and 26500 ug*h/mL.",
        "The IV mouse PK model Singh_2016_stamulumab_mouse gives about 3200",
        "ug*h/mL per mg/kg/week at steady state."
      ),
      source_name = "AUC0-168hr"
    )
  )

  population <- list(
    species = "mouse (SCID)",
    n_subjects = 40,
    n_studies = 1,
    age_range = NA,
    weight_range = NA,
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state = "Healthy female SCID mice (muscle-growth efficacy study)",
    dose_range = "PBS vehicle or 1, 2.5, 5 and 10 mg/kg/week IP for 12 weeks, with a double first dose (n = 8 per group)",
    regions = NA,
    notes = paste(
      "Singh 2016 Methods 'Mouse efficacy studies' and Results 'Muscle growth",
      "in mice': muscles collected and weighed at the end of 12 weeks of",
      "treatment. Five groups of eight mice (vehicle plus four doses) are",
      "implied by the design; the paper does not print the total. Table 2 'CV%'",
      "is the estimation precision of each parameter, not between-animal",
      "variability."
    )
  )

  ini({
    # Singh 2016 Table 2, 'Mouse' rows. Table 2 labels the E0 unit 'Gram' for
    # all three muscles; the EDL value 16.7 can only be milligrams (a mouse EDL
    # weighs about 10 mg) and the Figure 3a axis spans 16-22, so it is taken
    # here as mg. Gastrocnemius and quadriceps are in g.
    e0_edl <- 16.7; label("Vehicle EDL muscle weight (mg)") # Table 2, Mouse EDL E0 = 16.7 'Gram' (CV 1%), read as mg
    emax_edl <- 29.9; label("Maximum EDL weight increase over vehicle (%)") # Table 2, Mouse EDL Emax = 29.9% (CV 5%)
    lauc50_edl <- log(1152.6); label("Steady-state weekly AUC giving half the maximum EDL effect (ug*h/mL)") # Table 2, Mouse EDL AUC50 = 1,152.6 ug*h/mL (CV 21%)

    e0_gastro <- 0.24; label("Vehicle gastrocnemius muscle weight (g)") # Table 2, Mouse Gastro E0 = 0.24 Gram (CV 2%)
    emax_gastro <- 23.2; label("Maximum gastrocnemius weight increase over vehicle (%)") # Table 2, Mouse Gastro Emax = 23.2% (CV 13%)
    lauc50_gastro <- log(1560.5); label("Steady-state weekly AUC giving half the maximum gastrocnemius effect (ug*h/mL)") # Table 2, Mouse Gastro AUC50 = 1,560.5 ug*h/mL (CV 50%)

    e0_quad <- 0.34; label("Vehicle quadriceps muscle weight (g)") # Table 2, Mouse Quad E0 = 0.34 Gram (CV 3%)
    emax_quad <- 27.2; label("Maximum quadriceps weight increase over vehicle (%)") # Table 2, Mouse Quad Emax = 27.2% (CV 15%)
    lauc50_quad <- log(1762.6); label("Steady-state weekly AUC giving half the maximum quadriceps effect (ug*h/mL)") # Table 2, Mouse Quad AUC50 = 1,762.6 ug*h/mL (CV 54%)
  })

  model({
    # Table 2 footnote: 'Direct effect Emax model ... effect of steady-state
    # exposures (AUC0-168h) on % increase in muscle mass'; Emax is the
    # 'maximum % change over vehicle', so it multiplies the vehicle weight E0.
    muscle_edl <- e0_edl * (1 + emax_edl / 100 * AUC_STAMULUMAB / (exp(lauc50_edl) + AUC_STAMULUMAB))
    muscle_gastro <- e0_gastro * (1 + emax_gastro / 100 * AUC_STAMULUMAB / (exp(lauc50_gastro) + AUC_STAMULUMAB))
    muscle_quad <- e0_quad * (1 + emax_quad / 100 * AUC_STAMULUMAB / (exp(lauc50_quad) + AUC_STAMULUMAB))
  })
}

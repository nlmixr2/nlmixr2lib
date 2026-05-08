Archary_2018_lopinavir <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral lopinavir/ritonavir in severely malnourished HIV-infected children, with FFM allometric scaling and a linear total-cholesterol effect on apparent clearance (Archary 2018)."
  reference <- "Archary M, McIlleron H, Bobat R, La Russa P, Sibaya T, Wiesner L, Hennig S. Population Pharmacokinetics of Lopinavir in Severely Malnourished HIV Infected Children and the Effect on Treatment Outcomes. Pediatr Infect Dis J. 2018;37(4):349-355. doi:10.1097/INF.0000000000001867"
  vignette <- "Archary_2018_lopinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling reference 5.6 kg (Archary 2018 study-cohort mean baseline FFM 5.1-5.5 kg per Table 1). FFM is computed from total body weight, height, and sex per Al-Sallami et al. Clin Pharmacokinet 2015;54(11):1169-1178 (paper reference 20). Allometrically scaled FFM gave a better model fit (delta-OFV = -8.7) than total body weight (delta-OFV = -6.6) per paper Results page 6.",
      source_name        = "FFM"
    ),
    TCHOL = list(
      description        = "Total serum cholesterol",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear effect on apparent CL/F centered at 3 mmol/L (Archary 2018 Table 1 baseline mean 2.7-2.9 mmol/L). Treated as a surrogate for nutritional / hepatic-function recovery in this severely malnourished pediatric cohort. Paper text describes the effect as 20.7% increase in F (paper Results page 6); the published equation in Table 2 footnote applies the multiplier directly to apparent CL/F -- this model implements the equation as printed in Table 2 (see Assumptions and deviations in the validation vignette).",
      source_name        = "CHOL"
    )
  )

  population <- list(
    n_subjects     = 62L,
    n_studies      = 1L,
    n_observations = 502L,
    age_range      = "1 month to 12 years (eligibility); enrolled cohort 0.1-3.9 years",
    age_median     = "0.9 years",
    weight_range   = "approx 4-12 kg (WHO weight-band dosing range 3-19.9 kg eligible)",
    weight_mean    = "6.5 kg (early-ART arm) / 6.6 kg (delayed-ART arm) per Table 1",
    sex_female_pct = 43,
    race_ethnicity = "African (South African pediatric cohort, Durban / KwaZulu-Natal)",
    disease_state  = "HIV infection with severe acute malnutrition (weight-for-length Z-score < -3, mid-upper arm circumference < 115 mm, or peripheral edema). 20 of 62 patients on rifampicin-based anti-tuberculosis treatment received super-boosted LPV/rtv (LPV:rtv 1:1 ratio).",
    dose_range     = "Oral LPV/rtv (liquid) per WHO weight-band dosage charts. Children 3-5.9 kg received approximately 20/25 mg/kg LPV/rtv twice daily; 14-19.9 kg received approximately 11.7/12 mg/kg twice daily. Super-boosted regimen used in TB co-treated patients.",
    regions        = "South Africa (King Edward VIII Hospital, Durban)",
    notes          = "MATCH study (Malnutrition and ART Timing in Children with HIV; Clinical trial registry PACTR21609001751384). 56 patients had Day-14 sampling. Day 1 sampling schedule: 1.3-1.8 h, 3-4 h, 5-7 h, 8-10 h post-dose. Day 13 trough (pre-dose). Day 14 schedule: 30 min pre-dose plus 1.3-1.8 h, 3-4 h, 5-7 h, 8-10 h post-dose. LPV LoQ 0.0195 ug/mL. Patients randomized to early (within 14 days of admission) vs delayed (until nutritional recovery to WHZ -2 plus 14 days from admission) ART initiation."
  )

  ini({
    # Structural parameters from Archary 2018 Table 2 (final model, second column "Parameter estimates").
    # Reported in units of L/h/5.6 kg and L/5.6 kg (the FFM reference); the paper writes
    # the structural formulas with the (FFM/5.6) ratio so the typical-value parameters
    # apply directly at FFM = 5.6 kg.
    lcl <- log(3.1);  label("Apparent clearance at FFM=5.6 kg, TCHOL=3 mmol/L (CL/F, L/h)")  # Table 2 row 1: CL/F = 3.1 L/h/5.6 kg
    lvc <- log(9.6);  label("Apparent volume of distribution at FFM=5.6 kg (Vd/F, L)")        # Table 2 row 2: Vd/F = 9.6 L/5.6 kg
    lka <- log(0.39); label("First-order absorption rate constant (1/h)")                     # Table 2 row 3: ka = 0.39 /h

    # Allometric exponents fixed per Archary 2018 Methods page 4 ("Allometric exponents
    # were fixed to 0.75 for CL/F and 1 for Vd/F, when testing use of weight or fat-free
    # mass (FFM) to account for size").
    e_ffm_cl <- fixed(0.75); label("Allometric exponent on CL with FFM (unitless, fixed)")  # Methods page 4
    e_ffm_vc <- fixed(1.00); label("Allometric exponent on Vd with FFM (unitless, fixed)")  # Methods page 4

    # Cholesterol linear effect on apparent CL/F (paper Table 2 footnote equation):
    #   CL = 3.1 * (FFM/5.6)^0.75 * (1 + 0.207 * (CHOL - 3))
    # Paper text describes this as "20.7% increase in F per 1 mmol/L above 3 mmol/L"
    # (Results page 6); see vignette Assumptions and deviations for direction note.
    e_tchol_cl <- 0.207; label("Linear coefficient on (TCHOL - 3 mmol/L) for apparent CL/F (per mmol/L)")  # Table 2 footnote

    # Inter-individual variability. Archary 2018 reports the only IIV term on relative
    # bioavailability F (logit-transformed, IIV = 69.5% per Table 2 row 4 / Results page
    # 5: "The IIV estimated for F consequently is reflective of variability for apparent
    # CL/F and Vd/F estimates"). Because Cmax / AUC depend on CL/F = CL_true/F and
    # Vd/F = Vd_true/F, F-IIV propagates to apparent CL/F. This model represents the
    # paper's IIV as IIV on log(CL/F): omega^2 = log(1 + 0.695^2) per the standard
    # log-normal CV-to-variance conversion (vignette documents the simplification).
    etalcl ~ 0.395  # log(1 + 0.695^2) = 0.3946; CV = 69.5% per Table 2

    # Residual error. Archary 2018 reports a piecewise proportional RUV (Table 2 rows 6-7
    # / Results page 5-6): 37.7% within the first 5 h post-dose and 27.2% at 5+ h (with
    # an additional 15.5% BSV on the RUV magnitude). This model uses a single
    # proportional error fixed to the larger 37.7% value as a conservative envelope; the
    # vignette documents the simplification.
    propSd <- 0.377; label("Proportional residual error (CV, fraction)")  # Table 2 row 6 (<5 h post-dose)
  })

  model({
    # Individual PK parameters. Apparent CL/F is allometrically scaled with FFM
    # (reference 5.6 kg) and modified linearly by total cholesterol per the equation
    # in Archary 2018 Table 2 footnote.
    cl <- exp(lcl + etalcl) * (FFM / 5.6)^e_ffm_cl * (1 + e_tchol_cl * (TCHOL - 3))
    vc <- exp(lvc) * (FFM / 5.6)^e_ffm_vc
    ka <- exp(lka)

    # Micro-constant
    kel <- cl / vc

    # ODE system (1-compartment with first-order absorption from depot)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation and proportional residual error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

Bista_2015_fentanyl <- function() {
  description <- "One-compartment population PK model for transdermal fentanyl (Durogesic patch) in adult cancer patients with first-order absorption from the patch and allometric body-weight scaling on CL/F and V/F (Bista 2015)"
  reference <- paste(
    "Bista SR, Haywood A, Hardy J, Norris R, Hennig S.",
    "Exposure to fentanyl after transdermal patch administration",
    "for cancer pain management.",
    "Manuscript dated 2015 provided by the senior author (S. Hennig);",
    "published-journal citation / DOI not on the manuscript copy used",
    "for extraction."
  )
  vignette <- "Bista_2015_fentanyl"
  units <- list(time = "hr", dosing = "ug", concentration = "ug/L") # Methods + Tables 1 and 3: dose in ug/h, plasma concentration in ug/L

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric a priori scaling per Bista 2015 Table 2 footnote: CL/F = theta1 * (WT/70)^0.75 and V/F = theta2 * (WT/70). Reference weight 70 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 56L,                                # Bista 2015 Abstract / Table 1
    n_studies      = 1L,                                 # Single-centre observational study (Brisbane 2011-2014)
    age_range      = "39-90 years",                      # Bista 2015 Table 1
    age_median     = "69.5 years",                       # Bista 2015 Table 1
    weight_range   = "41.8-110.0 kg",                    # Bista 2015 Table 1
    weight_median  = "71.5 kg",                          # Bista 2015 Table 1
    sex_female_pct = 39.3,                               # 22 female / 56 total per Bista 2015 Table 1 (34/22 male/female)
    race_ethnicity = NULL,                               # Not reported in Bista 2015 Table 1
    disease_state  = "Adults with advanced malignant disease receiving Durogesic transdermal fentanyl matrix patches for cancer-pain management; common diagnoses included ovaries, prostate, breast, cervix, lung and bone (Bista 2015 Table 1).",
    dose_range     = "12-200 ug/h transdermal Durogesic patch (median 50 ug/h); patches replaced every 72 h; treatment duration variable, time since last patch change at sampling 0.5-77 h.",
    regions        = "Australia (single tertiary cancer centre, Brisbane).",
    sampling       = "Sparse opportunistic sampling: 163 plasma samples from 56 patients, median 2 samples/patient (range 1-10) on median 2 occasions/patient (range 1-10).",
    co_medication  = "Concomitant CYP3A4/3A5 inhibitors and inducers were collected; 24 patients on inducer only, 2 on inhibitor only, 5 on both; chloramphenicol (n=1), diltiazem (n=1), fluconazole (n=2), tamoxifen (n=2), fluoxetine (n=1) as inhibitors and dexamethasone (n=28) as inducer (Bista 2015 Table 1, Methods).",
    notes          = "BSA 0.06-2.26 m^2 (median 1.80); BMI 15-41 kg/m^2 (median 24.8); creatinine clearance 25.4-216.0 mL/min (median 99.0); ALT 8-153 IU/L (median 35.0); AST 15.0-169.0 IU/L (median 40.7); ALP 24.1-2246 IU/L (median 138). Pain score 0-10 (median 2); patch adhesion 0-3 (median 0). None of the tested covariates beyond a priori weight scaling were retained in the final model."
  )

  ini({
    # Structural PK parameters - reference weight 70 kg (Bista 2015 Table 2 footnote)
    lka  <- log(0.013); label("Absorption rate constant (1/h)")                # Bista 2015 Table 3 final model: ka = 0.013 (RSE 21.1%); 90% bootstrap CI 0.008-0.018
    lcl  <- log(122);   label("Apparent clearance CL/F at 70 kg (L/h)")        # Bista 2015 Table 3 final model: CL/F = 122 L/h/70kg (RSE 9.4%); 90% bootstrap CI 104.9-142.7
    lvc  <- fix(log(350)); label("Apparent volume of distribution V/F at 70 kg (L)") # Bista 2015 Table 3 + Table 2 footnote: V/F fixed to 350 L/70kg (Janssen Durogesic Product Information, ref 26 of Bista 2015)

    # Allometric scaling exponents - fixed a priori per Bista 2015 Table 2 footnote
    e_wt_cl <- fix(0.75); label("Body-weight allometric exponent on CL/F")     # Bista 2015 Table 2 footnote: CL/F = theta1 * (WT/70)^0.75 (a priori, theoretical allometric)
    e_wt_vc <- fix(1);    label("Body-weight allometric exponent on V/F")      # Bista 2015 Table 2 footnote: V/F = theta2 * (WT/70) (a priori, linear weight)

    # IIV - Bista 2015 Table 3 reports CV% on the log-normal scale
    # omega^2 = log(1 + CV^2) (per the standard log-normal CV-to-variance identity)
    etalcl ~ log(1 + 0.385^2) # Bista 2015 Table 3 final model: BSV CL = 38.5% CV (RSE 19.5%); 90% bootstrap CI 24.9-50.0%

    # Residual error - Bista 2015 Table 3 final model: proportional 36.3%
    propSd <- 0.363; label("Proportional residual error (fraction)") # Bista 2015 Table 3 final model: proportional error 36.3% (RSE 10.2%); 90% bootstrap CI 21.9-40.8%
  })
  model({
    # Individual structural parameters with allometric weight scaling
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl   # Bista 2015 Table 2 footnote: CL/F = theta1 * (WT/70)^0.75
    vc <- exp(lvc) * (WT / 70)^e_wt_vc            # Bista 2015 Table 2 footnote: V/F = theta2 * (WT/70); V/F fixed, no IIV

    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}

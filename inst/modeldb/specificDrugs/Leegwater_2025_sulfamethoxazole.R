Leegwater_2025_sulfamethoxazole <- function() {
  description <- paste(
    "Integrated parent-metabolite population PK model for sulfamethoxazole",
    "and its major metabolite N-acetyl sulfamethoxazole in hospitalized",
    "adults treated with oral or intravenous cotrimoxazole",
    "(trimethoprim/sulfamethoxazole), built from routine",
    "therapeutic-drug-monitoring data in three Dutch university medical",
    "centers (Leegwater 2025). One compartment for each analyte with",
    "first-order absorption into the sulfamethoxazole compartment; the",
    "metabolite formation clearance is fixed at 0.4 times the",
    "sulfamethoxazole elimination clearance because no urine data were",
    "available to identify it. BSA-normalized eGFR enters both clearances as",
    "a power function, much more steeply for the renally cleared metabolite",
    "(exponent 0.797) than for the largely hepatically metabolized parent",
    "(0.27), which is why N-acetyl sulfamethoxazole accumulates in renal",
    "impairment. Continuous renal replacement therapy replaces the eGFR term",
    "outright, raising sulfamethoxazole clearance 2.2-fold while LOWERING",
    "metabolite clearance to 0.68-fold. Oral bioavailability was estimated at",
    "~100% and then fixed to 1, so clearances and volumes are apparent values",
    "that apply to both routes. Dose amounts are the SULFAMETHOXAZOLE",
    "component of the combination product (cotrimoxazole 1,920 mg delivers",
    "1,600 mg sulfamethoxazole).",
    sep = " "
  )
  reference <- paste(
    "Leegwater E, Baidjoe L, Wilms EB, Visser LG, Touw DJT, de Winter BCM,",
    "de Boer MGJ, van Paassen J, van den Berg CHSB, van Prehn J,",
    "van Gelder T, Moes DJAR. Population Pharmacokinetics of",
    "Trimethoprim/Sulfamethoxazole: Dosage Optimization for Patients with",
    "Renal Insufficiency or Receiving Continuous Renal Replacement Therapy.",
    "Clin Pharmacol Ther. 2025;117(1):184-192. doi:10.1002/cpt.3421.",
    "Fixed effects, between-subject variability and residual error from",
    "Table 3; the rate-constant structure (K20 / K23 / K30), the piecewise",
    "CRRT-vs-eGFR clearance branches and the OMEGA variances from the",
    "sulfamethoxazole NONMEM control stream reproduced in the Supplement",
    "('Supplement NONMEM code', ADVAN5 TRANS1). The control stream also",
    "corrects a transcription error in the Table 3 footnote -- see the",
    "e_rrt_crrt_status_cl comment below and the vignette Errata.",
    sep = " "
  )
  vignette <- "Leegwater_2025_cotrimoxazole"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate calculated with the CKD-EPI equation and reported BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters both apparent clearances as a power function normalized to",
        "68 mL/min/1.73 m^2 (Table 3 footnote a; control stream",
        "'(EGFR/68)**THETA(11)' for the parent and '(EGFR/68)**THETA(10)'",
        "for the metabolite). The normalizer 68 is what the model uses; the",
        "cohort median eGFR quoted in the Results is 70 and the Table 1 mean",
        "is 70.8.",
        "The exponents differ sharply by analyte -- 0.27 for sulfamethoxazole",
        "versus 0.797 for N-acetyl sulfamethoxazole -- which the Discussion",
        "attributes to route: less than 30% of sulfamethoxazole is recovered",
        "unchanged in urine because its primary clearance route is hepatic",
        "metabolism, whereas the metabolite is mainly cleared renally.",
        "The eGFR term is switched OFF entirely for patients on CRRT, for",
        "both analytes. A CRCL value must nevertheless be supplied for every",
        "subject because the arithmetic switch in model() evaluates both",
        "branches; any positive placeholder works for CRRT subjects.",
        sep = " "
      ),
      source_name        = "EGFR"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Continuous renal replacement therapy during cotrimoxazole treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CRRT)",
      notes              = paste(
        "1 = concomitant CRRT, 0 = no CRRT. Called 'rel' in the NONMEM",
        "control stream and 'CRRT' in Table 3. Subject-level and time-fixed:",
        "patients treated with intermittent hemodialysis or ECMO were",
        "excluded, so the flag never changes within a subject. 18 of the 168",
        "subjects (10.7%) were on CRRT (Table 1).",
        "The effect is directionally OPPOSITE for the two analytes -- CRRT",
        "multiplies sulfamethoxazole clearance by 2.2 but multiplies",
        "metabolite clearance by 0.683. The Discussion explains this as loss",
        "of tubular reabsorption: sulfamethoxazole is extensively reabsorbed",
        "in the native kidney and is not reabsorbed once it is filtered into",
        "the ultrafiltrate, whereas N-acetyl sulfamethoxazole is reabsorbed",
        "much less and so is little affected.",
        "In both branches CRRT replaces the eGFR power term rather than",
        "multiplying on top of it.",
        sep = " "
      ),
      source_name        = "rel"
    )
  )

  compartmentData <- list(
    depot         = list(analyte = "sulfamethoxazole", units = "mg", specimen = "administration site", verified = TRUE),
    central       = list(analyte = "sulfamethoxazole", units = "mg", specimen = "plasma", verified = TRUE),
    central_nasmx = list(analyte = "N-acetyl sulfamethoxazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 168,
    n_studies      = 1,
    age_mean       = "58.2 years (SD 15.2); median 58 years",
    weight_mean    = "76.7 kg (SD 16.2)",
    sex_female_pct = 35.7,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Hospitalized adults (>= 18 years) treated with therapeutic doses of",
      "oral or intravenous cotrimoxazole, sampled as part of routine",
      "therapeutic drug monitoring. Indications include Pneumocystis",
      "jirovecii pneumonia and other infections requiring high-dose",
      "cotrimoxazole. 20.8% were solid-organ transplant recipients, 13.7%",
      "had a malignancy, 11.3% were living with HIV, 7.7% had a stem-cell",
      "transplant and 36.9% received corticosteroids.",
      sep = " "
    ),
    renal_function = paste(
      "eGFR (CKD-EPI) mean 70.8 mL/min/1.73 m^2 (SD 33.2), median 70; serum",
      "creatinine mean 119.7 umol/L (SD 91.8). 18 of 168 (10.7%) received",
      "concomitant CRRT. Patients on intermittent hemodialysis or ECMO were",
      "excluded.",
      sep = " "
    ),
    dose_range     = paste(
      "Routine care rather than protocol-assigned. Daily cotrimoxazole",
      "starting doses: <= 960 mg 6.0%, 1,920-2,400 mg 16.1%, 2,880 mg 8.9%,",
      "3,840-4,800 mg 20.8%, 5,760 mg 48.2%. 44.6% started oral, 55.4%",
      "intravenous.",
      sep = " "
    ),
    regions        = "The Netherlands (Leiden University Medical Center n = 16, Erasmus MC n = 116, University Medical Center Groningen n = 36).",
    n_observations = "348 paired sulfamethoxazole and N-acetyl sulfamethoxazole plasma concentrations from 168 patients, peaks and troughs. Sulfamethoxazole ranged 2-380 mg/L and N-acetyl sulfamethoxazole 2-173.4 mg/L.",
    notes          = paste(
      "Retrospective multicenter observational cohort, January 2016 to",
      "December 2021 (Methods, 'Study design' and 'Participants and data",
      "collection'); demographics from Table 1, left-hand column. The",
      "Discussion states this is the only population PK analysis to date",
      "reporting N-acetyl sulfamethoxazole pharmacokinetics. The",
      "trimethoprim model of the same paper was fitted to a smaller",
      "52-patient subset -- see modellib('Leegwater_2025_trimethoprim').",
      sep = " "
    )
  )

  ini({
    # Sulfamethoxazole (parent) structural parameters - Leegwater 2025
    # Table 3, cross-checked against the $THETA block of the supplement
    # control stream (all values there are the final estimates entered as
    # FIX for simulation).
    lka <- log(0.978); label("First-order absorption rate constant for sulfamethoxazole (1/h)")                     # Table 3: 0.978 1/h, RSE 45%, bootstrap 95% CI 0.25-8.61; control stream $THETA(3) "0.978 FIX"
    lcl <- log(0.97); label("Apparent sulfamethoxazole elimination clearance at eGFR 68 mL/min/1.73 m^2 without CRRT (L/h)") # Table 3: 0.97 L/h, RSE 4%; control stream $THETA(1) "0.97 FIX". NOTE this is the K20 (elimination) arm only -- the metabolite formation arm K23 is an ADDITIONAL 0.4*CL, so total sulfamethoxazole elimination is 1.4*CL. See f_clform_nasmx below.
    lvc <- log(37.0); label("Apparent sulfamethoxazole central volume of distribution (L)")                             # Table 3: 37.0 L, RSE 6%, bootstrap 95% CI 21.5-43.3; control stream $THETA(2) "37 FIX"
    lfdepot <- fixed(log(1)); label("Oral bioavailability (unitless)")                                               # Table 3: "Biological availability 1 Fixed"; control stream $THETA(5) "1 FIX". Results: estimated first, found ~100%, then fixed to 100%.

    # Sulfamethoxazole covariate effects. Mutually exclusive branches, not
    # multiplicative layers: control stream $PK evaluates
    # IF(rel==0) TVCL = THETA(1)*((EGFR/68)**THETA(11)) and
    # IF(rel==1) TVCL = THETA(1)*THETA(12).
    e_crcl_cl <- 0.27; label("Power exponent on (CRCL/68) for sulfamethoxazole CL in patients not receiving CRRT (unitless)") # Table 3 "eGFR on CL": 0.27, RSE 20%, bootstrap 95% CI 0.16-0.38; control stream $THETA(11) "0.27 FIX"
    e_rrt_crrt_status_cl <- 2.2; label("Multiplicative factor on sulfamethoxazole CL for patients receiving CRRT (unitless)") # Table 3 "CRRT on CL": 2.2, RSE 10%, bootstrap 95% CI 1.8-2.6; control stream $THETA(12) "2.2 FIX". The Table 3 footnote prints the CRRT branch as "1.34 x 2.2", which is a transcription error: 1.34 is the METABOLITE's typical clearance. The control stream multiplies THETA(1) = 0.97, the Results state CRRT clearance is "2.2 times higher" than at the median eGFR, and Table 4's "200% for sulfamethoxazole" dose requirement is 0.97*2.2/0.978 = 2.2-fold, not the 3.0-fold that 1.34*2.2 would give.

    # N-acetyl sulfamethoxazole (metabolite) parameters. The "_nasmx" suffix
    # is registered as a metabolite-suffix entry in
    # inst/references/compartment-names.md.
    f_clform_nasmx <- fixed(0.4); label("N-acetyl sulfamethoxazole formation clearance as a multiple of sulfamethoxazole CL (unitless)") # Table 3 "Conversion parent metabolite: 0.4 x CL"; control stream "K23 = 0.4*CL/V2". Methods: "No urine concentrations were available for sulfamethoxazole or N-acetyl sulfamethoxazole; therefore, 40% of the total sulfamethoxazole clearance was estimated to be converted to N-acetyl sulfamethoxazole" (citing Kucers' the Use of Antibiotics, ref. 19). Deliberately NOT named fm_nasmx: K23 sits alongside K20 rather than inside it, so 0.4 is the formation clearance relative to the ELIMINATION clearance, and the implied share of TOTAL clearance is 0.4/1.4 = 0.286. See the vignette Errata.
    lcl_nasmx <- log(1.34); label("Apparent N-acetyl sulfamethoxazole clearance at eGFR 68 mL/min/1.73 m^2 without CRRT (L/h)") # Table 3: 1.34 L/h, RSE 4%, bootstrap 95% CI 1.21-1.43; control stream $THETA(4) "1.34 FIX"
    lvc_nasmx <- log(3.98); label("Apparent N-acetyl sulfamethoxazole central volume of distribution (L)")              # Table 3: 3.98 L, RSE 24%, bootstrap 95% CI 1.55-6.33; control stream $THETA(7) "3.98 FIX"
    e_crcl_cl_nasmx <- 0.797; label("Power exponent on (CRCL/68) for N-acetyl sulfamethoxazole CL in patients not receiving CRRT (unitless)") # Table 3 "eGFR on CL" (metabolite block): 0.797, RSE 7%, bootstrap 95% CI 0.65-0.92; control stream $THETA(10) "0.797 FIX"
    e_rrt_crrt_status_cl_nasmx <- 0.683; label("Multiplicative factor on N-acetyl sulfamethoxazole CL for patients receiving CRRT (unitless)") # Table 3 "CRRT on CL" (metabolite block): 0.683, RSE 11%, bootstrap 95% CI 0.54-0.86; control stream $THETA(13) "0.683 FIX"; footnote a "in case of CRRT: 1.34 x 0.683"

    # Between-subject variability. The supplement's $OMEGA block carries
    # variances on the log scale (every parameter enters as TV*EXP(ETA)) and
    # Table 3 reports each as a percentage that is exactly sqrt(omega^2):
    # sqrt(0.132) = 36.3%, sqrt(0.396) = 62.9%, sqrt(0.166) = 40.7%. The
    # percentages are therefore omega itself, and ini() takes omega^2
    # directly from the control stream. IIV on ka, on F and on the
    # metabolite volume was fixed to zero ($OMEGA(3), $OMEGA(5), $OMEGA(7)
    # all "0 FIX"), so those parameters carry no eta.
    etalcl       ~ 0.132                                                                                               # control stream $OMEGA(1) '0.132 FIX'; Table 3 'IIV CL (%)' 36.3, RSE 8%, bootstrap 95% CI 29.6-42.4
    etalvc       ~ 0.396                                                                                               # control stream $OMEGA(2) '0.396 FIX'; Table 3 'IIV Vd (%)' 62.9, RSE 11%, bootstrap 95% CI 37.4-79.4
    etalcl_nasmx ~ 0.166                                                                                               # control stream $OMEGA(4) '0.166 FIX'; Table 3 'IIV CL (%)' (metabolite block) 40.7, RSE 8%, bootstrap 95% CI 34.5-48.3

    # Residual error. Control stream $ERROR, with $SIGMA 1 FIX:
    # Y = IPRED + W*EPS(1)*IPRED with W = THETA(8) for CMT 2 and
    # W = THETA(9) for CMT 3, i.e. a pure proportional error per analyte.
    propSd       <- 0.181; label("Sulfamethoxazole proportional residual error SD (fraction)")                          # Table 3 "Proportional error / Sulfamethoxazole": 0.181, RSE 8%, bootstrap 95% CI 0.144-0.208; control stream $THETA(8) "0.181 FIX"
    propSd_nasmx <- 0.201; label("N-acetyl sulfamethoxazole proportional residual error SD (fraction)")                 # Table 3 "Proportional error / N-acetyl sulfamethoxazole": 0.201, RSE 8%, bootstrap 95% CI 0.155-0.230; control stream $THETA(9) "0.201 FIX"
  })

  model({
    # Normalizing eGFR for both power terms (Table 3 footnote a and the
    # control stream both use 68, not the cohort median of 70).
    ref_crcl <- 68

    ka <- exp(lka)

    # Apparent clearances. The NONMEM IF branches are written here as an
    # arithmetic switch on the binary covariate, with the eGFR power term
    # multiplied by (1 - RRT_CRRT_STATUS) inside the sum rather than raised
    # to it so a CRRT subject never evaluates a fractional power of a
    # possibly-zero eGFR:
    #   RRT_CRRT_STATUS = 0 -> typical value * (CRCL/68)^exponent
    #   RRT_CRRT_STATUS = 1 -> typical value * CRRT factor
    cl <- exp(lcl + etalcl) *
      ((1 - RRT_CRRT_STATUS) * (CRCL / ref_crcl)^e_crcl_cl +
         RRT_CRRT_STATUS * e_rrt_crrt_status_cl)
    cl_nasmx <- exp(lcl_nasmx + etalcl_nasmx) *
      ((1 - RRT_CRRT_STATUS) * (CRCL / ref_crcl)^e_crcl_cl_nasmx +
         RRT_CRRT_STATUS * e_rrt_crrt_status_cl_nasmx)

    vc       <- exp(lvc + etalvc)
    vc_nasmx <- exp(lvc_nasmx)

    # Micro-constants, named for the ADVAN5 rate constants of the control
    # stream: K12 = KA, K20 = CL/V2, K23 = 0.4*CL/V2, K30 = CM/V3.
    #
    # K20 and K23 are BOTH first-order losses from the sulfamethoxazole
    # compartment, so total sulfamethoxazole elimination is
    # (1 + f_clform_nasmx) * cl / vc, i.e. 1.4 * cl / vc. The Methods
    # sentence describing 0.4 as a share of "total" clearance is looser than
    # what the control stream implements; the vignette Errata reproduces the
    # paper's own target-attainment percentages to confirm the 1.4 factor.
    kel       <- cl / vc                          # K20
    kform     <- f_clform_nasmx * cl / vc         # K23
    kel_nasmx <- cl_nasmx / vc_nasmx              # K30

    # The metabolite transfer is a plain amount-for-amount transfer, as
    # ADVAN5's K23 is: the source model applies no molar-mass correction
    # between sulfamethoxazole (253.3 g/mol) and N-acetyl sulfamethoxazole
    # (295.3 g/mol), and the apparent metabolite volume absorbs it.
    d/dt(depot)         <- -ka * depot
    d/dt(central)       <-  ka * depot - kel * central - kform * central
    d/dt(central_nasmx) <-  kform * central - kel_nasmx * central_nasmx

    f(depot) <- exp(lfdepot)

    Cc       <- central / vc
    Cc_nasmx <- central_nasmx / vc_nasmx

    Cc       ~ prop(propSd)
    Cc_nasmx ~ prop(propSd_nasmx)
  })
}

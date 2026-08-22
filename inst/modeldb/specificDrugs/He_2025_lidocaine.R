He_2025_lidocaine <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous lidocaine with two",
    "sequential first-order metabolite compartments (monoethylglycinexylidide",
    "[MEGX] and glycinexylidide [GX]) in Chinese adults undergoing partial",
    "hepatectomy. Lidocaine is given as a short loading infusion at anaesthesia",
    "induction followed by a continuous intraoperative infusion. Total lidocaine",
    "clearance CL is split so that a fraction fm_megx forms MEGX; MEGX is",
    "eliminated with clearance cl_megx, of which a fraction fm_gx forms GX; GX",
    "is eliminated with clearance cl_gx. Both metabolite compartments share the",
    "lidocaine central volume vc, which the authors fixed by assumption because",
    "no metabolite disposition data were available. Covariates retained after",
    "stepwise forward-addition / backward-elimination are tumour size on",
    "lidocaine clearance, total administered lidocaine dose on the",
    "lidocaine-to-MEGX fraction metabolised, on the peripheral volume and on GX",
    "clearance, and total body weight on the lidocaine-to-MEGX fraction",
    "metabolised. The paper reports covariate coefficients as bare power",
    "exponents without stating the normalisation constants; the packaged model",
    "normalises each covariate by its Table 1 study-typical value (see vignette",
    "Errata). The paper's variability columns are headed '%CV' but hold",
    "variances (see vignette Errata for the objective-function falsifier).",
    "Inter-individual variability on the MEGX-to-GX fraction and on MEGX",
    "clearance was fixed to zero by the authors and is therefore absent here."
  )
  reference <- paste(
    "He C, Qi X, Liu Y, Jin Y, Zhang M, Zhang Y, Fu L, Zheng L, Tu F, Wang Z.",
    "Optimizing Lidocaine Dosing in Hepatectomy Patients: A Population",
    "Pharmacokinetic Study of Active Metabolites.",
    "Drug Des Devel Ther. 2025;19:6255-6268. doi:10.2147/DDDT.S485389.",
    "Supplementary material (Dovepress supplementary file 485389,",
    "Tables S1-S4) supplies the covariate-search and sensitivity-analysis",
    "tables used for the variance-scale falsifier.",
    sep = " "
  )
  vignette <- "He_2025_lidocaine"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    TUMSZ = list(
      description        = "Preoperative liver tumour size, measured by imaging as the largest lesion diameter. Enters lidocaine clearance as a power term (TUMSZ / 50 mm)^e_tumsz_cl with a negative exponent, so clearance falls as tumour size rises; the authors read this as an indirect surrogate for the remaining functional liver mass (Discussion, p. 6266).",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column `SIZE`, reported in cm (Table 1 median 5 cm, range 1.5-15 cm). Converted to the canonical mm unit on ingestion (median 50 mm, range 15-150 mm) and the model's normalisation constant scaled to match, so (TUMSZ / 50)^exponent is numerically identical to (SIZE_cm / 5)^exponent. The normalisation constant is NOT stated in the paper; 50 mm is the Table 1 median. Single largest-lesion diameter, not a RECIST sum of diameters, hence TUMSZ rather than TUM_SLD.",
      source_name        = "SIZE"
    ),
    DOSE_LIDOCAINE_MG = list(
      description        = "Total lidocaine dose administered to the subject over the whole treatment episode (loading infusion plus continuous intraoperative infusion), in mg. Time-fixed per subject. Enters three parameters as power terms in (DOSE_LIDOCAINE_MG / 312 mg): the lidocaine-to-MEGX fraction metabolised, the peripheral volume, and GX clearance.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column `DOSE` (Table S2 abbreviation list: 'DOSE, Total doses of lidocaine'). A drug-specific member of the `DOSE_<drug>_<units>` family is used rather than the bare `DOSE` canonical because rxode2's event-table translator (`etTrans()`) consumes a column literally named `DOSE` and never exposes it to `model()`; a model reading this covariate from an event table fails at solve time with 'The following parameter(s) are required for solving: DOSE'. Confirmed for this model. The normalisation constant is NOT stated in the paper; 312 mg is derived from Table 1 as the mean short-infusion (loading) dose 86.07 mg plus the median continuous infusion rate 57.97 mg/h times the mean continuous infusion time 3.90 h. See vignette Errata for the reproduction of the paper's observed lidocaine and metabolite concentrations that supports this value.",
      source_name        = "DOSE"
    ),
    WT = list(
      description        = "Total body weight at baseline. Enters the lidocaine-to-MEGX fraction metabolised as a power term (WT / 60 kg)^e_wt_cl_fm_megx with a negative exponent.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column `TBW` (Table S2 abbreviation list: 'TBW, Total body weight'). Table 1 reports mean 60.32 kg (SD 8.27). The normalisation constant is NOT stated in the paper; 60 kg is used, matching both the Table 1 mean and the 'typical 60 kg, 55-year-old individual' the authors simulate in Results / Simulation and Weighted Lidocaine Exposure. Distinct from `IBW` (ideal body weight), which the paper also records and uses to express mg/kg dosing but which did NOT survive covariate selection.",
      source_name        = "TBW"
    )
  )

  covariatesDataExcluded <- list(
    IBW = list(
      description = "Ideal body weight. Screened on lidocaine clearance and on the central and peripheral volumes (Table S2 models 5, 35, 44) and carried through forward addition (models 48, 67, 82, 96, 107) but removed at backward elimination (Table S2 model 112, dOFV +6.532, not significant at p < 0.01).",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained in the final model; no point estimate is reported for it."
    ),
    HEIGHT = list(
      description = "Body height. Screened on lidocaine clearance and peripheral volume (Table S2 models 3, 42, 46, 65, 80, 94); never reached the final model.",
      units       = "m",
      type        = "continuous",
      notes       = "Screened but not retained; source column `HEI`."
    ),
    HGB = list(
      description = "Haemoglobin concentration. Screened on lidocaine clearance, MEGX clearance, the fraction metabolised, and both volumes (Table S2 models 4, 16, 28, 34, 43, 47, 53, 66, 71, 81, 86, 95, 106, 110); never reached the final model.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase activity. Screened on inter-compartmental clearance, GX clearance and both volumes (Table S2 models 9, 19, 32; Table S3 sensitivity models); never significant.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    AGE = list(
      description = "Age at surgery. Screened on the central volume (Table S2 model 31, dOFV -2.63, not significant).",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    APTT = list(
      description = "Activated partial thromboplastin time. Screened on MEGX clearance and the fraction metabolised (Table S2 models 12, 25, 50, 58, 68, 83, 97, 108) and entered forward addition at model 111, but removed at backward elimination (model 114, dOFV +4.624, not significant at p < 0.01).",
      units       = "s",
      type        = "continuous",
      notes       = "Screened but not retained; no point estimate is reported for it."
    ),
    INR = list(
      description = "International normalised ratio of prothrombin time. Screened on MEGX and GX clearance (Table S2 models 17, 22, 54, 56, 72, 74, 87, 88, 101); never reached the final model.",
      units       = "(ratio)",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    PT = list(
      description = "Prothrombin time. Screened on MEGX and GX clearance (Table S2 models 18, 23, 57, 75, 89); never reached the final model.",
      units       = "s",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    WBC = list(
      description = "White blood cell count. Screened on inter-compartmental clearance (Table S2 model 11, dOFV -2.215, not significant).",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    CIRR = list(
      description = "Presence or absence of liver cirrhosis (Child-Pugh grading). Screened on inter-compartmental clearance, MEGX clearance, GX clearance, the central volume and the peripheral volume (Table S2 models 10, 13, 20, 33, 39, 51, 63, 69, 84, 98, 109); never reached the final model.",
      units       = "(binary indicator)",
      type        = "categorical",
      notes       = "Screened but not retained. 32 of 33 graded patients were Child-Pugh A and 1 was Child-Pugh B (Table 1), so the split carries very little information."
    ),
    TYPE = list(
      description = "Extent of liver resection, segment (26 patients) versus lobe (9 patients) resection. Screened on lidocaine clearance, GX clearance and both volumes (Table S2 models 8, 24, 38; Table S3 sensitivity models); never significant.",
      units       = "(binary indicator)",
      type        = "categorical",
      notes       = "Screened but not retained."
    ),
    DURA = list(
      description = "Total duration of Pringle manoeuvres during the operation. Screened on lidocaine clearance, MEGX clearance, the fraction metabolised and the peripheral volume (Table S2 models 2, 15, 27, 41, 60, 64, 77, 91, 102); never reached the final model.",
      units       = "h",
      type        = "continuous",
      notes       = "Screened but not retained. The paper nevertheless fixes it at 0.63 h when describing the 'typical' simulated individual (Results / Simulation and Weighted Lidocaine Exposure); because it is not in the final model it has no effect on any prediction."
    )
  )

  compartmentData <- list(
    central      = list(analyte = "lidocaine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1  = list(analyte = "lidocaine", units = "mg", specimen = "plasma", verified = TRUE),
    central_megx = list(analyte = "monoethylglycinexylidide (MEGX)", units = "mg", specimen = "plasma", verified = TRUE),
    central_gx   = list(analyte = "glycinexylidide (GX)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 35L,
    n_studies      = 1L,
    age_range      = "20-65 years",
    age_median     = "54 years (Table 1; the Results text states 53 years -- the source is internally inconsistent by one year, and neither value enters the model since age was screened but not retained)",
    weight_range   = "not reported as a range; mean 60.32 kg (SD 8.27)",
    weight_median  = "60.32 kg (mean)",
    sex_female_pct = 20,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Partial hepatectomy for an initial diagnosis of liver tumour: 27 hepatocellular carcinoma, 5 cholangiocarcinoma, 2 benign focal hyperplasia, 1 metastatic liver cancer. Child-Pugh A in 32 of 33 graded patients, Child-Pugh B in 1. Lobe resection in 9 (26%), segment resection in 26 (74%). Median tumour size 5 cm (range 1.5-15 cm); median duration of Pringle manoeuvres 0.62 h (range 0.20-1.84 h). Severe hepatic insufficiency (total bilirubin > 2 mg/dL), severe renal dysfunction (GFR < 30 mL/min/1.73 m^2), severe cardiac disease, ideal body weight < 40 kg and any long-term co-medication affecting lidocaine metabolism were exclusion criteria.",
    dose_range     = "Intravenous loading infusion of 1.5 mg/kg ideal body weight given over more than 10 min at anaesthesia induction (mean 86.07 mg, SD 9.57), followed by a continuous intraoperative infusion targeting 1.0 mg/kg/h ideal body weight (median achieved rate 0.99 mg/kg/h, range 0.62-1.54; median 57.97 mg/h, range 38.83-100.00) for a mean of 3.90 h (SD 1.62).",
    regions        = "China (single centre, Affiliated Hospital of North Sichuan Medical College, Nanchong, Sichuan)",
    bmi_range      = "18-30 kg/m^2 by inclusion criterion; observed mean 22.87 kg/m^2 (SD 2.67)",
    notes          = "Demographics from Table 1. Single-centre, prospective, open-label study conducted January-December 2021; Chinese Clinical Trial Registry ChiCTR2100042730. Plasma lidocaine, MEGX and GX were measured by UPLC-MS/MS with LLOQs of 10, 2 and 2 ng/mL and calibration ranges of 10-5000, 2-1000 and 2-500 ng/mL respectively. Samples were drawn at baseline, 0.5 h and every 1 h after the start of surgery, and at 0, 0.5, 1, 2, 4, 8 and 12 h after the end of surgery. Estimation was FOCE with interaction in NONMEM 7.2.0. Table S2 traces the stepwise search from a base model at OFV 10208.285 to the last forward-addition model 111 at OFV 10144.856; backward elimination (Table S2 models 112-118) then dropped IBW on CL and APTT on CLE, and the OFV of the resulting final model is not tabulated. The five covariate effects reported in Table 2 -- SIZE on CL, DOSE on CL_FM, TBW on CL_FM, DOSE on V2 and DOSE on CLG -- are what this model encodes."
  )

  ini({
    # ------------------------------------------------------------------
    # Lidocaine structural disposition. Table 2 'Estimate (95% CI)'
    # column. Two-compartment model with first-order elimination from the
    # central compartment (Results / Population Pharmacokinetic Modeling;
    # schematic Figure 2).
    lcl <- log(26.1)
    label("Log typical lidocaine clearance CL at TUMSZ = 50 mm (L/h)")  # Table 2 row 'CL (L/h)' = 26.1 (95% CI 17.59, 34.61; RSE 16.6%)
    lvc <- log(8.73)
    label("Log lidocaine apparent central volume V1 (L)")  # Table 2 row 'V1 (L)' = 8.73 (95% CI 5.34, 12.12; RSE 19.8%)
    lvp <- log(63.6)
    label("Log lidocaine apparent peripheral volume V2 at DOSE_LIDOCAINE_MG = 312 mg (L)")  # Table 2 row 'V2 (L)' = 63.6 (95% CI 43.22, 83.94; RSE 16.4%)
    lq <- log(41.0)
    label("Log lidocaine inter-compartmental clearance CLD (L/h)")  # Table 2 row 'CLD (L/h)' = 41.0 (95% CI 17.87, 64.13; RSE 28.8%)

    # ------------------------------------------------------------------
    # Metabolite chain. Figure 2 draws CL_FM as the lidocaine -> MEGX
    # branch of the lidocaine clearance and CLE_FM as the MEGX -> GX
    # branch of the MEGX clearance; Table 2 defines both as "the
    # proportion of clearance for fraction metabolism", so each is a
    # unitless fraction multiplying the clearance it splits, not a
    # clearance in its own right.
    #
    # Confirmed arithmetically by the Discussion (p. 6266): "the
    # estimated rate of lidocaine conversion to MEGX was 0.475 L/h in
    # this study" = CL_FM * CL = 0.0182 * 26.1 = 0.4750 L/h.
    #
    # Encoded on the log scale (not the logit scale) because the paper's
    # covariate model is multiplicative-power and its inter-individual
    # variability is exponential (Methods: theta_i = theta_TV * exp(eta_i)),
    # both of which are linear on the log scale.
    lfm_megx <- log(0.0182)
    label("Log fraction of lidocaine clearance forming MEGX at DOSE_LIDOCAINE_MG = 312 mg and WT = 60 kg (unitless)")  # Table 2 row 'CL FM' = 0.0182 (95% CI 0.004, 0.032; RSE 39.6%)
    lcl_megx <- log(1.41)
    label("Log MEGX elimination clearance CLE (L/h)")  # Table 2 row 'CLE (L/h)' = 1.41 (95% CI 0.795, 2.025; RSE 22.3%)
    lfm_gx <- log(0.897)
    label("Log fraction of MEGX clearance forming GX (unitless)")  # Table 2 row 'CLE FM' = 0.897 (95% CI 0.794, 1; RSE 5.9%)
    lcl_gx <- log(4.77)
    label("Log GX elimination clearance CLG at DOSE_LIDOCAINE_MG = 312 mg (L/h)")  # Table 2 row 'CLG (L/h)' = 4.77 (95% CI 1.54, 8; RSE 34.6%)

    # ------------------------------------------------------------------
    # Covariate effects. Every retained covariate coefficient in Table 2
    # is a dimensionless power exponent. The paper prints the exponents
    # but never prints the covariate equations or the normalisation
    # constants; the packaged model applies the standard power form
    # parameter = typical * (covariate / reference)^exponent with each
    # reference taken from Table 1 (tumour size median 5 cm = 50 mm;
    # total body weight mean 60.32 kg, rounded to the 60 kg the authors
    # themselves use for their 'typical' simulated individual; total dose
    # 86.07 + 57.97 * 3.90 = 312 mg). See vignette Errata.
    e_tumsz_cl <- -0.382
    label("Power exponent on (TUMSZ / 50 mm) for lidocaine clearance (unitless)")  # Table 2 row 'The effect of SIZE on CL' = -0.382 (95% CI -0.717, -0.047; RSE 44.8%)
    e_dose_fm_megx <- 0.669
    label("Power exponent on (DOSE_LIDOCAINE_MG / 312 mg) for the lidocaine-to-MEGX fraction metabolised (unitless)")  # Table 2 row 'The effect of DOSE on CL FM' = 0.669 (95% CI 0.161, 1.177; RSE 38.7%)
    e_wt_fm_megx <- -1.09
    label("Power exponent on (WT / 60 kg) for the lidocaine-to-MEGX fraction metabolised (unitless)")  # Table 2 row 'The effect of TBW on CL FM' = -1.09 (95% CI -2.66, 0.48; RSE 73.5%)
    e_dose_vp <- 1.27
    label("Power exponent on (DOSE_LIDOCAINE_MG / 312 mg) for the lidocaine peripheral volume (unitless)")  # Table 2 row 'The effect of DOSE on V2' = 1.27 (95% CI 0.315, 2.225; RSE 38.3%)
    e_dose_cl_gx <- 1
    label("Power exponent on (DOSE_LIDOCAINE_MG / 312 mg) for GX clearance (unitless)")  # Table 2 row 'The effect of DOSE on CLG' = 1 (95% CI -0.196, 2.196; RSE 61%)

    # ------------------------------------------------------------------
    # Inter-individual variability. Methods state a log-normal
    # exponential model, theta_i = theta_TV * exp(eta_i).
    #
    # SCALE: Table 2's block is headed 'Inter-Individual Variability
    # (%CV)' but the numbers are omega^2 * 100, i.e. VARIANCES, not
    # percent CVs. The falsifier is the objective-function change of a
    # covariate addition, which for one random effect over N subjects is
    # dOFV = N * log(omega^2_after / omega^2_before) with N = 35:
    #   SIZE on CL   Table S3 13.1 -> 9.55 : variance reading gives
    #                35*log(13.1/9.55)  = -11.06 vs Table S2 model 6
    #                dOFV -10.926;  %CV reading gives -21.99.
    #   SIZE on FM   Table S3 11.63 -> 10.1: variance reading gives
    #                35*log(11.63/10.1) =  -4.94 vs Table S2 model 29
    #                dOFV  -4.525;  %CV reading gives  -9.82.
    #   TYPE on CL   Table S3 13.1 -> 12.6 : variance reading gives
    #                35*log(13.1/12.6)  =  -1.36 vs Table S2 model 8
    #                dOFV  -1.016;  %CV reading gives  -2.66.
    # The variance reading matches all three; the %CV reading is off by
    # roughly a factor of two in each. Values below are therefore
    # Table 2 / 100. See vignette Errata.
    etalcl      ~ 0.0916  # Table 2 IIV row 'CL'    = 9.16 (RSE 57.4%); omega^2 = 9.16/100; equals 30.9% CV
    etalvc      ~ 0.041   # Table 2 IIV row 'V1'    = 4.1  (RSE 44.6%); omega^2 = 4.1/100;  equals 20.5% CV
    etalfm_megx ~ 0.0908  # Table 2 IIV row 'CL FM' = 9.08 (RSE 43.2%); omega^2 = 9.08/100; equals 30.8% CV
    etalvp      ~ 0.194   # Table 2 IIV row 'V2'    = 19.4 (RSE 49.8%); omega^2 = 19.4/100; equals 46.3% CV
    etalq       ~ 0.523   # Table 2 IIV row 'CL D'  = 52.3 (RSE 47.4%); omega^2 = 52.3/100; equals 84.7% CV
    etalcl_gx   ~ 0.16    # Table 2 IIV row 'CLG'   = 16   (RSE 44.9%); omega^2 = 16/100;   equals 41.7% CV

    # No inter-individual variability on the MEGX-to-GX fraction
    # metabolised or on MEGX elimination clearance: Table 2 reports '0*'
    # for both with the footnote '*Fixed to 0'. Omitted rather than
    # written as `~ fixed(0)` because a zero-variance diagonal makes
    # OMEGA singular and breaks the Cholesky sampler used by rxSolve.

    # ------------------------------------------------------------------
    # Residual error. Methods specify a proportional model,
    # Y_ij = F_ij * (1 + eps_ij), separately for each of the three
    # analytes. Table 2's residual rows are labelled 'sigma 2' (a
    # variance) under the same mis-headed '(%CV)' column, and take the
    # same /100 scale as the IIV block above; the standard deviation is
    # therefore sqrt(value / 100).
    propSd <- sqrt(0.054)
    label("Proportional residual error SD on lidocaine concentration (fraction)")  # Table 2 row 'sigma 2 Lidocaine' = 5.4 (RSE 11.6%); variance 0.054 -> SD 23.2%
    propSd_megx <- sqrt(0.036)
    label("Proportional residual error SD on MEGX concentration (fraction)")  # Table 2 row 'sigma 2 MEGX' = 3.6 (RSE 7%); variance 0.036 -> SD 19.0%
    propSd_gx <- sqrt(0.11)
    label("Proportional residual error SD on GX concentration (fraction)")  # Table 2 row 'sigma 2 GX' = 11 (RSE 19.3%); variance 0.11 -> SD 33.2%
  })

  model({
    # 1. Individual parameters. Power covariate models on the normalised
    #    covariates; exponential (log-normal) inter-individual variability.
    cl <- exp(lcl + e_tumsz_cl * log(TUMSZ / 50) + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + e_dose_vp * log(DOSE_LIDOCAINE_MG / 312) + etalvp)
    q  <- exp(lq + etalq)

    # Fraction of lidocaine clearance routed to MEGX, and fraction of
    # MEGX clearance routed to GX. fm_gx carries neither covariates nor
    # IIV in the source.
    fm_megx <- exp(lfm_megx + e_dose_fm_megx * log(DOSE_LIDOCAINE_MG / 312) +
                     e_wt_fm_megx * log(WT / 60) + etalfm_megx)
    fm_gx   <- exp(lfm_gx)

    # Metabolite elimination clearances. cl_megx carries no covariates
    # and no IIV in the source.
    cl_megx <- exp(lcl_megx)
    cl_gx   <- exp(lcl_gx + e_dose_cl_gx * log(DOSE_LIDOCAINE_MG / 312) + etalcl_gx)

    # 2. Micro-constants.
    kel     <- cl / vc
    k12     <- q / vc
    k21     <- q / vp
    kel_megx <- cl_megx / vc
    kel_gx   <- cl_gx / vc

    # 3. ODE system (Figure 2). Lidocaine is dosed into `central`.
    #    Total lidocaine elimination is cl; the fraction fm_megx of that
    #    flux appears as MEGX. Total MEGX elimination is cl_megx; the
    #    fraction fm_gx of that flux appears as GX. Both metabolite
    #    compartments are scaled by the lidocaine central volume vc --
    #    Methods: "The distribution volume of metabolites was set equal
    #    to that of central lidocaine", Table 2 footnote: "Assuming that
    #    the apparent distribution volume (V1) of MEGX and GX were equal
    #    to central compartment (lidocaine)".
    d/dt(central)      <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(central_megx) <-  fm_megx * kel * central - kel_megx * central_megx
    d/dt(central_gx)   <-  fm_gx * kel_megx * central_megx - kel_gx * central_gx

    # 4. Observations. All three analytes are read out of the shared
    #    central volume. The metabolite states carry lidocaine-equivalent
    #    mass: the paper reports no molar conversion between lidocaine
    #    (MW 234.3), MEGX (MW 206.3) and GX (MW 178.2), and sums the three
    #    concentrations directly in its weighted-exposure metric
    #    (lidocaine + MEGX + 0.25 * GX), so the fractions metabolised are
    #    apparent mass fractions. See vignette Errata.
    Cc      <- central / vc
    Cc_megx <- central_megx / vc
    Cc_gx   <- central_gx / vc

    Cc      ~ prop(propSd)
    Cc_megx ~ prop(propSd_megx)
    Cc_gx   ~ prop(propSd_gx)
  })
}

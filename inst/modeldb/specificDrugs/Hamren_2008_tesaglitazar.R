Hamren_2008_tesaglitazar <- function() {
  description <- "Mechanistic parent + acyl-glucuronide population PK model for tesaglitazar (a dual PPAR alpha/gamma agonist) in 41 adult subjects with varying degrees of renal function (Hamren 2008). Parent tesaglitazar follows a two-compartment disposition with first-order oral absorption (ka fixed at 1.5 1/h, F fixed at 1); renal clearance CLrt = 0.027 L/h directs parent to a cumulative urine compartment, and metabolic clearance CLmt = 1.91 L/h generates the acyl glucuronide metabolite. The metabolite follows a one-compartment disposition (Vcm = 8.5 L) with saturable Michaelis-Menten renal clearance (Vmax = 0.188 umol/h, Km = 0.041 umol/L) routing to a cumulative urine compartment, linear non-renal clearance (CLnrm = 1.2 L/h), and biliary excretion (kbm = 11.7 1/h) into a paper-specific gut compartment. The gut compartment releases interconverted parent tesaglitazar back into the parent central compartment at rate kicv = 0.79 1/h, completing the futile-cycle interconversion loop that the source paper proposes as the mechanism for increased tesaglitazar exposure in renal-impairment subjects. Covariates: BSA-normalized renal function CRCL (iohexol-clearance-measured GFR, mL/min/1.73 m^2; linear centered slope on CLrt and direct linear normalised scaling on metabolite Vmax), per-subject free fraction FU (% by ultrafiltration; linear centered slope on CLmt), sex SEXF (women have 31% lower CLrt than men), concomitant probenecid CONMED_PROBENECID (75% reduction of both CLrt and metabolite Vmax), and body weight WT (shared centered linear slope on Vct and Vpt). Concentrations are molar (umol/L) and amounts are molar (umol) throughout to match the Michaelis-Menten parameterisation of the acyl-glucuronide renal elimination; the user converts mg-of-tesaglitazar doses to umol using the molecular weight of 408.45 g/mol (1 mg = 2.45 umol)."
  reference   <- paste(
    "Hamren B, Ericsson H, Samuelsson O, Karlsson MO.",
    "Mechanistic modelling of tesaglitazar pharmacokinetic data in subjects",
    "with various degrees of renal function -- evidence of interconversion.",
    "Br J Clin Pharmacol. 2008;65(6):855-863.",
    "doi:10.1111/j.1365-2125.2008.03110.x.",
    sep = " "
  )
  vignette <- "Hamren_2008_tesaglitazar"
  units    <- list(time = "hour", dosing = "umol", concentration = "umol/L")

  paper_specific_compartments <- c("gut_gluc")

  covariateData <- list(
    CRCL = list(
      description        = "BSA-normalized glomerular filtration rate measured by plasma iohexol clearance (mL/min/1.73 m^2). Iohexol is an exogenous contrast agent cleared exclusively by glomerular filtration; its plasma clearance is the clinical gold standard for measured GFR. Distinct from a creatinine-based estimate (eGFR or measured CrCl); the canonical CRCL register covers both creatinine-based and tracer-measured renal function as of the 2026-06-17 register update.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hamren 2008 Methods: 'During the study (days 1 and 42), plasma iohexol clearance (CLiohexol) was determined and used as a measure of renal function.' Reference 76 mL/min/1.73 m^2 is the pooled-cohort median observed CLiohexol (subjects with renal insufficiency median 32 mL/min/1.73 m^2 [range 16-94]; healthy controls median 90 mL/min/1.73 m^2 [range 75-120], Table 1). Enters the model two ways inside model(): (a) a linear centered slope on parent CLrt with 1 + 0.0099 * (CRCL - 76) (Table 2 'GFR on CLrt 0.99 %/(ml min-1 1.73m-2)'); (b) a direct linear normalised scaling on metabolite Vmax with (CRCL / 76) (Table 2 footnote 'Centred on a male subject with GFR 76 ml min-1 1.73m-2'). Per the standing tracer-measured-GFR-policy decision recorded in this task's sidecar response 001 (2026-06-17), CLiohexol is encoded under the broadened CRCL canonical rather than a separate GFR_TRACER canonical.",
      source_name        = "CLiohexol"
    ),
    FU = list(
      description        = "Per-subject measured fraction unbound of tesaglitazar in plasma at day 42, by ultrafiltration (%).",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hamren 2008 Results 'Pharmacokinetic analysis': 'Median (and range) of tesaglitazar fraction unbound at day 42 was 0.11% (0.08-0.2) and 0.09% (0.06-0.12) for subjects with IRF and controls, respectively. Fraction unbound data were missing for three individuals; for these, the median value (0.1%) was imputed.' Reference 0.1% is used in the linear centered effect 1 + 5.55 * (FU - 0.1) on the parent metabolic clearance CLmt. The realised cohort fu range 0.06-0.2% drives CLmt from approximately 1.5 L/h (FU = 0.06%) to 3.0 L/h (FU = 0.2%) at the typical individual (Table 2 footnote 'fu varied between 0.06 and 0.2%, hence mean CLmt varies from 1.5 to 3.0 l h-1').",
      source_name        = "f_u"
    ),
    SEXF = list(
      description        = "Sex indicator: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Hamren 2008 Table 1: 17 male / 6 female in the renal-insufficiency cohort and 10 male / 8 female in the control cohort (cohort-pooled 27 male / 14 female; sex_female_pct 34.1). Used as a multiplicative covariate on CLrt: 1 + (-0.31) * SEXF so that women have 31% lower renal clearance of tesaglitazar than men at the same iohexol-measured GFR (Table 2 'Gender on CLrt (% difference in female vs male) -31%').",
      source_name        = "SEX"
    ),
    CONMED_PROBENECID = list(
      description        = "Concomitant probenecid indicator: 1 = subject received probenecid co-administration around the modeled tesaglitazar dose, 0 = no probenecid.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no probenecid)",
      notes              = "Hamren 2008 Methods: 'For exploratory purposes, probenecid was administered to four subjects in the control group in order to investigate the likelihood of a drug-drug interaction with tesaglitazar. The first oral probenecid dose (500 mg) was given in the evening the day before start of tesaglitazar treatment, the second dose at start of tesaglitazar dosing and the third, on the same day, in the evening.' Used as a multiplicative covariate on the parent renal clearance CLrt and on the metabolite saturable renal Vmax: 1 + (-0.75) * CONMED_PROBENECID reduces both by 75% during probenecid intake (Table 2 'Probenecid on CLrt and Vmax (% difference during probenecid intake) -75%'). The published model used a single shared theta estimated jointly across both target parameters; this packaged encoding carries two separate effect parameters e_probenecid_cl_renal and e_probenecid_vmax_gluc initialised to the same -0.75 value to keep each covariate-parameter pair canonically named (see vignette Assumptions and deviations). Treated as time-fixed at the subject level here because the three probenecid doses around the first tesaglitazar dose were the only co-administration period; CLrt and Vmax effects are applied throughout for the four affected subjects to reproduce the published static-effect estimate.",
      source_name        = "PROBENECID"
    ),
    WT = list(
      description        = "Body weight (baseline; treated as time-fixed in the source paper).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hamren 2008 Table 1: median 84 kg (range 62-107) in the renal-insufficiency cohort and median 75 kg (range 60-97) in the control cohort. Reference 80 kg per the Table 2 footnote 'Centred on a subject with a body weight of 80 kg'. Linear centered effect 1 + 0.0094 * (WT - 80) shared on Vct and Vpt (Table 2 'Body weight on Vct and Vpt (% kg-1) 0.94'). Encoded as a single shared parameter e_wt_vc_vp per the canonical shared-exponent convention.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    n_observations = "707 plasma tesaglitazar + 323 plasma acyl glucuronide + 163 urine tesaglitazar + 163 urine acyl glucuronide concentration-time records (Results: Pharmacokinetic analysis).",
    age_range      = "34-78 years (renal-insufficiency cohort 34-78; healthy controls 41-73)",
    age_median     = "55 years (renal-insufficiency cohort), 53 years (healthy controls), Table 1",
    weight_range   = "60-107 kg (renal-insufficiency cohort 62-107 kg; healthy controls 60-97 kg)",
    weight_median  = "84 kg (renal-insufficiency cohort), 75 kg (healthy controls), Table 1",
    sex_female_pct = 34.1,
    race_ethnicity = "Not reported in source paper; recruitment at Swedish nephrology and clinical pharmacology sites (study SHSBC-0007, AstraZeneca R&D Molndal + Sahlgrenska University Hospital, Gothenburg).",
    disease_state  = "Adults with varying degrees of renal function. Renal-insufficiency cohort (n = 23, Part II only): mild GFR 51-80, moderate GFR 31-50, severe GFR 10-30 mL/min/1.73 m^2; no subjects on dialysis. Healthy controls (n = 18) matched for age and sex with renal-insufficiency subjects. Part I pilot study (n = 6 subjects with moderate or severe renal insufficiency) was pooled into the Part II analysis. Four healthy-control subjects also received oral probenecid (500 mg x 3 doses around the first tesaglitazar dose) as a probenecid-tesaglitazar drug-drug-interaction probe.",
    renal_function = "Assessed by plasma iohexol clearance (CLiohexol) on days 1 and 42; median 32 mL/min/1.73 m^2 (range 16-94) in the renal-insufficiency cohort, median 90 mL/min/1.73 m^2 (range 75-120) in the healthy-control cohort; pooled-cohort median 76 mL/min/1.73 m^2 used as the modelling reference.",
    dose_range     = "Tesaglitazar 0.5 mg orally once daily for 7 days (Part I pilot; 6 subjects with moderate or severe renal insufficiency) AND tesaglitazar 1 mg orally once daily for 42 +/- 3 days (Part II; 23 subjects with mild / moderate / severe renal insufficiency plus 18 healthy controls). 1 mg tesaglitazar approx 2.45 umol; 0.5 mg approx 1.22 umol (tesaglitazar molecular weight 408.45 g/mol).",
    regions        = "Sweden (open, stratified two-centre study SHSBC-0007).",
    sampling       = "Part I: predose and 1, 4, 24 h postdose on days 1 and 7. Part II: predose and 1, 2, 4, 12, 24 h postdose on days 1 and 42, trough samples on days 14 and 28, and four additional samples on days 3 +/- 1, 6 +/- 1, 9 +/- 1, 21 +/- 2 after the last dose intake. Urine: 0-6, 6-12, 12-24 h on day 1 and 0-24 h on day 42 of Part II.",
    notes          = "Patient characteristics from Table 1. One subject discontinued the study on day 28; data collected up to discontinuation were retained. Fraction-unbound data were missing for three individuals; the cohort median 0.1% was imputed for these. Probenecid was given to four control-cohort subjects as an exploratory drug-drug-interaction probe (Methods: 'Study design and population')."
  )

  ini({
    # All parameter estimates from Hamren 2008 Table 2 'Parameter estimates and
    # associated relative standard errors (RSE) from the final population PK
    # model of tesaglitazar and the acyl glucuronide'.
    #
    # Concentration / amount units throughout this file are molar (umol for
    # state amounts; umol/L for plasma concentrations) to match the
    # Michaelis-Menten parameterisation of the acyl-glucuronide renal
    # elimination. The user converts the mg-of-tesaglitazar dose to umol at
    # data assembly via dose_umol = dose_mg / 408.45 * 1000 (tesaglitazar MW
    # 408.45 g/mol; chemical structure depicted in Figure 1).

    # Parent tesaglitazar -- structural parameters (reference covariates:
    # CRCL = 76 mL/min/1.73 m^2, WT = 80 kg, no probenecid, male, FU = 0.1%)
    lka       <- fixed(log(1.5));  label("Absorption rate constant ka (1/h); fixed at 1.5 from prior PK knowledge of tesaglitazar")                # Table 2 ka = 1.5 (fixed); Discussion: 'The absorption rate constant of tesaglitazar was fixed in the present analysis due to limited plasma data early after dose administration.'
    lfdepot   <- fixed(log(1));    label("Bioavailability (fraction); fixed at 1 because oral absorption is complete")                              # Methods 'Pharmacokinetic modelling': 'Since the oral absorption of tesaglitazar is known to be complete, the absolute bioavailability (F) was set to 1.'
    lcl_renal <- log(0.027);       label("Renal clearance CLrt of tesaglitazar (L/h) at reference covariates")                                      # Table 2 CLrt = 0.027 (RSE 4.8%)
    lcl_met   <- log(1.91);        label("Metabolic clearance CLmt of tesaglitazar (L/h, parent -> acyl glucuronide formation) at reference covariates")  # Table 2 CLmt = 1.91 (RSE 8.2%)
    lq        <- log(0.22);        label("Inter-compartmental clearance Qt of tesaglitazar central <-> peripheral1 (L/h)")                          # Table 2 Qt = 0.22 (RSE 13%)
    lvc       <- log(2.1);         label("Central volume of distribution Vct of tesaglitazar (L) at reference WT 80 kg")                            # Table 2 Vct = 2.1 (RSE 8.1%)
    lvp       <- log(5.1);         label("Peripheral volume of distribution Vpt of tesaglitazar (L) at reference WT 80 kg")                         # Table 2 Vpt = 5.1 (RSE 3.2%)

    # Acyl glucuronide metabolite -- structural parameters (1-compartment
    # disposition; saturable Michaelis-Menten renal elimination + linear
    # non-renal clearance + first-order biliary excretion; all centered on
    # a male subject with CRCL = 76 mL/min/1.73 m^2 and not on probenecid).
    lvmax_gluc      <- log(0.188); label("Saturable renal elimination Vmax of acyl glucuronide (umol/h) at reference CRCL 76")                      # Table 2 Vmax = 0.188 (RSE 23%)
    lkm_gluc        <- log(0.041); label("Michaelis-Menten constant Km of acyl glucuronide renal elimination (umol/L)")                             # Table 2 Km = 0.041 (RSE 29%)
    lcl_nonren_gluc <- log(1.2);   label("Non-renal clearance CLnrm of acyl glucuronide (L/h, linear, loss from circulation)")                      # Table 2 CLnrm = 1.2 (RSE 11%)
    lkbm            <- log(11.7);  label("Biliary excretion rate constant kbm of acyl glucuronide (1/h, central_gluc -> gut_gluc)")                 # Table 2 kbm = 11.7 (RSE 5.9%)
    lkicv           <- log(0.79);  label("Interconversion rate constant kicv (1/h, gut hydrolysis + parent reabsorption back into central)")        # Table 2 'kint = 0.79' (RSE 7.3%); paper-named 'kint' renamed to canonical kicv to disambiguate from the TMDD-canonical kint (see references/parameter-names.md)
    lvc_gluc        <- log(8.5);   label("Central volume of distribution Vcm of acyl glucuronide (L)")                                              # Table 2 Vcm = 8.5 (RSE 17%)

    # Covariate effects (Table 2 'Significant covariate effects (multiplicative covariate models)')
    e_fu_cl_met            <- 5.55;    label("Linear centered slope of FU (%) on CLmt; CLmt = TVCLmt * (1 + e_fu_cl_met * (FU - 0.1))")             # Table 2 'fu on CLmt 555 %/unit fu' (RSE 23%); 555% per unit fu maps to coefficient 5.55 in the multiplicative form 1 + slope * (FU - 0.1)
    e_crcl_cl_renal        <- 0.0099;  label("Linear centered slope of CRCL on CLrt; CLrt = TVCLrt * (1 + e_crcl_cl_renal * (CRCL - 76))")          # Table 2 'GFR on CLrt 0.99 %/(ml min-1 1.73m-2)' (RSE 7.7%); 0.99% per mL/min/1.73 m^2 -> 0.0099 in the multiplicative form
    e_sexf_cl_renal        <- -0.31;   label("Sex effect on CLrt; CLrt multiplied by (1 + e_sexf_cl_renal * SEXF) (female -31%)")                   # Table 2 'Gender on CLrt (% difference in female vs male) -31%' (RSE 17%)
    e_probenecid_cl_renal  <- -0.75;   label("Probenecid effect on CLrt; CLrt multiplied by (1 + e_probenecid_cl_renal * CONMED_PROBENECID) (-75%)") # Table 2 'Probenecid on CLrt and Vmax (% difference during probenecid intake) -75%' (RSE 3.9%); shared estimate, encoded here and below as two separate parameters
    e_probenecid_vmax_gluc <- -0.75;   label("Probenecid effect on metabolite Vmax; Vmax multiplied by (1 + e_probenecid_vmax_gluc * CONMED_PROBENECID) (-75%)")  # Table 2 (same shared estimate as e_probenecid_cl_renal); see vignette Assumptions and deviations for the shared-vs-split parameterisation rationale
    e_wt_vc_vp             <- 0.0094;  label("Shared linear centered slope of WT (kg) on Vct and Vpt; multiplied by (1 + e_wt_vc_vp * (WT - 80))")  # Table 2 'Body weight on Vct and Vpt (% kg-1) 0.94' (RSE 33%); 0.94% per kg -> 0.0094 in the multiplicative form

    # Inter-individual variability (Table 2 parametric estimate column).
    # omega^2 = log(1 + (CV/100)^2) for log-normal IIV. The source paper
    # reports both a parametric and a nonparametric variance-covariance
    # matrix; the parametric diagonal CV% values are carried here. The
    # final published fit used the nonparametric matrix to improve
    # predictive performance (see vignette Assumptions and deviations).
    etalcl_met         ~ 0.0319   # Table 2 parametric CLmt   18% CV -> log(1 + 0.18^2) = 0.0319
    etalcl_renal       ~ 0.0392   # Table 2 parametric CLrt   20% CV -> log(1 + 0.20^2) = 0.0392
    etalvc             ~ 0.1416   # Table 2 parametric Vct    39% CV -> log(1 + 0.39^2) = 0.1416
    etalvp             ~ 0.0253   # Table 2 parametric Vpt    16% CV -> log(1 + 0.16^2) = 0.0253
    etalkm_gluc        ~ 0.1283   # Table 2 parametric Km     37% CV -> log(1 + 0.37^2) = 0.1283
    etalcl_nonren_gluc ~ 0.2561   # Table 2 parametric CLnrm  54% CV -> log(1 + 0.54^2) = 0.2561

    # Residual error (Table 2 'Residual error estimates (proportional models)').
    # The source paper used untransformed plasma and urine data with
    # proportional residual error for each analyte / matrix pair.
    propSd           <- 0.097; label("Proportional residual SD on plasma tesaglitazar Cc (fraction)")                       # Table 2 tesaglitazar in plasma 9.7% (RSE 7.9%)
    propSd_urineTesa <- 0.37;  label("Proportional residual SD on cumulative urine tesaglitazar urineTesa (fraction)")      # Table 2 tesaglitazar in urine   37% (RSE 6.1%)
    propSd_gluc      <- 0.18;  label("Proportional residual SD on plasma acyl glucuronide Cc_gluc (fraction)")              # Table 2 metabolite in plasma    18% (RSE 5.9%)
    propSd_urineGluc <- 0.26;  label("Proportional residual SD on cumulative urine acyl glucuronide urineGluc (fraction)")  # Table 2 metabolite in urine     26% (RSE 9.8%)
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters (log-normal IIV; covariate effects on CLrt,
    # CLmt, Vmax, Vct, Vpt).
    # ------------------------------------------------------------------
    ka     <- exp(lka)
    fdepot <- exp(lfdepot)

    cl_renal <- exp(lcl_renal + etalcl_renal) *
                (1 + e_crcl_cl_renal * (CRCL - 76)) *
                (1 + e_sexf_cl_renal * SEXF) *
                (1 + e_probenecid_cl_renal * CONMED_PROBENECID)
    cl_met   <- exp(lcl_met + etalcl_met) *
                (1 + e_fu_cl_met * (FU - 0.1))
    q        <- exp(lq)
    vc       <- exp(lvc + etalvc) *
                (1 + e_wt_vc_vp * (WT - 80))
    vp       <- exp(lvp + etalvp) *
                (1 + e_wt_vc_vp * (WT - 80))

    # Metabolite parameters (M-M renal Vmax scaled directly by CRCL / 76
    # per the paper text 'direct linear relationship between Vmax and
    # CLiohexol, normalized to the median observed renal function (76 ml
    # min-1 1.73 m-2)').
    vmax           <- exp(lvmax_gluc) *
                      (CRCL / 76) *
                      (1 + e_probenecid_vmax_gluc * CONMED_PROBENECID)
    km             <- exp(lkm_gluc + etalkm_gluc)
    cl_nonren_gluc <- exp(lcl_nonren_gluc + etalcl_nonren_gluc)
    kbm            <- exp(lkbm)
    kicv           <- exp(lkicv)
    vc_gluc        <- exp(lvc_gluc)

    # ------------------------------------------------------------------
    # Plasma concentrations (umol/L) used in M-M renal elimination of
    # the acyl glucuronide and as observation outputs.
    # ------------------------------------------------------------------
    Cc      <- central      / vc        # tesaglitazar plasma (umol/L)
    Cc_gluc <- central_gluc / vc_gluc   # acyl glucuronide plasma (umol/L)

    # ------------------------------------------------------------------
    # Micro-rate constants and fluxes (1/h or umol/h as noted).
    # ------------------------------------------------------------------
    kel_renal  <- cl_renal       / vc        # parent first-order CLrt rate (1/h)
    kel_met    <- cl_met         / vc        # parent first-order CLmt formation rate (1/h)
    k12        <- q              / vc        # parent central -> peripheral1 (1/h)
    k21        <- q              / vp        # parent peripheral1 -> central (1/h)
    kel_nonren <- cl_nonren_gluc / vc_gluc   # metabolite first-order non-renal CL rate (1/h)
    flux_mm    <- vmax * Cc_gluc / (km + Cc_gluc)  # metabolite M-M renal elimination (umol/h)

    # ------------------------------------------------------------------
    # ODE system. State variables (units in umol unless noted):
    #   depot          tesaglitazar in absorption depot
    #   central        tesaglitazar in central plasma
    #   peripheral1    tesaglitazar in peripheral 1
    #   central_gluc   acyl glucuronide in central plasma
    #   gut_gluc       acyl glucuronide accumulated in gut (paper-specific compartment for interconversion)
    #   urine          cumulative tesaglitazar excreted in urine
    #   urine_gluc     cumulative acyl glucuronide excreted in urine
    # Mass balance: the metabolite-formation flux (CLmt * Cc) leaves the
    # parent central compartment and arrives 1:1 in molar terms in
    # central_gluc; the interconversion flux (kicv * gut_gluc) leaves
    # gut_gluc as metabolite-equivalents and arrives 1:1 in molar terms
    # in central as parent tesaglitazar (the gut hydrolysis converts the
    # acyl glucuronide back to aglycone).
    # ------------------------------------------------------------------
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot -
                           kel_renal * central -
                           kel_met   * central -
                           k12       * central +
                           k21       * peripheral1 +
                           kicv      * gut_gluc
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(central_gluc) <-  kel_met    * central -
                           flux_mm -
                           kel_nonren * central_gluc -
                           kbm        * central_gluc
    d/dt(gut_gluc)     <-  kbm * central_gluc - kicv * gut_gluc
    d/dt(urine)        <-  kel_renal * central
    d/dt(urine_gluc)   <-  flux_mm

    # Oral bioavailability anchor (paper fixes F = 1).
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # Observations and residual errors. Cumulative urine outputs are
    # reported as masses (umol) directly from the cumulative compartments.
    # ------------------------------------------------------------------
    urineTesa <- urine
    urineGluc <- urine_gluc

    Cc        ~ prop(propSd)
    Cc_gluc   ~ prop(propSd_gluc)
    urineTesa ~ prop(propSd_urineTesa)
    urineGluc ~ prop(propSd_urineGluc)
  })
}

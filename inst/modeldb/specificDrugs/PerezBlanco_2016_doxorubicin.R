PerezBlanco_2016_doxorubicin <- function() {
  description <- paste(
    "Joint population PK model for IV doxorubicin (DOX, 3-compartment) and",
    "its active C-13 alcohol metabolite doxorubicinol (DOXol, 2-compartment)",
    "in adult patients with non-Hodgkin's lymphoma receiving R-CHOP",
    "chemotherapy (Perez-Blanco 2016). DOX was administered as a 0.5-h IV",
    "infusion at the protocol dose of 50 mg/m^2. The fraction Fm = 0.22 of",
    "total DOX clearance is routed to DOXol formation; the remaining",
    "(1 - Fm) fraction represents non-DOXol elimination pathways. The five",
    "volumes of distribution (V1/V2/V3 for DOX and V4/V5 for DOXol) were",
    "held fixed during estimation: the DOX volumes to the Kontny 2013",
    "(doi:10.1007/s00280-013-2261-3) adult-reference values, and the DOXol",
    "volumes to the values obtained from a sensitivity analysis carried out",
    "for that purpose. No covariates were retained in the final model;",
    "bilirubin and AST showed an influence on CL and CLm but the OFV",
    "decrease was not statistically significant. Residual variability is",
    "proportional for both DOX and DOXol. This was the first published",
    "two-compartment DOXol popPK model."
  )
  reference <- paste(
    "Perez-Blanco JS, Santos-Buelga D, Fernandez de Gatta MM,",
    "Hernandez-Rivas JM, Martin A, Garcia MJ.",
    "Population pharmacokinetics of doxorubicin and doxorubicinol in",
    "patients diagnosed with non-Hodgkin's lymphoma.",
    "Br J Clin Pharmacol. 2016;82(6):1517-1527.",
    "doi:10.1111/bcp.13070",
    sep = " "
  )
  vignette <- "PerezBlanco_2016_doxorubicin"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age in years (mean 66, range 26-84).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the stepwise covariate model on PK parameters but not retained -- did not show a significant influence on DOX CL (Perez-Blanco 2016 Discussion). Tabulated in Table 1."
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened on PK parameters but not retained. Cohort split 22F / 23M."
    ),
    WT = list(
      description        = "Body weight in kg (mean 71, range 43-110).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on PK parameters but not retained. Standard allometric scaling on CL/V was not preserved in the final model -- the paper attributes the lack of WT effect to the high homogeneity of the adult NHL population."
    ),
    HT = list(
      description        = "Height in m (mean 1.64, range 1.43-1.92).",
      units              = "m",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on PK parameters but not retained."
    ),
    BSA = list(
      description        = "Body surface area in m^2 (mean 1.8, range 1.3-2.3).",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on PK parameters but not retained. Used clinically for dose calculation (50 mg/m^2) but did not enter the structural model."
    ),
    BMI = list(
      description        = "Body mass index in kg/m^2 (mean 26.5, range 19.9-37.6).",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on PK parameters but not retained."
    ),
    LBW = list(
      description        = "Lean body weight in kg (mean 47.9, range 28.7-69.5).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed using the Janmahasatian formula; screened on PK parameters but not retained."
    ),
    CLCR = list(
      description        = "Creatinine clearance in mL/min (mean 91, range 40-201).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Estimated using the Cockcroft-Gault formula; screened on PK parameters but not retained."
    ),
    AST = list(
      description        = "Aspartate aminotransferase in IU/L (mean 25, range 12-64).",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and CLm; showed an influence on both, but the OFV decrease was not statistically significant and AST was not retained in the final model (Perez-Blanco 2016 Results, 'Final popPK model')."
    ),
    ALT = list(
      description        = "Alanine aminotransferase in IU/L (mean 23, range 7-88).",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on PK parameters but not retained."
    ),
    BILI = list(
      description        = "Total bilirubin in mg/dL (mean 0.44, range 0.10-0.70).",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and CLm; showed an influence on both, but the OFV decrease was not statistically significant and bilirubin was not retained in the final model."
    ),
    ECOG = list(
      description        = "Eastern Cooperative Oncology Group performance status.",
      units              = "(ordinal)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened on PK parameters but not retained."
    ),
    IPI = list(
      description        = "International Prognostic Index for NHL.",
      units              = "(ordinal)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened on PK parameters but not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 2L,
    n_observations = "125 DOX and 120 DOXol plasma concentrations; one outlier subject (|CWRES| > 4) was removed from the final n = 44 fit.",
    age_range      = "26-84 years (mean 66, SD 15; Table 1).",
    weight_range   = "43-110 kg (mean 71, SD 12; Table 1).",
    height_range   = "1.43-1.92 m (mean 1.64; Table 1).",
    bsa_range      = "1.3-2.3 m^2 (mean 1.8; Table 1).",
    sex_female_pct = 48.9,
    race_ethnicity = c(White = 100),
    disease_state  = "Adults with non-Hodgkin's lymphoma (diffuse large B-cell lymphoma n = 36, Burkitt-like lymphoma n = 2, follicular lymphoma n = 2, other n = 5) receiving R-CHOP every 21 days for six cycles. All patients had normal hepatic, renal, and cardiac function.",
    dose_range     = "Doxorubicin 50 mg/m^2 protocol; administered as a 0.5-h IV infusion (range 0.2-1.3 h; Table 1). Per-subject dose 53-130 mg (mean 89 mg).",
    regions        = "Six Spanish hospitals; 30 patients from the University Hospital of Salamanca, 15 patients from the GEL-R-COMP-2013 multi-centre trial (NCT02012088).",
    notes          = "Sparse PK sampling at 0, 30, 90 and 180 min after the end of the DOX infusion, selected by D-optimal design in WinNonlin 5.3. UHPLC-fluorescence assay with DOX LLOQ 8 ng/mL and DOXol LLOQ 3 ng/mL. Population PK fit in NONMEM 7.3 with FOCE, supported by PsN 4.4.0 and Xpose 4.5.0. Patients treated between June 2009 and June 2015."
  )

  ini({
    # ------------------------------------------------------------------
    # DOXORUBICIN STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Table 2 'Final model (n = 44)' column. The source uses an NM-TRAN
    # parameterisation in apparent clearances and volumes; the values
    # below are the back-transformed point estimates.
    #
    # Five of the eleven PK fixed-effects are FIXED:
    #   V1, V2, V3 -- carried from Kontny 2013 [doi:10.1007/s00280-013-2261-3]
    #     adult-reference DOX volumes (paper Results 'Final popPK model').
    #   V4, V5     -- fixed to the values obtained in a sensitivity analysis
    #     carried out on the DOXol volumes for this purpose (same paragraph).
    # The remaining six (CL, Q2, Q3, Fm, CLm, Q5) are estimated; RSE
    # values for the estimated parameters are tabulated in Table 2.

    lcl       <- log(62.4) ; label("Doxorubicin total clearance CL (L/h); includes the Fm-fraction routed to doxorubicinol formation and the (1 - Fm)-fraction non-DOXol elimination") # Table 2 final-model CL = 62.4 L/h, RSE 11.5%
    lvc       <- fixed(log(17.7))   ; label("Doxorubicin central volume V1 (L); FIXED from Kontny 2013") # Table 2 final-model V1 = 17.7 L FIX
    lq        <- log(50.7) ; label("Doxorubicin inter-compartmental clearance Q to peripheral1 (L/h); paper Q2") # Table 2 final-model Q2 = 50.7 L/h, RSE 18.4%
    lvp       <- fixed(log(1830))   ; label("Doxorubicin peripheral1 volume Vp (L); paper V2 FIXED from Kontny 2013") # Table 2 final-model V2 = 1830 L FIX
    lq2       <- log(28.4) ; label("Doxorubicin inter-compartmental clearance Q2 to peripheral2 (L/h); paper Q3") # Table 2 final-model Q3 = 28.4 L/h, RSE 13.5%
    lvp2      <- fixed(log(71))     ; label("Doxorubicin peripheral2 volume Vp2 (L); paper V3 FIXED from Kontny 2013") # Table 2 final-model V3 = 71 L FIX

    # ------------------------------------------------------------------
    # FRACTION METABOLISED TO DOXORUBICINOL
    # ------------------------------------------------------------------
    # Fm splits the total DOX clearance into the DOXol-formation
    # pathway (rate = Fm * CL/Vc * central) and the residual non-DOXol
    # elimination pathway (rate = (1 - Fm) * CL/Vc * central). The
    # paper's AUC equation confirms this parameterisation:
    #   AUC_DOX   = Dose / CL
    #   AUC_DOXol = Dose * Fm / CLm    (Methods 'Drug exposure and haematological toxicity')

    lfm <- log(0.22) ; label("Fraction of doxorubicin clearance routed to doxorubicinol formation (Fm, unitless)") # Table 2 final-model Fm = 0.22, RSE 14.7%

    # ------------------------------------------------------------------
    # DOXORUBICINOL STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # DOXol is described by a two-compartment model with first-order
    # distribution and elimination (paper Results 'Final popPK model').
    # The current paper is the first published two-compartment DOXol
    # popPK model; prior studies used a one-compartment DOXol.

    lcl_doxol <- log(26.8) ; label("Doxorubicinol total clearance CLm (L/h)") # Table 2 final-model CLm = 26.8 L/h, RSE 42.9%
    lvc_doxol <- fixed(log(79.8))   ; label("Doxorubicinol central volume V4 (L); FIXED from sensitivity analysis") # Table 2 final-model V4 = 79.8 L FIX
    lq_doxol  <- log(424)  ; label("Doxorubicinol inter-compartmental clearance Q5 to peripheral1_doxol (L/h)") # Table 2 final-model Q5 = 424 L/h, RSE 18.0%
    lvp_doxol <- fixed(log(653))    ; label("Doxorubicinol peripheral1 volume V5 (L); FIXED from sensitivity analysis") # Table 2 final-model V5 = 653 L FIX

    # ------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # Exponential IIV (paper Results 'Final popPK model'). Table 2
    # reports IIV as %CV via footnote 'The values of random-effect
    # parameters have been expressed as %'. NONMEM exponential IIV
    # relates CV% to the log-scale variance via
    #     omega^2 = log(1 + (CV/100)^2)
    # Three IIVs are FIXED at the model-building values to avoid
    # over-parameterisation (paper Results 'Final popPK model'):
    #   IIV(Q2), IIV(Q3), IIV(CLm).
    # No correlation between random-effect parameters was identified.

    etalcl       ~ 0.051112                                                                           # IIV(CL)  = 22.9% CV -> omega^2 = log(1 + 0.229^2) = 0.051112
    etalq        ~ fixed(0.344214)                                                                    # IIV(Q2)  = 64.1% CV -> omega^2 = log(1 + 0.641^2) = 0.344214 (FIXED at model-building value)
    etalq2       ~ fixed(0.076520)                                                                    # IIV(Q3)  = 28.2% CV -> omega^2 = log(1 + 0.282^2) = 0.076520 (FIXED at model-building value)
    etalfm       ~ 0.160322                                                                           # IIV(Fm)  = 41.7% CV -> omega^2 = log(1 + 0.417^2) = 0.160322
    etalcl_doxol ~ fixed(0.201130)                                                                    # IIV(CLm) = 47.2% CV -> omega^2 = log(1 + 0.472^2) = 0.201130 (FIXED at model-building value)
    etalq_doxol  ~ 0.297821                                                                           # IIV(Q5)  = 58.9% CV -> omega^2 = log(1 + 0.589^2) = 0.297821

    # ------------------------------------------------------------------
    # RESIDUAL ERROR (PROPORTIONAL)
    # ------------------------------------------------------------------
    # Proportional residual error for both DOX and DOXol on the linear
    # concentration scale (paper Results 'Final popPK model':
    # 'Residual variability was modelled using a proportional error
    # term for both entities'). Table 2 reports the residual SD as a
    # percentage (37.1% for DOX, 32.1% for DOXol) with RSE 8.3% and
    # 10.4% respectively.

    propSd       <- 0.371 ; label("Doxorubicin proportional residual SD on linear concentration (fraction)")   # Table 2 'IIVResidual (%)' = 37.1% (proportional residual SD; sigma in NONMEM Y = F*(1+EPS))
    propSd_doxol <- 0.321 ; label("Doxorubicinol proportional residual SD on linear concentration (fraction)") # Table 2 'IIVResidual,m (%)' = 32.1%
  })

  model({
    # Individual PK parameters.
    cl        <- exp(lcl       + etalcl)
    vc        <- exp(lvc)
    q         <- exp(lq        + etalq)
    vp        <- exp(lvp)
    q2        <- exp(lq2       + etalq2)
    vp2       <- exp(lvp2)
    fm        <- exp(lfm       + etalfm)
    cl_doxol  <- exp(lcl_doxol + etalcl_doxol)
    vc_doxol  <- exp(lvc_doxol)
    q_doxol   <- exp(lq_doxol  + etalq_doxol)
    vp_doxol  <- exp(lvp_doxol)

    # Doxorubicin disposition (3-compartment IV). Total DOX clearance
    # cl includes both the doxorubicinol-formation pathway (rate
    # = fm * cl / vc * central) and the residual non-DOXol pathway
    # (rate = (1 - fm) * cl / vc * central); the parent-loss term
    # therefore uses cl / vc directly, and the metabolite-input term
    # below uses fm * cl / vc separately.
    d/dt(central)     <- -cl / vc * central -
                          q  / vc * central + q  / vp  * peripheral1 -
                          q2 / vc * central + q2 / vp2 * peripheral2
    d/dt(peripheral1) <-  q  / vc * central - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc * central - q2 / vp2 * peripheral2

    # Doxorubicinol disposition (2-compartment). Metabolite input
    # is fm * cl / vc * central; the rest of the DOX clearance does
    # not appear here because it represents elimination to other
    # (non-DOXol) products.
    d/dt(central_doxol)     <-  fm * cl / vc * central -
                                cl_doxol / vc_doxol * central_doxol -
                                q_doxol  / vc_doxol * central_doxol +
                                q_doxol  / vp_doxol * peripheral1_doxol
    d/dt(peripheral1_doxol) <-  q_doxol  / vc_doxol * central_doxol -
                                q_doxol  / vp_doxol * peripheral1_doxol

    # Plasma concentrations in mg/L (dose in mg, volumes in L).
    # 1 mg/L = 1 ug/mL = 1000 ng/mL; the paper's bioanalytical LLOQs
    # of 8 ng/mL (DOX) and 3 ng/mL (DOXol) correspond to 0.008 mg/L
    # and 0.003 mg/L respectively.
    Cc       <- central       / vc
    Cc_doxol <- central_doxol / vc_doxol

    Cc       ~ prop(propSd)
    Cc_doxol ~ prop(propSd_doxol)
  })
}

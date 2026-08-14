Baiardi_2025_dalbavancin <- function() {
  description <- "Two-compartment population PK model with first-order elimination for intravenous dalbavancin in adult patients receiving multidose (off-label, long-term) regimens for difficult-to-treat Gram-positive infections with creatinine clearance above 30 mL/min. Total body weight is scaled by allometry with exponents fixed a priori at 0.75 on the clearances (CL, Q) and 1 on the volumes (V1, V2); no other covariate improved the fit, including creatinine clearance. Unbound concentration is returned as Ccu using the reported 7 percent free fraction, because the paper's PK/PD target is 100 percent fT above 4x MIC (total dalbavancin above 14.29 mg/L at the 0.25 mg/L EUCAST/USCAST breakpoint)."
  reference <- paste(
    "Baiardi G, Cameran Caviglia M, Boni S, Di Paolo A, Marini V, Cangemi G,",
    "Cafaro A, Pontali E, Mattioli F. Multidose Dalbavancin Population",
    "Pharmacokinetic Analysis for Prolonged Target Attainment in Patients",
    "Requiring Long-Term Treatment. Antibiotics. 2025;14(2):190.",
    "doi:10.3390/antibiotics14020190.",
    sep = " "
  )
  vignette <- "Baiardi_2025_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Baiardi 2025 Results",
        "section 2.2: 'Body weight scaled by allometric principles was set a",
        "priori as a covariate with fixed exponents of 0.75 on clearances (CL,",
        "Q) and of 1 on volumes of distribution (V1, V2).' The paper does not",
        "print the reference (standard) weight; 70 kg is used here because the",
        "three allometry references it cites for this step (Holford 1996 'A Size",
        "Standard for Pharmacokinetics'; Anderson & Holford 2008; Holford &",
        "Anderson 2017) all define the 70 kg size standard. The cohort median of",
        "72 kg (Table 1) is within 3% of 70 kg and is numerically",
        "indistinguishable at the resolution of any published result; see the",
        "validation vignette Assumptions and deviations section. Median 72 kg,",
        "range 44-179 kg (Table 1); the Monte Carlo simulations spanned 40-200",
        "kg in three weight bands (40-80, 80-120, 120-200 kg)."
      ),
      source_name        = "Weight"
    )
  )

  # Covariates screened by Baiardi 2025 but NOT retained in the final model.
  # Documentation only -- none of these is referenced in model().
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened on CL by both linear and non-linear relationships and NOT",
        "retained: Baiardi 2025 Results 2.2 reports the addition of CLcr on CL",
        "changed the objective function value by only 0.02 units, far below the",
        "3.84-unit forward-inclusion threshold. The paper attributes this to",
        "dalbavancin's mixed renal/non-renal elimination and to the cohort",
        "having CLcr > 30 mL/min throughout. Median 66.2 mL/min, range 31.7-283",
        "mL/min (Table 1). The packaged model therefore applies only to patients",
        "with CLcr > 30 mL/min; dose reduction is required below that threshold."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (centered on the cohort median) but not retained; Baiardi 2025 Methods 4.2 and Results 2.2. Median 72 years, range 26-97 years (Table 1)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened (centered on the cohort median) but not retained; Baiardi 2025 Methods 4.2 and Results 2.2. Median 168 cm, range 140-195 cm (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened (centered on the cohort median) but not retained; Baiardi 2025 Methods 4.2 and Results 2.2. Median 33.1 g/L, range 21.4-39.6 g/L (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened (centered on the cohort median) but not retained; Baiardi 2025 Methods 4.2 and Results 2.2. Median 1 mg/dL, range 0.5-1.7 mg/dL (Table 1)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened (centered on the cohort median; printed as 'PCR' in Baiardi 2025 Methods 4.2) but not retained. Median 8.4 mg/L, range 0.3-183.7 mg/L (Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    n_observations = 195L,
    age_range      = "26-97 years",
    age_median     = "72 years",
    weight_range   = "44-179 kg",
    weight_median  = "72 kg",
    sex_female_pct = 30,
    race_ethnicity = "Not reported; single-centre Italian cohort.",
    disease_state  = paste(
      "Adults with documented or suspected Gram-positive infection who failed",
      "primary antimicrobial therapy and were prescribed multidose dalbavancin",
      "as rescue treatment. Sites of infection (Table 1): ABSSSI (10),",
      "spondylodiscitis (9), septic arthritis (5), endocarditis (2), chronic",
      "knee bursitis (1), osteomyelitis (1), spinal implant infection (1),",
      "prostatic abscess (1). Isolates included MSSA (4), MRSE (5), MRSA (5),",
      "S. agalactiae, S. epidermidis, S. sanguinis and S. lugdunensis (1 each);",
      "12 had no growth or unavailable culture results."
    ),
    dose_range     = "At least two 1500 mg intravenous doses (30 min infusion) 7-14 days apart, at the discretion of the prescribing physician; median 3 doses per patient, range 2-10 (Table 1).",
    regions        = "Italy (Ente Ospedaliero Ospedali Galliera, Genoa)",
    renal_function = "Creatinine clearance (Cockcroft-Gault) always above 30 mL/min: median 66.2 mL/min, range 31.7-283 mL/min (Table 1). The model is not applicable to severe renal impairment (CLcr < 30 mL/min), for which the label requires dose reduction.",
    notes          = paste(
      "Prospective/retrospective single-centre TDM study (DALT DRUM, Liguria",
      "Territorial Ethics Committee 82/2023) run February 2023 to February 2024;",
      "baseline demographics in Baiardi 2025 Table 1. Sampling was Ctrough plus",
      "end-of-infusion Cmax at each administration, with additional sparse",
      "samples where feasible (median 5.5 TDM instances per patient, range 2-20).",
      "195 concentrations entered the analysis; five points with CWRES > 3 were",
      "excluded during model building. Total plasma dalbavancin was measured by",
      "LC-MS/MS over 0.66-800 mg/L (Cafaro 2024). Estimation was FOCEI in NONMEM",
      "7.5. The published Monte Carlo target-attainment analysis simulated 100,000",
      "profiles from a 3000-subject virtual population (1000 per weight band,",
      "40-200 kg) with covariates harvested from NHANES."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Baiardi 2025 Table 2 ("Final model
    # parameters estimate"), reported for a 70 kg reference subject (see
    # the WT entry in covariateData for the reference-weight rationale).
    # Table 2 prints CL, V1, Q and V2 on the linear scale; they are
    # log-transformed here per the nlmixr2lib convention.
    # ------------------------------------------------------------------
    lcl <- log(0.0273); label("Clearance CL (L/h)")                        # Baiardi 2025 Table 2: CL = 0.0273 L/h (RSE 5.1%; bootstrap median 0.0272, CI 0.0251-0.0294)
    lvc <- log(3.6);    label("Central volume of distribution V1 (L)")     # Baiardi 2025 Table 2: V1 = 3.6 L (RSE 3.7%; bootstrap median 3.6, CI 3.4-3.8)
    lq  <- log(0.0225); label("Intercompartmental clearance Q (L/h)")      # Baiardi 2025 Table 2: Q = 0.0225 L/h (RSE 28.4%; bootstrap median 0.0223, CI 0.0175-0.02936)
    lvp <- log(6.4);    label("Peripheral volume of distribution V2 (L)")  # Baiardi 2025 Table 2: V2 = 6.4 L (RSE 11.9%; bootstrap median 6.5, CI 5.6-7.5)

    # ------------------------------------------------------------------
    # Allometric exponents. Both were set a priori and held fixed, not
    # estimated: Baiardi 2025 Results 2.2 "Body weight scaled by
    # allometric principles [21-23] was set a priori as a covariate with
    # fixed exponents of 0.75 on clearances (CL, Q) and of 1 on volumes
    # of distribution (V1, V2)." Neither appears in Table 2, consistent
    # with being fixed rather than estimated.
    # ------------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of body weight on CL and Q (unitless)")   # Baiardi 2025 Results 2.2 (fixed a priori)
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent of body weight on V1 and V2 (unitless)")  # Baiardi 2025 Results 2.2 (fixed a priori)

    # ------------------------------------------------------------------
    # Plasma unbound fraction. Not estimated in the popPK fit; used by
    # the paper to convert total to free dalbavancin for the
    # 100% fT > 4x MIC target. Baiardi 2025 Discussion prints the
    # conversion as a display equation (p. 9): "If f DAL = 0.07 x
    # TotalDAL and MIC = 0.25 mg/L; TotalDAL > 4 x 0.25 mg/L / 0.07;
    # TotalDAL > 14.29 mg/L". Note that this display equation is dropped
    # by the automated markdown trim of the PDF -- 0.07 appears verbatim
    # only in the typeset equation, not in the surrounding prose, which
    # states it as "a 7% free fraction of DAL in plasma".
    # ------------------------------------------------------------------
    fu <- fixed(0.07); label("Dalbavancin plasma unbound fraction (unitless)")  # Baiardi 2025 Discussion display equation p. 9: f DAL = 0.07 x TotalDAL; consistent with the 93% protein binding cited in the Introduction

    # ------------------------------------------------------------------
    # Inter-individual variability. Baiardi 2025 Table 2 reports the
    # random effects as CV% under the heading "Random Effects (CV%)".
    # The equivalent log-normal variance is omega^2 = log(CV^2 + 1).
    # Eta shrinkage from Table 2 is quoted alongside each entry.
    # ------------------------------------------------------------------
    etalcl ~ 0.047265  # log(0.220^2 + 1); Baiardi 2025 Table 2 omega CL = 22%   (RSE 13.6%; bootstrap 21.4%, CI 16.0-26.3; eta shrinkage 8%)
    etalvc ~ 0.029490  # log(0.173^2 + 1); Baiardi 2025 Table 2 omega V1 = 17.3% (RSE 17.4%; bootstrap 16.9%, CI 11.0-21.2; eta shrinkage 17%)
    etalq  ~ 0.271919  # log(0.559^2 + 1); Baiardi 2025 Table 2 omega Q  = 55.9% (RSE 21.2%; bootstrap 53.0%, CI 31.3-71.1; eta shrinkage 30%)
    etalvp ~ 0.086729  # log(0.301^2 + 1); Baiardi 2025 Table 2 omega V2 = 30.1% (RSE 15.9%; bootstrap 28.6%, CI 19.4-37.4; eta shrinkage 34%)

    # ------------------------------------------------------------------
    # Residual variability. Table 2 reports a single proportional term
    # "b" with no additive component.
    # ------------------------------------------------------------------
    propSd <- 0.144; label("Proportional residual error (fraction)")  # Baiardi 2025 Table 2: b (residual error, proportional) = 0.144 (RSE 10%; bootstrap 0.144, CI 0.126-0.166)
  })

  model({
    # 1. Individual PK parameters. Allometric scaling on total body
    #    weight about the 70 kg size standard, exponents fixed a priori.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system -- two compartments with first-order elimination from
    #    the central compartment. Dalbavancin is given intravenously, so
    #    there is no absorption compartment; doses enter `central`
    #    directly (as a 30 min infusion in the source study).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 4. Observations. Cc is total plasma dalbavancin (the quantity the
    #    LC-MS/MS assay and the model were fit to) and Ccu is the free
    #    concentration that drives the paper's 100% fT > 4x MIC target.
    Cc  <- central / vc
    Ccu <- Cc * fu

    Cc ~ prop(propSd)
  })
}

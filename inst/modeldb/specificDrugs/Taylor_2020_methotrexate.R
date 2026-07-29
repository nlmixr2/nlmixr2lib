Taylor_2020_methotrexate <- function() {
  description <- "Three-compartment population PK model for intravenous high-dose methotrexate (5 or 8 g/m^2 over 24 h IV infusion) in pediatric NOPHO ALL2000 / ALL2008 patients with acute lymphoblastic leukemia; BSA-normalized PK parameters (reference 1.73 m^2) and a time-varying serum creatinine power effect on clearance (reference 29 umol/L) implemented as the default population PK model behind the MTXPK.org clinical decision support tool (Taylor 2020)"
  reference <- paste(
    "Taylor ZL, Mizuno T, Punt NC, Baskaran B, Navarro Sainz A, Shuman W,",
    "Felicelli N, Vinks AA, Heldrup J, Ramsey LB. (2020). MTXPK.org:",
    "A Clinical Decision Support Tool Evaluating High-Dose Methotrexate",
    "Pharmacokinetics to Inform Post-Infusion Care and Use of Glucarpidase.",
    "Clin Pharmacol Ther 108(3):635-643.",
    "doi:10.1002/cpt.1957.",
    sep = " "
  )
  vignette <- "Taylor_2020_methotrexate"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear normalization to 1.73 m^2 applied to all PK parameters (CL, V1, Q2, V2, Q3, V3) per Taylor 2020 Table 2 step 2 and Table 3 column headers (parameters reported as L/h/1.73 m^2 and L/1.73 m^2). The paper does not state the BSA computation formula (DuBois / Mosteller / Haycock); record 'unspecified' when assembling event data.",
      source_name        = "BSA"
    ),
    CREAT = list(
      description        = "Serum creatinine (time-varying)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying through the course of MTX infusion; recorded prior to each MTX cycle and at least daily during cycles. Power effect on CL only with reference 29 umol/L (the NOPHO population median; equivalent to 0.33 mg/dL) and exponent -0.247 per Taylor 2020 Table 3 and the covariate-model equation on page 638.",
      source_name        = "SCr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 772L,
    n_studies      = 2L,
    n_courses      = 4986L,
    n_observations = 31672L,
    age_range      = "1-18.83 years (median 4)",
    age_median     = "4 years",
    weight_range   = "7.2-105 kg (median 17.8)",
    weight_median  = "17.8 kg",
    bsa_range      = "0.4-2.31 m^2 (median 0.745)",
    bsa_median     = "0.745 m^2",
    sex_female_pct = 43,
    race_ethnicity = "Primarily European (Nordic / Baltic); race not reported in the source dataset and not used as a covariate.",
    disease_state  = "Philadelphia-chromosome-negative acute lymphoblastic leukemia treated under the NOPHO ALL2000 and ALL2008 protocols.",
    dose_range     = "5 or 8 g/m^2 IV over a 24 h infusion (Taylor 2020 Methods, range 0.6-10.1 g/m^2); 6-8 courses per patient with folinic acid rescue starting 36 or 42 h after the start of infusion.",
    regions        = "Denmark (27%), Finland (10%), Norway (8%), Sweden (55%).",
    notes          = "Patient demographics from Taylor 2020 Table 1. The NOPHO database excluded patients who received glucarpidase because the MTX immunoassay has interference from DAMPA (the glucarpidase cleavage product). 48 patients and 679 concentrations were excluded for missing dosing information; an additional 753 concentrations were excluded by the investigators as implausible because they caused NONMEM minimization errors. 5,535 of the 31,672 concentrations were recorded >= 96 h after the start of infusion (delayed-clearance follow-up). Finnish patients had ~26% faster clearance than Swedish / Danish / Norwegian patients in a post-hoc one-way ANOVA, attributed to earlier pre-hydration timing; country was not retained as a covariate in the final model."
  )

  ini({
    # Structural PK parameters (BSA-normalized typical values for a reference
    # subject with BSA = 1.73 m^2 and SCr = 29 umol/L). Values from Taylor 2020
    # Table 3 "Final model parameter estimates" / column "Mean". The model is
    # parameterized on the absolute scale per 1.73 m^2 (L/h or L); individual
    # parameters scale linearly with BSA inside model().
    lcl  <- log(11)    ; label("Clearance (L/h, normalized to 1.73 m^2)")            # Taylor 2020 Table 3: CL 11 L/h/1.73 m^2 (RSE 0.7%)
    lvc  <- log(16.5)  ; label("Central volume (V1, L, normalized to 1.73 m^2)")     # Taylor 2020 Table 3: V1 16.5 L/1.73 m^2 (RSE 5.2%)
    lq   <- log(0.602) ; label("Intercompartmental clearance to vascular peripheral (Q2, L/h, normalized to 1.73 m^2)")   # Taylor 2020 Table 3: Q2 0.602 L/h/1.73 m^2 (RSE 4.1%)
    lvp  <- log(4.55)  ; label("Vascular peripheral volume (V2, L, normalized to 1.73 m^2)")                                # Taylor 2020 Table 3: V2 4.55 L/1.73 m^2 (RSE 3.4%)
    lq2  <- log(0.111) ; label("Intercompartmental clearance to non-vascular peripheral (Q3, L/h, normalized to 1.73 m^2)") # Taylor 2020 Table 3: Q3 0.111 L/h/1.73 m^2 (RSE 2%)
    lvp2 <- log(13.1)  ; label("Non-vascular peripheral volume (V3, L, normalized to 1.73 m^2)")                            # Taylor 2020 Table 3: V3 13.1 L/1.73 m^2 (RSE 5%)

    # Covariate effect: power exponent of SCr / 29 on CL only.
    e_creat_cl <- -0.247; label("Power exponent of serum creatinine (SCr / 29) on CL (unitless)")   # Taylor 2020 Table 3: SCr -0.247 (RSE 5.7%); covariate equation on page 638

    # Inter-individual variability (log-normal). Taylor 2020 Table 3 reports
    # the "Interindividual variability" column as the variance (omega^2) on the
    # log-transformed parameter (consistent with NONMEM $OMEGA on a LOG-THETA
    # parameterization), so the value is used directly. V1 and Q2 IIVs
    # approached zero and were fixed at 0 in the final model (Taylor 2020
    # Results, page 638; Table 3 entries "0 FIX").
    etalcl  ~ 0.08          # Taylor 2020 Table 3: IIV(CL) variance 0.08 (RSE 4.7%)
    etalvc  ~ fixed(0)      # Taylor 2020 Table 3: IIV(V1) 0 FIX (variance fixed at zero in the final model)
    etalq   ~ fixed(0)      # Taylor 2020 Table 3: IIV(Q2) 0 FIX (variance fixed at zero in the final model)
    etalvp  ~ 0.12          # Taylor 2020 Table 3: IIV(V2) variance 0.12 (RSE 4.3%)
    etalq2  ~ 0.13          # Taylor 2020 Table 3: IIV(Q3) variance 0.13 (RSE 6.7%)
    etalvp2 ~ 0.10          # Taylor 2020 Table 3: IIV(V3) variance 0.10 (RSE 14.4%)

    # Residual error. Taylor 2020 describes an "exponential residual error
    # model" with log-transformation of both sides (Methods page 637 and
    # Supplementary Appendix S1) but does not report the numeric variance
    # estimate (Table 3 has no sigma row). A small placeholder is supplied
    # so the observation Cc has a valid residual model in nlmixr2; users
    # who need a realistic stochastic prediction interval should substitute
    # a paper-class typical value (e.g. 20-30% CV) and document the choice.
    # Precedent: Jelliffe_2014_digoxin.R applies the same placeholder pattern.
    propSd <- fixed(0.01); label("Proportional residual error (placeholder; not reported in source)")   # not reported in Taylor 2020 -- see vignette Assumptions and deviations
  })

  model({
    # Individual PK parameters. All six structural parameters scale linearly
    # with BSA (normalized to 1.73 m^2). Only CL carries the serum-creatinine
    # power effect; SCr enters as a time-varying covariate. V1 and Q2 have no
    # IIV (fixed at 0 above) so their etas contribute zero by construction.
    cl  <- exp(lcl  + etalcl)  * (BSA / 1.73) * (CREAT / 29)^e_creat_cl
    vc  <- exp(lvc  + etalvc)  * (BSA / 1.73)
    q   <- exp(lq   + etalq)   * (BSA / 1.73)
    vp  <- exp(lvp  + etalvp)  * (BSA / 1.73)
    q2  <- exp(lq2  + etalq2)  * (BSA / 1.73)
    vp2 <- exp(lvp2 + etalvp2) * (BSA / 1.73)

    # Three-compartment IV PK with first-order elimination from the central
    # compartment (NONMEM ADVAN11-style mass balance). Dose is delivered to
    # `central` over the 24 h infusion via the event-table rate/dur column.
    d/dt(central)     <-  q  / vp  * peripheral1 + q2 / vp2 * peripheral2 -
                          (cl + q + q2) / vc * central
    d/dt(peripheral1) <-  q  / vc  * central     -  q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc  * central     -  q2 / vp2 * peripheral2

    # Plasma concentration: dose in mg / volume in L gives mg/L (= ug/mL).
    # Taylor 2020 reports MTX concentrations in umol/L; users comparing to
    # the paper's figures or the MTXPK.org outputs should multiply Cc by
    # 1000 / 454.44 (methotrexate MW, g/mol) to obtain umol/L.
    Cc <- central / vc

    # Proportional residual error on the linear concentration scale; in the
    # source paper the exponential residual model maps to a log-additive form
    # under log-transformation of both sides.
    Cc ~ prop(propSd)
  })
}

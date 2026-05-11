Attarwala_2023_mRNA3927 <- function() {
  description <- "Preclinical (mouse, rat, cynomolgus monkey; allometrically scalable to humans). Translational semi-mechanistic PK and PK/PD model for mRNA-3927, an LNP-encapsulated dual mRNA encoding propionyl-CoA carboxylase (PCC) subunits PCCA and PCCB. PK: 3-compartment plasma1-tissue-plasma2 redistribution (V shared between the two plasma compartments; V and V2 fixed at the mouse reference and scaled allometrically) with body-weight allometric scaling of clearances (mouse reference 0.025 kg; estimated exponents cla on CL12/CL32 and clb on CL23/CL20). PD: liver PCC protein 2-compartment indirect-response model driven by an effect compartment linked to plasma mRNA, with synthesis linear in effect-compartment mRNA concentration and first-order degradation. Three downstream biomarkers (2-methylcitrate, 3-hydroxypropionate, C3/C2 carnitine ratio) follow direct sigmoidal Imax suppression by liver PCC protein with Imax fixed at 0.999."
  reference <- "Attarwala H, Lumley M, Liang M, Ivaturi V, Senn J. Translational Pharmacokinetic/Pharmacodynamic Model for mRNA-3927, an Investigational Therapeutic for the Treatment of Propionic Acidemia. Nucleic Acid Ther. 2023;33(2):141-147. doi:10.1089/nat.2022.0036"
  vignette <- "Attarwala_2023_mRNA3927"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (used for allometric scaling across species).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight is 0.025 kg (PCC-deficient mouse strain 4011020 (A138T) PCCA -/-, the species the PK + PD model was anchored to). Scales all PK clearances and volumes: volume exponent fixed at 1, clearance exponents estimated (cla = 0.631 on CL12 and CL32; clb = 1.10 on CL23 and CL20). The paper uses this same allometric scheme to extrapolate from mouse / rat / cyno PK directly to humans (Discussion, p. 145). PCC protein and downstream biomarker PD parameters are NOT body-weight-scaled in the published model; their fits are mouse-only and the human extrapolation re-uses them unchanged.",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "preclinical: mouse (PCC-deficient strain 4011020 (A138T) PCCA -/-; N = 111) + rat (juvenile Sprague-Dawley; N = 19) + cynomolgus monkey (N = 16). Model is allometrically scalable to humans (the paper's intended use; Discussion p. 145) by setting WT to the human value and using the same fixed structural parameters and allometric exponents.",
    n_subjects     = 146L,
    n_studies      = 3L,
    n_mice         = 111L,
    n_rats         = 19L,
    n_monkeys      = 16L,
    age_range      = "preclinical: juvenile rats and monkeys; mice not specified.",
    weight_range   = "mouse ~25 g (reference), rat ~200 g (juvenile Sprague-Dawley), cynomolgus monkey ~3 kg (paper does not report exact body weights; the allometric scaling uses BW / 0.025 with body weight supplied per subject).",
    sex_female_pct = NA_real_,
    disease_state  = "Propionic acidemia (PCC deficiency) in the mouse cohort; healthy juvenile Sprague-Dawley rats and cynomolgus monkeys.",
    dose_range     = "mouse 0.2, 0.5, 1, 2 mg/kg IV bolus (single dose) and 0.5, 2 mg/kg IV bolus q3W x 4 (multiple dose); rat 1, 3, 9 mg/kg IV bolus q2W (3 doses); monkey 1, 3, 5 mg/kg IV bolus q2W (3 doses).",
    regions        = "Preclinical study at Moderna, Inc. (Cambridge, MA, USA); rat / monkey work conducted under Moderna IACUC.",
    notes          = "PK was fit simultaneously across the three species using allometric scaling on body weight; PD (PCC protein and 2-MC / 3-HP / C3/C2 biomarkers) was fit on PCC-deficient mice only. The published model is the back-end of the first-in-human dose-selection analysis for Phase 1/2 study NCT04159103 in patients (>= 1 year of age) with propionic acidemia. Estimation method: FOCEI (Pumas 2.0) except 3-HP which used Laplacian. Per-subject body weight is required to recover species-specific parameter values; supply WT = 0.025 (mouse), 0.20 (rat), 3.0 (monkey), or 70 (adult human) for typical-individual simulations."
  )

  ini({
    # =========================================================================
    # PK model for plasma PCCA/B mRNA (Table 1, Supp Table S1 part A).
    # 3-compartment redistribution: A1 = plasma 1 (input here), A2 = tissue,
    # A3 = plasma 2; both plasma compartments share volume V; tissue has
    # its own volume V2 and is the elimination compartment.
    # Observed plasma concentration is C13 = C1 + C3.
    # All values typical-value estimates for a 0.025 kg PCC-deficient mouse
    # before WT scaling.
    # =========================================================================
    lcl12 <- log(19.7);            label("CL12, plasma1 -> tissue clearance (mL/h, mouse reference)")   # Table 1, RSE 23.1%
    lcl23 <- log(0.215);           label("CL23, tissue -> plasma2 clearance (mL/h, mouse reference)")   # Table 1, RSE 35.9%
    lcl32 <- log(2.96);            label("CL32, plasma2 -> tissue clearance (mL/h, mouse reference)")   # Table 1, RSE 48.2%
    lcl20 <- log(0.136);           label("CL20, tissue elimination clearance (mL/h, mouse reference)")  # Table 1, RSE 3.78%
    lvc   <- fixed(log(2.67));     label("V, distribution volume of plasma1 (= plasma2) (mL, mouse)")   # Table 1 footnote f, fixed; both plasma cmts share V
    lvp   <- fixed(log(0.961));    label("V2, distribution volume of tissue (mL, mouse)")               # Table 1 footnote f, fixed

    # Allometric exponents on body weight (paper symbols clalpha, clbeta).
    # cla (clalpha) applies to plasma -> tissue clearances (CL12 and CL32);
    # clb (clbeta)  applies to tissue -> plasma + tissue-elimination (CL23 and CL20).
    # Volume exponents are fixed at 1 per Supp Table S1.
    e_wt_cla <- 0.631;             label("Allometric exponent on CL12 and CL32 (clalpha, unitless)")    # Table 1, RSE 9.2%
    e_wt_clb <- 1.10;              label("Allometric exponent on CL23 and CL20 (clbeta, unitless)")     # Table 1, RSE 7.55%

    # IIV: only CL32 has IIV in the published PK model (Table 1).
    # omega^2 = log(1 + CV^2) ; CV% = 52.7 -> log(1 + 0.527^2) = 0.24510.
    etalcl32 ~ 0.24510   # Table 1, IIV on CL32 = 52.7 percent CV

    # Proportional residual error on plasma mRNA (Table 1, "Proportional residual error (%)" 37.5, RSE 11.0%).
    propSd <- 0.375;             label("Proportional residual SD on plasma mRNA (fraction)")          # Table 1

    # =========================================================================
    # PCC protein PK/PD model (Table 1 second block, Supp Table S1 part B).
    # Two-compartment indirect response with an effect compartment driven by
    # plasma mRNA C13. PCC central protein receives a synthesis input that is
    # LINEAR in the effect-compartment mRNA concentration (rate = Slope * Ce),
    # exchanges with a peripheral PCC compartment via kq, and decays first-order
    # with kdeg. See vignette Assumptions section: the supplement's d(PCC)/dt
    # equation is written "Ksyn * Ce" which is internally inconsistent with the
    # paper's main text ("linear function of plasma PCCA/B mRNA concentration",
    # p.143) and with the Slope units (mg/g/h)/(mg/mL); we implement the
    # linear-in-Ce form intended by the main text.
    # =========================================================================
    lke0    <- log(0.242);   label("ke0, equilibrium rate constant for effect compartment (1/h)")  # Table 1, RSE 3.35%
    lslope  <- log(55.6);    label("Slope, PCC synthesis per unit effect-compartment mRNA ((mg/g/h)/(mg/mL))") # Table 1, RSE 1.93%
    lkdeg   <- log(0.0075);  label("kdeg, PCC protein degradation rate (1/h)")                     # Table 1, RSE 2.82%
    lkq     <- log(0.00474); label("kq, intercompartmental rate constant for PCC protein (1/h)")   # Table 1, RSE 6.7%

    # Proportional residual error on liver PCC protein (Table 1, 30.5%, RSE 0.02%).
    propSd_PCC <- 0.305;            label("Proportional residual SD on liver PCC protein (fraction)") # Table 1

    # =========================================================================
    # Biomarker PD (Table 1 third block, Supp Table S1 part C).
    # Each biomarker (2-MC, 3-HP, C3/C2) is modeled as:
    #   biomarker = E0 + base * (1 - Imax * PCC / (IC50 + PCC))
    # where Imax is fixed at 0.999 (paper footnote m), E0 is the level not
    # amenable to PCC-mediated suppression, and base is the level that IS
    # amenable to suppression. PCC here is the central PCC protein
    # concentration in liver (mass-per-tissue units).
    # =========================================================================
    imax <- fixed(0.999);          label("Maximum inhibition fraction (fixed)")                    # Table 1 footnote m

    # ---- 2-methylcitrate (2-MC) ----
    le0_mc2    <- log(1.67);   label("E0 for 2-MC (umol/L)")                                       # Table 1, RSE 4.46%
    lbase_mc2  <- log(2.93);   label("Suppressible baseline for 2-MC (umol/L)")                    # Table 1, RSE 6.66%
    lic50_mc2  <- log(21.0);   label("IC50 of PCC for 2-MC suppression (PCC mass per liver tissue)") # Table 1, RSE 8.96%
    propSd_MC2  <- 0.234;       label("Proportional residual SD on 2-MC (fraction)")                # Table 1, RSE 3.84%
    # IIV (Table 1, "IIV (%)" column for the 2-MC block: 30.0 on E0, 45.4 on base).
    etale0_mc2   ~ 0.08618   # log(1 + 0.300^2)
    etalbase_mc2 ~ 0.18743   # log(1 + 0.454^2)

    # ---- 3-hydroxypropionate (3-HP) ----
    le0_hp3    <- log(0.0987); label("E0 for 3-HP (umol/L)")                                       # Table 1, RSE 6.74%
    lbase_hp3  <- log(60.1);   label("Suppressible baseline for 3-HP (umol/L)")                    # Table 1, RSE 14.4%
    lic50_hp3  <- log(37.5);   label("IC50 of PCC for 3-HP suppression (PCC mass per liver tissue)") # Table 1, RSE 12.2%
    propSd_HP3  <- 0.348;       label("Proportional residual SD on 3-HP (fraction)")                # Table 1, RSE 5.56%
    # IIV (Table 1, "IIV (%)" column for the 3-HP block: 51.1 on E0, 73.8 on base).
    etale0_hp3   ~ 0.23199   # log(1 + 0.511^2)
    etalbase_hp3 ~ 0.43485   # log(1 + 0.738^2)

    # ---- C3/C2 carnitine ratio ----
    le0_c3c2   <- log(0.347);  label("E0 for C3/C2 carnitine ratio (umol/L)")                      # Table 1, RSE 6.4%
    lbase_c3c2 <- log(2.31);   label("Suppressible baseline for C3/C2 ratio (umol/L)")             # Table 1, RSE 6.61%
    lic50_c3c2 <- log(32.1);   label("IC50 of PCC for C3/C2 suppression (PCC mass per liver tissue)") # Table 1, RSE 11.7%
    propSd_C3C2 <- 0.360;       label("Proportional residual SD on C3/C2 ratio (fraction)")         # Table 1, RSE 4.79%
    # IIV (Table 1, "IIV (%)" column for the C3/C2 block: 51.9 on E0, 56.4 on base).
    etale0_c3c2   ~ 0.23851  # log(1 + 0.519^2)
    etalbase_c3c2 ~ 0.27618  # log(1 + 0.564^2)
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Body-weight allometric scaling factors (Supp Table S1 part A).
    #    Reference body weight = 0.025 kg (mouse).
    # ----------------------------------------------------------------------
    wt_ratio <- WT / 0.025
    sf_cla   <- wt_ratio^e_wt_cla   # CL12, CL32
    sf_clb   <- wt_ratio^e_wt_clb   # CL23, CL20
    sf_v     <- wt_ratio            # V, V2 (exponent fixed at 1)

    # ----------------------------------------------------------------------
    # 2. Individual PK parameters.
    # ----------------------------------------------------------------------
    cl12 <- exp(lcl12)            * sf_cla
    cl23 <- exp(lcl23)            * sf_clb
    cl32 <- exp(lcl32 + etalcl32) * sf_cla
    cl20 <- exp(lcl20)            * sf_clb
    vc   <- exp(lvc)              * sf_v
    vp   <- exp(lvp)              * sf_v

    # PD parameters (mouse-only fit; not WT-scaled in the published model).
    ke0   <- exp(lke0)
    slope <- exp(lslope)
    kdeg  <- exp(lkdeg)
    kq    <- exp(lkq)

    e0_mc2     <- exp(le0_mc2     + etale0_mc2)
    base_mc2   <- exp(lbase_mc2   + etalbase_mc2)
    ic50_mc2   <- exp(lic50_mc2)

    e0_hp3     <- exp(le0_hp3     + etale0_hp3)
    base_hp3   <- exp(lbase_hp3   + etalbase_hp3)
    ic50_hp3   <- exp(lic50_hp3)

    e0_c3c2    <- exp(le0_c3c2    + etale0_c3c2)
    base_c3c2  <- exp(lbase_c3c2  + etalbase_c3c2)
    ic50_c3c2  <- exp(lic50_c3c2)

    # ----------------------------------------------------------------------
    # 3. Plasma / tissue concentrations from compartment amounts.
    #    central = plasma 1 (A1), peripheral1 = tissue (A2),
    #    peripheral2 = plasma 2 (A3). Both plasma compartments share V.
    # ----------------------------------------------------------------------
    C1  <- central     / vc
    C2  <- peripheral1 / vp
    C3  <- peripheral2 / vc
    C13 <- C1 + C3   # observed plasma mRNA concentration

    # ----------------------------------------------------------------------
    # 4. PK ODE system (Supp Table S1 part A).
    # ----------------------------------------------------------------------
    d/dt(central)     <- -cl12 * C1
    d/dt(peripheral1) <-  cl12 * C1 - cl23 * C2 + cl32 * C3 - cl20 * C2
    d/dt(peripheral2) <-  cl23 * C2 - cl32 * C3

    # ----------------------------------------------------------------------
    # 5. PD effect compartment + PCC protein 2-compartment indirect response
    #    (Supp Table S1 part B). The effect-compartment state is carried as
    #    a concentration (paper's Ce), driven toward C13 with rate ke0.
    #    PCC synthesis is linear in Ce (rate = slope * Ce); see Assumptions
    #    in the vignette for the supplement-equation typo discussion.
    # ----------------------------------------------------------------------
    d/dt(effect) <- ke0 * (C13 - effect)
    d/dt(pcc)    <- slope * effect + kq * (pcc_p - pcc) - kdeg * pcc
    d/dt(pcc_p)  <- kq * (pcc - pcc_p)

    # ----------------------------------------------------------------------
    # 6. Biomarker outputs (Supp Table S1 part C).
    #    biomarker = E0 + base * (1 - Imax * PCC / (IC50 + PCC))
    # ----------------------------------------------------------------------
    pcc_inh <- imax * pcc / (ic50_mc2 + pcc)
    MC2     <- e0_mc2  + base_mc2  * (1 - pcc_inh)

    pcc_inh_hp <- imax * pcc / (ic50_hp3 + pcc)
    HP3        <- e0_hp3 + base_hp3 * (1 - pcc_inh_hp)

    pcc_inh_c3 <- imax * pcc / (ic50_c3c2 + pcc)
    C3C2       <- e0_c3c2 + base_c3c2 * (1 - pcc_inh_c3)

    # ----------------------------------------------------------------------
    # 7. Observation variables and residual error.
    # ----------------------------------------------------------------------
    Cc  <- C13
    PCC <- pcc

    Cc  ~ prop(propSd)
    PCC ~ prop(propSd_PCC)
    MC2 ~ prop(propSd_MC2)
    HP3 ~ prop(propSd_HP3)
    C3C2 ~ prop(propSd_C3C2)
  })
}

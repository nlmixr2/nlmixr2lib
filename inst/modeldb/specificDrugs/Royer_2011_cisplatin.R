Royer_2011_cisplatin <- function() {
  description <- paste(
    "Four-compartment population PK model for platinum after cisplatin",
    "intraperitoneal perioperative chemotherapy (PIPC), with or without",
    "intraperitoneal epinephrine, in adult women with recurrent epithelial",
    "ovarian cancer (Royer 2011). Compartments are peritoneum (paper-",
    "mechanistic IP dosing / sampling site), central (serum ultrafilterable",
    "platinum), peripheral1, and bound (paper-mechanistic protein-bound",
    "plasma platinum). Ultrafiltered platinum transfers peritoneum -> central",
    "at rate IPCL * (peritoneum / IPV), distributes central <-> peripheral at",
    "rate constants k12 = 0.632 /h and k21 = 0.0425 /h (paper notation k23 /",
    "k32), and is eliminated from central at rate CL * Cc. Protein-bound",
    "platinum is formed from central via a Michaelis-Menten binding term",
    "modulated by baseline serum total protein (Vmax * Cc / (KM + Cc) * TPRO,",
    "per Supplementary Data S1) and eliminated first-order at rate kB.",
    "Epinephrine decreases IPCL by 53.1% and increases central volume V by",
    "80.5% (multiplicative fractional coefficients on the typical values).",
    "Ultrafiltered and protein-bound platinum are fitted simultaneously to 316",
    "peritoneum + 577 unbound plasma + 577 bound plasma observations from 55",
    "patients (26 with epinephrine, 29 without)."
  )
  reference <- paste0(
    "Royer B, Kalbacher E, Onteniente S, Jullien V, Montange D, Piedoux S, ",
    "Thiery-Vuillemin A, Delroeux D, Pili-Floury S, Guardiola E, Combe M, ",
    "Muret P, Nerich V, Heyd B, Chauffert B, Kantelip JP, Pivot X. ",
    "Intraperitoneal clearance as a potential biomarker of cisplatin after ",
    "intraperitoneal perioperative chemotherapy: a population pharmacokinetic ",
    "study. Br J Cancer. 2012;106(3):460-467. doi:10.1038/bjc.2011.557"
  )
  vignette <- "Royer_2011_cisplatin"

  # Peritoneum (intraperitoneal cavity) and bound (protein-bound plasma
  # platinum) are paper-mechanistic states not in the canonical register.
  # peritoneum receives the IP dose and delivers unbound platinum to central at
  # IPCL * (peritoneum / IPV); bound is formed via a Michaelis-Menten binding
  # term on central unbound platinum scaled by serum total protein and
  # eliminated first-order (Royer 2011 Figure 1 + Supplementary Data S1). The
  # bound compartment carries a concentration (mg/L) directly, following the
  # Urien 2004 convention that Royer 2011 explicitly cites as the source of
  # this protein-binding formulation.
  paper_specific_compartments <- c("peritoneum", "bound")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_EPI = list(
      description        = "Presence of epinephrine co-administered in the intraperitoneal cisplatin bath",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = cisplatin alone",
      notes              = "1 = epinephrine co-administered in the IP bath (Royer 2011: 1 mg/L n=11, 2 mg/L n=12, or 3 mg/L n=3, dichotomised in the model because plasma Pt concentrations were similar across dose levels per Guardiola 2010). 0 = cisplatin alone (n=29). Enters IPCL and V as multiplicative fractional coefficients: IPCL_i = exp(lipcl + etalipcl) * (1 + e_epi_ipcl * CONMED_EPI); V_i = exp(lvc + etalvc) * (1 + e_epi_vc * CONMED_EPI). The Royer 2011 Table 2 notation `y1 + y2 * EPI` is a proportional-change parameterisation: y2 = -0.531 gives a 53.1% decrease in IPCL, and y4 = +0.805 gives an 80.5% increase in V (both percentages are quoted verbatim in Royer 2011 Results 'Covariates').",
      source_name        = "EPI"
    ),
    TPRO = list(
      description        = "Baseline total serum protein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Median 34.0 g/L (range 16-56) across the pooled 55-patient cohort per Royer 2011 Table 1. Enters the protein-bound Pt formation flux multiplicatively (Supplementary Data S1 equation: PtB = (Vmax * Cc / (Cc + KM)) * Prot). The paper reports TPRO as the baseline serum total protein and uses it as a time-invariant per-subject covariate.",
      source_name        = "PROT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at surgery",
      units       = "years",
      type        = "continuous",
      notes       = "Tested in the covariate screening (Royer 2011 Methods 'Population pharmacokinetic analysis' and Results 'Covariates'); did not survive the joint OFV/IIV-reduction retention criterion and was not retained in the final model."
    ),
    WT = list(
      description = "Actual body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened, not retained (Royer 2011 Methods / Results 'Covariates')."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened, not retained."
    ),
    BSA = list(
      description = "Body surface area (Du Bois and Du Bois 1916 formula)",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened, not retained. BSA is used only downstream of the fit in the interstitial-Pt-penetration equation (3x0 = 3*D*BSA/IPCL) reported in Royer 2011 Discussion, not as a covariate on any structural PK parameter."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened, not retained."
    ),
    SCR = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened, not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault 1976)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened, not retained. Downstream ROC-curve analysis in Royer 2011 evaluated whether the fitted CL discriminates renal toxicity (AUC of ROC 0.514, non-discriminating; Table 3), but CL itself was not modelled as a function of CRCL."
    ),
    PRIP = list(
      description = "Intraperitoneal total protein concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened, not retained. Distinct from TPRO (serum total protein) which is retained on the protein-bound-Pt formation flux."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 1L,
    n_observations = "1470 platinum measurements (316 intraperitoneal + 577 ultrafiltered plasma + 577 protein-bound plasma)",
    age_range      = "25.5-75.0 years (median 58.3)",
    weight_range   = "49-85 kg (median 60.5)",
    height_range   = "150-178 cm (median 161.5)",
    bmi_range      = "17.9-31.2 kg/m^2 (median 23.4)",
    bsa_range      = "1.42-2.01 m^2 (median 1.63; Du Bois-Du Bois formula)",
    lbm_range      = "36.6-56.2 kg (median 43.6; Du Bois formula per Royer 2011 Table 1 footnote a)",
    scr_range      = "32-88 umol/L (median 58.0)",
    crcl_range     = "58.2-182.8 mL/min (median 94.0; Cockcroft-Gault)",
    tpro_range     = "16-56 g/L (median 34.0)",
    sex_female_pct = 100,
    disease_state  = "Recurrent epithelial ovarian cancer confined to the peritoneal cavity (no extra-peritoneal disease), with progression at least 6 months after first-line intravenous platinum-based chemotherapy. Patients had WHO performance status 0-1, life expectancy at least 3 months, and normal haematological / renal / hepatic function at baseline. Cardiac pathology was an exclusion criterion (given the anticipated cardiovascular effects of intraperitoneal epinephrine).",
    dose_range     = "Perioperative intraperitoneal chemotherapy: 3 L physiological saline bath containing cisplatin (typically 60-70 mg/m^2 as elemental cisplatin) delivered during optimal cytoreductive surgery. Two consecutive 1-h baths per procedure; the second bath was shortened to 45 min for the last n=11 epinephrine patients (see notes). Epinephrine (n=26) was given at 1 mg/L (n=11), 2 mg/L (n=12), or 3 mg/L (n=3) in the IP bath; n=29 patients received cisplatin alone.",
    regions        = "France (CHU Besancon multi-centre with CHU Amiens; French phase I trial infrastructure)",
    notes          = "Peritoneal and blood samples at 1, 30, 59 min after the start of each 1-h bath (5, 25, 44 min for the shortened 45-min bath). Additional blood samples at 4, 6, 8, 16, 24 h after IPC; the 16-h samples were discarded for the last n=11 patients because they were inconvenient (0300 h) and not informative. Concomitant IV fluid (Ca-glucuronate, Mg, KCl, NaCl in normal saline) administered for renal-toxicity prevention. Platinum assayed by flameless atomic absorption spectrophotometry (Varian SpectrAA 220Z Zeeman graphite furnace) on ultrafiltered and total plasma fractions. NONMEM VI.2 double precision, FOCE with INTERACTION. 500-resampling nonparametric bootstrap confirmed parameter stability (Royer 2011 Table 2 bootstrap columns)."
  )

  ini({
    # Structural parameters -- Royer 2011 Table 2 'Original data set estimate'
    # column, final model. Typical values hold at CONMED_EPI = 0 (no
    # epinephrine); epinephrine effects enter as multiplicative fractional
    # coefficients on IPCL and V (see e_epi_ipcl / e_epi_vc below).
    lipv  <- log(3.10)                ; label("Intraperitoneal volume of distribution IPV (L)")                                            # Royer 2011 Table 2: IPV = 3.10 L (RSE 3.1%)
    lipcl <- log(4.66)                ; label("Intraperitoneal clearance IPCL intercept y1 (L/h; typical value at CONMED_EPI = 0)")         # Royer 2011 Table 2: IPCL y1 = 4.66 L/h (RSE 4.3%)
    lcl   <- log(9.63)                ; label("Central (serum) clearance CL (L/h)")                                                        # Royer 2011 Table 2: CL = 9.63 L/h (RSE 5.2%)
    lvc   <- log(21.4)                ; label("Central volume V intercept y3 (L; typical value at CONMED_EPI = 0)")                         # Royer 2011 Table 2: V y3 = 21.4 L (RSE 5.9%)
    lk12  <- log(0.632)               ; label("Central-to-peripheral1 rate constant k12 (1/h; paper notation k23)")                        # Royer 2011 Table 2: k23 = 0.632 /h (RSE 3.7%)
    lk21  <- log(0.0425)              ; label("Peripheral1-to-central rate constant k21 (1/h; paper notation k32)")                        # Royer 2011 Table 2: k32 = 0.0425 /h (RSE 4.4%)
    lvmax <- log(0.0123)              ; label("Michaelis-Menten Vmax of protein-bound Pt formation flux (see vignette Errata for units)")   # Royer 2011 Table 2: Vmax = 0.0123 (paper table units 'mg/L'; see vignette Assumptions and deviations)
    lkm   <- log(2.00)                ; label("Michaelis-Menten constant KM (mg/L)")                                                       # Royer 2011 Table 2: KM = 2.00 mg/L (RSE 10.1%)
    lkb   <- log(0.382)               ; label("First-order elimination rate constant of protein-bound Pt kB (1/h; see vignette Errata)")   # Royer 2011 Table 2: kB = 0.382 (paper table units 'l h-1'; treated as 1/h for dimensional consistency, see vignette)

    # Epinephrine effect (multiplicative fractional coefficient). The Royer
    # 2011 Table 2 notation `y1 + y2 * EPI` is a proportional-change form:
    # y2 = -0.531 gives a 53.1% decrease in IPCL when EPI = 1, and y4 = 0.805
    # gives an 80.5% increase in V when EPI = 1. Both percentages are quoted
    # verbatim in Royer 2011 Results 'Covariates' and used here to confirm
    # the multiplicative rather than additive interpretation.
    e_epi_ipcl <- -0.531              ; label("Fractional epinephrine effect on IPCL (unitless; -0.531 = 53.1% decrease at CONMED_EPI = 1)") # Royer 2011 Table 2: y2 = -0.531 (RSE 6.9%)
    e_epi_vc   <-  0.805              ; label("Fractional epinephrine effect on V (unitless; +0.805 = 80.5% increase at CONMED_EPI = 1)")   # Royer 2011 Table 2: y4 = 0.805 (RSE 26.3%)

    # Interindividual variability (log-normal; omega^2 = log(1 + CV^2)). Royer
    # 2011 Table 2 reports IIV as %CV; IIV on k12, k21, KM, and kB could not
    # be estimated (Royer 2011 Results 'Patient population and structural PK
    # model'). Correlation between V (lvc) and Vmax (lvmax) is r = -0.124 per
    # Table 2 row 'r V/Vmax'; covariance = -0.124 * sqrt(var_V * var_Vmax) =
    # -0.124 * sqrt(0.06735 * 0.23361) = -0.01555.
    etalipv  ~ 0.03808                                                                                                                     # log(1 + 0.197^2); Royer 2011 Table 2: IIV_IPV  = 19.7% CV (RSE 25.9%)
    etalipcl ~ 0.04853                                                                                                                     # log(1 + 0.223^2); Royer 2011 Table 2: IIV_IPCL = 22.3% CV (RSE 19.0%)
    etalcl   ~ 0.14444                                                                                                                     # log(1 + 0.394^2); Royer 2011 Table 2: IIV_CL   = 39.4% CV (RSE 28.2%)
    etalvc + etalvmax ~ c(0.06735,
                          -0.01555, 0.23361)                                                                                               # var_V = log(1+0.264^2) = 0.06735; cov = -0.01555; var_Vmax = log(1+0.513^2) = 0.23361; Table 2 rows IIV_V (26.4% CV), IIV_Bmax (51.3% CV), r_V/Vmax (-0.124)

    # Residual error. Royer 2011 Table 2: one shared proportional error (e1 =
    # 17.8% CV) applied to all three outputs (peritoneum, unbound plasma,
    # bound plasma) and one additive error (e2 = 0.098 mg/L) applied to
    # peritoneum concentrations only (Royer 2011 Results paragraph following
    # Table 2: "The error model includes both proportional and additive
    # models, but the latter was only applied to IP concentrations."). The
    # bare `propSd` parameter is applied to the canonical parent output `Cc`
    # (unbound plasma), and per-output `propSd_<output>` aliases (both set
    # equal to `propSd`) are applied to Cip and Cbound so the model preserves
    # the paper's shared-value semantics while satisfying the
    # `checkModelConventions()` per-output residual-name pattern.
    propSd         <- 0.178          ; label("Proportional residual error CV, unbound plasma (Cc)")                                       # Royer 2011 Table 2: e1 = 17.8% CV (RSE 9.6%)
    propSd_Cip     <- 0.178          ; label("Proportional residual error CV, peritoneum concentration (shared value with propSd)")        # Royer 2011 Table 2: e1 (shared)
    propSd_Cbound  <- 0.178          ; label("Proportional residual error CV, protein-bound plasma Pt (shared value with propSd)")         # Royer 2011 Table 2: e1 (shared)
    addSd_Cip      <- 0.098          ; label("Additive residual SD, peritoneum concentration (mg/L)")                                     # Royer 2011 Table 2: e2 = 0.098 mg/L (RSE 7.0%)
  })

  model({
    # Individual pharmacokinetic parameters. CONMED_EPI (0/1) enters IPCL and
    # V as multiplicative fractional coefficients (see e_epi_ipcl /
    # e_epi_vc). k12, k21, KM, and kB have no interindividual variability
    # (fixed across the population per Royer 2011 Results 'Patient
    # population and structural PK model').
    ipv   <- exp(lipv  + etalipv)
    ipcl  <- exp(lipcl + etalipcl) * (1 + e_epi_ipcl * CONMED_EPI)
    cl    <- exp(lcl   + etalcl)
    vc    <- exp(lvc   + etalvc)   * (1 + e_epi_vc   * CONMED_EPI)
    k12   <- exp(lk12)
    k21   <- exp(lk21)
    vmax  <- exp(lvmax + etalvmax)
    km    <- exp(lkm)
    kb    <- exp(lkb)

    # Observation-variable concentrations (mg/L). Cip is peritoneum-cavity
    # cisplatin concentration; Cc is central unbound-plasma platinum
    # concentration; Cbound is protein-bound plasma platinum concentration.
    # The `bound` compartment carries concentration units directly (see the
    # paper_specific_compartments note above and Urien 2004 which Royer 2011
    # cites).
    Cip    <- peritoneum / ipv
    Cc     <- central    / vc
    Cbound <- bound

    # ODE system -- Royer 2011 Figure 1 scheme. Ultrafiltered Pt transfers
    # from peritoneum to central at flux IPCL * Cip, distributes central <->
    # peripheral1 at rate constants k12 and k21, and is eliminated from
    # central at flux CL * Cc. Protein-bound Pt is formed from central Cc by
    # the Michaelis-Menten binding term (Vmax * Cc / (KM + Cc)) scaled by
    # baseline serum total protein TPRO per Supplementary Data S1 equation
    # `PtB = (Vmax * PtUf / (PtUf + KM)) * Prot`, and eliminated first-order
    # at rate kB. Following Urien 2004 (explicitly cited in Royer 2011), the
    # bound-formation flux is NOT subtracted from central: the paper's CL is
    # the apparent total systemic clearance of Uf platinum and already
    # encompasses all elimination pathways (including covalent protein
    # binding) as an apparent parameter.
    d/dt(peritoneum)  <- -ipcl * Cip
    d/dt(central)     <-  ipcl * Cip - cl * Cc - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12  * central - k21 * peripheral1
    d/dt(bound)       <-  (vmax * Cc / (km + Cc)) * TPRO - kb * bound

    # Observation model. e1 (proportional, 17.8% CV) is shared across the
    # three outputs (Cip, Cc, Cbound); e2 (additive, 0.098 mg/L) is applied
    # to peritoneum concentrations only.
    Cip    ~ prop(propSd_Cip) + add(addSd_Cip)
    Cc     ~ prop(propSd)
    Cbound ~ prop(propSd_Cbound)
  })
}

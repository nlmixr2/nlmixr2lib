Sampson_2014_azithromycin <- function() {
  description <- "Four-compartment mamillary population PK model for oral azithromycin simultaneously describing concentrations in whole blood, peripheral blood mononuclear cells (PBMCs), and polymorphonuclear cells (PMNs) in healthy adults (Sampson 2014). First-order absorption with lag; unidirectional flow from central to PBMC and to PMN compartments; bidirectional flow between central and a peripheral tissue compartment; elimination from central, PBMC, and PMN compartments. The observed whole-blood concentration is a weighted sum of plasma, PBMC, and PMN concentrations."
  reference <- "Sampson MR, Dumitrescu TP, Brouwer KLR, Schmith VD. Population pharmacokinetics of azithromycin in whole blood, peripheral blood mononuclear cells, and polymorphonuclear cells in healthy adults. CPT Pharmacometrics Syst Pharmacol. 2014;3(3):e103. doi:10.1038/psp.2013.80"
  vignette <- "Sampson_2014_azithromycin"
  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "mg/L"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    age_range      = "21-63 years",
    age_median     = "48.5 years",
    weight_range   = ">=50 kg (inclusion criterion)",
    bmi_range      = "21.4-28.2 kg/m^2",
    sex_female_pct = 40,
    race_ethnicity = c(Caucasian = 80),
    disease_state  = "Healthy adults",
    dose_range     = "250 mg (n=10) or 1,000 mg (n=10) oral single dose",
    regions        = "United Kingdom (Cambridge, single site)",
    trial_id       = "NCT01416350",
    notes          = "Single-dose study; 269 blood, 227 PBMC, and 239 PMN observations collected predose and at 1, 2, 3, 4, 6, 9, 12, 16, 24, 48, 96, 144, 240, 336, and 504 h postdose (PBMC/PMN measured at all timepoints except 3 h and 240 h). Baseline demographics from Results section paragraph 1; sampling schedule from Methods (Study design). No covariates were retained in the final model."
  )

  ini({
    # Structural parameters - Sampson 2014 Table 1 (model estimate column).
    # Central compartment (Comp1): plasma. Whole-blood observations are
    # modelled as a weighted sum of plasma, PBMC, and PMN concentrations.
    # V1/F and clearances are apparent (oral, F unidentifiable). Cellular
    # volumes V2/F and V3/F are small because cellular concentrations are
    # >100-fold higher than blood/plasma. Cellular parameter names use
    # the no-underscore form (lqpbmc, lqpmn, lclpbmc, lclpmn) so the
    # shared eta on CL12 and CL13 (etalqpbmc_qpmn) is recognised by
    # checkModelConventions() as a valid shared-eta suffix.
    ltlag    <- log(0.41);     label("Absorption lag time (h)")                                    # Table 1
    lka      <- log(0.53);     label("First-order absorption rate constant (1/h)")                 # Table 1
    lvc      <- log(336);      label("Apparent central (plasma) volume V1/F (L)")                  # Table 1
    lvpbmc   <- log(0.62);     label("Apparent PBMC volume V2/F (L)")                              # Table 1
    lvpmn    <- log(2.96);     label("Apparent PMN volume V3/F (L)")                               # Table 1
    lvp      <- log(4597);     label("Apparent peripheral tissue volume V4/F (L)")                 # Table 1
    lqpbmc   <- log(9.0);      label("Apparent intercompartmental clearance central -> PBMC CL12/F (L/h)")  # Table 1
    lqpmn    <- log(26.7);     label("Apparent intercompartmental clearance central -> PMN CL13/F (L/h)")   # Table 1
    lq       <- log(73.2);     label("Apparent intercompartmental clearance central -> peripheral tissue CL14/F (L/h)")  # Table 1
    lcl      <- log(67.3);     label("Apparent central elimination clearance CL1/F (L/h)")         # Table 1
    lclpbmc  <- log(0.0091);   label("Apparent PBMC elimination clearance CL2/F (L/h)")            # Table 1
    lclpmn   <- log(0.026);    label("Apparent PMN elimination clearance CL3/F (L/h)")             # Table 1

    # Blood observation mixing coefficients (unitless): [blood] = a*Cplasma + b*Cpbmc + c*Cpmn.
    # A and B were estimated empirically; C was fixed at B/1,000 in the final model.
    a_blood  <- 0.51;          label("Blood mixing coefficient for plasma (unitless)")            # Table 1
    b_blood  <- 0.0016;        label("Blood mixing coefficient for PBMC (unitless)")              # Table 1
    # CL41 = CL14/2 was a fixed structural relationship in the source model
    # (Table 1 footnote). Implemented inside model() as q41 = q/2.

    # Interindividual variability - Sampson 2014 Table 1 (CV%, exponential model)
    # Conversion: omega^2 = log(CV^2 + 1) for exponential (log-normal) IIV.
    etalka         ~ 0.1554   # CV 41%:  log(0.41^2 + 1)   (Table 1, eta_Ka)
    etalvc         ~ 0.9117   # CV 122%: log(1.22^2 + 1)   (Table 1, eta_V1/F)
    etalvpbmc      ~ 0.2313   # CV 51%:  log(0.51^2 + 1)   (Table 1, eta_V2/F)
    etalvpmn       ~ 0.2476   # CV 53%:  log(0.53^2 + 1)   (Table 1, eta_V3/F)
    etalcl         ~ 0.8329   # CV 114%: log(1.14^2 + 1)   (Table 1, eta_CL1/F)
    etalqpbmc_qpmn ~ 0.4463   # CV 75%:  log(0.75^2 + 1)   (Table 1, shared eta on CL12 and CL13)

    # Residual error - Sampson 2014 Table 1 (proportional, CV%)
    propSd_Cblood <- 0.47;  label("Proportional residual error, whole-blood concentration (fraction)")  # Table 1
    propSd_Cpbmc  <- 0.74;  label("Proportional residual error, PBMC concentration (fraction)")         # Table 1
    propSd_Cpmn   <- 0.64;  label("Proportional residual error, PMN concentration (fraction)")          # Table 1
  })

  model({
    # Individual parameters
    tlag    <- exp(ltlag)
    ka      <- exp(lka      + etalka)
    vc      <- exp(lvc      + etalvc)
    vpbmc   <- exp(lvpbmc   + etalvpbmc)
    vpmn    <- exp(lvpmn    + etalvpmn)
    vp      <- exp(lvp)
    cl      <- exp(lcl      + etalcl)
    clpbmc  <- exp(lclpbmc)
    clpmn   <- exp(lclpmn)
    # Shared IIV applied to both qpbmc (CL12) and qpmn (CL13); ref. Table 1
    # footnote: "eta_CL12-CL13/F is the shared eta estimate for CL12 and CL13".
    qpbmc   <- exp(lqpbmc  + etalqpbmc_qpmn)
    qpmn    <- exp(lqpmn   + etalqpbmc_qpmn)
    q       <- exp(lq)
    # Fixed structural relationship: CL41 = CL14 / 2 (Table 1 footnote).
    q41     <- q / 2
    # Fixed structural relationship: C = B / 1000 (Table 1 footnote).
    c_blood <- b_blood / 1000

    # Micro-constants
    k_pbmc_out <- clpbmc / vpbmc
    k_pmn_out  <- clpmn  / vpmn
    k_per_out  <- q41    / vp

    # Four-compartment mamillary ODE system (Figure 2):
    # depot -> central (ka, with tlag);
    # unidirectional central -> pbmc (qpbmc) and central -> pmn (qpmn);
    # bidirectional central <-> peripheral1 (q in, q41 = q/2 back);
    # elimination from central (cl), pbmc (clpbmc), and pmn (clpmn).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl + qpbmc + qpmn + q) * central / vc +
                          k_per_out * peripheral1
    d/dt(pbmc)        <-  qpbmc * central / vc - k_pbmc_out * pbmc
    d/dt(pmn)         <-  qpmn  * central / vc - k_pmn_out  * pmn
    d/dt(peripheral1) <-  q     * central / vc - k_per_out  * peripheral1

    # Absorption lag time on the depot compartment
    lag(depot) <- tlag

    # Observed concentrations
    Cplasma <- central / vc
    Cpbmc   <- pbmc    / vpbmc
    Cpmn    <- pmn     / vpmn
    # Whole-blood concentration is the weighted sum of plasma, PBMC, and PMN
    # concentrations (paper: [blood] = A*[Comp1] + B*[Comp2] + C*[Comp3], where
    # A, B were estimated and C = B/1000 was fixed). Reported in mg/L.
    Cblood  <- a_blood * Cplasma + b_blood * Cpbmc + c_blood * Cpmn

    Cblood ~ prop(propSd_Cblood)
    Cpbmc  ~ prop(propSd_Cpbmc)
    Cpmn   ~ prop(propSd_Cpmn)
  })
}

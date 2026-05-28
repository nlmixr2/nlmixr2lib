Johnson_2011_olanzapine_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Mechanism-based hybrid physiology-based population PK-PD model for",
    "olanzapine and striatal dopamine D2 receptor occupancy (D2RO) in",
    "rats (Johnson 2011). Plasma PK is a 2-compartment model fitted",
    "across IP, SC, and IV routes in Wistar / Sprague-Dawley rats",
    "(single dose 0.01-40 mg/kg, pooled from 12 studies, n = 283);",
    "the absorption rate constant was not estimable, so all routes",
    "deposit drug directly into the central compartment, and the",
    "intraperitoneal bioavailability FIP is estimated (about 64%) with",
    "an 87% CV log-normal IIV. SC and IV bioavailability are fixed at 1.",
    "The brain submodel adds a brain-vascular compartment (Vbv, fed by",
    "cerebral blood flow CLbv from systemic central) and a brain-",
    "extravascular compartment (Vbev, fed across the BBB by an estimated",
    "clearance CLbev applied to the unbound concentration on each side",
    "via fixed fu_plasma and fu_brain).",
    "D2 receptor occupancy in striatum is the reduced model published by",
    "the authors (Bmax dropped per their sensitivity analysis):",
    "dRO/dt = kon * Cfree_bev * (1 - RO) - koff * RO, with kon = koff/Kd",
    "and Cfree_bev = fu_brain * (Cbev in nM), so the binding kinetics are",
    "driven by the free brain-extravascular concentration converted to",
    "nM via the olanzapine molecular weight (312.43 g/mol). All",
    "structural parameters are body-weight-normalised (per kg)."
  )
  reference <- paste(
    "Johnson M, Kozielska M, Pilla Reddy V, Vermeulen A, Li C, Grimwood S,",
    "de Greef R, Groothuis GMM, Danhof M, Proost JH. Mechanism-Based",
    "Pharmacokinetic-Pharmacodynamic Modeling of the Dopamine D2 Receptor",
    "Occupancy of Olanzapine in Rats.",
    "Pharm Res. 2011;28(10):2490-2504.",
    "doi:10.1007/s11095-011-0477-7.",
    sep = " "
  )
  vignette <- "Johnson_2011_olanzapine_rat"
  units <- list(
    time          = "h",
    dosing        = "mg/kg",
    concentration = "mg/L"
  )

  covariateData <- list(
    ROUTE_IP = list(
      description        = paste(
        "Intraperitoneal-route indicator: 1 = the dose record is",
        "intraperitoneal (IP), 0 = subcutaneous (SC) or intravenous (IV)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (SC or IV; bioavailability fixed at 1)",
      notes              = paste(
        "Per-dose-record covariate. Selects the intraperitoneal",
        "bioavailability FIP (about 0.636 with 87% CV log-normal IIV)",
        "when ROUTE_IP = 1; the encoding",
        "f(central) <- exp(ROUTE_IP * (lfip + etalfip))",
        "collapses to F = 1 when ROUTE_IP = 0 because exp(0) = 1, so SC",
        "and IV doses inherit complete bioavailability without IIV on F.",
        "All routes deposit drug directly into the central compartment",
        "because the absorption rate constant was not estimable from the",
        "Johnson 2011 dataset (Results, Population Plasma Pharmacokinetics)."
      ),
      source_name        = paste(
        "(implicit; the per-study route is reported in Table I of the",
        "source paper)"
      )
    )
  )

  population <- list(
    species        = "rat (Wistar or Sprague-Dawley)",
    n_subjects     = 283L,
    n_studies      = 12L,
    age_range      = "Adult (specific age not reported in the source publication)",
    weight_range   = "Not reported individually (per-kg structural parameters)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy",
    dose_range     = paste(
      "Single dose 0.01-40 mg/kg via intraperitoneal (IP), subcutaneous",
      "(SC), or intravenous (IV) routes. Per Table I: IP studies 1-6b",
      "(0.01-30 mg/kg); SC studies 7-11 (0.04-40 mg/kg); IV study 12",
      "(2.5 mg/kg)."
    ),
    regions        = "Pooled industrial datasets from three sponsors",
    notes          = paste(
      "Pooled dataset of plasma concentration (PC), total brain",
      "concentration (BC) and striatal dopamine D2 receptor occupancy",
      "(RO) measurements from 12 single-dose studies contributed to the",
      "TI Pharma mechanism-based PK-PD platform by Janssen Research and",
      "Development (Belgium), Merck Sharp and Dohme (The Netherlands)",
      "and Pfizer Global Research and Development (USA). Strains were",
      "Wistar or Sprague-Dawley; the source paper does not break down",
      "demographics by study (sex, weight, age). D2RO was measured by",
      "in vivo binding with [3H]raclopride (studies 2, 5, 6a, 6b, 7, 9,",
      "10) or ex vivo binding with [125I]sulpiride (studies 8, 11). See",
      "Table I and Methods, Data Management."
    )
  )

  ini({
    # ----- Plasma PK -- Johnson 2011 Table III (original dataset point estimates) -----
    # Methods, Modeling Tools + Population Pharmacokinetic Analysis: NONMEM
    # FOCE; ADVAN9 user-defined ODEs; dose enters central directly because
    # ka was not estimable (Results, Population Plasma Pharmacokinetics).
    lvc  <- log(4.22)
    label("Plasma central volume of distribution V1 (L/kg)")          # Table III: V1 = 4.22 (RSE 9%)

    lcl  <- log(3.21)
    label("Plasma systemic clearance CL (L/h/kg)")                    # Table III: CL = 3.21 (RSE 9%)

    lvp  <- log(2.23)
    label("Plasma peripheral volume of distribution V2 (L/kg)")       # Table III: V2 = 2.23 (RSE 14%)

    lq   <- log(1.70)
    label("Plasma inter-compartmental clearance Q (L/h/kg)")          # Table III: Q = 1.70 (RSE 30%)

    lfip <- log(0.636)
    label("Log bioavailability for intraperitoneal route FIP (unitless)") # Table III: FIP = 0.636 (RSE 12%); FSC = FIV = 1 (fixed per Results)

    # ----- Brain PK + PD -- Johnson 2011 Table V (reduced PBPKPD final estimates) -----
    # The final model is the *reduced* PBPKPD: Bmax is dropped from the full
    # model per the Sensitivity Analysis (Results, Hybrid Physiology-Based
    # PK-PD Model). The remaining estimated brain parameters are CLbev, Kd
    # and koff; kon is derived as koff/Kd.
    lcl_bev <- log(0.394)
    label("Brain-extravascular clearance CLbev (L/h/kg)")             # Table V: CLbev = 0.394 (RSE 15%)

    lkd     <- log(14.7)
    label("D2-receptor equilibrium dissociation constant Kd (nM)")    # Table V: Kd = 14.7 (RSE 8%)

    lkoff   <- log(2.62)
    label("D2-receptor dissociation rate constant koff (1/h)")        # Table V: koff = 2.62 (RSE 24%)

    # ----- Fixed physiological / binding parameters -----
    # Methods, Hybrid Physiology-Based PK-PD Model:
    #   "The volumes of brain-vascular (Vbv) and brain-extravascular (Vbev)
    #    compartments were assumed to be equal to the physiological values
    #    in the rat: 0.00024 L/kg and 0.00656 L/kg, respectively (18).
    #    The clearance from the brain-vascular compartment (CLbv) was
    #    assumed to be equal to the cerebral blood flow in rats, which is
    #    0.312 L/h/kg (18).  ... The unbound fraction of olanzapine in
    #    plasma (fu_plasma) and brain (fu_brain) was fixed to the values
    #    obtained from literature: 0.23 and 0.034, respectively (19)."
    lvbv      <- fixed(log(0.00024))
    label("Brain-vascular volume Vbv (L/kg; rat physiology, fixed)")          # Methods: Vbv = 0.00024 L/kg (ref. 18)

    lvbev     <- fixed(log(0.00656))
    label("Brain-extravascular volume Vbev (L/kg; rat physiology, fixed)")    # Methods: Vbev = 0.00656 L/kg (ref. 18)

    lcl_bv    <- fixed(log(0.312))
    label("Brain-vascular clearance CLbv = cerebral blood flow (L/h/kg; rat physiology, fixed)") # Methods: CLbv = 0.312 L/h/kg (ref. 18)

    fu_plasma <- fixed(0.23)
    label("Unbound fraction of olanzapine in plasma fu_plasma (rat; fixed)")  # Methods: fu_plasma = 0.23 (ref. 19)

    fu_brain  <- fixed(0.034)
    label("Unbound fraction of olanzapine in brain fu_brain (rat; fixed)")    # Methods: fu_brain = 0.034 (ref. 19)

    # ----- IIV (Johnson 2011 Table III; Methods convention) -----
    # Methods: "Inter-animal variability is expressed as percent coefficient
    # of variation which is the square root of omega^2 * 100", i.e.
    # omega^2 = (CV/100)^2 (the small-CV approximation -- distinct from the
    # exact log-normal formula log(1 + CV^2)). The encoded variances follow
    # the paper's convention literally so that the omega^2 round-trip to the
    # paper's reported 56% / 87% CV.
    #   IAV-CL  56% CV  -> omega^2 = 0.56^2 = 0.3136
    #   IAV-F1  87% CV  -> omega^2 = 0.87^2 = 0.7569
    # Brain PK and PD parameters have no IIV ("Due to the scarcity of data,
    # no inter-animal variability was assumed in the PBPKPD model").
    etalcl  ~ 0.3136
    # ^ Table III: IAV-CL = 56% CV (RSE 19%)

    etalfip ~ 0.7569
    # ^ Table III: IAV-F1 (IP only; SC/IV F=1 fixed) = 87% CV (RSE 18%)

    # ----- Residual error -----
    # Methods, Modeling Tools:
    #   ln(Yobs) = ln(Ypred) + eps  for plasma + brain olanzapine
    #     (log-additive error; here mapped to nlmixr2 prop(propSd) which is
    #     the small-CV approximation of the log-normal error structure
    #     -- the encoded SD equals the log-additive sigma published)
    #   Yobs    = Ypred + eps        for D2RO (additive on the 0-1 D2RO scale)
    # Plasma propSd = 0.141 matches Table III "Proportional error 0.141".
    # Brain  propSd_Cbrain = 0.479 matches Table V "Proportional error (BC)
    #   0.479"; sqrt-of-variance is NOT taken because the published value
    #   already matches the Discussion narrative "high residual variability
    #   (48%) for the brain concentrations".
    # D2RO addSd_D2RO = 0.136 matches Table V "Additive Error (D2RO) 0.136".
    propSd        <- 0.141
    label("Plasma olanzapine proportional residual SD (fraction; log-additive)")  # Table III: 0.141 (RSE 7%)

    propSd_Cbrain <- 0.479
    label("Brain olanzapine proportional residual SD (fraction; log-additive)")   # Table V:  0.479 (RSE 6%); Discussion confirms 48% residual variability

    addSd_D2RO    <- 0.136
    label("D2 receptor occupancy additive residual SD (fraction, 0-1 scale)")     # Table V:  0.136 (RSE 5%)
  })

  model({
    # ----- Constants -----
    # Olanzapine molecular weight (g/mol) for the mg/L -> nM concentration
    # conversion used by the receptor-binding ODE. Appendix 1 writes the
    # conversion as CEV [nM] = (A4/V4) / (MW/1000), which is dimensionally
    # consistent only when A4/V4 is in ug/L (i.e., NONMEM running with
    # amounts in ug and volumes in L); in this file amounts are in mg and
    # volumes in L, so the equivalent dimensionally-correct factor is
    # 1e6 / MW (1 mg/L = 1e6 nmol/L when divided by MW in g/mol). The
    # published kon = koff / Kd = 0.178 nM^-1 h^-1 is preserved.
    mw_olanzapine <- 312.43

    # ----- Individual structural parameters -----
    vc      <- exp(lvc)
    cl      <- exp(lcl + etalcl)
    vp      <- exp(lvp)
    q       <- exp(lq)
    vbv     <- exp(lvbv)
    vbev    <- exp(lvbev)
    cl_bv   <- exp(lcl_bv)
    cl_bev  <- exp(lcl_bev)
    kd      <- exp(lkd)
    koff    <- exp(lkoff)
    kon     <- koff / kd                       # nM^-1 h^-1; derived (Methods)

    # ----- Working concentrations -----
    Cp           <- central             / vc   # plasma central       (mg/L)
    Ct           <- peripheral1         / vp   # plasma peripheral    (mg/L)
    Cbv          <- brain_vascular      / vbv  # brain vascular       (mg/L)
    Cbev         <- brain_extravascular / vbev # brain extravascular  (mg/L)
    Cbev_nM      <- Cbev * 1e6 / mw_olanzapine     # mg/L -> nM
    Cbev_free_nM <- fu_brain * Cbev_nM             # free brain conc (nM)

    # ----- ODEs (Appendix 1, reduced model) -----
    # Central: linear CL elimination + Q exchange with peripheral + cerebral
    #   blood flow exchange with brain vascular (CLbv carries both forward
    #   and reverse flux because CLbv represents the bidirectional cerebral
    #   blood flow rather than a one-way clearance).
    d/dt(central)             <- -cl    * Cp -
                                  q     * (Cp  - Ct) -
                                  cl_bv * (Cp  - Cbv)
    d/dt(peripheral1)         <-  q     * (Cp  - Ct)

    # Brain-vascular: cerebral-blood-flow in/out + BBB-clearance to brain
    #   extravascular (CLbev applied to the *unbound* concentration on each
    #   side because "only the unbound olanzapine in this intravascular
    #   compartment crosses the BBB", Methods).
    d/dt(brain_vascular)      <-  cl_bv  * (Cp - Cbv) -
                                  cl_bev * (fu_plasma * Cbv - fu_brain * Cbev)

    # Brain-extravascular: receives the BBB flux only.
    d/dt(brain_extravascular) <-  cl_bev * (fu_plasma * Cbv - fu_brain * Cbev)

    # D2 receptor occupancy (fractional, 0-1). Reduced model: the
    # striatum-free compartment is dropped (the binding to D2 receptors is
    # assumed not to affect the free drug concentration in the brain) so
    # the receptor-occupancy ODE is driven by the free brain-extravascular
    # concentration directly.
    d/dt(effect)              <-  kon  * Cbev_free_nM * (1 - effect) -
                                  koff * effect

    # ----- Bioavailability (route-dependent) -----
    # The encoding f(central) = exp(ROUTE_IP * (lfip + etalfip)) collapses
    # to 1 when ROUTE_IP = 0 (SC or IV) and to exp(lfip + etalfip) when
    # ROUTE_IP = 1 (IP). Drug enters central directly because ka was not
    # estimable.
    f(central) <- exp(ROUTE_IP * (lfip + etalfip))

    # ----- Observations -----
    Cc     <- Cp        # plasma olanzapine (mg/L)
    Cbrain <- Cbev      # brain (extravascular) olanzapine (mg/L)
    D2RO   <- effect    # fractional D2 receptor occupancy (0-1)

    Cc     ~ prop(propSd)
    Cbrain ~ prop(propSd_Cbrain)
    D2RO   ~ add(addSd_D2RO)
  })
}

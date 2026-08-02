Henin_2012_diclofenac <- function() {
  description <- paste(
    "Semi-mechanistic Gastro-Intestinal Transit Time (GITT) absorption",
    "model applied to enteric-coated diclofenac (Henin 2012, The AAPS",
    "Journal 14:155-163). The intact tablet moves through the stomach",
    "(no drug release due to enteric coating) into the proximal small",
    "intestine, distal small intestine, and colon; drug release and",
    "absorption are collapsed into a single first-order rate per GI",
    "region because the drug release and absorption phenomena were",
    "indistinguishable in the observed plasma-only dataset. The",
    "effective absorption rate from the intact tablet is a step-",
    "function-weighted combination KA_PSI (proximal SI, 1.06 1/h),",
    "KA_DSI (distal SI, 8.64 1/h), and KA_Col (colon, 1 1/h FIXED)",
    "gated by tablet-position STEP functions with per-subject",
    "inflection points IP_APSI (gastric emptying), IP_PSI_DSI, and",
    "IP_DSI_C (population distributions from Table II). Total absorbed",
    "fraction is limited by the gut-wall bioavailability FA = 0.61",
    "(Table III, IIV 8% CV). Disposition is a three-compartment model",
    "with population parameters allometrically scaled to a 70 kg",
    "reference (CL = 16.5 L/h/70kg^0.75, V1 = 3.68 L/70kg, Q2 = 1.75",
    "L/h/70kg^0.75, V2 = 7.48 L/70kg, Q3 = 7.21 L/h/70kg^0.75, V3 =",
    "3.79 L/70kg). Combined additive (10 nmol/L) + proportional (12.6%)",
    "residual error. This extraction encodes only the 'no return to",
    "fundus' subpopulation and treats the fundus + antrum transit as a",
    "single stomach-to-PSI event, appropriate for an enteric-coated",
    "tablet whose contents do not release in the stomach."
  )
  reference <- paste(
    "Henin E, Bergstrand M, Standing JF, Karlsson MO. A mechanism-Based",
    "Approach for Absorption Modeling: The Gastro-Intestinal Transit",
    "Time (GITT) Model. The AAPS Journal. 2012;14(2):155-163.",
    "doi:10.1208/s12248-012-9324-y."
  )
  vignette <- "Henin_2012_gastrointestinal_transit_time_absorption"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight; used for allometric scaling of the disposition parameters to a 70 kg reference.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Continuous covariate. Reference weight is 70 kg (Table III uses per-70kg allometric normalization). Adult healthy volunteers in the bioequivalence source study.",
      source_name        = "WT"
    ),
    IP_APSI = list(
      description        = "Individual inflection-point time (h) for the enteric-coated tablet transit from stomach to proximal small intestine (STEP value 0.5). Population mean stomach transit time ~ 2 h (paper Results, 'Application of GITT Model to Diclofenac Data': range 1.5-3 h across the studied individuals). Fasted-condition MRT was used because the bioequivalence source study was fasted; the paper reports the population MRT for stomach (fundus + antrum) as a single value of about 2 h.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject covariate. In the Henin 2012 diclofenac fit the stomach residence variability was estimated from the diclofenac plasma data (paper Results reports 19% shrinkage on the stomach-residence random effect). The published population MRT was 2 h; the vignette shows an example of sampling IP_APSI from a log-normal (MRT_stomach = 2 h, sigma corresponding to ~ 100% CV per Table II combined stomach entry).",
      source_name        = "IP_APSI"
    ),
    IP_PSI_DSI = list(
      description        = "Individual inflection-point time (h) for tablet transition from proximal to distal small intestine (STEP value 0.5). Sampled as MRT_psi * exp(eta) with MRT_psi = 1.17 h and eta ~ N(0, VRT_psi) with VRT_psi = 1.37 h^2 per Table II (fasted / fixed from Bergstrand 2009 upstream Markov fit).",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject covariate. Note: IP_PSI_DSI is relative to time 0 (dose), so the effective PSI residence time is IP_PSI_DSI - IP_APSI. Table II VRT = 1.37 gives CV = 50% on the residence-time distribution.",
      source_name        = "IP_PSI_DSI"
    ),
    IP_DSI_C = list(
      description        = "Individual inflection-point time (h) for tablet transition from distal small intestine to colon (STEP value 0.5). Sampled as MRT_dsi * exp(eta) with MRT_dsi = 1.22 h and eta ~ N(0, VRT_dsi) with VRT_dsi = 1.48 h^2 per Table II.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject covariate. IP_DSI_C is time 0-referenced; effective DSI residence time is IP_DSI_C - IP_PSI_DSI. Table II CV = 58%.",
      source_name        = "IP_DSI_C"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    disease_state  = "Healthy adult volunteers",
    dose_range     = "Single 50 mg oral enteric-coated diclofenac tablet under fasting conditions.",
    regions        = "Ireland (Irish Medicines Board approved; unpublished study by Rosemont Pharmaceuticals Ltd, UK).",
    notes          = "Bioequivalence sub-study with three formulations (enteric-coated tablet, soluble tablet, suspension); the GITT model was applied to the enteric-coated arm only. Disposition parameters were fixed from a separate IV pediatric study (Korpela 1990; 10 children age 4-6 y receiving 0.5 mg/kg IV over 5-15 min) with allometric extrapolation to 70 kg adult weight. Absorption parameters (KA_PSI, KA_DSI, KA_Col, FA) were estimated from the adult bioequivalence data using the GITT step-function tablet-transit model. LOQ = 10 ng/mL (33.767 nmol/L); 58% of enteric-coated observations were below LOQ and were handled via the M3 likelihood-based method (paper Methods)."
  )

  ini({
    # ------------------------------------------------------------------
    # Region-specific first-order absorption rates from the intact
    # tablet (Table III). Table III labels KA in nmol/h but the units
    # are dimensionally 1/h (first-order rate constant); the label is
    # a paper artifact and is treated as h^-1 here.
    # ------------------------------------------------------------------
    lka_psi     <- log(1.06)  ; label("Log absorption rate from tablet in proximal SI KA_PSI (1/h)")               # Table III: KA_PSI = 1.06 1/h (RSE 21%)
    lka_dsi     <- log(8.64)  ; label("Log absorption rate from tablet in distal SI KA_DSI (1/h)")                 # Table III: KA_DSI = 8.64 1/h (RSE 258%; weak identifiability)
    lka_col     <- fixed(log(1))  ; label("Log absorption rate from tablet in colon KA_Col (1/h)")          # Table III: KA_Col = 1 1/h FIXED

    # Gut-wall bioavailability, logit-parameterized.
    logitfa     <- log(0.61 / (1 - 0.61))  ; label("Logit gut-wall bioavailability FA (unitless)")                 # Table III: FA = 0.61 (RSE 15%); logit(0.61) = 0.4479

    # ------------------------------------------------------------------
    # Disposition (Table III). Reference weight 70 kg per Table III
    # column headers 'L/h/70 kg^0.75' and 'L/70 kg'. Standard allometric
    # exponents (0.75 for clearances, 1 for volumes) are FIXED.
    # ------------------------------------------------------------------
    lcl         <- log(16.5)  ; label("Log clearance CL at 70 kg reference (L/h)")                                  # Table III: CL = 16.5 L/h/70kg^0.75 (RSE 5%)
    lvc         <- log(3.68)  ; label("Log central volume V1 at 70 kg reference (L)")                               # Table III: V1 = 3.68 L/70kg (RSE 17%)
    lq          <- log(1.75)  ; label("Log first inter-compartmental clearance Q2 at 70 kg reference (L/h)")        # Table III: Q2 = 1.75 L/h/70kg^0.75 (RSE 21%)
    lvp         <- log(7.48)  ; label("Log first peripheral volume V2 at 70 kg reference (L)")                      # Table III: V2 = 7.48 L/70kg (RSE 16%)
    lq2         <- log(7.21)  ; label("Log second inter-compartmental clearance Q3 at 70 kg reference (L/h)")       # Table III: Q3 = 7.21 L/h/70kg^0.75 (RSE 3%)
    lvp2        <- log(3.79)  ; label("Log second peripheral volume V3 at 70 kg reference (L)")                     # Table III: V3 = 3.79 L/70kg (RSE 7%)

    e_wt_cl     <- fixed(0.75); label("Allometric exponent on clearances")                                  # Table III headers 'L/h/70 kg^0.75' -- standard allometric scaling
    allo_v      <- fixed(1)   ; label("Allometric exponent on volumes")                                     # Table III headers 'L/70 kg' -- standard allometric scaling

    # ------------------------------------------------------------------
    # STEP function sigmoidicity (paper Equation 1; FIXED). Value not
    # reported in the paper; SIG = 20 is used to give a sharp transition
    # between GI regions (see vignette Errata).
    # ------------------------------------------------------------------
    sig         <- fixed(20)  ; label("STEP function sigmoidicity SIG (unitless;; not reported)")             # Paper Eq 1; SIG value not reported (see vignette Errata)

    # ------------------------------------------------------------------
    # Combined additive + proportional residual error (Table III).
    # Paper reports additive on the nmol/L scale (10 nmol/L, matching
    # LOQ 33.767 nmol/L = 10 ng/mL); converted to ng/mL via diclofenac
    # molecular weight 296.15 g/mol: 10 nmol/L * 296.15 / 1000 = 2.9615
    # ng/mL.
    # ------------------------------------------------------------------
    addSd       <- 2.9615     ; label("Additive residual SD (ng/mL; = 10 nmol/L via MW 296.15)")                    # Table III: additive = 10 nmol/L (RSE 0.3%); converted to 2.9615 ng/mL via diclofenac MW = 296.15 g/mol
    propSd      <- 0.126      ; label("Proportional residual SD (fraction)")                                        # Table III: proportional = 12.6% (RSE 21%)

    # ------------------------------------------------------------------
    # Inter-individual variability (Table III). Absorption IIVs are
    # very large (KA_PSI 109% CV, KA_DSI 53% CV) reflecting substantial
    # between-subject variability in regional absorption efficiency.
    # ------------------------------------------------------------------
    etalka_psi   ~ log(1.09^2 + 1)          # Table III: KA_PSI IIV = 109% CV (RSE 370%)
    etalka_dsi   ~ log(0.53^2 + 1)          # Table III: KA_DSI IIV = 53% CV (RSE 227%)
    etalogitfa   ~ log(0.08^2 + 1)          # Table III: FA IIV = 8% CV (RSE 13%); on logit scale approximation
    etalcl       ~ log(0.169^2 + 1)         # Table III: CL IIV = 16.9% (RSE 9%)
    etalvc       ~ log(0.201^2 + 1)         # Table III: V1 IIV = 20.1% (RSE 18%)
    etalq        ~ log(0.375^2 + 1)         # Table III: Q2 IIV = 37.5% (RSE 32%)
    etalvp       ~ log(0.327^2 + 1)         # Table III: V2 IIV = 32.7% (RSE 75%)
    etalq2       ~ log(0.531^2 + 1)         # Table III: Q3 IIV = 53.1% (RSE 37%)
    etalvp2      ~ log(0.354^2 + 1)         # Table III: V3 IIV = 35.4% (RSE 24%)
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters
    # ------------------------------------------------------------------
    ka_psi      <- exp(lka_psi + etalka_psi)
    ka_dsi      <- exp(lka_dsi + etalka_dsi)
    ka_col      <- exp(lka_col)

    logitfa_i   <- logitfa + etalogitfa
    fa          <- 1 / (1 + exp(-logitfa_i))

    # Allometric scaling to 70 kg reference.
    cl          <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc          <- exp(lvc + etalvc) * (WT / 70)^allo_v
    q           <- exp(lq  + etalq)  * (WT / 70)^e_wt_cl
    vp          <- exp(lvp + etalvp) * (WT / 70)^allo_v
    q2          <- exp(lq2 + etalq2) * (WT / 70)^e_wt_cl
    vp2         <- exp(lvp2 + etalvp2) * (WT / 70)^allo_v

    # ------------------------------------------------------------------
    # Tablet-position STEP indicators (Equation 1). No-return
    # subpopulation only. For the enteric-coated tablet the stomach
    # transit is a single event (IP_APSI); IP_FA is not used because
    # no drug is released in the stomach.
    # ------------------------------------------------------------------
    s_apsi      <- 1 / (1 + exp(-sig * (t - IP_APSI)))
    s_psidsi    <- 1 / (1 + exp(-sig * (t - IP_PSI_DSI)))
    s_dsic      <- 1 / (1 + exp(-sig * (t - IP_DSI_C)))

    pos_stomach <- 1 - s_apsi
    pos_psi     <- s_apsi - s_psidsi
    pos_dsi     <- s_psidsi - s_dsic
    pos_colon   <- s_dsic

    # ------------------------------------------------------------------
    # Effective absorption rate from the intact tablet: sum of region-
    # specific KA rates weighted by tablet-in-region indicator. The
    # tablet in the stomach releases nothing (enteric coating).
    # ------------------------------------------------------------------
    ka_eff      <- ka_psi * pos_psi + ka_dsi * pos_dsi + ka_col * pos_colon

    # ------------------------------------------------------------------
    # ODE system. Depot is the intact enteric-coated tablet drug
    # amount; only FA fraction of the absorbed amount reaches central.
    # The remaining (1 - FA) is treated as a bioavailability loss.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka_eff * depot
    d/dt(central)     <-  fa * ka_eff * depot - (cl / vc) * central -
                          (q  / vc) * central + (q  / vp) * peripheral1 -
                          (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <-  (q  / vc) * central - (q  / vp) * peripheral1
    d/dt(peripheral2) <-  (q2 / vc) * central - (q2 / vp2) * peripheral2

    # Observed plasma concentration. Central compartment holds mg
    # amount when dosing is in mg; central / vc gives mg/L; multiplied
    # by 1000 to yield ng/mL (matches paper LOQ 10 ng/mL = 33.767
    # nmol/L via MW 296.15 g/mol).
    Cc          <- central / vc * 1000

    Cc         ~ add(addSd) + prop(propSd)
  })
}

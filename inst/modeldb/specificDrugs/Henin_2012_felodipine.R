Henin_2012_felodipine <- function() {
  description <- paste(
    "Semi-mechanistic Gastro-Intestinal Transit Time (GITT) absorption",
    "model applied to extended-release felodipine (Henin 2012; the AAPS",
    "Journal 14:155-163). The tablet moves through five GI regions",
    "(fundus -> antrum -> proximal small intestine -> distal small",
    "intestine -> colon) governed by sigmoid step functions with",
    "per-subject inflection points IP_FA, IP_APSI, IP_PSI_DSI, IP_DSI_C",
    "supplied as covariates (population distributions from Table II).",
    "Zero-order drug release from the tablet into the current GI region",
    "at rates D1 (fundus 0.68 mg/h), D2/D3 (antrum + PSI 1.91 mg/h),",
    "D4/D5 (DSI + colon 1.16 mg/h); dissolved drug flows downstream",
    "via K23 (fundus -> antrum, 0.43 1/h) and K34 (antrum -> PSI,",
    "3.48 1/h fasted or 0.81 1/h fed) with both K23 and K34 increased",
    "by 5 1/h once the tablet reaches the small intestine (Table I",
    "footnote a). Dissolved drug in PSI, DSI, and colon is absorbed to",
    "a semi-physiological liver compartment (Vliver 0.0143 L/kg, QH 3.5",
    "L/h/kg^0.75) at rates K47 = K57 = 2.87 1/h and K67 = 1.15 1/h with",
    "gut-wall bioavailability FA (0.23 fasted or 0.39 fed). Hepatic",
    "extraction ratio EH = 0.50 and blood flow QH set the first-pass",
    "loss; systemic distribution is a 3-compartment model (Vcentral,",
    "Qperiph1/Vperiph1, Qperiph2/Vperiph2). This extraction encodes only",
    "the 'no return to fundus' subpopulation (7/12 felodipine subjects",
    "in the paper); the 3-component antrum-fundus-return mixture is",
    "documented in the vignette Assumptions and deviations."
  )
  reference <- paste(
    "Henin E, Bergstrand M, Standing JF, Karlsson MO. A mechanism-Based",
    "Approach for Absorption Modeling: The Gastro-Intestinal Transit",
    "Time (GITT) Model. The AAPS Journal. 2012;14(2):155-163.",
    "doi:10.1208/s12248-012-9324-y."
  )
  vignette <- "Henin_2012_gastrointestinal_transit_time_absorption"
  paper_specific_compartments <- c("fundus", "antrum", "proximal_si", "distal_si", "colon")

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "felodipine", units = "mg", specimen = "administration site", verified = FALSE),
    fundus      = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE),
    antrum      = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE),
    proximal_si = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE),
    distal_si   = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE),
    colon       = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE),
    liver       = list(analyte = "felodipine", units = "mg", specimen = "tissue", verified = FALSE),
    central     = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "felodipine", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight; used for allometric scaling of liver blood flow (QH = 3.5 * WT^0.75) and the fixed per-kg liver volume (Vliver = 0.0143 * WT).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Continuous covariate. The paper does not state a reference weight for the felodipine cohort (n = 6 healthy volunteers); the allometric parameters (QH and Vliver) are expressed per kg^0.75 and per kg respectively so any positive WT is admissible. See vignette Errata for the paper's population summary.",
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed-vs-fasted at-dosing indicator. 1 = fed condition (postprandial high-fat meal per Weitschies 2005 study design), 0 = fasted (overnight).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Per-dose indicator. Enters the model via three fed/fasted parameter switches (Table I): K34 (antrum -> PSI dissolved-drug transfer) = 3.48 h^-1 fasted / 0.81 h^-1 fed; FA (gut-wall bioavailability fraction) = 0.23 fasted / 0.39 fed; and additionally the tablet-transit inflection points IP_FA (fundus -> antrum) and IP_APSI (antrum -> PSI) are longer in the fed state per Table II. Set FED per observation record so the model interpolates correctly.",
      source_name        = "FED"
    ),
    IP_FA = list(
      description        = "Individual inflection-point time (h) at which the sigmoid step for tablet transition from fundus to antrum equals 0.5. Sampled from the fixed population distribution IP_FA = MRT_fundus * exp(eta_fundus) with MRT_fundus = 0.4 h (fasted) or 1.04 h (fed) and eta_fundus ~ N(0, VRT_fundus) with VRT_fundus = 0.46 (fasted) or 1.09 (fed) h^2 per Table II (values fixed from the upstream Bergstrand 2009 Markov-chain fit).",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject covariate. Population MRT and VRT are fasted/fed-dependent; the vignette shows an example of sampling IP_FA and its siblings from these distributions.",
      source_name        = "IP_FA"
    ),
    IP_APSI = list(
      description        = "Individual inflection-point time (h) for tablet transition from antrum to proximal small intestine (STEP value 0.5). Sampled as MRT_antrum * exp(eta_antrum) with MRT_antrum = 0.32 h (fasted) or 1.58 h (fed) and eta ~ N(0, VRT_antrum) with VRT_antrum = 0.15 (fasted) or 2.50 (fed) h^2 per Table II.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject covariate. Fed-vs-fasted stratification of the population distribution is required.",
      source_name        = "IP_APSI"
    ),
    IP_PSI_DSI = list(
      description        = "Individual inflection-point time (h) for tablet transition from proximal to distal small intestine (STEP value 0.5). Sampled as MRT_psi * exp(eta_psi) with MRT_psi = 1.17 h and eta ~ N(0, VRT_psi) with VRT_psi = 1.37 h^2 per Table II.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject covariate. Not fed/fasted dependent per Table II.",
      source_name        = "IP_PSI_DSI"
    ),
    IP_DSI_C = list(
      description        = "Individual inflection-point time (h) for tablet transition from distal small intestine to colon (STEP value 0.5). Sampled as MRT_dsi * exp(eta_dsi) with MRT_dsi = 1.22 h and eta ~ N(0, VRT_dsi) with VRT_dsi = 1.48 h^2 per Table II.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject covariate. Not fed/fasted dependent per Table II.",
      source_name        = "IP_DSI_C"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 6L,
    n_studies      = 1L,
    disease_state  = "Healthy volunteers",
    dose_range     = "Single 10 mg extended-release felodipine tablet, oral, under fasting or fed (postprandial high-fat meal) conditions in a crossover design.",
    regions        = "Charite, Berlin, Germany",
    notes          = "Cross-over MMM study by Weitschies et al. (2005; Ref 16). Six healthy volunteers received magnetically labelled 10 mg extended-release felodipine tablets under fasting and fed conditions with simultaneous MMM tracking of tablet GI position, drug release, and plasma sampling. Extended-release formulation; 11 mm tablet diameter. Felodipine parameters extracted from Table I; GI transit residence-time distribution (fixed as prior information) from Table II; both were originally reported by Bergstrand et al. 2009 (CPT 86:77-83)."
  )

  ini({
    # ------------------------------------------------------------------
    # Zero-order drug release rates from the tablet (Table I). Only
    # three unique typical values: D1 (fundus), D2 = D3 (antrum and
    # proximal small intestine share 1.91 mg/h), D4 = D5 (distal small
    # intestine and colon share 1.16 mg/h). IIV on D1 only (9.2% CV).
    # ------------------------------------------------------------------
    ld_fundus     <- log(0.68)  ; label("Log zero-order tablet release rate in fundus D1 (mg/h)")                       # Table I: D1 = 0.68 mg/h (RSE 5.1%)
    ld_antrum_psi <- log(1.91)  ; label("Log zero-order tablet release rate in antrum + proximal SI D2 = D3 (mg/h)")   # Table I: D2 = D3 = 1.91 mg/h (RSE 5.0%)
    ld_dsi_colon  <- log(1.16)  ; label("Log zero-order tablet release rate in distal SI + colon D4 = D5 (mg/h)")      # Table I: D4 = D5 = 1.16 mg/h (RSE 3.3%)

    # ------------------------------------------------------------------
    # Dissolved-substance distribution rate constants (Table I). K23
    # (fundus -> antrum) and K34 (antrum -> PSI). Both are increased by
    # 5 1/h once the tablet reaches the small intestine (Table I
    # footnote a). K34 has separate fasted / fed typical values.
    # ------------------------------------------------------------------
    lk23         <- log(0.43)  ; label("Log dissolved-drug transfer rate fundus -> antrum K23 (1/h)")                  # Table I: K23 = 0.43 1/h (RSE 9.3%)
    lk34_fasted  <- log(3.48)  ; label("Log dissolved-drug transfer rate antrum -> PSI K34 (fasted, 1/h)")             # Table I: K34 fasted = 3.48 1/h (RSE 7.4%)
    lk34_fed     <- log(0.81)  ; label("Log dissolved-drug transfer rate antrum -> PSI K34 (fed, 1/h)")                # Table I: K34 fed = 0.81 1/h (RSE 5.9%)
    k_si_bump    <- fixed(5)   ; label("Additive increase in K23 and K34 after tablet reaches SI (1/h)")        # Table I footnote a: rate increased by 5 h^-1 after tablet movement to small intestine

    # ------------------------------------------------------------------
    # Absorption first-order rate constants from PSI, DSI, colon to
    # liver (Table I). K47 = K57 = 2.87 1/h (PSI and DSI share), K67 =
    # 1.15 1/h (colon).
    # ------------------------------------------------------------------
    lk47_57      <- log(2.87)  ; label("Log absorption rate constant PSI -> liver (K47 = K57, 1/h)")                    # Table I: K47 = K57 = 2.87 1/h (RSE 8.4%)
    lk67         <- log(1.15)  ; label("Log absorption rate constant colon -> liver K67 (1/h)")                         # Table I: K67 = 1.15 1/h (RSE 18.8%)

    # ------------------------------------------------------------------
    # Gut-wall bioavailability fraction FA (Table I). Separate fasted /
    # fed typical values with IIV in each state. Parameterised on the
    # logit scale so FA_i stays in (0, 1) under the log-normal-like eta.
    # ------------------------------------------------------------------
    logitfa_fasted <- log(0.23 / (1 - 0.23))  ; label("Logit gut-wall bioavailability FA (fasted)")  # Table I: FA fasted = 0.23 (RSE 12.9%); logit(0.23) = -1.208
    logitfa_fed    <- log(0.39 / (1 - 0.39))  ; label("Logit gut-wall bioavailability FA (fed)")     # Table I: FA fed = 0.39 (RSE 13.8%); logit(0.39) = -0.447

    # ------------------------------------------------------------------
    # Hepatic first-pass parameters (Table I). Hepatic extraction ratio
    # EH is parameterised on the logit scale to keep individual EH_i in
    # (0, 1); the linear-scale typical is 0.50. Liver blood flow QH is
    # allometrically scaled per WT^0.75; liver volume Vliver is per WT.
    # Both allometric anchors are held FIXED per Table I.
    # ------------------------------------------------------------------
    logiteh      <- 0         ; label("Logit hepatic extraction ratio EH (logit(0.50) = 0)")                            # Table I: EH = 0.50 (RSE 3.5%); logit(0.5) = 0
    qh_per_kg075 <- fixed(3.5); label("Liver blood flow allometric constant QH (L/h/kg^0.75)")                   # Table I: QH = 3.5 L/h/kg^0.75 FIXED
    vliver_perkg <- fixed(0.0143); label("Liver volume per kg body weight Vliver (L/kg)")                         # Table I: Vliver = 0.0143 L/kg FIXED

    # ------------------------------------------------------------------
    # Systemic disposition (Table I). Central and two peripheral
    # compartments; log-transformed volumes and clearances. Table I
    # reports absolute values (not allometrically-scaled per kg^0.75)
    # for a 70-kg cohort; the model preserves this by keeping vc / vp /
    # q as WT-independent typical values. Downstream users who scale to
    # non-70-kg populations may want to add an allometric layer.
    # ------------------------------------------------------------------
    lvc          <- log(20.4) ; label("Log central volume of distribution Vcentral (L)")                                # Table I: Vcentral = 20.4 L (RSE 5.0%)
    lq           <- log(174)  ; label("Log inter-compartmental clearance Qperiph1 (L/h)")                               # Table I: Qperiph1 = 174 L/h (RSE 3.1%)
    lvp          <- log(88)   ; label("Log first peripheral volume Vperiph1 (L)")                                       # Table I: Vperiph1 = 88 L (RSE 3.0%)
    lq2          <- log(21.9) ; label("Log second inter-compartmental clearance Qperiph2 (L/h)")                        # Table I: Qperiph2 = 21.9 L/h (RSE 7.3%)
    lvp2         <- log(166)  ; label("Log second peripheral volume Vperiph2 (L)")                                      # Table I: Vperiph2 = 166 L (RSE 4.1%)

    # ------------------------------------------------------------------
    # Sigmoidicity factor of the tablet-position STEP functions
    # (Equation 1 of the paper). The paper describes SIG qualitatively
    # ('the higher the SIG, the steeper the step function') but does
    # not report a numeric typical value in the main text or tables and
    # the NONMEM control-file supplement is not on disk. SIG = 20 is a
    # commonly used value that yields a near-instantaneous transition
    # (the 10-90 percentile window spans ~ 0.22 h at SIG = 20). Held
    # FIXED to preserve the paper's mechanistic intent of a sharp
    # step-function switch between GI regions.
    # ------------------------------------------------------------------
    sig          <- fixed(20) ; label("STEP-function sigmoidicity factor SIG (unitless;; not reported in paper)") # Paper Eq 1; specific SIG value not reported in Henin 2012 (see vignette Errata) -- FIXED at 20 per common practice

    # ------------------------------------------------------------------
    # Residual error (Table I): 23% proportional on plasma felodipine
    # concentration.
    # ------------------------------------------------------------------
    propSd       <- 0.23      ; label("Proportional residual error (fraction; Cc scale)")                                # Table I: Residual error = 23% (RSE 8.9%)

    # ------------------------------------------------------------------
    # Inter-individual variability (Table I). Values converted from the
    # reported CV% to omega^2 via log(CV^2 + 1) (log-normal parameter
    # variance on the log-transformed structural scale). FA IIV is on
    # the logit scale via the reported linear-scale CV; the operational
    # variance carried on logitfa is documented alongside each entry.
    # ------------------------------------------------------------------
    etald_fundus  ~ log(0.092^2 + 1)         # Table I: D1 IIV = 9.2% (RSE 30.8%)
    etalk34       ~ log(0.464^2 + 1)         # Table I: K34 fasted IIV = 46.4% (RSE 32.1%); shared across fasted (fed IIV not reported)
    etalogitfa_fasted ~ log(0.174^2 + 1)     # Table I: FA fasted IIV = 17.4% (RSE 15.7%); carried on logit scale approximation
    etalogitfa_fed    ~ log(0.123^2 + 1)     # Table I: FA fed IIV = 12.3% (RSE 14.2%); carried on logit scale approximation
    etalogiteh    ~ log(0.108^2 + 1)         # Table I: EH IIV = 10.8% (RSE 17.5%); carried on logit scale approximation
    etalvc        ~ log(0.479^2 + 1)         # Table I: Vcentral IIV = 47.9% (RSE 30.8%)
    etalq         ~ log(0.234^2 + 1)         # Table I: Qperiph1 IIV = 23.4% (RSE 22.7%)
    etalvp        ~ log(0.250^2 + 1)         # Table I: Vperiph1 IIV = 25.0% (RSE 23.6%)
    etalvp2       ~ log(0.184^2 + 1)         # Table I: Vperiph2 IIV = 18.4% (RSE 21.5%)
  })

  model({
    # ------------------------------------------------------------------
    # Individual structural parameters
    # ------------------------------------------------------------------
    d_fundus       <- exp(ld_fundus + etald_fundus)
    d_antrum_psi   <- exp(ld_antrum_psi)
    d_dsi_colon    <- exp(ld_dsi_colon)

    k23_base       <- exp(lk23)
    k34_typ        <- exp(lk34_fasted) * (1 - FED) + exp(lk34_fed) * FED
    k34_base       <- k34_typ * exp(etalk34)

    k47            <- exp(lk47_57)
    k57            <- exp(lk47_57)
    k67            <- exp(lk67)

    # Gut-wall bioavailability on logit scale. Per-condition eta is
    # additively applied on the logit scale; the fed / fasted branch
    # selects the typical value and the matching eta variance.
    logitfa_fasted_i <- logitfa_fasted + etalogitfa_fasted
    logitfa_fed_i    <- logitfa_fed    + etalogitfa_fed
    logitfa_i        <- logitfa_fasted_i * (1 - FED) + logitfa_fed_i * FED
    fa               <- 1 / (1 + exp(-logitfa_i))

    # Hepatic extraction ratio via logit-scale IIV.
    logiteh_i      <- logiteh + etalogiteh
    eh             <- 1 / (1 + exp(-logiteh_i))

    # Semi-physiological liver: allometric flow, per-kg volume.
    qh             <- qh_per_kg075 * WT^0.75
    vliver_i       <- vliver_perkg * WT

    # Systemic disposition.
    vc             <- exp(lvc + etalvc)
    q              <- exp(lq + etalq)
    vp             <- exp(lvp + etalvp)
    q2             <- exp(lq2)
    vp2            <- exp(lvp2 + etalvp2)

    # ------------------------------------------------------------------
    # Tablet-position STEP indicators (Equation 1). Compact numerically
    # stable form: STEP(t, IP, SIG) = 1 / (1 + exp(-SIG * (t - IP))).
    # No-return-to-fundus subpopulation only: fundus -> antrum -> PSI
    # -> DSI -> colon in monotonic order (see description).
    # ------------------------------------------------------------------
    s_fa           <- 1 / (1 + exp(-sig * (t - IP_FA)))
    s_apsi         <- 1 / (1 + exp(-sig * (t - IP_APSI)))
    s_psidsi       <- 1 / (1 + exp(-sig * (t - IP_PSI_DSI)))
    s_dsic         <- 1 / (1 + exp(-sig * (t - IP_DSI_C)))

    pos_fundus     <- 1 - s_fa
    pos_antrum     <- s_fa - s_apsi
    pos_psi        <- s_apsi - s_psidsi
    pos_dsi        <- s_psidsi - s_dsic
    pos_colon      <- s_dsic

    # Tablet-reaches-SI indicator for the K23 / K34 bump (footnote a).
    in_si          <- pos_psi + pos_dsi + pos_colon
    k23_dyn        <- k23_base + k_si_bump * in_si
    k34_dyn        <- k34_base + k_si_bump * in_si

    # ------------------------------------------------------------------
    # Zero-order tablet release (only while drug remains in tablet).
    # Smooth non-negative gate keeps depot from crossing zero during
    # integration: depot_active ~ 1 while depot > 0.1 mg, ramps to 0
    # over a ~0.1 mg window. Preserves the paper's zero-order-until-
    # exhausted release intent without introducing an if / else branch.
    # ------------------------------------------------------------------
    depot_active   <- 1 / (1 + exp(-50 * (depot - 0.1)))
    release_fundus <- d_fundus     * pos_fundus * depot_active
    release_antrum <- d_antrum_psi * pos_antrum * depot_active
    release_psi    <- d_antrum_psi * pos_psi    * depot_active
    release_dsi    <- d_dsi_colon  * pos_dsi    * depot_active
    release_colon  <- d_dsi_colon  * pos_colon  * depot_active
    release_total  <- release_fundus + release_antrum + release_psi + release_dsi + release_colon

    # ------------------------------------------------------------------
    # Dissolved-drug and semi-physiological-liver mass balance.
    # Dissolved drug in fundus / antrum only travels downstream via
    # K23_dyn / K34_dyn; dissolved drug in PSI / DSI / colon leaves the
    # gut via K47 / K57 / K67 with only FA fraction reaching the liver
    # (the remaining (1 - FA) fraction is treated as a gut-wall
    # bioavailability loss per the paper's Figure 1 caption). The liver
    # follows the well-stirred model with blood flow qh; central losses
    # to liver at qh represent hepatic recirculation and (1 - eh) *
    # qh * liver / vliver_i returns to central, with eh * qh * liver /
    # vliver_i eliminated hepatically.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -release_total
    d/dt(fundus)      <-  release_fundus - k23_dyn * fundus
    d/dt(antrum)      <-  release_antrum + k23_dyn * fundus - k34_dyn * antrum
    d/dt(proximal_si) <-  release_psi    + k34_dyn * antrum - k47 * proximal_si
    d/dt(distal_si)   <-  release_dsi                        - k57 * distal_si
    d/dt(colon)       <-  release_colon                      - k67 * colon
    d/dt(liver)       <-  fa * (k47 * proximal_si + k57 * distal_si + k67 * colon) +
                          qh * central / vc - qh * liver / vliver_i
    d/dt(central)     <-  qh * (1 - eh) * liver / vliver_i - qh * central / vc -
                          q  * central / vc + q  * peripheral1 / vp -
                          q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1) <-  q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2

    # Observed plasma concentration. Central compartment holds mg
    # amount when dosing is in mg; central / vc gives mg/L, multiplied
    # by 1000 to yield ng/mL. Felodipine molecular weight is 384.25
    # g/mol for users who want to convert to nmol/L (1 ng/mL felodipine
    # = 1000 / 384.25 = 2.602 nmol/L).
    Cc             <- central / vc * 1000

    Cc            ~ prop(propSd)
  })
}

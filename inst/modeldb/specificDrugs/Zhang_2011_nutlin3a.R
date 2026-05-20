Zhang_2011_nutlin3a <- function() {
  description <- "Preclinical (mouse). Whole-body PBPK model for nutlin-3a (MDM2 inhibitor) in adult C57BL/6 mice after intravenous and oral administration (Zhang et al. 2011, DMD). Thirteen physiological tissue compartments (adipose, adrenal gland, bone marrow, brain, intestine + lumen, liver, lung, muscle, retina, spleen, vitreous, residual diffusion-limited tissue + residual vascular space) with arterial and venous blood pools (75/25 split of total blood volume). Perfusion-limited tissues use a partition coefficient K_i; the eye is modelled as retina + vitreous coupled by a permeability-surface-area product PA_VIT; the residual compartment is diffusion-limited (5% vascular space, 95% tissue, coupled by PA_RES). Elimination is combined linear (hepatic, k_e) and saturable Michaelis-Menten (arterial, V_max/K_m). Oral absorption is first-order from an intestinal lumen depot. Plasma protein binding is reported (B_max = 286 uM, K_A = 0.085 1/uM, Langmuir form) but is only used for an unbound-concentration derivation that the paper applies to tissue exposure / IC50 comparisons, not to the elimination ODEs; the ODE system operates on total concentrations. The model is intended for typical-value simulation."
  reference <- "Zhang F, Tagen M, Throm S, Mallari J, Miller L, Guy RK, Dyer MA, Williams RT, Roussel MF, Nemeth K, Zhu F, Zhang J, Lu M, Panetta JC, Boulos N, Stewart CF. Whole-body physiologically based pharmacokinetic model for nutlin-3a in mice after intravenous and oral administration. Drug Metab Dispos. 2011;39(1):15-21. doi:10.1124/dmd.110.035915."
  vignette <- "Zhang_2011_nutlin3a"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    weight = "kg"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Mouse body weight. Used to scale all tissue volumes linearly from the %-of-body-weight values in Table 1 (assumed tissue density 1 kg/L). The paper's adult C57BL/6 mice were ~25 g (0.025 kg); the model file defaults to WT = 0.025 if users do not supply it via the event table or rxSolve(params = c(WT = ...)). The Brown et al. 1997 physiological tissue blood flows tabulated in mL/h were converted to L/h.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "mouse (adult C57BL/6, male and female)",
    n_subjects     = 355L,
    n_studies      = 2L,
    age_range      = "adult",
    weight_range   = "~25 g (typical adult C57BL/6)",
    sex_female_pct = 4.8,
    disease_state  = "Healthy, non-tumour-bearing adult mice. The PBPK model was fit to pooled plasma + tissue PK from non-tumour-bearing animals.",
    dose_range     = "Intravenous: 10 and 20 mg/kg single bolus (tail vein). Oral: 50, 100, and 200 mg/kg single gavage. Model simulations explored 50, 100, 200, and 400 mg/kg QD and BID.",
    regions        = "St. Jude Children's Research Hospital (Memphis, TN, USA)",
    notes          = "Study 1: 145 mice (128 male + 17 female) split into oral 100 mg/kg, oral 200 mg/kg, and IV 10 mg/kg arms with destructive sampling for tissue distribution. Study 2: 210 adult male mice in oral 50, 100 mg/kg and IV 10, 20 mg/kg arms with serial plasma sampling. Tissues sampled: plasma, liver, spleen, intestine, muscle, lung, adipose, bone marrow, adrenal gland, brain, retina, vitreous (the union across the two studies)."
  )

  ini({
    # Estimated structural PBPK parameters (Table 2). Units are not stated in
    # the paper; the values are encoded in mass-based units (mg, mg/L, L/h)
    # consistent with the mg/kg dosing regimen and the Km / Vmax magnitudes:
    # at C >> Km the saturable arm reaches V_max ~ 0.0287 mg/h, which clears a
    # ~0.5 mg bolus in ~17 h; below ~10 uM (5.8 mg/L) the apparent clearance
    # V_max / K_m ~ 0.57 L/h gives a half-life of ~1.5 min, matching the
    # paper's description of "very rapid" clearance at low concentration. See
    # the vignette Assumptions and deviations section for the unit decision.
    lka      <- log(0.409)   ; label("First-order oral absorption rate (1/h)")             # Table 2 (ka = 0.409)
    lke      <- log(0.0160)  ; label("Hepatic linear elimination clearance ke (L/h), acts on arterial concentration per liver ODE") # Table 2 (Ke = 0.0160)
    lvmax    <- log(0.0287)  ; label("Saturable elimination Vmax in arterial blood (mg/h)") # Table 2 (Vmax = 0.0287)
    lkm      <- log(0.050)   ; label("Michaelis constant for arterial saturable elimination (mg/L)") # Table 2 (Km = 0.050)

    # Tissue:plasma partition coefficients (Table 2). Unitless ratios. Held as
    # estimated point values in the paper; encoded here as log-transformed
    # constants so users can perturb via rxSolve(params=) if desired.
    lk_adi   <- log(1.61)    ; label("Adipose tissue:plasma partition coefficient (unitless)")  # Table 2 (K_ADI)
    lk_adr   <- log(2.05)    ; label("Adrenal gland tissue:plasma partition coefficient (unitless)") # Table 2 (K_ADR)
    lk_bra   <- log(0.055)   ; label("Brain tissue:plasma partition coefficient (unitless)")    # Table 2 (K_BRA)
    lk_int   <- log(12.2)    ; label("Intestine tissue:plasma partition coefficient (unitless)") # Table 2 (K_INT)
    lk_liv   <- log(7.54)    ; label("Liver tissue:plasma partition coefficient (unitless)")    # Table 2 (K_LIV)
    lk_lun   <- log(1.78)    ; label("Lung tissue:plasma partition coefficient (unitless)")     # Table 2 (K_LUN)
    lk_mus   <- log(2.08)    ; label("Muscle tissue:plasma partition coefficient (unitless)")   # Table 2 (K_MUS)
    lk_ret   <- log(4.01)    ; label("Retina tissue:plasma partition coefficient (unitless)")   # Table 2 (K_RET)
    lk_spl   <- log(2.72)    ; label("Spleen tissue:plasma partition coefficient (unitless)")   # Table 2 (K_SPL)
    lk_vit   <- log(0.012)   ; label("Vitreous tissue:plasma partition coefficient (unitless)") # Table 2 (K_VIT)
    lk_mrw   <- fixed(log(1.0)) ; label("Bone marrow tissue:plasma partition coefficient (unitless); not reported in Table 2, fixed at 1.0") # see vignette Errata
    lk_res   <- log(4.8)     ; label("Residual diffusion-limited tissue:plasma partition coefficient (unitless)") # Table 2 (K_RES)

    # Permeability-surface area products for diffusion-limited compartments
    # (Table 2). Encoded in L/h to match the L/h flow units used in model().
    lpa_vit  <- log(0.0036)  ; label("Vitreous permeability-surface area product PA_VIT (L/h)")  # Table 2 (PA_VIT)
    lpa_res  <- log(0.0048)  ; label("Residual permeability-surface area product PA_RES (L/h)")  # Table 2 (PA_RES)

    # In-vitro derived (Results, p.18): blood/plasma concentration ratio. Used
    # to map venous-blood concentration to a plasma concentration for the Cc
    # observation. Plasma protein binding (Bmax, KA) is not in the ODE system;
    # see model() for the unbound-fraction derivation it powers.
    lbp      <- fixed(log(0.70)); label("Blood:plasma concentration ratio (unitless)")  # Results, Fig 2A (B/P = 0.70)

    # In-vitro Langmuir plasma protein binding (Results, p.18). Bmax in uM,
    # KA in 1/uM. Used to compute the unbound plasma concentration C_pu and
    # unbound fraction fub reported alongside total plasma C. NOT used in the
    # elimination ODEs.
    bmax     <- fixed(286)   ; label("Plasma protein binding sites Bmax (uM)")                 # Results (Bmax = 286)
    ka_bind  <- fixed(0.085) ; label("Plasma protein binding association constant KA (1/uM)")  # Results (KA = 0.085)
    mw       <- fixed(581.5) ; label("Nutlin-3a molecular weight (g/mol)")                     # PubChem CID 11433190 (581.5 g/mol)

    # IIV -- Table 2 reports CV% values. omega^2 on the log scale is
    # log(CV^2 + 1) for log-normal random effects.
    etalka   ~ 0.0926   # log(0.312^2 + 1); Table 2 (IIV ka = 31.2%)
    etalke   ~ 0.00408  # log(0.064^2 + 1); Table 2 (IIV ke =  6.4%)
    etalvmax ~ 0.153    # log(0.406^2 + 1); Table 2 (IIV Vmax = 40.6%)

    # Residual error (Table 2): 35.6% -- encoded as proportional SD on the
    # linear scale. The paper does not state whether this is additive on
    # log-scale (= proportional on linear scale) or strictly proportional;
    # 35.6% is consistent with either interpretation for a fitted model.
    propSd   <- 0.356   ; label("Proportional residual error (fraction)")  # Table 2 (Residual error = 35.6%)
  })

  model({
    # ------------------------------------------------------------------
    # Body weight (kg) and assumed tissue density (1 kg/L). All tissue
    # volumes in Table 1 are reported as % body weight (mass); we convert
    # to L assuming unit density. The Q_B column of Table 1 is in mL/h and
    # is converted to L/h here.
    # ------------------------------------------------------------------
    ka      <- exp(lka   + etalka)
    ke_cl   <- exp(lke   + etalke)
    vmax    <- exp(lvmax + etalvmax)
    km      <- exp(lkm)

    k_adi   <- exp(lk_adi)
    k_adr   <- exp(lk_adr)
    k_bra   <- exp(lk_bra)
    k_int   <- exp(lk_int)
    k_liv   <- exp(lk_liv)
    k_lun   <- exp(lk_lun)
    k_mus   <- exp(lk_mus)
    k_ret   <- exp(lk_ret)
    k_spl   <- exp(lk_spl)
    k_vit   <- exp(lk_vit)
    k_mrw   <- exp(lk_mrw)
    k_res   <- exp(lk_res)

    pa_vit  <- exp(lpa_vit)
    pa_res  <- exp(lpa_res)
    bp      <- exp(lbp)

    # Body weight in kg (mouse default 0.025 kg = 25 g). Users override
    # via the event-table covariate column WT or rxSolve(params=c(WT=...)).
    bw <- WT

    # Tissue volumes (L) from Table 1 fractions x bw (kg) x 1 L/kg.
    v_blo <- 0.049   * bw
    v_adi <- 0.068   * bw
    v_adr <- 0.00048 * bw
    v_mrw <- 0.058   * bw
    v_bra <- 0.0165  * bw
    v_int <- 0.0362  * bw
    v_liv <- 0.0549  * bw
    v_lun <- 0.0073  * bw
    v_mus <- 0.384   * bw
    v_ret <- 0.0004  * bw
    v_spl <- 0.0035  * bw
    v_vit <- 0.00035 * bw
    v_res_total <- 0.299 * bw

    # Venous / arterial blood pools -- text: "the volume of the venous pool
    # was fixed to 75% of the total blood volume."
    v_ven <- 0.75 * v_blo
    v_art <- 0.25 * v_blo

    # Residual diffusion-limited compartment -- text: "the vascular space
    # assumed to be 5% of the residual volume."
    v_res_vasc <- 0.05 * v_res_total
    v_res_tis  <- 0.95 * v_res_total

    # Blood flows (L/h) from Table 1 mL/h x 1e-3.
    q_blo <- 0.839
    q_adi <- 0.0587
    q_adr <- 0.00252
    q_mrw <- 0.0923
    q_bra <- 0.0277
    q_int <- 0.109
    q_ha  <- 0.0168     # hepatic artery (Table 1 row "Liver" Q_B = 16.8 mL/h)
    q_lun <- 0.839      # = cardiac output
    q_mus <- 0.133
    q_ret <- 0.00316
    q_spl <- 0.00948
    q_res <- 0.2569

    # Liver total blood flow Q_LIV = hepatic artery + splanchnic returns
    # (spleen + intestine), per text "Liver blood flow (QLIV) was the sum
    # of the blood flow to the hepatic artery, spleen, and liver".
    q_liv <- q_ha + q_spl + q_int

    # Concentrations (mg/L) in each tissue from amount (mg) / volume (L).
    c_ven <- venous   / v_ven
    c_art <- arterial / v_art
    c_lun <- lung     / v_lun
    c_adi <- adipose  / v_adi
    c_adr <- adrenal  / v_adr
    c_mrw <- bonemarrow / v_mrw
    c_bra <- brain    / v_bra
    c_int <- intestine / v_int
    c_liv <- liver    / v_liv
    c_mus <- muscle   / v_mus
    c_ret <- retina   / v_ret
    c_spl <- spleen   / v_spl
    c_vit <- vitreous / v_vit
    c_res_vasc <- res_vasc / v_res_vasc
    c_res_tis  <- res_tis  / v_res_tis

    # ------------------------------------------------------------------
    # ODEs
    # depot -- intestinal lumen for oral dosing. First-order absorption
    # into the intestine tissue compartment ("Absorption from the lumen
    # was assumed to be linear based on an absorption rate constant ka").
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot

    # Venous blood pool -- sum of venous returns from every tissue (except
    # lung, which receives from venous), plus the spleen and intestine,
    # whose returns feed the liver and then exit via the liver. The
    # spleen and intestine therefore appear in the liver inflow, not in
    # the venous pool inflow.
    d/dt(venous) <-
        q_adi * (c_adi / k_adi) +
        q_adr * (c_adr / k_adr) +
        q_mrw * (c_mrw / k_mrw) +
        q_bra * (c_bra / k_bra) +
        q_mus * (c_mus / k_mus) +
        q_ret * (c_ret / k_ret) +
        q_liv * (c_liv / k_liv) +
        q_res * c_res_vasc -
        q_blo * c_ven

    # Lung -- receives all cardiac output from the venous pool; output to
    # arterial. "The lungs received all venous input, and the arterial
    # input was the output from the lungs."
    d/dt(lung) <- q_lun * (c_ven - c_lun / k_lun)

    # Arterial blood pool -- receives lung output; contains saturable
    # elimination. "V_ART x dA_ART/dt = Q_BLO x (C_LUN/K_LUN - C_ART)
    # - V_max x C_ART / (K_m + C_ART)".
    d/dt(arterial) <- q_blo * (c_lun / k_lun - c_art) -
                      vmax * c_art / (km + c_art)

    # Perfusion-limited tissues. dA_i/dt = Q_i x (C_ART - C_i/K_i).
    d/dt(adipose)    <- q_adi * (c_art - c_adi / k_adi)
    d/dt(adrenal)    <- q_adr * (c_art - c_adr / k_adr)
    d/dt(bonemarrow) <- q_mrw * (c_art - c_mrw / k_mrw)
    d/dt(brain)      <- q_bra * (c_art - c_bra / k_bra)
    d/dt(muscle)     <- q_mus * (c_art - c_mus / k_mus)
    d/dt(spleen)     <- q_spl * (c_art - c_spl / k_spl)

    # Intestine tissue -- perfusion-limited inflow from arterial, outflow
    # to liver via portal vein, plus first-order absorption input from
    # the depot (lumen): ka x A_depot, per "Absorption from the lumen
    # was assumed to be linear..." and the Figure 1 schematic.
    d/dt(intestine) <- ka * depot + q_int * (c_art - c_int / k_int)

    # Liver -- hepatic artery inflow + portal returns from spleen and
    # intestine, outflow proportional to combined liver flow, plus a
    # linear elimination term. The paper writes ke x C_ART for the
    # elimination term; the literal form is preserved here. See vignette
    # Assumptions and deviations.
    d/dt(liver) <-
        q_ha  * c_art +
        q_spl * (c_spl / k_spl) +
        q_int * (c_int / k_int) -
        q_liv * (c_liv / k_liv) -
        ke_cl * c_art

    # Retina -- perfusion-limited inflow from arterial, plus diffusion to
    # vitreous (no direct blood flow to vitreous). "Input into the
    # vitreous was by diffusion from the retina."
    d/dt(retina)   <- q_ret * (c_art - c_ret / k_ret) -
                      pa_vit * (c_ret - c_vit / k_vit)
    d/dt(vitreous) <- pa_vit * (c_ret - c_vit / k_vit)

    # Residual diffusion-limited -- vascular subspace coupled to tissue
    # subspace via permeability product. "Modeling this compartment as
    # perfusion-limited did not adequately describe the multiexponential
    # profile of nutlin-3a. Therefore, the residual compartment was
    # modeled as diffusion-limited..."
    d/dt(res_vasc) <- q_res * (c_art - c_res_vasc) -
                      pa_res * (c_res_vasc - c_res_tis / k_res)
    d/dt(res_tis)  <- pa_res * (c_res_vasc - c_res_tis / k_res)

    # ------------------------------------------------------------------
    # Plasma protein binding (Results, p.18). Used to derive an unbound
    # fraction (fub) and unbound plasma concentration from the total
    # plasma concentration for downstream IC50 comparisons. The
    # elimination ODEs above are on total concentrations; fub is
    # reported as an auxiliary output. Both bmax and ka_bind are in uM
    # units; convert plasma concentration mg/L -> uM via mw (g/mol).
    # ------------------------------------------------------------------
    Cc_uM   <- (c_ven / bp) * 1000 / mw
    cpu_uM  <- (-(1 + bmax - ka_bind * Cc_uM) +
                  sqrt((1 + bmax - ka_bind * Cc_uM) ^ 2 + 4 * ka_bind * Cc_uM)) /
                 (2 * ka_bind)
    fub     <- cpu_uM / (Cc_uM + 1e-12)
    Cc_unbound <- (c_ven / bp) * fub

    # ------------------------------------------------------------------
    # Observation -- plasma concentration in mg/L. C_plasma = C_venous_blood
    # / (B/P). The proportional residual error from Table 2 is applied
    # here.
    # ------------------------------------------------------------------
    Cc <- c_ven / bp
    Cc ~ prop(propSd)
  })
}

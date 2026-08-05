Ai_2024_ractopamine_goat_pbpk <- function() {
  description <- "PBPK (whole-body, hybrid flow- and membrane-limited; goat). Eleven-tissue physiologically based pharmacokinetic model for the beta-adrenoceptor agonist ractopamine in Liaoning cashmere goats after repeated oral gavage, comprising gastric contents, intestinal contents, liver, spleen, kidney, heart, lung, muscle, fat, brain, a lumped rest-of-body compartment, arterial plasma, venous plasma and a urinary-excretion sink. Liver, spleen, kidney, heart and lung are perfusion (flow) limited; muscle, fat, brain and the rest of the body are membrane (permeability) limited and each carry a vascular plasma sub-compartment plus a tissue sub-compartment separated by a permeability-area product. Absorption is first-order gastric emptying (Kst) into the gut lumen followed by first-order uptake (Ka) into the liver in competition with first-order fecal loss of unabsorbed drug (Kgut); elimination is hepatic (Clhe) plus renal (Clre), both acting on the unbound fraction. Built to predict edible-tissue residues and withdrawal times against Codex Alimentarius maximum residue limits (Ai 2024)."
  reference   <- "Ai J, Gao Y, Yang F, Zhao Z, Dong J, Wang J, Fu S, Ma Y, Gu X. Development and application of a physiologically-based pharmacokinetic model for ractopamine in goats. Front Vet Sci. 2024;11:1399043. doi:10.3389/fvets.2024.1399043"
  vignette    <- "Ai_2024_ractopamine_goat_pbpk"
  units       <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Every state below integrates a drug AMOUNT in ug,
  # exactly as the supplementary acslXtreme DERIVATIVE block does; the
  # per-tissue concentrations are formed algebraically in model().
  compartmentData <- list(
    # The stomach and gut states hold luminal drug (gavage contents and
    # intestinal contents), not a sampled biological matrix; "administration
    # site" is the closest member of the specimen vocabulary.
    stomach     = list(analyte = "ractopamine", units = "ug", specimen = "administration site", verified = TRUE),
    a_gut       = list(analyte = "ractopamine", units = "ug", specimen = "administration site", verified = TRUE),
    liver       = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    spleen      = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    kidney      = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    heart       = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    lung        = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    vp_muscle   = list(analyte = "ractopamine", units = "ug", specimen = "plasma", verified = TRUE),
    int_muscle  = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    vp_adipose  = list(analyte = "ractopamine", units = "ug", specimen = "plasma", verified = TRUE),
    int_adipose = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    vp_brain    = list(analyte = "ractopamine", units = "ug", specimen = "plasma", verified = TRUE),
    int_brain   = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    vp_other    = list(analyte = "ractopamine", units = "ug", specimen = "plasma", verified = TRUE),
    int_other   = list(analyte = "ractopamine", units = "ug", specimen = "tissue", verified = TRUE),
    arterial    = list(analyte = "ractopamine", units = "ug", specimen = "plasma", verified = TRUE),
    venous      = list(analyte = "ractopamine", units = "ug", specimen = "plasma", verified = TRUE),
    urine       = list(analyte = "ractopamine", units = "ug", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Scales every compartment volume (Vxx = Vcxx * WT), the cardiac output allometrically (Qtot = QAR * WT^0.74 * (1 - pcv)), and converts both clearances from L/h/kg to L/h. The supplementary acslX code sets `constant bw = 30`, the mean weight of the 27 male Liaoning cashmere goats (30 +/- 5 kg; Ai 2024 Section 2.2). Body weight was held constant over the 28-day dosing period and the subsequent depletion phase.",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "goat (Liaoning cashmere goat, Capra hircus; male)",
    n_subjects     = 27L,
    n_studies      = 1L,
    age_range      = "10 months",
    weight_range   = "30 +/- 5 kg",
    sex_female_pct = 0,
    disease_state  = "Healthy; one-week acclimation on a drug-free diet before dosing",
    dose_range     = "Single oral gavage 1 mg/kg BW and single intravenous 1 mg/kg BW (n = 6, crossover with a 15-day washout) for the pharmacokinetic study; continuous oral gavage 1 mg/kg BW per day for 28 days for the residue-depletion study (n = 21, including 3 controls)",
    regions        = "China (Institute of Feed Research, Chinese Academy of Agricultural Sciences, Beijing)",
    n_observations = NA_integer_,
    notes          = "Three of the 27 goats were sacrificed to measure the organ-weight fractions of Ai 2024 Table 3. Tissue blood flows are literature values for sheep and goats (Ai 2024 refs 17 and 28). The residue-depletion concentration-time data come from two earlier publications by the same laboratory (Ai 2024 refs 24 and 25, Zhao 2017 and Zhao 2019). Partition coefficients, both clearances, Ka, Kst and Kgut were optimised in the acslXtreme OptStatModule against the measured data (Ai 2024 Table 4). No individual-level fitting was performed, so the model carries no inter-individual random effects; Ai 2024 propagated parameter uncertainty instead through a 500-iteration Monte Carlo analysis over the sensitive parameters to derive the withdrawal times of Ai 2024 Table 7."
  )

  ini({
    # ---------------------------------------------------------------------
    # Tissue:plasma partition coefficients (Pxx, unitless), log-transformed
    # following the `lk_<organ>` convention of the sibling veterinary
    # residue-PBPK model Yang_2023_diclazuril_chicken_pbpk.R. All were
    # optimised in the acslXtreme OptStatModule and are FIXED at the
    # published final values. Ai 2024 Table 4 prints these rounded to four
    # decimals; the full-precision values below are the `constant` values
    # of the supplementary acslX INITIAL block, which is the code that
    # generated the published predictions.
    # ---------------------------------------------------------------------
    lk_liv <- fixed(log(2.558408));  label("Liver:plasma partition coefficient (Pli, unitless)")        # Ai 2024 Table 4 (2.5584) / supplement `Pli=2.558408`
    lk_kid <- fixed(log(1.77341));   label("Kidney:plasma partition coefficient (Pki, unitless)")       # Ai 2024 Table 4 (1.7734) / supplement `Pki=1.77341`
    lk_lun <- fixed(log(1.533326));  label("Lung:plasma partition coefficient (Plu, unitless)")         # Ai 2024 Table 4 (1.5333) / supplement `Plu=1.533326`
    lk_spl <- fixed(log(1.004675));  label("Spleen:plasma partition coefficient (Psp, unitless)")       # Ai 2024 Table 4 (1.0047) / supplement `Psp=1.004675`
    lk_hrt <- fixed(log(1.383948));  label("Heart:plasma partition coefficient (Phe, unitless)")        # Ai 2024 Table 4 (1.3839) / supplement `Phe=1.383948`
    lk_mus <- fixed(log(1.068578));  label("Muscle:plasma partition coefficient (Pmu, unitless)")       # Ai 2024 Table 4 (1.0686) / supplement `Pmu=1.068578`
    lk_adi <- fixed(log(0.7526104)); label("Fat:plasma partition coefficient (Pfa, unitless)")          # Ai 2024 Table 4 (0.7526) / supplement `Pfa=0.7526104`
    lk_bra <- fixed(log(0.6896452)); label("Brain:plasma partition coefficient (Pbr, unitless)")        # Ai 2024 Table 4 (0.6896) / supplement `Pbr=0.6896452`
    lk_res <- fixed(log(9.088835));  label("Rest-of-body:plasma partition coefficient (Pre, unitless)") # Ai 2024 Table 4 (9.0888) / supplement `Pre=9.088835`

    # ---------------------------------------------------------------------
    # Membrane permeability coefficients (PPxx, unitless) for the four
    # membrane-limited tissues. Each is multiplied by that tissue's plasma
    # flow to give the permeability-area product PAxx in L/h (supplement
    # `PAmu=PPmu*Qmu` etc.), so PPxx is the permeability expressed as a
    # fraction of tissue plasma flow. Optimised and FIXED.
    # ---------------------------------------------------------------------
    lpp_mus <- fixed(log(0.02712495));  label("Muscle membrane permeability coefficient (PPmu, fraction of muscle plasma flow)") # Ai 2024 Table 4 Ppmu (0.0271) / supplement `PPmu=0.02712495`
    lpp_adi <- fixed(log(0.00537883));  label("Fat membrane permeability coefficient (PPfa, fraction of fat plasma flow)")       # Ai 2024 Table 4 Ppfa (0.0054) / supplement `PPfa=0.00537883`
    lpp_bra <- fixed(log(0.006772557)); label("Brain membrane permeability coefficient (PPbr, fraction of brain plasma flow)")   # Ai 2024 Table 4 Ppbr (0.0068) / supplement `PPbr=0.006772557`
    lpp_res <- fixed(log(0.002063587)); label("Rest-of-body membrane permeability coefficient (PPre, fraction of rest plasma flow)") # Ai 2024 Table 4 Ppre (0.0021) / supplement `PPre=0.002063587`

    # ---------------------------------------------------------------------
    # Absorption. Gastric emptying (Kst) moves drug from the gastric
    # contents to the intestinal contents; from there uptake into the liver
    # competes with fecal loss of the unabsorbed fraction (Kgut). Ai 2024
    # Section 2.5 gives the bioavailability as F = Ka / (Ka + Kgut) and the
    # supplement computes exactly that (`F=ka/(ka+kint)`) -- the paper's
    # Table 4 prints this rate constant as "K int" while Section 2.5 calls
    # it Kgut; they are the same quantity (see the vignette Errata).
    # ---------------------------------------------------------------------
    lkst  <- fixed(log(0.09103672)); label("First-order gastric emptying rate constant (Kst, 1/h)")                        # Ai 2024 Table 4 (0.0910) / supplement `kst=0.09103672`
    lka   <- fixed(log(0.9861192));  label("First-order absorption rate constant from intestinal contents (Ka, 1/h)")      # Ai 2024 Table 4 (0.9861) / supplement `ka=0.9861192`
    lkgut <- fixed(log(0.9016421));  label("First-order fecal excretion rate constant for unabsorbed drug (Kgut, 1/h)")    # Ai 2024 Table 4 "K int" (0.9016) / supplement `kint=0.9016421`

    # ---------------------------------------------------------------------
    # Elimination and plasma protein binding. Both clearances act on the
    # unbound concentration in the eliminating organ (supplement
    # `reli=CCLhe*Cli/Pli*pfree`, `reki=CCLre*Cki/Pki*pfree`). Ai 2024
    # Table 4 reports Clre at its initial value with a standard deviation
    # of exactly 0, i.e. the optimiser never moved it -- see the vignette
    # Errata for what that implies about the renal route.
    # ---------------------------------------------------------------------
    lclhe  <- fixed(log(0.06243522));   label("Hepatic clearance per unit body weight (Clhe, L/h/kg)")      # Ai 2024 Table 4 (0.0624) / supplement `clhe=0.06243522`
    lclre  <- fixed(log(0.0001090409)); label("Renal clearance per unit body weight (Clre, L/h/kg)")        # Ai 2024 Table 4 (0.0001) / supplement `clre=0.0001090409`
    lpbind <- fixed(log(0.1651634));    label("Plasma protein binding ratio (pbind, fraction bound)")       # supplement `pbind=0.1651634` (not tabulated in the paper; appears in Ai 2024 Table 6 as a sensitivity parameter)

    # ---------------------------------------------------------------------
    # Residual error. Ai 2024 is a deterministic PBPK model fitted to
    # pooled tissue-residue means; it reports no residual-error model
    # (validation used linear regression of observed on predicted plus
    # residual plots, Ai 2024 Figures 6 and 7). The placeholders below
    # exist only so the model is a syntactically complete nlmixr2 object
    # for forward simulation; they are NOT paper-derived. Same convention
    # as the sibling Yang_2023_diclazuril_chicken_pbpk.R.
    # ---------------------------------------------------------------------
    propSd          <- fixed(0.10); label("Proportional residual error placeholder, venous plasma (fraction)") # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Cliver   <- fixed(0.10); label("Proportional residual error placeholder, liver (fraction)")         # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Ckidney  <- fixed(0.10); label("Proportional residual error placeholder, kidney (fraction)")        # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Cmuscle  <- fixed(0.10); label("Proportional residual error placeholder, muscle (fraction)")        # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Cadipose <- fixed(0.10); label("Proportional residual error placeholder, fat (fraction)")           # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Clung    <- fixed(0.10); label("Proportional residual error placeholder, lung (fraction)")          # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Cspleen  <- fixed(0.10); label("Proportional residual error placeholder, spleen (fraction)")        # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Cheart   <- fixed(0.10); label("Proportional residual error placeholder, heart (fraction)")         # not reported in Ai 2024; placeholder for syntactic completeness only
    propSd_Cbrain   <- fixed(0.10); label("Proportional residual error placeholder, brain (fraction)")         # not reported in Ai 2024; placeholder for syntactic completeness only
  })

  model({
    # ==================================================================
    # Goat physiology. Organ-weight fractions were measured by Ai 2024 in
    # three sacrificed goats (Table 3); tissue blood flows are literature
    # values for sheep and goats (Table 2, refs 17 and 28). These are
    # physiology rather than fitted quantities, so they are carried as
    # traceable literals here rather than as ini() parameters - the same
    # convention as Yang_2023_diclazuril_chicken_pbpk.R.
    #
    # Every value below is the `constant` value of the supplementary
    # acslX INITIAL block, which is the code that generated the published
    # predictions. Where the printed tables disagree with the code
    # (spleen and brain weight fractions, the rest-of-body balances, and
    # the cardiac-output exponent) the code is authoritative and the
    # discrepancy is documented in the vignette Errata.
    # ==================================================================
    qar         <- 6.9      # cardiac output coefficient, L/h/kg^0.74 of whole blood; Ai 2024 Section 2.9 ("6.9 L/h/kg") / supplement `QAR=6.9`
    pcv         <- 0.29     # packed cell volume, converts blood volume and flow to plasma; supplement `pcv=0.29`

    vc_liver    <- 0.0129   # liver weight / BW;          Ai 2024 Table 3 (1.29%)  / supplement `Vcli=0.0129`
    vc_kidney   <- 0.0031   # kidney weight / BW;         Ai 2024 Table 3 (0.31%)  / supplement `Vcki=0.0031`
    vc_lung     <- 0.0078   # lung weight / BW;           Ai 2024 Table 3 (0.78%)  / supplement `Vclu=0.0078`
    vc_muscle   <- 0.3527   # muscle weight / BW;         Ai 2024 Table 3 (35.27%) / supplement `Vcmu=0.3527`
    vc_adipose  <- 0.0274   # fat weight / BW;            Ai 2024 Table 3 (2.74%)  / supplement `Vcfa=0.0274`
    vc_heart    <- 0.0035   # heart weight / BW;          Ai 2024 Table 3 (0.35%)  / supplement `Vche=0.0035`
    vc_spleen   <- 0.0022   # spleen weight / BW;         supplement `Vcsp=0.0022` (Ai 2024 Table 3 prints 0.28%)
    vc_brain    <- 0.005    # brain weight / BW;          supplement `Vcbr=0.005`  (Ai 2024 Table 3 prints 0.32%)
    vc_arterial <- 0.0188   # arterial blood vol / BW;    supplement `Vcab=0.0188` (not tabulated in the paper)
    vc_venous   <- 0.0376   # venous blood vol / BW;      supplement `Vcvb=0.0376` (not tabulated in the paper)

    qc_liver    <- 0.4832   # liver plasma flow / Qtot (hepatic artery + portal vein); Ai 2024 Table 2 (48.32%) / supplement `Qcli=0.4832`
    qc_kidney   <- 0.1705   # kidney flow / Qtot;         Ai 2024 Table 2 (17.05%) / supplement `Qcki=0.1705`
    qc_muscle   <- 0.14     # muscle flow / Qtot;         Ai 2024 Table 2 (14.00%) / supplement `Qcmu=0.14`
    qc_adipose  <- 0.085    # fat flow / Qtot;            Ai 2024 Table 2 (8.50%)  / supplement `Qcfa=0.085`
    qc_heart    <- 0.0498   # heart flow / Qtot;          Ai 2024 Table 2 (4.98%)  / supplement `Qche=0.0498`
    qc_spleen   <- 0.04     # spleen flow / Qtot;         supplement `Qcsp=0.04` (not tabulated in the paper)
    qc_brain    <- 0.02     # brain flow / Qtot;          supplement `Qcbr=0.02`  (not tabulated in the paper)
    qc_lung     <- 1.0      # lung flow / Qtot = 1 (the lung sees the whole cardiac output); supplement `Qclu=1`

    vf_muscle   <- 0.01     # fraction of the muscle that is vascular plasma;      supplement `Vfmu=0.01`
    vf_adipose  <- 0.005    # fraction of the fat that is vascular plasma;         supplement `Vffa=0.005`
    vf_brain    <- 0.01     # fraction of the brain that is vascular plasma;       supplement `Vfbr=0.01`
    vf_other    <- 0.02     # fraction of the rest of body that is vascular plasma; supplement `Vfre=0.02`

    # ------------------------------------------------------------------
    # Drug-specific parameters, back-transformed from the log scale.
    # ------------------------------------------------------------------
    kst    <- exp(lkst)
    ka     <- exp(lka)
    kgut   <- exp(lkgut)
    clhe   <- exp(lclhe)
    clre   <- exp(lclre)
    pbind  <- exp(lpbind)
    k_liv  <- exp(lk_liv)
    k_kid  <- exp(lk_kid)
    k_lun  <- exp(lk_lun)
    k_spl  <- exp(lk_spl)
    k_hrt  <- exp(lk_hrt)
    k_mus  <- exp(lk_mus)
    k_adi  <- exp(lk_adi)
    k_bra  <- exp(lk_bra)
    k_res  <- exp(lk_res)
    pp_mus <- exp(lpp_mus)
    pp_adi <- exp(lpp_adi)
    pp_bra <- exp(lpp_bra)
    pp_res <- exp(lpp_res)

    fu     <- 1 - pbind                    # unbound fraction; supplement `pfree=1-pbind`
    f_oral <- ka / (ka + kgut)              # bioavailability; Ai 2024 Section 2.5 / supplement `F=ka/(ka+kint)`

    # ------------------------------------------------------------------
    # Compartment volumes (L, or kg taken as L at unit tissue density).
    # The rest-of-body fraction is the balance of body weight after every
    # named tissue INCLUDING arterial and venous blood, exactly as the
    # supplement computes it (`Vcre=1-(Vcli+Vcki+Vclu+Vcmu+Vcfa+Vche+
    # Vcab+Vcvb+Vcsp+Vcbr)`); this evaluates to 0.529 rather than the
    # 0.5866 printed in Ai 2024 Table 3, which omits the blood volumes.
    # ------------------------------------------------------------------
    vc_other    <- 1 - (vc_liver + vc_kidney + vc_lung + vc_muscle +
                        vc_adipose + vc_heart + vc_arterial + vc_venous +
                        vc_spleen + vc_brain)
    v_liver     <- vc_liver   * WT
    v_kidney    <- vc_kidney  * WT
    v_lung      <- vc_lung    * WT
    v_muscle    <- vc_muscle  * WT
    v_adipose   <- vc_adipose * WT
    v_heart     <- vc_heart   * WT
    v_spleen    <- vc_spleen  * WT
    v_brain     <- vc_brain   * WT
    v_other     <- vc_other   * WT

    # Arterial and venous PLASMA volumes: the blood volume less the
    # packed cell volume (supplement `Vap=Vab*(1-pcv)`).
    v_arterial  <- vc_arterial * WT * (1 - pcv)
    v_venous    <- vc_venous   * WT * (1 - pcv)

    # Membrane-limited tissues split into a vascular plasma sub-volume and
    # a tissue sub-volume (supplement `Vmub=Vmu*Vfmu; Vmut=Vmu-Vmub`).
    v_vp_muscle   <- v_muscle  * vf_muscle
    v_int_muscle  <- v_muscle  - v_vp_muscle
    v_vp_adipose  <- v_adipose * vf_adipose
    v_int_adipose <- v_adipose - v_vp_adipose
    v_vp_brain    <- v_brain   * vf_brain
    v_int_brain   <- v_brain   - v_vp_brain
    v_vp_other    <- v_other   * vf_other
    v_int_other   <- v_other   - v_vp_other

    # ------------------------------------------------------------------
    # Plasma flows (L/h). Cardiac output is allometrically scaled and
    # converted from blood to plasma (supplement
    # `QTOT=QAR*bw**0.74*(1-pcv)`). The rest-of-body flow is the balance
    # after the named tissues EXCLUDING the spleen, whose venous outflow
    # drains into the liver and is therefore already inside Qcli
    # (supplement `Qcre=1-(Qcli+Qcki+Qcmu+Qcfa+Qche+Qcbr)`); this
    # evaluates to 0.0515 rather than the 0.0715 printed in Ai 2024
    # Table 2, and is the value for which the venous-return flows sum
    # exactly to the cardiac output.
    # ------------------------------------------------------------------
    qc_other   <- 1 - (qc_liver + qc_kidney + qc_muscle + qc_adipose +
                       qc_heart + qc_brain)
    q_total    <- qar * WT^0.74 * (1 - pcv)
    q_liver    <- q_total * qc_liver
    q_kidney   <- q_total * qc_kidney
    q_muscle   <- q_total * qc_muscle
    q_adipose  <- q_total * qc_adipose
    q_heart    <- q_total * qc_heart
    q_spleen   <- q_total * qc_spleen
    q_brain    <- q_total * qc_brain
    q_lung     <- q_total * qc_lung
    q_other    <- q_total * qc_other

    # Permeability-area products (L/h); supplement `PAmu=PPmu*Qmu` etc.
    pa_muscle  <- pp_mus * q_muscle
    pa_adipose <- pp_adi * q_adipose
    pa_brain   <- pp_bra * q_brain
    pa_other   <- pp_res * q_other

    # Clearances scaled to the animal (L/h); supplement `CCLhe=clhe*bw`.
    cl_hepatic <- clhe * WT
    cl_renal   <- clre * WT

    # ------------------------------------------------------------------
    # Compartment concentrations (ug/L, or ug/kg at unit tissue density).
    # For the membrane-limited tissues the whole-tissue concentration
    # reported by the paper is the summed amount over the whole organ
    # volume (supplement `amu=amub+amut; Cmu=amu/Vmu`).
    # ------------------------------------------------------------------
    Cc            <- venous   / v_venous          # venous plasma: the sampled matrix
    c_arterial    <- arterial / v_arterial
    Cliver        <- liver    / v_liver
    Ckidney       <- kidney   / v_kidney
    Clung         <- lung     / v_lung
    Cspleen       <- spleen   / v_spleen
    Cheart        <- heart    / v_heart
    c_vp_muscle   <- vp_muscle   / v_vp_muscle
    c_int_muscle  <- int_muscle  / v_int_muscle
    c_vp_adipose  <- vp_adipose  / v_vp_adipose
    c_int_adipose <- int_adipose / v_int_adipose
    c_vp_brain    <- vp_brain    / v_vp_brain
    c_int_brain   <- int_brain   / v_int_brain
    c_vp_other    <- vp_other    / v_vp_other
    c_int_other   <- int_other   / v_int_other
    Cmuscle       <- (vp_muscle  + int_muscle)  / v_muscle
    Cadipose      <- (vp_adipose + int_adipose) / v_adipose
    Cbrain        <- (vp_brain   + int_brain)   / v_brain
    c_other       <- (vp_other   + int_other)   / v_other

    # ------------------------------------------------------------------
    # Mass-balance ODEs (Ai 2024 Table 1; supplementary acslX DERIVATIVE
    # block, which is authoritative -- several rows of the printed
    # Table 1 are garbled, see the vignette Errata).
    #
    # Dosing: the daily oral gavage dose (ug) is delivered to `stomach`.
    # In the acslX code this is a 0.001 h zero-order pulse repeated every
    # 24 h; here it is carried by the event table, so the schedule is a
    # property of the data, not of the model.
    # ------------------------------------------------------------------
    d/dt(stomach) <- -kst * stomach                                        # supplement `rsto=Rdosepo-kst*asto`
    d/dt(a_gut)   <- kst * stomach - kgut * a_gut - ka * f_oral * a_gut     # supplement `rgut=kst*asto-kint*agut-ka*F*agut`

    # Liver: hepatic-artery-plus-portal inflow, splenic venous inflow,
    # the absorbed oral dose (first pass), and hepatic elimination of the
    # unbound drug.
    d/dt(liver)   <- (q_liver - q_spleen) * c_arterial +
                     q_spleen * Cspleen / k_spl -
                     q_liver * Cliver / k_liv +
                     ka * f_oral * a_gut -
                     cl_hepatic * Cliver / k_liv * fu                      # supplement `rali=...`
    d/dt(spleen)  <- q_spleen * (c_arterial - Cspleen / k_spl)             # supplement `rasp=Qsp*(Cap-Csp/Psp)`
    d/dt(kidney)  <- q_kidney * (c_arterial - Ckidney / k_kid) -
                     cl_renal * Ckidney / k_kid * fu                       # supplement `raki=...`
    d/dt(urine)   <- cl_renal * Ckidney / k_kid * fu                       # supplement `reki`, integrated as `eki`
    d/dt(heart)   <- q_heart * (c_arterial - Cheart / k_hrt)               # supplement `rahe=Qhe*(Cap-Che/Phe)`

    # Membrane-limited tissues. The vascular sub-compartment exchanges
    # with arterial plasma by flow and with the tissue sub-compartment
    # across the membrane; at equilibrium C_tissue = P * C_vascular.
    d/dt(vp_muscle)   <- q_muscle * (c_arterial - c_vp_muscle) +
                         pa_muscle * c_int_muscle / k_mus -
                         pa_muscle * c_vp_muscle                           # supplement `ramub=...`
    d/dt(int_muscle)  <- -pa_muscle * c_int_muscle / k_mus +
                         pa_muscle * c_vp_muscle                           # supplement `ramut=...`
    d/dt(vp_adipose)  <- q_adipose * (c_arterial - c_vp_adipose) +
                         pa_adipose * c_int_adipose / k_adi -
                         pa_adipose * c_vp_adipose                         # supplement `rafab=...`
    d/dt(int_adipose) <- -pa_adipose * c_int_adipose / k_adi +
                         pa_adipose * c_vp_adipose                         # supplement `rafat=...`
    d/dt(vp_brain)    <- q_brain * (c_arterial - c_vp_brain) +
                         pa_brain * c_int_brain / k_bra -
                         pa_brain * c_vp_brain                             # supplement `rabrb=...`
    d/dt(int_brain)   <- -pa_brain * c_int_brain / k_bra +
                         pa_brain * c_vp_brain                             # supplement `rabrt=...`
    d/dt(vp_other)    <- q_other * (c_arterial - c_vp_other) +
                         pa_other * c_int_other / k_res -
                         pa_other * c_vp_other                             # supplement `rareb=...`
    d/dt(int_other)   <- -pa_other * c_int_other / k_res +
                         pa_other * c_vp_other                             # supplement `raret=...`

    # Lung sees the whole venous return; arterial plasma is its outflow.
    d/dt(lung)     <- q_lung * (Cc - Clung / k_lun)                        # supplement `ralu=Qlu*(Cvp-Clu/Plu)`
    d/dt(arterial) <- q_lung * (Clung / k_lun - c_arterial)                # supplement `raap=Qlu*(Clu/Plu-Cap)`

    # Venous plasma collects every tissue's venous outflow. The
    # membrane-limited tissues return their vascular sub-compartment
    # concentration directly (it is already plasma); the flow-limited
    # tissues return C/P. The spleen does not appear because it drains
    # into the liver.
    d/dt(venous)   <- q_other   * c_vp_other +
                      q_brain   * c_vp_brain +
                      q_adipose * c_vp_adipose +
                      q_heart   * Cheart  / k_hrt +
                      q_muscle  * c_vp_muscle +
                      q_kidney  * Ckidney / k_kid +
                      q_liver   * Cliver  / k_liv -
                      q_lung    * Cc                                       # supplement `ravp=...`

    # ------------------------------------------------------------------
    # Observations. Venous plasma plus the eight tissues Ai 2024 compared
    # observed against predicted (Table S4 / Figures 5-7). Muscle, fat,
    # liver and kidney are the four Codex Alimentarius MRL tissues from
    # which the withdrawal times of Ai 2024 Table 7 are derived.
    # ------------------------------------------------------------------
    Cc       ~ prop(propSd)
    Cliver   ~ prop(propSd_Cliver)
    Ckidney  ~ prop(propSd_Ckidney)
    Cmuscle  ~ prop(propSd_Cmuscle)
    Cadipose ~ prop(propSd_Cadipose)
    Clung    ~ prop(propSd_Clung)
    Cspleen  ~ prop(propSd_Cspleen)
    Cheart   ~ prop(propSd_Cheart)
    Cbrain   ~ prop(propSd_Cbrain)
  })
}

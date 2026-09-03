Jones_2013_perfusionLimited_pbpk <- function() {
  description <- paste(
    "Methodology reference (whole-body PBPK, perfusion-rate-limited, 16",
    "ODEs; Berkeley Madonna). The generic human whole-body PBPK model",
    "deposited as the Supplementary Data of the Jones and Rowland-Yeo",
    "2013 CPT:PSP tutorial on basic concepts in PBPK modelling, offered",
    "there as the worked alternative to running the same structure on a",
    "commercial platform. Thirteen perfusion-rate-limited tissues",
    "(adipose, bone, brain, gut, heart, kidney, liver, lung, muscle,",
    "skin, spleen, testes and a lumped rest-of-body) are strung between",
    "arterial and venous blood pools; gut and spleen drain into the",
    "liver through the portal vein and the liver is entered additionally",
    "by the hepatic artery, so the whole splanchnic bed leaves through",
    "the hepatic vein. Elimination is hepatic metabolism, scaled from",
    "human liver microsomal intrinsic clearance through MPPGL and the",
    "microsomal unbound fraction, plus an optional renal arm; both are",
    "driven by an unbound tissue concentration. Oral dosing is a",
    "first-order depot emptying into the gut, so a first-pass extraction",
    "arises structurally rather than through a bioavailability term.",
    "There is no drug: every tissue-to-plasma partition coefficient, the",
    "blood-to-plasma ratio and both unbound fractions are set to 1 by the",
    "authors, which makes the model a neutral template into which a",
    "compound's own physicochemically predicted Kp vector, binding data",
    "and intrinsic clearance are substituted. The physiology is a single",
    "70 kg reference adult and both published fraction columns balance",
    "exactly: the fifteen fractional tissue volumes sum to 1.000000 and",
    "the arterial-side fractional blood flows sum to 1.000000.",
    "Deterministic typical-value model: the source reports no",
    "between-subject variability and no residual-error model.",
    sep = " "
  )
  reference <- paste(
    "Jones HM, Rowland-Yeo K. Basic concepts in physiologically based",
    "pharmacokinetic modeling in drug discovery and development. CPT",
    "Pharmacometrics Syst Pharmacol. 2013;2(8):e63.",
    "doi:10.1038/psp.2013.41. PMID 23945604. PMCID PMC3828005.",
    "The model is transcribed in full from the Supplementary Data file",
    "(psp201341x1.doc), sections {HUMAN PBPK MODEL}, {SPECIES SPECIFIC",
    "PARAMETERS}, {COMPOUND SPECIFIC PARAMETERS} and {DOSING}; the",
    "main-text narrative of the perfusion-limited mass balance is",
    "Equations 1 and 2 and the microsomal clearance scaling is Equation",
    "3. See the vignette Errata for the one divergence between the",
    "main-text equation and the deposited code.",
    sep = " "
  )
  vignette <- "Jones_2013_perfusionLimited_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Both routes the deposited code initialises: an oral dose into the gut
  # lumen (INIT D = PODOSE) and an IV bolus into venous blood
  # (INIT Ave = IVDOSE). `venous` is outside the auto-detected
  # depot/central set, so the dosing compartments are declared explicitly.
  dosing <- c("depot", "venous")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: each state is named by the
  # deposited code's own inline comment on its `INIT A..` line.
  compartmentData <- list(
    depot    = list(analyte = "generic small molecule", units = "mg", specimen = "administration site", verified = TRUE),
    adipose  = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    bone     = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    brain    = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    gut      = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    heart    = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    kidney   = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    liver    = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    lung     = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    muscle   = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    skin     = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    spleen   = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    testes   = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    venous   = list(analyte = "generic small molecule", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "generic small molecule", units = "mg", specimen = "whole blood", verified = TRUE),
    other    = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE)
  )

  # No covariates: a single 70 kg reference adult with no demographic or
  # laboratory inputs. Body weight is carried as the system constant `bw`
  # rather than as a covariate because the deposited code hard-codes it.
  covariateData <- list()

  population <- list(
    species       = "None (methodology/tutorial paper; human reference physiology carrying a hypothetical generic compound, no patients and no fitted estimates).",
    n_subjects    = 0L,
    n_studies     = 0L,
    disease_state = "N/A (deterministic forward simulation; not a fit of any real molecule).",
    dose_range    = "Supplementary Data {DOSING}: PODOSE = 100 mg oral bolus, IVDOSE = 0 mg. The 100 mg oral dose matches the proposed efficacious human dose simulated for the tutorial's illustrative compound X (main text, Example 1).",
    regions       = "N/A",
    scope_note    = paste(
      "Filed under inst/modeldb/pharmacokinetics/ (not specificDrugs/) because there is no",
      "drug, following the operator-ratified Beal_2001_iv1cmt_bql precedent for methodology",
      "papers whose only pharmacokinetic content is an author-supplied hypothetical compound.",
      "The file stem uses a structural descriptor in the slot where a drug name would sit,",
      "matching the sibling Gaohua_2023_permeabilityLimited_pbpk. The paper itself is a",
      "tutorial and its two worked examples (an anonymised compound X and repaglinide) were",
      "run on the Simcyp Population-Based Simulator; neither is extractable, because in both",
      "cases the Kp vector, the ADAM absorption model and, for repaglinide, the",
      "permeability-limited liver sub-compartment geometry are vendor-internal and unpublished.",
      "This model file is the deposited Berkeley Madonna code, which is complete and",
      "platform-free: every constant it consumes is printed in the supplement.",
      sep = " "
    ),
    notes         = paste(
      "The compound is deliberately neutral. Supplementary Data {COMPOUND SPECIFIC PARAMETERS}",
      "sets all thirteen Kp values, the blood-to-plasma ratio, the plasma unbound fraction and",
      "the microsomal unbound fraction to exactly 1, the absorption rate constant and absorbed",
      "fraction to 1, renal clearance to 0, and human liver microsomal intrinsic clearance to",
      "10 uL/min/mg. Substituting a real compound's values is the intended use. Because Kp = 1",
      "and BP = 1 throughout, the model's steady-state volume of distribution is the whole",
      "70 kg body volume by construction, and hepatic extraction is governed entirely by the",
      "ratio of the scaled metabolic clearance to hepatic blood flow.",
      sep = " "
    )
  )

  ini({
    # ============ System-specific parameters (human) ============
    # Supplementary Data {SPECIES SPECIFIC PARAMETERS}. A tissue density
    # of 1 kg/L is implied throughout: volumes in L are computed directly
    # as bw * fractional volume.
    bw <- fixed(70);     label("Body weight (kg)")                    # suppl {SPECIES SPECIFIC PARAMETERS} `BW = 70 ; BW (kg)`
    co <- fixed(108.33); label("Cardiac output (mL/s)")               # suppl {Total tissue blood flows - L/hr} `CO = 108.33 ; cardiac output (ml/s)`

    # Fractional tissue volumes (fraction of body weight). The fifteen
    # compartment fractions sum to exactly 1.000000; fvol_plasma and
    # fvol_erythrocytes are the composition of blood and are therefore
    # NOT part of that sum.
    fvol_adipose      <- fixed(0.213);    label("Fractional volume, adipose (L/kg body weight)")        # suppl `FVad = 0.213 ; adipose`
    fvol_bone         <- fixed(0.085629); label("Fractional volume, bone (L/kg body weight)")           # suppl `FVbo = 0.085629 ; bone`
    fvol_brain        <- fixed(0.02);     label("Fractional volume, brain (L/kg body weight)")          # suppl `FVbr = 0.02 ; brain`
    fvol_gut          <- fixed(0.0171);   label("Fractional volume, gut (L/kg body weight)")            # suppl `FVgu = 0.0171 ; gut`
    fvol_heart        <- fixed(0.0047);   label("Fractional volume, heart (L/kg body weight)")          # suppl `FVhe = 0.0047 ; heart`
    fvol_kidney       <- fixed(0.0044);   label("Fractional volume, kidney (L/kg body weight)")         # suppl `FVki = 0.0044 ; kidney`
    fvol_liver        <- fixed(0.021);    label("Fractional volume, liver (L/kg body weight)")          # suppl `FVli = 0.021 ; liver`
    fvol_lung         <- fixed(0.0076);   label("Fractional volume, lung (L/kg body weight)")           # suppl `FVlu = 0.0076 ; lung`
    fvol_muscle       <- fixed(0.4);      label("Fractional volume, muscle (L/kg body weight)")         # suppl `FVmu = 0.4 ; muscle`
    fvol_skin         <- fixed(0.0371);   label("Fractional volume, skin (L/kg body weight)")           # suppl `FVsk = 0.0371 ; skin`
    fvol_spleen       <- fixed(0.0026);   label("Fractional volume, spleen (L/kg body weight)")         # suppl `FVsp = 0.0026 ; spleen`
    fvol_testes       <- fixed(0.01);     label("Fractional volume, testes (L/kg body weight)")         # suppl `FVte = 0.01 ; testes`
    fvol_venous       <- fixed(0.0514);   label("Fractional volume, venous blood (L/kg body weight)")   # suppl `FVve = 0.0514 ; venous`
    fvol_arterial     <- fixed(0.0257);   label("Fractional volume, arterial blood (L/kg body weight)") # suppl `FVar = 0.0257 ; arterial`
    fvol_other        <- fixed(0.099771); label("Fractional volume, rest of body (L/kg body weight)")   # suppl `FVre = 0.099771 ; rest of body`
    fvol_plasma       <- fixed(0.0424);   label("Fractional volume, plasma (L/kg body weight)")         # suppl `FVpl = 0.0424 ; plasma`
    fvol_erythrocytes <- fixed(0.0347);   label("Fractional volume, erythrocytes (L/kg body weight)")   # suppl `FVrb = 0.0347 ; erythrocytes`

    # Fractional tissue blood flows (fraction of cardiac output). The
    # arterial-side shares sum to exactly 1.000000 once the hepatic
    # artery is resolved as fq_hepatic - fq_gut - fq_spleen, because
    # fq_hepatic is the venous-side (total hepatic outflow) share.
    fq_adipose <- fixed(0.05);     label("Fractional blood flow, adipose (fraction of cardiac output)")            # suppl `FQad = 0.05 ; adipose`
    fq_bone    <- fixed(0.05);     label("Fractional blood flow, bone (fraction of cardiac output)")               # suppl `FQbo = 0.05 ; bone`
    fq_brain   <- fixed(0.12);     label("Fractional blood flow, brain (fraction of cardiac output)")              # suppl `FQbr = 0.12 ; brain`
    fq_gut     <- fixed(0.146462); label("Fractional blood flow, gut (fraction of cardiac output)")                # suppl `FQgu = 0.146462 ; gut`
    fq_heart   <- fixed(0.04);     label("Fractional blood flow, heart (fraction of cardiac output)")              # suppl `FQhe = 0.04 ; heart`
    fq_kidney  <- fixed(0.19);     label("Fractional blood flow, kidney (fraction of cardiac output)")             # suppl `FQki = 0.19 ; kidney`
    fq_hepatic <- fixed(0.215385); label("Fractional blood flow, total hepatic outflow (fraction of cardiac output)")  # suppl `FQh = 0.215385 ; hepatic (venous side)`
    fq_lung    <- fixed(1);        label("Fractional blood flow, lung (fraction of cardiac output)")               # suppl `FQlu = 1 ; lung`
    fq_muscle  <- fixed(0.17);     label("Fractional blood flow, muscle (fraction of cardiac output)")             # suppl `FQmu = 0.17 ; muscle`
    fq_skin    <- fixed(0.05);     label("Fractional blood flow, skin (fraction of cardiac output)")               # suppl `FQsk = 0.05 ; skin`
    fq_spleen  <- fixed(0.017231); label("Fractional blood flow, spleen (fraction of cardiac output)")             # suppl `FQsp = 0.017231 ; spleen`
    fq_testes  <- fixed(0.01076);  label("Fractional blood flow, testes (fraction of cardiac output)")             # suppl `FQte = 0.01076 ; testes`
    fq_other   <- fixed(0.103855); label("Fractional blood flow, rest of body (fraction of cardiac output)")       # suppl `FQre = 0.103855 ; rest of body`

    # ============ Compound-specific parameters ============
    # Supplementary Data {COMPOUND SPECIFIC PARAMETERS}. Every value here
    # is a neutral placeholder chosen by the authors so that the template
    # runs; they are not estimates of anything and are the values a user
    # replaces with a real compound's data.
    kp_adipose <- fixed(1); label("Tissue-to-plasma partition coefficient, adipose")        # suppl `Kpad = 1 ; adipose`
    kp_bone    <- fixed(1); label("Tissue-to-plasma partition coefficient, bone")           # suppl `Kpbo = 1 ; bone`
    kp_brain   <- fixed(1); label("Tissue-to-plasma partition coefficient, brain")          # suppl `Kpbr = 1 ; brain`
    kp_gut     <- fixed(1); label("Tissue-to-plasma partition coefficient, gut")            # suppl `Kpgu = 1 ; gut`
    kp_heart   <- fixed(1); label("Tissue-to-plasma partition coefficient, heart")          # suppl `Kphe = 1 ; heart`
    kp_kidney  <- fixed(1); label("Tissue-to-plasma partition coefficient, kidney")         # suppl `Kpki = 1 ; kidney`
    kp_liver   <- fixed(1); label("Tissue-to-plasma partition coefficient, liver")          # suppl `Kpli = 1 ; liver`
    kp_lung    <- fixed(1); label("Tissue-to-plasma partition coefficient, lung")           # suppl `Kplu = 1 ; lung`
    kp_muscle  <- fixed(1); label("Tissue-to-plasma partition coefficient, muscle")         # suppl `Kpmu = 1 ; muscle`
    kp_skin    <- fixed(1); label("Tissue-to-plasma partition coefficient, skin")           # suppl `Kpsk = 1 ; skin`
    kp_spleen  <- fixed(1); label("Tissue-to-plasma partition coefficient, spleen")         # suppl `Kpsp = 1 ; spleen`
    kp_testes  <- fixed(1); label("Tissue-to-plasma partition coefficient, testes")         # suppl `Kpte = 1 ; testes`
    kp_other   <- fixed(1); label("Tissue-to-plasma partition coefficient, rest of body")   # suppl `Kpre = 1 ; rest of body`

    fup   <- fixed(1); label("Fraction unbound in plasma")                                  # suppl {In vitro binding data} `fup = 1 ; fraction unbound in plasma`
    bp    <- fixed(1); label("Blood-to-plasma concentration ratio")                          # suppl {In vitro binding data} `BP = 1 ; blood to plasma ratio`
    fumic <- fixed(1); label("Fraction unbound in liver microsomes")                         # suppl {In vitro binding data} `fumic = 1 ; fraction unbound in microsomes`

    mppgl    <- fixed(45); label("Microsomal protein per gram of liver (mg/g)")              # suppl {Clearance calculations} `MPPGL = 45 ; mg microsomal protein per g liver`
    clint    <- fixed(10); label("Apparent human liver microsomal intrinsic clearance (uL/min/mg protein)")  # suppl {Clearances} `HLM_CLint = 10 ; HLM CLint apparent (ul/min/mg)`
    cl_renal <- fixed(0);  label("Renal clearance acting on unbound kidney concentration (L/h)")             # suppl {Clearances} `CLrenal = 0 ; CLint renal (L/hr)`

    ka      <- fixed(1); label("First-order absorption rate constant (1/h)")                 # suppl {Absorption} `Ka = 1 ; Ka (hr-1)`
    fdepot  <- fixed(1); label("Fraction of the depot dose absorbed")                        # suppl {Absorption} `F  = 1 ; fraction absorbed`
  })

  model({
    # ================= Derived system physiology =================
    # Cardiac output, mL/s -> L/h.
    qc <- co / 1000 * 60 * 60  # suppl `QC = CO/1000*60*60 ; cardiac output (L/hr)`

    # Tissue volumes (L). suppl {Total tissue volumes - L}: Vx = BW*FVx.
    v_adipose      <- bw * fvol_adipose
    v_bone         <- bw * fvol_bone
    v_brain        <- bw * fvol_brain
    v_gut          <- bw * fvol_gut
    v_heart        <- bw * fvol_heart
    v_kidney       <- bw * fvol_kidney
    v_liver        <- bw * fvol_liver
    v_lung         <- bw * fvol_lung
    v_muscle       <- bw * fvol_muscle
    v_skin         <- bw * fvol_skin
    v_spleen       <- bw * fvol_spleen
    v_testes       <- bw * fvol_testes
    v_venous       <- bw * fvol_venous
    v_arterial     <- bw * fvol_arterial
    v_other        <- bw * fvol_other
    v_plasma       <- bw * fvol_plasma
    v_erythrocytes <- bw * fvol_erythrocytes

    # Venous and arterial plasma volumes. Carried because the deposited
    # code computes them; the perfusion-limited mass balance below works
    # in whole-blood volumes and does not consume them.
    v_plasma_venous   <- v_plasma * v_venous / (v_venous + v_arterial)    # suppl `Vplas_ven = Vpl*Vve/(Vve + Var)`
    v_plasma_arterial <- v_plasma * v_arterial / (v_venous + v_arterial)  # suppl `Vplas_art = Vpl*Var/(Vve + Var)`

    # Tissue blood flows (L/h). suppl {Total tissue blood flows - L/hr}.
    q_adipose <- qc * fq_adipose
    q_bone    <- qc * fq_bone
    q_brain   <- qc * fq_brain
    q_gut     <- qc * fq_gut
    q_heart   <- qc * fq_heart
    q_kidney  <- qc * fq_kidney
    q_hepatic <- qc * fq_hepatic
    q_lung    <- qc * fq_lung
    q_muscle  <- qc * fq_muscle
    q_skin    <- qc * fq_skin
    q_spleen  <- qc * fq_spleen
    q_testes  <- qc * fq_testes
    q_other   <- qc * fq_other

    # Hepatic artery is the balance of total hepatic outflow after the
    # two portal tributaries. suppl `Qha = Qh - Qgu - Qsp`.
    q_ha <- q_hepatic - q_gut - q_spleen

    # ================= Concentrations (mg/L) =================
    # suppl {Calculation of total concentrations - mg/L}.
    c_adipose  <- adipose / v_adipose
    c_bone     <- bone / v_bone
    c_brain    <- brain / v_brain
    c_gut      <- gut / v_gut
    c_heart    <- heart / v_heart
    c_kidney   <- kidney / v_kidney
    c_liver    <- liver / v_liver
    c_lung     <- lung / v_lung
    c_muscle   <- muscle / v_muscle
    c_skin     <- skin / v_skin
    c_spleen   <- spleen / v_spleen
    c_testes   <- testes / v_testes
    c_venous   <- venous / v_venous
    c_arterial <- arterial / v_arterial
    c_other    <- other / v_other

    # Venous plasma concentration - the model's observable.
    c_venous_plasma <- c_venous / bp  # suppl `Cplasmavenous = Cvenous/BP ; venous plasma`

    # Free concentrations driving the two elimination arms. Transcribed
    # exactly as the deposited code writes them (total tissue
    # concentration times fup). See vignette Errata: main-text Eq 2
    # drives elimination with the free concentration in the emergent
    # venous blood, C_T/(Kp/BP)*fup, which coincides with this only
    # because the template sets Kp = BP = 1.
    cu_liver  <- c_liver * fup    # suppl {Calculation of free concentrations - mg/L} `Cliverfree = Cliver*fup`
    cu_kidney <- c_kidney * fup   # suppl `Ckidneyfree = Ckidney*fup`

    # ================= Clearance and absorption =================
    # Microsomal scaling, main-text Eq 3. The factor 60/1000 is
    # (1000 g liver per L) * (60 min per h) / (1e6 uL per L).
    cl_met <- (clint / fumic) * mppgl * v_liver * 60 / 1000  # suppl `CLmet = (HLM_CLint/fumic)*MPPGL*Vli*60/1000 ; CLint scaled (L/hr)`

    absorption <- ka * depot * fdepot  # suppl `Absorption = Ka*D*F`

    # Blood returning to the venous pool. Gut and spleen are absent
    # because they drain into the liver; the liver returns the whole
    # splanchnic bed at q_hepatic. suppl `Venous = ...`.
    venous_return <-
      q_adipose * (c_adipose / kp_adipose * bp) +
      q_bone * (c_bone / kp_bone * bp) +
      q_brain * (c_brain / kp_brain * bp) +
      q_heart * (c_heart / kp_heart * bp) +
      q_kidney * (c_kidney / kp_kidney * bp) +
      q_hepatic * (c_liver / kp_liver * bp) +
      q_muscle * (c_muscle / kp_muscle * bp) +
      q_skin * (c_skin / kp_skin * bp) +
      q_testes * (c_testes / kp_testes * bp) +
      q_other * (c_other / kp_other * bp)

    # ================= Mass balance (mg/h) =================
    # suppl {Differential equations - mg/hr}. Non-eliminating tissues
    # follow main-text Eq 1; liver and kidney follow Eq 2.
    d/dt(depot)    <- -absorption                                                      # suppl `d/dt(D) = - Absorption ; oral dosing`
    d/dt(adipose)  <- q_adipose * (c_arterial - c_adipose / kp_adipose * bp)           # suppl `d/dt(Aad) = Qad*(Carterial - Cadipose/Kpad*BP)`
    d/dt(bone)     <- q_bone * (c_arterial - c_bone / kp_bone * bp)                    # suppl `d/dt(Abo) = Qbo*(Carterial - Cbone/Kpbo*BP)`
    d/dt(brain)    <- q_brain * (c_arterial - c_brain / kp_brain * bp)                 # suppl `d/dt(Abr) = Qbr*(Carterial - Cbrain/Kpbr*BP)`
    d/dt(gut)      <- absorption + q_gut * (c_arterial - c_gut / kp_gut * bp)          # suppl `d/dt(Agu) = Absorption + Qgu*(Carterial - Cgut/Kpgu*BP)`
    d/dt(heart)    <- q_heart * (c_arterial - c_heart / kp_heart * bp)                 # suppl `d/dt(Ahe) = Qhe*(Carterial - Cheart/Kphe*BP)`
    d/dt(kidney)   <- q_kidney * (c_arterial - c_kidney / kp_kidney * bp) -
      cl_renal * cu_kidney                                                             # suppl `d/dt(Aki) = Qki*(Carterial - Ckidney/Kpki*BP) - CLrenal*Ckidneyfree`
    d/dt(liver)    <- q_ha * c_arterial +
      q_gut * (c_gut / kp_gut * bp) +
      q_spleen * (c_spleen / kp_spleen * bp) -
      q_hepatic * (c_liver / kp_liver * bp) -
      cu_liver * cl_met                                                                # suppl `d/dt(Ali)= Qha*Carterial + Qgu*(Cgut/Kpgu*BP) + Qsp*(Cspleen/Kpsp*BP) - Qh*(Cliver/Kpli*BP) - Cliverfree*CLmet`
    d/dt(lung)     <- q_lung * c_venous - q_lung * (c_lung / kp_lung * bp)             # suppl `d/dt(Alu) = Qlu*Cvenous - Qlu*(Clung/Kplu*BP)`
    d/dt(muscle)   <- q_muscle * (c_arterial - c_muscle / kp_muscle * bp)              # suppl `d/dt(Amu) = Qmu*(Carterial - Cmuscle/Kpmu*BP)`
    d/dt(skin)     <- q_skin * (c_arterial - c_skin / kp_skin * bp)                    # suppl `d/dt(Ask) = Qsk*(Carterial - Cskin/Kpsk*BP)`
    d/dt(spleen)   <- q_spleen * (c_arterial - c_spleen / kp_spleen * bp)              # suppl `d/dt(Asp) = Qsp*(Carterial - Cspleen/Kpsp*BP)`
    d/dt(testes)   <- q_testes * (c_arterial - c_testes / kp_testes * bp)              # suppl `d/dt(Ate) = Qte*(Carterial - Ctestes/Kpte*BP)`
    d/dt(venous)   <- venous_return - q_lung * c_venous                                # suppl `d/dt(Ave) = Venous - Qlu*Cvenous`
    d/dt(arterial) <- q_lung * (c_lung / kp_lung * bp) - q_lung * c_arterial           # suppl `d/dt(Aar) = Qlu*(Clung/Kplu*BP) - Qlu*Carterial`
    d/dt(other)    <- q_other * (c_arterial - c_other / kp_other * bp)                 # suppl `d/dt(Are) = Qre*(Carterial - Crest/Kpre*BP)`

    # ================= Observation =================
    # Venous plasma concentration. The source reports no residual-error
    # model: this is a deterministic forward-simulation template.
    Cc <- c_venous_plasma
  })
}

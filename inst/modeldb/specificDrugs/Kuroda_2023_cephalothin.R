Kuroda_2023_cephalothin <- function() {
  description <- paste(
    "Veterinary (Thoroughbred horse). Three-compartment population PK model for the",
    "first-generation cephalosporin cephalothin (CET) in Thoroughbred horses, fitted",
    "jointly to single-dose intramuscular (11 mg/kg bwt, 8 horses, this study) and",
    "single-dose intravenous (22 mg/kg bwt, 12 horses, Kuroda 2021 Equine Vet J",
    "53:1239-1249) plasma data by nonlinear mixed-effects modelling in Phoenix",
    "WinNonlin 8.3. Absorption from the intramuscular site is first order (Kabs)",
    "with a bioavailability factor F; disposition is a central compartment (V1)",
    "connected to a rapidly equilibrating peripheral compartment (V2, distribution",
    "clearance CL2) and a slowly equilibrating peripheral compartment (V3,",
    "distribution clearance CL3), with elimination clearance CL from V1. Every",
    "structural parameter is reported per kilogram of body weight in Kuroda 2023",
    "Table 1, so the model scales each of them linearly with WT (exponent 1 fixed);",
    "doses are therefore given in mg and concentrations come out in ug/mL. The",
    "unbound concentration Cu = fu * Cc (fu = 0.8, taken from Ambrose 2007 and used",
    "unchanged by Kuroda 2023) is the quantity that drives the paper's PK/PD target,",
    "fT > MIC for 40% of the dosing interval, and its Monte Carlo probability of",
    "target attainment. Kuroda 2023 states that between-subject variability was",
    "described by an exponential model on the structural parameters but reports NO",
    "variance estimates anywhere (the Table 1 CV% column is bootstrap precision of",
    "the typical value, not BSV), so every eta is encoded as fixed(0) and the model",
    "simulates typical-value profiles unless the user supplies variances. The",
    "residual model is combined proportional plus additive; the additive component",
    "was estimated separately per route and the packaged value is the intramuscular",
    "one (see the addSd comment for the intravenous value)."
  )
  reference <- paste(
    "Kuroda T, Minamijima Y, Niwa H, Mita H, Tamura N, Fukuda K, Ohta M. (2023).",
    "Pharmacokinetic/pharmacodynamic analysis of cephalothin after intramuscular",
    "administration in Thoroughbred horses.",
    "Journal of Equine Science 34(4):111-114.",
    "doi:10.1294/jes.34.111.",
    sep = " "
  )
  vignette <- "Kuroda_2023_cephalothin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "cephalothin", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "cephalothin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "cephalothin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "cephalothin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed adult-horse body weight. Kuroda 2023 Table 1 reports every structural parameter per kilogram (l/kg, l/kg/hr), so WT enters as a linear multiplier (exponent 1) on cl, vc, q, vp, q2 and vp2 rather than as an estimated allometric term. Study horses weighed 490-570 kg.",
      source_name        = "bwt (body weight; Kuroda 2023 Methods paragraph 1)"
    )
  )

  population <- list(
    species        = "horse (Thoroughbred)",
    n_subjects     = 20L,
    n_studies      = 2L,
    age_range      = "4-9 years (intramuscular cohort)",
    weight_range   = "490-570 kg (intramuscular cohort)",
    sex_female_pct = 50,
    disease_state  = "healthy",
    dose_range     = "11 mg/kg bwt cephalothin as a single intramuscular bolus (<30 sec) into the right lateral neck, dissolved in 25 mL sterile physiological saline (this study); pooled with 22 mg/kg bwt single intravenous doses from Kuroda 2021",
    sampling       = "Plasma at 0, 5, 10, 20, 30 and 45 min and 1, 2, 3, 4, 6, 8 and 12 hr after the intramuscular dose; 10 mL blood from the left jugular vein via a 16 G catheter, LC-MS/MS assay (Nexera X2 / QTRAP 4500) with a limit of quantification of 0.1 ug/mL",
    regions        = "Japan (Equine Research Institute, Japan Racing Association, Tochigi)",
    organism       = "Streptococcus equi subsp. zooepidemicus (MIC90 0.12 mg/L) and Staphylococcus aureus (MIC90 0.5 mg/L), equine isolates (McGorum and Pirie 2010, reference 7 of Kuroda 2023)",
    notes          = "Eight Thoroughbred horses (four stallions and four mares) received the intramuscular dose; the model was fitted after aggregating those data with intravenous data from 12 different horses given 22 mg/kg (Kuroda 2021). Kuroda 2023 Table 1's caption says 'six horses' for the intramuscular arm, which conflicts with the Abstract, Methods and the Fig. 2 caption; all three state eight, and Fig. 2 plots eight profiles, so eight is used here. Ethics approvals 22-7 and 23-5 (Animal Care and Use Committee of the Equine Research Institute, Japan Racing Association). No neck pain, diarrhea or other side effects were observed."
  )

  ini({
    # =====================================================================
    # Structural disposition parameters -- Kuroda 2023 Table 1
    # =====================================================================
    # Table 1 reports every primary structural parameter normalized to body
    # weight (l/kg for volumes, l/kg/hr for clearances). They are carried
    # here on that per-kilogram basis and multiplied by WT (exponent 1) in
    # model(), which is exactly what the per-kilogram normalization means.
    # Naming map (Kuroda 2023 Fig. 1 / Table 1 -> nlmixr2lib canonical):
    #   V   -> vc            (central compartment, V1 in the Fig. 1 legend)
    #   V2  -> vp            (rapidly equilibrating peripheral, via Cl 12)
    #   V3  -> vp2           (slowly equilibrating peripheral,  via Cl 13)
    #   CL  -> cl            (plasma clearance out of V1)
    #   CL2 -> q             (V1 <-> V2 distribution clearance, Cl 12)
    #   CL3 -> q2            (V1 <-> V3 distribution clearance, Cl 13)
    lvc <- log(0.083); label("Weight-normalized central volume of distribution V1 (L/kg)")                       # Table 1: V = 0.083 l/kg (CV% 8.1; 2.5-97.5% 0.071-0.099)
    lvp <- log(0.060); label("Weight-normalized rapidly equilibrating peripheral volume V2 (L/kg)")              # Table 1: V2 = 0.060 l/kg (CV% 16.5; 2.5-97.5% 0.040-0.082)
    lvp2 <- log(0.054); label("Weight-normalized slowly equilibrating peripheral volume V3 (L/kg)")              # Table 1: V3 = 0.054 l/kg (CV% 71.1; 2.5-97.5% 0.039-0.057)
    lcl <- log(0.597); label("Weight-normalized plasma clearance CL (L/kg/h)")                                  # Table 1: CL = 0.597 l/kg/hr (CV% 4.3; 2.5-97.5% 0.522-0.597)
    lq <- log(0.106); label("Weight-normalized distribution clearance to peripheral1 CL2 (L/kg/h)")             # Table 1: CL2 = 0.106 l/kg/hr (CV% 14.5; 2.5-97.5% 0.098-0.156)
    lq2 <- log(0.018); label("Weight-normalized distribution clearance to peripheral2 CL3 (L/kg/h)")            # Table 1: CL3 = 0.018 l/kg/hr (CV% 18.1; 2.5-97.5% 0.016-0.030)

    # =====================================================================
    # Intramuscular absorption -- Kuroda 2023 Table 1
    # =====================================================================
    # Kabs and F were added to the disposition model for the intramuscular
    # arm only (Methods paragraph 2). Kabs is consistent with the reported
    # secondary parameter: log(2)/1.070 = 0.648 hr, and Table 1 reports an
    # absorption half-life of 0.65 hr.
    lka <- log(1.070); label("First-order intramuscular absorption rate constant Kabs (1/h)")                   # Table 1: Kabs = 1.070 1/hr (CV% 24.6); cross-checks against the secondary parameter absorption half-life = 0.65 hr
    lfdepot <- log(0.997); label("Intramuscular bioavailability F (fraction)")                                    # Table 1: F = 99.7% (CV% 9.8; 2.5-97.5% 67.7-99.8%); Results text repeats "bioavailability ... 99.7%"

    # =====================================================================
    # Plasma protein binding
    # =====================================================================
    # Not estimated by Kuroda 2023; taken unchanged from the cited source and
    # used to convert total plasma concentration to the free concentration
    # that drives the fT > MIC PK/PD target, hence fixed().
    fu <- fixed(0.8); label("Fraction unbound in plasma (assumed, not estimated)")                         # Methods paragraph 2: "Using this model and the free fraction value of 0.8 reported in a previous study [1]" (Ambrose 2007, Clin Infect Dis 44:79-86)

    # =====================================================================
    # Between-subject variability
    # =====================================================================
    # Kuroda 2023 Methods paragraph 2: "Interindividual variability was
    # assumed to obey a log-normal distribution, and for a given structural
    # parameter, the between-subject variability (BSV) was described using an
    # exponential model." NO variance, SD or CV of any eta is reported
    # anywhere in the paper -- the Table 1 CV% column is the bootstrap
    # precision of the TYPICAL VALUE (n = 50 replicates), not BSV, as the
    # table title states explicitly. Per the standing policy for unreported
    # IIV with structural values present, every eta is encoded as fixed(0)
    # rather than invented, so the packaged model simulates typical-value
    # profiles. Users with variance estimates can unfix any of these. The
    # paper does not say which subset of the structural parameters actually
    # carried BSV, so all eight estimated structural parameters get an eta
    # placeholder. See vignette Errata.
    etalvc ~ fixed(0)
    etalvp ~ fixed(0)
    etalvp2 ~ fixed(0)
    etalcl ~ fixed(0)
    etalq ~ fixed(0)
    etalq2 ~ fixed(0)
    etalka ~ fixed(0)
    etalfdepot ~ fixed(0)

    # =====================================================================
    # Residual unexplained variability -- Kuroda 2023 Table 1
    # =====================================================================
    # Methods paragraph 2: "The residual model was an additive plus
    # multiplicative (proportional) model", i.e. the Phoenix WinNonlin form
    # Cobs = C * (1 + CMultStdev * eps1) + stdev * eps2. A single
    # proportional term was estimated across both routes; the additive term
    # was estimated separately for the intravenous (stdev0) and intramuscular
    # (stdev1) datasets.
    propSd <- 0.117; label("Proportional residual error (fraction)")                                              # Table 1: CMultStdev (residual, proportional, IV) = 0.117 scalar (CV% 9.5; 2.5-97.5% 0.081-0.125)
    # The packaged additive term is the INTRAMUSCULAR one, matching this
    # model file's dosing route and the paper's subject. To simulate the
    # intravenous arm (dose straight into `central`), replace it with the
    # intravenous estimate stdev0 = 0.008 ug/mL (Table 1, CV% 28.8).
    # UNITS: Table 1 labels stdev0/stdev1 "ug/l". That label is a typo --
    # the assay's limit of quantification is 0.1 ug/mL and Fig. 2 plots
    # observed concentrations in ug/mL, so an additive SD of 0.109 ug/L
    # would sit three orders of magnitude below the quantification limit.
    # The values are in the assay's concentration units, ug/mL (= mg/L),
    # which is also what Phoenix's stdev parameter carries. See vignette
    # Errata.
    addSd <- 0.109; label("Additive residual error, intramuscular (ug/mL)")                                       # Table 1: stdev1 (residual, additive, IM) = 0.109 (CV% 50.0; 2.5-97.5% 0.002-0.109)
  })

  model({
    # 1. Individual parameters. Kuroda 2023 Table 1 is reported per kilogram
    #    of body weight, so every volume and clearance scales linearly with
    #    WT (exponent 1). Kabs and F are not weight-normalized in Table 1
    #    (1/hr and %) and are used as-is.
    vc <- exp(lvc + etalvc) * WT
    vp <- exp(lvp + etalvp) * WT
    vp2 <- exp(lvp2 + etalvp2) * WT
    cl <- exp(lcl + etalcl) * WT
    q <- exp(lq + etalq) * WT
    q2 <- exp(lq2 + etalq2) * WT
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system (Kuroda 2023 Fig. 1). The intramuscular dose enters the
    #    absorption compartment `depot` (Aa in Fig. 1) and reaches the
    #    central compartment V1 with bioavailability F at first-order rate
    #    Kabs; an intravenous dose is given directly into `central`.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Bioavailability on the intramuscular absorption compartment
    f(depot) <- fdepot

    # 5. Observations. Cc is total plasma cephalothin (mg/L = ug/mL: dose in
    #    mg, volume in L). Cu is the free (unbound) concentration that drives
    #    the paper's fT > MIC PK/PD target and its probability of target
    #    attainment; it is a deterministic transform of Cc, is never measured
    #    directly, and so carries no residual error of its own.
    Cc <- central / vc
    Cu <- fu * Cc

    Cc ~ prop(propSd) + add(addSd)
  })
}

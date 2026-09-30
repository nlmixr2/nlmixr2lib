Yang_2020_enrofloxacin_rainbowTrout_pbpk <- function() {
  description <- paste(
    "Veterinary (rainbow trout). PBPK (whole-body, flow-limited, acslXtreme",
    "3.0.2.1) for enrofloxacin (ENR) and its metabolite ciprofloxacin (CIP)",
    "in rainbow trout, with water temperature driving cardiac output (Yang",
    "et al. 2020, Front Vet Sci 7:608348). The ENR sub-model has stomach",
    "contents, intestinal contents (the absorption depot), gill, liver, gut,",
    "kidney, muscle, skin, rest of body, venous and arterial blood, plus the",
    "culture water; the CIP sub-model has gill, liver, kidney, muscle, skin,",
    "rest of body, venous and arterial blood. All cardiac output passes",
    "through the gill between venous and arterial blood. Three routes are",
    "encoded simultaneously: oral (dose into the stomach, gastric emptying,",
    "first-pass absorption into the liver in competition with faecal loss",
    "to the water), intravenous (dose into venous blood) and immersion bath",
    "(bath drug infused into the water compartment and taken up across the",
    "gill into venous blood). ENR is cleared by renal excretion, by loss",
    "from the gill to the water and by hepatic conversion to CIP; CIP is",
    "cleared renally. Cardiac output is the linear function of water",
    "temperature CO = 3.95 * Temp - 12.9 mL/min/kg, so the model applies",
    "only above about 3.3 degC (the authors validated it at 5-17 degC).",
    "Deterministic structure with no residual-error model; the eight",
    "parameters of the authors' 500-iteration Monte Carlo withdrawal-interval",
    "analysis carry between-fish variability as normal distributions",
    "truncated at the mean +/- 1 SD (Supplementary Table S1)."
  )
  reference <- paste(
    "Yang F, Yang F, Wang D, Zhang C-S, Wang H, Song Z-W, Shao H-T, Zhang M,",
    "Yu M-L, Zheng Y. Development and application of a water temperature",
    "related physiologically based pharmacokinetic model for enrofloxacin",
    "and its metabolite ciprofloxacin in rainbow trout. Front Vet Sci.",
    "2020;7:608348. doi:10.3389/fvets.2020.608348.",
    "Equations transcribed from the acslXtreme code in Supplementary",
    "Material section 1 (Data Sheet 1); parameter values from Tables 3-4 and",
    "the same code; Monte Carlo distributions from Supplementary Table S1.",
    sep = " "
  )
  vignette <- "Yang_2020_enrofloxacin_rainbowTrout_pbpk"

  # Doses in mg (the paper prescribes mg/kg; multiply by WT). An immersion
  # bath is dosed into `a_water` as bath concentration (ppm = mg/L) times the
  # water volume `vwater`, infused over the bath duration. States hold
  # amounts in mg and volumes are in L (tissue density 1 kg/L), so
  # amount / volume is mg/L = ug/mL = ug/g, the paper's reporting units.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  # Three dosing routes, none of which is `depot` or `central`.
  dosing <- c("stomach", "a_venous", "a_water")

  # First-paper, fish-specific states. A canonical compartment needs a
  # second paper before it is registered, so these are declared here.
  paper_specific_compartments <- c(
    "a_water",
    "a_degraded",
    "a_gill",
    "a_gill_cipro"
  )

  compartmentData <- list(
    a_water = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    a_degraded = list(analyte = "enrofloxacin", units = "mg", specimen = "not applicable", verified = TRUE),
    stomach = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    a_gill = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_gut = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_metabolized = list(analyte = "enrofloxacin", units = "mg", specimen = "not applicable", verified = TRUE),
    a_kidney = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_skin = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_remainder = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_venous = list(analyte = "enrofloxacin", units = "mg", specimen = "whole blood", verified = TRUE),
    a_arterial = list(analyte = "enrofloxacin", units = "mg", specimen = "whole blood", verified = TRUE),
    a_liver_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "urine", verified = TRUE),
    a_muscle_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_skin_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_remainder_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_gill_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_venous_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "whole blood", verified = TRUE),
    a_arterial_cipro = list(analyte = "ciprofloxacin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      reference_value = 0.05,
      notes = paste(
        "Every tissue volume is a fixed fraction of body weight, cardiac",
        "output and both renal clearances are per-kg, and oral and IV doses",
        "are per-kg, so ENR and CIP concentrations after oral or IV dosing do",
        "not depend on WT. After an immersion bath they do: the amount taken",
        "up across the gill depends on the bath, not on the fish, so",
        "concentrations scale as 1 / WT. The code default is 0.05 kg; Yang",
        "2020 re-set it to the mean weight of each validation study (Table",
        "1: 0.1, 0.15, 0.204, 0.45 and 0.05 kg). reference_value is the",
        "code default."
      ),
      source_name = "bw"
    ),
    BODYTEMP = list(
      description = "Water temperature, taken as the body temperature of the fish",
      units = "degC",
      type = "continuous",
      reference_category = NULL,
      reference_value = 10,
      notes = paste(
        "Yang 2020 uses the rearing-water temperature `Temp`. Rainbow trout",
        "are heterothermic, so their body temperature is the water",
        "temperature, which is why this is carried as BODYTEMP. It enters",
        "only through cardiac output, CO (mL/min/kg) = 3.95 * Temp - 12.9",
        "(Yang 2020 Materials and Methods, citing their ref 42), and so",
        "scales every tissue blood flow. CO is zero at 3.27 degC and",
        "negative below it, so the model cannot be used below about 3.3",
        "degC. It was calibrated and validated at 5-17 degC. The relationship",
        "is not centred; reference_value = 10 degC (the middle of the",
        "withdrawal-interval temperatures 5, 10 and 16 degC) is only a",
        "default simulation value."
      ),
      source_name = "Temp"
    )
  )

  population <- list(
    species = "rainbow trout (Oncorhynchus mykiss)",
    n_subjects = NA_integer_,
    n_studies = 5L,
    age_range = NA_character_,
    weight_range = "about 50-450 g mean body weight across the five source studies (Yang 2020 Table 1)",
    sex_female_pct = NA_real_,
    disease_state = "healthy farmed fish",
    dose_range = "ENR IV 5 and 10 mg/kg; oral 5, 10, 30 and 50 mg/kg single dose; immersion bath 20 ppm for 2.5 h and 100 ppm for 0.5 h (Yang 2020 Table 1)",
    regions = "literature data (studies listed in Yang 2020 Table 1); model developed in China",
    water_temperature_range = "5-17 degC",
    notes = paste(
      "No new animal experiment. Yang 2020 digitised published plasma,",
      "serum and tissue concentrations of ENR and CIP from five rainbow",
      "trout studies (Table 1, refs 3, 13, 37, 38, 39), used some of the data",
      "to optimise the unknown parameters one route at a time by Nelder-Mead",
      "in acslXtreme (Table 4) and the rest to check the fits by eye, by",
      "linear regression and by mean absolute percentage error (Table 5).",
      "Tissue volumes and blood flows come from earlier fish PBPK models",
      "(rainbow trout, crucian carp and grass carp). The ENR partition",
      "coefficients were computed by the area method from one study (ref",
      "38); CIP partition coefficients for liver, kidney, muscle, skin and",
      "gut were taken from a human model. The fitted model was used to",
      "predict ENR + CIP withdrawal intervals for muscle and skin at 5, 10",
      "and 16 degC by 500-iteration Monte Carlo simulation (Table 6,",
      "Figure 8)."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Oral absorption (Yang 2020 Table 4, optimisation step 2, fitted to
    # the ref 3 plasma and tissue data; also the SI code).
    # ---------------------------------------------------------------
    lkst <- fixed(log(0.175))
    label("Gastric emptying rate constant Kst, stomach to intestinal contents (1/h)") # Table 4 (Kst final 0.175 /h); SI code `Kst=0.175`
    lka <- fixed(log(0.052))
    label("Oral absorption rate constant KaPO from intestinal contents to liver (1/h)") # Table 4 (KaPO final 0.052 /h); SI code `KaPO=0.052`
    lkfec <- fixed(log(0.605))
    label("Faecal elimination rate constant Kguc for unabsorbed ENR, intestinal contents to water (1/h)") # Table 4 (Kguc final 0.605 /h); SI code `Kguc=0.605`
    lfdepot <- fixed(log(0.6613))
    label("Oral bioavailability F (fraction)") # Materials and Methods 'the bioavailability was 66.13% (3)'; SI code `F=0.6613`; Table S1 (66.13%)

    # ---------------------------------------------------------------
    # Immersion bath and water (Table 4, optimisation step 3, fitted to
    # the ref 38 100 ppm / 0.5 h bath data; also the SI code).
    # ---------------------------------------------------------------
    lka_gill <- fixed(log(1.103))
    label("Branchial uptake rate constant KaIB, water to venous blood (1/h)") # Table 4 (KaIB final 1.103 /h); SI code `KaIB=1.103`
    lk_gill_water <- fixed(log(0.061))
    label("Rate constant Kgw for ENR loss from the gill to the water (1/h)") # Table 4 (Kgw final 0.061 /h); SI code `Kgw=0.061`
    lkdeg_water <- fixed(log(12003.31))
    label("Degradation rate constant Kde of ENR in the water (1/h)") # Table 4 (Kde final 12003.310 /h); SI code `Kde=12003.31`
    vwater <- fixed(40)
    label("Volume of culture water Vwater (L)") # SI code `constant Vwater=40 ! L`; the dose into the water is ppm * Vwater

    # ---------------------------------------------------------------
    # Elimination (Table 4, optimisation steps 1 and 4; also the SI code).
    # ---------------------------------------------------------------
    lkmet_cipro <- fixed(log(0.0725))
    label("Rate constant Kf for hepatic conversion of ENR to CIP (1/h)") # Table 4 (Kf final 0.0725 /h); SI code `Kf=0.0725`
    lcl_renal <- fixed(log(0.058))
    label("ENR renal clearance per kg body weight Clre (L/h/kg)") # Table 4 (Clre final 0.058 L/kg/h); SI code `Clreperkg=0.058`
    lcl_renal_cipro <- fixed(log(116.14))
    label("CIP renal clearance per kg body weight Clmre (L/h/kg)") # Table 4 (Clmre final 116.14 L/kg/h); SI code `Clmreperkg=116.14`

    # ---------------------------------------------------------------
    # ENR tissue:plasma partition coefficients (Table 3 column 'P x';
    # Pgi and Pr optimised in Table 4 step 1; the rest by the area
    # method). Also the SI code `constant Pgi=3.46,Pl=4.9,...`.
    # ---------------------------------------------------------------
    lkp_gill <- fixed(log(3.46))
    label("ENR gill partition coefficient Pgi (unitless)") # Table 3 (gill 3.46); Table 4 (Pgi final 3.46)
    lkp_liver <- fixed(log(4.9))
    label("ENR liver partition coefficient Pl (unitless)") # Table 3 (liver 4.9)
    lkp_kidney <- fixed(log(11.53))
    label("ENR kidney partition coefficient Pk (unitless)") # Table 3 (kidney 11.53)
    lkp_muscle <- fixed(log(2.83))
    label("ENR muscle partition coefficient Pm (unitless)") # Table 3 (muscle 2.83)
    lkp_skin <- fixed(log(7.38))
    label("ENR skin partition coefficient Ps (unitless)") # Table 3 (skin 7.38)
    lkp_gut <- fixed(log(4.88))
    label("ENR gut partition coefficient Pg (unitless)") # Table 3 (gut 4.88)
    lkp_remainder <- fixed(log(0.13))
    label("ENR rest-of-body partition coefficient Pr (unitless)") # Table 3 (rest 0.13); Table 4 (Pr final 0.13)

    # ---------------------------------------------------------------
    # CIP tissue:plasma partition coefficients (Table 3 column 'P mx';
    # Pmgi and Pmr optimised in Table 4 step 4, the rest taken from a
    # human model). Table 3 also lists Pmg = 3.39 for gut, but the CIP
    # sub-model has no gut compartment and the SI code never uses Pmg.
    # ---------------------------------------------------------------
    lkp_gill_cipro <- fixed(log(2.45))
    label("CIP gill partition coefficient Pmgi (unitless)") # Table 3 (gill 2.45); Table 4 (Pmgi final 2.45)
    lkp_liver_cipro <- fixed(log(3.67))
    label("CIP liver partition coefficient Pml (unitless)") # Table 3 (liver 3.67)
    lkp_kidney_cipro <- fixed(log(8.2))
    label("CIP kidney partition coefficient Pmk (unitless)") # Table 3 (kidney 8.2)
    lkp_muscle_cipro <- fixed(log(1.6))
    label("CIP muscle partition coefficient Pmm (unitless)") # Table 3 (muscle 1.6)
    lkp_skin_cipro <- fixed(log(0.718))
    label("CIP skin partition coefficient Pms (unitless)") # Table 3 (skin 0.718)
    lkp_remainder_cipro <- fixed(log(0.15))
    label("CIP rest-of-body partition coefficient Pmr (unitless)") # Table 3 (rest 0.15); Table 4 (Pmr final 0.15)

    # ---------------------------------------------------------------
    # The three physiological constants that the Monte Carlo analysis
    # varied (Table S1). Every other volume and flow fraction is a
    # constant in model().
    # ---------------------------------------------------------------
    fq_kidney <- fixed(0.056)
    label("Kidney blood flow Qck (fraction of cardiac output)") # Table 3 (kidney 0.056); SI code `Qck=0.056`; Table S1 (mean 0.056)
    fvol_muscle <- fixed(0.66)
    label("Muscle volume Vcm (fraction of body weight)") # Table 3 (muscle 0.66); SI code `Vcm=0.66`; Table S1 (mean 0.66)
    fvol_liver <- fixed(0.0126)
    label("Liver volume Vcl (fraction of body weight)") # Table 3 (liver 0.0126); SI code `Vcl=0.0126`. Table S1 prints mean 0.029, the liver BLOOD-FLOW fraction; see vignette Errata

    # ---------------------------------------------------------------
    # Between-fish variability of the Monte Carlo analysis (Yang 2020
    # 'Monte Carlo Analysis' and Table S1). Each parameter is normal,
    # with the mean equal to its model value and SD = 10% of the mean
    # (Kf: SD 0.02 /h), truncated at mean - SD and mean + SD. Each eta
    # below is a standard-normal deviate (variance fixed at 1) that
    # model() maps onto that truncated normal exactly, so these are NOT
    # log-scale variances. Water temperature was varied too, but its
    # distribution is not reported, so it is not encoded.
    # ---------------------------------------------------------------
    etafq_kidney ~ fixed(1) # Table S1 Qck: mean 0.056, SD 0.0056, min 0.0504, max 0.0616
    etafvol_muscle ~ fixed(1) # Table S1 Vcm: mean 0.66, SD 0.066, min 0.594, max 0.726
    etafvol_liver ~ fixed(1) # Table S1 Vcl: SD 10% of the mean (the table prints mean 0.029, SD 0.0029)
    etalkp_liver ~ fixed(1) # Table S1 Pl: mean 4.9, SD 0.49, min 4.41, max 5.39
    etalkp_muscle ~ fixed(1) # Table S1 Pm: mean 2.83, SD 0.283, min 2.547, max 3.113
    etalkp_muscle_cipro ~ fixed(1) # Table S1 Pmm: mean 1.6, SD 0.16, min 1.44, max 1.76
    etalfdepot ~ fixed(1) # Table S1 F: mean 66.13%, SD 6.613%, min 59.517%, max 72.743%
    etalkmet_cipro ~ fixed(1) # Table S1 Kf: mean 0.0725, SD 0.02, min 0.0525, max 0.0925
  })

  model({
    # =================================================================
    # Monte Carlo transform. `probitInv(eta, lo, hi)` maps a standard
    # normal deviate onto (lo, hi) through the normal CDF; with lo and
    # hi set to pnorm(-1) and pnorm(1), `probit()` of the result is a
    # standard normal truncated to (-1, 1). Parameter = mean + SD * z
    # then reproduces the Table S1 truncated normals exactly. With all
    # etas at zero, z = 0 and each parameter equals its model value.
    # =================================================================
    z_fq_kidney <- probit(probitInv(etafq_kidney, 0.158655253931457, 0.841344746068543))
    z_fvol_muscle <- probit(probitInv(etafvol_muscle, 0.158655253931457, 0.841344746068543))
    z_fvol_liver <- probit(probitInv(etafvol_liver, 0.158655253931457, 0.841344746068543))
    z_kp_liver <- probit(probitInv(etalkp_liver, 0.158655253931457, 0.841344746068543))
    z_kp_muscle <- probit(probitInv(etalkp_muscle, 0.158655253931457, 0.841344746068543))
    z_kp_muscle_cipro <- probit(probitInv(etalkp_muscle_cipro, 0.158655253931457, 0.841344746068543))
    z_fdepot <- probit(probitInv(etalfdepot, 0.158655253931457, 0.841344746068543))
    z_kmet_cipro <- probit(probitInv(etalkmet_cipro, 0.158655253931457, 0.841344746068543))

    # =================================================================
    # Individual parameters
    # =================================================================
    kst <- exp(lkst)
    ka <- exp(lka)
    kfec <- exp(lkfec)
    fdepot <- exp(lfdepot) * (1 + 0.1 * z_fdepot) # Table S1 F, SD 10% of mean
    ka_gill <- exp(lka_gill)
    k_gill_water <- exp(lk_gill_water)
    kdeg_water <- exp(lkdeg_water)
    kmet_cipro <- exp(lkmet_cipro) + 0.02 * z_kmet_cipro # Table S1 Kf, SD 0.02 /h
    cl_renal <- exp(lcl_renal) * WT # SI code `Clre=Clreperkg*bw` (L/h)
    cl_renal_cipro <- exp(lcl_renal_cipro) * WT # SI code `Clmre=Clmreperkg*bw` (L/h)

    kp_gill <- exp(lkp_gill)
    kp_liver <- exp(lkp_liver) * (1 + 0.1 * z_kp_liver) # Table S1 Pl, SD 10% of mean
    kp_kidney <- exp(lkp_kidney)
    kp_muscle <- exp(lkp_muscle) * (1 + 0.1 * z_kp_muscle) # Table S1 Pm, SD 10% of mean
    kp_skin <- exp(lkp_skin)
    kp_gut <- exp(lkp_gut)
    kp_remainder <- exp(lkp_remainder)

    kp_gill_cipro <- exp(lkp_gill_cipro)
    kp_liver_cipro <- exp(lkp_liver_cipro)
    kp_kidney_cipro <- exp(lkp_kidney_cipro)
    kp_muscle_cipro <- exp(lkp_muscle_cipro) * (1 + 0.1 * z_kp_muscle_cipro) # Table S1 Pmm, SD 10% of mean
    kp_skin_cipro <- exp(lkp_skin_cipro)
    kp_remainder_cipro <- exp(lkp_remainder_cipro)

    # =================================================================
    # Trout physiology (Table 3 and the SI code INITIAL section)
    # =================================================================
    # Cardiac output from water temperature, SI code
    # `CO=3.95*Temp-12.9` (mL/min/kg) and `Qtot=CO/1000*60*bw` (L/h).
    co <- 3.95 * BODYTEMP - 12.9
    qc <- co / 1000 * 60 * WT

    # Blood-flow fractions of cardiac output (Table 3). The rest of body
    # is the complement, SI code `Qcr=1-(Qck+Qcm+Qcs+Qcl+Qcg)`, which
    # evaluates to 0.1081 as printed in Table 3. The gill receives the
    # whole cardiac output (Qcgi = 1).
    fq_kidney_i <- fq_kidney * (1 + 0.1 * z_fq_kidney) # Table S1 Qck, SD 10% of mean
    fq_muscle <- 0.6 # Table 3 (muscle 0.6); SI code `Qcm=0.6`
    fq_skin <- 0.053 # Table 3 (skin 0.053); SI code `Qcs=0.053`
    fq_gut <- 0.1539 # Table 3 (gut 0.1539); SI code `Qcg=0.1539`
    fq_liver <- 0.029 # Table 3 (liver 0.029); SI code `Qcl=0.029`
    fq_remainder <- 1 - (fq_kidney_i + fq_muscle + fq_skin + fq_liver + fq_gut)

    q_kidney <- fq_kidney_i * qc
    q_muscle <- fq_muscle * qc
    q_skin <- fq_skin * qc
    q_gut <- fq_gut * qc
    q_liver <- fq_liver * qc
    q_remainder <- fq_remainder * qc

    # Volume fractions of body weight (Table 3). The rest of body is the
    # complement, SI code `Vcr=1-(Vcgi+Vck+Vcm+Vcs+Vcg+Vcl+Vcvb+Vcab)`,
    # which evaluates to 0.02079 as printed in Table 3. It is formed from
    # the TYPICAL muscle and liver fractions, so the Monte Carlo draws on
    # those two do not reach it: the complement is only 0.021, and a draw
    # of Vcm more than 0.32 SD above its mean would otherwise make it
    # negative and the rest-of-body state unstable. See vignette Errata.
    fvol_muscle_i <- fvol_muscle * (1 + 0.1 * z_fvol_muscle) # Table S1 Vcm, SD 10% of mean
    fvol_liver_i <- fvol_liver * (1 + 0.1 * z_fvol_liver) # Table S1 Vcl, SD 10% of mean
    fvol_gill <- 0.039 # Table 3 (gill 0.039); SI code `Vcgi=0.039`
    fvol_kidney <- 0.00841 # Table 3 (kidney 0.00841); SI code `Vck=0.00841`
    fvol_skin <- 0.1 # Table 3 (skin 0.1); SI code `Vcs=0.1`
    fvol_gut <- 0.0852 # Table 3 (gut 0.0852); SI code `Vcg=0.0852`
    fvol_venous <- 0.059 # Table 3 (venous blood 0.059); SI code `Vcvb=0.059`
    fvol_arterial <- 0.015 # Table 3 (arterial blood 0.015); SI code `Vcab=0.015`
    fvol_remainder <- 1 - (fvol_gill + fvol_kidney + fvol_muscle + fvol_skin +
      fvol_gut + fvol_liver + fvol_venous + fvol_arterial)

    v_gill <- fvol_gill * WT
    v_kidney <- fvol_kidney * WT
    v_muscle <- fvol_muscle_i * WT
    v_skin <- fvol_skin * WT
    v_gut <- fvol_gut * WT
    v_liver <- fvol_liver_i * WT
    v_venous <- fvol_venous * WT
    v_arterial <- fvol_arterial * WT
    v_remainder <- fvol_remainder * WT
    # The CIP sub-model has no gut compartment. SI code `Vmr=Vr+Vg`: the
    # CIP rest of body takes the gut tissue volume as well.
    v_remainder_cipro <- v_remainder + v_gut

    # Plasma and serum fractions of blood. SI code `Vvp=Vvb*(1-pcv)` with
    # haematocrit pcv = 0.304 and `Vvs=Vvb*(1-pre)` with pre = 0.555
    # (Materials and Methods: 'The pre and pcv values were set as 55.5
    # and 30.4%, respectively (42)'). Drug is assumed not to enter blood
    # cells, so blood drug amount / plasma (or serum) volume is the
    # plasma (or serum) concentration.
    hct <- 0.304
    f_serum <- 1 - 0.555

    # Molecular weights (g/mol), SI code `MWENR=359.4, MWCIP=331.347`.
    mw_enr <- 359.4
    mw_cipro <- 331.347

    # =================================================================
    # Concentrations (mg/L). SI code DERIVATIVE section.
    # =================================================================
    c_gill <- a_gill / v_gill
    c_gut <- a_gut / v_gut
    c_liver <- a_liver / v_liver
    c_kidney <- a_kidney / v_kidney
    c_muscle <- a_muscle / v_muscle
    c_skin <- a_skin / v_skin
    c_remainder <- a_remainder / v_remainder
    c_venous <- a_venous / v_venous
    c_arterial <- a_arterial / v_arterial

    c_gill_cipro <- a_gill_cipro / v_gill
    c_liver_cipro <- a_liver_cipro / v_liver
    c_kidney_cipro <- a_kidney_cipro / v_kidney
    c_muscle_cipro <- a_muscle_cipro / v_muscle
    c_skin_cipro <- a_skin_cipro / v_skin
    c_remainder_cipro <- a_remainder_cipro / v_remainder_cipro
    c_venous_cipro <- a_venous_cipro / v_venous
    c_arterial_cipro <- a_arterial_cipro / v_arterial

    # =================================================================
    # Fluxes (mg/h)
    # =================================================================
    r_absorb <- fdepot * ka * depot # SI code `F*KaPO*aic`
    r_feces <- kfec * depot * (1 - fdepot) # SI code `Reli1=Kguc*aic*(1-F)`
    r_met <- kmet_cipro * a_liver # SI code `Reli2=Kf*al`
    r_renal <- cl_renal * c_kidney / kp_kidney # SI code `Reli3=Clre*Ck/Pk`
    r_gill_water <- k_gill_water * a_gill # SI code `Kgw*agi`
    r_uptake <- ka_gill * a_water # SI code `KaIB*aw`
    r_degrade <- kdeg_water * a_water # SI code `Rade=Kde*aw`
    r_formation <- r_met / mw_enr * mw_cipro # SI code `Rbiotrans=Kf*al/MWENR*MWCIP`
    r_renal_cipro <- cl_renal_cipro * c_kidney_cipro / kp_kidney_cipro # SI code `relim=Clmre*Cmk/Pmk`

    # =================================================================
    # ENR ODEs (SI code; Yang 2020 Table 2 prints a simplified version,
    # see the vignette Errata). The immersion-bath dose is infused into
    # `a_water`, the oral dose is given into `stomach` and the IV dose
    # into `a_venous`.
    # =================================================================
    # Water. Receives faecal, renal and branchial ENR and loses ENR to
    # uptake across the gill and to degradation. SI code
    # `raw=-KaIB*aw+Kgw*agi+Reli1+Reli3+Rdoseib-Rade`.
    d/dt(a_water) <- -r_uptake + r_gill_water + r_feces + r_renal - r_degrade
    d/dt(a_degraded) <- r_degrade # SI code `MASSOUT=integ(Rade,0)+outbio`

    d/dt(stomach) <- -kst * stomach # SI code `rastc=Rdosepo-Kst*astc`
    d/dt(depot) <- kst * stomach - r_absorb - r_feces # SI code `raic=Kst*astc-F*KaPO*aic-Reli1`

    # All cardiac output passes through the gill. SI code
    # `ragi=Qtot*(Cvb-Cgi/Pgi)-Kgw*agi`.
    d/dt(a_gill) <- qc * (c_venous - c_gill / kp_gill) - r_gill_water

    d/dt(a_gut) <- q_gut * (c_arterial - c_gut / kp_gut) # SI code `rag=Qg*(Cab-Cg/Pg)`

    # Liver: hepatic-artery inflow, gut venous inflow, absorbed oral dose;
    # outflow at (Ql + Qg). SI code
    # `ral=F*KaPO*aic+Ql*Cab+Qg*Cg/Pg-(Ql+Qg)*Cl/Pl-Reli2`.
    d/dt(a_liver) <- r_absorb + q_liver * c_arterial + q_gut * c_gut / kp_gut -
      (q_liver + q_gut) * c_liver / kp_liver - r_met
    d/dt(a_metabolized) <- r_met # SI code `outbio=integ(Reli2,0)`

    d/dt(a_kidney) <- q_kidney * (c_arterial - c_kidney / kp_kidney) - r_renal # SI code `rak=Qk*(Cab-Ck/Pk)-Reli3`
    d/dt(a_muscle) <- q_muscle * (c_arterial - c_muscle / kp_muscle) # SI code `ram=Qm*(Cab-Cm/Pm)`
    d/dt(a_skin) <- q_skin * (c_arterial - c_skin / kp_skin) # SI code `ras=Qs*(Cab-Cs/Ps)`
    d/dt(a_remainder) <- q_remainder * (c_arterial - c_remainder / kp_remainder) # SI code `rar=Qr*(Cab-Cr/Pr)`

    # SI code `ravb=Rdoseiv+KaIB*aw+Qr*Cr/Pr+Qs*Cs/Ps+Qm*Cm/Pm+Qk*Ck/Pk+
    # (Ql+Qg)*Cl/Pl-Qtot*Cvb`. Bath uptake enters venous blood.
    d/dt(a_venous) <- r_uptake + q_remainder * c_remainder / kp_remainder +
      q_skin * c_skin / kp_skin + q_muscle * c_muscle / kp_muscle +
      q_kidney * c_kidney / kp_kidney + (q_liver + q_gut) * c_liver / kp_liver -
      qc * c_venous
    d/dt(a_arterial) <- qc * (c_gill / kp_gill - c_arterial) # SI code `raab=Qtot*(Cgi/Pgi-Cab)`

    # =================================================================
    # CIP ODEs (SI code). CIP is formed in the liver at the ENR
    # conversion rate, converted from mg ENR to mg CIP by molecular
    # weight, and cleared renally. There is no CIP gut compartment:
    # gut blood flow reaches the liver straight from arterial blood.
    # =================================================================
    d/dt(a_liver_cipro) <- r_formation + (q_liver + q_gut) *
      (c_arterial_cipro - c_liver_cipro / kp_liver_cipro) # SI code `raml=Rbiotrans+(Ql+Qg)*(Cmab-Cml/Pml)`
    d/dt(a_kidney_cipro) <- q_kidney * (c_arterial_cipro - c_kidney_cipro / kp_kidney_cipro) -
      r_renal_cipro # SI code `ramk=Qk*(Cmab-Cmk/Pmk)-relim`
    d/dt(a_urine_cipro) <- r_renal_cipro # SI code `mout=integ(relim,0)`
    d/dt(a_muscle_cipro) <- q_muscle * (c_arterial_cipro - c_muscle_cipro / kp_muscle_cipro) # SI code `ramm`
    d/dt(a_skin_cipro) <- q_skin * (c_arterial_cipro - c_skin_cipro / kp_skin_cipro) # SI code `rams`
    d/dt(a_remainder_cipro) <- q_remainder * (c_arterial_cipro - c_remainder_cipro / kp_remainder_cipro) # SI code `ramr`, with `Cmr=amr/Vmr`
    d/dt(a_gill_cipro) <- qc * (c_venous_cipro - c_gill_cipro / kp_gill_cipro) # SI code `ramgi=Qtot*(Cmvb-Cmgi/Pmgi)`
    d/dt(a_venous_cipro) <- q_remainder * c_remainder_cipro / kp_remainder_cipro +
      q_skin * c_skin_cipro / kp_skin_cipro + q_muscle * c_muscle_cipro / kp_muscle_cipro +
      q_kidney * c_kidney_cipro / kp_kidney_cipro +
      (q_liver + q_gut) * c_liver_cipro / kp_liver_cipro - qc * c_venous_cipro # SI code `ramvb`
    d/dt(a_arterial_cipro) <- qc * (c_gill_cipro / kp_gill_cipro - c_arterial_cipro) # SI code `ramab=Qtot*(Cmgi/Pmgi-Cmab)`

    # =================================================================
    # Outputs (ug/mL for plasma and serum, ug/g for tissues). The SI code
    # forms plasma and serum from VENOUS blood (`Cvp=avb/Vvp`,
    # `Cvs=avb/Vvs`), which is what Yang 2020 compared with the data.
    # =================================================================
    Cc <- a_venous / (v_venous * (1 - hct)) # ENR plasma; SI code `Cvp`
    Cserum <- a_venous / (v_venous * f_serum) # ENR serum; SI code `Cvs`
    Cliver <- c_liver
    Ckidney <- c_kidney
    Cmuscle <- c_muscle
    Cskin <- c_skin
    Cgut <- c_gut
    Cgill <- c_gill
    Cwater <- a_water / vwater # SI code `Cwater=aw/Vwater`

    Cc_cipro <- a_venous_cipro / (v_venous * (1 - hct)) # CIP plasma; SI code `Cmvp`
    Cserum_cipro <- a_venous_cipro / (v_venous * f_serum) # CIP serum; SI code `Cmvs`
    Cliver_cipro <- c_liver_cipro
    Ckidney_cipro <- c_kidney_cipro
    Cmuscle_cipro <- c_muscle_cipro
    Cskin_cipro <- c_skin_cipro
    Cgill_cipro <- c_gill_cipro

    # Marker residue for the withdrawal interval: ENR + CIP in muscle and
    # in skin (Yang 2020 'Withdrawal Interval Estimation'; MRL 0.1 ug/g).
    Cmuscle_total <- Cmuscle + Cmuscle_cipro
    Cskin_total <- Cskin + Cskin_cipro
  })
}

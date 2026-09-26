Garessus_2019_isoniazid_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, coupled mother-infant dyad, 26 ODEs, R deSolve lsoda).",
    "Isoniazid transfer from a lactating mother to her breastfed newborn",
    "(Garessus 2019, Front Pharmacol). Two complete flow-limited whole-body",
    "models run side by side and are coupled only through breast milk: the",
    "mother carries ten perfusion-limited well-stirred tissue compartments",
    "plus blood, an oral depot, and a breast-milk compartment perfused at",
    "breast-tissue blood flow; the infant carries the same ten tissues plus",
    "blood and a depot, with no milk compartment of its own. At every feed the",
    "mother's milk compartment is flushed into the infant's depot, which is",
    "the infant's only route of exposure. Oral drug is delivered from the",
    "depot straight into the liver, so hepatic extraction acts as first pass;",
    "the modelled bioavailability that results is 86 percent (fast) or 93",
    "percent (slow). The spleen drains portally into the liver rather than",
    "into blood. Elimination is a single hepatic clearance applied to the",
    "TOTAL liver concentration and deliberately NOT divided by the liver",
    "partition coefficient, because the source clearances are apparent values",
    "fitted to plasma that never accounted for unbound fraction; this is why",
    "the published local sensitivity analysis finds plasma AUC equally",
    "sensitive to clearance and to kp_liver. NAT2 acetylator status is",
    "modelled as two discrete clearance values per person, selected",
    "independently for mother and infant by NAT2_SLOW and NAT2_SLOW_INFANT,",
    "so all four dyad phenotype pairs can be simulated. Physiology is fixed",
    "ICRP 2002 reference data for an adult woman and a newborn, not scaled by",
    "body weight. Deterministic simulation model: the paper's reported",
    "intervals come from the confidence interval of the two clearance",
    "estimates, not from a random-effects model, so there is no",
    "interindividual variability and propSd is fixed at 0."
  )
  reference <- paste(
    "Garessus EDG, Mielke H, Gundert-Remy U. Exposure of infants to isoniazid",
    "via breast milk after maternal drug intake of recommended doses is",
    "clinically insignificant irrespective of metaboliser status.",
    "A physiologically-based pharmacokinetic (PBPK) modelling approach to",
    "estimate drug exposure of infants via breast-feeding.",
    "Front Pharmacol. 2019;10:5. doi:10.3389/fphar.2019.00005.",
    "Equation 1 and Table 1 give the flow-limited distribution equation and",
    "every physiological and drug-specific constant. The complete deSolve",
    "model script is Supplementary Data Sheet 1",
    "('430130_resubmitted_Garessus_Supplementary_Material_RScript.R'), which",
    "is the authoritative source for the ODE topology, the placement of the",
    "clearance term, and the breastfeeding gate.",
    "Clearances are inherited from Wilkins JJ et al. Br J Clin Pharmacol.",
    "2011;72(1):51-62 (mother) and Rey E et al. Fundam Clin Pharmacol.",
    "2001;15:355-359 (infant); ka from Wilkins 2011; the milk:plasma",
    "partition coefficient from Singh N et al. Br J Clin Pharmacol.",
    "2008;65(3):418-422."
  )
  vignette <- "Garessus_2019_isoniazid_breast_milk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: every state below was traced to the
  # named compartment of the Supplementary Data Sheet 1 model script
  # (mother A1-A13, infant A1C-A13C).
  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    lung = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    milk = list(analyte = "isoniazid", units = "mg", specimen = "milk", verified = TRUE),
    blood = list(analyte = "isoniazid", units = "mg", specimen = "whole blood", verified = TRUE),
    infant_depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    infant_lung = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_brain = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_heart = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_liver = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_spleen = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_kidney = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_adipose = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_skin = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_muscle = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_bone = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE),
    infant_blood = list(analyte = "isoniazid", units = "mg", specimen = "whole blood", verified = TRUE),
    infant_a_oral = list(analyte = "isoniazid", units = "mg", specimen = "not applicable", verified = TRUE),
    a_metabolized = list(analyte = "isoniazid", units = "mg", specimen = "not applicable", verified = TRUE),
    infant_a_metabolized = list(analyte = "isoniazid", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    NAT2_SLOW = list(
      description = "NAT2 slow-acetylator phenotype of the MOTHER (1 = slow metaboliser, 0 = fast)",
      units = "binary",
      type = "categorical",
      reference_category = "0 (fast metaboliser)",
      source_name = "metaboliser status (mother)",
      notes = paste(
        "Selects the mother's hepatic clearance between the two discrete",
        "values fitted by Wilkins 2011 in South African tuberculosis patients",
        "and quoted in the Metabolism and Excretion section: 21.6 L/h for",
        "fast metabolisers and 9.7 L/h for slow metabolisers. The paper",
        "simulates only these two cases; it fits no continuous NAT2 covariate",
        "relationship and reports no intermediate-acetylator value, so the",
        "companion NAT2_RAPID canonical is deliberately not carried here.",
        "Set NAT2_SLOW = 1 to simulate a slow-metabolising mother."
      )
    ),
    NAT2_SLOW_INFANT = list(
      description = "NAT2 slow-acetylator phenotype of the breastfed INFANT (1 = slow, 0 = fast)",
      units = "binary",
      type = "categorical",
      reference_category = "0 (fast metaboliser)",
      source_name = "metaboliser status (child)",
      notes = paste(
        "Dyad-partner covariate: describes the breastfed infant, not the",
        "modelled mother, and so takes the _INFANT suffix of the mother-infant",
        "dyad partner namespace alongside WT_INFANT and AGE_INFANT. Selects",
        "the infant's hepatic clearance between the two values observed by",
        "Rey 2001 in children under 6 months weighing 4 kg: 2.55 L/h (fast)",
        "and 0.76 L/h (slow). Mother and infant phenotypes are independent, so",
        "the four dyad pairs of Figures 4 and 5 are simulated by crossing this",
        "covariate with NAT2_SLOW."
      )
    )
  )

  # Screened / described by the paper but never entered into the model: the
  # infant physiology is fixed ICRP newborn reference data, NOT scaled by
  # body weight, so the 4 kg infant weight appears only when the authors
  # convert the absolute infant dose (mg/day) to a per-kg dose (mg/kg/day)
  # in the Table 2 and Table 3 footnotes.
  covariatesDataExcluded <- list(
    WT_INFANT = list(
      description = "Breastfed-infant body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Fixed at 4 kg throughout (Metabolism and Excretion: the Rey 2001",
        "clearances are 'observed in 3 children aged less than 6 months ...",
        "in infants weighing 4 kg'). Used only to report the external infant",
        "dose per kilogram in the Table 2 / Table 3 footnotes c and d; no",
        "volume, flow or clearance in the model is scaled by it."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 0,
    disease_state = paste(
      "Simulated lactating women receiving isoniazid monotherapy for",
      "drug-susceptible tuberculosis, and their exclusively breastfed",
      "newborns. No subject was dosed for this analysis: the model is a",
      "forward simulation built entirely from published physiology and from",
      "clearance and absorption estimates fitted elsewhere."
    ),
    mother = paste(
      "ICRP 2002 reference adult woman, 20-50 years. Organ masses converted",
      "to volumes at a density of 1 kg/L; cardiac output 5.9 L/min = 354 L/h.",
      "Breast-milk volume fixed at 0.1134 L (Kent 2006). The ICRP data",
      "describe a mainly United States population and the authors note it is",
      "leaner-than-average female tuberculosis patients they intend to",
      "represent."
    ),
    infant_partner = paste(
      "ICRP 2002 reference newborn weighing 4 kg; cardiac output 0.6 L/min =",
      "36 L/h. Organ blood flows use the adult female percentages of cardiac",
      "output applied to the newborn cardiac output, as the authors state",
      "explicitly. Exclusively breastfed, feeding every 2 h."
    ),
    dose_range = paste(
      "Maternal oral isoniazid 300 mg once daily (simulation 1) or 900 mg",
      "every 3 days (simulation 2), the two highest doses recommended for",
      "breastfeeding mothers by CDC 2016 and Nahid 2016. Validation",
      "simulations additionally used 200 mg (Lass and Bunger 1953) and a",
      "40 mg direct infant dose (10 mg/kg, Rey 2001)."
    ),
    feeding_pattern = paste(
      "Breastfeeding every 2 h, the first feed 2 h after the maternal dose",
      "(12 feeds per 24 h). The authors chose this because breastfeeding",
      "mothers are advised to take the drug immediately after a feed."
    ),
    regions = "Model physiology from ICRP 2002; maternal clearances from a South African tuberculosis cohort (Wilkins 2011); infant clearances from a French paediatric cohort (Rey 2001).",
    notes = paste(
      "Validation was visual only, against four clinical studies of maternal",
      "plasma and breast-milk isoniazid (Lass and Bunger 1953, Ricci and",
      "Copaitich 1954, Berlin and Lee 1979, Singh 2007) and against the",
      "paediatric plasma data of Rey 2001; no goodness-of-fit statistic is",
      "reported. The intervals throughout Tables 2 and 3 are produced by",
      "re-running the deterministic model at the lower limit, mean and upper",
      "limit of the clearance 95% confidence intervals, so they are clearance",
      "uncertainty bands and NOT population prediction intervals."
    )
  )

  ini({
    # ---- Absorption ------------------------------------------------------
    # Not estimated here; taken from the Wilkins 2011 popPK fit and assumed
    # identical in mother and infant 'in the absence of infant-specific data'.
    lka <- fixed(log(1.82))
    label("First-order absorption rate constant ka, mother and infant (1/h)")
    # Absorption section: 'a first-order absorption rate constant (ka) of
    # 1.82 [/h] that was clinically determined (Wilikins et al., 2011)';
    # supplementary R script 'ka <- 1.82'

    # ---- Hepatic clearance, mother --------------------------------------
    # Apparent clearances fitted in adults by Wilkins 2011. Applied to the
    # TOTAL liver concentration, see the model() block.
    lcl_fast <- fixed(log(21.6))
    label("Hepatic clearance of the mother, NAT2 fast metaboliser (L/h)")
    # Metabolism and Excretion: 'for the adult population were: 21.6
    # (18.9-28.2) L/h ... for fast'; script 'CL_fast <- c(18.9,21.6,28.2)'
    lcl_slow <- fixed(log(9.7))
    label("Hepatic clearance of the mother, NAT2 slow metaboliser (L/h)")
    # Metabolism and Excretion: '9.7 (9.61-10.7) L/h, for ... slow
    # metabolisers'; script 'CL_slow <- c(9.61,9.7,10.7)'

    # ---- Hepatic clearance, breastfed infant -----------------------------
    lcl_fast_infant <- fixed(log(2.55))
    label("Hepatic clearance of the breastfed infant, NAT2 fast metaboliser (L/h)")
    # Metabolism and Excretion: 'In the infant population, the clearance
    # values observed were: 2.55 (0.77-4.32) L/h ... in infants weighing
    # 4 kg' (Rey 2001); script 'CL_fast_c <- c(0.77,2.55,4.32)'
    lcl_slow_infant <- fixed(log(0.76))
    label("Hepatic clearance of the breastfed infant, NAT2 slow metaboliser (L/h)")
    # Metabolism and Excretion: '0.76 (0.69-0.82) L/h' (Rey 2001);
    # script 'CL_slow_c <- c(0.69,0.76,0.82)'

    # ---- Breastfeeding pattern and milk-to-infant transfer ---------------
    feed_n <- fixed(12)
    label("Number of breastfeeds per 24 h, giving a 24/12 = 2 h feeding cycle (feeds/day)")
    # Breast Milk-Drinking Behaviour of the Infant: 'milk drinking was
    # modelled to take place every 2 h'; Table 2 footnote a 'n = 12;
    # drinking every 2 h in 24 h'
    feed_first <- fixed(2)
    label("Time of the first breastfeed after the maternal dose (h)")
    # Breast Milk-Drinking Behaviour of the Infant: 'first drug intake would
    # thus take place about 2 h after oral dosing of the mother'; script
    # dosing at t = 8 h and first milk pulse at t = 10 h
    feed_window <- fixed(0.01)
    label("Duration of each breastfeeding window (h; 0.01 h = 36 s)")
    # Supplementary R script 'duration_milkintake <- 0.01  #36sec'
    kmilkinf <- fixed(100)
    label("First-order transfer rate from the maternal milk compartment to the infant depot during a feed (1/h)")
    # Supplementary R script: dMilk_intake = A12 * gate / duration_milkintake,
    # i.e. a rate constant of 1/0.01 = 100 /h while the gate is open

    # ---- Residual error --------------------------------------------------
    # The paper reports no residual-error model and no interindividual
    # variability: its intervals are the clearance confidence intervals
    # propagated through a deterministic solve. Fixing the residual SDs at 0
    # keeps this a forward simulation rather than inventing a variance the
    # source does not report. See the vignette Assumptions and deviations.
    propSd <- fixed(0)
    label("Proportional residual error, maternal plasma (fraction; not reported by the source)")
    propSd_Cmilk <- fixed(0)
    label("Proportional residual error, breast milk (fraction; not reported by the source)")
    propSd_Cinfant <- fixed(0)
    label("Proportional residual error, infant plasma (fraction; not reported by the source)")
  })

  model({
    # =====================================================================
    # Physiological constants. ICRP 2002 reference values as tabulated in
    # Table 1 and set literally in Supplementary Data Sheet 1. They are
    # fixed data, not parameters: nothing here is scaled by body weight, so
    # the model describes one reference woman and one reference newborn.
    # =====================================================================

    # ---- Mother: volumes (L), Table 1 column 'V [L] Mother' -------------
    # Organ masses (g) converted at an assumed density of 1 kg/L.
    v_lung <- 0.42 # script 'Vlung <- 420*0.001'
    v_brain <- 1.3 # script 'Vbrain <- 1300*0.001'
    v_heart <- 0.25 # script 'Vheart <- 250*0.001'
    v_liver <- 1.4 # script 'Vliver <- 1400*0.001'
    v_spleen <- 0.13 # script 'Vspleen <- 130*0.001'
    v_kidney <- 0.275 # script 'Vkidneys <- 275*0.001'
    v_adipose <- 22.5 # script 'Vadiposetissue <- 22500*0.001'
    v_skin <- 2.3 # script 'Vskin <- 2300*0.001'
    v_muscle <- 17.5 # script 'Vskeletalmuscle <- 17500*0.001'
    v_bone <- 7.8 # script 'Vbone <- 7800*0.001'
    v_blood <- 3.9 # script 'Vblood <- 3900*0.001'
    v_milk <- 0.1134 # Table 1 'Breast milk 0.1134'; Kent 2006

    # ---- Mother: blood flows (L/h) as fractions of cardiac output -------
    # Table 1 column 'Q [L/h] Mother'. The percentages are the ICRP adult
    # female distribution; cardiac output 5.9 L/min = 354 L/h.
    cardiac_output <- 5.9 * 60 # script 'CardiacOutput <- 5.9*60' = 354 L/h
    q_lung <- 0.025 * cardiac_output # 8.85; Table 1 note a, bronchial arteries/veins
    q_brain <- 0.120 * cardiac_output # 42.48
    q_heart <- 0.050 * cardiac_output # 17.70; Table 1 note b, coronary flow.
    # The script defines Qheart <- CardiacOutput but never uses it: the heart
    # TISSUE compartment is perfused by Qcoronary, which is what Table 1 lists.
    q_spleen <- 0.030 * cardiac_output # 10.62
    q_kidney <- 0.170 * cardiac_output # 60.18
    q_adipose <- 0.085 * cardiac_output # 30.09
    q_skin <- 0.050 * cardiac_output # 17.70
    q_muscle <- 0.120 * cardiac_output # 42.48
    q_bone <- 0.050 * cardiac_output # 17.70
    q_milk <- 0.004 * cardiac_output # 1.416; 'assuming breast milk blood flow
    # to be identical to breast tissue blood flow'; script 'Qbreast'
    # Dual hepatic supply: arterial 20.5% + portal 6.5% = 27% leaves the
    # liver, but the splanchnic share that first passes through the spleen
    # (3%) arrives as spleen efflux rather than from blood.
    q_liver_in <- (0.205 + 0.065 - 0.030) * cardiac_output # 84.96; script 'Qliver_nospleen'
    q_liver_out <- 0.270 * cardiac_output # 95.58; Table 1 'Liver 95.58'; script 'Qliver_total'

    # ---- Infant: volumes (L), Table 1 column 'V [L] Infant' -------------
    v_lung_infant <- 0.03 # script 'Vlung_c <- 30*0.001'
    v_brain_infant <- 0.38 # script 'Vbrain_c <- 380*0.001'
    v_heart_infant <- 0.02 # script 'Vheart_c <- 20*0.001'
    v_liver_infant <- 0.13 # script 'Vliver_c <- 130*0.001'
    v_spleen_infant <- 0.0095 # script 'Vspleen_c <- 9.5*0.001'
    v_kidney_infant <- 0.025 # script 'Vkidneys_c <- 25*0.001'
    v_adipose_infant <- 0.93 # script 'Vadiposetissue_c <- 930*0.001'
    v_skin_infant <- 0.175 # script 'Vskin_c <- 175*0.001'
    v_muscle_infant <- 0.8 # script 'Vskeletalmuscle_c <- 800*0.001'
    v_bone_infant <- 0.37 # script 'Vbone_c <- 370*0.001'
    v_blood_infant <- 0.27 # script 'Vblood_c <- 270*0.001'

    # ---- Infant: blood flows (L/h), Table 1 column 'Q [L/h] Infant' -----
    # Same percentages of cardiac output as the mother, applied to the
    # newborn cardiac output of 0.6 L/min = 36 L/h.
    cardiac_output_infant <- 0.6 * 60 # script 'CardiacOutput_c <- 0.6*60' = 36 L/h
    q_lung_infant <- 0.025 * cardiac_output_infant # 0.90
    q_brain_infant <- 0.120 * cardiac_output_infant # 4.32
    q_heart_infant <- 0.050 * cardiac_output_infant # 1.80; coronary
    q_spleen_infant <- 0.030 * cardiac_output_infant # 1.08
    q_kidney_infant <- 0.170 * cardiac_output_infant # 6.12
    q_adipose_infant <- 0.085 * cardiac_output_infant # 3.06
    q_skin_infant <- 0.050 * cardiac_output_infant # 1.80
    q_muscle_infant <- 0.120 * cardiac_output_infant # 4.32
    q_bone_infant <- 0.050 * cardiac_output_infant # 1.80
    q_liver_in_infant <- (0.205 + 0.065 - 0.030) * cardiac_output_infant # 8.64
    q_liver_out_infant <- 0.270 * cardiac_output_infant # 9.72; Table 1 'Liver 9.72'

    # ---- Tissue:blood partition coefficients, Table 1 column 'PC' -------
    # Calculated by the Schmitt 2008 method from pKa 1.82, logKow, fu 0.9 and
    # a plasma water fraction of 0.935; shared by mother and infant, which is
    # how the supplementary script uses them (one PC object, both models).
    kp_lung <- 0.79
    kp_brain <- 0.73
    kp_heart <- 0.70
    kp_liver <- 0.70
    kp_spleen <- 0.75
    kp_kidney <- 0.74
    kp_adipose <- 0.15
    kp_skin <- 0.62
    kp_muscle <- 0.74
    kp_bone <- 0.33
    # Measured, not calculated: Schmitt's data collection has no breast-tissue
    # composition, so the authors used the average milk:plasma AUC ratio of
    # Singh 2007 instead.
    kp_milk <- 0.89

    # =====================================================================
    # Individual parameters
    # =====================================================================
    ka <- exp(lka)
    # Two discrete clearances per person, selected by acetylator phenotype.
    cl <- exp(lcl_fast * (1 - NAT2_SLOW) + lcl_slow * NAT2_SLOW)
    cl_infant <- exp(
      lcl_fast_infant * (1 - NAT2_SLOW_INFANT) +
        lcl_slow_infant * NAT2_SLOW_INFANT
    )

    # =====================================================================
    # Breastfeeding gate
    # ---------------------------------------------------------------------
    # The source script hardcodes ~90 hyperbolic-tangent pulses, one per
    # feed, each a smoothed rectangle of width feed_window opening at
    # feed_first + k * feed_cycle. That is reproduced here as a single
    # periodic pulse so the schedule is set by feed_n / feed_first rather
    # than by the length of a literal sum.
    #
    # feed_phase is built so the pulse sits at the MIDDLE of the cycle; the
    # modulo wrap is then half a cycle away from the pulse and cannot clip
    # its tanh tails. feed_on suppresses the spurious pulses that the
    # periodic construction would otherwise place before the first feed, and
    # its switch is likewise placed half a cycle before it.
    # =====================================================================
    feed_cycle <- 24 / feed_n
    feed_rel <- t - feed_first + feed_cycle / 2
    feed_phase <- feed_rel - feed_cycle * floor(feed_rel / feed_cycle)
    feed_pulse <- (tanh(100 * (feed_phase - feed_cycle / 2)) -
      tanh(100 * (feed_phase - feed_cycle / 2 - feed_window))) / 2
    # Steepness 100 and the /2 are the source script's `anaus` function
    # verbatim: (tanh(100*(t-t0)) - tanh(100*(t-t1)))/2.
    feed_on <- (1 + tanh(10 * (t - feed_first + feed_cycle / 2))) / 2
    sqwMilkToInfant <- feed_pulse * feed_on

    # =====================================================================
    # Concentrations
    # =====================================================================
    Cblood <- blood / v_blood
    Cblood_infant <- infant_blood / v_blood_infant
    # Flow-limited well-stirred efflux concentrations, i.e. the concentration
    # of blood leaving each organ: Eq. 1's (A_organ/V_organ)/PC_organ:blood.
    cv_lung <- (lung / v_lung) / kp_lung
    cv_brain <- (brain / v_brain) / kp_brain
    cv_heart <- (heart / v_heart) / kp_heart
    cv_liver <- (liver / v_liver) / kp_liver
    cv_spleen <- (spleen / v_spleen) / kp_spleen
    cv_kidney <- (kidney / v_kidney) / kp_kidney
    cv_adipose <- (adipose / v_adipose) / kp_adipose
    cv_skin <- (skin / v_skin) / kp_skin
    cv_muscle <- (muscle / v_muscle) / kp_muscle
    cv_bone <- (bone / v_bone) / kp_bone
    cv_milk <- (milk / v_milk) / kp_milk

    cv_lung_infant <- (infant_lung / v_lung_infant) / kp_lung
    cv_brain_infant <- (infant_brain / v_brain_infant) / kp_brain
    cv_heart_infant <- (infant_heart / v_heart_infant) / kp_heart
    cv_liver_infant <- (infant_liver / v_liver_infant) / kp_liver
    cv_spleen_infant <- (infant_spleen / v_spleen_infant) / kp_spleen
    cv_kidney_infant <- (infant_kidney / v_kidney_infant) / kp_kidney
    cv_adipose_infant <- (infant_adipose / v_adipose_infant) / kp_adipose
    cv_skin_infant <- (infant_skin / v_skin_infant) / kp_skin
    cv_muscle_infant <- (infant_muscle / v_muscle_infant) / kp_muscle
    cv_bone_infant <- (infant_bone / v_bone_infant) / kp_bone

    # Amount leaving the milk compartment into the infant during a feed.
    milk_to_infant <- kmilkinf * milk * sqwMilkToInfant

    # =====================================================================
    # Mother. Compartment order fixes the rxode2 slot numbering, so depot is
    # slot 1 and an oral dose needs no cmt bookkeeping beyond cmt = 'depot'.
    # =====================================================================
    d/dt(depot) <- -ka * depot # script 'dA1 <- -A1*ka + dMed_intake'; the
    # maternal intake pulse is replaced by a native rxode2 dose event, which
    # is the exact instantaneous limit of the script's 0.36-second pulse.

    d/dt(lung) <- Cblood * q_lung - cv_lung * q_lung
    d/dt(brain) <- Cblood * q_brain - cv_brain * q_brain
    d/dt(heart) <- Cblood * q_heart - cv_heart * q_heart
    d/dt(spleen) <- Cblood * q_spleen - cv_spleen * q_spleen
    d/dt(kidney) <- Cblood * q_kidney - cv_kidney * q_kidney
    d/dt(adipose) <- Cblood * q_adipose - cv_adipose * q_adipose
    d/dt(skin) <- Cblood * q_skin - cv_skin * q_skin
    d/dt(muscle) <- Cblood * q_muscle - cv_muscle * q_muscle
    d/dt(bone) <- Cblood * q_bone - cv_bone * q_bone

    # Liver. Four terms: the whole absorbed oral dose (so hepatic extraction
    # is the first pass), portal inflow carrying the spleen's efflux,
    # arterial + portal inflow from blood, and total efflux back to blood.
    #
    # The clearance term multiplies the TOTAL liver concentration
    # liver/v_liver and is NOT divided by kp_liver. This is the source
    # script verbatim -- 'dA5 <- ... - (CL_fast)*(A5/Vliver)' -- and it is
    # deliberate: 'this model did not restrict clearance to unbound drug
    # molecules, as clearance estimates were based on apparent measurements
    # that did not account for unbound drug fractions'. It is also what makes
    # plasma AUC scale as 1/(cl * kp_liver), which is exactly what the
    # published local sensitivity analysis reports (a 1% rise in clearance
    # drops AUC by 0.9988%, a 1% rise in kp_liver by 0.9987%). Dividing by
    # kp_liver here would make AUC = dose/cl and destroy that result.
    d/dt(liver) <- ka * depot +
      cv_spleen * q_spleen +
      Cblood * q_liver_in -
      cv_liver * q_liver_out -
      cl * (liver / v_liver)

    # Breast milk: perfused at breast-tissue blood flow, drained by feeding.
    d/dt(milk) <- Cblood * q_milk - cv_milk * q_milk - milk_to_infant

    # Blood collects every organ's efflux EXCEPT the spleen's, which the
    # liver already collected portally, and delivers the matching inflows.
    d/dt(blood) <- cv_lung * q_lung +
      cv_brain * q_brain +
      cv_heart * q_heart +
      cv_liver * q_liver_out +
      cv_kidney * q_kidney +
      cv_adipose * q_adipose +
      cv_skin * q_skin +
      cv_muscle * q_muscle +
      cv_bone * q_bone +
      cv_milk * q_milk -
      Cblood * (q_lung + q_brain + q_heart + q_liver_in + q_spleen +
        q_kidney + q_adipose + q_skin + q_muscle + q_bone + q_milk)

    # =====================================================================
    # Breastfed infant. Structurally identical, minus the milk compartment;
    # its depot is filled by feeding instead of by a dose event.
    # =====================================================================
    d/dt(infant_depot) <- -ka * infant_depot + milk_to_infant

    d/dt(infant_lung) <- Cblood_infant * q_lung_infant - cv_lung_infant * q_lung_infant
    d/dt(infant_brain) <- Cblood_infant * q_brain_infant - cv_brain_infant * q_brain_infant
    d/dt(infant_heart) <- Cblood_infant * q_heart_infant - cv_heart_infant * q_heart_infant
    d/dt(infant_spleen) <- Cblood_infant * q_spleen_infant - cv_spleen_infant * q_spleen_infant
    d/dt(infant_kidney) <- Cblood_infant * q_kidney_infant - cv_kidney_infant * q_kidney_infant
    d/dt(infant_adipose) <- Cblood_infant * q_adipose_infant - cv_adipose_infant * q_adipose_infant
    d/dt(infant_skin) <- Cblood_infant * q_skin_infant - cv_skin_infant * q_skin_infant
    d/dt(infant_muscle) <- Cblood_infant * q_muscle_infant - cv_muscle_infant * q_muscle_infant
    d/dt(infant_bone) <- Cblood_infant * q_bone_infant - cv_bone_infant * q_bone_infant

    d/dt(infant_liver) <- ka * infant_depot +
      cv_spleen_infant * q_spleen_infant +
      Cblood_infant * q_liver_in_infant -
      cv_liver_infant * q_liver_out_infant -
      cl_infant * (infant_liver / v_liver_infant)

    d/dt(infant_blood) <- cv_lung_infant * q_lung_infant +
      cv_brain_infant * q_brain_infant +
      cv_heart_infant * q_heart_infant +
      cv_liver_infant * q_liver_out_infant +
      cv_kidney_infant * q_kidney_infant +
      cv_adipose_infant * q_adipose_infant +
      cv_skin_infant * q_skin_infant +
      cv_muscle_infant * q_muscle_infant +
      cv_bone_infant * q_bone_infant -
      Cblood_infant * (q_lung_infant + q_brain_infant + q_heart_infant +
        q_liver_in_infant + q_spleen_infant + q_kidney_infant +
        q_adipose_infant + q_skin_infant + q_muscle_infant + q_bone_infant)

    # ---- Process accumulators -------------------------------------------
    # Cumulative isoniazid actually ingested by the infant via milk. This is
    # the mechanistically consistent external infant dose. It is NOT the
    # number the paper headlines as 'oral dose [mg/d]': that figure is a
    # separate bookkeeping sum of the milk AMOUNT standing in the compartment
    # at each of the 12 feed times, which assumes each feed empties the
    # compartment completely. Because milk keeps being perfused from blood
    # during the 36-second feed window while the compartment is draining, the
    # two differ; see the vignette Assumptions and deviations.
    d/dt(infant_a_oral) <- milk_to_infant

    # Cumulative hepatic elimination. The source script carries these as
    # 'dMetab_m' and 'dMetab_c' for exactly one reason: they are what closes
    # the mass balance, which the script then checks with its MASS term.
    d/dt(a_metabolized) <- cl * (liver / v_liver) # script 'dMetab_m'
    d/dt(infant_a_metabolized) <- cl_infant * (infant_liver / v_liver_infant) # script 'dMetab_c'

    # =====================================================================
    # Observations
    # =====================================================================
    Cc <- Cblood # maternal plasma, Figure 2 panel 1 / Table 2 'Mother Plasma'
    Cmilk <- cv_milk * kp_milk # = milk / v_milk, Table 2 'Breast milk'
    Cinfant <- Cblood_infant # infant plasma, Table 2 'Infant Plasma'

    Cc ~ prop(propSd)
    Cmilk ~ prop(propSd_Cmilk)
    Cinfant ~ prop(propSd_Cinfant)
  })
}

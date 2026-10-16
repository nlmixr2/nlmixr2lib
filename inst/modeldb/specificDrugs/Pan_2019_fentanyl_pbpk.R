Pan_2019_fentanyl_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, linear, perfusion-limited). Original 16-state fentanyl",
    "PBPK system of Bjorkman et al. 1994 for a 73 kg reference man, as coded",
    "by Pan and Duffull 2019 to serve as the full model for automated proper",
    "lumping. Venous and arterial blood, eight single-compartment tissues",
    "(brain, heart, kidney, muscle, skin, adipose, carcass as `other`), a",
    "two-compartment liver (`liver`, `liver_deep`) that is the only site of",
    "elimination, a two-compartment pancreas + spleen tissue (`spleen`,",
    "`spleen_deep`; volumes, flows and partition coefficients are the means of",
    "the pancreas and spleen values) and a three-compartment gut (`gut`,",
    "`gut_deep`, `gut_deep2`). Gut and pancreas + spleen drain into the liver.",
    "Every state is a CONCENTRATION (ug/L) and the system is the authors'",
    "deposited micro-rate-constant matrix K transcribed term by term. No",
    "between-subject variability and no residual error: a deterministic",
    "system. The two 4-state lumped reductions of this system are",
    "Pan_2019_fentanyl_lumped_a and Pan_2019_fentanyl_lumped_b."
  )
  reference <- paste(
    "Pan S, Duffull SB (2019). Automated proper lumping for simplification of",
    "linear physiologically based pharmacokinetic systems. J Pharmacokinet",
    "Pharmacodyn 46(4):361-370. doi:10.1007/s10928-019-09644-5.",
    "Parameter values and the 16 x 16 micro-rate-constant matrix from",
    "Supplementary Material 1 Sections 1A-1B and the deposited MATLAB files",
    "(Model_parameter_values.m, Runfile_Autolumping_FenPBPK.m). Physiology",
    "originally from Bjorkman S, Stanski DR, Verotta D, Harashima H (1994)",
    "Comparative tissue concentration profiles of fentanyl and alfentanil in",
    "humans predicted from tissue/blood partition data obtained in rats.",
    "Anesthesiology 80:380-393 (cited as ref 15)."
  )
  vignette <- "Pan_2019_fentanyl_lumping"
  units <- list(time = "min", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    venous = list(analyte = "fentanyl", units = "ug/L", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "fentanyl", units = "ug/L", specimen = "whole blood", verified = TRUE),
    brain = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    other = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    gut_deep = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    gut_deep2 = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    spleen_deep = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE),
    liver_deep = list(analyte = "fentanyl", units = "ug/L", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1L,
    age_range = "adult (reference man)",
    weight_range = "73 kg (reference man)",
    disease_state = "Not applicable: a physiology-based reference individual, not a fitted cohort.",
    dose_range = "750 ug intravenous bolus (Pan 2019); 750 ug over 5 min in the original Bjorkman 1994 study",
    notes = paste(
      "A single deterministic reference individual: Supplementary Material 1",
      "Section 1A gives WT = 73 kg ('average weight for man') and a cardiac",
      "output of 6.7 L/min. Organ volumes, blood flows and tissue:blood",
      "partition coefficients are those Bjorkman 1994 derived from rat",
      "tissue/blood partition data scaled to man. Tissue weights in kg are",
      "used as volumes in L (unit density)."
    )
  )

  ini({
    # All values: Pan 2019 Supplementary Material 1 Section 1A ('Adopted from
    # Bjorkman et al. 1994'), identical to the deposited Model_parameter_values.m.
    # Organ blood flows are printed as flow per kg of tissue times the tissue
    # volume (e.g. 'Qbra = 0.56*Vbra'), so the per-kg flow is carried here and
    # multiplied out in model(). The lung parameters of Section 1A are not
    # carried: the deposited K matrix sets the lung concentration equal to the
    # venous concentration and never uses them.
    qc <- fixed(6.7); label("Cardiac output (L/min)") # S1A 'Qcad = 6.7'
    v_venous <- fixed(3.9); label("Venous blood volume (L)") # S1A 'Vven = 3.9'
    v_arterial <- fixed(2.0); label("Arterial blood volume (L)") # S1A 'Vart = 2.0'

    v_brain <- fixed(1.4); label("Brain volume (L)") # S1A 'Vbra = 1.4'
    qn_brain <- fixed(0.56); label("Brain blood flow per kg tissue (L/min/kg)") # S1A 'Qbra = 0.56*Vbra'
    lkp_brain <- fixed(log(4)); label("Brain:blood partition coefficient (unitless)") # S1A 'KPbra = 4'

    v_heart <- fixed(0.33); label("Heart volume (L)") # S1A 'Vhea = 0.33'
    qn_heart <- fixed(0.8); label("Heart blood flow per kg tissue (L/min/kg)") # S1A 'Qhea = 0.8*Vhea'
    lkp_heart <- fixed(log(5.1)); label("Heart:blood partition coefficient (unitless)") # S1A 'KPhea = 5.1'

    v_kidney <- fixed(0.31); label("Kidney volume (L)") # S1A 'Vkid = 0.31'
    qn_kidney <- fixed(4); label("Kidney blood flow per kg tissue (L/min/kg)") # S1A 'Qkid = 4*Vkid'
    lkp_kidney <- fixed(log(13.7)); label("Kidney:blood partition coefficient (unitless)") # S1A 'KPkid = 13.7'

    v_liver <- fixed(1.8); label("Liver volume (L)") # S1A 'Vliv = 1.8'
    qn_liver <- fixed(0.96); label("Liver blood flow per kg tissue (L/min/kg)") # S1A 'Qliv = 0.96*Vliv'
    lkp_liver <- fixed(log(4.3)); label("Liver:blood partition coefficient (unitless)") # S1A 'KPliv = 4.3'

    v_muscle <- fixed(30); label("Muscle volume (L)") # S1A 'Vmus = 30'
    qn_muscle <- fixed(0.038); label("Muscle blood flow per kg tissue (L/min/kg)") # S1A 'Qmus = 0.038*Vmus'
    lkp_muscle <- fixed(log(3.5)); label("Muscle:blood partition coefficient (unitless)") # S1A 'KPmus = 3.5'

    v_skin <- fixed(3.3); label("Skin volume (L)") # S1A 'Vski = 3.3'
    qn_skin <- fixed(0.12); label("Skin blood flow per kg tissue (L/min/kg)") # S1A 'Qski = 0.12*Vski'
    lkp_skin <- fixed(log(2.3)); label("Skin:blood partition coefficient (unitless)") # S1A 'KPski = 2.3'

    v_adipose <- fixed(12.5); label("Adipose volume (L)") # S1A 'Vfat = 12.5'
    qn_adipose <- fixed(0.028); label("Adipose blood flow per kg tissue (L/min/kg)") # S1A 'Qfat = 0.028*Vfat'
    lkp_adipose <- fixed(log(30.3)); label("Adipose:blood partition coefficient (unitless)") # S1A 'KPfat = 30.3'

    v_other <- fixed(15.5); label("Carcass volume (L)") # S1A 'Vcar = 15.5'
    qn_other <- fixed(0.04); label("Carcass blood flow per kg tissue (L/min/kg)") # S1A 'Qcar = 0.04*Vcar'
    lkp_other <- fixed(log(1)); label("Carcass:blood partition coefficient, assumed (unitless)") # S1A 'KPcar = 1; (not available)'

    v_gut <- fixed(1.2); label("Gut volume, also the first gut compartment volume (L)") # S1A 'Vgut = 1.2', 'Vgut_1 = Vgut'
    qn_gut <- fixed(0.83); label("Gut blood flow per kg tissue (L/min/kg)") # S1A 'Qgut = 0.83*Vgut'
    lkp_gut <- fixed(log(9)); label("Gut:blood partition coefficient (unitless)") # S1A 'KPgut = 9'

    # Pancreas + spleen: Bjorkman 1994 pancreas and spleen values averaged.
    v_spleen <- fixed(0.14); label("Pancreas + spleen volume, mean of the two organs (L)") # S1A 'Vpas = (0.18+0.10)/2'
    qn_spleen <- fixed(0.9); label("Pancreas + spleen blood flow per kg tissue, mean of the two organs (L/min/kg)") # S1A 'Qpas = (0.6+1.2)/2*Vpas'
    lkp_spleen <- fixed(log(27.7)); label("Pancreas + spleen:blood partition coefficient, mean of the two organs (unitless)") # S1A 'KPpas = (24.1+31.3)/2'

    # Intra-organ sub-model constants, per kg of the organ (Bjorkman 1994 p 388
    # reports them in mL/min per g and mL per g; Section 1A converts to L/min
    # and L by multiplying by the organ volume in kg).
    cln_spleen_deep <- fixed(2.21); label("Pancreas + spleen shallow-deep exchange clearance per kg tissue (L/min/kg)") # S1A 'CL12_pas = 2.21*Vpas'
    vn_spleen <- fixed(0.082); label("Pancreas + spleen first (perfused) compartment volume per kg tissue (L/kg)") # S1A 'Vpas_1 = 0.082*Vpas'
    vn_spleen_deep <- fixed(22.3); label("Pancreas + spleen deep compartment volume per kg tissue (L/kg)") # S1A 'Vpas_2 = 22.3*Vpas'
    cln_liver_deep <- fixed(0.35); label("Liver shallow-deep exchange clearance per kg tissue (L/min/kg)") # S1A 'CL12_hep = 0.35*Vliv'
    cln_liver <- fixed(0.095); label("Hepatic elimination clearance per kg liver (L/min/kg)") # S1A 'CL10_hep = 0.095*Vliv' (95 mL/min per kg)
    vn_liver_deep <- fixed(8.7); label("Liver deep compartment volume per kg tissue (L/kg)") # S1A 'Vhep_2 = 8.7*Vliv'
    cln_gut_deep <- fixed(0.89); label("Gut first-to-second compartment exchange clearance per kg tissue (L/min/kg)") # S1A 'CL12_gut = 0.89*Vgut'
    cln_gut_deep2 <- fixed(0.03); label("Gut second-to-third compartment exchange clearance per kg tissue (L/min/kg)") # S1A 'CL23_gut = 0.03*Vgut'
    vn_gut_deep <- fixed(6.4); label("Gut second compartment volume per kg tissue (L/kg)") # S1A 'Vgut_2 = 6.4*Vgut'
    vn_gut_deep2 <- fixed(3.26); label("Gut third compartment volume per kg tissue (L/kg)") # S1A 'Vgut_3 = 3.26*Vgut'
  })

  model({
    kp_brain <- exp(lkp_brain)
    kp_heart <- exp(lkp_heart)
    kp_kidney <- exp(lkp_kidney)
    kp_liver <- exp(lkp_liver)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_adipose <- exp(lkp_adipose)
    kp_other <- exp(lkp_other)
    kp_gut <- exp(lkp_gut)
    kp_spleen <- exp(lkp_spleen)

    # Blood flows (L/min) and tissue efflux rate constants KT = Q / (V * KP)
    # (1/min), as in Section 1A.
    q_brain <- qn_brain * v_brain
    q_heart <- qn_heart * v_heart
    q_kidney <- qn_kidney * v_kidney
    q_liver <- qn_liver * v_liver
    q_muscle <- qn_muscle * v_muscle
    q_skin <- qn_skin * v_skin
    q_adipose <- qn_adipose * v_adipose
    q_other <- qn_other * v_other
    q_gut <- qn_gut * v_gut
    q_spleen <- qn_spleen * v_spleen

    kt_brain <- q_brain / (v_brain * kp_brain)
    kt_heart <- q_heart / (v_heart * kp_heart)
    kt_kidney <- q_kidney / (v_kidney * kp_kidney)
    kt_liver <- q_liver / (v_liver * kp_liver)
    kt_muscle <- q_muscle / (v_muscle * kp_muscle)
    kt_skin <- q_skin / (v_skin * kp_skin)
    kt_adipose <- q_adipose / (v_adipose * kp_adipose)
    kt_other <- q_other / (v_other * kp_other)
    kt_gut <- q_gut / (v_gut * kp_gut)
    kt_spleen <- q_spleen / (v_spleen * kp_spleen)

    # Intra-organ sub-models (L/min and L).
    cl_spleen_deep <- cln_spleen_deep * v_spleen
    v_spleen1 <- vn_spleen * v_spleen
    v_spleen_deep <- vn_spleen_deep * v_spleen
    cl_liver_deep <- cln_liver_deep * v_liver
    cl_liver <- cln_liver * v_liver
    v_liver_deep <- vn_liver_deep * v_liver
    cl_gut_deep <- cln_gut_deep * v_gut
    cl_gut_deep2 <- cln_gut_deep2 * v_gut
    v_gut_deep <- vn_gut_deep * v_gut
    v_gut_deep2 <- vn_gut_deep2 * v_gut

    # ODEs: rows 1-16 of the Section 1B matrix K, in its state order. States
    # are concentrations (ug/L). Lung is not a state: the deposited code sets
    # Clun = Cven. Gut and pancreas + spleen drain only to the liver, and the
    # liver receives no arterial input in K.
    d/dt(venous) <- (-qc * venous +
      kt_brain * v_brain * brain + kt_heart * v_heart * heart +
      kt_kidney * v_kidney * kidney + kt_liver * v_liver * liver +
      kt_muscle * v_muscle * muscle + kt_skin * v_skin * skin +
      kt_adipose * v_adipose * adipose + kt_other * v_other * other) / v_venous
    d/dt(arterial) <- qc / v_arterial * venous - qc / v_arterial * arterial
    d/dt(brain) <- q_brain / v_brain * arterial - kt_brain * brain
    d/dt(heart) <- q_heart / v_heart * arterial - kt_heart * heart
    d/dt(kidney) <- q_kidney / v_kidney * arterial - kt_kidney * kidney
    d/dt(liver) <- (-(cl_liver_deep + cl_liver) / v_liver - kt_liver) * liver +
      kt_gut * v_gut / v_liver * gut + kt_spleen * v_spleen / v_liver * spleen +
      cl_liver_deep / v_liver * liver_deep
    d/dt(muscle) <- q_muscle / v_muscle * arterial - kt_muscle * muscle
    d/dt(skin) <- q_skin / v_skin * arterial - kt_skin * skin
    d/dt(adipose) <- q_adipose / v_adipose * arterial - kt_adipose * adipose
    d/dt(other) <- q_other / v_other * arterial - kt_other * other
    d/dt(gut) <- q_gut / v_gut * arterial - (q_gut + cl_gut_deep) / v_gut * gut +
      cl_gut_deep / v_gut * gut_deep
    d/dt(gut_deep) <- cl_gut_deep / v_gut_deep * gut -
      (cl_gut_deep + cl_gut_deep2) / v_gut_deep * gut_deep +
      cl_gut_deep2 / v_gut_deep * gut_deep2
    d/dt(gut_deep2) <- cl_gut_deep2 / v_gut_deep2 * gut_deep - cl_gut_deep2 / v_gut_deep2 * gut_deep2
    d/dt(spleen) <- q_spleen / v_spleen1 * arterial - (q_spleen + cl_spleen_deep) / v_spleen1 * spleen +
      cl_spleen_deep / v_spleen1 * spleen_deep
    d/dt(spleen_deep) <- cl_spleen_deep / v_spleen_deep * spleen - cl_spleen_deep / v_spleen_deep * spleen_deep
    d/dt(liver_deep) <- cl_liver_deep / v_liver_deep * liver - cl_liver_deep / v_liver_deep * liver_deep

    # A bolus amount (ug) into venous enters as concentration Dose / Vven
    # (Runfile 'IC = [Dose/Vven; zeros(15,1)]').
    f(venous) <- 1 / v_venous

    # Observation: arterial blood concentration (Runfile 'output = [2]').
    Cc <- arterial
  })
}

Carter_2020_pitavastatin_eltrombopag_pbpk <- function() {
  description <- paste(
    "PBPK (semimechanistic, liver-focused; deSolve). Pitavastatin (victim) and",
    "eltrombopag (perpetrator) after oral dosing in simulated healthy adults,",
    "with mutual competitive inhibition of saturable hepatic uptake. Each drug",
    "has six states: stomach and gastrointestinal tract (amounts), and liver",
    "extracellular space, liver, gallbladder and central compartment",
    "(concentrations; whole blood for pitavastatin, plasma for eltrombopag).",
    "Unbound drug in the liver extracellular space enters the liver by",
    "Michaelis-Menten uptake, whose Km is multiplied by (1 + C_other / Km_other)",
    "with C_other the other drug's extracellular concentration, plus passive",
    "diffusion; unbound liver drug leaves by passive efflux, metabolism and",
    "biliary excretion into the gallbladder, which empties into the GI tract",
    "(enterohepatic recycling). The central compartment exchanges with the",
    "extracellular space through hepatic blood flow; pitavastatin also has a",
    "urinary clearance. Uptake and passive constants are scaled from the",
    "in vitro micro-rate-constant model (Carter_2020_pitavastatin_eltrombopag_invitro).",
    "Values and structure follow the authors' deSolve code. Between-subject",
    "variability is the authors' Monte-Carlo sampling: every parameter is",
    "normally distributed with a 10-30 percent CV; no residual error."
  )
  reference <- paste(
    "Carter SJ, Chouhan B, Sharma P, Chappell MJ. Prediction of Clinical",
    "Transporter-Mediated Drug-Drug Interactions via Comeasurement of",
    "Pitavastatin and Eltrombopag in Human Hepatocyte Models. CPT",
    "Pharmacometrics Syst Pharmacol. 2020;9(4):211-221.",
    "doi:10.1002/psp4.12505 (PMID 32142598, PMC7179958).",
    "ODEs: Supplemental Table S3 (Eqs. S41-S51) and the authors' deSolve",
    "script 'PBPK_Model.r' (Model Codes supplement), which is followed where",
    "the two differ. Parameter values: Table 2 ('PBPK model values') and the",
    "script's parameter block."
  )
  vignette <- "Carter_2020_pitavastatin_eltrombopag"
  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  dosing <- c("stomach", "stomach_eltrombopag")

  # Pitavastatin (the victim) keeps the bare state names and eltrombopag (the
  # perpetrator) carries the `_eltrombopag` suffix. The liver is
  # permeability-limited: `is_liver` is the liver extracellular space (V_ext.H)
  # and `int_liver` the liver cells (V_L). The liver, extracellular,
  # gallbladder and central states hold concentrations (ng/mL); the stomach and
  # gut-lumen states hold amounts (ng).
  paper_specific_compartments <- c(
    "stomach_eltrombopag",
    "gut_lumen_eltrombopag",
    "is_liver_eltrombopag",
    "int_liver_eltrombopag",
    "gallbladder_eltrombopag",
    "central_eltrombopag"
  )

  compartmentData <- list(
    stomach = list(analyte = "pitavastatin", units = "ng", specimen = "administration site", verified = TRUE),
    gut_lumen = list(analyte = "pitavastatin", units = "ng", specimen = "administration site", verified = TRUE),
    is_liver = list(analyte = "pitavastatin", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver = list(analyte = "pitavastatin", units = "ng/mL", specimen = "tissue", verified = TRUE),
    gallbladder = list(analyte = "pitavastatin", units = "ng/mL", specimen = "bile", verified = TRUE),
    central = list(analyte = "pitavastatin", units = "ng/mL", specimen = "whole blood", verified = TRUE),
    stomach_eltrombopag = list(
      analyte = "eltrombopag",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    gut_lumen_eltrombopag = list(
      analyte = "eltrombopag",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    is_liver_eltrombopag = list(analyte = "eltrombopag", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver_eltrombopag = list(analyte = "eltrombopag", units = "ng/mL", specimen = "tissue", verified = TRUE),
    gallbladder_eltrombopag = list(analyte = "eltrombopag", units = "ng/mL", specimen = "bile", verified = TRUE),
    central_eltrombopag = list(analyte = "eltrombopag", units = "ng/mL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 0L,
    disease_state = "healthy volunteers (simulated)",
    weight_range = "83 kg reference adult male (main text); Table S3 title states 70 kg",
    dose_range = "Single oral 1 mg pitavastatin and/or 75 mg eltrombopag",
    notes = paste(
      "Not fitted to clinical data: the model is built bottom-up from the",
      "in vitro hepatocyte estimates (scaled to a whole liver with 139e6",
      "hepatocytes per g liver and a 1695.6 g liver) and literature",
      "physiology, then compared with published clinical profiles - 1 mg",
      "pitavastatin from Prueksaritanont 2014 (n = 8, converted to blood",
      "with blood:plasma 0.425) and 75 mg eltrombopag from Deng 2011. The",
      "paper's simulations used 100 Monte-Carlo subjects. Eltrombopag",
      "Fa*Fg and the eltrombopag fraction unbound (0.005) were chosen by",
      "visual fit to the clinical eltrombopag data."
    )
  )

  ini({
    # Absorption and enterohepatic recycling (shared by both drugs)
    lkst <- fixed(log(0.1)); label("Gastric emptying rate constant, k_ge (1/min)") # Table 2, k_ge 0.1 (inverse of a 10-min gastric emptying time); script ktransP_v mean 0.1
    lkehc <- fixed(log(0.0618)); label("Gallbladder emptying (bile flow) rate constant, k_bile (1/min)") # Table 2 physiological values, k_bile 0.062; script BiTrans_v mean 0.0618

    # Pitavastatin
    lka <- fixed(log(0.07)); label("Pitavastatin absorption rate constant from the GI tract, K_a.P (1/min)") # Table 2, K_a.X for pitavastatin 0.07 (Eq. 5); script kaP_v mean 0.07
    lvmax <- fixed(log(29998636)); label("Pitavastatin whole-liver maximum uptake rate, V_max.up.P (ng/min)") # Table 2, V_max.up.X for pitavastatin 3.0e7 (Eq. 6); script VmP_v mean 29998636
    lkm <- fixed(log(9525)); label("Pitavastatin uptake Km, also its inhibition constant against eltrombopag uptake, K_m.up.P (ng/mL)") # Table 2, K_m.up.X for pitavastatin 9525 (Eq. 7); script KmP_v mean 9525
    lps_in <- fixed(log(117.8)); label("Pitavastatin passive diffusion clearance into the liver, P_dif.P (mL/min)") # Table 2, P_dif.X for pitavastatin 118 (Eq. 8); script PdPi_v mean 117.8
    lps_out <- fixed(log(104.2)); label("Pitavastatin passive diffusion clearance out of the liver, P_def.P (mL/min)") # Table 2, P_def.X for pitavastatin 104 (Eq. 9); script PdePi_v mean 104.2
    lclint_met <- fixed(log(93)); label("Pitavastatin hepatic metabolic clearance of unbound liver drug, CL_met.P (mL/min)") # Table 2, CL_met.X for pitavastatin 93 (Eq. 10); script CL_MePi_v mean 93
    lclint_bile <- fixed(log(165)); label("Pitavastatin biliary clearance of unbound liver drug, CL_bi.P (mL/min)") # Table 2, CL_bi.X for pitavastatin 165 (Eq. 12); script CL_BiPi_v mean 165
    lclurine <- fixed(log(11.5)); label("Pitavastatin urinary clearance parameter, CL_ur.P (mL/min)") # Table 2, CL_ur.P 12 (Eq. 14); script CLur_v mean 11.5
    lfu_blood <- fixed(log(0.009)); label("Pitavastatin fraction unbound in blood, fu_bl.P (fraction)") # Table 2, f_u.bl.P 0.009 (fu,pl 0.004 / Bl:Pl 0.425, Eq. 15); script Cfpp_v mean 0.009
    lfu_liver <- fixed(log(0.02)); label("Pitavastatin fraction unbound in the liver, fu_L.P (fraction)") # script fTP_v mean 0.02 ('fu, tissue from RED expt 2'); Table 2 prints 0.026 (Eq. 16)
    lvc <- fixed(log(5820)); label("Pitavastatin central (total blood) volume, V_c.P (mL)") # Table 2, V_c.X for pitavastatin 5820 (total blood volume); script VcP_v mean 5820

    # Eltrombopag
    lka_eltrombopag <- fixed(log(11)); label("Eltrombopag absorption rate constant from the GI tract, K_a.E (1/min)") # Table 2, K_a.X for eltrombopag 11 (footnote a: visual fit to tmax); script kaE_v mean 11
    lvmax_eltrombopag <- fixed(log(5104539)); label("Eltrombopag whole-liver maximum uptake rate, V_max.up.E (ng/min)") # Table 2, V_max.up.X for eltrombopag 5.1e6 (Eq. 6); script VmE_v mean 5104539
    lkm_eltrombopag <- fixed(log(3138)); label("Eltrombopag uptake Km, also its inhibition constant against pitavastatin uptake, K_m.up.E (ng/mL)") # Table 2, K_m.up.X for eltrombopag 3138 (Eq. 7); script KmE_v mean 3138
    lps_in_eltrombopag <- fixed(log(11784)); label("Eltrombopag passive diffusion clearance into the liver, P_dif.E (mL/min)") # script PdE_v mean 11784 (0.05 * 139 * 1695.6, Eq. 8); Table 2 prints 1.2e5
    lps_out_eltrombopag <- fixed(log(307)); label("Eltrombopag passive diffusion clearance out of the liver, P_def.E (mL/min)") # Table 2, P_def.X for eltrombopag 307 (Eq. 9); script PdeE_v mean 307
    lclint_met_eltrombopag <- fixed(log(6.8)); label("Eltrombopag hepatic metabolic clearance of unbound liver drug, CL_met.E (mL/min)") # Table 2, CL_met.X for eltrombopag 6.8 (Eq. 11); script CL_MeE_v mean 6.8
    lclint_bile_eltrombopag <- fixed(log(13)); label("Eltrombopag biliary clearance of unbound liver drug, CL_bi.E (mL/min)") # Table 2, CL_bi.X for eltrombopag 13 (Eq. 13); script CL_BiE_v mean 13
    lfu_eltrombopag <- fixed(log(0.005)); label("Eltrombopag fraction unbound in plasma, fu_pl.E (fraction)") # Table 2, f_u.pl.E 0.005 (best visual fit of 0.01, 0.002 and 0.005); script CfpE_v mean 0.005
    lfu_liver_eltrombopag <- fixed(log(0.4)); label("Eltrombopag fraction unbound in the liver, fu_L.E (fraction)") # Table 2, f_u.L.E 0.4 (footnote b: visual fit to data); script fTE_v mean 0.4
    lvc_eltrombopag <- fixed(log(2940)); label("Eltrombopag central (plasma) volume, V_c.E (mL)") # Table 2, V_c.X for eltrombopag 2940; script VcE_v mean 2940

    # Physiology
    lv_liver_ext <- fixed(log(556)); label("Liver extracellular space volume, V_ext.H (mL)") # script Vext_H_v mean 556 (Watanabe 2009); Table 2 prints 469
    lv_liver <- fixed(log(1570)); label("Liver volume, V_L (mL)") # Table 2 physiological values, V_L 1570; script VH_v mean 1570
    lv_gallbladder <- fixed(log(36)); label("Gallbladder volume, V_GaBl (mL)") # Table 2 physiological values, V_GaBl 36; script VGaBl_v mean 36
    lq_liver <- fixed(log(1320)); label("Hepatic blood flow, Q_H (mL/min)") # Table 2 physiological values, Q_H 1320; script Qh_v mean 1320
    lq_kidney <- fixed(log(1170)); label("Kidney blood flow, Q_K (mL/min)") # Table 2 physiological values, Q_K 1170; script Qk_v mean 1170

    # Monte-Carlo between-subject variability of the script: each parameter
    # is drawn as theta * (1 + eta), eta normal with variance CV^2.
    etalkst ~ fixed(0.04) # script ktransP_v sd 0.2 * mean
    etalkehc ~ fixed(0.04026) # script BiTrans_v sd 0.2 * 0.062 with mean 0.0618, CV 0.2006
    etalka ~ fixed(0.04) # script kaP_v sd 0.2 * mean
    etalvmax ~ fixed(0.09) # script VmP_v sd 0.3 * mean
    etalkm ~ fixed(0.09) # script KmP_v sd 0.3 * mean
    etalps_in ~ fixed(0.09) # script PdPi_v sd 0.3 * mean
    etalps_out ~ fixed(0.09) # script PdePi_v sd 0.3 * mean
    etalclint_met ~ fixed(0.09) # script CL_MePi_v sd 0.3 * mean
    etalclint_bile ~ fixed(0.04) # script CL_BiPi_v sd 0.2 * mean
    etalclurine ~ fixed(0.04) # script CLur_v sd 0.2 * mean
    etalfu_blood ~ fixed(0.04) # script Cfpp_v sd 0.2 * mean
    etalfu_liver ~ fixed(0.09) # script fTP_v sd 0.3 * mean
    etalvc ~ fixed(0.01) # script VcP_v sd 0.1 * mean
    etalka_eltrombopag ~ fixed(0.04) # script kaE_v sd 0.2 * mean
    etalvmax_eltrombopag ~ fixed(0.09) # script VmE_v sd 0.3 * mean
    etalkm_eltrombopag ~ fixed(0.09) # script KmE_v sd 0.3 * mean
    etalps_in_eltrombopag ~ fixed(0.09) # script PdE_v sd 0.3 * mean
    etalps_out_eltrombopag ~ fixed(0.09) # script PdeE_v sd 0.3 * mean
    etalclint_met_eltrombopag ~ fixed(0.04) # script CL_MeE_v sd 0.2 * mean
    etalclint_bile_eltrombopag ~ fixed(0.04) # script CL_BiE_v sd 0.2 * mean
    etalfu_eltrombopag ~ fixed(0.01) # script CfpE_v sd 0.1 * mean
    etalfu_liver_eltrombopag ~ fixed(0.04) # script fTE_v sd 0.2 * mean
    etalvc_eltrombopag ~ fixed(0.04) # script VcE_v sd 0.2 * mean
    etalv_liver_ext ~ fixed(0.04) # script Vext_H_v sd 0.2 * mean
    etalv_liver ~ fixed(0.04) # script VH_v sd 0.2 * mean
    etalv_gallbladder ~ fixed(0.04) # script VGaBl_v sd 0.2 * mean
    etalq_liver ~ fixed(0.04) # script Qh_v sd 0.2 * mean
    etalq_kidney ~ fixed(0.04) # script Qk_v sd 0.2 * mean
  })

  model({
    # Individual parameters: normal (not log-normal) Monte-Carlo draws, as in
    # the authors' rnorm() sampling.
    kst <- exp(lkst) * (1 + etalkst)
    kehc <- exp(lkehc) * (1 + etalkehc)
    ka <- exp(lka) * (1 + etalka)
    vmax <- exp(lvmax) * (1 + etalvmax)
    km <- exp(lkm) * (1 + etalkm)
    ps_in <- exp(lps_in) * (1 + etalps_in)
    ps_out <- exp(lps_out) * (1 + etalps_out)
    clint_met <- exp(lclint_met) * (1 + etalclint_met)
    clint_bile <- exp(lclint_bile) * (1 + etalclint_bile)
    clurine <- exp(lclurine) * (1 + etalclurine)
    fu_blood <- exp(lfu_blood) * (1 + etalfu_blood)
    fu_liver <- exp(lfu_liver) * (1 + etalfu_liver)
    vc <- exp(lvc) * (1 + etalvc)
    ka_eltrombopag <- exp(lka_eltrombopag) * (1 + etalka_eltrombopag)
    vmax_eltrombopag <- exp(lvmax_eltrombopag) * (1 + etalvmax_eltrombopag)
    km_eltrombopag <- exp(lkm_eltrombopag) * (1 + etalkm_eltrombopag)
    ps_in_eltrombopag <- exp(lps_in_eltrombopag) * (1 + etalps_in_eltrombopag)
    ps_out_eltrombopag <- exp(lps_out_eltrombopag) * (1 + etalps_out_eltrombopag)
    clint_met_eltrombopag <- exp(lclint_met_eltrombopag) * (1 + etalclint_met_eltrombopag)
    clint_bile_eltrombopag <- exp(lclint_bile_eltrombopag) * (1 + etalclint_bile_eltrombopag)
    fu_eltrombopag <- exp(lfu_eltrombopag) * (1 + etalfu_eltrombopag)
    fu_liver_eltrombopag <- exp(lfu_liver_eltrombopag) * (1 + etalfu_liver_eltrombopag)
    vc_eltrombopag <- exp(lvc_eltrombopag) * (1 + etalvc_eltrombopag)
    v_liver_ext <- exp(lv_liver_ext) * (1 + etalv_liver_ext)
    v_liver <- exp(lv_liver) * (1 + etalv_liver)
    v_gallbladder <- exp(lv_gallbladder) * (1 + etalv_gallbladder)
    q_liver <- exp(lq_liver) * (1 + etalq_liver)
    q_kidney <- exp(lq_kidney) * (1 + etalq_kidney)

    # Hepatic uptake clearances (mL/min) with mutual competitive inhibition.
    # The inhibitor term uses the other drug's TOTAL extracellular
    # concentration, as in Eqs. S43/S48 and the script.
    clup <- vmax / (km * (1 + is_liver_eltrombopag / km_eltrombopag) + fu_blood * is_liver) + ps_in
    clup_eltrombopag <- vmax_eltrombopag / (km_eltrombopag * (1 + is_liver / km) + fu_eltrombopag * is_liver_eltrombopag) + ps_in_eltrombopag

    # Pitavastatin, Eqs. S41-S46 as coded in the script
    d/dt(stomach) <- -kst * stomach
    d/dt(gut_lumen) <- kst * stomach - ka * gut_lumen + kehc * gallbladder
    d/dt(is_liver) <- (ka * gut_lumen - fu_blood * is_liver * clup + q_liver * fu_blood * (central - is_liver) + ps_out * fu_liver * int_liver) / v_liver_ext
    d/dt(int_liver) <- (fu_blood * is_liver * clup - (clint_met + clint_bile + ps_out) * fu_liver * int_liver) / v_liver
    d/dt(gallbladder) <- clint_bile * fu_liver * int_liver / v_gallbladder - kehc * gallbladder
    # Urinary loss as coded in the script: Q_K * CL_ur * fu_bl * C (Eq. S46
    # prints Q_K * CL_urine.P * S6 without fu_bl).
    d/dt(central) <- (q_liver * fu_blood * (is_liver - central) - q_kidney * clurine * fu_blood * central) / vc

    # Eltrombopag, Eqs. S47-S51 as coded in the script (which passes the
    # eltrombopag dose through the stomach with the shared k_ge)
    d/dt(stomach_eltrombopag) <- -kst * stomach_eltrombopag
    d/dt(gut_lumen_eltrombopag) <- kst * stomach_eltrombopag - ka_eltrombopag * gut_lumen_eltrombopag + kehc * gallbladder_eltrombopag
    d/dt(is_liver_eltrombopag) <- (ka_eltrombopag * gut_lumen_eltrombopag - fu_eltrombopag * is_liver_eltrombopag * clup_eltrombopag + q_liver * fu_eltrombopag * (central_eltrombopag - is_liver_eltrombopag) + ps_out_eltrombopag * fu_liver_eltrombopag * int_liver_eltrombopag) / v_liver_ext
    d/dt(int_liver_eltrombopag) <- (fu_eltrombopag * is_liver_eltrombopag * clup_eltrombopag - (clint_met_eltrombopag + clint_bile_eltrombopag + ps_out_eltrombopag) * fu_liver_eltrombopag * int_liver_eltrombopag) / v_liver
    d/dt(gallbladder_eltrombopag) <- clint_bile_eltrombopag * fu_liver_eltrombopag * int_liver_eltrombopag / v_gallbladder - kehc * gallbladder_eltrombopag
    d/dt(central_eltrombopag) <- q_liver * fu_eltrombopag * (is_liver_eltrombopag - central_eltrombopag) / vc_eltrombopag

    # Outputs (ng/mL): pitavastatin in whole blood, eltrombopag in plasma
    Cc <- central
    Cc_eltrombopag <- central_eltrombopag
  })
}

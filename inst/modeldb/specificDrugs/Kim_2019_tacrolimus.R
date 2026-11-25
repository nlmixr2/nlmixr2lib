Kim_2019_tacrolimus <- function() {
  description <- "Integrated population PK model of tacrolimus (TAC) and mycophenolate mofetil (MMF) with the TAC-MMF drug-drug interaction in stable Korean kidney transplant recipients (Kim 2019, final model). TAC follows a two-compartment model with first-order absorption and a lag time; apparent oral clearance (CL/F) is 1.49-fold higher in CYP3A5 expressers and is suppressed by plasma mycophenolic acid (MPA) through an inverse-exponential interaction (CL/F = 21.9 / exp(0.06*Cmpa) * 1.49^CYP3A5). MPA follows a two-compartment model with first-order absorption and is eliminated only by metabolism to MPAG (7-O-glucuronide; 85%) and AcMPAG (acyl glucuronide; 15%). MPAG undergoes enterohepatic recirculation through a gallbladder compartment that empties into the MPA absorption compartment during a meal-triggered window. The SLCO1B3 rs4149117 T allele increases the MPAG volume 1.2-fold and the UGT2B7 rs7439366 T allele lowers the AcMPAG formation rate constant to 0.812-fold. Tacrolimus concentrations are in ng/mL (whole blood); MPA, MPAG and AcMPAG are in ug/mL (plasma)."
  reference <- paste(
    "Kim JH, Han N, Kim MG, Kim YW, Jang H, Yun H-Y, Yu M-Y, Kim I-W, Kim YS, Oh JM.",
    "Model based development of tacrolimus dosing algorithm considering CYP3A5",
    "genotypes and mycophenolate mofetil drug interaction in stable kidney",
    "transplant recipients. Sci Rep. 2019;9:11740. doi:10.1038/s41598-019-47876-0.",
    sep = " "
  )
  vignette <- "Kim_2019_tacrolimus"
  # Oral tacrolimus into depot; oral mycophenolate mofetil into depot_mpa
  dosing <- c("depot", "depot_mpa")
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (tacrolimus); ug/mL (MPA, MPAG, AcMPAG)"
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "tissue", verified = TRUE),
    depot_mpa = list(
      analyte = "mycophenolate mofetil / mycophenolic acid (MPA)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_mpa = list(analyte = "mycophenolic acid (MPA)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_mpa = list(analyte = "mycophenolic acid (MPA)", units = "mg", specimen = "tissue", verified = TRUE),
    central_mpag = list(
      analyte = "mycophenolic acid 7-O-glucuronide (MPAG)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    gallbladder = list(
      analyte = "mycophenolic acid 7-O-glucuronide (MPAG)",
      units = "mg",
      specimen = "bile",
      verified = TRUE
    ),
    central_acmpag = list(
      analyte = "mycophenolic acid acyl glucuronide (AcMPAG)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status: 1 = CYP3A5 expresser (CYP3A5*1/*1 or *1/*3 at rs776746), 0 = nonexpresser (CYP3A5*3/*3).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5 nonexpresser, *3/*3; 19 of 32 patients in the Kim 2019 cohort)",
      notes = "Time-fixed germline genotype; 13 of 32 (40.6%) patients were expressers (Table 1). Multiplicative power-form effect on tacrolimus CL/F: cl = 21.9 * 1.49^CYP3A5_EXPR (Eq. 1; Table 2).",
      source_name = "CYP3A5"
    ),
    SNP_SLCO1B3_RS4149117_T_CARRIER = list(
      description = "SLCO1B3 334T>G (rs4149117) T-allele carrier: 1 = carries at least one T allele, 0 = GG.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (GG genotype; 17 of 32 patients)",
      notes = "Time-fixed germline genotype; 15 of 32 (46.9%) patients were T carriers (Table 1). Multiplicative power-form effect on the MPAG volume: vc_mpag = 1.76 * 1.2^SNP_SLCO1B3_RS4149117_T_CARRIER (Table 2 'SLCO1B3 on V MPAG').",
      source_name = "SLCO1B3"
    ),
    SNP_UGT2B7_RS7439366_T_CARRIER = list(
      description = "UGT2B7 802C>T (rs7439366) T-allele carrier: 1 = carries at least one T allele, 0 = CC.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CC genotype; 17 of 32 patients)",
      notes = "Time-fixed germline genotype; 15 of 32 (46.9%) patients were T carriers (Table 1). Multiplicative power-form effect on the MPA -> AcMPAG formation rate constant: k59 = k59_typ * 0.812^SNP_UGT2B7_RS7439366_T_CARRIER (Table 2 'UGT2B7 on k59'). The paper writes the variant both as 802T>C (Results) and 802C>T (Methods, Discussion); the carrier allele (T) and the reference genotype (CC) are stated consistently.",
      source_name = "UGT2B7"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    n_observations = "Up to 5 samples per patient (median 4; pre-dose and 1, 2, 3 and 4 h post-dose) for each of TAC, MPA, MPAG and AcMPAG.",
    age_range = "20-70 years",
    age_median = "52 years (Table 1)",
    weight_range = "43.9-102.4 kg",
    weight_median = "62.7 kg (Table 1)",
    height_range = "151-180 cm (median 165.8)",
    sex_female_pct = 37.5,
    race_ethnicity = "100% Korean.",
    disease_state = "Stable kidney transplant recipients at least 6 months post-transplant (median 5.7 years, range 0.6-10.4), on maintenance tacrolimus (Prograf) + mycophenolate mofetil (CellCept) + prednisolone.",
    dose_range = "Tacrolimus 1-6 mg/day (median 2 mg/day) and mycophenolate mofetil 500-2,000 mg/day (median 1,000 mg/day), oral, twice daily, at a stable dose for at least 2 weeks; prednisolone (or equivalent) 2.5-5 mg/day.",
    regions = "Republic of Korea (Seoul National University Hospital).",
    renal_function = "MDRD eGFR 35.0-85.7 mL/min/1.73 m2 (median 59.5); serum creatinine 0.69-1.66 mg/dL (median 1.22).",
    hematocrit = "32.3-53.8% (median 43.6).",
    genotypes = "CYP3A5 expressers 13/32 (40.6%); SLCO1B3 rs4149117 T carriers 15/32 (46.9%); UGT2B7 rs7439366 T carriers 15/32 (46.9%) (Table 1).",
    notes = "Single-centre steady-state sampling study (NCT02808065). Baseline demographics from Table 1. Absorption, lag-time and enterohepatic-circulation parameters, and the AcMPAG volume and elimination rate constant, were fixed rather than estimated because the 4-h sampling window could not inform them."
  )

  ini({
    # Final-model estimates from Kim 2019 Table 2 ('Final model' columns).
    # IIV is reported as %CV and converted with omega^2 = log(1 + CV^2).
    # Rate-constant numbering follows Figure 2 (1 = TAC GI, 2 = TAC central,
    # 3 = TAC peripheral, 4 = MPA GI, 5 = MPA central, 6 = MPA peripheral,
    # 7 = MPAG, 8 = gallbladder, 9 = AcMPAG).

    # Tacrolimus: 2-compartment, first-order absorption with lag time
    lka   <- fixed(log(1.78)); label("Tacrolimus absorption rate constant Ka (1/h)")          # Table 2 Ka TAC = 1.78 FIX
    lcl   <- log(21.9);        label("Tacrolimus apparent oral clearance CL/F (L/h)")         # Table 2 CL TAC = 21.9 L/h; Eq. 1
    lvc   <- log(103);         label("Tacrolimus apparent central volume V/F (L)")            # Table 2 V TAC = 103 L
    lk23  <- log(0.803);       label("Tacrolimus central-to-peripheral rate constant k23 (1/h)") # Table 2 k23 = 0.803 1/h
    lk32  <- log(0.583);       label("Tacrolimus peripheral-to-central rate constant k32 (1/h)") # Table 2 k32 = 0.583 1/h
    ltlag <- fixed(log(0.59)); label("Tacrolimus absorption lag time (h)")                    # Table 2 Lag time = 0.59 FIX

    e_cyp3a5_cl <- 1.49; label("CYP3A5 expresser multiplier on tacrolimus CL/F (unitless)") # Table 2 CYP3A5 on CL TAC = 1.49; Eq. 1
    # MPA plasma concentration (ug/mL) lowers TAC CL/F as 1/exp(slope * Cmpa).
    e_cmpa_cl <- 0.06; label("Inverse-exponential slope of MPA concentration on tacrolimus CL/F (mL/ug)") # Table 2 Interaction = 0.06; Eq. 1

    # Mycophenolic acid: 2-compartment, first-order absorption; dose is the
    # administered MMF amount (mg), so CL/F and V/F absorb F and MMF -> MPA.
    lka_mpa <- fixed(log(2.29)); label("MPA absorption rate constant Ka (1/h)")            # Table 2 Ka MPA = 2.29 FIX
    lcl_mpa <- log(3.27);        label("MPA apparent oral clearance CL/F (L/h)")           # Table 2 CL MPA = 3.27 L/h
    lvc_mpa <- log(23.2);        label("MPA apparent central volume V/F (L)")              # Table 2 V MPA = 23.2 L
    lk56    <- log(1.75);        label("MPA central-to-peripheral rate constant k56 (1/h)") # Table 2 k56 = 1.75 1/h
    lk65    <- log(0.0089);      label("MPA peripheral-to-central rate constant k65 (1/h)") # Table 2 k65 = 0.0089 1/h
    f_mpag  <- fixed(0.85);      label("Fraction of MPA metabolism forming MPAG (fMPA)")   # Table 2 f MPA = 0.85 FIX
    e_ugt2b7_k59 <- 0.812;       label("UGT2B7 rs7439366 T-carrier multiplier on AcMPAG formation rate constant k59 (unitless)") # Table 2 UGT2B7 on k59 = 0.812

    # MPAG with enterohepatic recirculation through the gallbladder
    lvc_mpag <- log(1.76);        label("MPAG apparent volume V/F (L)")                      # Table 2 V MPAG = 1.76 L
    e_slco1b3_vc_mpag <- 1.2;     label("SLCO1B3 rs4149117 T-carrier multiplier on MPAG volume (unitless)") # Table 2 SLCO1B3 on V MPAG = 1.2
    lk70     <- log(0.103);       label("MPAG elimination rate constant k70 (1/h)")          # Table 2 k70 = 0.103 1/h
    lehc     <- fixed(log(0.367)); label("Fraction of MPAG elimination via the bile (EHC)")  # Table 2 EHC = 0.367 FIX
    lk84     <- fixed(log(18.4)); label("Gallbladder emptying rate constant k84 (1/h)")      # Table 2 k84 = 18.4 FIX
    lmtime1  <- fixed(log(7.96)); label("Gallbladder emptying start MTIME1 (h after the dose)") # Table 2 MTIME1 = 7.96 FIX
    lmtime2  <- fixed(log(1));    label("Gallbladder emptying duration MTIME2 (h)")          # Table 2 MTIME2 = 1 FIX

    # AcMPAG: 1-compartment
    lvc_acmpag <- fixed(log(13.9)); label("AcMPAG apparent volume V/F (L)")                 # Table 2 V AcMPAG = 13.9 FIX
    lk90       <- fixed(log(0.407)); label("AcMPAG elimination rate constant k90 (1/h)")   # Table 2 k90 = 0.407 FIX

    # IIV (exponential); Table 2 final-model IIV CV%
    etalcl     ~ 0.05783 # Table 2 IIV CL TAC = 24.4% CV
    etalvc     ~ 0.1882  # Table 2 IIV V TAC = 45.5% CV
    etalcl_mpa ~ 0.06444 # Table 2 IIV CL MPA = 25.8% CV
    etalvc_mpa ~ 0.1689  # Table 2 IIV V MPA = 42.9% CV

    # Residual error: TAC and AcMPAG combined; MPA and MPAG proportional
    addSd         <- 0.585;  label("Additive residual error, tacrolimus (ng/mL)")    # Table 2 sigma add TAC = 0.585
    propSd        <- 0.211;  label("Proportional residual error, tacrolimus (fraction)") # Table 2 sigma prop TAC = 0.211
    propSd_mpa    <- 0.515;  label("Proportional residual error, MPA (fraction)")    # Table 2 sigma prop MPA = 0.515
    propSd_mpag   <- 0.216;  label("Proportional residual error, MPAG (fraction)")   # Table 2 sigma prop MPAG = 0.216
    addSd_acmpag  <- 0.0647; label("Additive residual error, AcMPAG (ug/mL)")        # Table 2 sigma add AcMPAG = 0.0647
    propSd_acmpag <- 0.268;  label("Proportional residual error, AcMPAG (fraction)") # Table 2 sigma prop AcMPAG = 0.268
  })

  model({
    ka   <- exp(lka)
    vc   <- exp(lvc + etalvc)
    k23  <- exp(lk23)
    k32  <- exp(lk32)
    tlag <- exp(ltlag)

    ka_mpa <- exp(lka_mpa)
    cl_mpa <- exp(lcl_mpa + etalcl_mpa)
    vc_mpa <- exp(lvc_mpa + etalvc_mpa)
    k56    <- exp(lk56)
    k65    <- exp(lk65)

    vc_mpag <- exp(lvc_mpag) * e_slco1b3_vc_mpag^SNP_SLCO1B3_RS4149117_T_CARRIER
    k70     <- exp(lk70)
    ehc     <- exp(lehc)
    k84     <- exp(lk84)
    mtime1  <- exp(lmtime1)
    mtime2  <- exp(lmtime2)

    vc_acmpag <- exp(lvc_acmpag)
    k90       <- exp(lk90)

    # MPA plasma concentration in mg/L = ug/mL, the unit of the interaction slope
    Cmpa <- central_mpa / vc_mpa

    # Eq. 1: CL/F of TAC = 21.9 * 1/exp(0.06 * Cmpa) * 1.49^CYP3A5
    cl  <- exp(lcl + etalcl) / exp(e_cmpa_cl * Cmpa) * e_cyp3a5_cl^CYP3A5_EXPR
    kel <- cl / vc

    # MPA leaves the central compartment only by metabolism (Figure 2), so
    # k57 + k59 = CL/F MPA / V/F MPA, split by fMPA. The UGT2B7 effect acts on
    # k59 alone. EHC = k78 / (k70 + k78) gives k78 = EHC * k70 / (1 - EHC).
    ke_mpa <- cl_mpa / vc_mpa
    k57    <- f_mpag * ke_mpa
    k59    <- (1 - f_mpag) * ke_mpa * e_ugt2b7_k59^SNP_UGT2B7_RS7439366_T_CARRIER
    k78    <- ehc * k70 / (1 - ehc)

    # Gallbladder empties into the MPA GI compartment from MTIME1 to
    # MTIME1 + MTIME2 h after each MMF dose. The window is anchored to the most
    # recent MMF dose so it recurs every dosing interval; before any MMF dose
    # tad() is undefined and the comparison is false, so no emptying occurs.
    tad_mmf  <- tad(depot_mpa)
    empty_on <- 0
    if (tad_mmf >= mtime1 && tad_mmf <= mtime1 + mtime2) {
      empty_on <- 1
    }
    k84_eff <- k84 * empty_on

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - kel * central - k23 * central + k32 * peripheral1
    d/dt(peripheral1) <- k23 * central - k32 * peripheral1
    lag(depot)        <- tlag

    d/dt(depot_mpa)       <- -ka_mpa * depot_mpa + k84_eff * gallbladder
    d/dt(central_mpa)     <- ka_mpa * depot_mpa - k56 * central_mpa + k65 * peripheral1_mpa - k57 * central_mpa - k59 * central_mpa
    d/dt(peripheral1_mpa) <- k56 * central_mpa - k65 * peripheral1_mpa

    d/dt(central_mpag) <- k57 * central_mpa - k70 * central_mpag - k78 * central_mpag
    d/dt(gallbladder)  <- k78 * central_mpag - k84_eff * gallbladder

    d/dt(central_acmpag) <- k59 * central_mpa - k90 * central_acmpag

    # Tacrolimus: mg/L * 1000 = ng/mL. MPA, MPAG, AcMPAG: mg/L = ug/mL.
    Cc        <- 1000 * central / vc
    Cc_mpa    <- central_mpa / vc_mpa
    Cc_mpag   <- central_mpag / vc_mpag
    Cc_acmpag <- central_acmpag / vc_acmpag

    Cc        ~ add(addSd) + prop(propSd)
    Cc_mpa    ~ prop(propSd_mpa)
    Cc_mpag   ~ prop(propSd_mpag)
    Cc_acmpag ~ add(addSd_acmpag) + prop(propSd_acmpag)
  })
}

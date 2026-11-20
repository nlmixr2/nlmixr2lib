Hughes_2025_clesacostat <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption for oral clesacostat (PF-05221304, a liver-targeted acetyl-CoA carboxylase inhibitor) in healthy participants, participants with hepatic impairment, and adults with non-alcoholic fatty liver disease (NAFLD) or presumed non-alcoholic steatohepatitis (NASH), alone or co-administered with ervogastat. Absorption is faster with fasting / fructose administration and with tablet formulation B; relative bioavailability is lower with fasting / fructose administration and rises stepwise at doses of 100, 200 and 240 mg. Apparent clearance scales allometrically with body weight, is higher with concomitant ervogastat, and is lower in NAFLD / presumed NASH and in SLCO1B1 intermediate or poor metabolizers. Interindividual variability on CL/F; combined additive and time-after-dose-dependent proportional residual error. This is the clesacostat exposure model whose individual average steady-state concentrations drive the Hughes 2025 liver-fat exposure-response model (Hughes_2025_ervogastat_clesacostat_liverfat)."
  reference <- paste(
    "Hughes JH, Amin NB, Wojciechowski J, Vourvahis M.",
    "Exposure-response modeling of liver fat imaging endpoints in",
    "non-alcoholic fatty liver disease populations administered ervogastat",
    "alone and co-administered with clesacostat.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):317-330.",
    "doi:10.1002/psp4.13275",
    sep = " "
  )
  vignette <- "Hughes_2025_ervogastat_clesacostat"
  paper_specific_residual_sds <- c("propSdMax", "propSdMin")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline body weight. Allometric power function referenced",
        "to 70 kg with exponents fixed to 0.75 on CL/F and Q/F and 1 on Vc/F",
        "and Vp/F (Supplement Table S4; Code S2 ALLMCL / ALLMVD)."
      ),
      source_name = "BWT"
    ),
    DOSE_CLESACOSTAT_MG = list(
      description = "Clesacostat amount per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source DOSE2 ('actual clesacostat dose'). Selects the stepwise",
        "relative-bioavailability factor: 1 below 100 mg, 1.110944 at",
        ">= 100 mg, 1.252650 at >= 200 mg and 1.288417 at >= 240 mg (Code S2",
        "DOSEFREL; each threshold overwrites the previous one, so the factors",
        "are not cumulative). Must equal the amt of the dose record."
      ),
      source_name = "DOSE2"
    ),
    CONMED_ERVOGASTAT = list(
      description = "Concomitant ervogastat indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (clesacostat monotherapy)",
      notes = paste(
        "1 = ervogastat co-administered (Code S2 'IF (DOSE1.GT.0)'), 0 =",
        "clesacostat alone. A treatment effect of +27.3884% on CL/F (Table S4",
        "'Effect of concomitant ervogastat on CL/F'); the Code S2 comment",
        "notes that an Emax relationship collapsed to an EC50 of about 0.001,",
        "hence the on/off form."
      ),
      source_name = "DOSE1"
    ),
    FED = list(
      description = "Fed state at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (fed: standard meal, high-fat meal, or dosed without regard to food)",
      notes = paste(
        "Code S2 FOOD: 0 fasting, 1 standard meal, 2 no regard to food, 3",
        "high-fat meal, 4 fructose. The fasting and fructose categories share",
        "one effect (Table S4 'Effect of fasting/fructose during",
        "administration'): +77.5485% on ka and -9.31887% on Frel. Encode FED",
        "= 0 for fasting and for fructose administration, FED = 1 for the",
        "three meal categories, which form the pooled reference. The model",
        "uses (1 - FED) as the fasting/fructose indicator. The reference",
        "population for ka and Frel in the Supplement (Figure S3 caption) is",
        "formulation A 'with food'."
      ),
      source_name = "FOOD"
    ),
    FORM_CLESACOSTAT_TABLETB = list(
      description = "Clesacostat tablet formulation B indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet formulation A)",
      notes = paste(
        "1 = tablet formulation B, alone or with the ervogastat",
        "immediate-release tablet (Code S2 FORM = 7 or 9), 0 = tablet",
        "formulation A (FORM = 6 or 8). Raises ka by 109.646% (Table S4",
        "'Effect of tablet formulation B on ka'). Formulation B was used in",
        "NCT03871439, NCT03597217 and NCT04399538 (Table S1)."
      ),
      source_name = "FORM"
    ),
    DIS_NAFLD = list(
      description = "Non-alcoholic fatty liver disease or presumed non-alcoholic steatohepatitis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant or participant in a hepatic-impairment cohort)",
      notes = paste(
        "1 = NAFLD or presumed NASH (Code S2 PTST = 4 or 5), 0 otherwise;",
        "-29.3157% on CL/F (Table S4). The two categories were combined into",
        "one estimate on the similarity of their bootstrap distributions",
        "(Supplement Results). Hepatic impairment was tested and judged not",
        "clinically important for clesacostat, so it has no effect here."
      ),
      source_name = "PTST"
    ),
    SLCO1B1_PM_IM = list(
      description = "Pooled SLCO1B1 intermediate-or-poor transporter-activity phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive SLCO1B1 activity, or phenotype not determined)",
      notes = paste(
        "Code S2 SLCO1B1: 0 not determined / wild type, 1-3 extensive, 4-7",
        "intermediate, 8-10 poor metabolizer; the effect applies when",
        "SLCO1B1 > 3, i.e. intermediate or poor (-9.99308% on CL/F, Table S4",
        "'Effect of SLCO1B1 IM/PM activity on CL/F'). Undetermined phenotype",
        "is pooled with extensive activity in the reference."
      ),
      source_name = "SLCO1B1"
    )
  )

  covariatesDataExcluded <- list(
    HEPIMP = list(
      description = "Any hepatic impairment indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Tested on CL/F in the full covariate model and classified as not clinically important; not retained (Supplement Results, Figure S3).",
      source_name = "PTST"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (reference 49 years) and classified as not clinically important; not retained (Supplement Results, Figure S3).",
      source_name = "Age"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "clesacostat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clesacostat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "clesacostat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 487L,
    n_studies = 8L,
    n_observations = 6722L,
    age_range = "19-72 years",
    age_median = "52 years",
    weight_range = "52-170 kg",
    weight_median = "90 kg",
    sex_female_pct = 39.6,
    race_ethnicity = c(White = 70.6, Black = 12.7, Asian = 11.5, Other = 5.1),
    disease_state = paste(
      "Healthy volunteers (24%), mild / moderate / severe hepatic impairment",
      "(1.2% each), NAFLD (27.1%) and presumed NASH (45.2%)"
    ),
    dose_range = paste(
      "1-240 mg single doses; 1-100 mg twice daily or 40-100 mg once daily;",
      "2-50 mg once daily for 16 weeks; 2.5-20 mg twice or once daily with",
      "ervogastat (Table S1)"
    ),
    regions = "United States; includes healthy Japanese participants (NCT03597217)",
    notes = paste(
      "Demographics from Hughes 2025 Table 1 (clesacostat popPK column).",
      "6928 concentrations measured, 6722 evaluable; the 3% below the 2.0",
      "ng/mL LLOQ were excluded (M1). SLCO1B1 phenotype: extensive 61.2%,",
      "intermediate 18.1%, poor 3.3%, not determined 17.5%. Estimation by",
      "FOCE-I in NONMEM 7.5 (Code S2, run100)."
    )
  )

  ini({
    # Structural parameters (Supplement Code S2 $THETA final estimates;
    # Table S4 prints the same values rounded). Reference: 70 kg, tablet
    # formulation A with food, dose below 100 mg, clesacostat alone, healthy
    # participant with extensive (or undetermined) SLCO1B1 activity.
    lka <- log(0.450538); label("First-order absorption rate constant (1/h)") # Code S2 THETA(1) TVKA = 0.450538; Table S4 'ka' = 0.451
    lcl <- log(1.10516); label("Apparent clearance CL/F (L/h)") # Code S2 THETA(2) TVCL = 1.10516; Table S4 'CL/F' = 1.11
    lvc <- log(10.5885); label("Apparent central volume Vc/F (L)") # Code S2 THETA(3) TVV2 = 10.5885; Table S4 'Vc/F' = 10.6
    lq <- log(0.541464); label("Apparent intercompartmental clearance Q/F (L/h)") # Code S2 THETA(4) TVQ = 0.541464; Table S4 'Q/F' = 0.541
    lvp <- log(6.56800); label("Apparent peripheral volume Vp/F (L)") # Code S2 THETA(5) TVV3 = 6.56800; Table S4 second 'Vc/F' row (the peripheral volume, mislabelled) = 6.57
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the reference condition (unitless)") # Code S2 POPF1 = 1*FOODFREL*DOSEFREL; structural anchor, no THETA

    # Allometric exponents, fixed (Table S4 'Fixed'; Code S2 '0.75 FIX', '1.00 FIX')
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/70 on CL/F and Q/F (unitless)") # Code S2 THETA(10) TVALLMCL = 0.75 FIX
    e_wt_vc <- fixed(1.00); label("Allometric exponent of WT/70 on Vc/F and Vp/F (unitless)") # Code S2 THETA(11) TVALLMVD = 1.00 FIX

    # Covariate effects (Code S2: categorical effects enter as 1 + theta)
    e_conmed_ervogastat_cl <- 0.273884; label("Concomitant ervogastat fractional change in CL/F (fraction)") # Code S2 THETA(12) TVDGATCL = 0.273884; Table S4 = 0.274
    e_fasted_ka <- 0.775485; label("Fasting / fructose administration fractional change in ka (fraction)") # Code S2 THETA(13) TVFOOD04KA = 0.775485; Table S4 = 0.775
    e_fasted_fdepot <- -0.0931887; label("Fasting / fructose administration fractional change in Frel (fraction)") # Code S2 THETA(14) TVFOOD04FREL = -0.0931887; Table S4 = -0.0932
    e_nafld_cl <- -0.293157; label("NAFLD or presumed NASH fractional change in CL/F (fraction)") # Code S2 THETA(15) TVPTST45CL = -0.293157; Table S4 = -0.293
    e_dose100_fdepot <- 0.110944; label("Fractional change in Frel at clesacostat doses >= 100 mg (fraction)") # Code S2 THETA(16) TVDOSE100FREL = 0.110944; Table S4 = 0.111
    e_dose200_fdepot <- 0.252650; label("Fractional change in Frel at clesacostat doses >= 200 mg (fraction)") # Code S2 THETA(17) TVDOSE200FREL = 0.252650; Table S4 = 0.253
    e_dose240_fdepot <- 0.288417; label("Fractional change in Frel at clesacostat doses >= 240 mg (fraction)") # Code S2 THETA(18) TVDOSE240FREL = 0.288417; Table S4 = 0.288
    e_tabletb_ka <- 1.09646; label("Tablet formulation B fractional change in ka (fraction)") # Code S2 THETA(19) TVFORMKA = 1.09646; Table S4 = 1.10
    e_slco1b1_cl <- -0.0999308; label("SLCO1B1 intermediate / poor activity fractional change in CL/F (fraction)") # Code S2 THETA(20) TVSLCO1B1CL = -0.0999308; Table S4 = -0.0999

    # IIV (Code S2 $OMEGA variance; Table S4 prints omega as SD 0.269)
    etalcl ~ 0.0737692 # Code S2 $OMEGA PPVCL = 0.0737692; Table S4 omega CL/F 0.269

    # Residual error: W = sqrt((IPRED*PRO)^2 + ADD^2), EPS variance 1 FIX,
    # with PRO = RUVMIN + (RUVMAX - RUVMIN) * exp(-RUVKT * TAD)
    propSdMax <- 1.21740; label("Proportional residual SD at the time of dosing (fraction)") # Code S2 THETA(6) TVRUVPROMAX = 1.21740; Table S4 'Maximum proportional residual error' = 1.22
    propSdMin <- 0.272273; label("Asymptotic proportional residual SD late after dosing (fraction)") # Code S2 THETA(7) TVRUVPROMIN = 0.272273; Table S4 'Minimum proportional residual error' = 0.272
    lkruv <- log(0.985922); label("Rate constant of the decline of the proportional residual SD with time after dose (1/h)") # Code S2 THETA(8) TVRUVKT = 0.985922; Table S4 'kRUV' = 0.986
    addSd <- 2.46300; label("Additive residual SD (ng/mL)") # Code S2 THETA(9) TVRUVADD = 2.46300; Table S4 = 2.46
  })

  model({
    # Allometric scaling referenced to 70 kg (Code S2 ALLMCL / ALLMVD)
    allomCl <- (WT / 70)^e_wt_cl
    allomV <- (WT / 70)^e_wt_vc

    # Fasting or fructose administration (Code S2 FOOD = 0 or 4)
    fasted <- 1 - FED

    # Stepwise dose effect on Frel; each threshold overwrites the previous
    # one (Code S2 DOSEFREL)
    dose_frel <- 1
    if (DOSE_CLESACOSTAT_MG >= 100) dose_frel <- 1 + e_dose100_fdepot
    if (DOSE_CLESACOSTAT_MG >= 200) dose_frel <- 1 + e_dose200_fdepot
    if (DOSE_CLESACOSTAT_MG >= 240) dose_frel <- 1 + e_dose240_fdepot

    fdepot <- exp(lfdepot) * (1 + e_fasted_fdepot * fasted) * dose_frel
    ka <- exp(lka) * (1 + e_fasted_ka * fasted) * (1 + e_tabletb_ka * FORM_CLESACOSTAT_TABLETB)
    cl <- exp(lcl + etalcl) * allomCl *
      (1 + e_conmed_ervogastat_cl * CONMED_ERVOGASTAT) *
      (1 + e_nafld_cl * DIS_NAFLD) *
      (1 + e_slco1b1_cl * SLCO1B1_PM_IM)
    vc <- exp(lvc) * allomV
    q <- exp(lq) * allomCl
    vp <- exp(lvp) * allomV

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Time-dependent proportional residual error (Karlsson 1995 form in
    # Code S2); TAD is 0 before the first dose in the control stream.
    tadRuv <- tad(depot)
    if (is.na(tadRuv)) tadRuv <- 0
    kruv <- exp(lkruv)
    propSdTad <- propSdMin + (propSdMax - propSdMin) * exp(-kruv * tadRuv)

    # Dose in mg, volume in L: x 1000 for ng/mL (Code S2 S2 = V2/1000)
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSdTad)
  })
}

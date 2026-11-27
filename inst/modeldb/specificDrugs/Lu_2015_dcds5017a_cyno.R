Lu_2015_dcds5017a_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Semi-mechanistic multiple-analyte population PK model for the anti-CD79b-MMAE antibody-drug conjugate DCDS5017A (Lu 2015): nine DAR species (DAR0-DAR8), each two-compartment, sharing Vc, Vp and CLd; per-species Michaelis-Menten proteolytic degradation and sequential one-MMAE-at-a-time deconjugation whose rate constant depends on DAR through a Weibull function; unconjugated MMAE two-compartment with input from both catabolic pathways. Outputs antibody-conjugated MMAE (acMMAE), total antibody (Tab) and unconjugated MMAE."
  reference <- "Lu D, Jin JY, Girish S, Agarwal P, Li D, Prabhu S, Dere R, Saad O, Nazzal D, Koppada N, Ramanujan S, Ng CM. Semi-mechanistic Multiple-Analyte Pharmacokinetic Model for an Antibody-Drug-Conjugate in Cynomolgus Monkeys. Pharm Res. 2015;32(6):1907-1919. doi:10.1007/s11095-014-1585-y. PMID 25467958; PMCID PMC4422865."
  vignette <- "Lu_2015_dcds5017a_cyno"
  dosing <- c(
    "dar8_central",
    "dar6_central",
    "dar4_central",
    "dar2_central",
    "dar1_central",
    "dar0_central",
    "central_mmae"
  )
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L (acMMAE, Cc); ug/mL (total antibody, Cc_tab); ng/mL (unconjugated MMAE, Cc_mmae)"
  )

  compartmentData <- list(
    dar8_central = list(analyte = "DCDS5017A DAR8 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar8_peripheral1 = list(analyte = "DCDS5017A DAR8 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar7_central = list(analyte = "DCDS5017A DAR7 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar7_peripheral1 = list(analyte = "DCDS5017A DAR7 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar6_central = list(analyte = "DCDS5017A DAR6 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar6_peripheral1 = list(analyte = "DCDS5017A DAR6 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar5_central = list(analyte = "DCDS5017A DAR5 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar5_peripheral1 = list(analyte = "DCDS5017A DAR5 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar4_central = list(analyte = "DCDS5017A DAR4 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar4_peripheral1 = list(analyte = "DCDS5017A DAR4 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar3_central = list(analyte = "DCDS5017A DAR3 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar3_peripheral1 = list(analyte = "DCDS5017A DAR3 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar2_central = list(analyte = "DCDS5017A DAR2 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar2_peripheral1 = list(analyte = "DCDS5017A DAR2 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar1_central = list(analyte = "DCDS5017A DAR1 species", units = "nmol", specimen = "plasma", verified = TRUE),
    dar1_peripheral1 = list(analyte = "DCDS5017A DAR1 species", units = "nmol", specimen = "tissue", verified = TRUE),
    dar0_central = list(
      analyte = "anti-CD79b antibody with no MMAE (DAR0)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    dar0_peripheral1 = list(
      analyte = "anti-CD79b antibody with no MMAE (DAR0)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    central_mmae = list(
      analyte = "unconjugated monomethyl auristatin E (MMAE)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_mmae = list(
      analyte = "unconjugated monomethyl auristatin E (MMAE)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 52L,
    n_studies = 3L,
    weight_range = "not reported; the paper assumes 2.5 kg as the typical body weight",
    weight_median = "2.5 kg (assumed typical value, Results)",
    sex_female_pct = NA_real_,
    disease_state = "Healthy (non-tumour-bearing) animals",
    dose_range = "DCDS5017A IV single dose 0.3, 1 or 3 mg/kg (n = 4 per dose) and 3 or 5 mg/kg IV q3w x 4 (n = 10 per dose); unconjugated MMAE IV single dose 0.03 or 0.063 mg/kg (n = 10 per dose)",
    n_observations = "447 total-antibody, 132 acMMAE and 454 unconjugated-MMAE concentrations",
    notes = "32 monkeys received DCDS5017A and 20 received unconjugated MMAE (Lu 2015 Methods 'Pharmacokinetic Data in Cynomolgus Monkeys' and Results). DCDS5017A is the cynomolgus-cross-reactive anti-CD79b surrogate ADC with an mc-vc-PABC-MMAE linker-drug and average DAR ~3.5. Fitted with the ITS-MCPEM algorithm in S-ADAPT 1.57 with below-LLOQ data handled as censored. The model has no covariates: volumes and clearances are absolute values for a monkey of about 2.5 kg, and doses in mg/kg must be converted to nmol with an assumed body weight."
  )

  ini({
    # ADC (all DAR species share Vc, Vp, CLd and the proteolytic Vmax / Km)
    lvmax <- log(0.0407)
    label("Maximum proteolytic degradation rate of each DAR species, Emax (nmol/h)") # Table I 'Emax (1/h)' 0.0407; the rate law Emax*X/(Km*Vc + X) (Eq. 1) requires amount/time
    lkm <- log(0.959)
    label("ADC concentration at half-maximal proteolytic degradation, Km (nmol/L)") # Table I 'Km (nmol/L)' 0.959
    lvc <- log(0.0996)
    label("ADC central volume, Vc (L)") # Table I 'Vc (L)' 0.0996
    lq <- log(0.00197)
    label("ADC intercompartmental clearance, CLd (L/h)") # Table I 'CLd (L/h)' 0.00197
    lvp <- log(0.183)
    label("ADC peripheral volume, Vp (L)") # Table I 'Vp (L)' 0.183

    # Weibull dependence of the per-MMAE deconjugation rate constant on DAR (Eq. 11)
    lkdec_wscale <- log(1.51)
    label("Weibull scale parameter lambda of the DAR-dependent deconjugation rate (unitless)") # Table I 'lambda [Weibull scale parameter]' 1.51
    lkdec_wshape <- log(1.10)
    label("Weibull shape parameter kappa of the DAR-dependent deconjugation rate (unitless)") # Table I 'kappa [Weibull shape parameter]' 1.10
    kdec_wamp <- fixed(1)
    label("Weibull amplitude K of the DAR-dependent deconjugation rate (1/h)") # Results text below Eq. 11: 'K is fixed to 1 h-1'

    # Fractions of released MMAE that reach the systemic MMAE central compartment
    lfm_proteolysis <- log(0.783)
    label("Fraction of MMAE released by proteolytic degradation that enters systemic circulation, Fmp (unitless)") # Table I 'Fmp' 0.783
    lfm_deconjugation <- log(0.0242)
    label("Fraction of MMAE released by deconjugation that enters systemic circulation, Fmd (unitless)") # Table I 'Fmd' 0.0242

    # Unconjugated MMAE
    lcl_mmae <- log(2.27)
    label("Unconjugated MMAE clearance, CLm (L/h)") # Table I 'CLm (L/h)' 2.27
    lvc_mmae <- log(1.00)
    label("Unconjugated MMAE central volume, Vm (L)") # Table I 'Vm (L)' 1.00
    lq_mmae <- log(23.0)
    label("Unconjugated MMAE intercompartmental clearance, CLmd (L/h)") # Table I 'CLmd (L/h)' 23.0
    lvp_mmae <- log(88.4)
    label("Unconjugated MMAE peripheral volume, Vmp (L)") # Table I 'Vmp (L)' 88.4

    # Composition of the dosing solution by hydrophobic interaction chromatography
    fdar0 <- fixed(0.05235)
    label("Fraction of the ADC dose present as DAR0 (unitless)") # Methods 'Model Development': DAR = 0, 5.235%
    fdar1 <- fixed(0.0061)
    label("Fraction of the ADC dose present as DAR1 (unitless)") # Methods 'Model Development': DAR = 1, 0.61%
    fdar2 <- fixed(0.3051)
    label("Fraction of the ADC dose present as DAR2 (unitless)") # Methods 'Model Development': DAR = 2, 30.51%
    fdar3 <- fixed(0)
    label("Fraction of the ADC dose present as DAR3 (unitless)") # Methods 'Model Development': no DAR = 3 listed
    fdar4 <- fixed(0.4835)
    label("Fraction of the ADC dose present as DAR4 (unitless)") # Methods 'Model Development': DAR = 4, 48.35%
    fdar5 <- fixed(0)
    label("Fraction of the ADC dose present as DAR5 (unitless)") # Methods 'Model Development': no DAR = 5 listed
    fdar6 <- fixed(0.12755)
    label("Fraction of the ADC dose present as DAR6 (unitless)") # Methods 'Model Development': DAR = 6, 12.755%
    fdar7 <- fixed(0)
    label("Fraction of the ADC dose present as DAR7 (unitless)") # Methods 'Model Development': no DAR = 7 listed
    fdar8 <- fixed(0.02545)
    label("Fraction of the ADC dose present as DAR8 (unitless)") # Methods 'Model Development': DAR = 8, 2.545%

    # Inter-individual variability: Table I column 'Inter-individual variability'
    # read as log-scale variances (omega^2); the parenthetical is the precision
    # of the estimate, as in the population-mean column.
    etalvmax ~ 0.0945 # Table I Emax IIV 0.0945
    etalvc ~ 0.0118 # Table I Vc IIV 0.0118
    etalvp ~ 0.230 # Table I Vp IIV 0.230
    etalfm_deconjugation ~ 1.78 # Table I Fmd IIV 1.78
    etalcl_mmae ~ 0.333 # Table I CLm IIV 0.333
    etalvp_mmae ~ 0.0670 # Table I Vmp IIV 0.0670

    # Proportional residual error (S-ADAPT SD as a fraction of the prediction)
    propSd <- 0.257
    label("Proportional residual SD of acMMAE (fraction)") # Table I 'sigma acMMAE' 0.257
    propSd_tab <- 0.194
    label("Proportional residual SD of total antibody (fraction)") # Table I 'sigma Ab [Residual error for Tab]' 0.194
    propSd_mmae <- 0.327
    label("Proportional residual SD of unconjugated MMAE (fraction)") # Table I 'sigma MMAE' 0.327
  })
  model({
    mw_tab <- 145239 # g/mol, DCDS5017A molecular weight (Methods 'Antibody Drug Conjugates')
    mw_mmae <- 718 # g/mol, MMAE molecular weight (Methods 'Antibody Drug Conjugates')

    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kdec_wscale <- exp(lkdec_wscale)
    kdec_wshape <- exp(lkdec_wshape)

    fm_proteolysis <- exp(lfm_proteolysis)
    fm_deconjugation <- exp(lfm_deconjugation + etalfm_deconjugation)

    cl_mmae <- exp(lcl_mmae + etalcl_mmae)
    vc_mmae <- exp(lvc_mmae)
    q_mmae <- exp(lq_mmae)
    vp_mmae <- exp(lvp_mmae + etalvp_mmae)

    # Per-MMAE deconjugation rate constant of the DARn species (Eq. 11),
    # evaluated at x = 8 - n. At x = 0 the Weibull density with kappa > 1 is
    # exactly zero, yet Fig. 8 and the Results give 8 * KADC8 = 3.52 1/h; the
    # DAR8 term is therefore evaluated at x = 0.01, the offset that reproduces
    # 3.52 1/h (8 * KADC8 = 3.51 1/h at the Table I values).
    x8 <- 0.01
    kdec8 <- kdec_wamp * kdec_wshape / kdec_wscale * (x8 / kdec_wscale)^(kdec_wshape - 1) * exp(-(x8 / kdec_wscale)^kdec_wshape)
    kdec7 <- kdec_wamp * kdec_wshape / kdec_wscale * (1 / kdec_wscale)^(kdec_wshape - 1) * exp(-(1 / kdec_wscale)^kdec_wshape)
    kdec6 <- kdec_wamp * kdec_wshape / kdec_wscale * (2 / kdec_wscale)^(kdec_wshape - 1) * exp(-(2 / kdec_wscale)^kdec_wshape)
    kdec5 <- kdec_wamp * kdec_wshape / kdec_wscale * (3 / kdec_wscale)^(kdec_wshape - 1) * exp(-(3 / kdec_wscale)^kdec_wshape)
    kdec4 <- kdec_wamp * kdec_wshape / kdec_wscale * (4 / kdec_wscale)^(kdec_wshape - 1) * exp(-(4 / kdec_wscale)^kdec_wshape)
    kdec3 <- kdec_wamp * kdec_wshape / kdec_wscale * (5 / kdec_wscale)^(kdec_wshape - 1) * exp(-(5 / kdec_wscale)^kdec_wshape)
    kdec2 <- kdec_wamp * kdec_wshape / kdec_wscale * (6 / kdec_wscale)^(kdec_wshape - 1) * exp(-(6 / kdec_wscale)^kdec_wshape)
    kdec1 <- kdec_wamp * kdec_wshape / kdec_wscale * (7 / kdec_wscale)^(kdec_wshape - 1) * exp(-(7 / kdec_wscale)^kdec_wshape)

    k12 <- q / vc
    k21 <- q / vp
    kmv <- km * vc # half-saturation amount (nmol)

    # Proteolytic degradation of each species saturates on that species' own
    # amount (Eq. 1, 3, 5).
    prot8 <- vmax * dar8_central / (kmv + dar8_central)
    prot7 <- vmax * dar7_central / (kmv + dar7_central)
    prot6 <- vmax * dar6_central / (kmv + dar6_central)
    prot5 <- vmax * dar5_central / (kmv + dar5_central)
    prot4 <- vmax * dar4_central / (kmv + dar4_central)
    prot3 <- vmax * dar3_central / (kmv + dar3_central)
    prot2 <- vmax * dar2_central / (kmv + dar2_central)
    prot1 <- vmax * dar1_central / (kmv + dar1_central)
    prot0 <- vmax * dar0_central / (kmv + dar0_central)

    # Deconjugation flux out of each species, in species amount per hour; the
    # factor n is the number of MMAE that can come off (Eq. 1, 3, 5).
    dec8 <- 8 * kdec8 * dar8_central
    dec7 <- 7 * kdec7 * dar7_central
    dec6 <- 6 * kdec6 * dar6_central
    dec5 <- 5 * kdec5 * dar5_central
    dec4 <- 4 * kdec4 * dar4_central
    dec3 <- 3 * kdec3 * dar3_central
    dec2 <- 2 * kdec2 * dar2_central
    dec1 <- 1 * kdec1 * dar1_central

    d/dt(dar8_central) <- -k12 * dar8_central + k21 * dar8_peripheral1 - prot8 - dec8
    d/dt(dar8_peripheral1) <- k12 * dar8_central - k21 * dar8_peripheral1
    d/dt(dar7_central) <- -k12 * dar7_central + k21 * dar7_peripheral1 - prot7 + dec8 - dec7
    d/dt(dar7_peripheral1) <- k12 * dar7_central - k21 * dar7_peripheral1
    d/dt(dar6_central) <- -k12 * dar6_central + k21 * dar6_peripheral1 - prot6 + dec7 - dec6
    d/dt(dar6_peripheral1) <- k12 * dar6_central - k21 * dar6_peripheral1
    d/dt(dar5_central) <- -k12 * dar5_central + k21 * dar5_peripheral1 - prot5 + dec6 - dec5
    d/dt(dar5_peripheral1) <- k12 * dar5_central - k21 * dar5_peripheral1
    d/dt(dar4_central) <- -k12 * dar4_central + k21 * dar4_peripheral1 - prot4 + dec5 - dec4
    d/dt(dar4_peripheral1) <- k12 * dar4_central - k21 * dar4_peripheral1
    d/dt(dar3_central) <- -k12 * dar3_central + k21 * dar3_peripheral1 - prot3 + dec4 - dec3
    d/dt(dar3_peripheral1) <- k12 * dar3_central - k21 * dar3_peripheral1
    d/dt(dar2_central) <- -k12 * dar2_central + k21 * dar2_peripheral1 - prot2 + dec3 - dec2
    d/dt(dar2_peripheral1) <- k12 * dar2_central - k21 * dar2_peripheral1
    d/dt(dar1_central) <- -k12 * dar1_central + k21 * dar1_peripheral1 - prot1 + dec2 - dec1
    d/dt(dar1_peripheral1) <- k12 * dar1_central - k21 * dar1_peripheral1
    d/dt(dar0_central) <- -k12 * dar0_central + k21 * dar0_peripheral1 - prot0 + dec1
    d/dt(dar0_peripheral1) <- k12 * dar0_central - k21 * dar0_peripheral1

    # MMAE released by each pathway (Eq. 9, 10): proteolysis of a DARn
    # species frees n MMAE, so DAR0 contributes none; deconjugation frees one
    # MMAE per event.
    amp <- 8 * prot8 + 7 * prot7 + 6 * prot6 + 5 * prot5 + 4 * prot4 + 3 * prot3 + 2 * prot2 + prot1
    amd <- dec8 + dec7 + dec6 + dec5 + dec4 + dec3 + dec2 + dec1

    d/dt(central_mmae) <- fm_proteolysis * amp + fm_deconjugation * amd -
      (cl_mmae + q_mmae) / vc_mmae * central_mmae + q_mmae / vp_mmae * peripheral1_mmae
    d/dt(peripheral1_mmae) <- q_mmae / vc_mmae * central_mmae - q_mmae / vp_mmae * peripheral1_mmae

    # Dose the total ADC amount (nmol of antibody) into every darN_central;
    # f() splits it by the dosing-solution composition.
    f(dar8_central) <- fdar8
    f(dar7_central) <- fdar7
    f(dar6_central) <- fdar6
    f(dar5_central) <- fdar5
    f(dar4_central) <- fdar4
    f(dar3_central) <- fdar3
    f(dar2_central) <- fdar2
    f(dar1_central) <- fdar1
    f(dar0_central) <- fdar0

    tab_central <- dar0_central + dar1_central + dar2_central + dar3_central + dar4_central +
      dar5_central + dar6_central + dar7_central + dar8_central
    acmmae_central <- dar1_central + 2 * dar2_central + 3 * dar3_central + 4 * dar4_central +
      5 * dar5_central + 6 * dar6_central + 7 * dar7_central + 8 * dar8_central

    Cc <- acmmae_central / vc
    Cc_tab <- tab_central / vc * mw_tab / 1e6
    Cc_mmae <- central_mmae / vc_mmae * mw_mmae / 1000

    Cc ~ prop(propSd)
    Cc_tab ~ prop(propSd_tab)
    Cc_mmae ~ prop(propSd_mmae)
  })
}

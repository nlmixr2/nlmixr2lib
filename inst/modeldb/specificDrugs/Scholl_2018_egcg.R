Scholl_2018_egcg <- function() {
  description <- "Two-compartment population PK model for epigallocatechin gallate (EGCG) after oral green tea extract in healthy adults (Scholl 2018), with sequential zero-order dissolution into the depot and lagged first-order absorption, relative bioavailability lowered in ABCC2 c.-24C>T (rs717620) C/C wild-type carriers, apparent clearance lowered linearly per SLCO1B1 c.521T>C (rs4149056) C allele, and combined proportional + additive residual error."
  reference <- paste(
    "Scholl C, Lepper A, Lehr T, Hanke N, Schneider KL, Brockmoller J, Seufferlein T, Stingl JC. (2018).",
    "Population nutrikinetics of green tea extract.",
    "PLoS One 13(2):e0193074. doi:10.1371/journal.pone.0193074"
  )
  vignette <- "Scholl_2018_green_tea_catechins"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "epigallocatechin gallate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "epigallocatechin gallate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "epigallocatechin gallate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SNP_ABCC2_RS717620_HET = list(
      description = "ABCC2 (MRP2) c.-24C>T (rs717620) heterozygote (C/T) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 coded MRP71T as the variant-allele count (C/C = 0, C/T = 1, T/T = 2). The EGCG model applies its effect only to the C/C wild-type stratum, i.e. when SNP_ABCC2_RS717620_HET = 0 and SNP_ABCC2_RS717620_HOM = 0: F1 is multiplied by 0.740 (26 percent lower than in T-allele carriers). Genotype frequencies C/C 58, C/T 22, T/T 3 of 83 genotyped (Table 4).",
      source_name = "MRP71T"
    ),
    SNP_ABCC2_RS717620_HOM = list(
      description = "ABCC2 (MRP2) c.-24C>T (rs717620) homozygous-variant (T/T) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Companion to SNP_ABCC2_RS717620_HET; derive from MRP71T == 2.",
      source_name = "MRP71T"
    ),
    SNP_SLCO1B1_RS4149056_HET = list(
      description = "SLCO1B1 (OATP1B1) c.521T>C (rs4149056) heterozygote (T/C) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Scholl 2018 coded OATP41T as the C-allele count (T/T = 0, C/T = 1, C/C = 2) and entered it linearly: CL/F * (1 + theta * OATP41T), theta = -0.152, i.e. -15 percent for C/T and -30 percent for C/C (Table 5). OATP41T = HET + 2 * HOM. Genotype frequencies T/T 54, C/T 27, C/C 3 (Table 4).",
      source_name = "OATP41T"
    ),
    SNP_SLCO1B1_RS4149056_HOM = list(
      description = "SLCO1B1 (OATP1B1) c.521T>C (rs4149056) homozygous-variant (C/C) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Companion to SNP_SLCO1B1_RS4149056_HET; derive from OATP41T == 2.",
      source_name = "OATP41T"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 84L,
    n_studies = 1L,
    age_range = "19-49 years",
    age_mean = "24.5 years (SD 3.9)",
    weight_mean = "68.0 kg (SD 11.7)",
    sex_female_pct = 64.3,
    race_ethnicity = "Western European ancestry (all participants)",
    disease_state = "Healthy volunteers (BMI 18-30 kg/m2; no acute or chronic disease; no regular medication except oral contraceptives).",
    dose_range = "Green tea extract capsules standardised to 150 mg EGCG (plus 18 mg EGC, 35 mg ECG, 10 mg EC) taken orally twice daily with a meal for 4 days, then one capsule after an overnight fast on day 5; PK sampled on day 5 pre-dose and 0.5-9 h post-dose.",
    regions = "Germany (University of Ulm)",
    notes = "Scholl 2018 Table 1 (84 analysed of 89 enrolled; 30 male, 54 female, 37 of the women on oral contraceptives). NCT01360320. The three catechin models (EGCG, EGC, ECG) were developed independently on the same subjects."
  )

  ini({
    # Structural parameters -- Scholl 2018 Table 3, EGCG column (final
    # estimates). Structure from the S6 File (Modeling Supplement) EGCG
    # NONMEM control stream: ADVAN6 with explicit $DES, zero-order input of
    # duration D1 into the depot (RATE = -2 in the data), first-order K12
    # from depot to central, lag time and F1 on the depot.
    lka <- log(1.86); label("First-order absorption rate constant ka, depot to central (1/h)") # Table 3 EGCG ka = 1.86 1/h (S6 File K12 = THETA(1))
    lcl <- log(148); label("Apparent clearance CL/F, SLCO1B1 rs4149056 T/T (L/h)") # Table 3 EGCG CL/F = 148 L/h
    lvc <- log(334); label("Apparent central volume Vcentral/F (L)") # Table 3 EGCG Vcentral/F = 334 L
    lq <- log(59.7); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 EGCG Q/F = 59.7 L/h
    lvp <- log(1260); label("Apparent peripheral volume Vperipheral/F (L)") # Table 3 EGCG Vperipheral/F = 1260 L
    ld1 <- log(0.301); label("Zero-order dissolution duration D1 into the depot (h)") # Table 3 EGCG D1 = 0.301 h
    ltlag <- log(0.327); label("Absorption lag time (h)") # Table 3 EGCG lag time = 0.327 h
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, ABCC2 rs717620 T-allele carriers (unitless)") # Table 3 EGCG F = 1.00 fixed; S6 File $THETA (1 FIX) ;F

    # Covariate effects -- Table 3 EGCG rows; implementation from the S6 File
    # EGCG $PK block.
    e_snp_slco1b1_rs4149056_cl <- -0.152; label("Fractional change in CL/F per SLCO1B1 rs4149056 C allele (unitless)") # Table 3 OATP1B1_rs4149056 = -0.152; S6 File CL = THETA(2)*(1+OATP41EFF*OATP41T)
    e_snp_abcc2_rs717620_cc_fdepot <- 0.740; label("F1 multiplier for ABCC2 rs717620 C/C wild-type carriers (unitless)") # Table 3 MRP2_rs717620_EGCG = 0.740; S6 File IF (MRP71T.EQ.0) MRP71EFF=THETA(10)

    # IIV -- Table 3 EGCG IIV rows (%CV). Encoded as omega^2 = (CV/100)^2,
    # the convention the paper states (Methods: the exponential-model omega
    # approximates the CV) and that the S6 File $OMEGA initial estimates of
    # the ECG model reproduce (see vignette, Assumptions and deviations).
    # Diagonal $OMEGA in the S6 File EGCG stream.
    etalka ~ 0.6241 # Table 3 EGCG IIV ka 79 %CV -> 0.79^2
    etalcl ~ 0.0625 # Table 3 EGCG IIV CL/F 25 %CV -> 0.25^2
    etald1 ~ 1.0609 # Table 3 EGCG IIV D1 103 %CV -> 1.03^2
    etalfdepot ~ 0.2116 # Table 3 EGCG IIV F 46 %CV -> 0.46^2
    etaltlag ~ 0.1024 # Table 3 EGCG IIV lag time 32 %CV -> 0.32^2

    # Residual error -- S6 File $ERROR Y = IPRED*(1+EPS(1)) + EPS(2).
    propSd <- 0.17; label("Proportional residual error (fraction)") # Table 3 EGCG prop. residual variability 17 %
    addSd <- 5; label("Additive residual error (ug/L)") # Table 3 EGCG add. residual variability 5 ug/L
  })

  model({
    # SLCO1B1 c.521T>C C-allele count (OATP41T in the S6 File).
    oatp41t <- SNP_SLCO1B1_RS4149056_HET + 2 * SNP_SLCO1B1_RS4149056_HOM
    # ABCC2 c.-24C>T C/C wild-type indicator (MRP71T == 0 in the S6 File).
    mrp71_cc <- 1 - SNP_ABCC2_RS717620_HET - SNP_ABCC2_RS717620_HOM

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (1 + e_snp_slco1b1_rs4149056_cl * oatp41t)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot + etalfdepot) * e_snp_abcc2_rs717620_cc_fdepot^mrp71_cc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Doses enter the depot as a zero-order infusion of duration D1 (dose
    # records need rate = -2), delayed by the lag time.
    dur(depot) <- d1
    alag(depot) <- tlag
    f(depot) <- fdepot

    # Amounts in mg, volumes in L; S2 = V2/1000 gives ug/L.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}

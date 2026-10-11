CherkaouiRbati_2023_dsm265 <- function() {
  description <- paste(
    "Semi-mechanistic chemoprophylaxis PK/PD (QSP) model of the antimalarial",
    "dihydroorotate-dehydrogenase (DHODH) inhibitor DSM265 against",
    "Plasmodium falciparum (Cherkaoui-Rbati 2023). The model mimics the",
    "parasite life cycle in two coupled pools: a liver (hepatic-schizont)",
    "stage and a blood stage, each growing exponentially and killed by an",
    "Emax/EC50 Hill turnover term driven by the DSM265 dried-blood-spot",
    "concentration. Merozoites are released from liver to blood through a",
    "sigmoidal maturation window centred at T50 ~ 6 days. PK is a",
    "two-compartment model with zero-order absorption, a lag time, linear",
    "elimination, and a dose-level effect on relative bioavailability;",
    "individual PK parameters feed the PD layer as drivers. Parameters were",
    "estimated sequentially (induced-blood-stage malaria, sporozoite human",
    "challenge, and published placebo-challenge studies) and assembled into",
    "the final chemoprophylaxis model. Liver- and blood-stage growth rates",
    "and EC50 are assumed 100% correlated within an individual, and liver",
    "Emax, Hill and turnover are fixed to the blood-stage values (not",
    "separately identifiable). Reproduces the Day-28 protection success",
    "rates in Cherkaoui-Rbati 2023 Figure 3.",
    sep = " "
  )
  reference <- paste(
    "Cherkaoui-Rbati MH, Andenmatten N, Burgert L, Egbelowo OF, Fendel R,",
    "Fornari C, Gabel M, Ward J, Moehrle JJ, Gobeau N.",
    "A pharmacokinetic-pharmacodynamic model for chemoprotective agents",
    "against malaria.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(1):50-61.",
    "doi:10.1002/psp4.12875.",
    sep = " "
  )
  vignette <- "CherkaouiRbati_2023_dsm265"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Life-cycle parasite pools (liver / blood stages) and their drug-kill
  # turnover states are paper-mechanistic, not canonical compartments.
  paper_specific_compartments <- c(
    "parasite_liver",
    "parasite_blood",
    "kill_liver",
    "kill_blood"
  )

  # The drug dose enters `central` (zero-order absorption); the sporozoite
  # inoculum enters `parasite_liver`.
  dosing <- c("central", "parasite_liver")

  compartmentData <- list(
    central = list(
      analyte = "DSM265",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "DSM265",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    parasite_liver = list(
      analyte = "Plasmodium falciparum parasites (liver / hepatic-schizont stage)",
      units = "parasites",
      specimen = "tissue",
      verified = TRUE
    ),
    kill_liver = list(
      analyte = "DSM265 drug-induced liver-stage kill rate (turnover state)",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    parasite_blood = list(
      analyte = "Plasmodium falciparum parasites (blood stage)",
      units = "parasites",
      specimen = "blood cell",
      verified = TRUE
    ),
    kill_blood = list(
      analyte = "DSM265 drug-induced blood-stage kill rate (turnover state)",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 190,
    n_studies = 15,
    disease_state = paste(
      "Healthy volunteers in induced-blood-stage malaria (IBSM) and",
      "sporozoite human-challenge (spz HuCh) studies, plus patients with",
      "uncomplicated Plasmodium falciparum malaria (phase IIa)"
    ),
    dose_range = "25-1200 mg DSM265 single oral dose; spz-HuCh used 400 mg",
    notes = paste(
      "Data pooled from 15 studies: 1 first-in-human single-ascending-dose",
      "(55), 2 IBSM (15), 1 phase IIa in Peru (24), 2 spz-HuCh (39,",
      "including 10 placebo), and 57 placebo subjects from 9 published",
      "sporozoite-challenge studies (Coffeng et al. 2017). See Table 1."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # PK: two-compartment, zero-order absorption with lag, linear elimination
    # (Cherkaoui-Rbati 2023 Table 2; final PK model). DSM265 dried-blood-spot
    # concentration in ug/mL. F anchored at 1 with a dose-level effect on
    # relative bioavailability.
    # -----------------------------------------------------------------------
    lcl <- log(0.476); label("Clearance (L/h)")                              # Table 2: CL = 0.476 L/h (RSE 2.71%)
    lvc <- log(8.58); label("Central volume (L)")                            # Table 2: Vc = 8.58 L (RSE 12.8%)
    lq <- log(37.1); label("Inter-compartmental clearance (L/h)")            # Table 2: Q1 = 37.1 L/h (RSE 4.61%)
    lvp <- log(57.3); label("Peripheral volume (L)")                         # Table 2: Vp,1 = 57.3 L (RSE 3.36%)
    ld1 <- log(2.85); label("Zero-order absorption duration (h)")            # Table 2: Tk,0 = 2.85 h (RSE 5.84%)
    ltlag <- log(0.147); label("Absorption lag time (h)")                    # Table 2: Tlag,1 = 0.147 h (RSE 12.7%)
    lfrel <- fixed(log(1)); label("Relative bioavailability (fraction)")     # Table 2: F = 1 (FIX)
    e_dose_frel <- -0.102; label("Dose-level effect on relative bioavailability (power on dose/400 mg)") # Table 2: beta_F,Dose = -0.102 (RSE 27%); transform log(dose/400)

    # -----------------------------------------------------------------------
    # PD: liver- and blood-stage parasite growth and drug kill
    # (Cherkaoui-Rbati 2023 Table 3; final chemoprophylaxis model). Growth
    # rates and EC50 are 100% correlated between stages (shared rho etas);
    # liver Emax/Hill/kin fixed to blood values (not separately identifiable).
    # -----------------------------------------------------------------------
    lkgrow_liver <- fixed(log(0.0716)); label("Liver-stage parasite growth rate (1/h)")      # Table 3: GR_L = 0.0716 1/h (fixed; derived from ~30000 merozoites per hepatocyte over ~6 days)
    lkgrow_blood <- fixed(log(0.0624)); label("Blood-stage parasite growth rate (1/h)")      # Table 3: GR_B = 0.0624 1/h (fixed from step 3; Table S2)
    lemax_liver <- fixed(log(0.205)); label("Liver-stage maximum kill rate (1/h)")           # Table 3: Emax,L = 0.205 1/h (fixed to Emax,B)
    lemax_blood <- fixed(log(0.205)); label("Blood-stage maximum kill rate (1/h)")           # Table 3: Emax,B = 0.205 1/h (fixed from step 2; Table S1)
    lec50_liver <- log(1.7); label("Liver-stage EC50, DBS concentration (ug/mL)")            # Table 3: EC50,L = 1.7 ug/mL (RSE 3.99%; the only typical estimated in the final model)
    lec50_blood <- fixed(log(1.25)); label("Blood-stage EC50, DBS concentration (ug/mL)")    # Table 3: EC50,B = 1.25 ug/mL (fixed from step 2; Table S1)
    lhill_liver <- fixed(log(10)); label("Liver-stage Hill coefficient (unitless)")          # Table 3: h_L = 10 (fixed to h_B)
    lhill_blood <- fixed(log(10)); label("Blood-stage Hill coefficient (unitless)")          # Table 3: h_B = 10 (fixed; chosen between 1 and 10 by AIC, Table S1)
    lkin <- fixed(log(0.0771)); label("Kill-effect turnover rate (1/h)")                     # Table 3: k_in = 0.0771 1/h (fixed from step 2; shared liver/blood)
    lklb <- fixed(log(6)); label("Merozoite liver-to-blood invasion rate (1/h)")             # Table 3: k_LB = 6 1/h (fixed; mean transfer time 10 min)
    lt50 <- fixed(log(144)); label("Time of half hepatocyte burst after infection (h)")      # Table 3: T50 = 144 h (fixed; midpoint of the 5-7 day liver stage)
    lsigtliver <- fixed(log(0.1)); label("Spread of the merozoite release window (h)")       # Table 3: sigma_LB = 0.1 h (fixed; >98% released within Day 6 +/- 30 min)
    lfinc <- fixed(log(0.00119)); label("Fraction of inoculated sporozoites invading hepatocytes (fraction)") # Table 3: F_inc = 0.00119 (fixed from step 4; Table S3)

    # Per-individual correlation scalings (typical 1): a single rho draw
    # scales both stage growth rates, another scales both stage EC50s, so
    # liver and blood are 100% correlated within an individual (Table 3:
    # rho_GRL,GRB = 1 and rho_EC50L,EC50B = 1).
    lrhogr <- fixed(log(1)); label("Growth-rate correlation scaling (fraction)")             # Table 3: rho_GRL,GRB = 1 (fixed); IIV carried on this scaling
    lrhoec50 <- fixed(log(1)); label("EC50 correlation scaling (fraction)")                  # Table 3: rho_EC50L,EC50B = 1 (fixed); IIV carried on this scaling

    # -----------------------------------------------------------------------
    # IIV (Table 2 for PK, Table 3 for PD; omegas reported as SD).
    # -----------------------------------------------------------------------
    etalcl ~ 0.0620    # Table 2: IIV CL = 0.249 SD; 0.249^2
    etalvc ~ 1.1236    # Table 2: IIV Vc = 1.06 SD; 1.06^2
    etalvp ~ 0.0790    # Table 2: IIV Vp,1 = 0.281 SD; 0.281^2
    etald1 ~ 0.3170    # Table 2: IIV Tk,0 = 0.563 SD; 0.563^2
    etaltlag ~ 0.5184  # Table 2: IIV Tlag,1 = 0.72 SD; 0.72^2
    etalfrel ~ 0.0437  # Table 2: IIV F = 0.209 SD; 0.209^2
    etalfinc ~ fixed(0.0202)   # Table 3: IIV F_inc = 0.142 SD, held from step 4; 0.142^2
    etalrhogr ~ fixed(0.0139)  # Table 3: IIV GR = 0.118 SD, held from step 3; 0.118^2
    etalrhoec50 ~ 0.1971       # Table 3: IIV EC50 = 0.444 SD; 0.444^2

    # -----------------------------------------------------------------------
    # Residual error.
    # -----------------------------------------------------------------------
    propSd <- 0.179; label("Proportional residual error on DSM265 concentration (fraction)") # Table 2: proportional error = 0.179 (RSE 2.1%)
    addSd <- 2.66; label("Additive residual error on log blood parasitemia (ln p/mL)")        # Table 3: additive error = 2.66 ln(p/mL) (RSE 11.7%)
  })

  model({
    # Fixed physiological constants (Cherkaoui-Rbati 2023 Methods):
    #   VBlood = total blood volume, 5 L = 5000 mL (standard human value).
    #   spz_per_bite = 640 sporozoites per mosquito bite (dose parasite_liver
    #     in sporozoites directly; for bite dosing supply bites * 640).
    VBlood <- 5000

    # 1. Individual PK parameters (Table 2). Relative bioavailability carries a
    #    dose-level power effect read from the central dose amount.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)
    frel <- exp(lfrel + etalfrel) * (podo(central) / 400)^e_dose_frel

    # 2. Individual PD parameters. A single rho draw scales both stage growth
    #    rates; a second scales both stage EC50s (100% within-subject
    #    correlation). Finc carries its own IIV.
    rhogr <- exp(lrhogr + etalrhogr)
    rhoec50 <- exp(lrhoec50 + etalrhoec50)
    kgrow_liver <- exp(lkgrow_liver) * rhogr
    kgrow_blood <- exp(lkgrow_blood) * rhogr
    ec50_liver <- exp(lec50_liver) * rhoec50
    ec50_blood <- exp(lec50_blood) * rhoec50
    emax_liver <- exp(lemax_liver)
    emax_blood <- exp(lemax_blood)
    hill_liver <- exp(lhill_liver)
    hill_blood <- exp(lhill_blood)
    kin <- exp(lkin)
    klb <- exp(lklb)
    t50 <- exp(lt50)
    sigtliver <- exp(lsigtliver)
    finc <- exp(lfinc + etalfinc)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system.
    #    PK (central DBS concentration in ug/mL).
    Cc <- central / vc
    # Positive-guarded concentration for the large-exponent Hill term.
    cceff <- Cc + 1e-9

    #    Sigmoidal fraction of hepatocytes burst since infection (infection is
    #    at absolute t = 0; drug is administered at negative times relative to
    #    infection, matching the source data convention).
    trlb <- expit((t - t50) / sigtliver)
    # Merozoite release flux from liver to blood (parasites/h).
    liver_to_blood <- klb * trlb * parasite_liver

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(parasite_liver) <- parasite_liver * (kgrow_liver - kill_liver) - liver_to_blood
    d/dt(kill_liver) <- kin * (emax_liver * cceff^hill_liver / (ec50_liver^hill_liver + cceff^hill_liver) - kill_liver)

    d/dt(parasite_blood) <- parasite_blood * (kgrow_blood - kill_blood) + liver_to_blood
    d/dt(kill_blood) <- kin * (emax_blood * cceff^hill_blood / (ec50_blood^hill_blood + cceff^hill_blood) - kill_blood)

    # 5. Dosing adjustments. Drug: zero-order duration into central with lag
    #    and dose-dependent relative bioavailability. Sporozoite inoculum:
    #    fraction Finc of the injected sporozoites invade hepatocytes.
    dur(central) <- d1
    alag(central) <- tlag
    f(central) <- frel
    f(parasite_liver) <- finc

    # 6. Observations: DSM265 concentration and log blood parasitemia
    #    (parasites per mL). The paper fits additive error on ln parasitemia.
    lparasitemia_blood <- log(parasite_blood / VBlood + 1e-9)

    Cc ~ prop(propSd)
    lparasitemia_blood ~ add(addSd)
  })
}

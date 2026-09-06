Li_2017_naproxen_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment population PK model for naproxen (NPX) after",
    "intraperitoneal dosing in male and female Lewis rats with collagen-induced arthritis",
    "(CIA, a model of rheumatoid arthritis) and in healthy controls (Li 2017). Absorption from",
    "the i.p. site is first order with a bioavailability fixed at 0.9 from literature i.v. rat",
    "data. The distinguishing feature is that every disposition process operates on UNBOUND",
    "drug while both compartments hold TOTAL drug: saturable albumin binding is solved",
    "algebraically at each time point, in plasma and again in the tissue interstitial fluid",
    "(ISF), and the resulting unbound concentrations drive elimination and distribution. The",
    "binding submodel is a Langmuir high-affinity site plus, after the paper's own",
    "Ka2 * Cup << 1 approximation, a linear low-affinity arm; both capacities are proportional",
    "to the albumin concentration, so the model reproduces the dose-dependent (nonlinear) PK",
    "of naproxen and the hypoalbuminaemia of arthritis from one mechanism. ISF albumin is a",
    "fixed fraction of plasma albumin (E/P), higher in arthritis because inflammation raises",
    "microvascular permeability. Binding constants were estimated from separate ultrafiltration",
    "data and fixed into the PK model; the PK data of all four groups were then fitted jointly",
    "by naive pooling in ADAPT 5, so the model carries no between-subject variability."
  )
  reference <- paste(
    "Li X, DuBois DC, Almon RR, Jusko WJ. Effect of Disease-Related Changes in Plasma Albumin",
    "on the Pharmacokinetics of Naproxen in Male and Female Arthritic Rats.",
    "Drug Metab Dispos. 2017;45(5):476-483. doi:10.1124/dmd.116.074500."
  )
  vignette <- "Li_2017_naproxen_rat"
  units <- list(time = "h", dosing = "ug/kg", concentration = "ug/mL")

  covariateData <- list(
    ALB = list(
      description        = "Plasma albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Li 2017 reports albumin in umol/L as `Pt` (Table 1, marked Fixed; measured by ELISA",
        "and shown in Fig. 2): 347 umol/L in CIA females, 550 in healthy females, 282 in CIA",
        "males and 422 in healthy males. Converted to the canonical SI unit g/L inside the",
        "model with albumin molecular weight 66500 g/mol (the conversion and molecular weight",
        "required by the ALB register entry, precedent `Fauchet_2015_lopinavir_unbound.R`), so",
        "the four group values are 23.1, 36.6, 18.8 and 28.1 g/L respectively. Albumin is the",
        "load-bearing covariate of this paper: it sets both the saturable binding capacity",
        "(n1 * Pt) and the linear non-saturable arm (n2 * Pt * Ka2) in plasma, and the same two",
        "quantities scaled by E/P in the tissue interstitial fluid."
      ),
      source_name        = "Pt"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Sex enters only through the plasma-protein-binding constants: Li 2017 Table 1",
        "estimates Ka1 and Ka2 separately in each of the four sex-by-disease groups, and the",
        "Discussion records that a single shared pair 'produced less satisfactory overall",
        "fittings'. Sex does NOT modify any disposition parameter -- Table 4 assigns CL, CLd",
        "and Vt by disease only, and the Results state that NCA 'did not indicate any",
        "significant differences in the PK parameters between female and male rats'."
      ),
      source_name        = "sex group"
    ),
    DIS_CIA = list(
      description        = "Collagen-induced arthritis indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy control rat)",
      notes              = paste(
        "1 = rat with collagen-induced arthritis studied at peak disease (day 16 post",
        "induction in females, day 21 in males); 0 = age- and sex-matched healthy control.",
        "Disease selects the unbound plasma clearance, the unbound distribution clearance and",
        "the peripheral volume (Li 2017 Table 4, which reports each of these as an Arthritic /",
        "Healthy pair), the ISF-to-plasma albumin ratio E/P (Discussion: 0.9 in CIA vs 0.5 in",
        "healthy rats), and -- jointly with SEXF -- the two albumin association constants of",
        "Table 1. The absorption rate constant and the central volume are shared across all",
        "groups."
      ),
      source_name        = "CIA"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "naproxen",
      units    = "ug/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "naproxen (total, bound plus unbound)",
      units    = "ug/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "naproxen (total, bound plus unbound)",
      units    = "ug/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species        = "rat (Lewis)",
    sex            = "male and female",
    n_subjects     = "3 CIA rats and 4 healthy rats sampled at each PK time point; 8 CIA and 8 healthy females plus 4 CIA and 4 healthy males for the albumin and protein-binding studies",
    n_studies      = 1L,
    age_range      = "5-8 weeks old at purchase, age-matched within each sex group at the time of the PK studies",
    weight_range   = "approximately 110-160 g (females) and 170-220 g (males)",
    disease_state  = "Collagen-induced arthritis (Chondrex protocol) studied at peak disease -- day 16 post induction in females, day 21 in males; approximately 80 percent of females and 60 percent of males developed arthritis in one or both hind paws. Healthy age-matched Lewis rats served as controls.",
    dose_range     = "Single intraperitoneal bolus of sodium naproxen in phosphate-buffered saline (pH 8), 1 mL/kg: 11, 27.5 or 55 mg/kg sodium naproxen (equivalent to 10, 25 or 50 mg/kg naproxen) in CIA rats and 55 mg/kg sodium naproxen (50 mg/kg naproxen) in healthy rats",
    regions        = "United States (State University of New York at Buffalo)",
    notes          = paste(
      "Serial saphenous-vein sampling at 15, 30 and 45 min and 1, 2, 4, 6, 9, 12 and 24 h",
      "post dose; naproxen quantified in plasma and in ultrafiltrate by LC-MS/MS. Plasma",
      "protein binding was measured by ultrafiltration over 2-500 ug/mL in pooled plasma from",
      "each group. All protein-binding and PK data were naive-pooled before analysis; the",
      "binding profiles were fitted first and the estimated binding constants were then fixed",
      "in the PK model. Fitting used the maximum likelihood algorithm in ADAPT 5. Because the",
      "analysis was naive-pooled, the published model has no between-subject variability."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Disposition of UNBOUND naproxen -- Li 2017 Table 4. Clearances are
    # in mL/h/kg and volumes in mL/kg, so amounts are per kg body weight
    # (ug/kg) and concentrations are ug/mL. Table 4 assigns CL, CLd and
    # Vt separately to arthritic and healthy rats and shares ka and Vp.
    # ------------------------------------------------------------------
    lka <- log(0.814); label("First-order absorption rate constant from the intraperitoneal site (1/h)")     # Table 4: ka = 0.814 1/h (CV 7.54%)

    lcl_cia     <- log(1370); label("Unbound plasma clearance in arthritic rats (mL/h/kg)")                  # Table 4: CL,Arthritic = 1370 mL/h/kg (CV 4.73%)
    lcl_healthy <- log(1879); label("Unbound plasma clearance in healthy rats (mL/h/kg)")                    # Table 4: CL,Healthy = 1879 mL/h/kg (CV 9.26%)

    lq_cia      <- log(647.2); label("Unbound distribution clearance in arthritic rats (mL/h/kg)")           # Table 4: CLd,Arthritic = 647.2 mL/h/kg (CV 18.61%)
    lq_healthy  <- log(1371);  label("Unbound distribution clearance in healthy rats (mL/h/kg)")             # Table 4: CLd,Healthy = 1371 mL/h/kg (CV 45.23%)

    lvp_cia     <- log(140.7); label("Peripheral (tissue interstitial fluid) distribution volume of total drug in arthritic rats (mL/kg)")  # Table 4: Vt,Arthritic = 140.7 mL/kg (CV 9.27%)
    lvp_healthy <- log(114.7); label("Peripheral (tissue interstitial fluid) distribution volume of total drug in healthy rats (mL/kg)")    # Table 4: Vt,Healthy = 114.7 mL/kg (CV 17.37%)

    # Vp was "difficult to estimate, perhaps owing to the i.p. doses, and
    # thus it was fixed to rat plasma volume to improve model stability"
    # (Discussion); Table 4 marks it Fixed with footnote a citing the
    # physiologic value of Shah and Betts (2012).
    lvc     <- fixed(log(32.36)); label("Central (plasma) distribution volume of total drug (mL/kg)")        # Table 4: Vp = 32.36 mL/kg (Fixed, rat plasma volume)

    # F was not estimated: "F is the bioavailability of the i.p. dose
    # calculated to be about 0.9 from literature-reported i.v. data in
    # rats (Lauroba et al., 1986)" (Methods, after eq. 7).
    lfdepot <- fixed(log(0.9));   label("Bioavailability of the intraperitoneal dose (fraction)")            # Methods after eq. 7: F = 0.9, fixed from Lauroba 1986

    # ------------------------------------------------------------------
    # Plasma albumin binding -- Li 2017 Table 1. These were estimated
    # from the ultrafiltration data in a separate step and then held
    # fixed in the PK model ("the protein-binding profiles were first
    # fitted and the estimated binding parameters were fixed in the PK
    # model", Methods), hence fixed() throughout.
    #
    # The two-class Langmuir model of eq. 1 is
    #   Cbp = n1*Pt*Ka1*Cup/(1 + Ka1*Cup) + n2*Pt*Ka2*Cup/(1 + Ka2*Cup)
    # and the paper's stated approximation Ka2*Cup << 1 collapses the
    # second class to a linear arm, giving the canonical saturable-plus-
    # linear shape with kaff = Ka1, bmax = n1*Pt (a concentration, in
    # umol/L) and kns = n2*Pt*Ka2 (unitless). Both capacities are built
    # in model() from the live ALB covariate so the paper's mechanism --
    # albumin drives the nonlinearity -- stays exercisable, which is why
    # only the association constants appear here.
    #
    # Ka1 and Ka2 were estimated separately in each of the four sex-by-
    # disease groups; the Discussion records that "attempted modeling of
    # all protein-binding data using single values for Ka1 and Ka2
    # produced less satisfactory overall fittings". Group suffixes are
    # cf/hf/cm/hm, the paper's own labels in the Fig. 2 caption for CIA
    # females, healthy females, CIA males and healthy males.
    # ------------------------------------------------------------------
    kaff_cf <- fixed(0.28); label("Association constant of the first (saturable) albumin binding site, CIA females (L/umol)")      # Table 1: Ka1 = 0.28 /uM (CV 3.53%)
    kaff_hf <- fixed(0.25); label("Association constant of the first (saturable) albumin binding site, healthy females (L/umol)")  # Table 1: Ka1 = 0.25 /uM (CV 3.33%)
    kaff_cm <- fixed(0.26); label("Association constant of the first (saturable) albumin binding site, CIA males (L/umol)")        # Table 1: Ka1 = 0.26 /uM (CV 4.00%)
    kaff_hm <- fixed(0.26); label("Association constant of the first (saturable) albumin binding site, healthy males (L/umol)")    # Table 1: Ka1 = 0.26 /uM (CV 1.50%)

    kaff2_cf <- fixed(0.0041); label("Association constant of the second (non-saturable) albumin binding site, CIA females (L/umol)")      # Table 1: Ka2 = 0.0041 /uM (CV 4.2%)
    kaff2_hf <- fixed(0.0043); label("Association constant of the second (non-saturable) albumin binding site, healthy females (L/umol)")  # Table 1: Ka2 = 0.0043 /uM (CV 4.35%)
    kaff2_cm <- fixed(0.0056); label("Association constant of the second (non-saturable) albumin binding site, CIA males (L/umol)")        # Table 1: Ka2 = 0.0056 /uM (CV 11.75%)
    kaff2_hm <- fixed(0.0054); label("Association constant of the second (non-saturable) albumin binding site, healthy males (L/umol)")    # Table 1: Ka2 = 0.0054 /uM (CV 2.65%)

    # E/P, the ratio of albumin concentration in interstitial fluid to
    # that in plasma. Not estimated: fixed from the literature, and
    # disease-specific because "protein concentrations in tissue
    # interstitial space are known to increase due to the increase of
    # microvascular permeability with inflammation" (Discussion).
    f_alb_isf_cia     <- fixed(0.9); label("Ratio of interstitial-fluid to plasma albumin concentration, arthritic rats (fraction)")  # Discussion: "the E/P value used for CIA rats was about 0.9"
    f_alb_isf_healthy <- fixed(0.5); label("Ratio of interstitial-fluid to plasma albumin concentration, healthy rats (fraction)")    # Discussion: "The average E/P based on the literature is about 0.5 for healthy rats"

    # ------------------------------------------------------------------
    # Residual error. Li 2017 eq. 8 gives the variance model
    # Vi = (s1 + s2*Yi)^2 -- an additive-plus-proportional error on the
    # standard-deviation scale of the model-predicted TOTAL plasma
    # concentration -- but the fitted s1 and s2 are not reported anywhere
    # in the paper: Table 4 lists only the eight structural rows and the
    # supplemental ADAPT code is publisher-hosted and unavailable. Both
    # terms are therefore fixed at zero, which reproduces the published
    # deterministic predictions exactly; see the vignette Errata.
    #
    # There is no between-subject variability to encode: the analysis was
    # naive-pooled ("All protein-binding and PK data were naive-pooled
    # before analysis", Methods), so the published model is a
    # typical-value mechanism by construction rather than by omission.
    # ------------------------------------------------------------------
    addSd  <- fixed(0); label("Additive residual error standard deviation on total plasma naproxen (ug/mL; s1, not reported)")   # eq. 8: Vi = (s1 + s2*Yi)^2; s1 not reported
    propSd <- fixed(0); label("Proportional residual error standard deviation on total plasma naproxen (fraction; s2, not reported)")  # eq. 8: Vi = (s1 + s2*Yi)^2; s2 not reported
  })

  model({
    # ------------------------------------------------------------------
    # 1. Group selection. Ka1 and Ka2 are indexed by sex and disease
    #    (Table 1); CL, CLd, Vt and E/P by disease alone (Table 4 and the
    #    Discussion). ka and Vp are shared by every group.
    # ------------------------------------------------------------------
    kaff <- SEXF * (DIS_CIA * kaff_cf + (1 - DIS_CIA) * kaff_hf) +
      (1 - SEXF) * (DIS_CIA * kaff_cm + (1 - DIS_CIA) * kaff_hm)
    kaff2 <- SEXF * (DIS_CIA * kaff2_cf + (1 - DIS_CIA) * kaff2_hf) +
      (1 - SEXF) * (DIS_CIA * kaff2_cm + (1 - DIS_CIA) * kaff2_hm)
    f_alb_isf <- DIS_CIA * f_alb_isf_cia + (1 - DIS_CIA) * f_alb_isf_healthy

    ka <- exp(lka)
    cl <- DIS_CIA * exp(lcl_cia) + (1 - DIS_CIA) * exp(lcl_healthy)
    q  <- DIS_CIA * exp(lq_cia) + (1 - DIS_CIA) * exp(lq_healthy)
    vp <- DIS_CIA * exp(lvp_cia) + (1 - DIS_CIA) * exp(lvp_healthy)
    vc <- exp(lvc)
    fdepot <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 2. Binding capacities. The canonical covariate ALB is carried in
    #    SI g/L and converted to the umol/L the binding constants of
    #    Table 1 are defined on, using albumin molecular weight
    #    66500 g/mol (ALB register entry; precedent
    #    Fauchet_2015_lopinavir_unbound.R). Table 1 fixes the site
    #    stoichiometries at n1 = 1 and n2 = 4, which appear below as the
    #    literal multipliers on the plasma albumin concentration.
    #
    #    In the interstitial fluid the same affinities and stoichio-
    #    metries are assumed but the albumin concentration is scaled by
    #    E/P: "described by eq. 4 with a difference in ISF and plasma
    #    protein concentrations, where Pt is multiplied by E/P" (Methods).
    # ------------------------------------------------------------------
    albuM <- ALB * 1e6 / 66500

    bmaxPlasma <- 1 * albuM
    knsPlasma  <- 4 * albuM * kaff2
    bmaxIsf    <- 1 * albuM * f_alb_isf
    knsIsf     <- 4 * albuM * f_alb_isf * kaff2

    # ------------------------------------------------------------------
    # 3. Unbound concentrations, solved from the total concentrations by
    #    the positive root of eq. 3 (given as eq. 4):
    #      a*Cu^2 + b*Cu - Ctot = 0
    #      a = Ka1*(1 + n2*Pt*Ka2)
    #      b = n1*Pt*Ka1 + n2*Pt*Ka2 - Ctot*Ka1 + 1
    #      Cu = (-b + sqrt(b^2 + 4*a*Ctot)) / (2*a)
    #    Because a > 0 and Ctot >= 0 the discriminant is never negative
    #    and Cu -> 0 as Ctot -> 0, so the root is numerically safe.
    #
    #    Naproxen molecular weight 230.26 g/mol converts the ug/mL of the
    #    compartments to the umol/L of the binding constants and back;
    #    it is corroborated by the paper's own LC-MS/MS transition,
    #    m/z 229.2 for the [M-H]- ion (Methods, Drug Analysis).
    # ------------------------------------------------------------------
    cpTotUM <- (central / vc) * 1000 / 230.26
    ctTotUM <- (peripheral1 / vp) * 1000 / 230.26

    aPlasma <- kaff * (1 + knsPlasma)
    bPlasma <- bmaxPlasma * kaff + knsPlasma - cpTotUM * kaff + 1
    cuPlasmaUM <- (-bPlasma + sqrt(bPlasma * bPlasma + 4 * aPlasma * cpTotUM)) / (2 * aPlasma)

    aIsf <- kaff * (1 + knsIsf)
    bIsf <- bmaxIsf * kaff + knsIsf - ctTotUM * kaff + 1
    cuIsfUM <- (-bIsf + sqrt(bIsf * bIsf + 4 * aIsf * ctTotUM)) / (2 * aIsf)

    cuPlasma <- cuPlasmaUM * 230.26 / 1000
    cuIsf    <- cuIsfUM * 230.26 / 1000

    # ------------------------------------------------------------------
    # 4. ODE system, eqs. 5-7. The compartments hold TOTAL drug while
    #    every flux is driven by the UNBOUND concentration, which is the
    #    "free hormone hypothesis" structure the paper builds on.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + q * (cuIsf - cuPlasma) - cl * cuPlasma
    d/dt(peripheral1) <-  q * (cuPlasma - cuIsf)

    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 5. Observation. The LC-MS/MS assay and the variance model of eq. 8
    #    are both on the TOTAL plasma concentration, so Cc is total drug;
    #    Cu exposes the unbound plasma concentration that Fig. 7 plots and
    #    that the Results quote as peak unbound concentrations.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cu <- cuPlasma

    Cc ~ add(addSd) + prop(propSd)
  })
}

# nlmixr2lib

# development version

* Add Yoo 2009 cilostazol ([doi:10.1111/j.1365-2125.2009.03558.x](https://doi.org/10.1111/j.1365-2125.2009.03558.x)) -- 104 healthy Korean male volunteers; CYP3A5 and CYP2C19 polymorphism effects on apparent oral clearance.
* Add Berges 2007 enoxaparin ([doi:10.1111/j.1365-2125.2007.02920.x](https://doi.org/10.1111/j.1365-2125.2007.02920.x)) -- elderly inpatients (>75 years) on prophylactic subcutaneous enoxaparin 4000 IU once daily for VTE prophylaxis (PROPHRE.75 study).
* Add Padoin 1998 cephalexin ([doi:10.1128/aac.42.6.1463](https://doi.org/10.1128/aac.42.6.1463)) -- preclinical male Wistar rats with a cephalexin / quinapril oral coadministration DDI (ratifies new `CONMED_QPRL_ORAL` covariate canonical).
* Add van Rongen 2018 metformin ([doi:10.1007/s40272-018-0293-1](https://doi.org/10.1007/s40272-018-0293-1)) -- overweight and obese Caucasian adolescents.
* Add Tarning 2012 artemether and dihydroartemisinin ([doi:10.1186/1475-2875-11-293](https://doi.org/10.1186/1475-2875-11-293)) -- pregnant women with uncomplicated Plasmodium falciparum malaria in Uganda (joint parent + DHA popPK with zero-order dissolution and 6-compartment transit absorption).
* Add Schoemaker 2017 brivaracetam ([doi:10.1007/s00228-017-2230-6](https://doi.org/10.1007/s00228-017-2230-6)) -- paediatric patients with epilepsy aged 1 month to 16 years; ratifies new `CONMED_PB`, `CONMED_CBZ`, and `CONMED_VPA` covariate canonicals.
* Add Nanga 2019 tacrolimus meta-model ([doi:10.1111/bcp.14110](https://doi.org/10.1111/bcp.14110)) -- pooled paediatric and adult solid-organ transplant recipients (n = 281; 201 liver and 80 kidney) on oral tacrolimus; ratifies new `TX_LIVER` and `FORM_SYRUP` covariate canonicals.
* Rewrite Fiedler-Kelly 2019 fremanezumab ([doi:10.1111/bcp.14096](https://doi.org/10.1111/bcp.14096)) to support both IV and SC routes -- the packaged model now carries the route-specific central volume (Vc,IV = 2.98 L FIXED vs Vc,SC = 1.88 L) and the route-specific residual-error structure (IV proportional-only vs SC additive+proportional) from Fiedler-Kelly 2019 Table 2, switched by the `ROUTE_IV` covariate column. Removes the previous spurious allometric weight effect on Vp (the source paper holds Vp FIXED with no weight effect).
* Add Yuan 2019 concizumab ([doi:10.1016/j.ejps.2019.105032](https://doi.org/10.1016/j.ejps.2019.105032)) -- healthy adult males.
* `_pkgdown.yml` is now refreshed in place using `# AUTOGEN:<name>:BEGIN`/`:END` signpost comments rather than a full `yaml::read_yaml` + `yaml::write_yaml` round-trip. The signpost-managed regions are the `specific_drugs` and `ddmore` navbar menus and the top-level `articles:` index; everything else in the file (template, navbar.structure, maintainer-added keys, comments, formatting) is preserved byte-for-byte. The modeldb gains two new cached columns -- `label` (human-readable navbar text like `Ustekinumab (Aguiar 2021)` or `DDMoRe: lidocaine`) and `category` (`specificDrugs` / `ddmore` / `other`) -- so the regeneration costs no extra work beyond the existing per-file md5 cache. Fixes the recurring pkgdown CI breakage where stale entries crept into the `articles:` index when the YAML round-trip captured a transient filesystem state.
* Add Quartino 2016 trastuzumab ([doi:10.1007/s00280-015-2922-5](https://doi.org/10.1007/s00280-015-2922-5)) -- women with HER2-positive early breast cancer (HannaH study; subcutaneous fixed 600 mg q3w and intravenous weight-based regimens).
* Add Diepstraten 2013 propofol ([doi:10.1038/psp.2013.47](https://doi.org/10.1038/psp.2013.47)) -- morbidly obese and nonobese adults, adolescents, and children (pooled meta-analysis across five previously published studies).
* List of models vignette is grouped into General / Specific Drugs / DDMoRe with capped-width tables and uses the navbar-style "Drug (Author Year)" name (`DDMoRe: <drug>` for the parameterless DDMoRe entries). Drug-specific articles in `vignettes/articles/` are now hidden from the auto-generated pkgdown Articles index and reachable only via the navbar dropdowns or direct URL. Validation vignettes that previously titled themselves with the machine basename now use the human form so the page title matches the link.
* Add Diep 2026 donidalorsen ([doi:10.1002/psp4.70206](https://doi.org/10.1002/psp4.70206)) -- pooled healthy volunteers and adult / adolescent patients with hereditary angioedema (ratifies new `DIS_HAE` covariate canonical).
* Add Attarwala 2023 mRNA-3927 ([doi:10.1089/nat.2022.0036](https://doi.org/10.1089/nat.2022.0036)) -- preclinical PCC-deficient mice, juvenile Sprague-Dawley rats, and cynomolgus monkeys (translational semi-mechanistic PK + PK/PD model for an LNP-encapsulated dual mRNA encoding PCCA / PCCB subunits of propionyl-CoA carboxylase; allometrically scalable to humans).
* Add Dunn 2025 tranexamic acid ([doi:10.1002/jcph.70031](https://doi.org/10.1002/jcph.70031)) -- pregnant individuals receiving IV, IM, or oral TXA for prevention or treatment of postpartum hemorrhage.
* Add Parhiz 2024 mRNA-LNP ([doi:10.1016/j.omtn.2024.102175](https://doi.org/10.1016/j.omtn.2024.102175)) -- preclinical mice (C57BL/6, 25 g); whole-body PBPK plus luciferase-expression model for systemically administered firefly-luciferase mRNA in lipid nanoparticles (bare-LNP homogenate-assay parameter set).
* Add Schoenmakers 2025 betamethasone ([doi:10.1002/bcp.70035](https://doi.org/10.1002/bcp.70035)) -- pregnant women with imminent preterm birth, including early-onset pre-eclampsia (ratifies new `DIS_EOPE` covariate canonical for early-onset pre-eclampsia indicator).
* Add Rambiritch 2016 glibenclamide ([doi:10.2147/CPAA.S102676](https://doi.org/10.2147/CPAA.S102676)) -- poorly controlled South African adults with type 2 diabetes mellitus.
* Add Lahu 2010 roflumilast ([doi:10.2165/11536600-000000000-00000](https://doi.org/10.2165/11536600-000000000-00000)) -- adult healthy volunteers and patients with moderate-to-severe COPD (joint parent + roflumilast N-oxide popPK; ratifies new `noxide` registered metabolite suffix and new `DIS_COPD` covariate canonical).
* Add Brill 2014 midazolam ([doi:10.1007/s40262-014-0166-x](https://doi.org/10.1007/s40262-014-0166-x)) -- 20 morbidly obese surgical patients and 12 non-obese healthy volunteers receiving semi-simultaneous oral and intravenous midazolam.
* Add Diao 2014 rFIXFc ([doi:10.1007/s40262-013-0129-7](https://doi.org/10.1007/s40262-013-0129-7)) -- adolescents and adults with severe to moderate haemophilia B (three-compartment IV popPK with estimated body-weight power exponents on CL and V1).
* Add Plock 2014 ferumoxytol ([doi:10.1007/s40262-014-0203-9](https://doi.org/10.1007/s40262-014-0203-9)) -- pooled adult healthy volunteers and CKD-stage-5D haemodialysis patients (two-compartment Michaelis-Menten IV PK with WT and sex covariates on V1; haemodialysis time-varying V1 mechanism documented in the vignette but not enabled in the model file).
* Add Qi 2014 sapropterin ([doi:10.1007/s40262-014-0196-4](https://doi.org/10.1007/s40262-014-0196-4)) -- pediatric and adult patients with phenylketonuria pooled from BioMarin studies PKU-015 (0-6 years) and PKU-004 (9-50 years).
* Add Retlich 2015 linagliptin ([doi:10.1007/s40262-014-0232-4](https://doi.org/10.1007/s40262-014-0232-4)) -- adults with type 2 diabetes mellitus (coupled two-compartment popPK with concentration-dependent DPP-4 binding in both compartments plus a sigmoid Emax PD model relating linagliptin to plasma DPP-4 activity; ratifies new `GGT`, `FPG`, `DPP4_BL_RFU`, `CONMED_METFORMIN`, and `FORM_LINAG_TAB1` covariate canonicals).
* Add Buil-Bruna 2015 lanreotide ([doi:10.1007/s40262-015-0329-4](https://doi.org/10.1007/s40262-015-0329-4)) -- adults with gastroenteropancreatic neuroendocrine tumors.
* Add van Rongen 2016 acetaminophen ([doi:10.1007/s40262-015-0357-0](https://doi.org/10.1007/s40262-015-0357-0)) -- morbidly obese and non-obese adults receiving 2 g IV acetaminophen (parent + glucuronide + sulphate + CYP2E1-oxidation metabolite popPK; ratifies new `cysmer` registered metabolite suffix for combined acetaminophen cysteine + mercapturate).
* Add Gupta 2017 ixazomib ([doi:10.1007/s40262-017-0526-4](https://doi.org/10.1007/s40262-017-0526-4)) -- adult patients with multiple myeloma, lymphoma, solid tumours, or light-chain amyloidosis (pooled phase I-III TOURMALINE-MM1 popPK; three-compartment oral + IV PK with BSA on V4 and time-after-dose-varying residual SD).
* Add Savic 2017 cladribine ([doi:10.1007/s40262-017-0516-6](https://doi.org/10.1007/s40262-017-0516-6)) -- adults with relapsing-remitting multiple sclerosis (three-compartment popPK with renal-clearance scaling on Cockcroft-Gault CRCL and a non-renal-clearance effect of concomitant subcutaneous IFN beta-1a; ratifies new `CONMED_IFNB1A` covariate canonical).
* Add Hendriksen 2013 artesunate ([doi:10.1038/clpt.2013.26](https://doi.org/10.1038/clpt.2013.26)) -- Tanzanian children (7 months to 11 years) with severe Plasmodium falciparum malaria receiving intramuscular artesunate.
* Add Ait-Oudhia 2012 canakinumab ([doi:10.1038/psp.2012.6](https://doi.org/10.1038/psp.2012.6)) -- adults with rheumatoid arthritis (integrated PK / quasi-equilibrium IL-1beta binding / CRP transduction / ACR latent-variable model with logit-transform mapping to ACR20/50/70 response probabilities).
* Add Wu 2012 sirolimus ([doi:10.1038/psp.2012.18](https://doi.org/10.1038/psp.2012.18)) -- adults with advanced solid tumors receiving 1-60 mg/week oral sirolimus (first popPK report describing the nonlinear absorption with a two-compartment + Michaelis-Menten saturable absorption structure; hematocrit power covariate on apparent oral clearance).
* Add de Hoogd 2017 morphine ([doi:10.1007/s40262-017-0544-2](https://doi.org/10.1007/s40262-017-0544-2)) -- 20 morbidly obese adults (post gastric-bypass) and 20 healthy adult volunteers (joint parent + M3G + M6G popPK).
* Add Doldan-Martelli 2013 EGF-IFN chimeric ligand ([doi:10.1038/psp.2013.2](https://doi.org/10.1038/psp.2013.2)) -- in-vitro Daudi / Daudi-EGFR cells (mechanistic two-receptor binding model for an EGF-IFNalpha-2a chimera with selectivity-by-affinity-mutant analysis).
* Add Tetschke 2018 erythropoiesis ([doi:10.3390/pr6090157](https://doi.org/10.3390/pr6090157)) -- healthy adult male volunteers undergoing a single 1-unit blood donation (three-compartment endogenous mixed-effects model for red blood cell regeneration after phlebotomy; ratifies new `THB_MASS` covariate canonical for total hemoglobin mass measured by CO-rebreathing).
* Add Ekhart 2008 carboplatin ([doi:10.1007/s00280-008-0856-x](https://doi.org/10.1007/s00280-008-0856-x)) -- adult cancer patients (underweight to obese) receiving combination chemotherapy.
* Add Aregbe 2012 alvespimycin / 17-DMAG ([doi:10.1007/s00280-012-1859-1](https://doi.org/10.1007/s00280-012-1859-1)) -- adult patients with advanced solid tumors.
* Add Ferron 2013 cabazitaxel ([doi:10.1007/s00280-012-2058-9](https://doi.org/10.1007/s00280-012-2058-9)) -- patients with advanced solid tumors (pooled five Phase I-III studies including TROPIC mCRPC).
* Add Panoilia 2015 bevacizumab ([doi:10.1007/s00280-015-2701-3](https://doi.org/10.1007/s00280-015-2701-3)) -- adults with stage IV metastatic colorectal cancer (TMDD QSS binding model for bevacizumab and free VEGF165 with fixed allometric body-weight scaling).
* Add Varatharajan 2016 daunorubicin ([doi:10.1007/s00280-016-3166-8](https://doi.org/10.1007/s00280-016-3166-8)) -- adult de novo acute myeloid leukaemia patients (joint daunorubicin + daunorubicinol popPK with two-compartment disposition for each component; ratifies new `dol` registered metabolite suffix).
* Add Ozawa 2007 docetaxel ([doi:10.1111/j.1349-7006.2007.00615.x](https://doi.org/10.1111/j.1349-7006.2007.00615.x)) -- Japanese adult cancer patients (62 subjects with breast, NSCLC, head-and-neck, or other solid tumours; coupled three-compartment IV PK and modified Friberg myelosuppression PD with a dexamethasone-induced ANC-bump compartment; AGP power-form effect on the drug-effect slope).
* Add Alsultan 2017 pyrazinamide ([doi:10.1128/AAC.02625-16](https://doi.org/10.1128/AAC.02625-16)) -- adults with drug-susceptible pulmonary tuberculosis (TBTC studies 27 and 28).
* Add Yukawa 1990 phenytoin ([doi:10.1248/cpb.38.1973](https://doi.org/10.1248/cpb.38.1973)) -- Japanese epileptic outpatients on chronic oral phenytoin maintenance therapy (ratifies `CONMED_AED`, `FORM_POWDER`, and `DOSE_PHT_MGKGD` covariate canonicals).
* Add Vezina 2010 valganciclovir ([doi:10.2147/CPAA.S8341](https://doi.org/10.2147/CPAA.S8341)) -- pediatric solid organ transplant recipients receiving valganciclovir prophylaxis for Epstein-Barr virus disease.
* Add Kunisawa 2014 olprinone ([doi:10.2147/CPAA.S50626](https://doi.org/10.2147/CPAA.S50626)) -- healthy adult Japanese male volunteers.
* Add Wang 2021 pertuzumab ([doi:10.1007/s00280-021-04296-0](https://doi.org/10.1007/s00280-021-04296-0)) -- women with HER2-positive early breast cancer (FeDeriCa study, NCT03493854; pooled IV and subcutaneous fixed-dose-combination pertuzumab+trastuzumab cohorts).
* Add Storset 2014 tacrolimus ([doi:10.1111/bcp.12361](https://doi.org/10.1111/bcp.12361)) -- adult kidney-transplant recipients (theory-based plasma-parameterised popPK with FFM allometry, CYP3A5-expresser effects, prednisolone-induced reduction in F, first-day-post-transplant F spike, and saturable haematocrit-dependent RBC binding linking plasma to whole-blood concentration).
* Add Overgaard 2019 semaglutide ([doi:10.1007/s13300-019-0581-y](https://doi.org/10.1007/s13300-019-0581-y)) -- pooled clinical pharmacology cohort of healthy volunteers and adults with type 2 diabetes.
* Add Rovei 1982 theophylline ([doi:10.1111/j.1365-2125.1982.tb02035.x](https://doi.org/10.1111/j.1365-2125.1982.tb02035.x)) -- healthy adult volunteers receiving single oral 125-500 mg theophylline tablets.
* Add Dunlap 2025 tacrolimus ([doi:10.1007/s40262-025-01529-w](https://doi.org/10.1007/s40262-025-01529-w)) -- adult allogeneic hematopoietic cell transplant recipients (CYP3A5 metabolizer phenotype and reduced-intensity conditioning effects on apparent clearance; ratifies new `HCT_COND_RIC` covariate canonical).
* Add Benkali 2010 tacrolimus (Clin Pharmacokinet 2010;49(10):683-92; DOI not on disk) -- stable adult renal transplant recipients switched to once-daily extended-release tacrolimus (Advagraf).
* Add Boer-Perez 2026 piperacillin ([doi:10.1128/aac.00998-25](https://doi.org/10.1128/aac.00998-25)) -- preterm and term neonates with severe infections.
* Add Jonckheere 2019 cefepime ([doi:10.1128/AAC.01552-19](https://doi.org/10.1128/AAC.01552-19)) -- critically ill ICU adults on continuous-infusion cefepime via target-controlled infusion.
* Add Wang 2015 rucaparib ([doi:10.1002/cpdd.176](https://doi.org/10.1002/cpdd.176)) -- adults with advanced solid tumors (Phase 1 first-in-patient study A4991002).
* Add Brown 2017 osimertinib ([doi:10.1111/bcp.13223](https://doi.org/10.1111/bcp.13223)) -- adults with advanced EGFR-mutation-positive non-small cell lung cancer pooled with healthy volunteers (joint osimertinib + AZ5104 metabolite popPK).
* Add Caldes 2009 ganciclovir ([doi:10.1128/AAC.00085-09](https://doi.org/10.1128/AAC.00085-09)) -- adult solid organ transplant recipients (kidney, liver, heart) with cytomegalovirus infection receiving IV ganciclovir followed by oral valganciclovir.
* Add Vezina 2014 valganciclovir ([doi:10.1111/bcp.12343](https://doi.org/10.1111/bcp.12343)) -- paediatric and adult solid organ transplant recipients receiving valganciclovir prophylaxis.
* Add Chen 2021 ganciclovir ([doi:10.1002/jcph.1735](https://doi.org/10.1002/jcph.1735)) -- adult Chinese renal allograft recipients receiving oral valganciclovir for CMV prophylaxis.
* Add Guo 2022 PF-06939999 ([doi:10.1002/psp4.12882](https://doi.org/10.1002/psp4.12882)) -- adults with advanced solid tumors (Phase I FIP study NCT03854227; coupled PK + plasma SDMA indirect-response + Friberg semi-mechanistic platelet model).
* Add Urien 2005 capecitabine ([doi:10.1007/s10928-005-0018-2](https://doi.org/10.1007/s10928-005-0018-2)) -- adults with metastatic cancer.
* Add Garg 2014 pertuzumab ([doi:10.1007/s00280-014-2560-3](https://doi.org/10.1007/s00280-014-2560-3)) -- patients with a variety of HER2-targeted solid tumors.
* Add Chien 2022 imatinib ([doi:10.1007/s00280-022-04454-y](https://doi.org/10.1007/s00280-022-04454-y)) -- healthy Caucasian volunteers receiving a single oral 400 mg dose.
* Add Akbar 2025 voriconazole ([doi:10.1371/journal.pone.0318883](https://doi.org/10.1371/journal.pone.0318883)) -- adult and pediatric Pakistani cancer patients receiving therapeutic drug monitoring of intravenous voriconazole.
* Add Ayyar 2024 givosiran ([doi:10.1016/j.xphs.2023.10.026](https://doi.org/10.1016/j.xphs.2023.10.026)) -- adults with acute hepatic porphyria (mechanistic translational PK with ASGPR-mediated hepatocyte uptake, endolysosomal kinetics, and cytoplasmic RISC loading; ratifies new `asn1` registered metabolite suffix for AS(N-1)3' GalNAc-siRNA).
* Add Salinger 2013 magnesium sulphate ([doi:10.1111/1471-0528.12222](https://doi.org/10.1111/1471-0528.12222)) -- pregnant women with pre-eclampsia receiving IV or IM MgSO4-7H2O for eclampsia prevention.
* Add Chakraborty 2012 canakinumab ([doi:10.2165/11599820-000000000-00000](https://doi.org/10.2165/11599820-000000000-00000)) -- adults and children with cryopyrin-associated periodic syndromes (CAPS) and four other inflammatory cohorts (population PK / IL-1b binding model with quasi-steady-state binding; ratifies new `il1b` registered metabolite suffix).
* Add Rodrigues 2017 oxcarbazepine ([doi:10.1111/bcp.13392](https://doi.org/10.1111/bcp.13392)) -- epileptic children aged 2-12 years (parent oxcarbazepine + active monohydroxy-derivative metabolite popPK with allometric weight scaling, MHD-to-OXC back-transformation, and concomitant enzyme-inducing antiepileptic drug effect; ratifies new `CONMED_EIAED` covariate canonical and `mhd` registered metabolite suffix).
* Add Li 2017 cediranib ([doi:10.1111/bcp.13266](https://doi.org/10.1111/bcp.13266)) -- adult cancer patients on oral cediranib (AZD2171).
* Add Tortorici 2017 alpha-1 proteinase inhibitor ([doi:10.1111/bcp.13358](https://doi.org/10.1111/bcp.13358)) -- adults with severe alpha-1 antitrypsin deficiency.
* Add Valitalo 2017 ketorolac ([doi:10.1111/bcp.13311](https://doi.org/10.1111/bcp.13311)) -- women at delivery, postpartum women, nonpregnant women, and men receiving single IV racemic ketorolac (joint R- / S-enantiomer popPK).
* Add Talke 2018 dexmedetomidine ([doi:10.1111/bcp.13571](https://doi.org/10.1111/bcp.13571)) -- healthy adult volunteers (3-cmt IV PK + effect-compartment sigmoid Emax PD for vasoconstriction).
* Add Huynh 2026 VRC07-523LS ([doi:10.1093/jac/dkaf449](https://doi.org/10.1093/jac/dkaf449)) -- HIV-exposed infants and healthy adults.
* Add Pu 2021 evinacumab ([doi:10.1002/psp4.12711](https://doi.org/10.1002/psp4.12711)) -- healthy volunteers and adult / pediatric patients with homozygous familial hypercholesterolemia.
* Add Zhang 2021 dupilumab ([doi:10.1002/psp4.12667](https://doi.org/10.1002/psp4.12667)) -- adult and adolescent patients with moderate-to-severe asthma (pooled with healthy adults).
* Add Hennig 2015 rifabutin ([doi:10.1128/AAC.01195-15](https://doi.org/10.1128/AAC.01195-15)) -- HIV-infected adults with pulmonary tuberculosis (joint parent + 25-O-desacetyl rifabutin metabolite popPK with allometric weight scaling, sex on V/F, and SLCO1B1 rs11045819 on bioavailability).
* Add Hennig 2013 tobramycin ([doi:10.1007/s40262-013-0036-y](https://doi.org/10.1007/s40262-013-0036-y)) -- adults and children with and without cystic fibrosis.
* Add Yassen 2025 asundexian ([doi:10.1002/psp4.70142](https://doi.org/10.1002/psp4.70142)) -- healthy volunteers and adults at risk for thromboembolic / cardiovascular events (PACIFIC-AF, PACIFIC-STROKE, PACIFIC-AMI Phase II + six Phase I studies; pooled n = 2914).
* Add Jelliffe 2014 digoxin ([doi:10.1097/FTD.0000000000000023](https://doi.org/10.1097/FTD.0000000000000023)) -- adults requiring digoxin therapy across normal-to-anephric renal function.
* Add Aksenov 2018 uricAcid ([doi:10.14814/phy2.13614](https://doi.org/10.14814/phy2.13614)) -- healthy adults and gout patients with hyperuricemia.
* Add Garmann 2017 BAY 81-8973 ([doi:10.1111/hae.13192](https://doi.org/10.1111/hae.13192)) -- patients aged 1-61 years with severe haemophilia A (LEOPOLD I/II/Kids).
* Add Shah 2012 mAb PBPK ([doi:10.1007/s10928-011-9232-2](https://doi.org/10.1007/s10928-011-9232-2)) -- platform PBPK model parameterized for the human (71 kg) reference subject; 15 anatomical tissues with FcRn-mediated recycling.
* Add Kunarajah 2017 doxorubicin ([doi:10.1007/s00280-017-3309-6](https://doi.org/10.1007/s00280-017-3309-6)) -- paediatric oncology patients (joint doxorubicin / doxorubicinol popPK with cardiac troponin I turnover sub-model; ratifies new `PRIOR_ANTHRACYCLINE_DOSE` covariate canonical and `doxol` registered metabolite suffix).
* Add Llanos-Paez 2017 gentamicin ([doi:10.1208/s12248-017-0173-6](https://doi.org/10.1208/s12248-017-0173-6)) -- pediatric oncology patients with febrile or fever-only neutropenia.
* Add Hennig 2015 phenytoin ([doi:10.1002/jcph.417](https://doi.org/10.1002/jcph.417)) -- critically ill children in a paediatric ICU (joint protein-unbound and protein-bound popPK with linear albumin partition).
* Add Bergmann 2014 tacrolimus ([doi:10.1097/FTD.0b013e31829f1ab8](https://doi.org/10.1097/FTD.0b013e31829f1ab8)) -- adult kidney transplant recipients.
* Add Harun 2019 cysticFibrosis ([doi:10.1136/thoraxjnl-2018-211550](https://doi.org/10.1136/thoraxjnl-2018-211550)) -- children with classic cystic fibrosis (sigmoid-Emax disease-progression model of FEV1% predicted from age 5 to 14 years).
* Add Hennig 2007 itraconazole ([doi:10.1111/j.1365-2125.2006.02778.x](https://doi.org/10.1111/j.1365-2125.2006.02778.x)) -- adult cystic fibrosis patients (D-optimal-designed cross-over single-dose 200 mg oral itraconazole with parallel parent / hydroxy-itraconazole metabolite popPK and per-formulation absorption / bioavailability).
* Add Hennig 2008 tobramycin ([doi:10.1111/j.1365-2125.2007.03045.x](https://doi.org/10.1111/j.1365-2125.2007.03045.x)) -- paediatric cystic fibrosis patients on once-daily intravenous tobramycin.
* Add Kirubakaran 2022 tacrolimus ([doi:10.1111/bcp.15566](https://doi.org/10.1111/bcp.15566)) -- adult heart transplant recipients with concomitant azole antifungal therapy.
* Add Hennig 2006 itraconazole ([doi:10.2165/00003088-200645110-00004](https://doi.org/10.2165/00003088-200645110-00004)) -- paediatric cystic-fibrosis and bone-marrow-transplant patients (parent + active metabolite popPK for oral itraconazole and hydroxy-itraconazole).
* Add Lawson 2022 busulfan ([doi:10.1002/psp4.12809](https://doi.org/10.1002/psp4.12809)) -- pediatric hematopoietic stem cell transplant recipients receiving once-daily IV busulfan.
* Add Archary 2019 abacavir ([doi:10.1111/bcp.13998](https://doi.org/10.1111/bcp.13998)) -- severely malnourished HIV-infected children.
* Add Bista 2015 fentanyl (manuscript, journal/DOI not on disk) -- adults with advanced cancer receiving Durogesic transdermal fentanyl matrix patches.
* Add Archary 2018 lopinavir ([doi:10.1097/INF.0000000000001867](https://doi.org/10.1097/INF.0000000000001867)) -- severely malnourished HIV-infected children (1-month to 12-year-olds) on twice-daily oral LPV/rtv with FFM allometric scaling and a total-cholesterol covariate effect on apparent clearance.
* Add Xu 2023 sabatolimab MBG453 ([doi:10.1002/psp4.12962](https://doi.org/10.1002/psp4.12962)) -- adults with advanced solid tumors or hematologic malignancies (AML, MDS, CMML).
* Add Goel 2016 sonidegib ([doi:10.1007/s00280-016-2982-1](https://doi.org/10.1007/s00280-016-2982-1)) -- healthy subjects and adults with advanced solid tumors or basal cell carcinoma.
* Add Tikiso 2021 abacavir ([doi:10.1111/bcp.14984](https://doi.org/10.1111/bcp.14984)) -- HIV-infected African children on abacavir-containing antiretroviral therapy.
* Add Stein 2019 tisagenlecleucel ([doi:10.1002/psp4.12388](https://doi.org/10.1002/psp4.12388)) -- pediatric and young adult patients with relapsed or refractory B-cell ALL (CAR-T cellular kinetic model: exponential expansion to Tmax followed by biexponential effector / memory-cell decline).
* Add van Wijk 2019 paracetamol ([doi:10.1038/s41598-019-38530-w](https://doi.org/10.1038/s41598-019-38530-w)) [DDMODEL00000294] -- preclinical zebrafish (Danio rerio) larvae 3-5 dpf under continuous 1 mM paracetamol bath exposure.
* Add Llanos-Paez 2017 gentamicin ([doi:10.1128/AAC.00205-17](https://doi.org/10.1128/AAC.00205-17)) -- pediatric oncology patients with febrile neutropenia.
* Add Llanos-Paez 2020 gentamicin ([doi:10.1128/AAC.01730-19](https://doi.org/10.1128/AAC.01730-19)) -- pediatric oncology and nononcology patients.
* Add Delattre 2010 amikacin ([doi:10.1097/FTD.0b013e3181f675c2](https://doi.org/10.1097/FTD.0b013e3181f675c2)) -- critically ill adults with severe sepsis or septic shock.
* Add Diep 2022 eplontersen ([doi:10.1111/bcp.15468](https://doi.org/10.1111/bcp.15468)) -- healthy volunteers in two phase 1 studies (two-compartment popPK with site-specific SC absorption and indirect-response PD on serum transthyretin).
* Add Laffont 2024 nalmefene ([doi:10.3389/fpsyt.2024.1399803](https://doi.org/10.3389/fpsyt.2024.1399803)) -- healthy adult volunteers receiving intranasal nalmefene HCl.
* Add Laffont 2024 naloxone ([doi:10.3389/fpsyt.2024.1399803](https://doi.org/10.3389/fpsyt.2024.1399803)) -- healthy adult volunteers receiving intranasal naloxone HCl.
* Add Cao 2013 mPBPK 12-mAb cohort ([doi:10.1007/s10928-013-9332-2](https://doi.org/10.1007/s10928-013-9332-2)) -- second-generation minimal PBPK structural model packaged as one entry per fit: 10 human-mAb entries (adecatumumab, mepolizumab, gevokizumab, GNbAC1, MEDI528, tefibazumab, PAmAb, PRO95780, siltuximab, visilizumab) under specificDrugs/, plus 2 preclinical mouse entries (mab7E3, mab8C2) under pharmacokinetics/.
* Add Sadouki 2025 meropenem / gentamicin / ciprofloxacin ([doi:10.1038/s41598-025-29354-y](https://doi.org/10.1038/s41598-025-29354-y)) -- in-vitro pharmacodynamic time-kill model for two- and three-way antibiotic combinations against Escherichia coli NCTC 12,241.
* Add Themans 2019 meropenem ([doi:10.1007/s40268-019-0268-x](https://doi.org/10.1007/s40268-019-0268-x)) [DDMODEL00000301] -- adults with severe pneumonia.
* Add Kovalenko 2016 dupilumab ([doi:10.1002/psp4.12136](https://doi.org/10.1002/psp4.12136)) [DDMODEL00000273] -- DDMORE-bundle replicate of the existing Kovalenko_2016_dupilumab specificDrugs entry; encodes the bundle's non-standard V2 weight covariate and `Output_simulated_*.lst` final estimates (no `Output_real_*.lst` shipped).
* Add Terranova 2018 paclitaxel ([doi:10.1016/j.jtbi.2018.04.012](https://doi.org/10.1016/j.jtbi.2018.04.012)) [DDMODEL00000274] -- preclinical xenograft mice (Dynamic Energy Budget tumor-growth-inhibition model coupling paclitaxel PK to tumor mass and host body weight / cachexia).
* Add Bajaj 2017 nivolumab ([doi:10.1002/psp4.12143](https://doi.org/10.1002/psp4.12143)) [DDMODEL00000284] -- adults with advanced solid tumors (DDMORE-source replicate of the paper-source Bajaj 2017 nivolumab popPK; time kept in hours to mirror the bundle).
* Add Birgersson 2019 artesunate ([doi:10.12688/wellcomeopenres.14849.2](https://doi.org/10.12688/wellcomeopenres.14849.2)) [DDMODEL00000297] -- pregnant and non-pregnant women with uncomplicated Plasmodium falciparum malaria in Burkina Faso (joint parent-metabolite popPK of artesunate and dihydroartemisinin with 3-compartment transit absorption and pregnancy / ALT / log-parasitaemia covariates).
* Add Dao 2020 sultiame ([doi:10.1002/prp2.558](https://doi.org/10.1002/prp2.558)) [DDMODEL00000298] -- healthy adult volunteers (4-compartment popPK with saturable plasma <-> erythrocyte binding plus cumulative urinary excretion).
* Add Conrado 2014 alzheimer ([doi:10.1007/s10928-014-9375-z](https://doi.org/10.1007/s10928-014-9375-z)) [DDMODEL00000290] -- adults with Alzheimer's disease (CAMD ADAS-Cog database; Richards three-parameter logistic disease-progression model with beta-regression residual).
* Add Schoemaker 2018 levetiracetam ([doi:10.1007/s40262-017-0597-2](https://doi.org/10.1007/s40262-017-0597-2)) [DDMODEL00000239] -- adults and children (4-16 years) with focal seizures (negative-binomial seizure-count PD model with mixture and Markovian dependence on previous-day count; LEV adult+pediatric fit used in the publication to scaffold a brivaracetam pediatric extrapolation).
* Add MPD6 Sutent sunitinib NSCLC PK/PD/tumor-growth model (no linked publication) [DDMODEL00000231] -- semi-mechanistic 15-state ODE model with parent + metabolite 2-compartment oral PK, four indirect-response PD biomarkers, sphere-volume tumor growth, and three resistance / memory chains; MDL-only deposit, validated by F.2 self-consistency only (no `.lst`, no companion paper, no simulated dataset).
* Add Leuppi-Taegtmeyer 2019 colistin ([doi:10.1128/AAC.01957-18](https://doi.org/10.1128/AAC.01957-18)) [DDMODEL00000295] -- critically ill adults receiving colistimethate sodium / colistin during continuous renal replacement therapy.
* Add Voller 2017 phenobarbital ([doi:10.1016/j.ejps.2017.05.026](https://doi.org/10.1016/j.ejps.2017.05.026)) [DDMODEL00000256] -- preterm and term newborns receiving a phenobarbital loading dose followed by oral maintenance.
* Add BAST 2017 PTTE four-event teaching library [DDMODEL00000243] -- 200 simulated patients, no linked publication; four parametric time-to-event hazard models packaged separately as `NA_NA_tte_gompertz` (Event 1, exponential / NEUT + AGE), `NA_NA_tte_gompertz_ev2` (Event 2, Gompertz / first-week AUC), `NA_NA_tte_lognormal` (Competing Event 1, log-normal / AGE), and `NA_NA_tte_loglogistic` (Competing Event 2, log-logistic / no covariate).
* Add Cook 2016 paracetamol ([doi:10.1007/s40262-016-0408-1](https://doi.org/10.1007/s40262-016-0408-1)) [DDMODEL00000271] -- term and preterm newborns receiving IV paracetamol (parent + glucuronide + sulphate plasma compartments with cumulative urinary excretion).
* Add Wang 2013 morphine ([doi:10.1007/s40261-013-0097-6](https://doi.org/10.1007/s40261-013-0097-6)) [DDMODEL00000269] -- neonates / infants / children / adolescents / adults across the entire paediatric age range (bodyweight-dependent allometric exponent on morphine clearance).
* Add Stevens 2012 remoxipride ([doi:10.1007/s10928-012-9262-4](https://doi.org/10.1007/s10928-012-9262-4)) [DDMODEL00000268] -- preclinical Wistar rats (mechanism-based PK/PD pool model for the prolactin response with brain-ECF-driven Emax stimulation and positive feedback on synthesis).
* Add Allegaert 2015 paracetamol ([doi:10.1186/s12871-015-0144-3](https://doi.org/10.1186/s12871-015-0144-3)) [DDMODEL00000267] -- young women across pregnancy / postpartum / non-pregnant-volunteer reproductive states with and without oral contraceptive use (8-compartment IV propacetamol PK with glucuronide and sulphate metabolites and cumulative-urine outputs).
* Update miridesap (CPHPC) [DDMODEL00000262] to the FINAL model from Sahota 2015 ([doi:10.1002/psp4.15](https://doi.org/10.1002/psp4.15)) -- healthy volunteers (CPH113776) and patients with systemic amyloidosis (CPH114527); two-compartment CPHPC PK plus two-compartment SAP turnover with CPHPC-SAP binding/internalization, now with the full final-model covariate panel (AMLIVER on Q4, AMLOAD on V4, SEXF on SAP baseline, CRCL on CL) and parameter values cross-checked against Sahota 2015 Table 2. Ratifies new `AMLOAD` and `AMLIVER` covariate canonicals and adds a `SEX 1=male/2=female` source-alias note to `SEXF`.
* Add Svensson 2018 rifampicin ([doi:10.1002/cpt.778](https://doi.org/10.1002/cpt.778)) [DDMODEL00000244] -- adult pulmonary tuberculosis patients on high-dose rifampicin (PanACEA HIGHRIF1 dose-escalation trial; auto-induction enzyme turnover with Michaelis-Menten clearance and dose-dependent bioavailability).
* Add Wilbaux 2015 prostate ([doi:10.1002/psp4.34](https://doi.org/10.1002/psp4.34)) [DDMODEL00000261] -- adults with metastatic castration-resistant prostate cancer on chemotherapy and/or hormonotherapy (joint K-PD model of CTC count and PSA kinetics with cell-lifespan delay and negative-binomial count likelihood).
* Add Clewe 2018 rifampicin ([doi:10.1093/jac/dkx380](https://doi.org/10.1093/jac/dkx380)) [DDMODEL00000259] -- in vitro M. tuberculosis B1585 time-kill under rifampicin + isoniazid + ethambutol triple-combination GPDI model (scenario 4).
* Add van Rongen 2018 midazolam ([doi:10.1007/s40262-017-0579-4](https://doi.org/10.1007/s40262-017-0579-4)) [DDMODEL00000250] -- 19 obese adolescents and 20 morbidly obese adults (CYP3A probe, oral + IV midazolam, sub-population-specific weight covariates on CL and Vp).
* Add Vet 2016 midazolam ([doi:10.1164/rccm.201510-2114OC](https://doi.org/10.1164/rccm.201510-2114OC)) [DDMODEL00000249] -- critically ill children in the paediatric ICU receiving continuous IV midazolam (per-stratum CL by number of failing organs, with CRP and body weight covariates).
* Add Knibbe 2009 morphine ([doi:10.2165/00003088-200948060-00003](https://doi.org/10.2165/00003088-200948060-00003)) [DDMODEL00000248] -- preterm neonates, term newborns, infants and toddlers <3 years (joint parent-metabolite popPK for morphine + M3G + M6G with body-weight allometric scaling and a postnatal-age-stratified glucuronidation step at PNA = 10 days).
* Add Zurlinden 2016 paracetamol ([doi:10.1007/s13318-015-0253-x](https://doi.org/10.1007/s13318-015-0253-x)) [DDMODEL00000237] -- healthy adults receiving a single 1000 mg oral dose (whole-body PBPK with Michaelis-Menten liver metabolism, cofactor depletion / resynthesis, and renal elimination of APAP and its conjugates AG and AS; first MCSim-source extraction in `inst/modeldb/ddmore/`).
* Add Zierhut 2008 osteoprotegerin ([doi:10.1007/s10928-008-9093-5](https://doi.org/10.1007/s10928-008-9093-5)) [DDMODEL00000233] -- healthy postmenopausal women receiving single IV or SC Fc-osteoprotegerin (two-peripheral-compartment PK with parallel linear and Michaelis-Menten elimination, logistic-style SC bioavailability F = FSC / (1 + FSC), and an indirect-response uNTX biomarker turnover model with route-conditional PK observation residual SD).
* Add Jonsson 2005 disufenton ([doi:10.2165/00003088-200544080-00007](https://doi.org/10.2165/00003088-200544080-00007)) [DDMODEL00000245] -- adults with acute ischaemic or haemorrhagic stroke receiving 72-h IV NXY-059 infusion across CRCL 20-143 mL/min.
* Add Wilkins 2008 rifampicin ([doi:10.1128/AAC.00461-07](https://doi.org/10.1128/AAC.00461-07)) [DDMODEL00000280] -- adult South African pulmonary tuberculosis patients.
* Add NA NA lidocaine ([DDMORE Foundation Model Repository: DDMODEL00000281](https://repository.ddmore.eu/model/DDMODEL00000281)) [DDMODEL00000281] -- population unspecified (no linked publication; lidocaine + MEGX + GX + 2,6-xylidide ADVAN5 parent-and-three-metabolites model).
* Add Laouenan 2015 ribavirin ([doi:10.1002/psp4.8](https://doi.org/10.1002/psp4.8)) [DDMODEL00000285] -- HCV genotype 1 cirrhotic patients on telaprevir- or boceprevir-based triple therapy (hemoglobin turnover PD model with upstream-PK ribavirin regressors).
* Add Clewe 2016 rifampicin ([doi:10.1093/jac/dkv416](https://doi.org/10.1093/jac/dkv416)) [DDMODEL00000240] -- in vitro Mycobacterium tuberculosis H37Rv multistate natural-growth scaffold (fast / slow / non-multiplying bacterial states; rifampicin exposure-response not encoded in the bundled scenario).
* Add Germovsek 2016 gentamicin ([doi:10.1128/AAC.00577-16](https://doi.org/10.1128/AAC.00577-16)) [DDMODEL00000238] -- neonates and infants receiving gentamicin.
* Add Novakovic 2017 cladribine ([doi:10.1208/s12248-016-9977-z](https://doi.org/10.1208/s12248-016-9977-z)) [DDMODEL00000223] -- adults with relapsing-remitting multiple sclerosis (8-item EDSS Item Response Theory disease-progression model with FREM covariates).
* Add Schindler 2016 sunitinib ([doi:10.1002/psp4.12057](https://doi.org/10.1002/psp4.12057)) [DDMODEL00000221] -- adults with imatinib-resistant or imatinib-intolerant advanced gastrointestinal stromal tumor.
* Add Jager 2011 gemtuzumab ([doi:10.1371/journal.pone.0024265](https://doi.org/10.1371/journal.pone.0024265)) [DDMODEL00000229] -- patients with acute myeloid leukemia (mechanism-based PKPD with explicit CD33 binding and leukemic-blast cell killing; DDMORE Monolix re-fit with added peripheral PK compartment).
* Add Khan 2015 ciprofloxacin ([doi:10.1093/jac/dkv233](https://doi.org/10.1093/jac/dkv233)) [DDMODEL00000225] -- in vitro time-kill experiments on E. coli K-12 wild-type and quinolone-resistant single-step mutants.
* Add Ter Heine 2014 tamoxifen ([doi:10.1111/bcp.12388](https://doi.org/10.1111/bcp.12388)) [DDMODEL00000212] -- adults with breast cancer on chronic 20 mg PO QD tamoxifen at steady state (joint parent-metabolite popPK with CYP2D6 and CYP3A4/5 covariates on endoxifen formation).
* Add Mohamed 2016 colistin + meropenem ([doi:10.1093/jac/dkv488](https://doi.org/10.1093/jac/dkv488)) [DDMODEL00000173] -- in vitro time-kill PK/PD against P. aeruginosa wild-type ATCC 27853 and meropenem-resistant ARU552.
* Add Bizzotto 2016 glucose ([doi:10.1152/ajpendo.00045.2016](https://doi.org/10.1152/ajpendo.00045.2016)) [DDMODEL00000227] -- adults across the glucose-tolerance spectrum (mechanistic glucose-tracer kinetics simulator).
* Add Netterberg 2017 docetaxel ([doi:10.1007/s00280-017-3366-x](https://doi.org/10.1007/s00280-017-3366-x)) [DDMODEL00000224] -- adults receiving docetaxel chemotherapy (Friberg-style myelosuppression PD model with Kloft 2006 parameter set).
* Add Zecchin 2016 survival ([doi:10.1111/bcp.12994](https://doi.org/10.1111/bcp.12994)) [DDMODEL00000218] -- adults with advanced epithelial ovarian cancer.
* Add Hansson 2013 sunitinib fatigue Markov + proportional-odds model ([doi:10.1038/psp.2013.62](https://doi.org/10.1038/psp.2013.62)) [DDMODEL00000222] -- adults with imatinib-resistant gastrointestinal stromal tumors.
* Add Svensson 2016 bedaquiline ([doi:10.1002/psp4.12147](https://doi.org/10.1002/psp4.12147)) [DDMODEL00000219] -- adults with multidrug-resistant tuberculosis (parent + N-desmethyl M2 metabolite, time-varying weight and albumin).
* Add Jonsson 2011 ethambutol ([doi:10.1128/AAC.00274-11](https://doi.org/10.1128/AAC.00274-11)) [DDMODEL00000220] -- adult South African pulmonary tuberculosis patients.
* Add Zecchin 2016 tumorovarian ([doi:10.1111/bcp.12994](https://doi.org/10.1111/bcp.12994)) [DDMODEL00000217] -- adults with advanced epithelial ovarian cancer receiving carboplatin monotherapy or carboplatin + gemcitabine combination chemotherapy.
* Add Girard 2012 pimasertib ([www.page-meeting.org/?abstract=2458](https://www.page-meeting.org/?abstract=2458)) [DDMODEL00000215] -- adults with advanced solid tumours and hematological malignancies (joint K-PD / cumulative-logit Markov ocular-AE-grade and Weibull-TTE dropout model).
* Add Lestini 2015 TGF-beta inhibitor ([doi:10.1007/s11095-015-1693-3](https://doi.org/10.1007/s11095-015-1693-3)) [DDMODEL00000192] -- simulated 50-subject oncology cohort (one-compartment first-order absorption PK + indirect-response biomarker turnover; first nlmixr2lib `inst/modeldb/ddmore/` entry).
* Add Li 2006 meropenem ([doi:10.1177/0091270006291035](https://doi.org/10.1177/0091270006291035)) [DDMODEL00000213] -- adult patients receiving meropenem.
* Add Friberg 2002 paclitaxel ([doi:10.1200/JCO.2002.02.140](https://doi.org/10.1200/JCO.2002.02.140)) [DDMODEL00000186] -- adult cancer patients receiving paclitaxel chemotherapy (semi-mechanistic myelosuppression PK/PD with leukocyte output).
* Add Hansson 2013 sunitinib ([doi:10.1038/psp.2013.61](https://doi.org/10.1038/psp.2013.61)) [DDMODEL00000197] -- adults with imatinib-resistant gastrointestinal stromal tumours (biomarker PD model for VEGF, sVEGFR-2, sVEGFR-3, sKIT).
* Add Henin 2009 capecitabine ([doi:10.1038/clpt.2008.220](https://doi.org/10.1038/clpt.2008.220)) [DDMODEL00000214] -- adults with metastatic colorectal or advanced/metastatic breast cancer (Markov-proportional-odds model for hand-and-foot syndrome grades 0-2 with K-PD capecitabine exposure).
* Add Plan 2012 pain ([doi:10.1038/clpt.2011.301](https://doi.org/10.1038/clpt.2011.301)) [DDMODEL00000194] -- adults in placebo arm of three Phase III neuropathic-pain trials (Markov Integer Model for daily 11-point Likert pain scores).
* Add Kloft 2004 sibrotuzumab ([doi:10.1023/B:DRUG.0000006173.72210.1c](https://doi.org/10.1023/B:DRUG.0000006173.72210.1c)) [DDMODEL00000195] -- adults with metastatic FAP-positive solid tumors.
* Change: Standardize Xie 2019 agomelatine parameter, compartment, IOV-eta, and residual-error naming. Compartments rename `DEPOT1` / `DEPOT2` / `LIVER` / `CENTPRNT` / `ALMTPERI` / `METB3OH` / `METB7DM` / `METB7DMPERI` to canonical `depot` / `depot2` / `liver` / `central` / `peripheral1` / `central_3oh` / `central_7dm` / `peripheral1_7dm`; the per-occasion IOV etas (formerly `e.IOV1`-`e.IOV5` and `eta17`-`eta35`) become descriptive `etaiov_<param>_<occ>` names paired with the structural parameter they apply to (`k13`, `alag2`, `k23`, `clint`, `fpop`); residual-error parameters `sdalmt` / `sd3oh` / `sd7dm` become canonical `addSd_lcalmt` / `addSd_lc3oh` / `addSd_lc7dm`; primary/secondary `1` suffixes are dropped where redundant (`fDepot1` -> `fDepot`, `alag1` -> `alag`, `ltvalag1` -> `ltvalag`, `etaltvalag1` -> `etaltvalag`, `F1` -> `fpop`, `etaF1` -> `etafpop`); per-metabolite molecular-weight ratios are renamed descriptively (`mpr1` / `mpr2` -> `mpr_3oh` / `mpr_7dm`). Convention infrastructure (`R/conventions.R`, `R/checkModelConventions.R`) extended to register the `liver` compartment, the `depot[0-9]+` numbered-depot pattern, the `3oh` / `7dm` metabolite suffixes, and the `etaiov_<param>_<occ>` IOV-eta naming pattern.
* Change: Standardize parameter, compartment, covariate-effect, and residual-error naming across the model library. Volumes use `vc` / `vp` / `vp2`; Michaelis-Menten Vmax uses `lvmax` / `vmax`; multi-component CL uses `lcl_ss` / `lcl_time`; shared exponents use `e_<cov>_<param1>_<param2>`; reversed-order covariate effects are flipped to `e_<cov>_<param>`. Parent-metabolite ADC models drop the parent `_adc` suffix and add canonical `<canonical>_<metab>` compartments (`central_mmae`, `central_dxd`, `central_sn38`, `central_tab`, `central_nab`, `central_dm4`, `central_medm4`); metabolite outputs become `Cc_<metab>`. Residual-error parameters now use the parameter-first form `propSd_<X>` / `addSd_<X>` (e.g. `propSd_dxd`, `addSd_tumorSize`) for every non-parent output, while the parent observation `Cc` keeps bare `propSd` / `addSd`. Convention infrastructure (`R/checkModelConventions.R`) extended to enforce the new canonical forms and flag deprecated names.
* Fix Robbie 2012 palivizumab population metadata after revisit against the full-text PDF: pediatric `n` corrected from 1,767 (derived) to 1,684 (paper PK dataset), `sex_female_pct` denominator from 1,660 to 1,684, and added a vignette note clarifying that the 2012 erratum's equation 10 correction (6,900 -> 16,900) is a NONMEM `$PRIOR` variance for adult Vp, not a final-model parameter. All `ini()` parameter values continue to match Table 2 exactly.
* Add de Vries Schultink 2020 zenocutuzumab (MCLA-128) ([doi:10.1007/s40262-020-00858-2](https://doi.org/10.1007/s40262-020-00858-2)) -- adults with advanced solid tumors (first bispecific antibody in nlmixr2lib).
* Add Almquist 2022 anifrolumab ([doi:10.1002/jcph.2055](https://doi.org/10.1002/jcph.2055)) -- adults with moderate-to-severe systemic lupus erythematosus and healthy volunteers (QSS-TMDD with dynamic IFNAR1 receptor pool and time-varying linear clearance).
* Fix Clegg 2024 nirsevimab vignette Figure 4: pass disjoint `id_offset` per cohort and carry `trial` via `rxSolve(keep = )` so the four panels no longer collapse into a single Frankenstein-subject simulation (predictions had been ~3x too high).
* Add Lu 2022 patritumab deruxtecan ([doi:10.1002/jcph.2137](https://doi.org/10.1002/jcph.2137)) -- adults with HER3-expressing advanced or metastatic solid tumors (NSCLC, breast cancer, colorectal cancer).
* Add Cheng 2026 immunoglobulin ([doi:10.1002/bcp.70420](https://doi.org/10.1002/bcp.70420)) -- children with primary immunodeficiency or secondary antibody deficiency receiving intravenous immunoglobulin replacement therapy.
* Add Frey 2010 tocilizumab ([doi:10.1177/0091270009350623](https://doi.org/10.1177/0091270009350623)) -- adults with moderate-to-severe rheumatoid arthritis.
* Add Hayashi 2007 omalizumab ([doi:10.1111/j.1365-2125.2006.02803.x](https://doi.org/10.1111/j.1365-2125.2006.02803.x)) -- Japanese adults with atopic asthma or seasonal allergic rhinitis.
* Add Hu 2014 bapineuzumab ([doi:10.1002/jcph.393](https://doi.org/10.1002/jcph.393)) -- adults with mild-to-moderate Alzheimer's disease.
* Add Huang 2017 VRC01 ([doi:10.1080/19420862.2017.1311435](https://doi.org/10.1080/19420862.2017.1311435)) -- HIV-uninfected healthy adults receiving IV or SC HIV-1 broadly neutralizing monoclonal antibody.
* Add Ide 2020 elotuzumab ([doi:10.1002/jcph.1698](https://doi.org/10.1002/jcph.1698)) -- Japanese and non-Japanese adults with multiple myeloma.
* Add Marquez-Megias 2023 adalimumab ([doi:10.3390/biomedicines11102822](https://doi.org/10.3390/biomedicines11102822)) -- adults with inflammatory bowel disease.
* Add Nestorov 2014 factor VIII Fc fusion protein ([doi:10.1002/cpdd.167](https://doi.org/10.1002/cpdd.167)) -- previously treated patients with severe hemophilia A.
* Add Perez-Ruixo 2025 posdinemab ([doi:10.1002/cpt.70173](https://doi.org/10.1002/cpt.70173)) -- healthy adults and adults with Alzheimer's disease.
* Add Petrov 2024 romiplostim ([doi:10.1002/cpdd.1494](https://doi.org/10.1002/cpdd.1494)) -- adults with chronic immune thrombocytopenia (ITP).
* Add Pouzin 2022 tusamitamab ravtansine ([doi:10.1007/s10928-021-09799-0](https://doi.org/10.1007/s10928-021-09799-0)) -- adults with advanced solid tumors expressing CEACAM5 (multi-analyte ADC + NAB + DM4 + MeDM4 model).
* Add Takeuchi 2023 ozoralizumab ([doi:10.1002/jcph.2380](https://doi.org/10.1002/jcph.2380)) -- Japanese adults with rheumatoid arthritis.
* Add Toukam 2025 BIIB107 ([doi:10.1002/jcph.70109](https://doi.org/10.1002/jcph.70109)) -- healthy adult volunteers (first-in-human study supporting MS dose optimization).
* Add Yao 2018 guselkumab ([doi:10.1002/jcph.1063](https://doi.org/10.1002/jcph.1063)) -- adults with moderate-to-severe plaque psoriasis.
* Add Frey 2013 tocilizumab DAS28 exposure-response ([doi:10.1177/0091270012437585](https://doi.org/10.1177/0091270012437585)) -- adults with moderate-to-severe rheumatoid arthritis (OPTION + TOWARD phase III pool).
* Add Hwang 2023 monalizumab ([doi:10.1002/jcph.2220](https://doi.org/10.1002/jcph.2220)) -- adults with advanced solid tumors or squamous cell carcinoma of the head and neck.
* Add Suri 2018 brentuximab vedotin ([doi:10.1002/cpt.1037](https://doi.org/10.1002/cpt.1037)) -- adults with CD30-positive Hodgkin lymphoma, systemic anaplastic large-cell lymphoma, mycosis fungoides, or primary cutaneous ALCL pooled across six clinical studies (including the phase III ALCANZA trial in CTCL).
* Add Zhong 2026 abatacept ([doi:10.1002/jcph.70156](https://doi.org/10.1002/jcph.70156)) -- pooled adults with rheumatoid arthritis, patients aged 2-17 years with polyarticular juvenile idiopathic arthritis, and patients aged 6+ years with hematologic malignancies receiving HLA-matched unrelated-donor HSCT.
* Add Mukai 2019 mogamulizumab ([doi:10.1002/jcph.1564](https://doi.org/10.1002/jcph.1564)) -- adults with cutaneous T-cell lymphoma or adult T-cell lymphoma.
* Add Lin 2024 casirivimab ([doi:10.1007/s11095-024-03764-5](https://doi.org/10.1007/s11095-024-03764-5)) -- pediatric and adult subjects, non-infected or with SARS-CoV-2 infection or household contacts.
* Add Faelens 2021 infliximab ([doi:10.3390/pharmaceutics13101623](https://doi.org/10.3390/pharmaceutics13101623)) -- adults with moderate-to-severe ulcerative colitis.
* Add Sathe 2024 sacituzumab govitecan ([doi:10.1007/s40262-024-01366-3](https://doi.org/10.1007/s40262-024-01366-3)) -- adults with metastatic triple-negative breast cancer and other solid tumors.
* Add Takahashi 2023 abatacept ([doi:10.1182/blood.2023020035](https://doi.org/10.1182/blood.2023020035)) -- pooled adult/pediatric rheumatoid-arthritis / pJIA patients and adult/pediatric hematopoietic-cell-transplant recipients in the ABA2 trial.
* Add Hong 2025 datopotamab deruxtecan ([doi:10.1002/psp4.70118](https://doi.org/10.1002/psp4.70118)) -- adults with advanced or metastatic NSCLC and breast cancer.
* Add Fau 2020 isatuximab ([doi:10.1002/psp4.12561](https://doi.org/10.1002/psp4.12561)) -- adults with relapsed/refractory multiple myeloma.
* Add Lu 2019 polatuzumab vedotin ([doi:10.1002/psp4.12482](https://doi.org/10.1002/psp4.12482)) -- adults with relapsed/refractory or previously untreated B-cell non-Hodgkin lymphoma.
* Add Okada 2025 rocatinlimab ([doi:10.1002/psp4.70069](https://doi.org/10.1002/psp4.70069)) -- adults with moderate-to-severe atopic dermatitis (pooled with plaque psoriasis, ulcerative colitis, and healthy volunteers).
* Add Zhou 2025 brentuximab vedotin ([doi:10.1002/cpt.3629](https://doi.org/10.1002/cpt.3629)) -- pediatric patients (5-18 years) with newly diagnosed Hodgkin lymphoma or relapsed/refractory systemic anaplastic large-cell lymphoma.
* Add Lin 2024 pozelimab ([doi:10.1007/s10928-024-09941-8](https://doi.org/10.1007/s10928-024-09941-8)) -- healthy adult volunteers, adults with paroxysmal nocturnal hemoglobinuria, and pediatric and adult patients with CHAPLE disease.
* Add Yin 2021 trastuzumab deruxtecan ([doi:10.1002/cpt.2096](https://doi.org/10.1002/cpt.2096)) -- adults with HER2-positive breast cancer or other HER2-expressing solid tumors.
* Add Hwang 2022 tremelimumab ([doi:10.1111/bcp.15622](https://doi.org/10.1111/bcp.15622)) -- adults with advanced solid tumours.
* Add Papathanasiou 2025 belantamab mafodotin ([doi:10.1007/s40262-025-01508-1](https://doi.org/10.1007/s40262-025-01508-1)) -- adults with relapsed/refractory multiple myeloma (ADC moiety only).
* Add Kuchimanchi 2024 dostarlimab ([doi:10.1111/bcp.16325](https://doi.org/10.1111/bcp.16325)) -- adults with primary advanced or recurrent endometrial cancer (RUBY Part 1) and advanced solid tumours (GARNET).
* Add Diao 2016 daclizumab CD25 occupancy ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) -- adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Diao 2016 daclizumab CD56 bright NK expansion ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) -- adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Diao 2016 daclizumab Treg reduction ([doi:10.1111/bcp.13051](https://doi.org/10.1111/bcp.13051)) -- adults with relapsing-remitting multiple sclerosis (PD model with Othman 2014 PK backbone).
* Add Fiedler-Kelly 2020 fremanezumab exposure-response ([doi:10.1111/head.13855](https://doi.org/10.1111/head.13855)) -- adults with episodic migraine and adults with chronic migraine (two PD-only models: `FiedlerKelly_2020_fremanezumab_em` and `FiedlerKelly_2020_fremanezumab_cm`).
* Add Koopman 2023 factor IX-Fc ([doi:10.1111/bcp.15881](https://doi.org/10.1111/bcp.15881)) -- children, adolescents, and adults with haemophilia B (real-world data including patients aged < 12 years).
* Add Le Tilly 2021 trastuzumab ([doi:10.1002/cpt.2188](https://doi.org/10.1002/cpt.2188)) -- adults with HER2+ breast cancer leptomeningeal carcinomatosis receiving intrathecal and intravenous trastuzumab.
* Add Wang 2024 sugemalimab ([doi:10.1111/bcp.16276](https://doi.org/10.1111/bcp.16276)) -- adults with advanced solid tumours or lymphomas across nine Phase I-III sugemalimab trials.
* Add Yang 2024 axatilimab ([doi:10.1002/cpt.3503](https://doi.org/10.1002/cpt.3503)) -- pooled healthy adults, adults with advanced solid tumors, and adults / children with chronic graft-versus-host disease.
* Add Hood 2021 MEDI7836 ([doi:10.3390/pharmaceutics13050613](https://doi.org/10.3390/pharmaceutics13050613)) -- healthy adult males in a first-in-human single-ascending-dose trial.
* Amend Castro-Suarez 2020 nimotuzumab: V1 decreased 53% for 50 mg cohort per Figure 4 visual inspection (direction not stated in paper text; corresponding author contacted).
* Add Castro-Suarez 2020 nimotuzumab ([doi:10.3390/pharmaceutics12121147](https://doi.org/10.3390/pharmaceutics12121147)) -- adults with autosomal dominant polycystic kidney disease.
* Add Yang 2021 cemiplimab ([doi:10.1007/s10928-021-09739-y](https://doi.org/10.1007/s10928-021-09739-y)) -- adults with advanced solid tumors including cutaneous squamous cell carcinoma.
* Add Papachristos 2020 bevacizumab ([doi:10.3390/ijms21113753](https://doi.org/10.3390/ijms21113753)) -- adults with metastatic colorectal cancer (three co-equal final models: descriptive PK, binding QSS TMDD, and PK/PD Imax inhibition of free VEGF-A).
* Add Ngo 2020 HL2351 ([doi:10.1002/psp4.12552](https://doi.org/10.1002/psp4.12552)) -- healthy adult Korean men.
* Add Wojciechowski 2022 domagrozumab ([doi:10.1111/cts.13418](https://doi.org/10.1111/cts.13418)) -- healthy adult volunteers and pediatric patients with Duchenne muscular dystrophy.
* Add Yu 2022 ofatumumab ([doi:10.1007/s40263-021-00895-w](https://doi.org/10.1007/s40263-021-00895-w)) -- adults with relapsing multiple sclerosis.
* Add Melhem 2022 dostarlimab ([doi:10.1111/bcp.15339](https://doi.org/10.1111/bcp.15339)) -- adults with advanced solid tumours.
* Add Brillac 2025 isatuximab ([doi:10.1007/s00280-025-04832-2](https://doi.org/10.1007/s00280-025-04832-2)) -- pediatric and adult patients with relapsed/refractory acute leukemias.
* Add Wu 2024 inotuzumab ozogamicin ([doi:10.1007/s40262-024-01386-z](https://doi.org/10.1007/s40262-024-01386-z)) -- pediatric and adult patients with relapsed/refractory B-cell precursor acute lymphoblastic leukemia and adults with B-cell non-Hodgkin's lymphoma.
* Fix vignettes: derive concentration units in `labs()` from model `$units` metadata; replace inline trapezoidal NCA with PKNCA; add PKNCA sections to nalmefene and clesrovimab vignettes; add PKNCA treatment grouping to benralizumab and durvalumab vignettes.
* Add Chen 2022 guselkumab ([doi:10.1111/cts.13197](https://doi.org/10.1111/cts.13197)) -- adults with active psoriatic arthritis (DISCOVER-1 and DISCOVER-2 phase 3 trials).
* Add `concentration` field to `units` metadata for `Wang_2017_benralizumab` and `Ogasawara_2020_durvalumab` models.
* Add Chen 2020 luspatercept ([doi:10.1002/psp4.12515](https://doi.org/10.1002/psp4.12515)) -- adults with anemia due to lower-risk myelodysplastic syndromes.
* Add Martinez 2019 alirocumab ([doi:10.1007/s40262-018-0669-y](https://doi.org/10.1007/s40262-018-0669-y)) -- healthy volunteers and adults with hypercholesterolemia.
* Add Li 2017 brentuximab vedotin ([doi:10.1002/jcph.920](https://doi.org/10.1002/jcph.920)) -- adults with relapsed/refractory CD30-expressing hematologic malignancies (Hodgkin lymphoma and systemic anaplastic large cell lymphoma).
* Add Quartino 2019 trastuzumab ([doi:10.1007/s00280-018-3728-z](https://doi.org/10.1007/s00280-018-3728-z)) -- adults with metastatic breast cancer, early breast cancer, advanced gastric cancer, or other solid tumors.
* Add Suleiman 2019 risankizumab ([doi:10.1007/s40262-019-00759-z](https://doi.org/10.1007/s40262-019-00759-z)) -- healthy volunteers and adults with moderate-to-severe plaque psoriasis.
* Add Berends 2019 infliximab ([doi:10.1007/s10928-019-09652-5](https://doi.org/10.1007/s10928-019-09652-5)) -- adults with moderate-to-severe ulcerative colitis.
* Add Nikanjam 2019 siltuximab ([doi:10.1007/s00280-019-03939-7](https://doi.org/10.1007/s00280-019-03939-7)) -- adults pooled across healthy volunteers and oncology cohorts including Castleman's disease and smoldering multiple myeloma.
* Add Sanghavi 2020 ipilimumab ([doi:10.1002/psp4.12477](https://doi.org/10.1002/psp4.12477)) -- adults with advanced solid tumors receiving ipilimumab alone or in combination with nivolumab.
* Add Zhang 2019 nivolumab ([doi:10.1002/psp4.12476](https://doi.org/10.1002/psp4.12476)) -- adults with advanced solid tumors, alone or in combination with ipilimumab or chemotherapy.
* Add Wang 2020 ontamalimab ([doi:10.1002/jcph.1590](https://doi.org/10.1002/jcph.1590)) -- adults with moderate-to-severe ulcerative colitis or Crohn's disease.
* Add Hanzel 2021 infliximab CT-P13 ([doi:10.1111/apt.16609](https://doi.org/10.1111/apt.16609)) -- adults with Crohn's disease or ulcerative colitis.
* Add Zhou 2021 belimumab ([doi:10.1007/s40268-021-00363-2](https://doi.org/10.1007/s40268-021-00363-2)) -- adult and pediatric patients with systemic lupus erythematosus (Chinese and non-Chinese).
* Add Aguiar 2021 ustekinumab ([doi:10.3390/pharmaceutics13101587](https://doi.org/10.3390/pharmaceutics13101587)) -- adults with active Crohn's disease.
* Add Li 2019 abatacept ([doi:10.1002/jcph.1308](https://doi.org/10.1002/jcph.1308)) -- adults with rheumatoid arthritis.
* Add Gandhi 2021 abatacept ([doi:10.1002/jcph.1781](https://doi.org/10.1002/jcph.1781)) -- pooled adults with rheumatoid arthritis and patients aged 2-17 years with polyarticular juvenile idiopathic arthritis.
* Add Mulyukov 2018 ranibizumab ([doi:10.1002/psp4.12322](https://doi.org/10.1002/psp4.12322)) -- anti-VEGF-naive adults with neovascular age-related macular degeneration.
* Add Bajaj 2017 nivolumab ([doi:10.1002/psp4.12143](https://doi.org/10.1002/psp4.12143)) -- patients with advanced solid tumors (melanoma, NSCLC, RCC, other).
* Add Kielbasa 2020 galcanezumab ([doi:10.1002/jcph.1511](https://doi.org/10.1002/jcph.1511)) -- healthy adults and adults with episodic or chronic migraine.
* Add Othman 2014 daclizumab ([doi:10.1007/s40262-014-0159-9](https://doi.org/10.1007/s40262-014-0159-9)) -- healthy adult volunteers.
* Add Timmermann 2019 brodalumab ([doi:10.1111/bcpt.13202](https://doi.org/10.1111/bcpt.13202)) -- adults with moderate-to-severe plaque psoriasis.
* Add Kuchimanchi 2018 evolocumab ([doi:10.1007/s10928-018-9592-y](https://doi.org/10.1007/s10928-018-9592-y)) -- healthy adults and patients with hypercholesterolemia (including heterozygous familial hypercholesterolemia and statin-intolerant cohorts).
* Add Mo 2018 olaratumab ([doi:10.1007/s40262-017-0562-0](https://doi.org/10.1007/s40262-017-0562-0)) -- patients with advanced or metastatic cancer (soft tissue sarcoma, NSCLC, GIST, glioblastoma multiforme).
* Add Kovalenko 2016 dupilumab ([doi:10.1002/psp4.12136](https://doi.org/10.1002/psp4.12136)) -- healthy volunteers and adults with moderate-to-severe atopic dermatitis.
* Add Narwal 2013 sifalimumab ([doi:10.1007/s40262-013-0085-2](https://doi.org/10.1007/s40262-013-0085-2)) -- adults with moderately active systemic lupus erythematosus.
* Add Valenzuela 2025 nipocalimab ([doi:10.1002/psp4.70109](https://doi.org/10.1002/psp4.70109)) -- healthy adults and adults with generalized myasthenia gravis.
* Add Hua 2015 anrukinzumab ([doi:10.1111/bcp.12589](https://doi.org/10.1111/bcp.12589)) -- healthy volunteers, asthma, and ulcerative colitis patients.
* Add Robbie 2012 palivizumab ([doi:10.1128/AAC.06446-11](https://doi.org/10.1128/AAC.06446-11)) -- adults and pediatric patients at high risk of RSV lower respiratory tract disease.
* Add Bender 2014 trastuzumab emtansine (T-DM1, ADC) ([doi:10.1208/s12248-014-9618-3](https://doi.org/10.1208/s12248-014-9618-3)) -- preclinical: rats + cynomolgus monkeys; two model variants (reduced + mechanistic DAR chain).
* Add Rosario 2015 vedolizumab ([doi:10.1111/apt.13243](https://doi.org/10.1111/apt.13243)) -- adults with moderately-to-severely active ulcerative colitis or Crohn's disease.
* Add Long 2017 necitumumab ([doi:10.1007/s40262-016-0452-x](https://doi.org/10.1007/s40262-016-0452-x)) -- cancer patients with advanced solid tumors.
* Add Wade 2015 certolizumab ([doi:10.1002/jcph.491](https://doi.org/10.1002/jcph.491)) -- adults with Crohn's disease.
* Add Lon 2013 abatacept ([doi:10.1007/s10928-013-9341-1](https://doi.org/10.1007/s10928-013-9341-1)) -- male Lewis rats with collagen-induced arthritis (preclinical).
* Add Lu 2014 trastuzumab emtansine ([doi:10.1007/s00280-014-2500-2](https://doi.org/10.1007/s00280-014-2500-2)) -- adults with HER2-positive locally advanced or metastatic breast cancer.
* Add Gupta 2016 amatuximab ([doi:10.1007/s00280-016-2984-z](https://doi.org/10.1007/s00280-016-2984-z)) -- adults with advanced mesothelin-expressing cancers including malignant pleural mesothelioma.
* Add Zheng 2016 sifalimumab ([doi:10.1111/bcp.12864](https://doi.org/10.1111/bcp.12864)) -- adults with systemic lupus erythematosus.
* Add Mould 2007 alemtuzumab ([doi:10.1111/j.1365-2125.2007.02914.x](https://doi.org/10.1111/j.1365-2125.2007.02914.x)) -- adults with B-cell chronic lymphocytic leukaemia.
* Add Farrell 2012 farletuzumab ([doi:10.1007/s00280-012-1959-y](https://doi.org/10.1007/s00280-012-1959-y)) -- women with advanced epithelial ovarian cancer.
* Add Xu 2011 sirukumab ([doi:10.1111/j.1365-2125.2011.03964.x](https://doi.org/10.1111/j.1365-2125.2011.03964.x)) -- healthy adult volunteers in a first-in-human study.
* Add Yamada 2025 zolbetuximab ([doi:10.1111/cts.70280](https://doi.org/10.1111/cts.70280)) -- patients with locally advanced unresectable or metastatic gastric/gastroesophageal junction adenocarcinoma.
* Add Jackson 2022 ixekizumab ([doi:10.1111/bcp.15034](https://doi.org/10.1111/bcp.15034)) -- paediatric patients with moderate-to-severe plaque psoriasis.
* Add Moein 2022 etrolizumab ([doi:10.1002/psp4.12846](https://doi.org/10.1002/psp4.12846)) -- patients with moderately-to-severely active ulcerative colitis.
* Add Chua 2025 mirikizumab ([doi:10.1111/cts.70320](https://doi.org/10.1111/cts.70320)) -- patients with moderately-to-severely active Crohn's disease.
* Add Ma 2020 sarilumab ANC ([doi:10.1007/s40262-020-00899-7](https://doi.org/10.1007/s40262-020-00899-7)) -- adults with rheumatoid arthritis.
* Add Ma 2020 sarilumab DAS28-CRP ([doi:10.1007/s40262-020-00899-7](https://doi.org/10.1007/s40262-020-00899-7)) -- adults with rheumatoid arthritis.
* Add Tiraboschi 2025 amlitelimab ([doi:10.1002/psp4.70121](https://doi.org/10.1002/psp4.70121)) -- adults with moderate-to-severe atopic dermatitis.
* Add Budha 2023 tislelizumab ([doi:10.1002/psp4.12880](https://doi.org/10.1002/psp4.12880)) -- patients with advanced tumors.
* Add Masters 2022 avelumab ([doi:10.1002/psp4.12771](https://doi.org/10.1002/psp4.12771)) -- patients with advanced solid tumors.
* Add Xu 2019 sarilumab ([doi:10.1007/s40262-019-00765-1](https://doi.org/10.1007/s40262-019-00765-1)) -- adults with rheumatoid arthritis.
* `checkModelConventions()` -- new function that reports deviations from the package's parameter-naming, covariate, compartment, and metadata conventions for a single model or the entire `modeldb`. Called automatically during `buildModelDb()` so convention drift surfaces at package-build time (existing grandfathered deviations continue to build). Canonical standards (PK parameter prefixes, `eta`-prefixed IIV, `propSd`/`addSd` residual error, canonical covariate column register with aliases, canonical compartment vocabulary) are codified in an internal `.nlmixr2libConventions` list that mirrors the `extract-literature-model` skill references. Addresses issue #39.
* Added canonical TMDD archetype models under `inst/modeldb/pharmacokinetics/`: `PK_1cmt_tmdd_full` (Mager & Jusko 2001), `PK_1cmt_tmdd_qss`, `PK_1cmt_tmdd_mm`, `PK_2cmt_tmdd_qss`, and `PK_2cmt_tmdd_mm` (Gibiansky et al. 2008 QSS and MM approximations). Built to the `extract-literature-model` skill conventions with `reference` / `units` / `population` metadata, per-parameter source-trace comments, and CL/V parameterization with `kel` derived inside `model()` so later transforms can re-parameterize. Replaces the 38 draft models from PR #60 (#60).
* Registered `target`, `complex`, and `total_target` as canonical TMDD compartment names in the `extract-literature-model` skill naming conventions (per @iamstein's proposal on PR #60).
* Added `vignettes/tmdd_archetypes.Rmd` comparing the five TMDD archetypes with typical-value trajectories of drug, free target, and complex, and a regime-convergence check showing QSS/MM collapsing onto the full Mager & Jusko 2001 model in the fast-binding regime of Gibiansky 2008.
* Add Thakre 2022 risankizumab ([doi:10.1007/s40744-022-00495-0](https://doi.org/10.1007/s40744-022-00495-0)) -- patients with active psoriatic arthritis.
* `addEta()`, `addResErr()`, `addDepot()`, `removeDepot()`, `addTransit()`, and `removeTransit()` now accept `model` as a deprecated alias for `ui` (issue #84). Passing `model = ...` emits a deprecation warning; passing both `ui` and `model` is an error.
* `addDepot()` and `addTransit()` now work correctly when `d/dt(central)` or `d/dt(depot)` appears at the beginning or end of the model block, or when transit-compartment ODEs and residual-error (`~`) specs are interleaved with assignment lines. The newly introduced helper and ODE lines are inserted immediately adjacent to the modified ODE so that the relative order of every pre-existing model line is preserved (#77, #78).
* Markov modeling creation functions including `createMarkovModel()` were added
* Add Kotani 2022 astegolimab ([doi:10.1002/jcph.2021](https://doi.org/10.1002/jcph.2021)) -- adults with severe asthma.
* Add Fasanmade 2009 infliximab ([doi:10.1007/s00228-009-0718-4](https://doi.org/10.1007/s00228-009-0718-4)) -- adults with moderately-to-severely active ulcerative colitis.
* Add Fiedler-Kelly 2019 fremanezumab ([doi:10.1111/bcp.14096](https://doi.org/10.1111/bcp.14096)) -- adults with chronic or episodic migraine.
* Add Hu 2026 clesrovimab ([doi:10.1002/cpt.70199](https://doi.org/10.1002/cpt.70199)) -- preterm and full-term infants.
* Add Clegg 2024 nirsevimab ([doi:10.1002/jcph.2401](https://doi.org/10.1002/jcph.2401)) -- preterm and term infants.
* Verified all published-literature specific-drug and mAb-consensus models against their source papers and fixed several parameter-encoding bugs that had been latent in the package since their original addition:
  - **CarlssonPetri 2021 liraglutide**: fixed categorical covariate encoding that was zeroing individual clearance for subjects not in the indexed group. `(1 - SEXF)^e_sex_cl` -> `e_sex_cl^(1 - SEXF)` (previously evaluated `0^1.12 = 0` for females); `CHILD^e_age_child_cl * ADOLESCENT^e_age_adolescent_cl` -> `e_age_child_cl^CHILD * e_age_adolescent_cl^ADOLESCENT` (previously evaluated `0^1.11 * 0^1.06 = 0` for adults). IIV rewritten as `omega^2 = log(1 + CV^2)` per Table 3's explicit `%CV = sqrt(exp(omega^2) - 1) * 100` footnote.
  - **Zhu 2017 lebrikizumab**: fixed IIV variance-covariance block that was storing `sqrt(variance)` (SDs) instead of variances/covariances from Table 3.
  - **Soehoel 2022 tralokinumab**: fixed IIV block (SDs -> variance-covariance matrix from Table 2 footnote `IIV = sqrt(exp(omega^2) - 1)`, correlation 0.61 applied on variances not SDs); corrected body-weight exponent on V2/V3 from `0.793` to Table 2's `0.783`.
  - **Kovalenko 2020 dupilumab**: squared the five IIV values so they store variances. Paper Methods explicitly defines `omega` as "the standard deviation [SD] of between-subject variability", and nlmixr2 `etaX ~ value` stores the variance (omega^2), so SDs needed squaring.
  - **Davda 2014 mAb consensus (PK_2cmt_mAb_Davda_2014)**: all parameters verified; no changes required.
* Filled in previously-TODO `population` metadata blocks for the five verified models above with demographics from each paper's Table 1.
* Added validation vignettes for the five verified models above (`CarlssonPetri_2021_liraglutide.Rmd`, `Zhu_2017_lebrikizumab.Rmd`, `Soehoel_2022_tralokinumab.Rmd`, `Kovalenko_2020_dupilumab.Rmd`, `PK_2cmt_mAb_Davda_2014.Rmd`). Each vignette follows the `extract-literature-model` skill conventions with a population description, source trace, virtual cohort, simulation, figure replication, and PKNCA-based NCA validation.
* Retrofit Cirincione 2017 exenatide model to the `extract-literature-model` skill conventions and fix parameter encoding bugs: `ka_max` corrected from `0.0813` to paper value `12.8` /hr, `Km` rescaled to ng/mL so units are consistent with `Cc = central / vc`, and IIV variances rewritten as `log(1 + CV^2)` rather than the `log(1 + CV)` shortcut. Replaces the character-valued `DVID` covariate with `STUDY1` / `STUDY5` binary indicators and adds a companion validation vignette with PKNCA checks against the paper's Figure 5 typical values.
* Fix endogenous-model bugs surfaced during a parameter audit:
  - **Charbonneau 2021 phenylalanine**: fixed `vd_phe` typo on the `f_gut_plasma` line (the undefined symbol made `rxSolve()` fail) -- should read `vd`; removed a stray `vd *` factor from the `daily_phe_intake` augmentation output so the reported value is in mg/day rather than (L/kg)*mg/day. All Table 4 parameter values verified against the authors' Zenodo companion notebook.
  - **Kim 2006 IgG**: corrected `V1` units label from `(mg/kg)` to `(mL/kg)` (paper Table 1); closed a missing `)` in the `ljmax` label; removed a redundant `igg_0 <- css` line. Parameter values verified against Table 1.
* Extended the `extract-literature-model` skill to cover endogenous, mechanistic, and turnover models: added naming conventions for Vmax/Km/baseline/rate parameters, a new `references/endogenous-validation.md` reference covering steady-state / perturbation-recovery / mass-balance / dimensional-analysis patterns, and an explicit dimensional-analysis item in `references/verification-checklist.md` (driven by the Kim units-label and Charbonneau `daily_phe_intake` bugs, which were only catchable by dimensional analysis).
* Retrofit the remaining specific-drug models and the Davda 2014 mAb 2-compartment model to the `extract-literature-model` skill conventions (no parameter-value changes). Structured `covariateData` with `description` / `units` / `type` / `reference_category` / `notes` / `source_name`, canonical `units` list (time/dosing/concentration), and `population` metadata blocks added (with TODO placeholders where not yet sourced). IIV etas renamed to `eta` + transformed-parameter name (e.g., `etacl` -> `etalcl`, `iiv_lka` -> `etalka`, `bsv_fpla_*` -> `etalfpla_*`). Race columns renamed to the `RACE_` prefix (`BLACK` -> `RACE_BLACK`, `ASIAN` -> `RACE_ASIAN`, `MULTIRACIAL` -> `RACE_MULTI`, `BLACK_OTH` -> `RACE_BLACK_OTH`, `ASIAN_AMIND_MULTI` -> `RACE_ASIAN_AMIND_MULTI`); `ADA` -> `ADA_POS`; `SEXM` -> `SEXF` (value inverted to keep the original effect magnitude). Touched vignettes updated to use the canonical column names. Models touched: `CarlssonPetri_2021_liraglutide`, `Clegg_2024_nirsevimab`, `Grimm_2023_gantenerumab`, `Grimm_2023_trontinemab`, `Hu_2026_clesrovimab`, `Kovalenko_2020_dupilumab`, `Kyhl_2016_nalmefene`, `PK_2cmt_mAb_Davda_2014`, `Soehoel_2022_tralokinumab`, `Xie_2019_agomelatine`, `Zhu_2017_lebrikizumab`.

# Version 0.3.2

* Add Kim 2006 model for IgG metabolism
* Add Xie 2019 agomelatine PK model
* Drop `qs` since it will be archived
* Update tumor growth inhibition models
* `addResErr()` now works with multiple-endpoint models
* Additional testing

# Version 0.3.1

* Bug fix for replacement of multiplicative expressions in `nlmixr2lib`
* phenylalanine_charbonneau_2021 had its net protein breakdown parameter corrected
* Kyhl_2016_nalmefene model was added

# Version 0.3.0

* Added ability to choose style type when modifying models.  Currently
  supported styles are: "camel" for `variablesLikeThis`, "snake" for
  `variables_like_this`, "dot" for `variables.like.this` and "blank"
  for `variableslikethis`.  This can be selected with
  `setCombineType()`.

* With the new combination style, you can change how `eta` variables
  are constructed with the `option(nlmixr2lib.etaCombineType="camel")`
  or whatever you wish it to the variable style to be.

* Added new model building framework for building models

  - **PK model building functions**

     - `addTransit()`/`removeTransit()` which were present before, but now modified and
       made a bit more robust, more closely matching literature method
       of transit compartments.

     - `addDepot()`/`removeDepot()` which were present before, but
       modified to be a bit more robust.

     - `addWeibullAbs()` which adds a Weibull absorption to a PK model

     - `convertMM()` converts linear elimination to Michaelis-Menten elimination

     - `transPK()` converts the `cl` style parameter transformations
       to various other PK transformations like `k`, `aob`, `alpha`,
       `k12`

  - **PD model building functions**

   - `addIndirectLin()` -- this adds an indirect effect model to a PK
     model that has a concentration `Cc` in the model.  This purposely
     uses a simple linear effect of `Cc*Ek` or `Cc*Ik` so it will be
     easy to parse and turn into other functional forms (like `Emax`
     or `Hill`).  If the PK model is not present it will use `Cc` as a
     covariate in a purely PD models.

   - `addIndirect()` -- this builds on `addIndirectLin()` and adds
     `Emax` or `Hill` models to a PK model. You can also set `imax=1`
     or `emax=1` to drop these parameters from being estimated in the
     model.  Additionally `hill=TRUE` will add a Hill coefficient to
     the sigmoid model.

   - `addEffectCmtLin()` -- this adds an effect compartment based on
     the `Cc` in the model.  The linear effect can be modified into
     other function forms.

   - `addDirectLin()` -- this adds a direct effect model based on the
     `Cc` in the model.

   - **Changing functional forms of Effect models**

     - `convertEmax()` changes linear effect models to Emax models

     - `convertEmaxHill()` changes linear effect models to Hill models

     - `convertQuad()` changes linear effect models to quadratic models

     - `convertLogLin()` changes linear effect models to log-linear models

   - **Changing functional forms of Baselines in non-indirect response models**

     - `addBaselineConst()` changes the zero baseline to a estimated
       constant

     - `addBaselineLin()` changes the zero baseline to a estimated
       constant and a linear constant with respect to `time`.

     - `addBaselineExp()` changes the zero baseline to a exponential
       decay with respect to time

     - `addBaseline1exp()` -- the baseline effect is changed from zero
       to to an exponential approaching to a constant (with respect to
       time).

   - **Changing model properties** (all use `addCmtProp()`)

      - `addBioavailability()` adds bioavailability property to a
        compartment

      - `addRate()` adds a modeled rate to a compartment

      - `addDur()` adds modeled duration to a compartment

      - `addIni()` adds an initial value to a compartment

      - `addLag()` adds a lag time to the a compartment

* Add Carlsson Petri (2021) liraglutide PK model
* Add Cirincione (2017) exenatide immediate-release PK model
* Add a variety of indirect response models
* Add a variety of tumor growth inhibition models and move all oncology models
  into a new model database directory
* Add a variety of double-absorption PK models
* `cp` and related `cpddSd` and `cppropSd` were renamed to `Cc`, `CcAddSd` and
  `CcPropSd` (fix #70).
* Multiple-endpoint models will have the `DV` column in the modeldb separated by
  commas.

# Version 0.2.0

* Work with the new `rxode2` version 2.0.12 `model()` and `ini()` assignment
  methods.
* Therapeutic-area specific models have begun being added.
* Models can now give the user some additional information load via the
  `message` meta-data.
* Models can now be in different directories.  The change is for ease of
  maintaining the library, it is not a change that affects users.
* A regression where `addEta()` did not change the parameter, related to a
  change in `rxode2`, was fixed.
* `addEta()` detects where to add etas more robustly when covariates are on the
  parameter.

## Models added

* Add Davda (2014) mAb consensus model
* Add Liu (2017) time-dependent clearance model based on nivolumab
* Add Kovalenko (2020) dupilumab PK model
* Add Soehoel (2022) tralokinumab PK model
* Add Zhu (2017) lebrikizumab PK model

# Version 0.1.0

* Initial version

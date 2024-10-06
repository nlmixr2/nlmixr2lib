# nlmixr2lib

# Version 0.2.0.9000

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

     - `convertMM()` converts linear elimination to Michelis-Menton elimination

     - `transPK()` converts the `cl` style parameter transformations
       to various other PK transformations like `k`, `aob`, `alpha`,
       `k12`

  - **PD model building functions**

   - `addIndirectLin()` -- this adds an indirect effect model to a PK
     model that has a concentration `Cc` in the model.  This purposely
     uses a simple linear effect of `Cc*Ek` or `Cc*Ik` so it will be
     easy to parse and turn into other functional forms (like `Emax`
     or `Hill`).  If the PK model is not present it will use `Cc` as a
     covarite in a purely PD models.

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
* Add a variety of indirect response models
* Add a variety of tumor growth inhibition models and move all oncology models
  into a new model database directory
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

# nlmixr2lib

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

* Add Kovalenko (2020) dupilumab PK model
* Add Davda (2014) mAb consensus model
* Add Soehoel (2022) tralokinumab PK model

# Version 0.1.0

* Initial version

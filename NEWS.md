# detectseparation 0.3

## New functionality

* New `detect_infinite_estimates()` method for the `glm` function:
  `detect_infinite_estimates()` detects infinite components in the
  maximum likelihood estimates of generalized linear models with
  binomial responses, using the methods of Schwendinger et
  al. (2021) for `"log"` links and the methods in Konis (2007) for
  '"logit"', '"probit"', and '"cauchit"' links. See
  `?detect_infinite_estimates` for details

* `detect_separation()` is now a wrapper of `detect_infinite_estimates`.

## Other improvements, updates and additions

* Major rewrite and enrichment of help pages and examples.

* New vignette "Detecting separation and infinite estimates in log-binomial regression".


# detectseparation 0.2

## Other improvements, updates and additions

* `detect_separation()` returns a warnings if the link function is not
  one of `logit`, `probit`, `cloglog`, `cauchit`.
  
* Documentation and citations updates.

* Added new tests.

# detectseparation 0.1

* First public release.

# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################

## Fill the DESCRIPTION ----
## Add meta data about your application
##
## /!\ Note: if you want to change the name of your app during development,
## either re-run this function, call golem::set_golem_name(), or don't forget
## to change the name in the app_sys() function in app_config.R /!\
##
golem::fill_desc(
  pkg_name = "PSM", # The Name of the package containing the App
  pkg_title = "Non-Linear Mixed-Effects Modelling using Stochastic Differential Equations", # The Title of the package containing the App
  pkg_description = "Functions for fitting linear and non-linear mixed-effects models using
        stochastic differential equations (SDEs). The package allows for any multivariate
        non-linear time-variant model to be specified, and it also handles multidimensional
        input, covariates, missing observations, and specification of dosage regimen.
        The provided pipeline relies on the coupling of the FOCE algorithm and Kalman filtering
        as outlined by Klim et al (2009) and has been validated against the proprietary
        software 'NONMEM' (Tornoe et al, 2005). Further functions are provided for finding
        smoothed estimates of model states and for simulation.", # The Description of the package containing the App
  author_first_name = "Robert", # Your First Name
  author_last_name = "Miller", # Your Last Name
  author_email = "r.miller@phb.de", # Your Email
  repo_url = "https://github.com/waidschrat/PSM.git" # The URL of the GitHub Repo (optional)
)

# Imports: deSolve, graphics, MASS, numDeriv, stats, ucminf
# Suggests: datasets
# Author: Stig Bousgaard Mortensen, Soeren Klim, Robert Miller
# Maintainer: Robert Miller <robert.miller@tu-dresden.de>
# License: GPL (>= 2)
# URL: http://www.imm.dtu.dk/psm

## Set {golem} options ----
?golem::set_golem_options()

golem::set_golem_version(
  version = "0.8-13",
  path = pkgload::pkg_path(),
  talkative = TRUE
)
#desc::desc_set_version("0.1", file = ".", normalize = FALSE)

## Create Common Files ----
## See ?usethis for more information
usethis::use_gpl_license()
usethis::use_readme_rmd(open = FALSE)
devtools::build_readme()
# Note that `contact` is required since usethis version 2.1.5
# If your {usethis} version is older, you can remove that param
# usethis::use_code_of_conduct(contact = "Golem User")

usethis::use_cran_badge()
usethis::use_lifecycle_badge("Stable")
usethis::use_coverage(type = "codecov", repo_spec = "waidschrat/bayesROE")
# usethis::use_news_md(open = FALSE)


## Use git ----
usethis::use_git()

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()

## Favicon ----
# If you want to change the favicon (default is golem's one)
golem::use_favicon() # path = "path/to/ico". Can be an online file.
# golem::remove_favicon() # Uncomment to remove the default favicon

## Add helper functions ----
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)

# You're now set! ----

# go to dev/02_dev.R
rstudioapi::navigateToFile("dev/02_dev.R")
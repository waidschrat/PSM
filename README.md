# PSM
PSM - Population Stochastic Modelling

This R package provides functions for fitting linear and non-linear mixed-effects models using
        stochastic differential equations (SDEs). The package allows for any multivariate
        non-linear time-variant model to be specified, and it also handles multidimensional
        input, covariates, missing observations, and specification of dosage regimen.
        The provided pipeline relies on the coupling of the FOCE algorithm and Kalman filtering
        as outlined by Klim et al (2009) and has been validated against the proprietary
        software 'NONMEM' (Tornoe et al, 2005). Further functions are provided for finding
        smoothed estimates of model states and for simulation.

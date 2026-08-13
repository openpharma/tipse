# tipse 2.0

This update contains code breaking changes, specifically:

-   Updated naming conventions of tipping point analysis
    -   Hazard inflation / deflation -\> hazard multiplication
    -   Deterministic sampling -\> landmark
    -   Random sampling -\> percentile
-   Allowed imputation in both arms
-   Added heatmap plot for tipping point results in case of two-dimensional tipping points
-   Added flexibility in plotting options
-   Exported a few internal functions for custom analysis pipeline

# tipse 1.2

-   Corrected degrees of freedom calculation when applying Rubin's rule using t-distribution
-   Corrected model-based tipping point analysis to utilize all the data during imputation model building
-   Updated random seed handling
-   Improved various documentations

# tipse 1.1

-   Allowed user to suppress assess_plausibility() print messages

# tipse 1.0

-   Initial CRAN submission.

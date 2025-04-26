# MPCH
Computationally Efficient Methods for Estimating Phenome-Wide Coheritability of Multi-Type Phenotypes Using Biobank Data

Step 1: Extract data. After quality control, we obtain 290 phenotypes.

Step 2: Estimate heritability. We first estimate heritability using family data, and then using all data. This step is accomplished by parallel computation.

Step 3: Summarize the heritability results.

Step 4: Compare heritability estimates with other methods.

Step 5: Estimate coheritability. We first estimate coheritability using family data, and then using all data. This step is accomplished by parallel computation.

Step 6: Summarize the coheritability results.

Step 7: Compare coheritability estimates with other methods.

## Supplementary data
Supp Data 1: Dictionary of data fields extracted from UK biobank data. Continuous, binary, ordinal, time-to-event phenotypes, and covariates are listed in separate sheets.

Supp Data 2: Summary statistics for the 290 phenotypes in all data and in family data. We report the complete-case rate, mean, standard deviation, and number of families with observations. Especially for time-to-event phenotypes, we report the event rate, mean follow-up time, standard deviation of the follow-up time, and the number of families with events.

Supp Data 3: Estimated heritability, family-shared environmental effect, fractions of genetic/environmental effects, number of iterations, and computation time. We compare MPCH, HEc and MPCH_NE that does not model the family-shared environmental effect.

Supp Data 4: Estimated coheritability, environmental correlation, and fractions of genetic/environmental effects for pairs of phenotypes.

Supp Data 5: Estimated coheritability by HEc for continuous phenotypes.

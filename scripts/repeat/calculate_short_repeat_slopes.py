import pandas as pd
import numpy as np
import statsmodels.api as sm

# Load data
df = pd.read_excel("short_repeat_data.xlsx")  # or .csv

# Log-transform E-values
E_vals = np.array([6, 3, 1, 1e-1, 1e-3, 1e-6])
logE = np.log10(E_vals)  # convert to log10 scale

results = []

# Loop through each species
for idx, row in df.iterrows():
    species = row['Species']
    y = row.iloc[1:].values.astype(float)  # select E-value columns
    X = sm.add_constant(logE)  # add intercept for linear regression
    model = sm.OLS(y, X).fit()  # fit OLS regression
    
    slope = model.params[1]     # regression slope
    se = model.bse[1]           # standard error of slope
    ci_lower = slope - 1.96*se  # 95% confidence interval lower bound
    ci_upper = slope + 1.96*se  # 95% confidence interval upper bound
    pval = model.pvalues[1]     # p-value for slope
    
    results.append([species, slope, ci_lower, ci_upper, pval])

# Save results
df_results = pd.DataFrame(results, columns=['Species','Slope','CI_lower','CI_upper','P_value'])
df_results.to_excel("short_repeat_slopes_CI.xlsx", index=False)

print(df_results)

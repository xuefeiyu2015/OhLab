function [chi_squared, df, p_value, is_significant] = chisq_prop_test(a, b, c, d)
%CHISQ_PROP_TEST Test the difference of two proportions using the chi-squared test.
%   Inputs:
%       a - number of successes in group 1
%       b - number of failures in group 1
%       c - number of successes in group 2
%       d - number of failures in group 2
%   Outputs:
%       chi_squared - the chi-squared test statistic
%       df - degrees of freedom
%       p_value - the p-value of the test
%       is_significant - true if the difference is statistically significant at the 0.05 level, false otherwise

% Define the observed frequencies for each category
observed = [a b; c d];

% Calculate the row and column totals
row_totals = sum(observed, 2);
col_totals = sum(observed, 1);

% Calculate the total sample size
n = sum(row_totals);

% Calculate the expected frequencies under the null hypothesis of no difference
expected = row_totals * col_totals / n;

% Calculate the chi-squared statistic
chi_squared = sum(sum((observed - expected).^2 ./ expected));

% Calculate the degrees of freedom
df = (size(observed, 1) - 1) * (size(observed, 2) - 1);

% Calculate the p-value
p_value = 1 - chi2cdf(chi_squared, df);

% Determine if the difference is statistically significant
is_significant = (p_value < 0.05);
end

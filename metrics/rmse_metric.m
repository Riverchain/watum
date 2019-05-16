function rmse = rmse_metric(observed, forecasted)
%%  RMSE is a statistical metric and has applications in numerical modeling.
% 
% RMSE is used to calculate the ’’coefficient of determination’’
% which is the square of the ’’Pearson product moment correlation
% coefficient’’ (RSqr; Pearson, 1896). This metric
% comprises the squared ratio of the combined dispersion of
% two series to the total dispersion of the observed and modelled
% series. It describes the proportion of the total statistical variance
% in the observed dataset that can be explained by the model. It
% ranges from 0.0 (poor model) to 1.0 (perfect model). It records
% as a ratio the level of overall agreement between the observed
% and modelled datasets; the equation, however, is based on a consideration
% of linear relationships and limited in that it standardizes
% to the observed and modelled means and variances. The
% metric is insensitive to additive and proportional differences between
% the observed and modelled datasets, such that high
% scores can be obtained, even if the simulated values are considerably
% different from the observed values in terms of magnitude
% and variability. Indeed, since quantification is restricted to
% a consideration of differences in dispersion, a solution with systematic
% errors that over-estimated or under-estimated on each
% occasion would still produce a good result even if all of the
% numbers were wrong. The model in such cases would exhibit
% serious flaws that should, but does not, preclude it from being
% assigned a ‘‘near perfect’’ score. The metric is also oversensitive
% to outliers and thus biased towards a consideration of extreme
% events such that the true overall relationship is
% obscured. The limitations of this metric and other correlationbased
% measures are well documented (e.g. Kessler and Neas,
% 1994; Legates and Davis, 1997; Legates and McCabe, 1999);
% it was, nevertheless, still ’’common practise’’ to use such measures
% in the 1990s (Chiew and McMahon, 1993). To redress
% such quandaries it is possible to make better use of additional
% material such as the intercept and gradient of the regression
% equation upon which this metric is based. For good agreement,
% the intercept should be close to zero, and the gradient can be
% used to provide a weighted version of this metric (wRSqr;
% Krause et al., 2005). This metric is a basic statistical method
% and its output can be tested for ‘‘statistical significance’’. Testing
% would involve the use of traditional parametric procedures
% and requirements not least of which are the assumptions of a bivariate
% normal distribution and a homoscedastic relationship.

if length(observed) ~= length(forecasted)
   error('inputs are not the same size')
end

n = length(observed);
rmse = sqrt((sum((observed-forecasted).^2))/n);

end


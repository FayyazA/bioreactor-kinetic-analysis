function AIC = aic_3(residuals,n,k)
sum_squares = residuals.^2;
v = [sum(sum_squares(1,:)),sum(sum_squares(2,:)),sum(sum_squares(3,:))];
AIC = [n * log(v(1)/n) + 2 * k,n * log(v(2)/n) + 2 * k,n * log(v(3)/n) + 2 * k];
end
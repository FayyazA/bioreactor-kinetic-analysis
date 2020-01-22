function Z = overall_impedance(R,C,f)
Xc = -1/(f*2*pi*C);
Z = sqrt(R^2 + Xc^2);
end
function my_Ans = voltage_signal_analysis(time_vector,voltage_vector,current_vector,current_frequency,resistor_w_value)
%frequencies_vector = [3,10,20,51,100,201,348,1010];
%resistances_vector = [169001,122323,73802,34398,17902,9316,5407,2727];
size(voltage_vector);
size(current_vector);
resistor_w_value;
bio_impedance_vector = voltage_vector./current_vector - 2 * resistor_w_value;
coefficients_vector = [10^3,current_frequency,0];
y = @(coefficients_vector,x)coefficients_vector(1)*cos(coefficients_vector(2)*x + coefficients_vector(3));
size(time_vector);
size(bio_impedance_vector);
size(y(time_vector,coefficients_vector));
[myFit,resnorm,residual] = lsqcurvefit(y,coefficients_vector,time_vector,bio_impedance_vector,[0,current_frequency - 0.00001,-pi],[10^9,current_frequency + 0.00001,pi],optimset('TolFun',1e-10));
disp(resnorm)
figure
plot(time_vector,myFit(1) * cos(myFit(2) * time_vector + myFit(3)),'g-',time_vector,bio_impedance_vector,'b-')
R = myFit(1) / sqrt(1 + (tan(myFit(3)) ^ 2));
Xc = R * tan(myFit(3));
my_Ans = [R,Xc];
end
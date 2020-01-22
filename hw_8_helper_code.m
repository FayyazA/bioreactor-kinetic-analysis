t = [0:0.00005:1];
frequencies_vector = [3,10,20,51,100,201,348,1010];
Rb = [10,100,10^3,10^4,10^5,10^6,10^6,10^6];
Xbc = [10^3,10^2,10^1,10^-1,10^-3,10^-5,10^-7,10^-9];
Xbc = -1 ./ (frequencies_vector .* Xbc * 2 * pi)
Relectrode = [169001,122323,73802,34398,17902,9316,5407,2727];
for j = 1:8
input_vector = 1*10^-6*cos(frequencies_vector(j) * 2 * pi * t);
resistor_w_value = Relectrode(j);
output_vector =  (sqrt(Rb(j)^2 + Xbc(j)^2) * cos(frequencies_vector(j)*2*pi*t + atan(Xbc(j)/Rb(j))) + 2 * Relectrode(j)) .* input_vector;
disp(voltage_signal_analysis(t,output_vector,input_vector,frequencies_vector(j)*2*pi,resistor_w_value))
figure
plot(t,output_vector);
xlabel('time')
ylabel('voltage')

end
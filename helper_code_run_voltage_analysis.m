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
output_vector = output_vector; %+ 10^-6 * randn(1,length(output_vector));
smoothed = smooth(output_vector)';
disp(voltage_signal_analysis(t,output_vector,input_vector,frequencies_vector(j)*2*pi,resistor_w_value))
omega = frequencies_vector(j);
output_vector_times_transfer_function = output_vector * (1000 * omega * 1i / (2000 * pi)) / (1 - omega^2/(20000*pi^2) + omega * 1i * (1/(2000*pi) + 1/(10*pi)));
denom = (1 - omega^2/(20000*pi^2))^2 + omega^2 * (1/(2000*pi) + 1/(10*pi))^2;
num_real = (omega^2 / (2000*pi)) * (1/(2000*pi) + 1/(10*pi));
num_imaginary = (1 - omega^2 / (20000*pi^2)) * (omega/(2000*pi));
transfer_function_mag = sqrt((1000*num_real/denom)^2 + (1000*num_imaginary/denom)^2);
transfer_function_phase = atan(num_imaginary/num_real);
transfer_function = transfer_function_mag * cos(omega*t + transfer_function_phase);
final_output_voltage = output_vector .* transfer_function;
%figure
%plot(t,output_vector)
%xlabel('time')
%coefficients_vector = [max(output_vector),2*pi,0,max(output_vector),2*pi,0,10];
%y = @(coefficients_vector,x)coefficients_vector(1)*cos(coefficients_vector(2)*x + coefficients_vector(3)) + coefficients_vector(4) * sin(coefficients_vector(5) * x + coefficients_vector(6)) + coefficients_vector(7)*cos(omega*2*pi*x);
%[myFit,resnorm,residual] = lsqcurvefit(y,coefficients_vector,t,output_vector,[-50,0,-pi,-50,0,-pi,-50],[50,10000*2*pi,pi,50,10000*2*pi,pi,50],optimset('TolFun',1e-15,'MaxFunEval',10^4,'MaxIter',10^4));
%disp(resnorm)
figure
plot(t,output_vector);%,'g-',t,myFit(1)*cos(myFit(2)*t+myFit(3)) + myFit(4)*cos(myFit(5)*t+myFit(6)) + myFit(7)*cos(omega*2*pi*t),'b-')
xlabel('time')
ylabel('voltage')

end
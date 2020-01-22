function F = myfun(data,transfer,frequency,vector,Rb,Cb,Rdl,Cdl,Rlyte,R,C)
w = frequency*2*pi;
Zo = R-1i/(w*C); %R* ( (1 + 1i * w * C * R) / (1i * w * C * R));
Zb = Rb-1i/(w*Cb); %Rb * ( (1 + 1i * w * Cb * Rb) / (1i * w * Cb * Rb));
Zel = Rlyte + (Rdl-1i*w*Cdl*Rdl^2)/(1+(w^2) * (Cdl^2) * (Rdl^2)); %(Rdl / (1 + 1i * w * Cdl * Rdl)) + Rlyte;
Top = (1/Zo + 1/(Zb+Zel))^-1;
Bot = (1/Zel + 1/Zo)^-1;
Impedance_waveform = Top + Bot;
time = data(:,1);
voltage = data(:,3);
[mag,phase,w0] = bode(transfer,frequency*2*pi);
y = @(cv,x)cv(1)*-cos(cv(2)*2*pi*x+cv(3));
[fit,resnorm,residuals] = lsqcurvefit(y,vector,time,voltage);
disp(resnorm)
figure
plot(time,voltage,'g-',time,fit(1)*-cos(fit(2)*2*pi*time+fit(3)));
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Output Data')
legend('Fitted Data','Raw Data')
figure
plot(time,data(:,2),'g-',time,data(:,2) * -0.35,'b-')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Input Data')
legend('Voltage on Bottom Side of Current Source', 'Voltage on Top Side of Current Source')
v_of_t = [fit(1) / mag, fit(3) - phase * 2 * pi /360];
z_from_data = [v_of_t(1) * 10^6,v_of_t(2)];
%Xc = -1/(frequency*2*pi*capacitance);
%Xc = 0;
%z_known = [sqrt(resistance^2 + Xc^2),atan(Xc/resistance)];
z_known = [abs(Impedance_waveform),angle(Impedance_waveform)];
disp(z_known)
disp(Impedance_waveform)
disp(z_from_data)
disp([mag,phase*2*pi/360])
disp(fit)
F = [z_known,z_from_data];
end
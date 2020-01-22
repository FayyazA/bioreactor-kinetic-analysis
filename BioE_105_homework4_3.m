%% Question 3
t = 0:1e-4:10;
one_khz = 2 * sin(2*pi*1*t) + 2;
data1 = sim_rc(one_khz,0.22*10^-6);
ten_khz = 2 * sin(2*pi*10*t) + 2;
data2 = sim_rc(ten_khz,0.22*10^-6);
plot(t,data1.Vout,'g-',t,data2.Vout,'b-');
    legend('1 Khz Input', '10 Khz Input')
    ylabel('Voltage across C1 (V)')
    xlabel('Time (ms)')
%% The 10 Khz signal overall has a much smoother signal with a smaller amplitude of oscillation that much more accurately follows the overall charging pattern of the capacitor.  The 1 Khz signal has a much larger amplitude of oscillation of about 0.25 V whereas the 10 Khz one has a largely negligible amplitude of about 0.01 - 0.05 V.
combined = one_khz + ten_khz;
data_combined = sim_rc(combined,0.22*10^-6);
figure
subplot(2,1,1)
plot(t,combined,'g-'); 
    legend('Input Signal');
    ylabel('Input Voltage (V)')
    xlabel('Time (ms)')
subplot(2,1,2)
plot(t,data_combined.Vout,'g-'); 
    legend('Vout_capacitor');
    ylabel('Voltage across C1 (V)')
    xlabel('Time (ms)')
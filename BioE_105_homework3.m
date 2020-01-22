function BioE_105_homework3()
t = 0:0.005:0.35;
signal = 10 * (1-exp(-t/0.1));
amplitudes = [0.25,1,3];
strings = ["Noise Amplitude = 0.25", "Noise Amplitude = 1", "Noise Amplitude = 3"];
figure
plot(t,signal)
title("Transient Analysis")
xlabel('Time (s)')
ylabel('Voltage (Volts)')
print(gcf,'-dtiff','-r300',strcat('TransientAnalysis',".tif"))
for i = 1:3
    signal_with_noise = signal + randn(1,71) * amplitudes(i)/3;
    params = [10,0.1];
    Vo = @(params,t)params(1)-params(1)*exp(-t/params(2));
    [x,resnorm,residuals] = lsqcurvefit(Vo,params,t,signal_with_noise);
    rmse = (sum(residuals.^2)/71)^0.5;
    fitted_signal = x(1) * (1-exp(-t/x(2)));
    figure
    plot(t,signal_with_noise,'k*--',t,fitted_signal,'r-')
    legend('Data','Fitted Graph')
    title(strings(i))
    xlabel('Time (s)')
    ylabel('Voltage (Volts)')
    print(gcf,'-dtiff','-r300',strcat(strings(i),".tif"))
    disp(x)
    disp(rmse)
end
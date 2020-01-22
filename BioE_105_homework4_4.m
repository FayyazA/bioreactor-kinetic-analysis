%% Question 4
t = 0:1e-4:10;
freq_vector = [1:10:1001];
table = zeros(length(freq_vector));
for i = 1:length(freq_vector)
    x = 2 * sin(2*pi*freq_vector(i)*t) + 2;
    data = sim_rc(x,3.18*10^-11);
    second_half = data.Vout(50000:length(data.Vout));
    amplitude = max(second_half) - min(second_half);
    table(i) = amplitude;
end
plot(log10(freq_vector),20*log10(table/4))
ylabel('Signal Amplitude (LogVolts)')
    xlabel('Frequency (LogHertz)')
%%This is a low-pass filter because the ouput signal amplitude drops off dramatically as the frequency of the input signal increases.  It allows the low frequency signals to pass through whereas the high frequencies don't, so it's a low-pass filter.  
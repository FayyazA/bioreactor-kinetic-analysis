%%Testing sim_rccopy response
t = 0:1e-4:10;
x = 2 * square(2*pi*1*t) + 2;
possible_C = [.00047, .001, .0022, .0047, .01, .022]*10^(-6);
for i = 1:6
    data2 = sim_rccopy(x,possible_C(i));
    figure
    plot(t,data2.Vout)
    ylabel('Voltage across C1 (V)')
    xlabel('Time (ms)')
    axis([0,10, -0.5,4.5])
end
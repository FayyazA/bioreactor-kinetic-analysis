%%Question 1a: Example where capacitor doesn't fully saturate
t = 0:1e-4:10;
x = 2 * square(2*pi*1*t) + 2;
data = sim_rc(x,0.022*10^-6)
figure
plot(t,data.Vout)
ylabel('Voltage across C1 (V)')
xlabel('Time (ms)')
axis([0,10, -0.5,4.5])
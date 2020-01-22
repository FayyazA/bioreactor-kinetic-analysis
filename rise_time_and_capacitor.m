function rt_and_c = rise_time_and_capacitor(signal,t)
%%Function To calculate Rise Time and Capacitor Value of Circuit
times = [0,0];
counter = 1;
for i = 1:10000
    if abs(signal(i) - max(signal) * 0.1) <= 0.005 * max(signal)
        times(1) = t(i);
    end
    if abs(signal(i) - max(signal) * 0.9) <= 0.005 * max(signal)
        times(2) = t(i);
        
    end
    
end

rise_time = abs(times(2) - times(1));
tau = rise_time * 10^-3/log(9);
capacitor = tau/(10^4);
rt_and_c = [rise_time,capacitor];
end
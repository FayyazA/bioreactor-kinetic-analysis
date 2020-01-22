%%Question 1b: Testing Rise Time function
t = 0:1e-4:10;
x = 2 * square(2*pi*1*t) + 2;
result_table = zeros(6,2);
possible_C = [.00047, .001, .0022, .0047, .01, .022]*10^(-6);
for i = 1:6
    data2 = sim_rc(x,possible_C(i));
    rt_and_c = rise_time_and_capacitor(data2.Vout,t);
    result_table(i,:) = rt_and_c;
end
result_table
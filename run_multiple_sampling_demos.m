function result_vector = run_multiple_sampling_demos()
frequencies = [1000,1250,1500,1750,2000,3000,4000,5000,7500,10000,20000,30000,40000,49950];
result_vector = zeros(length(frequencies));
residuals = zeros(length(frequencies));
for i = 1:length(frequencies)
    frequencies = [1000,1250,1500,1750,2000,3000,4000,5000,7500,10000,20000,30000,40000,49950];
    hertz = frequencies(i);
    outputs = sampling_demo_function_version(hertz);
    output_signal = outputs(2)
    data = output_signal{1,1};
    counter = 0;
    time_index = 1
    while counter < 3
        for j = 1:41
            if abs(data(j) - 0)  <= 0.005 & counter < 3
                counter = counter + 1
                time_index = j
            end
        end
    end
    period = time_index * 0.0005 - 0.0005;
    result_vector(i) = 1/period;
end
figure
plot(frequencies,result_vector(:,1));
title('Reconstructed Frequency vs. Input Frequency')
xlabel('Input Frequency (Hertz)')
ylabel('Reconstructed Frequency (Hertz)')
end
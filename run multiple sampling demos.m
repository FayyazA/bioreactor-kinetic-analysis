frequencies = [1000,1250,1500,1750,2000,3000,4000,5000,7500,10000,20000,30000,40000,49950];
for i = 1:length(frequencies)
    outputs = sampling_demo_function_version(frequencies(i));
    input_signal = outputs(1)
    output_signal = outputs(2)
    figure
    plot([0:0.0005:0.02],input_signal(1:41),'b-',[0:0.0005:0.02],output_signal(1:41),'g-')
end
function errors = calculate_errors(input)
mag_error = (input(3)-input(1))/ input(1);
phase_error = (input(4)-input(2))/ input(2);
errors = [mag_error,phase_error];
end
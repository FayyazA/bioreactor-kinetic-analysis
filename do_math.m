function my_ans = do_math(n)
syms T w
%my_ans = (1/(2*pi*1i*n^2)) *( exp(1i*n*pi/T) * (1i*n*pi/T - 1) - exp(-1i*n*pi/T) * (-1i*n*pi/T - 1)) * exp(-1i * w * n);
my_ans = (1/(2*pi*1i*T*n^2)) * (2 + 2*1i*n*pi) * cos(n*pi) * exp(-1i*w*n);
%simplify(my_ans,1000);
my_ans = vpa(my_ans,3);
end
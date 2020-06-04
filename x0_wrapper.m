%clear all;
close all;
R1P = 0.020055852;
kPL = 0.0187;
FP = 0.12795384;
R1Lin = 0.06091017;
kLP = 0.00021363;
kMCT4 = 0.59140871;
R1Lex = 0.027305961;
FL = 1.463854452;
kMCT1 = 0.008302768;
k = 3.862436329;
theta = 16.83617857;
gamma = 410169604.5;
%x0_vec = [0.019607843,0.020025031,0.020442219,0.020859408,0.021276596];
diff_in_fit = zeros(120,18);
my_vec = ["R1P","kPL","FP","R1Lin","kLP","kMCT4","R1Lex","FL","kMCT1","k","theta","gamma"];
table_data = zeros(1,18);
K = [-R1P-kPL-FP,kLP,0;kPL,-R1Lin-kLP-kMCT4,kMCT1;0,kMCT4,-R1Lex-kMCT1-FL];
A = [0,0,0,k,theta,gamma];
hk_1 = generate_synthetic(K,A);
count = 1;
for i = 1:length(my_vec)
    for p = 1:10
        hk_1_fits = run_multiple_files_modeling_synthetic_data_x0_testing(hk_1,"x0_test",thing(p,i),my_vec(i));
        hk_1_actual = [R1P,kPL,FP,R1Lin,kLP,kMCT4,R1Lex,FL,kMCT1,A(4),A(5),A(6),0,0,0,0,0,0];
        diff_in_fit(count,:) = hk_1_fits-hk_1_actual;
        count = count + 1;
    end
end
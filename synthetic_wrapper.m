% clear all;
% close all;
% R1P = 0.020781;
% kPL = 0.003479516;
% FP = 0.12115537;
% R1Lin = 0.035353019;
% kLP = 2.22045E-14;
% kMCT4 = 0.221997574;
% R1Lex = 0.026802025;
% FL = 0.119317026;
% kMCT1 = 2.22E-14;
diff_in_fit = zeros(60,18);
for i = 1:60
K = [-R1P(i)-kPL(i)-FP(i),kLP(i),0;kPL(i),-R1Lin(i)-kLP(i)-kMCT4(i),kMCT1(i);0,kMCT4(i),-R1Lex(i)-kMCT1(i)-FL(i)];
A = [0,0,0,k(i),theta(i),gamma(i)];
hk_1 = generate_synthetic(K,A);
hk_1_fits = run_multiple_files_modeling_synthetic_data(hk_1, "hk_1");
hk_1_actual = [R1P(i),kPL(i),FP(i),R1Lin(i),kLP(i),kMCT4(i),R1Lex(i),FL(i),kMCT1(i),A(4),A(5),A(6),0,0,0,0,0,0];
diff_in_fit(i,:) = hk_1_fits-hk_1_actual;
end
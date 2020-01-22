clear all
close all
file_import = xlsread("C:\Users\Owner\Documents\JW330.xlsx");


pyrflip = file_import(3,:);
lacflip = file_import(4,:);
TR = 3;
pyruvate = file_import(1,:);
lactate = file_import(2,:);
save("JW330",'pyrflip','lacflip','TR','pyruvate','lactate')
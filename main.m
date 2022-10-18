clear all;
close all;
clc;

type = 1;
prob = probSet(type);

soln = numericalMethods();

soln = soln.solveNumerical(prob);

err = errorCompute(prob, soln);
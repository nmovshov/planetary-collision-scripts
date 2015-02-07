%% Test drive bound_mass.m and bound_mass.py
clear
close all
clc

%% Data to work on
fnl = readtable('FNL_4.txt');
pos = fnl{:,1:3};
vel = fnl{:,4:6};
m   = fnl{:,end};

%% Compare algorithms in bound_mass.m
[M1, ind1] = bound_mass(pos,vel,m,'kory');
fprintf('Found %d bound nodes massing %g kg.\n',sum(ind1),M1);

[M2, ind2] = bound_mass(pos,vel,m,'jutzi');
fprintf('Found %d bound nodes massing %g kg.\n',sum(ind2),M2);

[M3, ind3] = bound_mass(pos,vel,m,'naor');
fprintf('Found %d bound nodes massing %g kg.\n',sum(ind3),M3);

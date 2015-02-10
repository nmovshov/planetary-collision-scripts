%% Generate synthetic particle distributions to test clumping algorithm
clear
close all
clc

%% A grid of identical stationary particles
x = linspace(-10,10,21);
[X, Y, Z] = meshgrid(x);
X = X(:); Y = Y(:); Z = Z(:);
VX = zeros(size(X));
VY = VX;
VZ = VX;
M = ones(size(X));
FNL_1 = table(X, Y, Z, VX, VY, VZ, M);
save synth_FNLs FNL_1
writetable(FNL_1,'FNL_1.txt')

%% A smaller grid within a grid
clear
x = linspace(-5,5,5);
[X, Y, Z] = meshgrid(x);
inX = X(:); inY = Y(:); inZ = Z(:);
inV = zeros(size(inX));
inM = ones(size(inX));
x = linspace(-30,30,11);
[X, Y, Z] = meshgrid(x);
outX = X(:); outY = Y(:); outZ = Z(:);
mask = outX>-20 & outX<20 & outY>-20 & outY<20 & outZ>-20 & outZ<20;
outX(mask) = [];
outY(mask) = [];
outZ(mask) = [];
outV = zeros(size(outX));
outM = ones(size(outX));
T1 = table(inX,inY,inZ,inV,inV,inV,inM,...
             'VariableNames',{'X','Y','Z','VX','VY','VZ','M'});
T2 = table(outX,outY,outZ,outV,outV,outV,outM,...
             'VariableNames',{'X','Y','Z','VX','VY','VZ','M'});
FNL_2 = [T1;T2];
save synth_FNLs FNL_2 -append
writetable(FNL_2,'FNL_2.txt')

%% Static grid within a flying grid
clear
x = linspace(-5,5,5);
[X, Y, Z] = meshgrid(x);
inX = X(:); inY = Y(:); inZ = Z(:);
inV = zeros(size(inX));
inM = ones(size(inX));
x = linspace(-30,30,11);
[X, Y, Z] = meshgrid(x);
outX = X(:); outY = Y(:); outZ = Z(:);
mask = outX>-20 & outX<20 & outY>-20 & outY<20 & outZ>-20 & outZ<20;
outX(mask) = [];
outY(mask) = [];
outZ(mask) = [];
outV = ones(size(outX));
outM = ones(size(outX));
T1 = table(inX,inY,inZ,inV,inV,inV,inM,...
             'VariableNames',{'X','Y','Z','VX','VY','VZ','M'});
T2 = table(outX,outY,outZ,0*outV,1*outV,0*outV,outM,...
             'VariableNames',{'X','Y','Z','VX','VY','VZ','M'});
FNL_3 = [T1;T2];
save synth_FNLs FNL_3 -append
writetable(FNL_3,'FNL_3.txt')

%% Slow grid within a fast grid
clear
offset = 100;
x = offset + linspace(-5,5,5);
[X, Y, Z] = meshgrid(x);
inX = X(:); inY = Y(:); inZ = Z(:);
inVX = 1 + randn(size(inX));
inVY = 1 + randn(size(inX));
inVZ = 1 + randn(size(inX));
inM = ones(size(inX));
x = offset + linspace(-30,30,11);
[X, Y, Z] = meshgrid(x);
outX = X(:); outY = Y(:); outZ = Z(:);
mask = outX>offset-20 & outX<offset+20 & outY>offset-20 & outY<offset+20 & outZ>offset-20 & outZ<offset+20;
outX(mask) = [];
outY(mask) = [];
outZ(mask) = [];
outVX = 10 + randn(size(outX));
outVY = 10 + randn(size(outX));
outVZ = 10 + randn(size(outX));
outM = ones(size(outX));
T1 = table(inX,inY,inZ,inVX,inVY,inVZ,inM,...
             'VariableNames',{'X','Y','Z','VX','VY','VZ','M'});
T2 = table(outX,outY,outZ,outVX,outVY,outVZ,outM,...
             'VariableNames',{'X','Y','Z','VX','VY','VZ','M'});
FNL_4 = [T1;T2];
save synth_FNLs FNL_4 -append
writetable(FNL_4,'FNL_4.txt')
FNL_4a = table2array(FNL_4);
save('FNL_4.fnl','FNL_4a','-ascii')

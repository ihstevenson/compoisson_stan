

% Need to install MatlabStan (dev-branch)...
% https://github.com/brian-lau/MatlabStan
% https://github.com/brian-lau/MatlabStan/wiki/Getting-Started
% and adjust paths

addpath(genpath('C:\Users\ian\Documents\MATLAB\MatlabProcessManager-master'))
addpath(genpath('C:\Users\ian\Documents\MATLAB\MatlabStan-dev'))
addpath('fastBSpline\')

model(1) = StanModel('file','stan\poisson_regression.stan');
model(1).compile();

model(2) = StanModel('file','stan\negbino_regression.stan');
model(2).compile();

model(3) = StanModel('file','stan\compoisson_regression.stan');
model(3).compile();

model(4) = StanModel('file','stan\gencount_regression.stan');
model(4).compile();

for i=1:length(model)
    model(i).working_dir = 'C:\Users\ian\Documents\MATLAB\';
end
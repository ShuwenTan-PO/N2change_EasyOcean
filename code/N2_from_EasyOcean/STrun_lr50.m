%
% sensitivity test, linear regression 50 m
%
clear all
close all
fl = 500; 
lr = 50;

A = {'75N', 'A02', 'A03', 'A05', 'A10', 'A12', 'A13', 'A16-A23', 'A20', 'A22', 'AR07E', 'AR07W'};
I = {'I01', 'I02', 'I03-I04', 'I05', 'I06S', 'I07', 'I08N', 'I08S-I09N', 'I09S', 'I10', 'IR06', 'IR06E', 'IR06-I10'};
P = {'P01', 'P02', 'P03', 'P04', 'P06', 'P09', 'P10', 'P11', 'P13', 'P14', 'P15', 'P16', 'P17', 'P17E', 'P18', 'P21'};
S = {'S04I', 'S04P', 'SR01', 'SR03', 'SR04'};
oceans = {'atlantic', 'indian', 'pacific', 'southern'};

reported_file = '/Users/stan/Desktop/WORK/DATA/GOSHIP/EasyOcean/reported/';
gridded_file = '/Users/stan/Desktop/WORK/DATA/GOSHIP/EasyOcean/gridded/';
easyocean_file = './';

% infos:
% A
A16A23.years = {'1988', '1993', '1995', '2003', '2011', '2013'}';
A16A23.ll_grid = [-72.5:0.1:63.4];
A16A23.pr_grid = [0:10:6500];

% S
S04P.years = {'1992', '2011', '2018'};
S04P.ll_grid = [168.0:0.1:287.1];
S04P.pr_grid = [0:10:6500];

D_pr = [];
% A
n = 8;
%     tic;
clear D_pr
sec = A{n};
disp(strcat('-----',sec,'-----'))
fname = lower(sec);
com = ['load ' reported_file oceans{1} '/' sec '/' fname '.mat'];
eval(com);
if n==1
    years = eval(['A' sec]).years;
    ll_grid = eval(['A' sec '.ll_grid']);
    pr_grid = eval(['A' sec '.pr_grid']);
elseif n==8
    years = eval([sec(1:3) sec(5:7)]).years;
    ll_grid = eval([sec(1:3) sec(5:7) '.ll_grid']);
    pr_grid = eval([sec(1:3) sec(5:7) '.pr_grid']);
else
    years = eval(sec).years;
    ll_grid = eval([sec '.ll_grid']);
    pr_grid = eval([sec '.pr_grid']);
end
for i = 1:length(years)      
    disp(['gridded - MAT, ', years{i}]);
    com = ['[s, m] = copyfile(''' easyocean_file sec '/' 'configuration_' years{i} '.m'', ''configuration.m'');'];
    eval(com);
    if s ~= 1, error(['copyfile ' num2str(i) ':', m]); end
    com = ['[D_pr(' num2str(i) ')' ,', pr_grid_lr] ', '= grid_data_pressure(D_reported(' num2str(i) '), ll_grid, pr_grid, ', num2str(fl), ',', num2str(lr) ');'];
    eval(com);
end
if ~exist(strcat([gridded_file oceans{1} '/' sec '/']), 'dir')
    mkdir(strcat([gridded_file oceans{1} '/' sec '/']))
end
com = ['save ' gridded_file oceans{1} '/' sec '/' fname '_lr50.mat' ' D_pr ll_grid pr_grid pr_grid_lr'];
eval (com);
%     toc;


% S
n = 2;
%     tic;
clear D_pr
sec = S{n};
disp(strcat('-----',sec,'-----'))
fname = lower(sec);
com = ['load ' reported_file oceans{4} '/' sec '/' fname '.mat'];
eval(com);
years = eval(sec).years;
ll_grid = eval([sec '.ll_grid']);
pr_grid = eval([sec '.pr_grid']);
for i = 1:length(years)      
    disp(['gridded - MAT, ', years{i}]);
    com = ['[s, m] = copyfile(''' easyocean_file sec '/' 'configuration_' years{i} '.m'', ''configuration.m'');'];
    eval(com);
    if s ~= 1, error(['copyfile ' num2str(i) ':', m]); end
    com = ['[D_pr(' num2str(i) ')' ,', pr_grid_lr] ', '= grid_data_pressure(D_reported(' num2str(i) '), ll_grid, pr_grid, ', num2str(fl), ',', num2str(lr) ');'];
    eval(com);
end
if ~exist(strcat([gridded_file oceans{4} '/' sec '/']), 'dir')
    mkdir(strcat([gridded_file oceans{4} '/' sec '/']))
end
com = ['save ' gridded_file oceans{4} '/' sec '/' fname '_lr50.mat' ' D_pr ll_grid pr_grid pr_grid_lr'];
eval (com);
%     toc;

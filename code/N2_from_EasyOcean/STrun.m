%
% do all sections
% use modeified grid_data_pressure, compute N2s, and save in .mat
%
clear all
close all
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
A75N.years = {'1994', '1995', '1997', '1998', '1999', '2000', '2001', '2006', '2016'};
A75N.ll_grid = [-16.0:0.1:18.0];
A75N.pr_grid = [0:10:6500];
A02.years = {'1994', '1997', '2017'};
A02.ll_grid = [309.9:0.1:349.4];
A02.pr_grid = [0:10:6500];
A03.years = {'1993'};
A03.ll_grid = [286.3:0.1:351.5];
A03.pr_grid = [0:10:6500];
A05.years = {'1992', '1998', '2004', '2010', '2011', '2015', '2020'};
A05.ll_grid = [279.9:0.1:346.8];
A05.pr_grid = [0:10:6500];
A10.years = {'1992', '2003', '2011'};
A10.ll_grid = [-53.5:0.1:15.1];
A10.pr_grid = [0:10:6500];
A12.years = {'1992', '1996', '1999', '2000', '2002', '2005', '2008a', '2008b', '2010', '2014'};
A12.ll_grid = [-69.8:0.1:-34.1];
A12.pr_grid = [0:10:6500];
A13.years = {'1983', '2010'};
A13.ll_grid = [-69.4:0.1:4.9];
A13.pr_grid = [0:10:6500];
A16A23.years = {'1988', '1993', '1995', '2003', '2011', '2013'}';
A16A23.ll_grid = [-72.5:0.1:63.4];
A16A23.pr_grid = [0:10:6500];
A20.years = {'1997', '2003', '2012', '2021'};
A20.ll_grid = [6.8:0.1:43.3];
A20.pr_grid = [0:10:6500];
A22.years = {'1997', '2003', '2012', '2021'};
A22.ll_grid = [11.0:0.1:40.1];
A22.pr_grid = [0:10:6500];
AR07E.years = {'1990', '1991a', '1991b', '1991c', '1992', '1994', '1995', ...
         '1996', '1997', '2000', '2005', '2014a', '2014b', '2015'};
AR07E.ll_grid = [314.3:0.1:346.0];
AR07E.pr_grid = [0:10:6500];
AR07W.years = {'1990', '1992', '1993', '1994', '1995', '1996', '1997', '1998', ...
         '2011b', '2012', '2013a', '2013b'};
AR07W.ll_grid = [54.9:0.1:60.6];
AR07W.pr_grid = [0:10:6500];

% I
I01.years = {'1995a', '1995b'};
I01.ll_grid = [50.5:0.1:97.6];
I01.pr_grid = [0:10:6500];
I02.years = {'1995', '2000'};
I02.ll_grid = [39.8:0.1:105.7];
I02.pr_grid = [0:10:6500];
I03I04.years = {'1995', '2003'};
I03I04.ll_grid = [35.3:0.1:115.8];
I03I04.pr_grid = [0:10:6500];
I05.years = {'1987', '1995', '2002', '2009'};
I05.ll_grid = [30.3:0.1:115.4];
I05.pr_grid = [0:10:6500];
I06S.years = {'1993', '1996', '2008', '2019'};
I06S.ll_grid = [-69.1:0.1:-33.1];
I06S.pr_grid = [0:10:6500];
I07.years = {'1995', '2018'};
I07.ll_grid = [-65.4:0.1:19.9];
I07.pr_grid = [0:10:6500];
I08N.years = {'1995', '2019'};
I08N.ll_grid = [-24.0:0.1:5.9];
I08N.pr_grid = [0:10:6500];
I08SI09N.years = {'1995', '2007', '2016'};
I08SI09N.ll_grid = [-65.9:0.1:19.8];
I08SI09N.pr_grid = [0:10:6500];
I09S.years = {'1995', '2004', '2012'};
I09S.ll_grid = [-65.4:0.1:-34.4];
I09S.pr_grid = [0:10:6500];
I10.years = {'1995', '2015'};
I10.ll_grid = [-24.4:0.1:-8.6];
I10.pr_grid = [0:10:6500];
IR06.years = {'1995a', '1995b', '2000'};
IR06.ll_grid = [-24.4:0.1:-8.9];
IR06.pr_grid = [0:10:6500];
IR06E.years = {'1989', '1992', '2000'};
IR06E.ll_grid = [-19.7:0.1:-8.8];
IR06E.pr_grid = [0:10:6500];
IR06I10.years = {'1995a', '1995b', '1995', '2000', '2015'};
IR06I10.ll_grid = [110.4:0.1:112.8];
IR06I10.pr_grid = [0:10:6500];

% P
P01.years = {'1985', '1999', '2007', '2014'};
P01.ll_grid = [145.4:0.1:235.1];
P01.pr_grid = [0:10:6500];
P02.years = {'1993', '2004', '2013'};
P02.ll_grid = [133.0:0.1:242.7];
P02.pr_grid = [0:10:6500];
P03.years = {'1985', '2005'};
P03.ll_grid = [124.9:0.1:242.7];
P03.pr_grid = [0:10:6500];
P04.years = {'1989'};
P04.ll_grid = [126.5:0.1:274.3];
P04.pr_grid = [0:10:6500];
P06.years = {'1992', '2003', '2010', '2017'}';
P06.ll_grid = [153.4:0.1:288.6];
P06.pr_grid = [0:10:6500];
P09.years = {'1994', '2010', '2016'};
P09.ll_grid = [-2.87:0.1:34.25];
P09.pr_grid = [0:10:6500];
P10.years = {'1993', '2005', '2011'};
P10.ll_grid = [-4.02:0.1:43.26];
P10.pr_grid = [0:10:6500];
P11.years = {'1993'};
P11.ll_grid = [-65.9:0.1:-11.7];
P11.pr_grid = [0:10:6500];
P13.years = {'1991', '1992', '2011'};
P13.ll_grid = [-8.0:0.1:52.1];
P13.pr_grid = [0:10:6500];
P14.years = {'1992', '2007'};
P14.ll_grid = [-66.02:0.1:59.1];
P14.pr_grid = [0:10:6500];
P15.years = {'1990', '2001', '2009', '2016'};
P15.ll_grid = [-74.6:0.1:54];
P15.pr_grid = [0:10:6500];
P16.years = {'1992', '2006', '2015'};
P16.ll_grid = [-70.1:0.1:56.3];
P16.pr_grid = [0:10:6500];
P17.years = {'1991', '2001'};
P17.ll_grid = [-62.5:0.1:55.9];
P17.pr_grid = [0:10:6500];
P17E.years = {'1992', '2017'};
P17E.ll_grid = [-66.4:0.1:-52.52];
P17E.pr_grid = [0:10:6500];
P18.years = {'1994', '2007', '2016'};
P18.ll_grid = [-70.1:0.1:22.9];
P18.pr_grid = [0:10:6500];
P21.years = {'1994', '2009'};
P21.ll_grid = [153.6:0.1:284.9];
P21.pr_grid = [0:10:6500];

% S
S04I.years = {'1994', '2012'};
S04I.ll_grid = [20.0:0.1:162.3];
S04I.pr_grid = [0:10:6500];
S04P.years = {'1992', '2011', '2018'};
S04P.ll_grid = [168.0:0.1:287.1];
S04P.pr_grid = [0:10:6500];
SR01.years = {'1993', '1994', '1996', '1997', '2000', '2001', '2002', '2003', '2009a', '2009b', '2011', '2015a', '2015b', '2016'};
SR01.ll_grid = [-61.1:0.1:-54.6];
SR01.pr_grid = [0:10:6500];
SR03.years = {'1991', '1993', '1994a', '1994b', '1995', '1996', '2001', '2008', '2011'};
SR03.ll_grid = [-65.9:0.1:-43.9];
SR03.pr_grid = [0:10:6500];
SR04.years = {'1989', '1990', '1992', '1996', '1998', '2005', '2008', '2010'};
SR04.ll_grid = [303.4:0.1:348.7];
SR04.pr_grid = [0:10:6500];

D_pr = [];
% A
for n = 1:length(A)
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
        com = ['[D_pr(' num2str(i) ')' ,', pr_grid_lr] ', '= grid_data_pressure(D_reported(' num2str(i) '), ll_grid, pr_grid);'];
        eval(com);
    end
    if ~exist(strcat([gridded_file oceans{1} '/' sec '/']), 'dir')
        mkdir(strcat([gridded_file oceans{1} '/' sec '/']))
    end
    com = ['save ' gridded_file oceans{1} '/' sec '/' fname '.mat' ' D_pr ll_grid pr_grid pr_grid_lr'];
    eval (com);
%     toc;
end

% I
for n = 1:length(I)
%     tic;
    clear D_pr
    sec = I{n};
    disp(strcat('-----',sec,'-----'))
    fname = lower(sec);
    com = ['load ' reported_file oceans{2} '/' sec '/' fname '.mat'];
    eval(com);
    if n==3
        years = eval([sec(1:3) sec(5:7)]).years;
        ll_grid = eval([sec(1:3) sec(5:7) '.ll_grid']);
        pr_grid = eval([sec(1:3) sec(5:7) '.pr_grid']);
    elseif n==8
        years = eval([sec(1:4) sec(6:9)]).years;
        ll_grid = eval([sec(1:4) sec(6:9) '.ll_grid']);
        pr_grid = eval([sec(1:4) sec(6:9) '.pr_grid']);
    elseif n==13
        years = eval([sec(1:4) sec(6:8)]).years;
        ll_grid = eval([sec(1:4) sec(6:8) '.ll_grid']);
        pr_grid = eval([sec(1:4) sec(6:8) '.pr_grid']);
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
        com = ['[D_pr(' num2str(i) ')' ,', pr_grid_lr] ', '= grid_data_pressure(D_reported(' num2str(i) '), ll_grid, pr_grid);'];
        eval(com);
    end
    if ~exist(strcat([gridded_file oceans{2} '/' sec '/']), 'dir')
        mkdir(strcat([gridded_file oceans{2} '/' sec '/']))
    end
    com = ['save ' gridded_file oceans{2} '/' sec '/' fname '.mat' ' D_pr ll_grid pr_grid pr_grid_lr'];
    eval (com);
%     toc;
end

% P
for n = 1:length(P)
%     tic;
    clear D_pr
    sec = P{n};
    disp(strcat('-----',sec,'-----'))
    fname = lower(sec);
    com = ['load ' reported_file oceans{3} '/' sec '/' fname '.mat'];
    eval(com);
    years = eval(sec).years;
    ll_grid = eval([sec '.ll_grid']);
    pr_grid = eval([sec '.pr_grid']);
    for i = 1:length(years)      
        disp(['gridded - MAT, ', years{i}]);
        com = ['[s, m] = copyfile(''' easyocean_file sec '/' 'configuration_' years{i} '.m'', ''configuration.m'');'];
        eval(com);
        if s ~= 1, error(['copyfile ' num2str(i) ':', m]); end
        com = ['[D_pr(' num2str(i) ')' ,', pr_grid_lr] ', '= grid_data_pressure(D_reported(' num2str(i) '), ll_grid, pr_grid);'];
        eval(com);
    end
    if ~exist(strcat([gridded_file oceans{3} '/' sec '/']), 'dir')
        mkdir(strcat([gridded_file oceans{3} '/' sec '/']))
    end
    com = ['save ' gridded_file oceans{3} '/' sec '/' fname '.mat' ' D_pr ll_grid pr_grid pr_grid_lr'];
    eval (com);
%     toc;
end

% S
for n = 1:length(S)
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
        com = ['[D_pr(' num2str(i) ')' ,', pr_grid_lr] ', '= grid_data_pressure(D_reported(' num2str(i) '), ll_grid, pr_grid);'];
        eval(com);
    end
    if ~exist(strcat([gridded_file oceans{4} '/' sec '/']), 'dir')
        mkdir(strcat([gridded_file oceans{4} '/' sec '/']))
    end
    com = ['save ' gridded_file oceans{4} '/' sec '/' fname '.mat' ' D_pr ll_grid pr_grid pr_grid_lr'];
    eval (com);
%     toc;
end

%     figure
%     semilogx(N2, -p_);
%     hold on
%     semilogx(N2_filter1, -p_filter1, "linewidth", 1.5);
%     hold on
%     semilogx(N2_filter2, -p_filter2, "linewidth", 1.5);
%     hold on
%     semilogx(N2_filter3, -p_filter3, "linewidth", 1.5);

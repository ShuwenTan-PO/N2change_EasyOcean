function [D_pressure, pr_grid_lr] = grid_data_pressure(D_reported, ll_grid, pr_grid, fl, lr)
%
% Horizontal & vertical interpolation of `reported_data` D_reported.
% lat/lon grid `ll_grid` and pressure_grid `pr_grid` are optional
%

%%%
%%% User defined interpolation functions in configuration.m
%%%

%%% add new features: compute sigma4, and N2 from
% 0) raw SA and CT; 
% 1) filtered SA and CT; 
% 2) filtered sigma4; 
% 3) linear regression in sigma4 layers (layer thickness = fl/10 by default)
% fl is the cut-off distance for a hanning filter in meters, default = 400 m
% lr is the layer thickness for computing N2 using the linear regression method, default = fl/10 m
%%% modified by S.Tan 04/11/2021, IDEO
% 4) linear regression in sigma_center layers (layer thickness = fl/10 by default)
% all except N2_fl3 and N2_fl4 are on pr_grid, the two outliers are on pr_grid_lr
%%% modified by S.Tan 03/15/2022, MorningSide Heights
% fix bug: 
% 1. horizontal objective mapping use maximum depth as the cutoff z, should be maximum pressure
% 2. after horizontal interpolation, get rid of aritificially filled near bottom data
% 3. fix maxpinterp instead of maxpiterp in hinterp_bylat and hinterp_bylon
%%% modified by S.Tan 03/16/2022, MorningSide Heights

configuration;
%%%

% default fl = 400 m, ll = fl/10 m
if nargin<4;fl=400;end
if nargin<5;lr=fl/10;end
   
pr_grid_lr = [pr_grid(1):lr:pr_grid(end)]';

% lat/lon and depth
stations = D_reported.Station;
nstn = length(stations);
[lats, lons, deps, tims, pre_max012, pre_max34] = deal(NaN(1, nstn));
for i = 1:nstn
    lats(i) = stations{i}.Lat;
    lons(i) = stations{i}.Lon;
    deps(i) = stations{i}.Depth;
    tims(i) = stations{i}.Time;
end

[idx, ll, isAtlantic] = sort_stations(lons, lats);
ls = ll(idx); % sort
if isAtlantic
    lons = ll;
end

ctdprs = D_reported.CTDprs;
ctdtem = D_reported.CTDtem;
ctdsal = D_reported.CTDsal;
ctdoxy = D_reported.CTDoxy;
ctdCT = D_reported.CTDCT;
ctdSA = D_reported.CTDSA;

% vertical interpolation
[ctdtem_v, ctdsal_v, ctdoxy_v] = deal(NaN(length(pr_grid), nstn));
[ctdCT_v, ctdSA_v] = deal(NaN(length(pr_grid), nstn));
for i = 1:nstn
    ctdtem_v(:,i) = vinterp_handle(ctdprs(:,i), ctdtem(:,i), pr_grid);
    ctdsal_v(:,i) = vinterp_handle(ctdprs(:,i), ctdsal(:,i), pr_grid);
    ctdoxy_v(:,i) = vinterp_handle(ctdprs(:,i), ctdoxy(:,i), pr_grid);
    ctdCT_v(:,i) = vinterp_handle(ctdprs(:,i), ctdCT(:,i), pr_grid);
    ctdSA_v(:,i) = vinterp_handle(ctdprs(:,i), ctdSA(:,i), pr_grid);
end
% compute sig4, n2 from raw SA & CT, and n2 from smoothed data
% then vertical interpolation
ctdsig4 = gsw_sigma4(ctdSA, ctdCT);
[ctdsig4_v, ctdN2_v, ctdN2_filter1_v, ctdN2_filter2_v] = deal(NaN(length(pr_grid), nstn));
[ctdN2_filter3_v, ctdN2_filter4_v] = deal(NaN(length(pr_grid_lr), nstn));
for i = 1:nstn
    % N2: raw SA and CT; 
    [N2, p_] = gsw_Nsquared(ctdSA(:,i), ctdCT(:,i), ctdprs(:,i), lats(i));
    % N2_filter1: filtered SA and CT; 
    ig = find(isfinite(ctdprs(:,i)) & isfinite(ctdSA(:,i)) & isfinite(ctdCT(:,i)) & isfinite(ctdsig4(:,i)));
    p = ctdprs(ig,i);
    t = ctdtem(ig,i);
    sa = ctdSA(ig,i);
    ct = ctdCT(ig,i);
    sig4 = ctdsig4(ig,i);
    dp = mode(diff(p));
    nfilt = round(fl / dp);
    if length(ig) > nfilt
        sa_ = filtend(sa, nfilt); % default hanning filter
        ct_ = filtend(ct, nfilt);
        sig4_ = filtend(sig4, nfilt);
    else
        sa_ = NaN(size(sa));
        ct_ = NaN(size(ct));
        sig4_ = NaN(size(sig4));
    end
    [N2_filter1, p_filter1] = gsw_Nsquared(sa_, ct_, p, lats(i));
    % N2_filter2: filtered sigma4; 
    g = gsw_grav(lats(i), p);
    b = -g.*(sig4_-mean(sig4_))./(mean(sig4_)+1000);    
    dep = gsw_z_from_p(p, lats(i));
    N2_filter2 = gradient(b, dep); p_filter2 = p;
    % N2_filter3: linear regression in sigma4 layers
    if length(p)>0
        p_layer = [p(1):lr:p(end)]';
        p_filter3 = nan(length(p_layer)-1, 1);
        N2_filter3 = nan(length(p_layer)-1, 1);
        for j=1:length(p_filter3)        
            mask = find(p>=p_layer(j) & p<p_layer(j+1));
            if length(mask) >= 2
                p_filter3(j) = mean(p(mask));
                g = gsw_grav(lats(i), p_filter3(j));       
                b = -g.*(sig4(mask)-mean(sig4))./(mean(sig4)+1000);    
                dep = gsw_z_from_p(p(mask), lats(i));
                dummy = polyfit(dep, b, 1);
                N2_filter3(j) = dummy(1);
            end
        end
    else
        p_filter3 = [];
        N2_filter3 = [];
    end
    % N2_filter4: linear regression in sigma_center layers
    if length(p)>0
        p_layer = [p(1):lr:p(end)]';
        p_filter4 = nan(length(p_layer)-1, 1);
        N2_filter4 = nan(length(p_layer)-1, 1);
        for j=1:length(p_filter4)        
            mask = find(p>=p_layer(j) & p<p_layer(j+1));
            if length(mask) >= 2
                p_filter4(j) = mean(p(mask));
                g = gsw_grav(lats(i), p_filter4(j));  
                sig = gsw_pot_rho_t_exact(sa(mask), t(mask), p(mask), p_filter4(j));
                b = -g.*(sig-mean(sig))./(mean(sig));    
                dep = gsw_z_from_p(p(mask), lats(i));
                dummy = polyfit(dep, b, 1);
                N2_filter4(j) = dummy(1);
            end 
        end
    else
        p_filter4 = [];
        N2_filter4 = [];
    end
    
    ctdsig4_v(:,i) = vinterp_handle(ctdprs(:,i), ctdsig4(:,i), pr_grid);
    ctdN2_v(:,i) = vinterp_handle(p_, N2, pr_grid);
    if length(p_filter1) < 2
        ctdN2_filter1_v(:,i) = vinterp(p_filter1, N2_filter1, pr_grid);
    else
        ctdN2_filter1_v(:,i) = vinterp_handle(p_filter1, N2_filter1, pr_grid);
    end
    if length(p_filter2) < 2
        ctdN2_filter2_v(:,i) = vinterp(p_filter2, N2_filter2, pr_grid);
    else
        ctdN2_filter2_v(:,i) = vinterp_handle(p_filter2, N2_filter2, pr_grid);
    end
    if length(p_filter3) > 2
        mask = find(isnan(p_filter3));
        if length(mask) > 0
            ctdN2_filter3_v(:,i) = interp1(p_filter3(~isnan(p_filter3)), N2_filter3(~isnan(p_filter3)), pr_grid_lr, 'linear');
            if mask(1) > 1 & mask(end) < length(p_filter3)
                mask_ = find(pr_grid_lr > p_filter3(mask(1)-1) & pr_grid_lr < p_filter3(mask(end)+1));
            elseif mask(1) == 1 & mask(end) < length(p_filter3)
                mask_ = find(pr_grid_lr > p_filter3(mask(1)) & pr_grid_lr < p_filter3(mask(end)+1));
            else
                mask_ = find(pr_grid_lr > p_filter3(mask(1)-1) & pr_grid_lr < p_filter3(mask(end)));
            end
            ctdN2_filter3_v(mask_,i) = nan;
        else
            ctdN2_filter3_v(:,i) = interp1(p_filter3, N2_filter3, pr_grid_lr, 'linear');
        end 
            
    end
    if length(p_filter4) > 2
        mask = find(isnan(p_filter4));
        if length(mask) > 0
            ctdN2_filter4_v(:,i) = interp1(p_filter4(~isnan(p_filter4)), N2_filter4(~isnan(p_filter4)), pr_grid_lr, 'linear'); 
            if mask(1) > 1 & mask(end) < length(p_filter4)
                mask_ = find(pr_grid_lr > p_filter4(mask(1)-1) & pr_grid_lr < p_filter4(mask(end)+1));
            elseif mask(1) == 1 & mask(end) < length(p_filter4)
                mask_ = find(pr_grid_lr > p_filter4(mask(1)) & pr_grid_lr < p_filter4(mask(end)+1));
            else
                mask_ = find(pr_grid_lr > p_filter4(mask(1)-1) & pr_grid_lr < p_filter4(mask(end)));
            end
            ctdN2_filter4_v(mask_,i) = nan;
        else 
            ctdN2_filter4_v(:,i) = interp1(p_filter4, N2_filter4, pr_grid_lr, 'linear'); 
        end 
    end
end

for i = 1:nstn
    a = pr_grid(find(~isnan(ctdN2_v(:,i))));
    if length(a)>0
        pre_max012(i) = a(end);
    else
        pre_max012(i) = deps(i);
    end
    a = pr_grid_lr(find(~isnan(ctdN2_filter4_v(:,i))));
    if length(a)>0
        pre_max34(i) = a(end);
    else
        pre_max34(i) = deps(i);
    end
end

% Interpolate chunk by chunk -- do not interpolate if more than MAX_SEPARATION deg apart
% (defined in configuration)
chunk = {};
len = 0;
idxhere = [1];
for i = 2:length(ls)
    if abs(ls(i) - ls(i-1))  < MAX_SEPARATION;
        idxhere = [idxhere, i];
        continue;
    end
    len = len + 1;
    chunk(1, len) = {idxhere};
    idxhere = [i];
end
chunk(1, len+1) = {idxhere};

[ctdtem_hv, ctdsal_hv, ctdoxy_hv] = deal(NaN(length(pr_grid), length(ll_grid)));
[ctdCT_hv, ctdSA_hv] = deal(NaN(length(pr_grid), length(ll_grid)));
[ctdsig4_hv, ctdN2_hv, ctdN2_filter1_hv, ctdN2_filter2_hv] = deal(NaN(length(pr_grid), length(ll_grid)));
[ctdN2_filter3_hv, ctdN2_filter4_hv] = deal(NaN(length(pr_grid_lr), length(ll_grid)));
maxp_ll_grid = deal(NaN(1, length(ll_grid)));
for i = 1:length(chunk)
    idxhere = chunk{i};
    lshere = ls(idxhere);
    ig = find(min(lshere) < ll_grid & ll_grid < max(lshere));
    if length(idxhere) < 2 || isempty(ig)
        continue;
    end
    if length(ig) > 1
        ctdtem_hv(:,ig) = hinterp_handle(ctdtem_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                         pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdsal_hv(:,ig) = hinterp_handle(ctdsal_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                         pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdoxy_hv(:,ig) = hinterp_handle(ctdoxy_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                         pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdCT_hv(:,ig) = hinterp_handle(ctdCT_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                        pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdSA_hv(:,ig) = hinterp_handle(ctdSA_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                        pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdsig4_hv(:,ig) = hinterp_handle(ctdsig4_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                        pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdN2_hv(:,ig) = hinterp_handle(ctdN2_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                        pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdN2_filter1_hv(:,ig) = hinterp_handle(ctdN2_filter1_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                                pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdN2_filter2_hv(:,ig) = hinterp_handle(ctdN2_filter2_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                                pre_max012(idxhere), pr_grid, ll_grid(ig));
        ctdN2_filter3_hv(:,ig) = hinterp_handle(ctdN2_filter3_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                                pre_max34(idxhere), pr_grid_lr, ll_grid(ig));
        [ctdN2_filter4_hv(:,ig), maxp_ll_grid(ig)] = hinterp_handle(ctdN2_filter4_v(:,idxhere), lons(idxhere), lats(idxhere), ...
                                                                    pre_max34(idxhere), pr_grid_lr, ll_grid(ig));
    end 
end

% Use nearest time to the grid point
ntime = interp1(ll, tims(idx), ll_grid, 'nearest');

D_pressure = struct('Station', {stations}, ...
                    'NTime', ntime, ...
                    'CTDtem', ctdtem_hv, ...
                    'CTDsal', ctdsal_hv, ...
                    'CTDoxy', ctdoxy_hv, ...
                    'CTDCT', ctdCT_hv, ...
                    'CTDSA', ctdSA_hv, ...
                    'CTDsig4', ctdsig4_hv, ...
                    'CTDN2', ctdN2_hv, ...
                    'CTDN2_f1', ctdN2_filter1_hv, ...
                    'CTDN2_f2', ctdN2_filter2_hv, ...
                    'CTDN2_f3', ctdN2_filter3_hv, ...
                    'CTDN2_f4', ctdN2_filter4_hv);
end
% 
% figure
% plot(ll_grid, -maxp_ll_grid, 'k');
% hold on
% plot(lons, -pre_max34, 'r')
% for i = 1:length(ll_grid)
%     a = -pr_grid_lr(find(~isnan(ctdN2_filter3_hv(:,i))));
%     if length(a)>0
%         hold on
%         plot(ll_grid(i), a(end), 'r.')
%     end
% end
% for i = 1:length(lons)
%     a = -pr_grid_lr(find(~isnan(ctdN2_filter3_v(:,i))));
%     if length(a)>0
%         hold on
%         plot(lons(i), a(end), 'b.')
%     end
% end
% 
% figure
% plot(ll_grid, -maxp_ll_grid, 'k');
% hold on
% plot(lons, -pre_max012, 'r')
% for i = 1:length(ll_grid)
%     a = -pr_grid(find(~isnan(ctdN2_filter2_hv(:,i))));
%     if length(a)>0
%         hold on
%         plot(ll_grid(i), a(end), 'r.')
%     end
% end
% for i = 1:length(lons)
%     a = -pr_grid(find(~isnan(ctdN2_filter2_v(:,i))));
%     if length(a)>0
%         hold on
%         plot(lons(i), a(end), 'b.')
%     end
% end
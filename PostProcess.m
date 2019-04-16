%% Generate URFs for the 7 refinements
% add the path to the NPSAT toolbox matlab commands
%addpath('/home/giorgk/CODES/NPSAT/matlab')
%
npsat_data = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/';
res_dir = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/Results/';
%%
for iref = 8:8
    WellURF = [];
    urf_files = dir([res_dir 'Ref' num2str(iref) 'TR/*.urfs']);
    for jj = 1:length(urf_files)
        display(['Refinement: ' num2str(iref)])
        if urf_files(jj,1).bytes > 0
            temp = readURFs([res_dir 'Ref' num2str(iref) 'TR/' urf_files(jj,1).name], []);
            WellURF = [WellURF;temp];
        end
    end
    %save([res_dir 'URFs_Ref_an' num2str(iref)], 'WellURF');
end
%% Import and process concentration input
% import the Conc540transientmatrix.dat file as matrix and arrange the rows
% rearrange the data
Conc = nan(153,137,540);
for ii = 1:size(Conc,3)
   Conc(:,:,ii) = Conc540transientmatrix(1+(ii-1)*153:ii*153,:);
end
%% version 2
% import the ConcbasedonyearlyRCH481.dat file as matrix and arrange the rows
% rearrange the data
Conc = nan(153,137,540);
for ii = 1:size(Conc,3)
   Conc(:,:,ii) = ConcbasedonyearlyRCH481(1+(ii-1)*153:ii*153,:);
end
%% plot concentration for any cell
figure(3);
hold on
for ii = 1:153
    for jj = 1:137
        plot(reshape(Conc(ii,jj,:),540,1))
    end
end
%%
Input_LF = zeros(153*137,45);
for ii = 1:153
    ii
    for jj = 1:137
        Input_LF(sub2ind([153 137],ii,jj),:) = reshape(Nload(ii,jj,:),45,1)';
    end
end
%%
prc = prctile(Input_LF,[10:10:90]);
%%
figure(20);plot(reshape(Conc(110,80,:),45,1))
%% define coordinate system
xcoord = @(ii)200 + (ii-1).*400;
ycoord = @(ii)61200 - (200 + (ii-1).*400);
%}
%% Load and preprocess Nloading scenario
load([npsat_data 'NloadingScenarios/300_N_concentration_Load']);
B = zeros(153,137,size(A,1)/153);
for ii = 1:size(A,1)/153
    B(:,:,ii) = A((ii-1)*153+1:ii*153,:);
end
%% Average loading every 12 months
Conc = zeros(153,137, size(B,3)/12);
for ii = 1:size(Conc,3)
    Conc(:,:,ii) = mean(B(:,:,(ii-1)*12+1:ii*12),3);
end
%% load and process the Nloading reduce Scenario
load([npsat_data 'NloadingScenarios/reduced_NLoadConc_input_45yrs'])
B = zeros(153,137,size(A,1)/153);
for ii = 1:size(A,1)/153
    B(:,:,ii) = A((ii-1)*153+1:ii*153,:);
end
%% Average loading every 12 months
Conc = zeros(153,137, size(B,3)/12);
for ii = 1:size(Conc,3)
    Conc(:,:,ii) = mean(B(:,:,(ii-1)*12+1:ii*12),3);
end
%% Convolute concentration with URFs
wells = readWells('input/well_data4.npsat');
exclID = find(wells(:,2)<=13000) - 1;

Nyears = 45;

for iref = 8:8
    iref
    load([res_dir 'URFs_Ref' num2str(iref) 'TR']);
%     dlt = [];
%     for jj = 1:length(WellURF)
%         if ~isempty(find(exclID == WellURF(jj,1).Eid, 1))
%             dlt = [dlt; jj];
%         end
%     end
%     WellURF(dlt,:) = [];
    Eid = [WellURF.Eid]';
    Sid = [WellURF.Sid]';
    v_cds = [WellURF.v_cds]';
    URFS = zeros(length(WellURF), 200);
    LF = zeros(length(WellURF), Nyears);
    % build loading functions
    for jj = 1:length(WellURF)
        p = WellURF(jj,1).p_lnd;
        [I, J] = findIJ_Modesto(p(1),p(2));
        lf = reshape(Conc(I,J,:), size(Conc,3), 1)';
        LF(jj,:) = lf(1:Nyears);
        URFS(jj,:) = WellURF(jj,1).URF;
    end
    LF(isinf(LF)) = 0;
    LF(isnan(LF)) = 0;
    % do the convolution
    BTC = ConvoluteURF(URFS, LF, 'cpp');
    % calculate average concentration per well
    wellids = unique(Eid);
    for jj = 1:length(wellids)
        ids = find(Eid == wellids(jj,1));
        w = v_cds(ids)./sum(v_cds(ids));
        btc = bsxfun(@times, BTC(ids,:), w);
        WellBTC(jj,1).Eid = wellids(jj,1);
        WellBTC(jj,1).BTC = sum(btc, 1);
    end
    AllRefRes(iref+1,1).ref = iref;
    AllRefRes(iref+1,1).btc = WellBTC;
end
%% Plot loading functions
LF_prc = prctile(LF, [10:10:90]);
figure(2);
plot(LF_prc')
%% Group all BTC into one variable
for iref = 9:9
    BTCALL = [];
    for jj = 1:length(AllRefRes(iref,1).btc)
        if any(isinf(AllRefRes(iref,1).btc(jj,1).BTC))
            continue
        end
        
        BTCALL = [BTCALL; AllRefRes(iref,1).btc(jj,1).BTC];
    end
    AllRefRes(iref,1).BTCALL = BTCALL;
end
%% 
save('WellBTCResults_v1','AllRefRes');
%% 
clf
prc = 50;
figure(3);hold on
for ii = 4:4
    plot(prctile(AllRefRes(ii,1).BTCALL,prc),'linewidth', 1.5)
    hold on
end
%% ovelay mt3d results
load('MT3D_Results.mat')
mt3d_disp = [];
for ii = 1:length(MT3DRes)
    mt3d_disp = [mt3d_disp; prctile(MT3DRes(ii,1).BTC, prc)];
end
plot(mt3d_disp');
%%
btc_prc = prctile(AllRefRes(ii,1).BTCALL, [10:10:90]);
figure(2);
plot(btc_prc')
%% plot the concentration around the wells
% copy paste to an empty variable wellIJ the row and columns from the well.xlsx
wellIJ = unique(wellIJ,'rows');
expnd = [-1 -1;-1 0;-1 1;0 -1; 0 1;-1 1;0 1;1 1];
moreWellCells = [];
for ii = 1:size(wellIJ,1)
    moreWellCells = [moreWellCells;wellIJ(ii,:)];
    for jj = 1:size(expnd,1)
        newij = wellIJ(ii,:) + expnd(jj,:);
        if newij(1) < 1 || newij(1) > 153
            continue;
        end
        if newij(2) < 1 || newij(2) > 137
            continue;
        end
        moreWellCells = [moreWellCells;newij];
    end
end
moreWellCells = unique(moreWellCells,'rows');
%% plot the concentration for the selected cells only
figure(10); hold on
temp = [];
for ii = 1:size(moreWellCells,1)
    temp =[temp;reshape(Conc(moreWellCells(ii,1),moreWellCells(ii,2),:),45,1)'];
end
temp_prc = prctile(temp, [10:10:90]);
%%
%plot(temp','linewidth',0.5, 'Color',[0.5 0.5 0.5])
%hold on
plot(temp_prc','linewidth', 2)
%}
%% Run the 7th refinement 
%{
iref = 7;
alphas = [0.1  0.15  0.2 0.25 0.3 0.4 0.6 1 1.3 1.6 1.9];
betas = 0.5:0.1:1.2;
for ia = 1:length(alphas)
    for ib = 1:length(betas)
        opt.alpha = alphas(ia);
        opt.beta = betas(ib);
        WellURF = [];
        urf_files = dir(['Results/Ref' num2str(iref) '/*.urfs']);
        for jj = 1:length(urf_files)
            if urf_files(jj,1).bytes > 0
                temp = readURFs(['Results/Ref' num2str(iref) '/' urf_files(jj,1).name], opt);
                WellURF = [WellURF;temp];
            end
        end
        save(['Results/URFs_Ref_a' num2str(ia) '_b' num2str(ib)], 'WellURF');
    end
end
%}
%% This is an attempt to make the above loop to run under Octave
% pkg load statistics %(for the nanmean function)
res_path = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/Results/';
iref = 7;
alphas = [0.1  0.15  0.2 0.25 0.3 0.4 0.6 1 1.3 1.6 1.9];
betas = 0.5:0.1:1.2;
% define istart iend on the workspace
for ii = istart:iend
    [ia, ib] = ind2sub([length(alphas) length(betas)], ii); 
    opt.alpha = alphas(ia);
    opt.beta = betas(ib);
    WellURF = [];
    urf_files = dir([res_path 'Ref' num2str(iref) '/*.urfs']);
    for jj = 1:length(urf_files)
        if urf_files(jj,1).bytes > 0
            temp = readURFs([urf_files(jj,1).folder '/' urf_files(jj,1).name], opt);
            WellURF = [WellURF;temp];
        end
    end
    save([res_path 'URFs_Ref_a' num2str(ia) '_b' num2str(ib)], 'WellURF');
end
%% 
URFS = nan(length(WellURF),200);
for ii = 1:length(WellURF)
    URFS(ii,:) = WellURF(ii).URF;
end
%%
plot(URFS(22001:44000,:)')
%% Distribution of lengths
clear hh
for iref = 0:8
    iref
    load([res_dir 'URFs_Ref_' num2str(iref)]);
    dlt = [];
    for jj = 1:length(WellURF)
        if ~isempty(find(exclID == WellURF(jj,1).Eid, 1))
            dlt = [dlt; jj];
        end
    end
    WellURF(dlt,:) = [];
    
    [xx, ff] = ecdf([WellURF.L]');
    hh(iref+1,1).xx = xx;
    hh(iref+1,1).ff = ff;
end
%%
for ii = 1:length(hh)
    semilogx(hh(ii,1).ff, hh(ii,1).xx)
    hold on
end


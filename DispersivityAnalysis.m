%% Convolute the dispersivity results with the loading functions
res_path = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/Results/SensitivityAnalysisData/';
load('N_mass_load.mat', 'Nload_conc')
Conc = Nload_conc;


for ii = 1:144
    ii
    load([res_path 'URFs_Ref_scenID_' num2str(ii) '.mat']);
    Eid = [WellURF.Eid]';
    Sid = [WellURF.Sid]';
    v_cds = [WellURF.v_cds]';
    URFS = zeros(length(WellURF), 200);
    LF = zeros(length(WellURF),45);
    % build loading functions
    for jj = 1:length(WellURF)
        p = WellURF(jj,1).p_lnd;
        [I, J] = findIJ_Modesto(p(1),p(2));
        lf = reshape(Conc(I,J,:),45,1)';
        LF(jj,:) = lf;
        URFS(jj,:) = WellURF(jj,1).URF;
    end
    URFS(URFS<0) = 0;
    LF(isinf(LF)) = 0;
    LF(isnan(LF)) = 0;
    BTC = ConvoluteURF(URFS, LF, 'cpp');
    wellids = unique(Eid);
    WellBTC = [];
    for jj = 1:length(wellids)
        ids = find(Eid == wellids(jj,1));
        w = v_cds(ids)./sum(v_cds(ids));
        btc = bsxfun(@times, BTC(ids,:), w);
        WellBTC = [WellBTC; sum(btc, 1)];
    end
    DispSensBTC(ii,1).wellids = wellids;
    DispSensBTC(ii,1).BTC = WellBTC;
end
%%
plot(prctile(DispSensBTC(92,1).BTC, [10:10:90])')
%%
meanbtc = [];
for ii = 1:144
    meanbtc = [meanbtc; prctile(DispSensBTC(ii,1).BTC, 50)];
end
plot(meanbtc')
%% Compare NPSAT with MT3D
load('MT3D_Results.mat')
load('DispSensBTC.mat')
%% 
prc = 10;
npsat_disp = [];
mt3d_disp = [];
for ii = 1:length(DispSensBTC)
    npsat_disp = [npsat_disp; prctile(DispSensBTC(ii,1).BTC, prc)];
end

for ii = 1:length(MT3DRes)
    mt3d_disp = [mt3d_disp; prctile(MT3DRes(ii,1).BTC, prc)];
end
clf
plot(npsat_disp', 'Color', [0.5, 0.5 0.5], 'linewidth', 0.5);
hold on
plot(mt3d_disp', 'linewidth', 1.5);

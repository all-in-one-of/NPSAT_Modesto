NPSAT_input_files_BAU = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/NPSAT_input_files_BAU/';
team_drive = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/';
%% Load results from MT3D results - different alpha
cnt = 1;
clear MT3DRes
for ii = [1 5 11 50 150]
    tbl = xlsread([NPSAT_input_files_BAU 'Results_MT3D_disp11.xlsx'], ['Disp ' num2str(ii)]);
    tbl(tbl(:,1) >= 121,:) = [];
    MT3DRes(cnt,1).disp = ii;
    MT3DRes(cnt,1).BTC = tbl(2:end,4:2:92);
    MT3DRes(cnt,:).RC = tbl(2:end,1:2);
    cnt = cnt + 1;
end
%% Load newest results
tbl = xlsread([NPSAT_input_files_BAU 'Results_MT3D_disp11.xlsx'], 'Sheet1');
%tbl(tbl(:,1) >= 121,:) = [];
MT3DRes.disp = 11;
MT3DRes.BTC = tbl(2:end,4:2:92);
MT3DRes.RC = tbl(2:end,1:2);
%% load 300 yr BCT from Mt3D
load([team_drive 'Results/BTC_Steady_yr2000_300yrs_MT3D'])
BTC_2000SS = BTC;
%% Compare 200 year predictions between Best of NPSAT and MT3D
load([team_drive 'Results/BTC_Steady_yr2000_300yrs_MT3D'])
BTC_2000SS = BTC;
load('WellBTCResults_v1.mat')
BTC_NPSAT = AllRefRes(9,1).BTCALL;

%load('MT3D_Results.mat')
%load('WellBTCResults_v1.mat')
load('WellData_Init.mat')
%wells = readWells('input/well_data4.npsat');
xcoord = @(ii)200 + (ii-1).*400;
ycoord = @(ii)61200 - (200 + (ii-1).*400);
%% Create the correspondance between NPSAT and MT3D well ids and 
for ii = 1:size(WellData,1)
    id = find(BTC_2000SS(:,1) == WellData(ii,1).I & BTC_2000SS(:,2) == WellData(ii,1).J);
    WellData(ii,1).mt3dID = id;
end
%% Exclude Wells
% Here we will exclude the wells that their source area involves the
% Constant head boundary. These are the wells with I >= 121
incID = find([WellData.I]' < 121 & [WellData.I]' > 0);
%% preproc for plots
colormat = [ ...
    166,206,227; ...
    31,120,180; ...
    178,223,138; ...
    51,160,44; ...
    251,154,153; ...
    227,26,28; ...
    253,191,111; ...
    255,127,0; ...
    202,178,214; ...
]/255;
colormat1 = [...
    228,26,28; ...
    55,126,184; ...
    77,175,74; ...
    152,78,163; ...
    255,127,0; ...
    255,255,51; ...
    166,86,40; ...
    247,129,191; ...
    153,153,153; ...
]/255;
colormat2 = [ ...
    255,247,236; ...
    254,232,200; ...
    253,212,158; ...
    253,187,132; ...
    252,141,89; ...
    239,101,72; ...
    215,48,31; ...
    179,0,0; ...
    127,0,0; ...
]/255;
colormat3 = [ ...
    102,194,165; ...
    252,141,98; ...
    141,160,203; ...
]/255;
%% plot all refinements vs Mt3DMS dispersivity 11
figure(10);clf
figure(10); hold on
prc = 75;
temp = [];
leg_text = [];
for ii = 1:9
    temp = [temp; prctile(AllRefRes(ii,1).BTCALL(incID,:),prc)];
    leg_text{ii} = ['# Ref: ' num2str(ii-1)];
end
h = plot(temp', '-','linewidth', 2);
set(h, {'color'}, num2cell(colormat,2));
plot(prctile(MT3DRes.BTC(incID,:), prc)', 'k','linewidth', 3);
leg_text{10} = 'MT3D';
xlabel('Time [years]');
ylabel('Concentration [M/L^3]')
lgd = legend(leg_text, 'Location','northwest');
lgd.NumColumns = 2;
grid on
%%
print([team_drive '/Figures/AllRef_Disp11_50perc'],'-dpng','-r600')
%% Compare percentiles between ref 8 and mt3dms
prc =[10:10:90];% [5 25 50 75 95];%[10:10:90];
yl = 14; %this is the y max
figure(10);clf
figure(10); hold on
npsat_prc = prctile(AllRefRes(9,1).BTCALL(incID,:), prc);
mt3d_prc = prctile(BTC_2000SS(incID,3:202), prc);
for ii = 1:length(prc)
    plot(npsat_prc(ii,:),'--', 'color', colormat(ii,:), 'linewidth', 3);
    plot(mt3d_prc(ii,:), 'color', colormat(ii,:), 'linewidth', 3);
    text(201,mt3d_prc(ii,end),['perc_{' num2str(prc(ii)) '}'],'fontsize', 14);
end
xlabel('Time [years]','fontsize', 14);
ylabel('Concentration [M/L^3]','fontsize', 14)
text(5, yl -1.5, 'Solid lines:  MT3D','fontsize', 16)
text(5, yl -0.5, 'Dashed lines: NPSAT (# Ref 8)','fontsize', 16)
ylim([0 yl]);
grid on
set(gca,'FontSize',12)
set(gca,'ytick',[0:2:yl])
%%
print([team_drive '/Figures/Ref8_Disp11_perc_v1'],'-dpng','-r600')
%% compare selected percentiles between NPSAT, steady, transient
% load the transiet simulation

load([team_drive 'Results/BTC_Transient_MT3D'])
BTC_TR = New_BTC;
%% plot the comparison
prc =[50];
figure(10);clf
figure(10); hold on
npsat_prc = prctile(AllRefRes(9,1).BTCALL(incID,1:45), prc);
mt3d_SSprc = prctile(BTC_2000SS(incID,3:47), prc);
mt3d_TRprc = prctile(BTC_TR(incID,3:47), prc);
plot(npsat_prc, 'color', colormat3(1,:), 'linewidth', 3, 'DisplayName', 'NPSAT');
plot(mt3d_SSprc, 'color', colormat3(2,:), 'linewidth', 3, 'DisplayName', 'MT3D (SS)');
plot(mt3d_TRprc, 'color', colormat3(3,:), 'linewidth', 3, 'DisplayName', 'MT3D (TR)');
xlabel('Time [years]','fontsize', 14);
ylabel('Concentration [M/L^3]','fontsize', 14)
title([num2str(prc) ' percentile'])
grid on
legend('Location', 'northwest')
set(gca,'FontSize',12)
%%
print([team_drive '/Figures/NPSAT_Mt3D_SS_TR_mean'],'-dpng','-r600')
%%
cnt = 0;
clf
yl = [0.5 3 9 13];
for ii = [5 25 75 95]
    cnt = cnt + 1;
    subplot(2,2,cnt); hold on
    npsat_prc = prctile(AllRefRes(9,1).BTCALL(incID,1:45), ii);
    mt3d_SSprc = prctile(BTC_2000SS(incID,3:47), ii);
    mt3d_TRprc = prctile(BTC_TR(incID,3:47), ii);
    plot(npsat_prc, 'color', colormat3(1,:), 'linewidth', 3, 'DisplayName', 'NPSAT');
    plot(mt3d_SSprc, 'color', colormat3(2,:), 'linewidth', 3, 'DisplayName', 'MT3D (SS)');
    plot(mt3d_TRprc, 'color', colormat3(3,:), 'linewidth', 3, 'DisplayName', 'MT3D (TR)');
    if cnt > 2
        xlabel('Time [years]','fontsize', 14);
    end
    if cnt == 1 || cnt == 3
        ylabel('Concentration [M/L^3]','fontsize', 14)
    end
    
    title([num2str(ii) ' percentile'])
    grid on
    legend('Location', 'northwest')
    set(gca,'FontSize',12)
    ylim([0 yl(cnt)]);
    subplot(2,2,cnt); hold off
end
%%
print([team_drive '/Figures/NPSAT_Mt3D_SS_TR_perc'],'-dpng','-r600')
%% Compare the differential change
load([team_drive 'Results/BTC_steady_reducedN_yr2000_45yrs'])
BTC_red_SS = BTC;
load([team_drive 'Results/BTC_Steady_yr2000_45years_MT3D'])
BTC_SS = BTC;
load([team_drive 'Results/BTC_Transient_MT3D'])
BTC_TR = New_BTC;
load([team_drive 'Results/BTC_transient_reducedN_45yrs'])
BTC_red_TR = New_BTC;
load('WellBTCResults_v1.mat')
npsat_btc = AllRefRes(9,1).BTCALL(incID,1:45);
load('WellBTCResults_v1_reduced.mat')
npsat_red_btc = AllRefRes(9,1).BTCALL(incID,1:45);
%%
prc = [50];
prc_btc_red_SS = prctile(BTC_red_SS(incID,3:end), prc);
prc_btc_SS = prctile(BTC_SS(incID,3:end), prc);
prc_btc_red_TR = prctile(BTC_red_TR(incID,3:end), prc);
prc_btc_TR = prctile(BTC_TR(incID,3:end), prc);
prc_btc_red_npsat = prctile(npsat_btc(incID,1:end), prc);
prc_btc_npsat = prctile(npsat_red_btc(incID,1:end), prc);
clf
hold on
plot(100*((prc_btc_red_SS - prc_btc_SS)./prc_btc_SS));
plot(100*((prc_btc_red_TR - prc_btc_TR)./prc_btc_SS));
plot(100*((prc_btc_red_npsat - prc_btc_npsat)./prc_btc_npsat));


%%
clf
hold on
plot(prctile(BTC_TR(incID, 3:end), [5 25 50 75 95])','b');
plot(prctile(BTC_red_TR(incID, 3:end), [5 25 50 75 95])','r');



%% plot fit vs screen length
screen_length_range = [WellData.SL]';
inc_test = [WellData.I]' < 121 & [WellData.I]' > 0;
screenDiscretization = linspace(min(screen_length_range), max(screen_length_range),100);
figure(10);clf
figure(10); hold on
iclr = 0;
ldg_txt = [];
for prc = 10:10:90
    iclr = iclr + 1;
%prc = 90;%[10:10:90];
cnt = 1;
clear MeanErr
for ii = 2:length(screenDiscretization)
    ids = find(inc_test & screen_length_range > screenDiscretization(ii-1) & screen_length_range < screenDiscretization(ii));
    if length(ids) < 30
        continue;
    end
    npsat_prc = prctile(AllRefRes(9,1).BTCALL(ids,:), prc);
    mt3d_prc = prctile(MT3DRes.BTC(ids,:), prc);
    R = corrcoef(npsat_prc, mt3d_prc);
    MeanErr(cnt,:) = [screenDiscretization(ii-1) screenDiscretization(ii) mean(mean(sqrt((npsat_prc - mt3d_prc).^2))) R(1,2)];
    cnt = cnt + 1;
end

stairs(mean(MeanErr(:,1:2),2), MeanErr(:,3), 'color', colormat2(iclr,:)  , 'linewidth', 3);
%for ii = 1:size(MeanErr)
%    plot(MeanErr(ii,1:2), MeanErr(ii,4)*ones(1,2), 'color', colormat(iclr,:)  , 'linewidth', 2)
ldg_txt{iclr} = ['p_{' num2str(prc) '}'];
%end
end
lgd = legend(ldg_txt, 'Location','southwest');
whitebg([.8 .8 .8])
grid on
xlabel('Screen length [m]')
ylabel('Correlation coefficient')
ylabel('Square error [M/L^3]')
title('NPSAT vs MT3D')
%%
print([team_drive '/Figures/SL_CorrCoef'],'-dpng','-r600')
%%
print([team_drive '/Figures/SL_SqErr'],'-dpng','-r600')
%% plot fit vs Depth
depth_range = [WellData.D_wt_top_mod]';
inc_test = [WellData.I]' < 121 & [WellData.I]' > 0;
DepthDiscr = linspace(min(depth_range), max(depth_range),100);
figure(10);clf
figure(10); hold on
iclr = 0;
ldg_txt = [];
for prc = 10:10:90
    iclr = iclr + 1;
%prc = 90;%[10:10:90];
cnt = 1;
clear MeanErr
for ii = 2:length(DepthDiscr)
    ids = find(inc_test & depth_range > DepthDiscr(ii-1) & screen_length_range < DepthDiscr(ii));
    if length(ids) < 30
        continue;
    end
    npsat_prc = prctile(AllRefRes(9,1).BTCALL(ids,:), prc);
    mt3d_prc = prctile(MT3DRes.BTC(ids,:), prc);
    R = corrcoef(npsat_prc, mt3d_prc);
    MeanErr(cnt,:) = [DepthDiscr(ii-1) DepthDiscr(ii) mean(mean(sqrt((npsat_prc - mt3d_prc).^2))) R(1,2)];
    cnt = cnt + 1;
end

stairs(mean(MeanErr(:,1:2),2), MeanErr(:,3), 'color', colormat2(iclr,:)  , 'linewidth', 3);
%for ii = 1:size(MeanErr)
%    plot(MeanErr(ii,1:2), MeanErr(ii,4)*ones(1,2), 'color', colormat(iclr,:)  , 'linewidth', 2)
ldg_txt{iclr} = ['p_{' num2str(prc) '}'];
%end
end
%lgd = legend(ldg_txt, 'Location','southwest');
lgd = legend(ldg_txt, 'Location','northwest');
whitebg([.8 .8 .8])
grid on
xlabel('Depth (water table to top screen) [m]')
%ylabel('Correlation coefficient')
ylabel('Square error [M/L^3]')
title('NPSAT vs MT3D')
%%
print([team_drive '/Figures/Depth_CorrCoef'],'-dpng','-r600')
%%
print([team_drive '/Figures/Depth_SqErr'],'-dpng','- r600')
function RunDispAnalysis(scenarioID)

pkg load statistics

scenarioID = scenarioID + 1;
iref = 7;

res_path = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/Results/';
% Scenario list
alphas = [0.1 0.125 0.15 0.175 0.2 0.25 0.3 0.4 0.5:0.2:2]';
betas = [0.4:0.1:1.2]';

[ia, ib] = ind2sub([length(alphas) length(betas)], scenarioID);
opt.alpha = alphas(ia);
opt.beta = betas(ib);

urf_files = dir([res_path 'Ref' num2str(iref) '/*.urfs']);
WellURF = [];
for jj = 1:length(urf_files)
    if urf_files(jj,1).bytes > 0
        temp = readURFs([urf_files(jj,1).folder '/' urf_files(jj,1).name], opt);
        WellURF = [WellURF;temp];
    end
end

save('-mat', [res_path 'URFs_Ref_scenID_' num2str(scenarioID) '.mat'], 'WellURF');


end


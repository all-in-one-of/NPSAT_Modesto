%% paths
NPSAT_input_files_BAU = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/NPSAT_input_files_BAU/';

%% Read input N loading
% The units are kg/hectare/yr
frmt = '%f';%[^\n\r]
for ii = 2:137
    frmt = [frmt '%f'];
end
frmt = [frmt '%[^\n\r]'];
Nload = nan(153, 137, 45);
for ii = 1:45
   fid = fopen([NPSAT_input_files_BAU 'N_mass_load/Nitrate_Mass_Load_' num2str(ii) '.dat'],'r');
   array = textscan(fid, frmt, 'Delimiter', ' ','MultipleDelimsAsOne', true);
   Nload(:,:,ii) = [array{1:end-1}];
   fclose(fid);
end
%% save 
save('N_mass_load', 'Nload');
%% plot Mass 
Nplot = nan(153*137,45);
cnt = 1;
for ii = 1:153
    for jj = 1:137
        nl = reshape(Nload(ii,jj,:),1,45);
        if any(nl > 0)
            Nplot(cnt,:) = nl;
            cnt = cnt + 1;
        end
    end
end
Nplot(cnt:end,:) = [];
%%
figure(4);plot(Nplot')
%% load Groundwater recharge for the steady steady
% The units are m^3/day
CBC = readModflowFlowdata([NPSAT_input_files_BAU 'MF_BAU_scheme9.crc']);
RCH_cbc = CBC(1,1).data;
CBC = readModflowFlowdata([NPSAT_input_files_BAU 'MF_BAU_scheme9.cbe']);
ET_cbc = CBC(1,1).data;
rch = sum(RCH_cbc,3);
et = sum(ET_cbc,3);
rch = rch - et;
%% Compute concentration
Nload = Nload./10000; %Kg/ha / year -> Kg/m^2 / year
rch = rch * 365; % m^3/ day -> m^3 /year
Nload = Nload .* (400*400); %Kg/m^2 -> Kg /year for each element
for ii = 1:size(Nload,3)
   Nload(:,:,ii) =  Nload(:,:,ii)./rch; % Kg/year -> Kg/m^3 / year
end
Nload = Nload.*1000; % Kg/m^3 / year -> mg/lt / year
%%
Nload_conc = Nload;
save('N_mass_load', 'Nload_conc','-append');
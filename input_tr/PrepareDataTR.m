%% Read wells from the 2000 run
wells = readWells('../input/well_data4.npsat');
%% calculate ratio
Qtarget = 3195370;
rat = Qtarget/abs(sum(wells(:,5)));
wells(:,5) = wells(:,5)*rat;
%% write new wells
fid = fopen('well_dataTR.npsat','w');
fprintf(fid, '%d\n', size(wells, 1));
fprintf(fid, '%0.3f %0.3f %0.3f %0.3f %0.3f\n', wells');
fclose(fid);
%% Modify Recharge
rch = read_Scattered('../input/rch_data.npsat', 2);
rch.v = rch.v*rat;
writeScatteredData('rch_dataTR.npsat', ...
                   struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'),...
                   [rch.p(:,1) rch.p(:,2) rch.v]);
%% Modify Stream recharge
streams = readStreams('../input/stream_data.npsat');
for ii = 1:length(streams)
    streams(ii,1).Q = rat*streams(ii,1).Q;
end
writeStreams('stream_dataTR.npsat',streams);
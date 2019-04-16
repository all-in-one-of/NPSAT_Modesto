%
%% Here I'll read Mehrdad's input files and convert them to npsat format
mehrdad_dir = '/home/giorgk/Documents/UCDAVIS/Mehrdad/NPSAT_input_files_BAU/';
mehrdad_dir = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/NPSAT_input_files_BAU/';
input_dir = 'input/';
output_dir = '/home/giorgk/Documents/UCDAVIS/Mehrdad/output/';
%% Top and bottom
% Read the xlc after converted to csv
raw_elev_data = csvread([mehrdad_dir 'bot_top_elv.csv'],1,0);
%%
xcoord = @(ii)200 + (ii-1).*400;
ycoord = @(ii)61200 - (200 + (ii-1).*400);
for ii = 1:16
    LAYELEV{ii,1} = raw_elev_data(raw_elev_data(:,3) == ii,:);
end
ftop = scatteredInterpolant(xcoord(LAYELEV{1,1}(:,2)), ycoord(LAYELEV{1,1}(:,1)), LAYELEV{1,1}(:,5),'linear','nearest');
fbot = scatteredInterpolant(xcoord(LAYELEV{1,1}(:,2)), ycoord(LAYELEV{1,1}(:,1)), LAYELEV{16,1}(:,4),'linear','nearest');

% extended grid
[Xgrid, Ygrid] = meshgrid(-200:400:55200+400, -200:400:61400);
Xgrid = reshape(Xgrid, size(Xgrid,1)*size(Xgrid,2), 1);
Ygrid = reshape(Ygrid, size(Ygrid,1)*size(Ygrid,2), 1);
top = ftop(Xgrid, Ygrid);
bot = fbot(Xgrid, Ygrid);

%% top elev
writeScatteredData([input_dir 'top_elev.npsat'], ...
                   struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'),...
                   [Xgrid Ygrid top]);
               
%% bottom elev
writeScatteredData([input_dir 'bot_elev.npsat'], ...
                   struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'),...
                   [Xgrid Ygrid bot]);

%% General Head boundary conditions
raw_GHB_data = csvread([mehrdad_dir 'GH_BC.csv'],1,0);
% get the mesh outline
outline = shaperead('ActiveModelArea');
%%
% split the outline to separate lines
LBC = [];
for ii = 1:length(outline(1,1).X)-2
    LBC = [LBC; outline(1,1).X(ii) outline(1,1).Y(ii) outline(1,1).X(ii+1) outline(1,1).Y(ii+1)];
end
LBC(2,:) = []; % THis is the east boundary that is no flow
%%
% concatenate the colinear consecutive segments
dlt = []; ids = []; val = [];
for ii = 1:size(LBC,1)-1
   if LBC(ii,3) == LBC(ii+1,1) && LBC(ii,4) == LBC(ii+1,4)
      if sum(LBC(ii,2) - [LBC(ii,4) LBC(ii+1,2) LBC(ii+1,4)]) == 0 || sum(LBC(ii,1) - [LBC(ii,3) LBC(ii+1,1) LBC(ii+1,3)])==0
          ids = [ids;ii];
          val = [val;LBC(ii,1) LBC(ii,2) LBC(ii+1,3) LBC(ii+1,4)];
          dlt = [dlt; ii+1];
      end
   end
end
LBC(ids,:) = val;
LBC(dlt,:) = [];
% rearrange the lines to be one continuous line
LBC = [LBC; LBC(1,:)];
LBC(1,:) = [];
%%
clf
hold on
for ii = 1:size(LBC,1)
   plot(LBC(ii,[1 3]), LBC(ii,[2 4])) 
end
plot(raw_GHB_data(:,1),raw_GHB_data(:,2),'xr')
%% 
% for each line find the GHB nodes that are closer than 200 m
for ii = 1:size(LBC,1)
    dst = Dist_Point_LineSegment(raw_GHB_data(:,1),raw_GHB_data(:,2), LBC(ii,:));
    ids = find(dst < 210);
    LBC_ids{ii,1} = ids;
end
% manulally fix some errors
LBC_ids{2,1}(raw_GHB_data(LBC_ids{2,1},1) > 13200,:) = [];
LBC_ids{3,1}(raw_GHB_data(LBC_ids{3,1},2) < 400,:) = [];
LBC_ids{162,1}(raw_GHB_data(LBC_ids{162,1},1) > 11200,:) = [];


%%
for ii = 1:size(LBC,1)
    [t, tmp] = projectPoints2Line(raw_GHB_data(LBC_ids{ii,1},[1 2]), LBC(ii,:));
    t_uniq = unique(t);
    V{ii,1} = nan(length(t_uniq), 1+1+2*15); % x v1 (z v) 16 times
    Ldist = pdist([LBC(ii,1:2);LBC(ii,3:4)]);
    V{ii,1}(:,1) = t_uniq.*Ldist;
    %V{ii,1}(1,1) = 0;
    %V{ii,1}(end,1) = Ldist;
    
    % for each point alont the line find the layer Z and value
    for jj = 1:length(t_uniq)
       id = find(t ==  t_uniq(jj));
       lzv = raw_GHB_data(LBC_ids{ii,1}(id),[6 3 7]);
       % find missing layer values
       miss_lay = setdiff([1:16]',lzv(:,1));
       v_miss = interp1(lzv(:,1),lzv(:,3), miss_lay,'nearest','extrap');
       z_miss = interp1(lzv(:,1),lzv(:,2), miss_lay,'linear','extrap');
       Vtable = nan(16,2);
       Vtable(lzv(:,1),1) = lzv(:,2);
       Vtable(lzv(:,1),2) = lzv(:,3);
       Vtable(miss_lay,1) = z_miss;
       Vtable(miss_lay,2) = v_miss;
       if any(diff(Vtable(:,1))>0)
           display(['Found errors on ' num2str(ii) ' at point ' num2str(jj)]);
       end
       V{ii,1}(jj,2:end) = [Vtable(1,2) reshape([(Vtable(1:end-1,1)+Vtable(2:end,1))/2 Vtable(2:end,2)]',1,30)];
    end
    V{ii,1} = [V{ii,1}(1,:); V{ii,1}; V{ii,1}(end,:)];
    V{ii,1}(1,1) = 0;
    V{ii,1}(end,1) = Ldist;
    if any(diff(V{ii,1}(:,1))<0)
        display(['Found errors on ' num2str(ii)]);
    end
end
%%
% average the corner points
for ii = 1:size(LBC,1)-1
    av = (V{ii,1}(end,2:end) + V{ii+1,1}(1,2:end))/2;
    V{ii,1}(end,2:end) = av;
    V{ii+1,1}(1,2:end) = av;
end
%%
% write finally boundary conditions to files
fid = fopen([input_dir 'BC_main.npsat'],'w');
fprintf(fid, '%d\n',size(LBC,1));
for ii = 1:size(LBC,1)
   fprintf(fid, 'EDGE 2 BC_files/lnbndr%d.npsat\n', ii);
   fprintf(fid, '%0.3f %0.3f\n', LBC(ii,1:2));
   fprintf(fid, '%0.3f %0.3f\n', LBC(ii,3:4));
   writeScatteredData([input_dir 'BC_files/lnbndr' num2str(ii) '.npsat'],...
       struct('PDIM', 1, 'TYPE', 'VERT', 'MODE', 'STRATIFIED'),...
       V{ii,1});
end
fclose(fid);



%% No flow boundary conditions
raw_noflow_data = csvread([mehrdad_dir 'no_flow.csv'],1,0);
% create a raster file with the active cells
active = ones(153, 137);
active(sub2ind(size(active),raw_noflow_data(:,1),raw_noflow_data(:,2))) = 0;
active(:,60:end) = 1;
WriteAscii4Raster('ActiveModelArea',active, 0, 0, 400, -9);
% write raster as grid
[Xgrid, Ygrid] = meshgrid(1:137, 1:153); 
ids = sub2ind(size(active),Ygrid, Xgrid);
ids(active == 0) = 0;
WriteAscii4Raster('ActiveModelCells',ids, 0, 0, 400, 0);
%
%% make mesh
% read the converted from raster ActiveModelCells shapefile and make a
% mesh file
%msh_shp = shaperead('ActiveModelCells');
p = []; msh = [];
for ii = 1:size(msh_shp,1)
    ii
    inds = [];
    pnts = [msh_shp(ii,1).X(1:4); msh_shp(ii,1).Y(1:4)]';
    for jj = 1:size(pnts,1)
        if isempty(p)
            p = pnts(jj,:);
            inds = [inds size(p,1)];
        else
            dst = sqrt((p(:,1) - pnts(jj,1)).^2 + (p(:,2) - pnts(jj,2)).^2);
            [cc, dd] = min(dst);
            if cc < 0.1
                inds = [inds dd];
            else
                p = [p; pnts(jj,:)];
                inds = [inds size(p,1)];
            end
        end
    end
    msh = [msh; inds];
end
% write mesh file
writeMeshfile([input_dir 'init_mesh.npsat'], p, msh);
%% Groundwater Recharge
% raw_rch_data = csvread([mehrdad_dir 'Rch_BAU.csv'],1,0);
% rch_nonzero = raw_rch_data(raw_rch_data(:,7)~=0,:);
CBC = readModflowFlowdata([mehrdad_dir 'MF_BAU_scheme9.crc']);
RCH_cbc = CBC(1,1).data;
CBC = readModflowFlowdata([mehrdad_dir 'MF_BAU_scheme9.cbe']);
ET_cbc = CBC(1,1).data;

rch = sum(RCH_cbc,3);
et = sum(ET_cbc,3);
rch = rch - et;
[Xgrid, Ygrid] = meshgrid(200:400:54800,200:400:61000);
Xgrid = reshape(Xgrid, 137*153,1);
Ygrid = reshape(flipud(Ygrid), 137*153,1);
rch = reshape(rch, 137*153,1)/(400*400);
frch = scatteredInterpolant(Xgrid, Ygrid, rch,'linear','nearest');
% extended grid
[Xgrid, Ygrid] = meshgrid(-200:400:55200+400, -200:400:61400);
Xgrid = reshape(Xgrid, size(Xgrid,1)*size(Xgrid,2), 1);
Ygrid = reshape(Ygrid, size(Ygrid,1)*size(Ygrid,2), 1);
rch = frch(Xgrid, Ygrid);
writeScatteredData([input_dir 'rch_data.npsat'], ...
                   struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'),...
                   [Xgrid Ygrid rch]);
%% K Porosity
% X Y Z row col lay Ss Por Kx Ky Kz
raw_K_data = csvread([mehrdad_dir 'Ss_K_n.csv'],1,0);
%%
% make a list of 2D points where the properties are defined
pnts = raw_K_data(raw_K_data(:,6) == 1,1:2);
% add an extra row or column outside the domain
id = find(pnts(:,2) == 200);
pnts = [pnts; pnts(id,1) -200*ones(length(id),1)];
id = find(pnts(:,2) == 61000);
pnts = [pnts; pnts(id,1) 61400*ones(length(id),1)];
id = find(pnts(:,1) == 200);
pnts = [pnts; -200*ones(length(id),1) pnts(id,2)];
id = find(pnts(:,1) == 54600);
pnts = [pnts; 55000*ones(length(id),1) pnts(id,2)];
%%
% create interpolants
KX = []; KZ = []; PR = []; EL = [];
for ii = 1:16
    vv = raw_K_data(raw_K_data(:,6) == ii,:);
    flay = scatteredInterpolant(vv(:,1), vv(:,2), vv(:,3),'linear','nearest');
    fKx = scatteredInterpolant(vv(:,1), vv(:,2), vv(:,10),'linear','nearest');
    fKz = scatteredInterpolant(vv(:,1), vv(:,2), vv(:,11),'linear','nearest');
    fPr = scatteredInterpolant(vv(:,1), vv(:,2), vv(:,8),'linear','nearest');
    EL = [EL flay(pnts(:,1), pnts(:,2))];
    KX = [KX fKx(pnts(:,1), pnts(:,2))];
    KZ = [KZ fKz(pnts(:,1), pnts(:,2))];
    PR = [PR fPr(pnts(:,1), pnts(:,2))];
end
EL = (EL(:,1:end-1) + EL(:,2:end))/2;
%%
DATA_Kx = [pnts KX(:,1)];
DATA_Kz = [pnts KZ(:,1)];
DATA_Pr = [pnts PR(:,1)];
for ii = 1:15
    DATA_Kx = [DATA_Kx EL(:,ii) KX(:,ii+1)];
    DATA_Kz = [DATA_Kz EL(:,ii) KZ(:,ii+1)];
    DATA_Pr = [DATA_Pr EL(:,ii) PR(:,ii+1)];
end
%%
% write the files
 writeScatteredData([input_dir 'KX_data.npsat'],...
       struct('PDIM', 2, 'TYPE', 'FULL', 'MODE', 'STRATIFIED'),...
       DATA_Kx);
writeScatteredData([input_dir 'KZ_data.npsat'],...
       struct('PDIM', 2, 'TYPE', 'FULL', 'MODE', 'STRATIFIED'),...
       DATA_Kz);
writeScatteredData([input_dir 'Por_data.npsat'],...
       struct('PDIM', 2, 'TYPE', 'FULL', 'MODE', 'STRATIFIED'),...
       DATA_Pr);
%% Find the elevation of the water table.
load('SteadyFlowHeadMatrixAll.mat')
SteadyFlowHead = zeros(153,137,16);
for ii = 1:16
    SteadyFlowHead(:,:,ii) = SteadyFlowHeadMatrixAll((ii-1)*153+1:ii*153,:);
end

%% Wells
% row,column,layer,x,y,pumpage (m3/d),Elv of bottom of screen (m),Elv of top of screen (m)
raw_well_data = csvread([mehrdad_dir 'well.csv'],1,0);
well_IJ = unique(raw_well_data(:,1:2),'rows');
clear well_data WellData
well_data = nan(size(well_IJ,1),5);
for ii = 1:size(well_IJ,1)
    % Find the heads in the steady state solution for modflow
    H = reshape(SteadyFlowHead(well_IJ(ii,1), well_IJ(ii,2),:),16,1);
    id = find(raw_well_data(:,1) == well_IJ(ii,1) & raw_well_data(:,2) == well_IJ(ii,2));
    well_data(ii,1) = raw_well_data(id(1),4);
    well_data(ii,2) = raw_well_data(id(1),5);
    well_data(ii,3) = max(raw_well_data(id,8));
    well_data(ii,4) = min(raw_well_data(id,7));
    well_data(ii,5) = sum(raw_well_data(id,6));
    
    WellData(ii,1).I = well_IJ(ii,1);
    WellData(ii,1).J = well_IJ(ii,2);
    WellData(ii,1).X = raw_well_data(id(1),4);
    WellData(ii,1).Y = raw_well_data(id(1),5);
    WellData(ii,1).Top_init = max(raw_well_data(id,8));
    WellData(ii,1).Bop_init = min(raw_well_data(id,7));
    WellData(ii,1).SL = WellData(ii,1).Top_init - WellData(ii,1).Bop_init;
    WellData(ii,1).WT = H(find(H > -9998,1));
    WellData(ii,1).D_wt_top = WellData(ii,1).WT - WellData(ii,1).Top_init;
    WellData(ii,1).Q = well_data(ii,5);
end
%% 
% Assign to the wells that have top of the screen higher that the water table
% a random depth based on the distribution of the other wells
D_wt_top = [WellData.D_wt_top]';
[d, fd] = ecdf(D_wt_top(D_wt_top > 0));
for ii = 1:length(WellData)
    if WellData(ii,1).D_wt_top > 0
        WellData(ii,1).D_wt_top_mod = WellData(ii,1).D_wt_top;
    else
        WellData(ii,1).D_wt_top_mod = interp1(d,fd,rand);
    end
end
%%
fid = fopen([input_dir 'well_data.npsat'],'w');
fprintf(fid, '%d\n', size(well_data, 1));
fprintf(fid, '%0.3f %0.3f %0.3f %0.3f %0.3f\n', well_data');
fclose(fid);
%% streams
% raw_stream_data = csvread([mehrdad_dir 'streams.csv'],1,0);
CBC = readModflowFlowdata([mehrdad_dir 'MF_BAU_scheme9.cs1']);
STRM_cbc = CBC(1,1).data;
strm = sum(STRM_cbc,3);
[Xgrid, Ygrid] = meshgrid(200:400:54800,200:400:61000);
Xgrid = reshape(Xgrid, 137*153,1);
Ygrid = reshape(flipud(Ygrid), 137*153,1);
strm = reshape(strm, 137*153,1);
% ids = find(strm~=0);
% strmID = strm;
% strmID(ids) = ids;
% WriteAscii4Raster('StreamBAU9',strmID, 0,0,400,0);
%%
strm_polys = shaperead('StreamPolys');
%%
clear polys
p = []; polys.Q = []; polys.poly = [];
cnt = 1;
for ii = 1:size(strm_polys, 1)
    if length(strm_polys(ii,1).X) == 6
        tmp = [mean(strm_polys(ii,1).X(1:end-2)) mean(strm_polys(ii,1).Y(1:end-2))];
        add_this = false;
        if isempty(p)
            add_this = true;
        else
            dst = sqrt((p(:,1) - tmp(1)).^2 + (p(:,2) - tmp(2)).^2);
            if min(dst) > 0.1
                add_this = true;
            end
        end
        
        if add_this
            dst = sqrt((Xgrid - tmp(1)).^2 + (Ygrid - tmp(2)).^2);
            [cc, dd ] = min(dst);
            if strm(dd) ~= 0
                p = [p;tmp];
                polys(cnt,1).Q = strm(dd)/(400*400);
                polys(cnt,1).poly = [strm_polys(ii,1).X(1:end-2)' strm_polys(ii,1).Y(1:end-2)'];
                cnt = cnt+1;
            end
        end
    end
end
writeStreams([input_dir 'stream_data.npsat'], polys);
%% Adjust the wells so that they are below the water table
% read the wells 
load('WellData_Init.mat', 'WellData')
% read the water table cloud points
WTC = ReadWaterTableCloudPoints([output_dir 'Modst_ref7_top_' num2str(4,'%03d') '_'], 32, true);
%% Create an interpolant;
Fnewtop = scatteredInterpolant(WTC(:,1), WTC(:,2), WTC(:,3));
%% Update well elevation
for ii = 1:length(WellData)
    xw = WellData(ii,1).X;
    yw = WellData(ii,1).Y;
    wtb = Fnewtop(xw, yw);
    d = WellData(ii,1).D_wt_top_mod;
    sl = WellData(ii,1).SL;
    top = wtb - d;
    bot = top - sl;
    wells(ii,:) = [xw yw top bot WellData(ii,1).Q];
end
%% write wells
fid = fopen([input_dir 'well_data4.npsat'],'w');
fprintf(fid, '%d\n', size(wells, 1));
fprintf(fid, '%0.3f %0.3f %0.3f %0.3f %0.3f\n', wells');
fclose(fid);
%% Write all water table cloud iterations
output_dir = '/media/giorgk/DATA/giorgk/Documents/NPSAT_Modesto/Results/tempOutput/';
for ii = 0:4
    WTC = ReadWaterTableCloudPoints([output_dir 'Modst_ref7_top_' num2str(ii,'%03d') '_'], 32, true);
end
%% Create an initial top based on a previous solution
top0 = read_Scattered([input_dir 'top_elev.npsat'],2);
% read the solution that we want to make it initial elevation
WTC = ReadWaterTableCloudPoints([output_dir 'Modst_ref7_top_' num2str(4,'%03d') '_'], 32, true);
Fnewtop = scatteredInterpolant(WTC(:,1), WTC(:,2), WTC(:,3));
elev_new = Fnewtop(top0.p(:,1), top0.p(:,2));
writeScatteredData([input_dir 'top_elev1.npsat'], ...
                   struct('PDIM', 2, 'TYPE', 'HOR', 'MODE', 'SIMPLE'),...
                   [top0.p(:,1), top0.p(:,2),elev_new]);

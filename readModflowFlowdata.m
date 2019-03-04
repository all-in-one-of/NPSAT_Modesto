function CBC = readModflowFlowdata(filename)
%http://gwlab.blogspot.gr/2012/08/reading-modflow-output-files-again.html

fid = fopen(filename, 'r');
k = 0;
while ~feof(fid)
    k = k+1;
    CBC(k,1).KSPT = fread(fid, 1, 'uint');
    CBC(k,1).KPER = fread(fid, 1, 'uint');
    CBC(k,1).DESC = fread(fid, 16, 'char');
    CBC(k,1).DESC = char(CBC(k,1).DESC');
    CBC(k,1).NCOL = fread(fid, 1, 'uint');
    CBC(k,1).NROW = fread(fid, 1, 'uint');
    CBC(k,1).NLAY = fread(fid, 1, 'uint');
    fprintf('Reading %s for time step %i, stress period %i\n', CBC(k,1).DESC, CBC(k,1).KSPT, CBC(k,1).KPER);
    
    if CBC(k,1).NLAY > 0
        % If NLAY is greater than zero, a 3D array of real numbers follows NLAY.  
        % The number of values is NCOL x NROW x NLAY. To read it, you can use 
        % a loop over the layers that contains a loop over rows that contains 
        % a loop over columns.
        CBC(k,1).data = read_modflow_array(fid, CBC(k,1).NROW, CBC(k,1).NCOL, CBC(k,1).NLAY);
    elseif CBC(k,1).NLAY <= 0
        error('negative NLAY cases not implemented yet')
        ITYPE = fread(fid, 1, 'uint');
        DELT = fread(fid, 1, 'float');
        PERTIM = fread(fid, 1, 'float');
        TOTIM = fread(fid, 1, 'float');
        if ITYPE == 0 || ITYPE == 1 
            data = read_modflow_array(fid, NROW, NCOL, abs(NLAY));
        elseif ITYPE == 3
            for i=1:NROW*NCOL
                data(i,1)=fread(fid,1,'uint');
            end
            
        end
        
        
    end
    
    
end


fclose(fid);

function data = read_modflow_array(fid, NR, NC, NL)
data = zeros(NR, NC, NL);
for il = 1:NL
    for ir = 1:NR
        for ic = 1:NC
            data(ir,ic, il)=fread(fid, 1, 'float');
        end
    end
end
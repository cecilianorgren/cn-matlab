function out = gmmCrossEventStatistics_v0(listFiles,varargin)

p = inputParser;
addParameter(p,'version',[]);
addParameter(p,'thrL1',0.05:0.05:0.3);
addParameter(p,'coldArgs',{'maxColdEV',500,'minVabs',50});
parse(p,varargin{:});

if ~iscell(p.Results.version), p.Results.version = {p.Results.version}; end

nVer = numel(p.Results.version);
out = struct;
out.nEvents = numel(listFiles);
out.nVDFs = 0;


% Initialize output structure
for iVer = 1:nVer
    switch lower(p.Results.version{iVer})
      case 'minbick'   
        out.M_bicmin = [];  
        out.minT_bicmin = [];  
      case 'chordbick'
        out.M_bicchord = [];  
      case 'percentage_cold'
        %out = struct('minBIC',[],'chordBIC',[],'L1',[]);
      
      case ''
    
    end
end


fprintf('iFile = ')
for iFile = 1:numel(listFiles)
  fprintf('%g,',iFile)
  fileName = [listFiles(iFile).folder filesep listFiles(iFile).name];
  GM = load(fileName,'GM_save').GM_save; 
  
  iPreDF = GM.time.tlim([GM.userdata.db_tFF_start GM.userdata.db_tDF]);
  iPostDF = GM.time.tlim([GM.userdata.db_tDF GM.userdata.db_tFF_stop]);

  GMpreDF = GM.sliceTime(iPreDF);
  GMpostDF = GM.sliceTime(iPostDF);

  out.nVDFs = out.nVDFs + GM.nT;
%  out.pre.nVDFs = out.pre.nVDFs + GMpreDF.nT;
%  out.post.nVDFs = out.post.nVDFs + GMpostDF.nT;
%  M_bicmin = accumData(GM,'M_bicmin',M_bicmin);


  for iVer = 1:nVer
    switch lower(p.Results.version{iVer})
      case 'minbick'
        G = GM.selectCol(GM.kneeK('method','min'));
        out.M_bicmin = [out.M_bicmin; G.nComp]; % append nComp
        data = G.sortBy('T').byNComp;
        for iM = 1:numel(data)
          minT = [data(iM).it data(iM).T(:,1)];
          nMinT = [data(iM).it data(iM).w(:,1)];
          if iM == 1 && iFile == 1
            out.minT_bicmin = minT;
            out.nMinT_bicmin = nMinT;
          else
            out.minT_bicmin = [out.minT_bicmin; minT];
            out.nMinT_bicmin = [out.nMinT_bicmin; nMinT];
          end
          [~,isort] = sort(out.minT_bicmin(:,1));
          out.minT_bicmin = out.minT_bicmin(isort,:);
          out.nMinT_bicmin = out.nMinT_bicmin(isort,:);
        end
      case 'chordbick'      
        G = GM.selectCol(GM.kneeK('method','chord'));
        out.M_bicchord = [out.M_bicchord; G.nComp]; % append nComp
        data = G.sortBy('T').byNComp;
        for iM = 1:numel(data)
          minT = [data(iM).it data(iM).T(:,1)];
          nMinT = [data(iM).it data(iM).w(:,1)];
          if iM == 1 && iFile == 1
            out.minT_bicchord = minT;
            out.nMinT_bicchord = nMinT;
          else
            out.minT_bicchord = [out.minT_bicchord; minT];
            out.nMinT_bicchord = [out.nMinT_bicchord; nMinT];
          end
          [~,isort] = sort(out.minT_bicchord(:,1));
          out.minT_bicchord = out.minT_bicchord(isort,:);
          out.nMinT_bicchord = out.nMinT_bicchord(isort,:);
        end
      case 'percentage_cold'
        %out = struct('minBIC',[],'chordBIC',[],'L1',[]);
      
      case ''
    
    end
  end
end
out.minT_bicmin = out.minT_bicmin(:,2);
out.minT_bicchord = out.minT_bicchord(:,2);
out.nMinT_bicmin = out.nMinT_bicmin(:,2);
out.nMinT_bicchord = out.nMinT_bicchord(:,2);
fprintf('\n',iFile)


  function out = accumData(in,version,GM)

    
  end
end
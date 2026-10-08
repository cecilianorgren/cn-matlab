files = dir('/Users/cecilianorgren/IRFU Dropbox/Cecilia Norgren/Days/survey_sorted_gender/');
files = dir('/Users/cecilianorgren/IRFU Dropbox/Cecilia Norgren/Days/survey_sorted_language/');
files(1:2) = []; % remove . and ..

nFiles = numel(files);
for iFile = 1:nFiles  
  file = [files(iFile).folder filesep files(iFile).name];
  data = read_edi_csv(file);

  
  question = data.question; 
  title_header = data.section;
  maxChars = 90;
  words = split(strtrim(string(question)));
  title_header = [title_header '-'  cat(1,words(1:2)')];
  title_header = char(join(title_header, ' '));
  words = words(3:end);
  lines = strings(0, 1);
  currentLine = "";

  for iWord = 1:numel(words)
    candidate = strtrim(currentLine + " " + words(iWord));
    if strlength(currentLine) > 0 && strlength(candidate) > maxChars
      lines(end+1) = currentLine;
      currentLine = words(iWord);
    else
      currentLine = candidate;
    end
  end

  if strlength(currentLine) > 0
    lines(end+1) = currentLine;
  end

  titlestring = cellstr(lines);

  groups = data.answers.Properties.VariableNames;
  groups = groups(2:end);

  answers = data.answers;
  nGroups = size(data.answers,2)-1;
  
  if istable(answers)
      options = string(answers{:,1});
      values = answers{:,2:end};    
      values(:,end+1) = sum(values,2); % add overall
  elseif isstruct(answers) && all(isfield(answers, {'option','count'}))
      options = string({answers.option});
      values = [answers.count];
  elseif isnumeric(answers)
      values = answers(:);
      options = string(1:numel(values));
  else
      options = unique(string(answers), 'stable');
      values = countcats(categorical(string(answers), options));
  end
  

  values_norm = 100*values./repmat(sum(values,1),[size(values,1) 1]);

  colors = [         0    0.4470    0.7410; ...
                0.8500    0.3250    0.0980; ...
                0.9290    0.6940    0.1250; ...
                0.4940    0.1840    0.5560; ...
                0.4660    0.6740    0.1880; ...
                0.3010    0.7450    0.9330; ...
                0.6350    0.0780    0.1840];
  %colors = [     0    0.4470    0.7410; ...
  %          0.8500    0.3250    0.0980; ...
  %          0.8 0.8 0.8];

  hca = subplot(2,1,1); % counts
  hca.ColorOrder = colors;
  b = barh(hca, values);
  hca.YTick = 1:numel(options);
  hca.YTickLabel = options;
  hca.Title.String = {title_header,titlestring{:}};
  hca.YDir = 'reverse';
  xlabel(hca, 'Count');
  legend(hca,{groups{:},'Overall'},'location','eastoutside','box','off')

  % Leave some room to the right for the labels 
  xmax = max(values(:)); 
  hca.XLim = [0 xmax * 1.12]; % Add values at the end of each bar 
  for j = 1:numel(b) 
    x = b(j).XEndPoints; 
    y = b(j).YEndPoints; 
    labels = compose('  %.1f%%', values_norm(:,j)); 
    text(hca, y, x, labels, ...
      'HorizontalAlignment', 'left', ... 
      'VerticalAlignment', 'middle', ... 
      'FontSize', 9); 
  end

  hca = subplot(2,1,2); % normalized
  b = barh(hca, values_norm);
  hca.YTick = 1:numel(options);
  hca.YTickLabel = options;
  hca.Title.String = {title_header,titlestring{:}};
  hca.YDir = 'reverse';
  xlabel(hca, 'Percentage');
  legend(hca,{groups{:},'Overall'},'location','eastoutside','box','off')
 
  % Leave some room to the right for the labels 
  hca.XLim = [0 max(values_norm(:)) * 1.12]; % Add percentages at the end of each bar 
  for j = 1:numel(b) 
    x = b(j).XEndPoints; 
    y = b(j).YEndPoints; 
    labels = compose('  %g', values(:,j)); 
    text(hca, y, x, labels, ... 
      'HorizontalAlignment', 'left', ... 
      'VerticalAlignment', 'middle', ... 
      'FontSize', 10); 
  end
  %pause

  h = findobj(gcf,'type','axes'); h = h(end:-1:1);
  for ip = 1:numel(h)
    h(ip).XGrid = 'on';
    h(ip).YGrid = 'on';
    h(ip).FontSize = 12;
    h(ip).Title.FontSize = 11;
    h(ip).ColorOrder = [colors(1:nGroups,:); 0.7 0.7 0.7];
    h(ip).Position(1) = 0.22;
    h(ip).Position(3) = 0.55;
  end

  str_print = [data.section '___' question];
  str_print = strrep(str_print,' ','_');
  str_print = strrep(str_print,'´','');
  str_print = char(str_print);
  try
    cn.print(str_print)
  catch
    str_print = str_print(1:50);
    cn.print(str_print)
  end
end
%disp(data.section)
%disp(data.question)
%disp(data.answers)

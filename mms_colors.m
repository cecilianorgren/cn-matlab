function colors = mms_colors(colorOption)
% MMS_COLORS 
%   colors = mms_colors('1234'); % mms1, mms2, mms3, mms4
%   colors = mms_colors('xyza'); % x, y, z, abs
%   colors = mms_colors('12yzab'); % mms1, mms2, y, z, abs, unknown (=yellow)
if strcmp(colorOption,'matlab')
  colors =   [   0    0.4470    0.7410;...
              0.8500    0.3250    0.0980;...
              0.9290    0.6940    0.1250;...
              0.4940    0.1840    0.5560;...
              0.4660    0.6740    0.1880;...
              0.3010    0.7450    0.9330;...
              0.6350    0.0780    0.1840];
elseif strcmp(colorOption,'blues')
  %colors = colormap('parula');
  colors = cn.cmap('bluepink');
  colors = colors(1:end/2,:);
  colors = colors(fix(linspace(1,size(colors,1),10)),:);
elseif 0 % official MMS colors, bad contrast
  colors = [];
  for ii = 1:numel(colorOption);    
    switch colorOption(ii)
      case '1' % mms1, black
        newColor = [0 0 0];
      case '2' % mms2, red
        newColor = [0.8 0.4 0];
      case '3' % mms3, bluish green
        newColor = [0 0.6 0.5];
      case '4' % mms4, sky blue
        newColor = [0.35 0.7 0.9];
      case 'x' % x
        newColor = [0.35 0.7 0.9];
      case 'y' % y
        newColor = [0 0.6 0.5];
      case 'z' % z
        newColor = [0.8 0.4 0];
      case 'a' % absolute value
        newColor = [0.3 0.3 0.3]; 
      otherwise 
        irf.log('warning',['Can''t recognize input ''' colorOption(ii) '''. Inserting a nice yellow instead.']);
        newColor = [0.8 0.6 0.0]; 
    end  
    colors(ii,1:3) = newColor;  
  end
else % better colors
  colors = [];
  for ii = 1:numel(colorOption);    
    switch colorOption(ii)
      case '1' % mms1, black
        newColor = [0 0 0];
      case '2' % mms2, red
        newColor = [1 0.2 0];
      case '3' % mms3, bluish green
        newColor = [0 0.8 0.0];
      case '4' % mms4, sky blue
        newColor = [0.1 0.4 1];
      case 'x' % x
        newColor = [0.1 0.4 1];
      case 'y' % y
        newColor = [0 0.7 0.0];
      case 'z' % z
        newColor = [1 0.2 0];
      case 'a' % absolute value
        newColor = [0.3 0.3 0.3]; 
      otherwise 
        irf.log('warning',['Can''t recognize input ''' colorOption(ii) '''. Inserting a nice yellow instead.']);
        newColor = [0.95 0.7 0.0]; 
    end  
    colors(ii,1:3) = newColor;  
  end
end
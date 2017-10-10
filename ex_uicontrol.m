function ex_uicontrol
    % Example code for uicontrol reference page

    % Create a figure and an axes to contain a 3-D surface plot.
    figure
    hax = axes('Units','pixels');
    surf(peaks)
    % Create a uicontrol object to let users change the colormap
    % with a pop-up menu. Supply a function handle as the object's 
    % Callback:
    uicontrol('Style', 'popup',...
           'String', 'hsv|hot|cool|gray',...
           'Position', [20 340 100 50],...
           'Callback', @setmap);       % Popup function handle callback
                                       % Implemented as a subfunction
    
    % Add a different uicontrol. Create a push button that clears 
    % the current axes when pressed. Position the button inside
    % the axes at the lower left. All uicontrols have default units
    % of pixels. In this example, the axes does as well.
    uicontrol('Style', 'pushbutton', 'String', 'Clear',...
        'Position', [20 20 50 20],...
        'Callback', 'cla');        % Pushbutton string callback
                                   % that calls a MATLAB function

    % Add a slider uicontrol to control the vertical scaling of the
    % surface object. Position it under the Clear button.
    uicontrol('Style', 'slider',...
        'Min',1,'Max',50,'Value',41,...
        'Position', [400 20 120 20],...
        'Callback', {@surfzlim,hax});   % Slider function handle callback
                                        % Implemented as a subfunction
   
    % Add a text uicontrol to label the slider.
    uicontrol('Style','text',...
        'Position',[400 45 120 20],...
        'String','Vertical Exaggeration')
end


function setmap(hObj,event) %#ok<INUSD>
    % Called when user activates popup menu 
    val = get(hObj,'Value');
    if val ==1
        colormap(jet)
    elseif val == 2
        colormap(hsv)
    elseif val == 3
        colormap(hot)
    elseif val == 4
        colormap(cool)
    elseif val == 5
        colormap(gray)
    end
end

function surfzlim(hObj,event,ax) %#ok<INUSL>
    % Called to set zlim of surface in figure axes
    % when user moves the slider control 
    val = 51 - get(hObj,'Value');
    zlim(ax,[-val val]);
end
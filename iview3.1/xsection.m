%% Add cross section plots to the image, created by imagesc
% Left click - activate/deactivate trace mode

% Author R. Rokitski, 2006
function    [] = xsection(option)
if ~exist('option'),
    if isempty(get(gcf,'WindowButtonUpFcn')),
        option = 'initialize';
    else
        option = '';
    end
end

switch option,
    case 'initialize',
        % save current image data
        % rearrange the figure to add subplots
        % plot xsections and xhairs

        % Identify active axes with image data
        H.figure  =   gcf;
        H.children      =   get(H.figure, 'Children');
        switch length(H.children),
            case    0,
                error('Figure is empty'); return
            case    1,
                if ~strfind(get(get(H.figure, 'Children'), 'Type'), 'axes'),
                    error('No valid axes in current figure'); return
                end
                H.axes      =   gca;
            otherwise
                axesInd =   find(strcmp(get(H.children, 'Type'), 'axes') == 1);
                if length(axesInd) ~= 1,
                    error('Figure should have only one subplot'); return
                end
                H.axes      =  H.children(axesInd);
        end

        H.lines  =   get(H.axes, 'Children');
        if length(H.lines) ~= 1,
            error('Axes should have only one object');return
        end

        H.image  =   get(H.lines);
        if ~sum(strcmp(fieldnames(H.image), 'CData')),
            error('No image in current axes'); return
        end

        %% Create GUI in the figure
        x = H.image.CData;
        figure(H.figure)
        AXim = subplot(4,4,[1 2 3 5 6 7 9 10 11]);
        imagesc(x); axis xy; %colorbar;

        AXvp = subplot(4,4,[4 8 12]);
        plot(x(:, fix(size(x, 2)/2)), 1:1:size(x, 1))
        set(AXvp, 'YLim', [1 size(x, 1)])
        set(AXvp, 'XAxisLocation', 'top')
        set(AXvp, 'YAxisLocation', 'right')
        AXhp = subplot(4,4,[13 14 15]);
        plot(x(fix(size(x, 1)/2),:))
        set(AXhp, 'XLim', [1 size(x, 2)])
        
        AXvpTI =   get(AXvp, 'Position');
        AXhpTI =   get(AXhp, 'Position');

        Font.FontName   = 'Helvetica';
        Font.FontUnits  = 'Pixels';
        Font.FontSize   = 18;
        Font.FontWeight = 'bold';
        Font.FontAngle  = 'normal';
        
        H.xvalue      = uicontrol(H.figure,'Style','edit','Units','Normalized', Font, ...
                                  'Position',[AXvpTI(1) AXhpTI(1)+AXhpTI(1)-0.05  AXvpTI(3), 0.05],...
                                  'BackGroundColor',[ 0 0 1],'ForeGroundColor',[ 1 1 1],...
                                  'Tag','XVALUE',...
                                  'TooltipString','X value (Read Only)',...
                                  'String',num2str(fix(size(x, 2)/2)));
        H.yvalue      = uicontrol(H.figure,'Style','edit','Units','Normalized', Font, ...
                                  'Position',[AXvpTI(1) AXhpTI(1)+AXhpTI(1)-0.05 - 0.075  AXvpTI(3), 0.05],...
                                  'BackGroundColor',[ 0 0 1],'ForeGroundColor',[ 1 1 1],...
                                  'Tag','YVALUE',...
                                  'TooltipString','Y value (Read Only)',...
                                  'String', num2str(fix(size(x, 1)/2)));
                          
        set(H.figure,'WindowButtonUpFcn','xsection(''toggle'');');
        set(H.figure, 'Toolbar', 'Figure');
    case 'trace',
        axesHandles =   get(gcf, 'Children');
        controlHandles =  axesHandles(find(strcmp(get(axesHandles, 'Type'), 'uicontrol') == 1));
        axesHandles = axesHandles(find(strcmp(get(axesHandles, 'Type'), 'axes') == 1));
        axisHandle  =   axesHandles(end);

        x_rng = get(axisHandle,'Xlim');
        y_rng = get(axisHandle,'Ylim');
        temp = get(axisHandle,'currentpoint');
        set(controlHandles(1), 'String', num2str(fix(temp(1,1))));
        set(controlHandles(2), 'String', num2str(fix(temp(1,2))));

        if temp(1,1) >= x_rng(1) && temp(1,1) <= x_rng(2) && ...
                temp(1,2) >= y_rng(1) && temp(1,2) <= y_rng(2),
            AllLines    =   get(axisHandle, 'Children');
            delete(AllLines(find(strcmp(get(AllLines, 'Type'), 'line')==1)));
            x = get(get(axisHandle, 'Children'), 'CData');
            line(x_rng, temp(:,2));
            line(temp(:,1), y_rng);

           % plot(axesHandles(2), x(:, fix(temp(1,1))), 1:1:size(x, 1))
           % slavas line
		   axes(axesHandles(2)) 
		   plot(x(:, fix(temp(1,1))), 1:1:size(x, 1))
            
			set(axesHandles(2), 'YLim', y_rng);
            set(axesHandles(2), 'XAxisLocation', 'top')
            set(axesHandles(2), 'YAxisLocation', 'right')
          
			axes(axesHandles(1)) 
			plot( x(fix(temp(1,2)), :))
            %% plot(axesHandles(1), x(fix(temp(1,2)), :)) 
			% slavas line
			
            set(axesHandles(1), 'XLim', x_rng);

        end

    case 'toggle',
        if isempty(get(gcf,'WindowButtonMotionFcn')),
            set(gcf,'WindowButtonMotionFcn','xsection(''trace'');');
        else
            set(gcf,'WindowButtonMotionFcn','');            
        end
    otherwise
end

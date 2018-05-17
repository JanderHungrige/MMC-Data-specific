function [fig, ha, el] = template_video_resp

% Create figure
fig = figure('Position', [0 0 1024 768], 'Color', [1 1 1]);%1320 715]); 800x600

% Create subplot
ha(1) = axes('Parent',fig, 'Position',[0.05 0.55 0.334659090909091+0.05 0.281721585811883+0.05],...
    'YTick',[],...
    'XTickLabel',{},...
    'XTick',[],...
    'DataAspectRatio',[1 1 1],...
    'Xcolor',[1 1 1],...
    'Ycolor',[1 1 1]);
%hold(ha(1),'all');

% Create subplot
ha(2) = axes('Parent',fig, 'Position',[0.05 0.1 0.334659090909091+0.05 0.281721585811883+0.05],...
    'YTick',[],...
    'XTick',[],...
    'DataAspectRatio',[1 1 1],...
    'Xcolor',[1 1 1],...
    'Ycolor',[1 1 1]);
title(ha(2),'Pressure distribution',...
    'FontName','Tw Cen MT',...
    'FontSize',15);

% Create axes
ha(3) = axes('Parent',fig,'YTick',[],'XTick',[],...
    'Position',[0.5 0.799264705882353 0.4046+0.05 0.1157]);
box(ha(3),'on');
title(ha(3),'Reference ECG',...
    'FontName','Tw Cen MT',...
    'FontSize',15);

% Create axes
ha(7) = axes('Parent',fig,'YTick',[],'XTick',[],...
    'Position',[0.5 0.399632352941176 0.4046+0.05 0.1157]);
box(ha(7),'on');
title(ha(7),'Contactless Respiration',...
    'FontName','Tw Cen MT',...
    'FontSize',15);

% Create axes
ha(5) = axes('Parent',fig,...
    'Position',[0.63 0.22 0.27466 0.09]);
box(ha(5),'on');
title(ha(5),'Heart Rate',...
    'FontName','Tw Cen MT',...
    'FontSize',15);

% Create axes
ha(6) = axes('Parent',fig,...
    'Position',[0.63 0.10 0.27466 0.09]);
box(ha(6),'on');

% Create axes
ha(4) = axes('Parent',fig,...
    'Position',[0.5 0.6 0.4046+0.05 0.1157]);
box(ha(4),'on');
title(ha(4),'Contactless ECG',...
    'FontName','Tw Cen MT',...
    'FontSize',15);

% FLOWER
% ------

% Sensor array
% ------------
offsetx=0.5;
offsety=-0.4;
% Create ellipse 1
el(1)=annotation(fig,'ellipse',...
    [0.0481701279630793+offsetx 0.64312598686944+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');

% Create ellipse 2
el(2)=annotation(fig,'ellipse',...
    [0.0619201279630794+offsetx 0.610030748774204+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');

% Create ellipse 3
el(3)=annotation(fig,'ellipse',...
    [0.0557891755821269+offsetx 0.573999002742458+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');

% Create ellipse 4
el(4)=annotation(fig,'ellipse',...
    [0.0318010803440316+offsetx 0.559792653536106+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');

% Create ellipse 5
el(5)=annotation(fig,'ellipse',...
    [0.00751536605831742+offsetx 0.576141859885315+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');

% Create ellipse 6
el(6)=annotation(fig,'ellipse',...
    [0.00221774701069838+offsetx 0.613681542424997+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');

% Create ellipse 7
el(7)=annotation(fig,'ellipse',...
    [0.0192415565345078+offsetx 0.647252970996424+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');

% Create ellipse 8
el(8)=annotation(fig,'ellipse',...
    [0.0313248898678412+offsetx 0.603919637663093+offsety 0.02 0.03],...
    'FaceColor',[0.8 0 0],...
    'Color','none');



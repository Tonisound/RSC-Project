function make_syncmovie(full_path,flag_save,time_step,time_window)
% Generating sync  movie of behavior and calcium imaging
% Input argument fullpath: path to folder
% flag_save true/false to save movie
% time_step: temporal resolution (seconds)
% time_window: temporal window (seconds)

% if nargin == 0
%     folder = '/Users/tonio/Desktop';
%     recording_name = 'RecExample_RSC';
%     full_path = fullfile(folder,recording_name);
% end
if nargin < 4
    time_window = 30;
end
if nargin < 3
    time_step = .1;
end
if nargin < 2
    flag_save = false;
end


% Parameters
theta_rot = -pi/4; % radians
hd_offset = 3*pi/16; % radians ;
flags = [1,1,1]; % tuning curves, place fields, movie

Info = get_sync_data(full_path);

save_dir = fullfile(full_path,'processed');
if flag_save && ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% Loading behavior data
v_tracking = [];
data_tracking = [];
labels_tracking = [];

if Info.Metadata.flag_video
    
    % Video
    v_tracking = VideoReader(fullfile(Info.video_tracking.folder,Info.video_tracking.name));
    % Data
    %     data_tracking = readmatrix(fullfile(Info.csv_tracking.folder,Info.csv_tracking.name),'NumHeaderLines',6);
    data_tracking = csvread(fullfile(Info.csv_tracking.folder,Info.csv_tracking.name),7);
    fid = fopen(fullfile(Info.csv_tracking.folder,Info.csv_tracking.name));
    for i=1:6
        hl1 = fgetl(fid);
    end
    hl2 = fgetl(fid);
    fclose(fid);
    % Labels
    temp1=regexp(hl1,',','split');
    temp2=regexp(hl2,',','split');
    for i=1:length(temp1)
        labels_tracking = [labels_tracking;{strcat('[vid] ',char(temp1(i)),'-',char(temp2(i)))}];
    end
    
    % Evening NumFrames and data points
    if v_tracking.NumberOfFrames<size(data_tracking,1)
        data_tracking = data_tracking(1:v_tracking.NumberOfFrames,:);
    end
    
    % Defining t_relative_tracking & t_true_tracking
    t_relative_tracking = data_tracking(:,2);
    t_true_tracking = Info.Timing.t_optitrack(1)+t_relative_tracking;
    
end

% Getting Tracking Position
x_position = rescale(data_tracking(:,strcmp(labels_tracking,'[vid]Position-X')),-1,1);
y_position = rescale(data_tracking(:,strcmp(labels_tracking,'[vid]Position-Y')),-1,1);
% z_position = rescale(data_tracking(:,strcmp(labels_tracking,'[vid]Position-Z')),-1,1);
% Applying Rotation
x_position_rot = rescale(x_position*cos(theta_rot)-y_position*sin(theta_rot),-1,1);
y_position_rot = rescale(x_position*sin(theta_rot)+y_position*cos(theta_rot),-1,1);
x_position = x_position_rot;
y_position = y_position_rot;

% Getting Tracking Rotation
x_rotation = data_tracking(:,strcmp(labels_tracking,'[vid]Rotation-X'));
y_rotation = data_tracking(:,strcmp(labels_tracking,'[vid]Rotation-Y'));
z_rotation = data_tracking(:,strcmp(labels_tracking,'[vid]Rotation-Z'));

% Lowpass filtering
fs = 1000;
fpass1 = 2;
[B,A]  = butter(1,fpass1/(fs/2),'low');
t = t_relative_tracking;
t_q = t_relative_tracking(1):(1/fs):t_relative_tracking(end);
% x_rotation
y = x_rotation;
y_q = interp1(t,y,t_q);
y_lp = filtfilt(B,A,y_q);
x_rotation_lp = interp1(t_q,y_lp,t);
% y_rotation
y = y_rotation;
y_q = interp1(t,y,t_q);
y_lp = filtfilt(B,A,y_q);
y_fa = correct_angle(y_q,t_q);
y_rotation_lp = interp1(t_q,y_lp,t);
y_rotation_fa = interp1(t_q,y_fa,t);
% z_rotation
y = z_rotation;
y_q = interp1(t,y,t_q);
y_lp = filtfilt(B,A,y_q);
z_rotation_lp = interp1(t_q,y_lp,t);

x_rotation_rad = x_rotation*(pi/180);
y_rotation_rad = y_rotation*(pi/180);
z_rotation_rad = z_rotation*(pi/180);

x_rotation_lp_rad = x_rotation_lp*(pi/180);
y_rotation_lp_rad = y_rotation_lp*(pi/180);
y_rotation_fa_rad = y_rotation_fa*(pi/180);
z_rotation_lp_rad = z_rotation_lp*(pi/180);


% figure;
% ax1 = subplot(311);
% hold(ax1,'on');
% plot(t_true_tracking,x_rotation)
% plot(t_true_tracking,x_rotation_lp,'o')
% ax2 = subplot(312);
% hold(ax2,'on');
% plot(t_true_tracking,y_rotation)
% plot(t_true_tracking,y_rotation_lp,'o')
% ax3 = subplot(313);
% hold(ax3,'on');
% plot(t_true_tracking,z_rotation)
% plot(t_true_tracking,z_rotation_lp,'o')
% linkaxes([ax1,ax2,ax3],'x');


% Loading calcium data
v_miniscope = [];
data_miniscope = [];
labels_miniscope_cells = [];
labels_miniscope_accel = [];

if Info.Metadata.flag_miniscope
    
    % Video
    v_miniscope = VideoReader(fullfile(Info.video_miniscope.folder,Info.video_miniscope.name));
    
    % Data
    % data_miniscope_accel = readmatrix(fullfile(Info.csv_miniscope_accel.folder,Info.csv_miniscope_accel.name),'NumHeaderLines',1);
    data_miniscope_accel = csvread(fullfile(Info.csv_miniscope_accel.folder,Info.csv_miniscope_accel.name),1);
    % data_miniscope_cells = readmatrix(fullfile(Info.csv_miniscope_cells.folder,Info.csv_miniscope_cells.name))';
    data_miniscope_cells = csvread(fullfile(Info.csv_miniscope_cells.folder,Info.csv_miniscope_cells.name))';
    data_miniscope = [data_miniscope_accel,data_miniscope_cells];
    fid = fopen(fullfile(Info.csv_miniscope_accel.folder,Info.csv_miniscope_accel.name));
    hl = fgetl(fid);
    fclose(fid);
    
    % Labels
    temp1=regexp(hl,',','split');
    for i=1:length(temp1)
        labels_miniscope_accel = [labels_miniscope_accel;{strcat('[min] ',char(temp1(i)))}];
    end
    for i=1:size(data_miniscope_cells,2)
        labels_miniscope_cells = [labels_miniscope_cells;{sprintf('Cell-%03d',i)}];
    end
    
    % Evening NumberOfFrames and data points
    if v_miniscope.NumberOfFrames<size(data_miniscope,1)
        data_miniscope = data_miniscope(1:v_miniscope.NumberOfFrames,:);
    end
    
    % Defining t_relative_miniscope & t_true_miniscope
    t_relative_miniscope = ((data_miniscope(:,1)-data_miniscope(1,1))/1000);
    t_true_miniscope  = Info.Timing.t_miniscope(1)+t_relative_miniscope;
end


% Normalizing miniscope data
labels_miniscope = [labels_miniscope_accel;labels_miniscope_cells];
for i =1:size(data_miniscope_cells,2)
    data_miniscope_cells(:,i) = rescale(data_miniscope_cells(:,i),0,1);
end


%% Tuning curves
if flags(1)
    
    n_cells = size(data_miniscope_cells,2);
    n_columns = 10;
    n_rows = ceil(n_cells/n_columns);
    
    hd_miniscope = interp1(t_true_tracking,y_rotation_rad,t_true_miniscope);
    hd_miniscope_lp = interp1(t_true_tracking,y_rotation_lp_rad,t_true_miniscope);
    hd_miniscope_fa = interp1(t_true_tracking,y_rotation_fa_rad,t_true_miniscope);
    bin_edges = -pi:(2*pi)/36:pi;
    
    bin_counts = zeros(length(bin_edges)-1,n_cells);
    for i = 1:size(bin_counts,1)
        index_keep = (hd_miniscope>=bin_edges(i)).*(hd_miniscope<bin_edges(i+1));
        bin_counts(i,:) = mean(data_miniscope_cells(index_keep==1,:));
    end
    bin_counts(isnan(bin_counts))=0;
    
    bin_counts_lp = zeros(length(bin_edges)-1,n_cells);
    for i = 1:size(bin_counts_lp,1)
        index_keep = (hd_miniscope_lp>=bin_edges(i)).*(hd_miniscope_lp<bin_edges(i+1));
        bin_counts_lp(i,:) = mean(data_miniscope_cells(index_keep==1,:));
    end
    bin_counts_lp(isnan(bin_counts_lp))=0;
    
    bin_counts_fa = zeros(length(bin_edges)-1,n_cells);
    for i = 1:size(bin_counts_fa,1)
        index_keep = (hd_miniscope_fa>=bin_edges(i)).*(hd_miniscope_fa<bin_edges(i+1));
        bin_counts_fa(i,:) = mean(data_miniscope_cells(index_keep==1,:),'omitnan');
    end
    bin_counts_fa(isnan(bin_counts_fa))=0;
    
    bin_centers = bin_edges(1:end-1)+.5*(bin_edges(2)-bin_edges(1));
    %     mvl_x=sum(repmat(cos(bin_centers)',[1 n_cells]).*bin_counts);
    %     mvl_y=sum(repmat(sin(bin_centers)',[1 n_cells]).*bin_counts);
    %     mvl = sqrt(mvl_x.^2+mvl_y.^2);
    all_mvl = [];
    all_mvl_lp = [];
    all_mvl_fa = [];
    all_pd = [];
    all_pd_lp = [];
    all_pd_fa = [];
    pimp_factor = 2;
    
    f1_1 = figure('Name','Raw tuning curves');
    for counter = 1:n_cells
        pos = get_position(n_rows,n_columns,counter);
        pax = polaraxes('Parent',f1_1,'Position',pos);
        hold(pax,'on');
        polarplot(hd_miniscope,data_miniscope_cells(:,counter),'Parent',pax,...
            'LineStyle','none','Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b');
        %         set(pax,'RTick',[],'RTickLabel',[],'ThetaTick',[],'ThetaTickLabel',[]);
        set(pax,'RTickLabel',[],'ThetaTickLabel',[]);
        
        this_bin_counts = bin_counts(:,counter);
        polarhistogram('BinEdges',bin_edges,'BinCounts',pimp_factor*this_bin_counts,'Parent',pax,...
            'FaceAlpha',.75,'FaceColor',[.5 .5 .5],'EdgeColor','k');
        
        this_mvl = sqrt((sum(cos(bin_centers)'.*this_bin_counts)).^2+(sum(sin(bin_centers)'.*this_bin_counts)).^2)/10;
        [~,ind_pd] = max(this_bin_counts);
        this_pd = bin_centers(ind_pd);
        polarplot([this_pd this_pd],[0 pimp_factor*this_mvl],'Parent',pax,...
            'LineStyle','-','Color','r','LineWidth',2);
        pax.Title.String = strcat(char(labels_miniscope_cells(counter)),sprintf('[MVL=%.2f]',this_mvl));
        
        all_mvl = [all_mvl;this_mvl];
        all_pd = [all_pd;this_pd];
    end
    
    f1_2 = figure('Name','Lowpass tuning curves');
    for counter = 1:n_cells
        pos = get_position(n_rows,n_columns,counter);
        pax = polaraxes('Parent',f1_2,'Position',pos);
        hold(pax,'on');
        polarplot(hd_miniscope_lp,data_miniscope_cells(:,counter),'Parent',pax,...
            'LineStyle','none','Marker','.','MarkerFaceColor','r','MarkerEdgeColor','r');
        pax.Title.String = labels_miniscope_cells(counter);
        set(pax,'RTickLabel',[],'ThetaTickLabel',[]);
        
        this_bin_counts= bin_counts_lp(:,counter);
        polarhistogram('BinEdges',bin_edges,'BinCounts',pimp_factor*this_bin_counts,'Parent',pax,...
            'FaceColor',[.5 .5 .5],'EdgeColor','k');
        
        this_mvl = sqrt((sum(cos(bin_centers)'.*this_bin_counts)).^2+(sum(sin(bin_centers)'.*this_bin_counts)).^2)/10;
        [~,ind_pd] = max(this_bin_counts);
        this_pd = bin_centers(ind_pd);
        polarplot([this_pd this_pd],[0 pimp_factor*this_mvl],'Parent',pax,...
            'LineStyle','-','Color','r','LineWidth',2);
        pax.Title.String = strcat(char(labels_miniscope_cells(counter)),sprintf('[MVL=%.2f]',this_mvl));
        
        all_mvl_lp = [all_mvl_lp;this_mvl];
        all_pd_lp = [all_pd_lp;this_pd];
    end
    
    f1_3 = figure('Name','Corrected angle tuning curves');
    for counter = 1:n_cells
        pos = get_position(n_rows,n_columns,counter);
        pax = polaraxes('Parent',f1_3,'Position',pos);
        hold(pax,'on');
        polarplot(hd_miniscope_fa,data_miniscope_cells(:,counter),'Parent',pax,...
            'LineStyle','none','Marker','.','MarkerFaceColor','g','MarkerEdgeColor','g');
        pax.Title.String = labels_miniscope_cells(counter);
        set(pax,'RTickLabel',[],'ThetaTickLabel',[]);
        
        this_bin_counts = bin_counts_fa(:,counter);
        polarhistogram('BinEdges',bin_edges,'BinCounts',pimp_factor*this_bin_counts,'Parent',pax,...
            'FaceColor',[.5 .5 .5],'EdgeColor','k');
        
        this_mvl = sqrt((sum(cos(bin_centers)'.*this_bin_counts)).^2+(sum(sin(bin_centers)'.*this_bin_counts)).^2)/10;
        [~,ind_pd] = max(this_bin_counts);
        this_pd = bin_centers(ind_pd);
        polarplot([this_pd this_pd],[0 pimp_factor*this_mvl],'Parent',pax,...
            'LineStyle','-','Color','r','LineWidth',2);
        pax.Title.String = strcat(char(labels_miniscope_cells(counter)),sprintf('[MVL=%.2f]',this_mvl));
        
        all_mvl_fa = [all_mvl_fa;this_mvl];
        all_pd_fa = [all_pd_fa;this_pd];
    end
    
    f1_4 = figure('Name','Ordered Bin Counts');
    ax = subplot(1,3,1);
    [~,ind_ordered] = sort(all_pd,'ascend');
    imagesc('XData',180*bin_centers/pi,'YData',1:n_cells,'CData',bin_counts(:,ind_ordered)','Parent',ax);
    ax.XLim = [-180 180];
    ax.YLim = [.5 n_cells+.5];
    ax.YTick = 1:n_cells;
    ax.YTickLabel = labels_miniscope_cells(ind_ordered);
    ax.Title.String = 'Raw Bin Counts';
    ax = subplot(1,3,2);
    [~,ind_ordered] = sort(all_pd_lp,'ascend');
    imagesc('XData',180*bin_centers/pi,'YData',1:n_cells,'CData',bin_counts_lp(:,ind_ordered)','Parent',ax);
    ax.XLim = [-180 180];
    ax.YLim = [.5 n_cells+.5];
    ax.YTick = 1:n_cells;
    ax.YTickLabel = labels_miniscope_cells(ind_ordered);
    ax.Title.String = 'Lowpass Bin Counts';
    ax = subplot(1,3,3);
    [~,ind_ordered] = sort(all_pd_fa,'ascend');
    imagesc('XData',180*bin_centers/pi,'YData',1:n_cells,'CData',bin_counts_fa(:,ind_ordered)','Parent',ax);
    ax.XLim = [-180 180];
    ax.YLim = [.5 n_cells+.5];
    ax.YTick = 1:n_cells;
    ax.YTickLabel = labels_miniscope_cells(ind_ordered);
    ax.Title.String = 'Corrected Bin Counts';
    
    if flag_save
        % Saving mode
        set([f1_1;f1_2;f1_3;f1_4],'Units','normalized','OuterPosition',[0 0 1 1]);
        saveas(f1_1,fullfile(save_dir,'tuning-curves.jpg'),'jpeg');
        saveas(f1_2,fullfile(save_dir,'tuning-curves-lowpass.jpg'),'jpeg');
        saveas(f1_3,fullfile(save_dir,'tuning-curves-corrected.jpg'),'jpeg');
        saveas(f1_4,fullfile(save_dir,'bin-counts-ordered.jpg'),'jpeg');
        fprintf('Tuning Curves Saved [%s].\n',save_dir);
        close([f1_1;f1_2;f1_3;f1_4]);
    end
end


%% Place fields
if flags(2)
    
    n_cells = size(data_miniscope_cells,2);
    n_columns = 10;
    n_rows = ceil(n_cells/n_columns);
    
    x_position_miniscope = rescale(interp1(t_true_tracking,x_position,t_true_miniscope),0,1);
    y_position_miniscope = rescale(interp1(t_true_tracking,y_position,t_true_miniscope),0,1);
    list_thresh = 0:.1:1;
    t_colors = get_colors(length(list_thresh)-1,'jet');
    
    % Computing maps
    x_step = .02;
    y_step = .02;
    x_position_grid = 0:x_step:1;
    y_position_grid = 0:y_step:1;
    heatmap = NaN(length(x_position_grid)-1,length(y_position_grid)-1,n_cells);
    densitymap = NaN(length(x_position_grid)-1,length(y_position_grid)-1);
    for i=1:length(x_position_grid)-1
        ind_keep_x = (x_position_miniscope>=x_position_grid(i)).*(x_position_miniscope<x_position_grid(i+1));
        for j=1:length(y_position_grid)-1
            ind_keep_y = (y_position_miniscope>=y_position_grid(j)).*(y_position_miniscope<y_position_grid(j+1));
            ind_keep = ind_keep_x.*ind_keep_y;
            if sum(ind_keep)>0
                temp = mean(data_miniscope_cells(ind_keep==1,:));
                heatmap(i,j,:) = permute(temp,[1 3 2]);
                densitymap(i,j) = sum(ind_keep)/length(ind_keep);
            end
        end
    end
    
    % Plotting Occupancy
    f2_1 = figure('Name','Occupancy Map');
    ax = axes('Parent',f2_1);
    hold(ax,'on');
    im = imagesc('XData',x_position_grid(1:end-1)+x_step/2,'YData',y_position_grid(1:end-1)+y_step/2,'CData',log(densitymap)','Parent',ax);
    im.AlphaData = ~isnan(im.CData);
    ax.Title.String = 'Occupancy Map';
    colorbar(ax);
    ax.XLim = [x_position_grid(1) x_position_grid(end)];
    ax.YLim = [y_position_grid(1) y_position_grid(end)];
    set(ax,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'XDir','reverse','YDir','reverse');
    l = line('XData',x_position_miniscope,'YData',y_position_miniscope,'Parent',ax,'LineStyle','-','Color','r',...
        'LineWidth',.1,'Marker','.','MarkerFaceColor','r','MarkerEdgeColor','r');
    axis(ax,'equal');
    l.Color(4)=.5;
    ax.Visible='off';
    ax.Title.Visible='on';
    
    % Plotting Trajectories
    f2_2 = figure('Name','Trajectories with activity');
    for counter = 1:n_cells
        ydata = rescale(data_miniscope_cells(:,counter),0,1);
        pos = get_position(n_rows,n_columns,counter,[.05 .05 .001;.05 .05 .01]);
        ax = axes('Parent',f2_2,'Position',pos);
        hold(ax,'on')
        l = line('XData',x_position_miniscope,'YData',y_position_miniscope,'Parent',ax,'LineStyle','-','Color',[.25 .25 .25],...
            'Marker','none','MarkerFaceColor','r','MarkerEdgeColor','r');
        l.Color(4)=.25;
        for k=1:length(list_thresh)-1
            ind_keep = (ydata>list_thresh(k)).*(ydata<list_thresh(k+1));
            line('XData',x_position_miniscope(ind_keep==1),'YData',y_position_miniscope(ind_keep==1),...
                'Parent',ax,'LineStyle','none','Color',[.25 .25 .25],...
                'Marker','.','MarkerFaceColor',t_colors(k,:),'MarkerEdgeColor',t_colors(k,:));
        end
        ax.Title.String = labels_miniscope_cells(counter);
        set(ax,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'XDir','reverse','YDir','reverse');
        axis(ax,'equal');
        ax.Visible='off';
        ax.Title.Visible='on';
    end
    
    % Plotting place fiels
    f2_3 = figure('Name','Place fields');
    for counter = 1:n_cells
        pos = get_position(n_rows,n_columns,counter,[.05 .05 .001;.05 .05 .01]);
        ax = axes('Parent',f2_3,'Position',pos);
        grid(ax,'on');
        im = imagesc('XData',x_position_grid(1:end-1)+x_step/2,'YData',y_position_grid(1:end-1)+y_step/2,'CData',heatmap(:,:,counter)','Parent',ax);
        im.AlphaData = ~isnan(im.CData);
        ax.Title.String = labels_miniscope_cells(counter);
        ax.XLim = [x_position_grid(1) x_position_grid(end)];
        ax.YLim = [y_position_grid(1) y_position_grid(end)];
        ax.CLim = [0 .5];
        %         colorbar(ax);
        set(ax,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[],'XDir','reverse','YDir','reverse');
        axis(ax,'equal');
        ax.Visible='off';
        ax.Title.Visible='on';
    end
    
    if flag_save
        % Saving mode
        set([f2_1;f2_2;f2_3],'Units','normalized','OuterPosition',[0 0 1 1]);
        saveas(f2_1,fullfile(save_dir,'occupancy-map.jpg'),'jpeg');
        saveas(f2_2,fullfile(save_dir,'place-fields-trajectories.jpg'),'jpeg');
        saveas(f2_3,fullfile(save_dir,'place-fields-binned.jpg'),'jpeg');
        fprintf('Place Fields Saved [%s].\n',save_dir);
        close([f2_1;f2_2;f2_3]);
    end
    
end


%% Sync movie
if flags(3)
    
    % Building absolute time vector
    % tts_true_tracking = datestr(t_true_tracking/(24*3600),'HH:MM:SS.FFF');
    % tts_true_miniscope = datestr(t_true_miniscope/(24*3600),'HH:MM:SS.FFF');
    tts_relative_tracking = datestr(t_relative_tracking/(24*3600),'HH:MM:SS.FFF');
    tts_relative_miniscope = datestr(t_relative_miniscope/(24*3600),'HH:MM:SS.FFF');
    
    t = (0:time_step:Info.Timing.t_ephys(end))';
%     t = (30:time_step:60)';
    
    tts_t = datestr(t/(24*3600),'HH:MM:SS.FFF');
    
    labels_full = [labels_tracking;labels_miniscope];
    flags_traces = [1*ones(size(labels_tracking));2*ones(size(labels_miniscope))];
    
    % Trace Selection
    [ind_selected,v] = listdlg('Name','Trace Selection','PromptString','Select traces to display',...
        'SelectionMode','multiple','ListString',labels_full,'ListSize',[300 500]);
    if v==0 || isempty(ind_selected)
        return;
    end
    %     ind_selected = 3:8;
    
    labels_selected = labels_full(ind_selected);
    flags_selected = flags_traces(ind_selected);
    n_traces = length(labels_selected);
    g_colors = get_colors(n_traces,'jet');
    
    
    f = figure('Units','normalized','Name',sprintf('[%s] Synchronized Movie',Info.Metadata.recording_name),'OuterPosition',[0 0 1 1],'Tag','MainFigure');
    t1 = uicontrol('Units','normalized','Style','text','Parent',f,'BackgroundColor',[.75 .75 .75],'FontSize',14,'Position',[0 .975 .1 .025],'Tag','Text1');
    
    ax1 = axes('Parent',f,'Position',[.05 .65 .3 .3],'Tag','Ax1');
    axis(ax1,'equal');
    ax11 = axes('Parent',f,'Position',[.05 .35 .3 .3],'Tag','Ax11','YDir','reverse','XDir','reverse');
    axis(ax11,'equal');
    
    ax12 = polaraxes('Parent',f,'Position',[0 .6 .1 .1],'Tag','Ax12','FontSize',8);
    ax13 = polaraxes('Parent',f,'Position',[0 .75 .1 .1],'Tag','Ax13','FontSize',8);
    %     ax13b = polaraxes('Parent',f,'Position',[.1 .875 .05 .05],'Tag','Ax13b','FontSize',8);
    ax14 = polaraxes('Parent',f,'Position',[0 .45 .1 .1],'Tag','Ax14','FontSize',8);
    hold(ax12,'on');
    hold(ax13,'on');
    %     hold(ax13b,'on');
    hold(ax14,'on');
    
    t2 = uicontrol('Units','normalized','Style','text','Parent',f,'BackgroundColor',[.75 .75 .75],'FontSize',14,'Position',[0 .925 .1 .025],'Tag','Text2');
    
    ax2 = axes('Parent',f,'Position',[.05 .05 .3 .3],'Tag','Ax2');
    axis(ax2,'equal');
    t3 = uicontrol('Units','normalized','Style','text','Parent',f,'BackgroundColor',[.75 .75 .75],'FontSize',14,'Position',[0 .325 .1 .025],'Tag','Text3');
    
    t100 = uicontrol(f,'Units','normalized','Style','text','String','','BackgroundColor','k','Position',[.95-(.4/6) .01 (.4/6) .005]);
    t101 = uicontrol(f,'Units','normalized','Style','text','String',sprintf('%d s',time_window/3),'FontSize',10,'Position',[.95-(.4/6) .015 (.4/6) .025]);
    
    
    all_axes_traces = [];
    n_iqr = 3;
    r_colors = .25*ones(3,3);
    
    for i=1:n_traces
        %     ax = axes('Parent',f,'Position',[.4 .05+.9*((i-1)/n_traces) .55 (.9/n_traces)-.1/n_traces],'Tag',sprintf('AxTrace%d',i));
        ax = axes('Parent',f,'Position',[.4 .05+.9*((n_traces-i)/n_traces) .55 (.9/n_traces)-.1/n_traces],'Tag',sprintf('AxTrace%d',i));
        all_axes_traces = [all_axes_traces;ax];
        cur_label = labels_selected(i);
        
        if flags_selected(i)==1
            % tracking
            ind_trace = strcmp(labels_tracking,cur_label);
            data_selected = data_tracking(:,ind_trace==1);
            
            if strcmp(cur_label,'[vid]Position-X')
                data_selected = x_position;
            end
            if strcmp(cur_label,'[vid]Position-Y')
                data_selected = y_position;
            end
            line('XData',t_true_tracking,'YData',data_selected,'Parent',ax,'Color',g_colors(i,:),'Linewidth',2);
            
            if strcmp(cur_label,'[vid]Rotation-X')
                r_colors(1,:) = g_colors(i,:);
            elseif strcmp(cur_label,'[vid]Rotation-Y')
                r_colors(2,:) = g_colors(i,:);
            elseif strcmp(cur_label,'[vid]Rotation-Z')
                r_colors(3,:) = g_colors(i,:);
            end
            
        elseif flags_selected(i)==2
            % miniscope
            ind_trace = strcmp(labels_miniscope,cur_label);
            data_selected = data_miniscope(:,ind_trace==1);
            line('XData',t_true_miniscope,'YData',data_selected,'Parent',ax,'Color',g_colors(i,:),'Linewidth',2);
        end
        
        %     if contains(cur_label,'Cell')
        %         ylim1=median(data_selected,'omitnan')-n_iqr*iqr(data_selected);
        %         ylim2=median(data_selected,'omitnan')+n_iqr*iqr(data_selected);
        %     else
        %         ylim1=min(data_selected,[],'omitnan');
        %         ylim2=max(data_selected,[],'omitnan');
        %     end
        
        if strcmp(cur_label,'[vid]Rotation-X')
            line('XData',t_true_tracking,'YData',x_rotation_lp,'Parent',ax,'Color','r','Linewidth',1);
        end
        if strcmp(cur_label,'[vid]Rotation-Y')
            %             line('XData',t_true_tracking,'YData',y_rotation_lp,'Parent',ax,'Color','r','Linewidth',1);
            line('XData',t_true_tracking,'YData',y_rotation_fa,'Parent',ax,'Color','r','Linewidth',1);
        end
        if strcmp(cur_label,'[vid]Rotation-Z')
            line('XData',t_true_tracking,'YData',z_rotation_lp,'Parent',ax,'Color','r','Linewidth',1);
        end
        
        ylim1=min(data_selected,[],'omitnan');
        ylim2=max(data_selected,[],'omitnan');
        ax.YLim = [ylim1,ylim2];
        line('XData',[NaN NaN],'YData',[ylim1 ylim2],'Parent',ax,'Color',[.5 .5 .5],'Tag','Cursor','Linewidth',2);
        
        
        %     set(ax,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
        set(ax,'XTick',[],'XTickLabel',[]);
        % ax.Visible = 'off';
        ax.YLabel.String = cur_label;
        ax.YLabel.Rotation = 0;
        ax.FontSize=10;
        
        % Tracker, trajectory, hd vector
        l_trajectory = line('XData',NaN,'YData',NaN,'Tag','Trajectory','Parent',ax11,...
            'Color',[.25 .25 .25],'LineWidth',.1,'MarkerFaceColor','k');
        l_tracker = line('XData',NaN,'YData',NaN,'Tag','Tracker','Parent',ax11,...
            'Marker','o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r');
        %         l_hdvector = line('XData',NaN,'YData',NaN,'Tag','HDVector','Parent',ax11,...
        %             'Color','r','LineWidth',1);
        %         l_hdvector_tip = line('XData',NaN,'YData',NaN,'Tag','HDVectorTip','Parent',ax11,...
        %             'Marker','^','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r');
        l_hdvector2 = line('XData',NaN,'YData',NaN,'Tag','HDVector','Parent',ax11,...
            'Color','r','LineWidth',1);
        l_hdvector_tip2 = line('XData',NaN,'YData',NaN,'Tag','HDVectorTip','Parent',ax11,...
            'Marker','^','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','r');
        
        set(ax11,'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
        ax11.XLabel.String = 'Position-X';
        ax11.YLabel.String = 'Position-Y';
        
        x_vector = polarplot([NaN NaN],[0 1],'Parent',ax12,'Color',r_colors(1,:),'LineWidth',2);
        x_vector_lp = polarplot([NaN NaN],[0 1],'Parent',ax12,'Color','r','LineWidth',1);
        y_vector = polarplot([NaN NaN],[0 1],'Parent',ax13,'Color',r_colors(2,:),'LineWidth',2);
        y_vector_fa = polarplot([NaN NaN],[0 1],'Parent',ax13,'Color','r','LineWidth',1);
        %         y_vector_lp = polarplot([NaN NaN],[0 1],'Parent',ax13,'Color','r','LineWidth',1);
        z_vector = polarplot([NaN NaN],[0 1],'Parent',ax14,'Color',r_colors(3,:),'LineWidth',2);
        z_vector_lp = polarplot([NaN NaN],[0 1],'Parent',ax14,'Color','r','LineWidth',1);
        
        set(ax12,'RTick',[],'RTickLabel',[],'ThetaDir','clockwise','ThetaZeroLocation','left');
        set(ax13,'RTick',[],'RTickLabel',[],'ThetaDir','clockwise','ThetaZeroLocation','left');
        set(ax14,'RTick',[],'RTickLabel',[],'ThetaDir','clockwise','ThetaZeroLocation','left');
        ax12.Title.String = 'Rotation-X';
        ax13.Title.String = 'Rotation-Y';
        ax14.Title.String = 'Rotation-Z';
        
        %         hold(ax13b,'on')
        %         y_vector2 = polarplot([NaN NaN],[0 1],'Parent',ax13b,'Color','r','LineWidth',1);
        %         y_vector3 = polarplot([NaN NaN],[0 1],'Parent',ax13b,'Color','b','LineWidth',1);
        %         set(ax13b,'RTick',[],'RTickLabel',[],'ThetaDir','clockwise');
    end
    
    
    if flag_save
        work_dir = fullfile(full_path,'processed','sync-frames');
        if isfolder(work_dir)
            rmdir(work_dir,'s');
        end
        mkdir(work_dir);
    end
    
    for i = 1:length(t)
        
        cur_t = t(i);
        t1.String = sprintf('%s (Intan)',tts_t(i,:));
        
        % Video Frame
        [delta_tracking,ind_frame_tracking] = min(abs((t_true_tracking-cur_t)));
        %     cla(ax1);
        %     cla(ax11);
        if delta_tracking < .5/v_tracking.FrameRate
            cur_frame_vid = v_tracking.read(ind_frame_tracking);
            delete(findobj(ax1,'Tag','Frame-tracking'));
            image(cur_frame_vid,'Parent',ax1,'Tag','Frame-tracking');
            ax1.Visible = 'off';
            axis(ax1,'equal');
            drawnow;
            
            t2.String = sprintf('%s (Optitrack)',tts_relative_tracking(ind_frame_tracking,:));
            set(l_trajectory,'XData',x_position(1:ind_frame_tracking),'YData',y_position(1:ind_frame_tracking));
            set(l_tracker,'XData',x_position(ind_frame_tracking),'YData',y_position(ind_frame_tracking));
            
            norm_hdvector = .25;
            %             theta = - y_rotation_lp_rad(ind_frame_tracking) + hd_offset;
            %             set(l_hdvector,'XData',[x_position(ind_frame_tracking),x_position(ind_frame_tracking)+norm_hdvector*cos(theta)],...
            %                 'YData',[y_position(ind_frame_tracking),y_position(ind_frame_tracking)+norm_hdvector*sin(theta)]);
            %             set(l_hdvector_tip,'XData',x_position(ind_frame_tracking)+norm_hdvector*cos(theta),'YData',y_position(ind_frame_tracking)+norm_hdvector*sin(theta));
            
            theta = - y_rotation_fa_rad(ind_frame_tracking) + hd_offset;
            set(l_hdvector2,'XData',[x_position(ind_frame_tracking),x_position(ind_frame_tracking)+norm_hdvector*cos(theta)],...
                'YData',[y_position(ind_frame_tracking),y_position(ind_frame_tracking)+norm_hdvector*sin(theta)]);
            set(l_hdvector_tip2,'XData',x_position(ind_frame_tracking)+norm_hdvector*cos(theta),'YData',y_position(ind_frame_tracking)+norm_hdvector*sin(theta));
            
            
            set(ax11,'XLim',[-1 1],'YLim',[-1 1]);
            set(x_vector,'ThetaData',[x_rotation_rad(ind_frame_tracking),x_rotation_rad(ind_frame_tracking)]);
            set(x_vector_lp,'ThetaData',[x_rotation_lp_rad(ind_frame_tracking),x_rotation_lp_rad(ind_frame_tracking)]);
            
            set(y_vector,'ThetaData',[y_rotation_rad(ind_frame_tracking),y_rotation_rad(ind_frame_tracking)]);
            %             set(y_vector_lp,'ThetaData',[y_rotation_lp_rad(ind_frame_tracking),y_rotation_lp_rad(ind_frame_tracking)]);
            set(y_vector_fa,'ThetaData',[y_rotation_fa_rad(ind_frame_tracking),y_rotation_fa_rad(ind_frame_tracking)]);
            
            set(z_vector,'ThetaData',[z_rotation_rad(ind_frame_tracking),z_rotation_rad(ind_frame_tracking)]);
            set(z_vector_lp,'ThetaData',[z_rotation_lp_rad(ind_frame_tracking),z_rotation_lp_rad(ind_frame_tracking)]);
            
            %             offset_rad = pi/4;
            %             set(y_vector2,'ThetaData',[offset_rad+y_rotation_lp_rad(ind_frame_tracking),pi+y_rotation_lp_rad(ind_frame_tracking)]);
            %             set(y_vector3,'ThetaData',[offset_rad+y_rotation_fa_rad(ind_frame_tracking),pi+y_rotation_fa_rad(ind_frame_tracking)]);
        end
        
        % Miniscope Frame
        [delta_miniscope,ind_frame_miniscope] = min(abs(t_true_miniscope-cur_t));
        %     cla(ax2);
        if delta_miniscope < 1/v_miniscope.FrameRate
            cur_frame_vid = v_miniscope.read(ind_frame_miniscope);
            delete(findobj(ax2,'Tag','Frame-miniscope'));
            image(cur_frame_vid,'Parent',ax2,'Tag','Frame-miniscope');
            ax2.Visible = 'off';
            axis(ax2,'equal');
            drawnow;
            
            t3.String = sprintf('%s (Miniscope)',tts_relative_miniscope(ind_frame_miniscope,:));
        end
        
        % Traces
        for k=1:length(all_axes_traces)
            ax = all_axes_traces(k);
            ax.XLim = [cur_t-time_window,cur_t+time_window];
            ax.YLabel.Units = 'normalized';
            ax.YLabel.Position = [-.05 0 0];
            l = findobj(ax,'Tag','Cursor');
            l.XData = [cur_t cur_t];
        end
        
        if flag_save
            % Saving mode
            pic_name = strcat(sprintf('Frame%05d.jpg',i));
            saveas(f,fullfile(work_dir,pic_name),'jpeg');
            fprintf('Frame Saved [%s].\n',pic_name);
        else
            % Display mode
            %         pause(.1);
        end
    end
    
    if flag_save
        save_video(work_dir,save_dir,'sync-movie',25);
    end
    
    close(f);
end

end

function Info = get_sync_data(full_path)

if nargin == 0
    folder = '/Users/tonio/Desktop';
    recording_name = 'RecExample_RSC';
    full_path = fullfile(folder,recording_name);
end

if ~exist(fullfile(full_path),'dir')
    errordlg(sprintf('Missing folder [%s]',full_path));
    return;
end


% Parameters
num_channels = 2;
sampling_freq = 2000;
% sampling_freq = 2048;
v_thresh_miniscope = 5;
v_thresh_optitrack = 5;


Info = struct('Timing',[],'Metadata',[],...
    'video_tracking',[],'csv_tracking',[],...
    'video_miniscope',[],'csv_miniscope_cells',[],'csv_miniscope_accel',[],...
    'ephys_amplifier',[],'ephys_lowpass',[],'ephys_time',[]);


% Opening analogin.dat
d_analogin = dir(fullfile(full_path,'analogin.dat'));
if isempty(d_analogin)
    errordlg(sprintf('Cannot synchronize: Missing analogin.dat [%s]',full_path));
    return;
else
    
    num_samples = d_analogin.bytes/(num_channels * 2); % uint16 = 2 bytes
    fid = fopen(fullfile(d_analogin.folder,d_analogin.name), 'r');
    v = fread(fid, [num_channels, num_samples], 'uint16');
    fclose(fid);
    v = (v-32768)*0.0003125;
    v_optitrack = v(1,:);
    v_miniscope = [0,.5*abs(diff(v(2,:)))];
    
    % Getting timestamps
    t_ephys = ((0:num_samples-1)/sampling_freq)';
    t_miniscope = t_ephys(v_miniscope>v_thresh_miniscope);
    t_optitrack = t_ephys(v_optitrack>v_thresh_optitrack);
    
    % Completing skipped frames
    
    Info.Timing.t_ephys = t_ephys;
    Info.Timing.duration_ephys = Info.Timing.t_ephys(end)-Info.Timing.t_ephys(1);
    Info.Timing.NumFrames_ephys = length(Info.Timing.t_ephys);
    
    Info.Timing.t_miniscope = t_miniscope;
    Info.Timing.duration_miniscope = Info.Timing.t_miniscope(end)-Info.Timing.t_miniscope(1);
    Info.Timing.NumFrames_miniscope = length(Info.Timing.t_miniscope);
    
    Info.Timing.t_optitrack = t_optitrack;
    Info.Timing.duration_optitrack = Info.Timing.t_optitrack(end)-Info.Timing.t_optitrack(1);
    Info.Timing.NumFrames_optitrack = length(Info.Timing.t_optitrack);
    
    Info.Timing.v_miniscope = v_miniscope;
    Info.Timing.v_optitrack = v_optitrack;
    
    %     figure;
    %     ax1=subplot(411);plot(v_miniscope);ax1.YLabel.String = 'v-miniscope;';
    %     ax2=subplot(412);plot(t_miniscope,[NaN,diff(t_miniscope)]);ax2.YLabel.String = 'diff-t-miniscope;';
    %     ax3=subplot(413);plot(v_optitrack);ax3.YLabel.String = 'v-optitrack;';
    %     ax4=subplot(414);plot(t_optitrack,[NaN,diff(t_optitrack)]);ax4.YLabel.String = 'diff-t_optitrack;';
    %     linkaxes([ax1;ax2;ax3;ax4],'x');
    
end

% Opening video data
d_video = dir(fullfile(full_path,'vid'));
if ~isempty(d_video)
    flag_video = true;
    d_video_avi = dir(fullfile(full_path,'vid','Take*.avi'));
    % Removing hidden files
    d_video_avi = d_video_avi(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_video_avi));
    if ~isempty(d_video_avi)
        Info.video_tracking = d_video_avi;
    end
    d_csv_tracking = dir(fullfile(full_path,'vid','Take*.csv'));
    % Removing hidden files
    d_csv_tracking = d_csv_tracking(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_csv_tracking));
    if ~isempty(d_csv_tracking)
        Info.csv_tracking = d_csv_tracking;
    end
else
    flag_video = false;
end


% Opening miniscope data
d_miniscope = dir(fullfile(full_path,'min'));
if ~isempty(d_miniscope)
    flag_miniscope = true;
    d_video_miniscope = dir(fullfile(full_path,'min','filtered.avi'));
    % Removing hidden files
    d_video_miniscope = d_video_miniscope(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_video_miniscope));
    if ~isempty(d_video_miniscope)
        Info.video_miniscope = d_video_miniscope;
    end
    d_csv_miniscope1 = dir(fullfile(full_path,'min','headOrientation.csv'));
    % Removing hidden files
    d_csv_miniscope1 = d_csv_miniscope1(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_csv_miniscope1));
    if ~isempty(d_csv_miniscope1)
        Info.csv_miniscope_accel = d_csv_miniscope1;
    end
    d_csv_miniscope2 = dir(fullfile(full_path,'min','*_C.csv'));
    % Removing hidden files
    d_csv_miniscope2 = d_csv_miniscope2(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_csv_miniscope2));
    if ~isempty(d_csv_miniscope2)
        Info.csv_miniscope_cells = d_csv_miniscope2;
    end
    
else
    flag_miniscope = false;
end

% Opening ephys data
d_ephys = dir(fullfile(full_path,'ephys'));
if ~isempty(d_ephys)
    flag_ephys = true;
    d_ephys_amplifier = dir(fullfile(full_path,'ephys','amplifier.dat'));
    % Removing hidden files
    d_ephys_amplifier = d_ephys_amplifier(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_ephys_amplifier));
    if ~isempty(d_ephys_amplifier)
        Info.ephys_amplifier = d_ephys_amplifier;
    end
    d_ephys_lowpass = dir(fullfile(full_path,'ephys','lowpass.dat'));
    % Removing hidden files
    d_ephys_lowpass = d_ephys_lowpass(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_ephys_lowpass));
    if ~isempty(d_ephys_lowpass)
        Info.ephys_lowpass = d_ephys_lowpass;
    end
    d_ephys_time = dir(fullfile(full_path,'ephys','time.dat'));
    % Removing hidden files
    d_ephys_time = d_ephys_time(arrayfun(@(x) ~strcmp(x.name(1),'.'),d_ephys_time));
    if ~isempty(d_ephys_time)
        Info.ephys_time = d_ephys_time;
    end
else
    flag_ephys = false;
end

Metadata.flag_video = flag_video;
Metadata.flag_miniscope = flag_miniscope;
Metadata.flag_ephys = flag_ephys;
Metadata.full_path = full_path;
temp = regexp(full_path,filesep,'split');
Metadata.recording_name = char(temp(end));
Metadata.folder = strrep(full_path,strcat(filesep,char(temp(end))),'');
Metadata.num_channels = num_channels;
Metadata.sampling_freq = sampling_freq;
Metadata.v_thresh_miniscope = v_thresh_miniscope;
Metadata.v_thresh_optitrack = v_thresh_optitrack;
Info.Metadata = Metadata;

end

function g_colors = get_colors(n_colors,cmap)
% Gives interpolated color list from specified colormap

if nargin < 2
    cmap = 'jet';
end
if nargin < 1
    n_colors = 64;
end

% g_colors = get(groot,'DefaultAxesColorOrder');
switch cmap
    case 'jet'
        f_colors = [0         0    0.5625;
            0         0    0.6250;
            0         0    0.6875;
            0         0    0.7500;
            0         0    0.8125;
            0         0    0.8750;
            0         0    0.9375;
            0         0    1.0000;
            0    0.0625    1.0000;
            0    0.1250    1.0000;
            0    0.1875    1.0000;
            0    0.2500    1.0000;
            0    0.3125    1.0000;
            0    0.3750    1.0000;
            0    0.4375    1.0000;
            0    0.5000    1.0000;
            0    0.5625    1.0000;
            0    0.6250    1.0000;
            0    0.6875    1.0000;
            0    0.7500    1.0000;
            0    0.8125    1.0000;
            0    0.8750    1.0000;
            0    0.9375    1.0000;
            0    1.0000    1.0000;
            0.0625    1.0000    0.9375;
            0.1250    1.0000    0.8750;
            0.1875    1.0000    0.8125;
            0.2500    1.0000    0.7500;
            0.3125    1.0000    0.6875;
            0.3750    1.0000    0.6250;
            0.4375    1.0000    0.5625;
            0.5000    1.0000    0.5000;
            0.5625    1.0000    0.4375;
            0.6250    1.0000    0.3750;
            0.6875    1.0000    0.3125;
            0.7500    1.0000    0.2500;
            0.8125    1.0000    0.1875;
            0.8750    1.0000    0.1250;
            0.9375    1.0000    0.0625;
            1.0000    1.0000         0;
            1.0000    0.9375         0;
            1.0000    0.8750         0;
            1.0000    0.8125         0;
            1.0000    0.7500         0;
            1.0000    0.6875         0;
            1.0000    0.6250         0;
            1.0000    0.5625         0;
            1.0000    0.5000         0;
            1.0000    0.4375         0;
            1.0000    0.3750         0;
            1.0000    0.3125         0;
            1.0000    0.2500         0;
            1.0000    0.1875         0;
            1.0000    0.1250         0;
            1.0000    0.0625         0;
            1.0000         0         0;
            0.9375         0         0;
            0.8750         0         0;
            0.8125         0         0;
            0.7500         0         0;
            0.6875         0         0;
            0.6250         0         0;
            0.5625         0         0;
            0.5000         0         0];
        
    case 'hot'
        f_colors = [0.0417         0         0;
            0.0833         0         0;
            0.1250         0         0;
            0.1667         0         0;
            0.2083         0         0;
            0.2500         0         0;
            0.2917         0         0;
            0.3333         0         0;
            0.3750         0         0;
            0.4167         0         0;
            0.4583         0         0;
            0.5000         0         0;
            0.5417         0         0;
            0.5833         0         0;
            0.6250         0         0;
            0.6667         0         0;
            0.7083         0         0;
            0.7500         0         0;
            0.7917         0         0;
            0.8333         0         0;
            0.8750         0         0;
            0.9167         0         0;
            0.9583         0         0;
            1.0000         0         0;
            1.0000    0.0417         0;
            1.0000    0.0833         0;
            1.0000    0.1250         0;
            1.0000    0.1667         0;
            1.0000    0.2083         0;
            1.0000    0.2500         0;
            1.0000    0.2917         0;
            1.0000    0.3333         0;
            1.0000    0.3750         0;
            1.0000    0.4167         0;
            1.0000    0.4583         0;
            1.0000    0.5000         0;
            1.0000    0.5417         0;
            1.0000    0.5833         0;
            1.0000    0.6250         0;
            1.0000    0.6667         0;
            1.0000    0.7083         0;
            1.0000    0.7500         0;
            1.0000    0.7917         0;
            1.0000    0.8333         0;
            1.0000    0.8750         0;
            1.0000    0.9167         0;
            1.0000    0.9583         0;
            1.0000    1.0000         0;
            1.0000    1.0000    0.0625;
            1.0000    1.0000    0.1250;
            1.0000    1.0000    0.1875;
            1.0000    1.0000    0.2500;
            1.0000    1.0000    0.3125;
            1.0000    1.0000    0.3750;
            1.0000    1.0000    0.4375;
            1.0000    1.0000    0.5000;
            1.0000    1.0000    0.5625;
            1.0000    1.0000    0.6250;
            1.0000    1.0000    0.6875;
            1.0000    1.0000    0.7500;
            1.0000    1.0000    0.8125;
            1.0000    1.0000    0.8750;
            1.0000    1.0000    0.9375;
            1.0000    1.0000    1.0000];
        
    case 'parula'
        f_colors = [0.2422    0.1504    0.6603;
            0.2504    0.1650    0.7076;
            0.2578    0.1818    0.7511;
            0.2647    0.1978    0.7952;
            0.2706    0.2147    0.8364;
            0.2751    0.2342    0.8710;
            0.2783    0.2559    0.8991;
            0.2803    0.2782    0.9221;
            0.2813    0.3006    0.9414;
            0.2810    0.3228    0.9579;
            0.2795    0.3447    0.9717;
            0.2760    0.3667    0.9829;
            0.2699    0.3892    0.9906;
            0.2602    0.4123    0.9952;
            0.2440    0.4358    0.9988
            0.2206    0.4603    0.9973;
            0.1963    0.4847    0.9892;
            0.1834    0.5074    0.9798;
            0.1786    0.5289    0.9682;
            0.1764    0.5499    0.9520;
            0.1687    0.5703    0.9359;
            0.1540    0.5902    0.9218;
            0.1460    0.6091    0.9079;
            0.1380    0.6276    0.8973;
            0.1248    0.6459    0.8883;
            0.1113    0.6635    0.8763;
            0.0952    0.6798    0.8598;
            0.0689    0.6948    0.8394;
            0.0297    0.7082    0.8163;
            0.0036    0.7203    0.7917;
            0.0067    0.7312    0.7660;
            0.0433    0.7411    0.7394;
            0.0964    0.7500    0.7120;
            0.1408    0.7584    0.6842;
            0.1717    0.7670    0.6554;
            0.1938    0.7758    0.6251;
            0.2161    0.7843    0.5923;
            0.2470    0.7918    0.5567;
            0.2906    0.7973    0.5188;
            0.3406    0.8008    0.4789;
            0.3909    0.8029    0.4354;
            0.4456    0.8024    0.3909;
            0.5044    0.7993    0.3480;
            0.5616    0.7942    0.3045;
            0.6174    0.7876    0.2612;
            0.6720    0.7793    0.2227;
            0.7242    0.7698    0.1910;
            0.7738    0.7598    0.1646;
            0.8203    0.7498    0.1535;
            0.8634    0.7406    0.1596;
            0.9035    0.7330    0.1774;
            0.9393    0.7288    0.2100;
            0.9728    0.7298    0.2394;
            0.9956    0.7434    0.2371;
            0.9970    0.7659    0.2199;
            0.9952    0.7893    0.2028;
            0.9892    0.8136    0.1885;
            0.9786    0.8386    0.1766;
            0.9676    0.8639    0.1643;
            0.9610    0.8890    0.1537;
            0.9597    0.9135    0.1423;
            0.9628    0.9373    0.1265;
            0.9691    0.9606    0.1064;
            0.9769    0.9839    0.0805];
        
    otherwise
        f_colors = [0         0         0;
            0.0159    0.0159    0.0159;
            0.0317    0.0317    0.0317;
            0.0476    0.0476    0.0476;
            0.0635    0.0635    0.0635;
            0.0794    0.0794    0.0794;
            0.0952    0.0952    0.0952;
            0.1111    0.1111    0.1111;
            0.1270    0.1270    0.1270;
            0.1429    0.1429    0.1429;
            0.1587    0.1587    0.1587;
            0.1746    0.1746    0.1746;
            0.1905    0.1905    0.1905;
            0.2063    0.2063    0.2063;
            0.2222    0.2222    0.2222;
            0.2381    0.2381    0.2381;
            0.2540    0.2540    0.2540;
            0.2698    0.2698    0.2698;
            0.2857    0.2857    0.2857;
            0.3016    0.3016    0.3016;
            0.3175    0.3175    0.3175;
            0.3333    0.3333    0.3333;
            0.3492    0.3492    0.3492;
            0.3651    0.3651    0.3651;
            0.3810    0.3810    0.3810;
            0.3968    0.3968    0.3968;
            0.4127    0.4127    0.4127;
            0.4286    0.4286    0.4286;
            0.4444    0.4444    0.4444;
            0.4603    0.4603    0.4603;
            0.4762    0.4762    0.4762;
            0.4921    0.4921    0.4921;
            0.5079    0.5079    0.5079;
            0.5238    0.5238    0.5238;
            0.5397    0.5397    0.5397;
            0.5556    0.5556    0.5556;
            0.5714    0.5714    0.5714;
            0.5873    0.5873    0.5873;
            0.6032    0.6032    0.6032;
            0.6190    0.6190    0.6190;
            0.6349    0.6349    0.6349;
            0.6508    0.6508    0.6508;
            0.6667    0.6667    0.6667;
            0.6825    0.6825    0.6825;
            0.6984    0.6984    0.6984;
            0.7143    0.7143    0.7143;
            0.7302    0.7302    0.7302;
            0.7460    0.7460    0.7460;
            0.7619    0.7619    0.7619;
            0.7778    0.7778    0.7778;
            0.7937    0.7937    0.7937;
            0.8095    0.8095    0.8095;
            0.8254    0.8254    0.8254;
            0.8413    0.8413    0.8413;
            0.8571    0.8571    0.8571;
            0.8730    0.8730    0.8730;
            0.8889    0.8889    0.8889;
            0.9048    0.9048    0.9048;
            0.9206    0.9206    0.9206;
            0.9365    0.9365    0.9365;
            0.9524    0.9524    0.9524;
            0.9683    0.9683    0.9683;
            0.9841    0.9841    0.9841;
            1.0000    1.0000    1.0000];
        
end

g_colors = (interp1(1:length(f_colors),f_colors,round(rescale(1:n_colors,1,length(f_colors)))));

end

function save_video(workingDir,savedir,video_name,video_quality)

if nargin < 4
    video_quality = 100;
end

yourfolder = dir(fullfile(workingDir,'*.jpg'));
% Removing hidden files
yourfolder = yourfolder(arrayfun(@(x) ~strcmp(x.name(1),'.'),yourfolder));

xlsfiles = {yourfolder.name};
[~,idx] = sort(xlsfiles);
new_folder = yourfolder(idx);
imageNames = {new_folder.name}';
 
try
    outputVideo = VideoWriter(fullfile(savedir,strcat(video_name,'.mp4')),'MPEG-4');
catch
    outputVideo = VideoWriter(fullfile(savedir,strcat(video_name,'.avi')),'Motion JPEG AVI');
end
outputVideo.FrameRate = 20;
outputVideo.Quality = video_quality;
% outputVideo.LosslessCompression = video_compression;
open(outputVideo);

% Writing Video + waitbar
h = waitbar(0,'Writing video file: 0.0 % completed.');
for ii = 1:length(imageNames)
    img = imread(fullfile(workingDir,imageNames{ii}));
    writeVideo(outputVideo,img);
    waitbar(ii/length(imageNames),h,sprintf('Writing video file: %.1f %% completed.',100*ii/length(imageNames)));
end
close(h);

fprintf('Video Saved [%s].\n',fullfile(savedir,video_name));
fprintf('[Working Directory: %s]\n',workingDir);
close(outputVideo);

end

function pos = get_position(n_rows,n_columns,counter,margins)
% Gives precise axes positions in multiple subplots

if nargin<4
    w_margin_1 = .05; % left margin
    w_margin_2 = .05; % right margin
    w_eps = .01;      % horizontal spacing
    h_margin_1 = .05; % bottom margin
    h_margin_2 = .05; % top margin
    h_eps = .01;      % vertical spacing
    
    margins = [w_margin_1,w_margin_2,w_eps;
        h_margin_1,h_margin_2,h_eps];
end

w_margin_1 = margins(1,1); % left margin
w_margin_2 = margins(1,2); % right margin
w_eps = margins(1,3);      % horizontal spacing
h_margin_1 = margins(2,1); % bottom margin
h_margin_2 = margins(2,2); % top margin
h_eps = margins(2,3);

pos1 = w_margin_1 + (mod(counter-1,n_columns)/n_columns)*(1-(w_margin_1+w_margin_2));
pos2 = 1 - h_margin_2 - (ceil(counter/n_columns)/n_rows)*(1-(h_margin_1+h_margin_2));
pos3 = ((1-(w_margin_1+w_margin_2))/n_columns) - w_eps;
pos4 = ((1-(h_margin_1+h_margin_2))/n_rows) - h_eps;
pos=[pos1,pos2,pos3,pos4];

end

function Y_filtered = correct_angle(Y,tq)
% low pass filters angle

thresh_angle = 5;
ind_remove = abs([0,diff(Y)])>thresh_angle;

Y_filtered=Y;
t_nn = tq(ind_remove==0);
y_nn = Y_filtered(ind_remove==0);
Y_filtered(ind_remove==1)=NaN;
% Y_filtered = interp1(t_nn,y_nn,tq);

end

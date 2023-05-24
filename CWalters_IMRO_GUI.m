% Custom Electrode Selection GUI for type 24 Neuropixels Probes
% Charlie Walters 12/4/22
% Parts of this code taken from SpikeGLX's IMRO scripts
% IMRO entries: (channel shank bank reference electrode)

% User inputs
ref_electrode = 0; % 0 is external, 1-4 are the tips
name_string = 'custom'; % will add the probe type and the ref electrode automatically
top_row = 250; % top row of sites you want to select from
bottom_row = 200; % bottom row of sites you want to select from
load_imro = 0; % if you set to 1, you will be asked to select the file in a ui

% Warnings
if ref_electrode>4
    error('Ref electrode out of bounds.')
    % Technically I think you can still pass the on-shank references to
    % SpikeGLX but it's not recommended.
end
if (top_row-bottom_row+1)<48
    error('Not enough rows displayed.')
end

% Fixed params
nameStr = sprintf('NPtype24_%s_ref%d', name_string, ref_electrode);
top_site = top_row*2+1;
bottom_site = bottom_row*2;
probe_type = 24;
n_channels = 384;
n_shanks = 4;
n_electrodes_per_shank = 1280;

% Calculate all the electrode,shank pairings that belong to a particular channel
channel_map = zeros(n_shanks,n_electrodes_per_shank);
bank_map = zeros(n_shanks,n_electrodes_per_shank);
for shankInd = 1:n_shanks
    for elecInd = 1:n_electrodes_per_shank
        [bank_map(shankInd,elecInd), channel_map(shankInd,elecInd)] = ElecToChan(shankInd-1,elecInd-1);
    end
end

% Re-organize channel map
channel_lookup = cell(1,n_channels);
for channelInd = 1:n_channels
    [shanks,electrodes] = find(channel_map==(channelInd-1));
    shanks = shanks-1; % find function produces indices that are 1-indexed
    electrodes = electrodes-1;
    channel_lookup{channelInd} = [shanks electrodes]; % either 13 or 14 electrodes per channel
end

active_electrodes = zeros(n_shanks,n_electrodes_per_shank);


%% Load IMRO table

if load_imro
    filename = uigetfile('*.imro');
    opts = delimitedTextImportOptions('Delimiter',{'(',')'},...
        'ConsecutiveDelimitersRule','join',...
        'LeadingDelimitersRule','ignore',...
        'TrailingDelimitersRule','ignore');
    imro = readtable(filename,opts);
    imro = imro(1,2:end);
    imro = table2array(imro);
    for entry = 1:n_channels
        data = imro{entry}(length(char(string(entry-1)))+2:end);
        shnk = str2num(data(1));
        chan = str2num(data(7:end));
        active_electrodes(shnk+1,chan+1) = 1; 
    end
    clear entry imro opts
    bottom_imro_site = find(sum(active_electrodes,1),1,'first');
    top_imro_site = find(sum(active_electrodes,1),1,'last');
    if bottom_imro_site < bottom_site || top_imro_site > top_site
        if rem(top_imro_site,2)==0
            top_imro_site = top_imro_site + 1;
        end
        if rem(bottom_imro_site,2)~=0
            bottom_imro_site = bottom_imro_site - 1;
        end

        top_site = top_imro_site;
        bottom_site = bottom_imro_site;
        
        top_row = (top_site-1)/2;
        bottom_row = (bottom_site)/2;
        warning('displaying more rows than entered')
    end
end


%% Initialize GUI

active_color = 'b';
n_disp_rows = top_row - bottom_row + 1;
n_buttons = n_disp_rows * n_shanks * 2;

figure("Name","GUI","Position",[0,0,400,800])
axis([0 13 0 top_row+2])
ylim([bottom_row-2, top_row+2])
yticks(bottom_row:2:top_row)
xlim([0,18])
xticks([])

% Make the probe/sites layout
row_height = 1;
shank_width = 2;
col_width = 1;
button_coords = zeros(n_shanks,n_electrodes_per_shank,4); % 4 coordinates per rectangle
x_position = 1;
rectangles = gobjects(n_shanks,n_electrodes_per_shank);

for iShank = 1:4
    for iSite = bottom_site+1:top_site+1
        row = floor((iSite-1)/2);
        if rem(iSite,2)==0
           col = 1;
        else
           col = 0;
        end
        % calculate the position coordinates
        button_coords(iShank,iSite,:) = [x_position+col,row-0.5,col_width,row_height];
        % plot the rectangle and store it in the array rectangles
        rectangles(iShank,iSite) = rectangle('Position',button_coords(iShank,iSite,:));
        if active_electrodes(iShank,iSite) == 1
            rectangles(iShank,iSite).FaceColor = active_color;
        end
    end
        x_position = x_position+3; % spacing between shanks
end
clear iShank iSite

% Make the 'finished selecting' button
finished_coords = [x_position,bottom_row,4,3];
finished_rectanlge = rectangle('Position',finished_coords,'FaceColor','g');
text(finished_coords(1)+0.1,finished_coords(2)+(finished_coords(4)/2),'print table')
channel_display = text(finished_coords(1),finished_coords(2)+5,string(sum(active_electrodes,'all'))+' sites selected');


%% While Loop

finished_selecting = 0;

while finished_selecting==0

    [clickx,clicky] = ginput(1);
    csite = 0;
    for iShank = 1:4
        for iSite = bottom_site+1:top_site+1
            if clickx>button_coords(iShank,iSite,1) && clickx<button_coords(iShank,iSite,1)+col_width && ...
                    clicky>button_coords(iShank,iSite,2) && clicky<button_coords(iShank,iSite,2)+row_height
                csite = [iShank,iSite];
                csite_truename = [iShank,iSite]-1;
            end
        end
    end
    clear iShank iSite

    if csite==0
        no_site_clicked = 1;
    else % the click was within a rectangle
        cchannel = channel_map(csite(1),csite(2)); % look up the channel the site corresponds to
        other_sites = channel_lookup{cchannel+1};
        if active_electrodes(csite(1),csite(2))==0 % if the site was not already selected
            for iOsite = 1:length(other_sites)
                % check whether other site is within bottom_site:top_site
                if other_sites(iOsite,2) >= bottom_site && other_sites(iOsite,2) <= top_site
                    % set other site to false
                    active_electrodes(other_sites(iOsite,1)+1,other_sites(iOsite,2)+1) = 0;
                    % de-highlight other site
                    rectangles(other_sites(iOsite,1)+1,other_sites(iOsite,2)+1).FaceColor = [0.8 0.8 0.8];
                end
            end
            clear iOsite

            % set clicked site to true
            active_electrodes(csite(1),csite(2)) = 1;
            % highlight clicked site
            rectangles(csite(1),csite(2)).FaceColor = active_color;
        else % if the site was already selected
            % set clicked site to false
            active_electrodes(csite(1),csite(2)) = 0;
            for iOsite = 1:length(other_sites)
                % check whether other site is within bottom_site:top_site
                if other_sites(iOsite,2) >= bottom_site && other_sites(iOsite,2) <= top_site
                    rectangles(other_sites(iOsite,1)+1,other_sites(iOsite,2)+1).FaceColor = [1 1 1];
                end
            end
            clear iOsite
        end
        channel_display.Color = 'w';
        channel_display = text(finished_coords(1),finished_coords(2)+5,string(sum(active_electrodes,'all'))+' sites selected');

    end

    % finished selecting
    if clickx>finished_coords(1) && clickx<finished_coords(1)+finished_coords(3) && ...
                    clicky>finished_coords(2) && clicky<finished_coords(2)+finished_coords(4)
        if sum(active_electrodes,'all')~=384
            warning('not enough sites selected')
        else
            finished_selecting = 1;
        end
    end        


end


%% Find info needed for IMRO table

imrows = nan(n_channels,5);
[imrows(:,2),imrows(:,5)] = find(active_electrodes);
imrows(:,2) = imrows(:,2)-1;
imrows(:,5) = imrows(:,5)-1;

for channel = 1:n_channels
    imrows(channel,1) = channel_map(imrows(channel,2)+1,imrows(channel,5)+1);
    imrows(channel,3) = bank_map(imrows(channel,2)+1,imrows(channel,5)+1);
end
clear channel

imrows = sortrows(imrows,1);
imrows(:,4) = ref_electrode;


%% Write the IMRO table

% open a new file wherever we are
fileName = [nameStr,'.imro'];
nmID = fopen(fileName,'w');

% imro table
% print first entry, specifying probe type and number of channels
fprintf(nmID,'(%d,%d)', 24, 384);
for channel = 1:n_channels
    fprintf(nmID,'(%d %d %d %d %d)', imrows(channel,1), imrows(channel,2), imrows(channel,3), imrows(channel,4), imrows(channel,5) );
end
clear channel
fprintf(nmID, '\n');
fclose(nmID);


%% Functions

function [bank, chan] = ElecToChan( shank, elecInd )

    % electrode to channel map
    elecMap = zeros(4,8,'single');
    elecMap(1,:) = [0,2,4,6,5,7,1,3];
    elecMap(2,:) = [1,3,5,7,4,6,0,2];
    elecMap(3,:) = [4,6,0,2,1,3,5,7];
    elecMap(4,:) = [5,7,1,3,0,2,4,6];
    
    bank = floor(elecInd/384);
    
    % which block within the bank?
    blockID = floor((elecInd - bank*384)/48);
    
    % which channel within the block?
    subBlockInd = mod((elecInd - bank*384), 48);
    
    % get the channel number for that shank, bank, and block combo
    chan = 48*elecMap(shank+1, blockID+1) + subBlockInd;

end
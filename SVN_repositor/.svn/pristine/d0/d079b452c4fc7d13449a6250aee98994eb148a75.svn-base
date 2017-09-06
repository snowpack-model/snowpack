function display_smet(fn)
% display_smet(filename)
% function to display smet meteo data
% input is filename
%
% optional requirement:
% 'datetickzoom' to make it pretty:
%       http://www.mathworks.com/matlabcentral/fileexchange/15029-datetickzoom-automatically-update-dateticks/content/datetickzoom.m
%
% 'suptitle' to have a title with filename:
%       http://www.mathworks.com/matlabcentral/fileexchange/45049-bland-altman-and-correlation-plot/content/BlandAltman/suptitle.m
%
% Author: A. Koehler, koehler@slf.ch
% Date: 30-01-2015

data = import_smet(fn);

if isempty(data)
    warning('struct is empty!');
    return;
end

% build figure with subplots
h = figure;
col = 2;
row = ceil(length(data.sensors)/col);

for nr = 1:length(data.sensors)
    ha(nr) = subplot(row, col, nr);
    exclude = data.(data.sensors{nr})~=data.nodata;
    hl(nr) = plot(ha(nr), data.time(exclude), data.(data.sensors{nr})(exclude));
    % color excludes red!
    if ~all(exclude)
        hold on;
        gaps = diff(exclude);
        gap_idx = find(diff(exclude));
        if gaps(gap_idx(1)) == 1 % missing values at start
            hle = plot(data.time([1 gap_idx(1)+1]), [data.(data.sensors{nr})(gap_idx(1)+1) data.(data.sensors{nr})(gap_idx(1)+1)], 'r');
            set(hle, 'LineWidth', 0.1);
            gap_idx = gap_idx(2:end);
        end
        if  gaps(gap_idx(end)) == -1 % missing values at end
            hle = plot(data.time([gap_idx(end) length(data.time)]), [data.(data.sensors{nr})(gap_idx(end)) data.(data.sensors{nr})(gap_idx(end))], 'r');
            set(hle, 'LineWidth', 0.1);
            gap_idx = gap_idx(1:end-1);
        end
        % now plotting the rest of the values in red
        if mod(length(gap_idx),2)==0 & ~isempty(gap_idx)
            for gap_nr = 1:2:length(gap_idx)
                idx = [ gap_idx([gap_nr]) gap_idx([gap_nr+1])+1 ];
                hle = plot(data.time( idx ), data.(data.sensors{nr})(idx) , 'r');
                set(hle, 'LineWidth', 0.1);
                hll = line([mean(data.time( idx ))*[1 1]], [min(data.(data.sensors{nr})(exclude)) max(data.(data.sensors{nr})(exclude))], 'LineWidth', 0.1, 'Color', 'r');
            end
        else
            warning('cannot color missing values... something is wrong!')
        end
    end
    
    set(hl(nr), 'LineWidth', 0.1);
    ylabel(data.sensors{nr});
    xlabel('time')
    if exist('datetickzoom')==2
        datetickzoom('x', 'mm/dd HH', 'keepticks')
    else
        datetick('x', 'mm/dd HH', 'keepticks')
    end
    axis tight
end
linkaxes(ha, 'x');
if exist('suptitle')==2
    suptitle(fn);
end
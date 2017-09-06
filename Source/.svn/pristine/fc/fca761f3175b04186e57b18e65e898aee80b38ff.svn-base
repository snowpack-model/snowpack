function data = import_smet(fn)
% data = import_smet(filename)
% function to import a smet weather station, give the filename as input
% return is struct with weather data
%
% Author: A. Koehler, koehler@slf.ch
% Date: 30-01-2015

data = struct([]);

if exist(fn, 'file') ~= 2
    warning(sprintf('file %s not found', fn))
    return;
end

% write header and find number header lines
fnid = fopen(fn);
headerline = fgetl(fnid);
lnr = 0;

while ischar(headerline)
    lnr = lnr + 1 ; % headerlines
    headerfield = regexp(headerline, '=.');
    if strcmpi(headerline,'[DATA]') % skip when data block comes
        break;
    end
    if headerfield
        fieldname = sscanf(headerline(1:headerfield-1), '%s%*c'); % ommit tailing spaces!
        fieldvalue = headerline(headerfield+2:end);
        if ~isempty(str2num(fieldvalue))
            fieldvalue = str2num(fieldvalue);
        end
        
        if isempty(data)
            data = struct(fieldname, fieldvalue);
        else
            data.(fieldname) = fieldvalue;
        end
    end   
    headerline = fgetl(fnid);
end
fclose(fnid);

meteo = importdata(fn, ' ', lnr);

% read timestamp
data.time = datenum(meteo.textdata(lnr + 1:end,1),'yyyy-mm-ddTHH:MM');

% assign data fields
sensors = textscan(data.fields, '%s'); % get sensors
sensors = sensors{1,1}(2:end); % remove timestamp

data.sensors = sensors;

for i = 1:length(sensors)
   data.(sensors{i}) = meteo.data(:,i);
   if all(size(data.(sensors{i})) ~= size(data.time))
       warning('data inconsistency!')
       data = struct([]);
       return;
   end
end

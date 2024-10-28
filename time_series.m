clear
close all

% generate the time series


read_dir = fullfile(pwd,'pics/raw');
list = dir(fullfile(read_dir,'*.jpg'));

expression = '(?<week>\D+)_(?<month>\D+)_(?<date>\d+)_(?<hour>\d+)-(?<minute>\d+)-(?<second>\d+)_(?<year>\d+).jpg';
timeInfo = struct('week', {}, 'month', {}, 'date', {}, 'hour', {}, 'minute', {}, 'second', {}, 'year', {});
time_matrix = zeros(length(list),6);
% year % month % date
% hour % minute % second
months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
expression2 = '(?<week>\D+)_(?<month>\D+)__(?<date>\d+)_(?<hour>\d+)-(?<minute>\d+)-(?<second>\d+)_(?<year>\d+).jpg';
timeInfo2 = struct('week', {}, 'month', {}, 'date', {}, 'hour', {}, 'minute', {}, 'second', {}, 'year', {});


for i = 1:length(list)
    try
        timeInfo(i) = regexp(list(i).name,expression,'names');
        timeInfo2(i) = regexp(list(i).name,expression2,'names');
    catch
    end
end


for i = 1:length(list)
    time_matrix(i,1) = str2double(timeInfo(i).year);
    for j = 1:12
        if strcmp(timeInfo(i).month, months{j})
            time_matrix(i,2) = j;
        end
    end
    if time_matrix(i,2) == 0
        for j = 1:12
            if strcmp(timeInfo2(i).month, months{j})
                time_matrix(i,2) = j;
            end
        end
    end
    time_matrix(i,3) = str2double(timeInfo(i).date);
    time_matrix(i,4) = str2double(timeInfo(i).hour);
    time_matrix(i,5) = str2double(timeInfo(i).minute);
    time_matrix(i,6) = str2double(timeInfo(i).second);
end

time_long = datetime(time_matrix);

[time_ascend, time_index] = sort(time_long,'ascend');
% time index helps we find the order from time_long
% i.e. 1st element from time_ascend can be found from time_long(time_index(1))

save('time_series.mat','time_long','time_matrix','time_ascend','time_index','list')

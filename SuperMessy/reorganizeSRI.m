%% Same as regular reorganize but added 'SRI' at end so won't get mixed up with other code


%The purpose of this Function is to reorganize the parcels 
%in the passed in array based on the Parc.key.
%%I you want to reorder the parcel in your own way you just need to follow
%these rules:
%   create a csv 
%   with numbers in the fisrt column that represent the new
%   order of the parcels
%   the first row should be a hearder with the labels to the rows
%

function newArray = reorganizeSRI(Array, ParcInfo)
    fid = fopen(fullfile(ParcInfo.key), 'r'); % open the parcel Key
    key = textscan(fid, repmat('%s',1,10), 'delimiter',',', 'CollectOutput',true); % scan place into a cell
    order = key{1}; %pull first colum
    fclose(fid);%close the csv
    newArray = zeros(ParcInfo.NP, size(Array, 2)); %build an empty array with the size of the the old array
    for i = 2:length(order)
        newArray(i-1, :) = Array(str2num(order{i}), :);
    end
   
end
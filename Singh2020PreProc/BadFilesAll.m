function[Bad]=BadFilesAll(Dir,MinSize)
tCD=cd;
%% Removes Server Jobs with output less than MinSize
%if Dir(end)~='\'
%    Dir=strcat(Dir,'\');
%end
G=GetFileNames(Dir);
Bad=[];
for i=1:numel(G)
    file=dir(fullfile(Dir,G{i}));
    fullfile(Dir,G{i})
    if isempty(file)||file.bytes<MinSize
        Bad=[Bad G(i)]; %#ok<*AGROW>
        %       strcat(Dir,G{i})
        delete(fullfile(Dir,G{i}));
    end
end
cd(tCD);
end

function [filenames] = GetFileNames(Inpath)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    list=dir(Inpath);  %get info of files/folders in current directory
    filenames={list(~[list.isdir]).name}; %Get only files, not folders
end


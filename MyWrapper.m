function Subjlist = MyWrapper(nSub, out_dir, in_dir, exclude_sublist)
%MYWRAPPER A simple wrapper to MyStartServerDT
%   This function compare the folders in in_dir and out_dir and
%   process |nSub| unprocessed subjects in ascending order of ID.
%   (Note: ascending order of STRINGS, not necessarily the same
%   as the order of numbers, though true for HCP).
%
%   We add an optional input |exclude_sublist| (path to a file)
%   to allow excluding some subjects. |exclude_sublist| should
%   be a text file with one subject ID in each line. Here we
%   simply use it to exclude the subjects already preprocessed
%   in (Singh et. al., 2020). You may want to remove this part
%   altogether.
%
%   Usually, each subject needs ~18min on CCPLINUX1.

% Compatiability issue
assert(exist('readlines', 'builtin'), ...
    'This script uses the function readlines(). Please modify the code or use R2020b+!')

% Add scripts to path
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'Singh2020PreProc'))

% Inputs
if nargin < 4
    exclude_sublist = fullfile(script_dir, 'old_subjlist.txt');
end
if nargin < 3 || isempty(in_dir)
    in_dir = '/net/10.27.136.121/hcpdb/packages/unzip/HCP_1200';
end
if nargin < 2 || isempty(out_dir)
    out_dir = fullfile(script_dir, 'GSR3');
end
if nargin < 1 || isempty(nSub)
    nSub = 10;
end

% Find unprocessed subjects
try
    allSub = dir(in_dir);
    allSub = {allSub(3:end).name};
catch
    allSub = {};
end
if isempty(allSub)  % Not sure why it couldn't be catched, but we do need this
    error('No data found in input dir! Could be file permission error.');
end
try
    oldSub = dir(fullfile(out_dir, 'Out'));
    oldSub = {oldSub(3:end).name};
catch  % Output dir has not been created yet or contains no subject
    oldSub = {};
end
if isfile(exclude_sublist)
    try
        exSub = cellstr(readlines(exclude_sublist));
    catch
        warning([exclude_sublist ' is not formatted correctly as a subject list!'])
        exSub = {};
    end
    oldSub = union(oldSub, exSub);
end
newSub = setdiff(allSub, oldSub);

% Processing
nSub = min(nSub, numel(newSub));
Subjlist = newSub(1:nSub)
t1 = tic;
MyStartServerDT(Subjlist, 3, 'y', out_dir, in_dir);
disp(['Elapsed time for all ' num2str(nSub) ' subjects:'])
toc(t1);

end
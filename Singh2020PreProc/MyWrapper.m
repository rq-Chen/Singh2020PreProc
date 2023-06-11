function Subjlist = MyWrapper(nSub, out_dir, in_dir)
%MYWRAPPER A simple wrapper to MyStartServerDT
%   This function compare the folders in in_dir and out_dir and
%   process |nSub| unprocessed subjects in ascending order of ID.
%   (Note: ascending order of STRINGS, not necessarily the same
%   as the order of numbers, though true for HCP).
%
%   Usually, each subject needs ~18min on CCPLINUX1.

% Inputs
if nargin < 3 || isempty(in_dir)
    in_dir = '/net/10.20.145.162/HCPpackages03/unzip/1200subject';
end
if nargin < 2 || isempty(out_dir)
    script_dir = fileparts(mfilename('fullpath'));
    out_dir = fullfile(script_dir, '..', 'GSR3');
end
if nargin < 1 || isempty(nSub)
    nSub = 10;
end

% Unprocessed subjects
try
    allSub = dir(in_dir);
    allSub = {allSub(3:end).name};
catch
    error('No data found in input dir! Could be file permission error.');
end
if isempty(allSub)  % Not sure why it couldn't be catched, but we do need this
    error('No data found in input dir! Could be file permission error.');
end
try
    oldSub = dir(fullfile(out_dir, 'Out'));
    oldSub = {oldSub(3:end).name};
catch
    oldSub = {};
end
newSub = setdiff(allSub, oldSub);

% Processing
nSub = min(nSub, numel(newSub));
Subjlist = newSub(1:nSub);
t1 = tic;
MyStartServerDT(Subjlist, 3, 'y', out_dir, in_dir);
disp(['Elapsed time for all ' num2str(nSub) ' subjects:'])
toc(t1);

end
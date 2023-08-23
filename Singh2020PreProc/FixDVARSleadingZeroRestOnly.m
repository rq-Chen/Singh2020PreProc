%% FixDVARSleadingZeroRestOnly.m - Add a leading zero to DVARS
%
% Just a quick fix. Add a line with a leading zero to the DVARS file.
% Then refilter the DVARS file and rerun Myc_fcprocess_HCP_dtseries.m
% and MyAllYeoParcellate.m
%
% Note: need to run this script after switching to its directory.

%% Setup

% Directory of this script
scriptdir = fileparts(mfilename('fullpath'));
oldpwd = pwd;
cd(scriptdir);

% Output directory
outdir = fullfile('..', 'GSR3');

% Filter
switches = SwitchSetup('HCP');
Filter = designfilt('bandstopfir', ...
    'FilterOrder', 40, ...
    'CutoffFrequency1', 0.06, ...
    'CutoffFrequency2', 0.14, ...
    'SampleRate', 1/switches.TR);


%% Add a leading zero to DVARS files and refilter

% Get the list of subjects
dirs = dir(fullfile(outdir, 'Out'));
dirs = {dirs(3:end-1).name};  % Exclude . and .. and 'outputs'

% Loop over subjects
for i = 1:numel(dirs)
    sessions = dir(fullfile(outdir, 'Out', dirs{i}, 'Results'));
    sessions = {sessions(3:end).name};  % Exclude . and ..
    for j = 1:numel(sessions)
        % DVARS file
        dvarsfile = fullfile(outdir, 'Out', dirs{i}, 'Results', sessions{j}, ...
            [dirs{i} '_' sessions{j} '_DVARS.txt']);
        dvarsfilteredfile = fullfile(outdir, 'Out', dirs{i}, 'Results', sessions{j}, ...
            [dirs{i} '_' sessions{j} '_DVARS_FILT.txt']);
        if ~isfile(dvarsfile)
            disp(['Cannot find ' dvarsfile])
            continue;
        end
        % Add a leading zero to the DVARS file
        dvars = [0; readmatrix(dvarsfile)];
        % Check if already processed
        if dvars(2) == 0
            disp(['Already processed ' dvarsfile])
            continue;
        end
        % Save the DVARS file
        writematrix(dvars, dvarsfile);
        % Refilter the DVARS file
        dvars = filtfilt(Filter, dvars);   
        % Save the filtered DVARS file
        writematrix(dvars, dvarsfilteredfile);
    end
end


%% Rerun Myc_fcprocess_HCP_dtseries.m and MyAllYeoParcellate.m

fieldtrip_dir = fullfile('utilities', 'fieldtrip');
addpath(fullfile(fieldtrip_dir))
addpath(fullfile(fieldtrip_dir, 'utilities'))
addpath(fullfile(fieldtrip_dir, 'external', 'freesurfer'))
addpath(fullfile(fieldtrip_dir, 'fileio', 'special')) % (some fieldtrip stuff used by Siegel..not mentioned, though)
addpath(fullfile(fieldtrip_dir, 'fileio'))

in_dir = '/net/10.27.136.121/hcpdb/packages/unzip/HCP_1200/';
dtseries_dir = fullfile(outdir, 'Out');

%Using HCP data or DMCC data?
Data.Type ='HCP';%HCP or DMCC
tseries = RestNameSetup(Data.Type);%sets up the tseries name used in either DMCC or HCP
switches = SwitchSetup(Data.Type);%sets up the parameters CHECK BEFORE RUNNING
switches.GSR = 3;

% Subject list
Subjlist = dirs;

for i = 1:length(Subjlist)
    fprintf('Processing subject %s', Subjlist{i});
    t0 = tic;
    Myc_fcprocess_HCP_dtseries(switches, Subjlist(i), in_dir, dtseries_dir, tseries)
    MyAllYeoParcellate(str2double(Subjlist(i)), switches.GSR, outdir, [], 'n');
    disp(['Elapsed time for subject ' Subjlist{i} ' :'])
    toc(t0)
end


%% Cleanup

cd(oldpwd);
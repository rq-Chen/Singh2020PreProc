%% Scrub Only
switches.AveRun = 1;
switches.GSR = 0;

switches.bpss = 1;
switches.bpss_lo = .08; % (recommend 0.08)
switches.bpss_hi = .009; % (recommend 0.009)

switches.skipfirst = 0; %skip this many frames at beginning of each run
switches.FDcut = 0.2; % threshold from petersen values to cut (recommend 0.04)
switches.DVARcut = 1.05; % percent above median to cutoff. zero means no cutoff (recommend 1.05)
switches.permutemask = 0; % randomly scrambles the tmask
switches.figs = 0; % generate figs to QC processing
switches.minframes = 400; % (recommend 400, unless you are afraid to throw away 1/4 subjects)
switches.TR = 0.720; % Acquisition rate in seconds
switches.FIX = 1;
switches.corrtype = 'pearson'; % 'pearson' or 'partial'
switches.runs = 1; % save run FC
switches.restart = 1; % save run FC
c_fcprocess_HCP_ptseries(switches)
c_fcprocess_HCP_dtseries(switches)
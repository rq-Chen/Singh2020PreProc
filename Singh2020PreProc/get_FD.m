function get_FD(sub, bold, movregs, out_dir)

% sub       : subject id (string)
% bold      : bold series name
% movregs   : movement regressors file
% out_dir   : location to write results

mov = dlmread(movregs); % or Movement_Regressors_dt?
Dtrans = mov(:,7:9); % these are the derivatives of translations
Drot = mov(:,10:12); % derivative of rotations, in degrees
% Drot units in HCP are degrees --> need conversion to radians,
% then multiply for a radius of 50mm (typical brain consensus for FDs)
Drot = (pi/180) * Drot * 50;
D = cat(2, Dtrans, Drot);

FD = sum(abs(D), 2);
disp('Writing to the Out Dir')
dlmwrite(fullfile(out_dir,[sub '_' bold '_FD.txt']), FD)
end

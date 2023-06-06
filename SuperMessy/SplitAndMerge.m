%This Function was built to correct for the the ordering of the parcells
%The ordering Used in Siegels paper and code is organized as Follows:
%Community1:
%   Left
%   Right
%Community2:
%   Left
%   Right
%The organization normally used in dlabels is the following
%Left:
%   Community1
%   Community2
%   etc
%Right:
%   Community1
%   Community2
%
%This function will split an array in two Left and Right the reorganize
%Right Now it just splits in half and interlives the parcels but that
%doesnt account for the different number of parcels on the left and right
%side may use raw_cii.bainordinate.parcellationlabel
function C = SplitAndMerge(Array)

s = size(Array,1);
half = ceil(s/2);
s1 = Array(1:half, :);
s2 = Array(half+1:s,:);
C = s1([1;1]*(1:size(s1,1)),:);
C(1:2:end,:) = s2;
end
%%function for getting th EpochIdx

function [EpochIdx,allbad] = getepochidx(i, datatoepoch)

for p = 1:size(datatoepoch,2)%for each response
    if i == p %stim == response 
        z = p+3;
    else
        z = p;
    end
    singlechanneldata = datatoepoch(:,z);
    EpochIdx = find((abs(diff(singlechanneldata))) > 100);
    if isempty(EpochIdx)
        continue;
    else
        %Do this 10 times to be thorough
        for count = 1:10
            for m = 2:(length(EpochIdx))
                D = EpochIdx(m)-EpochIdx(m-1);
                if abs(D) < 50
                    EpochIdx(m) = nan;
                else
                    continue;
                end 
            end
            Nanidx = isnan(EpochIdx);
            EpochIdx(Nanidx) = [];
            count = count+1;
        end
        if p == size(datatoepoch,2) %if it gets to the end and is still bad; just skip this stim channel
            allbad = 1;
            EpochIdx = zeros(10,10);
            return;
        elseif length(EpochIdx) < 21
            continue;%if fewer than 21 stims detected, it didnt actually pick up anything and should be excluded
        elseif length(EpochIdx) > 35 
            continue;
        else
            allbad = 0;
            return;
        end
    end
end
if ~exist('allbad','var')
    allbad = 2;
end

%This is calculating "PW ratios" and the like as metric of which Bode plot
%is best... 
function [hnorm,fpeak,freqs_log,singValues_log,pfRatio,peakFlag,bode_matrix,freq_matrix] = calc_BodePlotParams(A,B,C,fs,minSet,maxSet,nSteps)
        % Algorithm:
        %   find local peaks
        %   identify tallest local peak
        %   find minimum value before that value, DC* at w1
        %   find DC* value on right side of peak, at w2
        %   compute AmpRatio: ratio of peak magnitude to DC*
        %   compute freqRatio: ratio of log(w2)/log(w1)
        %   new PF ratio = AmpRatio/freqRatio

        % AS DEFAULTS: I set minSet to -1, maxSet to 3, and nSteps to 100.

        sys = ss(A,B,C,0,1/fs);%discrete (!!)  
        [hnorm,fpeak] = hinfnorm(sys); 
        [singValues,freqs] = sigma(sys); %,logspace(minSet,maxSet,nSteps)); %Just to see
        
        %make Bode matrix
%         bode_matrix = zeros(size(A,1),length(freqs));
%         freq_matrix = zeros(size(A,1),length(freqs));
        for i = 1:size(A,1) %every response
            [bode_mag,~,bode_freqs] = bode(sys(i));
            bode_matrix{i} = bode_mag(1,:);
            freq_matrix{i} = bode_freqs';
        end


        singValues_log = log(singValues);
        freqs_log = log(freqs);
        [peaks,locs] = findpeaks(singValues_log,freqs_log);
        peakFlag = 0;
        if ~isempty(peaks)
            [maxPeak,indMaxPeak] = max(peaks);
            maxPeak_freq = locs(indMaxPeak);
            indMaxPeak = find(freqs_log == maxPeak_freq);
            [minValLeft,~] = min(singValues_log(1:indMaxPeak));
            rightSingVec = singValues_log(indMaxPeak:end);
            [~,indMinValRight_preScale] = min(abs(rightSingVec-minValLeft));
            indMinValRight = indMinValRight_preScale+indMaxPeak-1;
            ampRatio = maxPeak - minValLeft;
            freqRatio = freqs_log(indMinValRight)-maxPeak_freq;
            if freqRatio == 0
                freqRatio = freqs_log(indMinValRight+1)-maxPeak_freq;
            end
            pfRatio = ampRatio./freqRatio;
        else
            peakFlag = 1;
            pfRatio = 0;
        end
end

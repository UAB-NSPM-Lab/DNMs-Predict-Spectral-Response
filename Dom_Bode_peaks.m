%Flattening Bode plots and selecting dominant resonant frequencies

clearvars
cd('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison')
load("PercentFreqMatchAllRand.mat")
pt_list = string(fieldnames(PercentFreqMatchAllRand));
for j = 1:length(pt_list)
    pt = pt_list(j);
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt));
    load("ArtifactualChannels.mat")
    load("base.mat")
    cd TFMsBode_figures2/ModelParams/
    dirNames = dir;
    for i = 1:length(dirNames)
        D = dirNames(i).name;
        if length(D) < 5
            continue;
        else
            stimname = erase(dirNames(i).name,'_modelParams_SEEG.mat');
            stimname = string(stimname); 
            var_name = append(stimname,'_modelParams_SEEG.mat');
            load(var_name);
            [Peaks] = peaks(stimname,A,B,C,fs,Arts,base); 
            AllPeaks_Bode.(stimname) = Peaks.(stimname); 
            clearvars Peaks
        end  
    end
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt,'/TFMsBode_figures2/Bode-CCSR Freq Comparison'));
    save("AllPeaks_Bode","AllPeaks_Bode",'-v7.3');
end 

function [Peaks] = peaks(stimname,A,B,C,fs,Arts,base) 
    %% For a specified stimulation channel pair, i
    tnames = base.bipolarresplabels;
    tidx = find((Arts.(stimname)==1));
%     tidx(1:142) = []; tidx = tidx-142; %P0120 only
    tnames(tidx) = [];
    t_end = length(B); %number of non-artifactual response channels 
    for t=1:t_end %non-artifactual response channels 
        tname = string(tnames(t));
        sys = ss(A,B,C(t,:),0,1/fs);
%         [mag,phase,wout] = bode(sys,logspace(-0.2018,3.8085,200)); %wout is omega (frequency) 
        [mag,phase,wout] = bode(sys,logspace(0,3.0103,200)); 
        magnitude = log10(squeeze(mag)); %gain manipulation; "flattening the Bode curve"
        [slope1, intercept1, slope2, intercept2, x_c] = fit_piecewise(magnitude);
        if x_c >= 200
            x_c = 200;
            %pause
%         elseif ~islogical(x_c)
%             pause
        elseif x_c<0
            Peaks.(stimname).(tname).adj_gain = nan; % 
            Peaks.(stimname).(tname).freq_Hz = nan; 
            continue
%             return
%             x_c = 0;
%             pause
        end
        x_c = round(x_c);% changed round to floor
        magnitude(1:x_c) = magnitude(1:x_c) - slope1*(1:x_c)'-intercept1;
        magnitude(x_c+1:200) = magnitude(x_c+1:200) - slope2*(x_c+1:200)'-intercept2; %is 200 aribrary? 
        semilogx(wout/2/pi,magnitude)
        hold on
        [peakValues, peakLocations] = findNPeaks(magnitude, 3);
        w_Hz = wout/2/pi; 
        pL_Hz = w_Hz(peakLocations);
        Peaks.(stimname).(tname).adj_gain = peakValues; % 
        Peaks.(stimname).(tname).freq_Hz = w_Hz(peakLocations); 
        semilogx(wout(peakLocations)/2/pi,peakValues,'r*')
        semilogx(wout/2/pi,log10(squeeze(mag)))
        grid on
        xlabel("Frequency [Hz]")
        ylabel("Magnitude [dB]")
        xline(5)
        xline(30)
        %pause
        hold off
    end

    function [slope1, intercept1, slope2, intercept2, x_c] = fit_piecewise(V)
        N = length(V);
        x = (1:N)'; 
        
        piecewiseModel = fittype(@(slope1, intercept1, x_c, slope2,  x) ...
            (x <= x_c) .* (slope1 * x + intercept1) + (x > x_c) .* (slope2 * x + x_c*(slope1-slope2)+intercept1), ...
            'independent', 'x', 'coefficients', {'slope1', 'intercept1', 'x_c', 'slope2'});
        
        % Initial parameter guesses
        x_c0 = round(N/2);
        slope1_0 = (V(x_c0) - V(1)) / (x_c0 - 1);
        intercept1_0 = V(1);
        slope2_0 = (V(end) - V(x_c0)) / (N - x_c0);
        
        % Fit the model
        fit_result = fit(x, V, piecewiseModel, 'StartPoint', [slope1_0, intercept1_0, x_c0, slope2_0]);
        
        % Extract fitted parameters
        slope1 = fit_result.slope1;
        intercept1 = fit_result.intercept1;
        slope2 = fit_result.slope2;
        x_c = fit_result.x_c;
        intercept2 = x_c*(slope1-slope2)+intercept1;
        if x_c < 0
            return;
%             pause
        end
    end
    
    function [peakValues, peakLocations] = findNPeaks(y, N)
        [pks, locs, ~, prom] = findpeaks(y);
        
        % Sort peaks by prominence
        [~, idx] = sort(prom, 'descend');
        
        % Select the top N peaks (or all available if less than N found)
        numPeaks = min(N, length(pks));
        peakValues = pks(idx(1:numPeaks));
        peakLocations = locs(idx(1:numPeaks));
        % PeakValues()
    end
end
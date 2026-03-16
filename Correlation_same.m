%Calculating the Bode-CCSR correlation for true (N1, N2, and Avg CCSR curves for SAME S-R pair) and comparing to Bode-random CCSR correlation

clearvars  
%get patient list
cd /data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output
keyword = 'P';
addpath(genpath('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output'));
fstruct = dir(append('*',keyword,'*'));
z = 0;
for i = 1:length(fstruct)
    if fstruct(i).isdir == 1
        cd(fstruct(i).name)
        keyword2 = 'TFMsBode_figures';
        fstruct2 = dir(append('*',keyword2,'*'));
        if ~isempty(fstruct2)
            z = z+1;
            Dir(z).names = cellstr(fstruct(i).name);
            cd ../
        else
            cd ../
            continue
        end
    else
        continue;
    end
end
%Load Things
cd /data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output;
load("fvec2.mat");
freqlim = find(fvec<55 & fvec>4); 

%Loop through all patients and calculate Bode-CCSR correlation and Bode-Rand CCSR correlation
for patients = 1:length(Dir)
    Pt = string(Dir(patients).names);
    %Load Things
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',Pt));
    load(append(Pt,'_short.mat'));
    load("CCSRs2.mat");
    load("ArtifactualChannels.mat");
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',Pt,'/TFMsBode_figures/ModelReconstructions'));
    load("CC_array.mat");
    cc_avg = nanmean(CC_array,1);
    remove_idx = find(isnan(cc_avg));
    CC_array(:,remove_idx) = [];
    load("StimNames.mat");
    fs = short.fs;
    CorrN1 = nan(size(CC_array,1),size(CC_array,2));
    CorrN2 = nan(size(CC_array,1),size(CC_array,2));
    CorrAvg = nan(size(CC_array,1),size(CC_array,2));
    CorrNrand = nan(size(CC_array,1),size(CC_array,2));
    %% Pick FM curves
    stims = string(fieldnames(Arts));
    nancount = 0;
    for j = 1:length(stims) %each stim for non-artifactual Channels
        S = stims(j);
        arts = Arts.(S);
        arts(remove_idx) = [];
        arts_idx = find(arts == 1); 
        cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',Pt,'/TFMsBode_figures/Bode')); 
        load(append((S),'_Bode_Matrices.mat'));
        for jj = 1:length(arts_idx)
            bode_1 = bode_matrix;
            bode_1(arts_idx(jj):end)=[];%first part
            freq_1 = freq_matrix;
            freq_1(arts_idx(jj):end) = [];
            bode_1(end+1) = {zeros(1,100)};%insert placeholder cell
            freq_1(end+1) = {zeros(1,100)};
            bode_2 = bode_matrix(arts_idx(jj):end);%second bit
            freq_2 = freq_matrix(arts_idx(jj):end);
            bode_matrix = [bode_1 bode_2];
            freq_matrix = [freq_1 freq_2];
        end%repeat until length of bode_matrix is same length as response channels
        N1_window_start = ceil(0.510*fs);
        N1_window_end = ceil(0.55*fs);
        N1_length = N1_window_end - N1_window_start;
        N2_window_start = ceil(0.55*fs);
        N2_window_end = ceil(fs);
        N2_length = N2_window_end - N2_window_start;
        %% Get S-R FM curves 
        resps = short.resplabels;
        resps(remove_idx) = [];
        art_Rs = string(resps((arts==1)));
        for k = 1:length(resps)
            R = string(resps(k));
            if ismember(R, art_Rs) %Artifactual channel means no Bode
                CorrN1(j,k) = nan;
                CorrN2(j,k) = nan;
                CorrAvg(j,k) = nan;
                CorrNrand(j,k) = nan;
                nancount = nancount+1;
                continue
            end
            ccsr = CCSRs.ccsr.(S).(R);
            r = find(ismember(string(short.resplabels),string(R)));
            ccep = short.average_waveforms.(S)(:,r);
            bode = bode_matrix{1,k};
            bode_freqs = freq_matrix{1,k};
            fs = short.fs;
            if all(bode) == 0%Artifactual channel means no Bode
                continue
            end
            %Situate Bode
            bode_freqs_Hz = bode_freqs./(2*pi);
            cut_idx = find(bode_freqs_Hz >= fvec(end) & bode_freqs_Hz<=fvec(1));
            cut_start = cut_idx(1)-1;
            cut_end = cut_idx(end) + 1;
            bode_short = bode(cut_start:cut_end);
            bode_freqs_short = bode_freqs_Hz(cut_start:cut_end);
            bode_interp = interp1(bode_freqs_short,bode_short,fvec);
            [~,N1_loc] = findpeaks(ccep(N1_window_start:N1_window_end),'MinPeakDistance',N1_length-1);
            [~,N2_loc] = findpeaks(ccep(N2_window_start:N2_window_end),'MinPeakDistance',N2_length-1);
            N1_peak = N1_window_start + N1_loc;
            N2_peak = N2_window_start+ N2_loc;
            if ~isempty(N1_peak)
                N1_vec = ccsr(:,N1_peak);
            else
                N1_vec = zeros(60,1);
            end
            if ~isempty(N2_peak)
                N2_vec = ccsr(:,N2_peak);
            else
                N2_vec = zeros(60,1);
            end
            N_vec_avg = mean(ccsr(:,N1_window_start:N2_window_end),2);
            rand_n = randi([N1_window_start N2_window_end]);
            N_vec_rand = ccsr(:,rand_n);
            CCSR_curves.(S).(R).N1 = N1_vec;
            CCSR_curves.(S).(R).N2 = N2_vec;
            CCSR_curves.(S).(R).Avg = N_vec_avg;
            CCSR_curves.(S).(R).Rand = N_vec_rand;
            %Correlated Curves?
            cor1 = corrcoef(N1_vec(freqlim),bode_interp(freqlim)); cor1 = cor1(1,2); 
            cor2 = corrcoef(N2_vec(freqlim),bode_interp(freqlim)); cor2 = cor2(1,2);
            coravg = corrcoef(N_vec_avg(freqlim),bode_interp(freqlim)); coravg = coravg(1,2); 
            corrand = corrcoef(N_vec_rand(freqlim),bode_interp(freqlim)); corrand = corrand(1,2);
            CorrN1(j,k) = cor1;
            CorrN2(j,k) = cor2;
            CorrAvg(j,k) = coravg;
            CorrNrand(j,k) = corrand;
        end
    end
    figure; 
    subplot(2,2,1); heatmap(CorrN1); clim([-1 1]); colorbar; title("N1-Bode"); grid off
    subplot(2,2,2); heatmap(CorrN2); clim([-1 1]); colorbar; title("N2-Bode");grid off
    subplot(2,2,3); heatmap(CorrAvg); clim([-1 1]); colorbar; title("Avg-Bode");grid off
    subplot(2,2,4); heatmap(CorrNrand); clim([-1 1]); colorbar; title("Rand-Bode");grid off
    sgtitle(append(Pt," Bode-CCSR Comparison - Same S and R"));
    colormap("jet")
    cd ../
    mkdir("Bode-CCSR Freq Comparison")
    cd("Bode-CCSR Freq Comparison")
    save("CC_array","CC_array")
    save("CorrNrand","CorrNrand")
    save("CorrAvg","CorrAvg")
    save("CorrN1","CorrN1")
    save("CorrN2","CorrN2")
    save("CCSR_curves", "CCSR_curves")
    cd ../
    disp(append("Done with ",Pt));
    clearvars -except Dir patients not_pts fvec freqlim
end

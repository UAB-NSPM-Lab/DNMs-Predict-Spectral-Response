%Calculating the Bode-CCSR correlation for true (N1, N2, and Avg CCSR curves for different S-R pairs) and comparing to Bode-random CCSR correlation 

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
    load("ArtifactualChannels.mat");
    load("base.mat");
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',Pt,'/TFMsBode_figures/ModelReconstructions'));
    load("CC_array.mat");
    cc_avg = nanmean(CC_array,1);
    remove_idx = find(isnan(cc_avg));
    CC_array(:,remove_idx) = [];
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',Pt,'/TFMsBode_figures/Bode-CCSR Freq Comparison'));
    load("CCSR_curves.mat");
    load("CorrAvg.mat");
    load("CorrN1.mat");
    load("CorrN2.mat");
    load("CorrNrand.mat");
    load("freqmatch_Bavg.mat");
    load("freqmatch_BN1.mat");
    load("freqmatch_BN2.mat");
    load("freqmatch_Brand.mat");
    load("PercentFreqMatch.mat");
    RandCorrN1 = nan(size(CC_array,1),size(CC_array,2),100);
    RandCorrN2 = nan(size(CC_array,1),size(CC_array,2),100);
    RandCorrAvg = nan(size(CC_array,1),size(CC_array,2),100);
    RandCorrNrand = nan(size(CC_array,1),size(CC_array,2),100);
    meanRandCorrN1 = nan(size(CC_array,1),size(CC_array,2));
    meanRandCorrN2 = nan(size(CC_array,1),size(CC_array,2));
    meanRandCorrAvg = nan(size(CC_array,1),size(CC_array,2));
    meanRandCorrNrand = nan(size(CC_array,1),size(CC_array,2));
    freqmatch_BN1 = nan(size(CC_array,1),size(CC_array,2));
    freqmatch_BN2 = nan(size(CC_array,1),size(CC_array,2));
    freqmatch_Bavg = nan(size(CC_array,1),size(CC_array,2));
    freqmatch_Brand = nan(size(CC_array,1),size(CC_array,2));
    meanFreqMatchRBN1 = nan(size(CC_array,1),size(CC_array,2));
    meanFreqMatchRBN2 = nan(size(CC_array,1),size(CC_array,2));
    meanFreqMatchRBavg = nan(size(CC_array,1),size(CC_array,2));
    meanFreqMatchRBrand = nan(size(CC_array,1),size(CC_array,2));
    %FM curves
    fs = base.fs;
    stims = string(fieldnames(CCSR_curves));
    Rand.Names = string(nan(1000,2));
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',Pt,'/TFMsBode_figures/Bode'))
    for i = 1:1000
        Rand.Names(i,1) = stims(randperm(length(stims),1));
        resps = string(base.bipolarresplabels(~Arts.(Rand.Names(i,1))));
        Rand.Names(i,2) = resps(randperm(length(resps),1));
        b_idx = find(ismember(resps,Rand.Names(i,2)));
        load(append((Rand.Names(i,1)),'_Bode_Matrices.mat'));
        Rand.Bode(i,1) = bode_matrix(b_idx);
        Rand.Freq(i,1) = freq_matrix(b_idx);
        clearvars bode_matrix freq_matrix
    end
    for j = 1:length(stims) %each stim for non-artifactual Channels
        S = string(stims(j));
        allRs = base.bipolarresplabels;
        allRs(remove_idx) = [];
        r = 0;
        for k = 1:length(allRs)
            if isnan(CorrAvg(j,k))
                continue;
            end
            r = r+1;
            R_list = base.bipolarresplabels(~Arts.(S));
            R = string(R_list(r));
            %ground truth
            ccsr_N1 = CCSR_curves.(S).(R).N1;
            ccsr_N2 = CCSR_curves.(S).(R).N2;
            ccsr_Avg = CCSR_curves.(S).(R).Avg;
            ccsr_Rand = CCSR_curves.(S).(R).Rand;
            for kk = 1:length(Rand.Bode)
                bode = Rand.Bode{kk,1};
                bode_freqs = Rand.Freq{kk,1};
                bode_freqs_Hz = bode_freqs./(2*pi);
                cut_idx = find(bode_freqs_Hz >= fvec(end) & bode_freqs_Hz<=fvec(1));
                cut_start = cut_idx(1)-1;
                cut_end = cut_idx(end) + 1;
                bode_short = bode(cut_start:cut_end);
                bode_freqs_short = bode_freqs_Hz(cut_start:cut_end);
                bode_interp = interp1(bode_freqs_short,bode_short,fvec);
                %corr
                cor1 = corrcoef(ccsr_N1(freqlim),bode_interp(freqlim)); cor1 = cor1(1,2); 
                cor2 = corrcoef(ccsr_N2(freqlim),bode_interp(freqlim)); cor2 = cor2(1,2);
                coravg = corrcoef(ccsr_Avg(freqlim),bode_interp(freqlim)); coravg = coravg(1,2); 
                corrand = corrcoef(ccsr_Rand(freqlim),bode_interp(freqlim)); corrand = corrand(1,2);
                RandCorrN1(j,k,kk) = cor1;
                RandCorrN2(j,k,kk) = cor2;
                RandCorrAvg(j,k,kk) = coravg;
                RandCorrNrand(j,k,kk) = corrand;
            end      
        end
    end
    meanRandCorrNrand = nanmean(RandCorrNrand,3);
    z_idx = find(meanRandCorrNrand==0); meanRandCorrNrand(z_idx) = nan;
    meanRandCorrN1 = nanmean(RandCorrN1,3);
    z_idx = find(meanRandCorrN1==0); meanRandCorrN1(z_idx) = nan;
    meanRandCorrN2 = nanmean(RandCorrN2,3);
    z_idx = find(meanRandCorrN2==0); meanRandCorrN2(z_idx) = nan;
    meanRandCorrAvg = nanmean(RandCorrAvg,3);
    z_idx = find(meanRandCorrAvg==0); meanRandCorrAvg(z_idx) = nan;
    meanFreqMatchRBN1 = nanmean(freqmatch_BN1,3);
    z_idx = find(meanFreqMatchRBN1==0); meanFreqMatchRBN1(z_idx) = nan;
    meanFreqMatchRBN2 = nanmean(freqmatch_BN2,3);
    z_idx = find(meanFreqMatchRBN2==0); meanFreqMatchRBN2(z_idx) = nan;
    meanFreqMatchRBavg = nanmean(freqmatch_Bavg,3);
    z_idx = find(meanFreqMatchRBavg==0); meanFreqMatchRBavg(z_idx) = nan;
    meanFreqMatchRBrand = nanmean(freqmatch_Brand,3);
    z_idx = find(meanFreqMatchRBrand==0); meanFreqMatchRBrand(z_idx) = nan;
    nancount = sum((isnan(meanFreqMatchRBN1)),"all");
    figure; 
    subplot(1,3,1); heatmap(meanRandCorrN1); clim([-1 1]); colorbar; title("N1-Bode");grid off
    subplot(1,3,2); heatmap(meanRandCorrN2); clim([-1 1]); colorbar; title("N2-Bode");grid off
    subplot(1,3,3); heatmap(meanRandCorrAvg); clim([-1 1]); colorbar; title("Avg-Bode");grid off
    sgtitle(append(Pt," Bode-CCSR Comparison - Rand S and R - Average"));
    colormap("jet")
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',Pt,'/TFMsBode_figures/Bode-CCSR Freq Comparison'));
    save("meanRandCorrNrand","meanRandCorrNrand");
    save("meanRandCorrAvg","meanRandCorrAvg");
    save("meanRandCorrN1","meanRandCorrN1");
    save("meanRandCorrN2","meanRandCorrN2");
    save("RandCorrNrand","RandCorrNrand",'-v7.3');
    save("RandCorrAvg","RandCorrAvg",'-v7.3');
    save("RandCorrN1","RandCorrN1",'-v7.3');
    save("RandCorrN2","RandCorrN2",'-v7.3');
    save("Rand","Rand",'-v7.3');
    cd ../
    disp(append("Done with ",Pt));
    clearvars -except Dir patients not_pts fvec freqlim
end

clearvars
cd('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison')
load("SplitCorr.mat")

pt_list = string(fieldnames(SplitCorr));
pt_list_ignore = 2;
pt_list(pt_list_ignore) = [];
for i = 1:length(pt_list)
    a = 0;
    pt = pt_list(i);
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt));
    load("base.mat");
    load("ArtifactualChannels.mat")
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt,'/TFMsBode_figures2/Bode-CCSR Freq Comparison'))
    load("ActualDomFreqs.mat");
    load("AllPeaks_Bode.mat");
    load("Rand.mat")
    s_list = string(fieldnames(Arts));
    r_list = base.bipolarresplabels;
    S = fieldnames(Arts);
    M_avg = zeros(length(S),length(Arts.(string(S(1)))));
    for m = 1:length(S)
        M_avg(m,:) = Arts.(string(S(m)));
    end
    F = find(mean(M_avg,1) == 1);
    r_list(F) = [];
    for h = 1:1000
        BS = Rand.Names(h,1);
        BR = Rand.Names(h,2);
        bode_freqs = round(AllPeaks_Bode.(BS).(BR).freq_Hz);
        for j = 1:length(s_list)
            S = s_list(j);
            r_cand_list = fieldnames(AllPeaks_Bode.(S));
            for k = 1:length(r_list)
                K = k*3-2;
                R = string(r_list(k));
                if ~ismember(R,r_cand_list)
                    RandFreqMatches.(pt).Avg{j,k} = "NaN";
                    avglayer(j,K) = nan;
                    avglayer(j,K+1) = nan;
                    avglayer(j,K+2) = nan;
                    RandFreqMatches.(pt).N1{j,k} = "NaN";
                    N1layer(j,K) = nan;
                    N1layer(j,K+1) = nan;
                    N1layer(j,K+2) = nan;
                    RandFreqMatches.(pt).N2{j,k} = "NaN";
                    N2layer(j,K) = nan;
                    N2layer(j,K+1) = nan;
                    N2layer(j,K+2) = nan;
                    RandFreqMatches.(pt).Baseline{j,k} = "NaN";
                    Baselinelayer(j,K) = nan;
                    Baselinelayer(j,K+1) = nan;
                    Baselinelayer(j,K+2) = nan;
                    continue;
                end
                if isnan(ActualDomFreqs.Avg{j,k})
                    RandFreqMatches.(pt).Avg{j,k} = "NaN";
                    avglayer(j,K) = nan;
                    avglayer(j,K+1) = nan;
                    avglayer(j,K+2) = nan;
                    RandFreqMatches.(pt).N1{j,k} = "NaN";
                    N1layer(j,K) = nan;
                    N1layer(j,K+1) = nan;
                    N1layer(j,K+2) = nan;
                    RandFreqMatches.(pt).N2{j,k} = "NaN";
                    N2layer(j,K) = nan;
                    N2layer(j,K+1) = nan;
                    N2layer(j,K+2) = nan;
                    RandFreqMatches.(pt).Baseline{j,k} = "NaN";
                    Baselinelayer(j,K) = nan;
                    Baselinelayer(j,K+1) = nan;
                    Baselinelayer(j,K+2) = nan;
                    continue;
                end
                a = a+1;
                avg_freqs = round(ActualDomFreqs.Avg{j,k});%in Hz
                n1_freqs = round(ActualDomFreqs.N1{j,k});
                n2_freqs = round(ActualDomFreqs.N2{j,k});
                baseline_freqs = round(ActualDomFreqs.Baseline{j,k});
    
                a_length.(pt)(a) = length(avg_freqs);
                n1_length.(pt)(a) = length(n1_freqs);
                n2_length.(pt)(a) = length(n2_freqs);
                n_length.(pt)(a) = length(baseline_freqs);
                
                I_avg = intersect(avg_freqs,bode_freqs);
                if ~isempty(I_avg)
                    RandFreqMatches.(pt).Avg{j,k} = I_avg;
                    avglayer(j,K) = I_avg(1);
                    if length(I_avg)>1
                        avglayer(j,K+1) = I_avg(2);
                    else
                        avglayer(j,K+1) = nan;
                    end
                    if length(I_avg)>2
                        avglayer(j,K+2)=I_avg(3);
                    else
                        avglayer(j,K+2) = nan;
                    end
                else
                    RandFreqMatches.(pt).Avg{j,k} = "NaN";
                    avglayer(j,K) = nan;
                    avglayer(j,K+1) = nan;
                    avglayer(j,K+2) = nan;
                    M_avg(j,k) = 0.5;
                end
                clearvars I_avg avg_freqs
                
                I_n1 = intersect(n1_freqs,bode_freqs);
                if ~isempty(I_n1)
                    RandFreqMatches.(pt).N1{j,k} = I_n1;
                    N1layer(j,K) = I_n1(1);
                    if length(I_n1)>1
                        N1layer(j,K+1) = I_n1(2);
                    else
                        N1layer(j,K+1) = nan;
                    end
                    if length(I_n1)>2
                        N1layer(j,K+2)=I_n1(3);
                    else
                        N1layer(j,K+2) = nan;
                    end
                else
                    RandFreqMatches.(pt).N1{j,k} = "NaN";
                    N1layer(j,K) = nan;
                    N1layer(j,K+1) = nan;
                    N1layer(j,K+2) = nan;
                    M_n1(j,k) = 0.5;
                end
                clearvars I_n1 n1_freqs
            
                I_n2 = intersect(n2_freqs,bode_freqs);
                if ~isempty(I_n2)
                    RandFreqMatches.(pt).N2{j,k} = I_n2;
                    N2layer(j,K) = I_n2(1);
                    if length(I_n2)>1
                        N2layer(j,K+1) = I_n2(2);
                    else
                        N2layer(j,K+1) = nan;
                    end
                    if length(I_n2)>2
                        N2layer(j,K+2)=I_n2(3);
                    else
                        N2layer(j,K+2) = nan;
                    end
                else
                    RandFreqMatches.(pt).N2{j,k} = "NaN";
                    N2layer(j,K) = nan;
                    N2layer(j,K+1) = nan;
                    N2layer(j,K+2) = nan;
                    M_n2(j,k) = 0.5;
                end
                clearvars I_n2 n2_freqs
    
                I_baseline = intersect(baseline_freqs,bode_freqs);
                if ~isempty(I_baseline)
                    RandFreqMatches.(pt).Baseline{j,k} = I_baseline;
                    Baselinelayer(j,K) = I_baseline(1);
                    if length(I_baseline)>1
                        Baselinelayer(j,K+1) = I_baseline(2);
                    else
                        Baselinelayer(j,K+1) = nan;
                    end
                    if length(I_baseline)>2
                        Baselinelayer(j,K+2)=I_baseline(3);
                    else
                        Baselinelayer(j,K+2) = nan;
                    end
                else
                    RandFreqMatches.(pt).Baseline{j,k} = "NaN";
                    Baselinelayer(j,K) = nan;
                    Baselinelayer(j,K+1) = nan;
                    Baselinelayer(j,K+2) = nan;
                    M_baseline(j,k) = 0.5;
                end
                clearvars I_baseline baseline_freqs
            end
        end
        Match_avg = find(M_avg==0);
        NonMatch_avg = find(M_avg==0.5);
        CXRrandom.(pt).avg(h) = length(Match_avg)/(length(Match_avg)+length(NonMatch_avg));
        Match_n1 = find(M_n1==0);
        NonMatch_n1 = find(M_n1==0.5);
        CXRrandom.(pt).n1(h) = length(Match_n1)/(length(Match_n1)+length(NonMatch_n1));
        Match_n2 = find(M_n2==0);
        NonMatch_n2 = find(M_n2==0.5);
        CXRrandom.(pt).n2(h) = length(Match_n2)/(length(Match_n2)+length(NonMatch_n2));
        Match_baseline = find(M_baseline==0);
        NonMatch_baseline = find(M_baseline==0.5);
        CXRrandom.(pt).baseline(h) = length(Match_baseline)/(length(Match_baseline)+length(NonMatch_baseline));
        clearvars bode_freqs Match_baseline Match_avg Match_n2 Match_rand Match_n1 NonMatch_baseline NonMatch_rand NonMatch_n2 NonMatch_n1 NonMatch_avg M_baseline M_rand M_n2 M_n1 M_avg
    end
    clearvars -except RandFreqMatches i pt_list CustomColormap CXRrandom pt
    save("CXRrandom.mat","CXRrandom",'-v7.3');
    save("RandFreqMatches.mat","RandFreqMatches",'-v7.3');
    disp(append("Done with ",(pt)));
    clearvars CXRrandom RandFreqMatches
end
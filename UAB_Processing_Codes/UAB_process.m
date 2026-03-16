function UAB_process(pt)

%paths
code_path = '/data/project/NSPMlab/hbriny99/UAB_Processing_Codes';
ccep_data_path = '/data/project/NSPMlab/hbriny99/UAB_CCEPdata';
out_path = '/data/project/NSPMlab/hbriny99/UAB_Output';

%read .edf
dirname = UAB_edfread(pt,ccep_data_path,out_path);

%preprocess
out = UAB_preprocess(dirname,ccep_data_path,out_path);

%remove artifacts
[out, base, Arts] = UAB_removeartifacts(out);

%filter
[filtered_data] = UAB_filter(out,Arts,base,pt,out_path);

%calculate average waveforms
[avg_waveforms] = UAB_avg_waveforms(pt,filtered_data,out_path);
UAB_plot_avg_waveforms(avg_waveforms)

%build models and Bode plots
UAB_Bode(pt, avg_waveforms, code_path,out_path);

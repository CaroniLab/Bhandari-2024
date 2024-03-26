clear
read_Intan_RHD2000_file

fileinfo = dir(fullfile(path, 'time.dat'));
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
fid = fopen(fullfile(path,'time.dat'), 'r');
t = fread(fid, num_samples, 'int32');
fclose(fid);
t = t / frequency_parameters.amplifier_sample_rate;

num_channels = length(board_adc_channels); % ADC input info from header file
fileinfo = dir(fullfile(path,'analogin.dat'));
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen(fullfile(path,'analogin.dat'), 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
v = v * 0.000050354;


ch =2;
fileinfo = dir(fullfile(path,'digitalin.dat'));
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen(fullfile(path,'digitalin.dat'), 'r');
digital_word = fread(fid, num_samples, 'uint16');
fclose(fid);
digital_input_ch = (bitand(digital_word, 2^ch) > 0);


ch_beam = 5;
fileinfo = dir(fullfile(path,'digitalin.dat'));
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen(fullfile(path,'digitalin.dat'), 'r');
digital_word = fread(fid, num_samples, 'uint16');
fclose(fid);
digital_input_ch_beam = (bitand(digital_word, 2^ch_beam) > 0);

writeNPY(v, fullfile(path,'analog.npy'))
writeNPY(digital_input_ch, fullfile(path,'camera.npy'))
writeNPY(t, fullfile(path,'time.npy'))
writeNPY(digital_input_ch_beam, fullfile(path,'beam.npy'))

%fileinfo = dir(fullfile(path,'amplifier.dat'));
%num_channels = length(amplifier_channels); % amplifier channel info from header file
%num_samples = fileinfo.bytes/(num_channels * 2);
%fid = fopen(fullfile(path,'amplifier.dat'), 'r');
%v_amp = fread(fid, [num_channels, num_samples], 'int16');
%fclose(fid);
%v = v * 0.195;


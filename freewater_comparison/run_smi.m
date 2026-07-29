addpath('/home/user/SMI'); addpath(pwd); addpath(fullfile(pwd,'stubs'));
warning('off','all'); pkg load statistics; pkg load image;
args = argv();
gen_and_fit_smi(str2double(args{1}), args{2});

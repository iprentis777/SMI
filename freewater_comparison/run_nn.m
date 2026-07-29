addpath('/home/user/SMI'); addpath(pwd); addpath(fullfile(pwd,'stubs'));
warning('off','all'); pkg load statistics; pkg load image;
args = argv(); sweep_nonneg(args{1});

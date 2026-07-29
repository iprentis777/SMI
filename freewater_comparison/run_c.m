addpath('/home/user/SMI'); addpath(pwd); addpath(fullfile(pwd,'stubs'));
warning('off','all'); pkg load statistics; pkg load image;
a = argv();
if strcmp(a{1},'cond'), gen_compartment(str2double(a{2}), a{3});
else, gen_brain(str2double(a{2}), a{3}, a{4}); end

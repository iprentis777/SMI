addpath('/home/user/SMI'); addpath(pwd); addpath(fullfile(pwd,'stubs'));
warning('off','all');
IO=binio();
bs=[0 0.5 1 2 3 4.5];
ks=[0.36 2.0 2.0 0.50 0.40; 0.60 2.0 2.0 0.50 0.00; 0.02 2.0 3.0 3.00 0.95];
R=zeros(size(ks,1),numel(bs),5);
for ik=1:size(ks,1)
  kv=[ks(ik,:) 1 1];
  for il=0:2:8
    R(ik,:,il/2+1)=SMI.RotInv_Kell_wFW_b_beta_TE_numerical(il,bs,ones(size(bs)),zeros(size(bs)),kv,3);
  end
end
IO.save('_kcheck',R); IO.save('_kcheck_b',bs); IO.save('_kcheck_k',ks);
printf('octave Kell written\n');

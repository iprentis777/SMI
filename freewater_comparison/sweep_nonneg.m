function sweep_nonneg(tag)
% sweep_tikhonov(tag)
%
% Re-fits the ALREADY GENERATED dwi_<tag> at a range of lambda_tikhonov, so
% every lambda sees byte-identical data and the only thing that varies is the
% regularizer. Motivated by REPORT_fODF_freewater.md section 4: SMI never
% resolves a 45 degree crossing, not even with zero noise, while CSD at b=3
% resolves it in 24/40 realisations. lambda_tikhonov was tuned for stability
% and has never been traded off against angular resolution.
%
% Writes smi_amp_lam<i>_<tag> plus the lambda values used.

more off
IO = binio();

LAMBDAS  = [0 1 3 10 30];  % lambda_nonneg; -1 means no constraint at all
LMAX_FIT = 6;
CS       = 1;
D_FW     = 3;

bvals     = IO.load('bvals')(:)';
bvecs     = IO.load('bvecs');
eval_dirs = IO.load('eval_dirs');
dwi       = IO.load(['dwi_' tag]);
GRID      = size(dwi); GRID = GRID(1:3);
NVOX      = prod(GRID);

% sigma is recoverable from the tag; keep it explicit rather than inferred
switch tag
    case 'snr30', sigma = 1/30;
    case 'snr15', sigma = 1/15;
    case 'clean', sigma = 1/100000;
    otherwise, error('unknown tag %s', tag);
end

Y = SMI.get_even_SH(eval_dirs, LMAX_FIT, CS);
L = repelem(0:2:LMAX_FIT, 2*(0:2:LMAX_FIT)+1)';

for il = 1:numel(LAMBDAS)
    lam = LAMBDAS(il);
    options = struct();
    options.b     = bvals;
    options.dirs  = bvecs;
    options.sigma = sigma*ones(GRID);
    options.mask  = true(GRID);
    options.compartments = {'IAS','EAS','FW'};
    options.NoiseBias    = 'Rician';
    options.Lmax         = [0 LMAX_FIT LMAX_FIT LMAX_FIT];
    options.CS_phase     = CS;
    options.D_FW         = D_FW;
    options.flag_fit_fODF = 1;
    options.fODF_regularization = struct('flag_nonneg',lam>0,'lambda_nonneg',max(lam,1e-6), ...
                                         'lambda_tikhonov',0.3);
    t0 = tic;
    out = SMI.fit(dwi, options);
    plm = reshape(out.plm, [NVOX size(out.plm,4)]);
    flm = [ones(NVOX,1) plm] .* repmat(sqrt((2*L+1)/(4*pi))', NVOX, 1);
    flm(~isfinite(flm)) = 0;
    IO.save(sprintf('smi_amp_nn%d_%s', il, tag), Y*flm');
    IO.save(sprintf('smi_pl_nn%d_%s', il, tag), ...
            reshape(out.pl, [NVOX size(out.pl,4)]));
    fprintf('lambda_nonneg %.2f  ->  %.1f s\n', lam, toc(t0));
end
IO.save(['nn_lambdas_' tag], LAMBDAS);
end

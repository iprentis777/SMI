function run_all_tests()
% RUN_ALL_TESTS Run the self-contained repository test scripts.
%
% From the repository root:
%   addpath('tests'); run_all_tests
%
% Each legacy test starts with `clear`, so it is executed in a small wrapper
% function to keep the runner's own state intact.

tests_dir = fileparts(mfilename('fullpath'));
test_files = {
    'test_fODF_outlier_cap.m'
    'test_SMI_outlier_cap.m'
    'test_SMI_response_helpers.m'
};

fprintf('Running %d self-contained SMI tests.\n', numel(test_files));
for it = 1:numel(test_files)
    fprintf('\n[%d/%d] %s\n', it, numel(test_files), test_files{it});
    run_test_file(fullfile(tests_dir, test_files{it}));
end
fprintf('\nAll self-contained SMI tests completed.\n');
end

function run_test_file(test_file)
run(test_file);
end

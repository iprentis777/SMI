# Tests

The scripts in this directory are self-contained and need no input data. Each
script resolves the repository root and shared `helpers/` directory from its
own location.

To run the full suite from the repository root:

```matlab
addpath('tests');
run_all_tests
```

To run individual checks from MATLAB:

```matlab
run('tests/test_fODF_outlier_cap.m')
run('tests/test_SMI_outlier_cap.m')
run('tests/test_SMI_response_helpers.m')
```

GNU Octave additionally needs the `statistics` and `image` packages. The test
scripts add the compatibility shims from `deconv_comparison/stubs/`
automatically.

The active simulation has one additional check, `deconv_comparison/test_csd_arms.m`,
kept with that package because it requires MRtrix3. It exercises the CSD arms
alone in ~2 s and is the regression test for the `-neg_lambda` bug.

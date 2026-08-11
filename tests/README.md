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

The active manuscript comparison has two additional checks in
`deconv_comparison/`: `check_manuscript_static.m` and `test_csd_arms.m`. They
are kept with that package because the latter requires MRtrix3.

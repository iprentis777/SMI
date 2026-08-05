# Tests

The scripts in this directory are self-contained and need no input data. Each
script resolves the repository root and shared `helpers/` directory from its
own location.

From MATLAB, run:

```matlab
run('tests/test_fODF_outlier_cap.m')
run('tests/test_SMI_outlier_cap.m')
run('tests/test_SMI_response_helpers.m')
```

GNU Octave additionally needs the `statistics` and `image` packages. The test
scripts add the compatibility shims from `deconv_comparison/stubs/`
automatically.

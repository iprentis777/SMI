# Examples

The examples fall into two groups: the original data-driven toolbox examples
and self-contained synthetic experiments added while developing the fODF work.

## Original toolbox example

[`example.m`](example.m) demonstrates four parameter-estimation workflows using
the three public SMI example datasets. Download the datasets from the
[SMI toolbox resource page](https://cai2r.net/resources/standard-model-of-diffusion-in-white-matter-the-smi-toolbox/)
and keep their directory names as `dataset_1`, `dataset_2`, and `dataset_3`.

Tell MATLAB where their parent directory is before running a section:

```matlab
setenv('SMI_EXAMPLE_DATA', 'C:\path\to\SMI_3datasets');
run('examples/example.m');
```

The example also uses the MATLAB NIfTI tools referenced in the root README.
Each `%% EXAMPLE` section performs its own path setup and can be run
independently.

## Self-contained examples

These scripts generate their own synthetic inputs and need no downloaded data:

| file | purpose |
|---|---|
| [`example_fODF_regularization.m`](example_fODF_regularization.m) | compare unregularized, Tikhonov, and non-negativity-constrained fODFs |
| [`example_fODF_regularization_sweep.m`](example_fODF_regularization_sweep.m) | sweep regularization settings and produce comparison figures |

The modulation and response-kernel viewers are retained as learning exercises,
not recommended pipeline steps. Their status and limitations are indexed in
[`Archive/README.md`](../Archive/README.md).

[`example_SMI_SSM.m`](example_SMI_SSM.m) is a specialized SSM example and is
not the primary first-run path for the toolbox.

# Historical patches

The numbered files in this directory preserve the measured changes as
mail-formatted patches. They duplicate changes already present in Git history
and are retained as project provenance, not as the recommended way to install
or use SMI.

Patch numbers record application order. When reconstructing the history, apply
them from lowest to highest with `git am --3way`. New contributors normally do
not need to read or apply these files; start with the root README and the
directory-specific documentation instead.

# A2MD Tests

I've realized that a big problem I have is that I have three
different libraries working in different examples. Therefore,
I decided to unify them into a common folder that will also be
version-controlled, so we can take account of the changes in the
examples.

## Contents

- gdb_000001: methane. It's here because it is easy, small and symmetric.
- gdb_000002: amonia. More or less the same.
- gdb_000003: water. Nothing to explain here
- gdb_000007: ethane. This is the simplest case of two-heavy atoms organic molecule.
- gdb_000037: random molecule with an eter group. I don't remember why it's here.
- gdb_000214: benzene. C6 symmetry and aromaticity. It's my favourite molecule.
- gdb_000397: propanotrio. I've never really understood why, but it has some strange
property that made all my previous fits crash, until I found out the importance of regularization.
It is also an example of symmetry.
- ani_000479: test to check issues of format in the new dataset

Within each folder you should be able to find .wfn files, .mol2 files, .csv files with values of
density and some picture. Try to keep this folder clean.

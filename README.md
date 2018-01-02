# DDF-ROM-POD
Data-Driven Filtered Reduced Order Modeling using POD

We implements data-driven filtered reduced order modeling (DDF-ROM) using
proper orthogonal decompositon (POD) in matlab. We calculated Snapshots, Mass,
Stiffness and other necessary (offline) matrices from DNS of a benchmark 2D channel flow past a cylinder problem. We upload sample offline matrices.

Here is the working procedure:

1. Run createROMDriver_35k.m. This will create ROM basis and ROM matrices, and save them in ROMtestSV35K_N16_166.mat. 


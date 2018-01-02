# DDF-ROM-POD
Data-Driven Filtered Reduced Order Modeling using POD

We implements data-driven filtered reduced order modeling (DDF-ROM) using
proper orthogonal decompositon (POD) in matlab. We calculated Snapshots, Mass,
Stiffness and other necessary (offline) matrices from DNS of a benchmark 2D channel flow past a cylinder problem. As the sizes of these offline matrices are huge, I uploaded sample offline matrices on my websites (these are all in snapshotData35Kdt002SV_Re100.mat file): http://www.math.vt.edu/people/jaman/Research.html .

Here is the working procedure:

1. Run createROMDriver_35k.m. This will create ROM basis and ROM matrices, and save them in ROMtestSV35K_N16_166.mat. 
2. Run createGsnap.m which will create the stress tensor \tau and save it as Gsnap_SV35K_r8_d16_N16_166.mat.
3. Run createABtilde_noconstraints.m to create the matrices \tilde{A} and \tilde{B} using truncated SVD. It saves the matrices as ABtilde_N16_r8_d16_166.mat.
4. Now we are ready to compute G-ROM and DDF-ROM solutions together using BDF2 (second-oder) time-stepping methods. Run ROMDriverPlot.m. For the first time step, we use backward-Euler (first order accurate) method. We compute lift, drag and energy for both RoMs and plot them together with results from DNS (i.e DNS projection onto ROM which is located in DNSProjectionMatrix_r8.mat file).


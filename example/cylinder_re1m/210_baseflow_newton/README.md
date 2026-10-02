# cylinder_re1m / 210_baseflow_newton

Probe. Newton fixed point (`userParam01 = 2`) started from the SFD field
`base_converged.f00001`. `userParam07 = 5`. Same mesh, dt, and k-tau
model as `../110_baseflow_sfd`. The first forward residual was 0.54. Five
Arnoldi columns left it there, rate 0.9999. The linear residual then jumped
to about 100 and stayed there. The next Newton residual was 0.542.
`nwt1cyl0.f00001` matches the SFD velocity exactly.
`newton_vs_sfd.png` is that comparison. This is not a converged base flow.

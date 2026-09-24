# MuonSystemArcVertexFinder

Displaced-vertex finding for long-lived-particle decays **behind the solenoid**, using
muon-system tracks only (e.g. `StandaloneMuonTracks`).

Author: Mahmoud Althakeel

## Method

The muon tracks' own parameters are defined at the IP; extrapolating them back through the
calorimeter and the field reversal at the solenoid gives metre-scale vertex errors. Instead,
each track is refitted **locally** from its muon-chamber hits and extrapolated back only to
the decay:

1. **Track model** - circle in xy + straight line in z, generalised least-squares fits in a
   frame attached to the innermost hit. The hit covariance includes the correlated multiple
   scattering in the iron yoke between stations.
2. **Error tube** - fit covariance propagated to the extrapolation point, plus multiple
   scattering in the calorimeter between that point and the first hit, plus a field-mismatch
   term (the circle is fitted in the chambers' field and extrapolated through the
   calorimeter region). The field term (0.5 T) is calibrated on truth so that the per-track
   pulls are ~1 at every extrapolation distance.
3. **Pair vertex** - Gauss-Newton chi2 fit, two residuals per track (xy distance
   perpendicular to the track, z at the same s): position, 3x3 covariance, chi2 (ndf = 1).
4. **Consistency cuts** - vertex upstream of both first hits, between the IP and the
   chambers, behind the solenoid, pair DCA, opposite charge.

## Running

    k4run runMuonSystemArcVertexFinder.py --inputFile tracks.root --outputFile vtx.root
    # --loose        : no DCA / chi2 / charge / region cut, all pairs kept (cut scans)
    # --noMaterial   : switch off the multiple-scattering part of the error tube
    # --fieldUncertainty <T> : field-mismatch term, 0 = off (default 0.5)
    # --noCurvature  : straight line instead of a circle in xy

Main properties: `CalorimeterX0` (90), `YokeX0` / `YokeThickness` (set per geometry: iron
between the stations and their radial span), `MomentumResolution` (0.20),
`FieldUncertainty` (0.5 T), `UseMaterialEffects`, `MaxPairDCA` (100 mm),
`RequireOppositeCharge`, `RequireBetweenIPAndHits`, `RequireBehindSolenoid`,
`MinUpstreamDistance`.

Vertex parameters: `[2, index1, index2, pairDCA, s1, s2, upstream, behindSolenoid,
chargeProduct, betweenIPAndHits]` (indices into the input collection; s = distance of the
vertex from each track's first hit, negative = upstream).

## Performance (HNL m_N = 20 GeV, |V_muN|^2 = 1e-10, decays at r_xy > 2.3 m)

Vertex cuts DCA < 300 mm, reconstructed r_xy > 2.3 m; Gaussian-core resolution:

| geometry | vertex found (given 2 opposite-charge tracks) | sigma x / y / z |
|---|---|---|
| 3L-30cm | 83% | 10.2 / 10.0 / 5.9 mm |
| 4L-50cm | 84% | 7.9 / 6.5 / 5.8 mm |
| 6L-30cm | 90% | 6.3 / 5.1 / 3.1 mm |

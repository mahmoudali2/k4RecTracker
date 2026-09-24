/*
 * Vertexing of long-lived-particle decays BEHIND the solenoid with muon-system
 * tracks only.
 *
 * Outside the coil the tracks cross only the calorimeter and the muon-system
 * yoke, and the field is weak. Rather than extrapolating a helix to a perigee at
 * the origin (metres away from the decay), each track is described LOCALLY from
 * its own muon-chamber hits and extrapolated back only as far as the decay:
 *
 *   track  : circle in xy + straight line z(s), generalised least-squares fits
 *            in a frame attached to the innermost hit; the hit covariance
 *            includes the correlated multiple scattering in the yoke between hits
 *   errors : fit covariance propagated to the extrapolation point, plus the
 *            multiple scattering in the material between that point and the
 *            first hit (dual-readout calorimeter)
 *   vertex : Gauss-Newton fit with two residuals per track (xy distance
 *            perpendicular to the track, z at the same s) -> position,
 *            3x3 covariance, chi2 with ndf = 2N - 3
 *
 * Units: mm, GeV, Tesla.
 * Author: Mahmoud Althakeel
 */
#pragma once

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVector3.h>

#include <optional>
#include <vector>

namespace k4vertexing::muonarc {

struct Config {
  // momentum used for the multiple-scattering estimate
  double bFieldForMomentum = 2.0;   // T, converts track omega [1/mm] into pT
  double momentumResolution = 0.20; // sigma(p)/p of the standalone muon system
  double minMomentum = 2.0;         // GeV, floor
  double hitResolution = 1.0;       // mm per coordinate
  // field mismatch between the chambers (where the circle is fitted) and the
  // calorimeter region it is extrapolated through; adds 0.5*0.3*dB/pT*s^2 to
  // the bending-plane error. 0.5 T calibrated on truth (flat per-track pulls)
  double fieldUncertainty = 0.5;    // T
  bool useMaterial = true;          // multiple scattering in yoke + calorimeter
  // material: dual-readout calorimeter between solenoid and muon system
  double caloX0 = 90.0;
  double caloInnerR = 2400.0, caloOuterR = 4500.0;
  double caloInnerZ = 2400.0, caloOuterZ = 4300.0;
  // material: muon-system yoke (between and beyond the stations)
  double yokeX0 = 40.0;
  double yokeThickness = 2000.0;
  bool fitCurvature = true;
  unsigned minHits = 3;
  // region and topology
  double solenoidR = 2300.0, solenoidZ = 2180.0;
  double upstreamTolerance = 100.0; // vertex may lie this far past a first hit
  int nSteps = 40;                  // steps of the material integrals
};

/// Radiation lengths per mm of path at a point (simple IDEA-like shells).
double x0PerMM(const TVector3& point, const Config& cfg);

class ArcTrack {
public:
  /// hits in mm; omega [1/mm] only sets the momentum for the scattering model
  bool build(std::vector<TVector3> hits, double omega, const Config& cfg);

  TVector3 point(double s) const;
  /// local coordinates (s, u) of a point in the frame at the first hit
  void local(const TVector3& v, double& s, double& u) const;
  /// u(s) of the fitted circle and its slope du/ds
  void uOfS(double s, double& u, double& slope) const;
  double zOfS(double s) const { return m_h0.Z() + m_pz[0] + m_pz[1] * s; }
  /// sigma of the trajectory perpendicular to it in xy, and in z, at local s
  void sigmas(double s, const Config& cfg, double& sigmaU, double& sigmaZ) const;

  const TVector3& firstHit() const { return m_h0; }
  double dirX() const { return m_t0[0]; }
  double dirY() const { return m_t0[1]; }
  double slopeZ() const { return m_pz[1]; }
  double momentum() const { return m_p; }
  double fitChi2() const { return m_chi2; }
  int nHits() const { return m_nHits; }

private:
  TVector3 m_h0;
  double m_t0[2] = {1, 0}, m_n0[2] = {0, 1};
  std::vector<double> m_pu; // a, b, (kappa)
  TMatrixDSym m_cu;
  double m_pz[2] = {0, 0};
  TMatrixDSym m_cz{2};
  double m_cosLambda = 1, m_p = 0, m_pt = 0, m_msFactor = 0, m_chi2 = 0;
  int m_nHits = 0;
};

struct VertexResult {
  TVector3 position;
  TMatrixDSym covariance{3};
  double chi2 = 0;
  int ndf = 0;
  double pairDCA = 0;            // distance between the two trajectories at the vertex
  std::vector<double> sFromFirst; // per track: local s of the vertex (negative = upstream)
  bool upstream = false;          // vertex before every track's first hit
  bool behindSolenoid = false;
  bool betweenIPAndHits = false;  // closer to the origin than every first hit, same side
};

std::optional<VertexResult> fitVertex(const std::vector<const ArcTrack*>& tracks, const Config& cfg,
                                      int maxIterations = 12);

} // namespace k4vertexing::muonarc

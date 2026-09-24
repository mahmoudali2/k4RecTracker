/*
 * Gaudi component: two-track vertexing of long-lived-particle decays BEHIND the
 * solenoid, using muon-system tracks only (e.g. StandaloneMuonTracks).
 *
 * Each track is refitted locally from its own muon-chamber hits (circle in xy,
 * line in z, with multiple scattering in the yoke) and extrapolated back only to
 * the decay, with the calorimeter scattering added to its error. Track pairs are
 * vertexed with a chi2 fit. See include/MuonArcVertexing.h.
 *
 * Author: Mahmoud Althakeel
 *
 * Default cuts = working point from the HNL (mN=20 GeV) study: highest efficiency
 * with < 0.1% fake vertices per event on prompt decays:
 *   MaxPairDCA = 100 mm, RequireOppositeCharge = true, RequireBetweenIPAndHits = true
 *   -> 67% efficiency for decays with r_xy > 2.3 m and >= 2 muon tracks,
 *      ~15-24 mm resolution (68%) per coordinate, 0.09% fakes per event.
 *
 * Input : edm4hep::TrackCollection  (tracks with muon-system hits attached)
 * Output: edm4hep::VertexCollection
 *   parameters = [ntracks(=2), index1, index2, pairDCA, s1, s2, upstream, behindSolenoid,
 *                 chargeProduct, betweenIPAndHits]
 *   (indices into the input collection; s = distance of the vertex from each
 *    track's first hit along the track, negative = upstream)
 */
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/VertexCollection.h>

#include "MuonArcVertexing.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace k4vertexing::muonarc;

struct MuonSystemArcVertexFinder final
    : k4FWCore::Transformer<edm4hep::VertexCollection(const edm4hep::TrackCollection&)> {

  MuonSystemArcVertexFinder(const std::string& name, ISvcLocator* serviceLocator)
      : Transformer(name, serviceLocator, {KeyValues("InputTracks", {"StandaloneMuonTracks"})},
                    {KeyValues("OutputVertices", {"MuonSystemVertices"})}) {}

  StatusCode initialize() override {
    const auto status = Transformer::initialize();
    if (!status.isSuccess())
      return status;
    m_cfg.bFieldForMomentum = m_bFieldForMomentum;
    m_cfg.momentumResolution = m_momentumResolution;
    m_cfg.minMomentum = m_minMomentum;
    m_cfg.hitResolution = m_hitResolution;
    m_cfg.fieldUncertainty = m_fieldUncertainty;
    m_cfg.useMaterial = m_useMaterial;
    m_cfg.caloX0 = m_caloX0;
    m_cfg.caloInnerR = m_caloInnerR;
    m_cfg.caloOuterR = m_caloOuterR;
    m_cfg.caloInnerZ = m_caloInnerZ;
    m_cfg.caloOuterZ = m_caloOuterZ;
    m_cfg.yokeX0 = m_yokeX0;
    m_cfg.yokeThickness = m_yokeThickness;
    m_cfg.fitCurvature = m_fitCurvature;
    m_cfg.minHits = m_minHits;
    m_cfg.solenoidR = m_solenoidR;
    m_cfg.solenoidZ = m_solenoidZ;
    m_cfg.upstreamTolerance = m_upstreamTolerance;
    if (m_caloOuterR <= m_caloInnerR || m_yokeThickness <= 0) {
      error() << "Inconsistent material description" << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "Muon-system arc vertexer: material " << (m_useMaterial ? "ON" : "OFF")
           << ", field uncertainty " << m_fieldUncertainty.value() << " T, calo " << m_caloX0.value() << " X0, yoke " << m_yokeX0.value()
           << " X0, sigma(p)/p = " << m_momentumResolution.value() << ", MaxPairDCA = " << m_maxPairDCA.value()
           << " mm, MaxChi2 = " << m_maxChi2.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  edm4hep::VertexCollection operator()(const edm4hep::TrackCollection& tracks) const override {
    edm4hep::VertexCollection vertices;

    // ---- local arc fit of every track from its muon hits -------------------
    std::vector<ArcTrack> arcs;
    std::vector<int> sourceIndex;
    std::vector<int> charge;
    int index = 0;
    for (const auto& track : tracks) {
      std::vector<TVector3> hits;
      for (const auto& hit : track.getTrackerHits())
        hits.emplace_back(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z);
      double omega = 0;
      for (const auto& state : track.getTrackStates())
        if (state.location == m_momentumStateLocation) {
          omega = state.omega;
          break;
        }
      ArcTrack arc;
      if (arc.build(hits, omega, m_cfg)) {
        arcs.push_back(std::move(arc));
        sourceIndex.push_back(index);
        // only the product of the two charges is used, so the sign convention
        // of the input state does not matter
        charge.push_back(omega > 0 ? 1 : (omega < 0 ? -1 : 0));
      }
      ++index;
    }
    debug() << "Arc-fitted " << arcs.size() << " of " << tracks.size() << " tracks" << endmsg;

    // ---- vertex every pair, apply the selection ----------------------------
    struct Candidate {
      VertexResult fit;
      std::size_t a, b;
    };
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < arcs.size(); ++i)
      for (std::size_t j = i + 1; j < arcs.size(); ++j) {
        auto fit = fitVertex({&arcs[i], &arcs[j]}, m_cfg);
        if (!fit)
          continue;
        if (m_requireUpstream && !fit->upstream)
          continue;
        if (m_requireBehindSolenoid && !fit->behindSolenoid)
          continue;
        if (m_requireBetween && !fit->betweenIPAndHits)
          continue;
        if (fit->pairDCA > m_maxPairDCA || fit->chi2 > m_maxChi2)
          continue;
        // vertex must lie at least MinUpstreamDistance before BOTH first hits:
        // rejects one muon split into two nearby tracks, which "meet" at the chambers
        if (*std::max_element(fit->sFromFirst.begin(), fit->sFromFirst.end()) > -m_minUpstreamDistance)
          continue;
        if (m_requireOppositeCharge && charge[i] * charge[j] >= 0)
          continue;
        candidates.push_back({std::move(*fit), i, j});
      }

    // ---- best chi2 first, each track used once ------------------------------
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& x, const Candidate& y) { return x.fit.chi2 < y.fit.chi2; });
    std::vector<bool> used(arcs.size(), false);
    for (const auto& cand : candidates) {
      if (!m_keepOverlapping && (used[cand.a] || used[cand.b]))
        continue;
      used[cand.a] = used[cand.b] = true;
      const auto& f = cand.fit;
      auto vertex = vertices.create();
      vertex.setPosition({static_cast<float>(f.position.X()), static_cast<float>(f.position.Y()),
                          static_cast<float>(f.position.Z())});
      const TMatrixDSym& c = f.covariance;
      vertex.setCovMatrix({static_cast<float>(c(0, 0)), static_cast<float>(c(0, 1)), static_cast<float>(c(1, 1)),
                           static_cast<float>(c(0, 2)), static_cast<float>(c(1, 2)), static_cast<float>(c(2, 2))});
      vertex.setChi2(static_cast<float>(f.chi2));
      vertex.setNdf(f.ndf);
      vertex.setPrimary(false);
      vertex.setAlgorithmType(m_algorithmType);
      for (float value : {2.f, static_cast<float>(sourceIndex[cand.a]), static_cast<float>(sourceIndex[cand.b]),
                          static_cast<float>(f.pairDCA), static_cast<float>(f.sFromFirst[0]),
                          static_cast<float>(f.sFromFirst[1]), f.upstream ? 1.f : 0.f,
                          f.behindSolenoid ? 1.f : 0.f, static_cast<float>(charge[cand.a] * charge[cand.b]),
                          f.betweenIPAndHits ? 1.f : 0.f})
        vertex.addToParameters(value);
      if (!m_allowMultipleVertices)
        break;
    }
    debug() << "Reconstructed " << vertices.size() << " muon-system vertices" << endmsg;
    return vertices;
  }

private:
  Config m_cfg;

  // --- momentum for the scattering model ---
  Gaudi::Property<int> m_momentumStateLocation{this, "MomentumStateLocation", edm4hep::TrackState::AtLastHit,
                                               "Track state whose omega gives the momentum"};
  Gaudi::Property<double> m_bFieldForMomentum{this, "BFieldForMomentum", 2.0,
                                              "Field [T] converting that state's omega into pT"};
  Gaudi::Property<double> m_momentumResolution{this, "MomentumResolution", 0.20,
                                               "sigma(p)/p of the standalone muon system"};
  Gaudi::Property<double> m_minMomentum{this, "MinMomentum", 2.0, "Momentum floor for scattering [GeV]"};
  Gaudi::Property<double> m_hitResolution{this, "HitResolution", 1.0, "Muon hit resolution [mm]"};
  Gaudi::Property<double> m_fieldUncertainty{
      this, "FieldUncertainty", 0.5,
      "Field mismatch [T] between chambers and calorimeter region (bending-plane error 0.5*0.3*dB/pT*s^2)"};
  // --- material ---
  Gaudi::Property<bool> m_useMaterial{this, "UseMaterialEffects", true,
                                      "Multiple scattering in yoke and calorimeter in the track errors"};
  Gaudi::Property<double> m_caloX0{this, "CalorimeterX0", 90.0,
                                   "Radiation lengths of the calorimeter between solenoid and muon system"};
  Gaudi::Property<double> m_caloInnerR{this, "CalorimeterInnerR", 2400.0, "[mm]"};
  Gaudi::Property<double> m_caloOuterR{this, "CalorimeterOuterR", 4500.0, "[mm]"};
  Gaudi::Property<double> m_caloInnerZ{this, "CalorimeterInnerZ", 2400.0, "[mm]"};
  Gaudi::Property<double> m_caloOuterZ{this, "CalorimeterOuterZ", 4300.0, "[mm]"};
  Gaudi::Property<double> m_yokeX0{this, "YokeX0", 40.0, "Radiation lengths of the muon-system yoke"};
  Gaudi::Property<double> m_yokeThickness{this, "YokeThickness", 2000.0,
                                          "Depth over which the yoke X0 is spread [mm]"};
  // --- track model ---
  Gaudi::Property<bool> m_fitCurvature{this, "FitCurvature", true, "Circle in xy (true) or straight line (false)"};
  Gaudi::Property<unsigned> m_minHits{this, "MinHits", 3, "Minimum muon hits per track"};
  // --- selection ---
  Gaudi::Property<bool> m_requireUpstream{this, "RequireUpstream", true,
                                          "Vertex must lie before every track's first hit"};
  Gaudi::Property<double> m_upstreamTolerance{this, "UpstreamTolerance", 100.0, "[mm]"};
  Gaudi::Property<bool> m_requireBehindSolenoid{this, "RequireBehindSolenoid", true,
                                                "Keep only vertices outside the solenoid volume"};
  Gaudi::Property<bool> m_requireBetween{this, "RequireBetweenIPAndHits", true,
                                         "Vertex closer to the origin than every first hit, on the same side"};
  Gaudi::Property<double> m_solenoidR{this, "SolenoidR", 2300.0, "[mm]"};
  Gaudi::Property<double> m_solenoidZ{this, "SolenoidHalfZ", 2180.0, "[mm]"};
  Gaudi::Property<double> m_maxPairDCA{this, "MaxPairDCA", 100.0,
                                       "Max distance between the two trajectories at the vertex [mm]"};
  Gaudi::Property<double> m_maxChi2{this, "MaxChi2", 1.0e9, "Max vertex chi2 (ndf = 1)"};
  Gaudi::Property<double> m_minUpstreamDistance{this, "MinUpstreamDistance", 0.0,
                                                "Vertex at least this far before both first hits [mm]"};
  Gaudi::Property<bool> m_requireOppositeCharge{this, "RequireOppositeCharge", true,
                                                "Keep only opposite-charge pairs (LLP -> mu+ mu- X)"};
  Gaudi::Property<bool> m_keepOverlapping{this, "KeepOverlappingPairs", false,
                                          "Write pairs sharing a track too (for cut studies)"};
  Gaudi::Property<bool> m_allowMultipleVertices{this, "AllowMultipleVertices", true,
                                                "Write every non-overlapping pair, not only the best"};
  Gaudi::Property<int> m_algorithmType{this, "AlgorithmType", 3, "edm4hep::Vertex algorithmType tag"};
};

DECLARE_COMPONENT(MuonSystemArcVertexFinder)

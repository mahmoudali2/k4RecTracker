# Author : Mahmoud Althakeel (mahmoud.althakeel@cern.ch)
#
# Steering file for running the standalone muon-system tracking algorithm
# (pattern recognition + optional GenFit fit over IDEA Muon-System hits).
# This script configures and runs the StandaloneMuonTracking Gaudi component.
#

import os
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import StandaloneMuonTracking
from Configurables import GeoSvc

# Geometry service - provides access to detector geometry
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"] + "/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

# Standalone muon tracking algorithm
standaloneMuonTracking = StandaloneMuonTracking()
standaloneMuonTracking.DetectorName = "Muon-System"           # Name of the detector to process
standaloneMuonTracking.ParticleType = "muon"                  # Particle type for material effects
standaloneMuonTracking.EncodingStringParameterName = "MuonSystemReadoutID"
standaloneMuonTracking.OutputLevel = INFO

# GenFit configuration
standaloneMuonTracking.UseGenFit = False                      # Enable GenFit track fitting
standaloneMuonTracking.MaxFitIterations = 100                 # Maximum iterations for the Kalman fit
# Inner propagation (disabled: solenoid boundary extrapolation not yet calibrated)
standaloneMuonTracking.DoInnerPropagation = False
standaloneMuonTracking.InnerPropTargetRadius = 100.0          # target radius in cm if enabled

# ── N-hit combinatorial strategy ─────────────────────────────────────────────
standaloneMuonTracking.MinTrackHits = 4                       # minimum hits to form a track
standaloneMuonTracking.MaxCombinatorialHits = 8               # cap on k in C(N,k) enumeration
standaloneMuonTracking.NHitMaxChi2NDF = 10.0                  # max chi2/NDF of N-hit circle fit
standaloneMuonTracking.HitIsolationCut = 0                    # min cross-region distance (cm); 0 to disable
standaloneMuonTracking.MinLayerSpan = 4                       # min distinct detector layers per combo
standaloneMuonTracking.MaxComboPT = 200.0                     # GeV: maximum combo pT (excludes combos above this threshold during search)
standaloneMuonTracking.OutlierSigma = 5.0                     # outlier removal threshold (σ units)
standaloneMuonTracking.SigmaHitDefault = 0.04                 # default hit resolution (cm = 0.4 mm)
standaloneMuonTracking.MaxOutlierIterations = 10              # max rounds of outlier removal
# Road cuts: guard against cross-track hit mixing
standaloneMuonTracking.MaxConsecDeltaPhi = 0.15               # max |Δφ| between consecutive layers (rad); -1 to disable
standaloneMuonTracking.MaxComboPhiSpread = 0.5                # max total φ spread across combo (rad); -1 to disable
# Distance cuts: guard against ghost combos spanning different detector regions
standaloneMuonTracking.MaxConsecutiveHitDistance = -1         # max 3D dist between adjacent (inner->outer) hits (cm); -1 to disable
standaloneMuonTracking.MaxPairHitDistance = 350.0             # max 3D dist between any two hits in a combo (cm); 300 accommodates barrel->endcap
# Hit quality pre-filter and proximity-based combo selection
standaloneMuonTracking.MaxHitEdepKeV = 0                      # reject hits with edep > 10 keV (delta-rays / hadronic secondaries)
standaloneMuonTracking.ProximityScoreWeight = 1.0             # weight for edep-homogeneity in combo ranking (0 to disable)
standaloneMuonTracking.EdepNormKeV = 3.0                      # normalisation: 2 keV avg deviation ~ +1 chi2 unit (for k>=4 combined score)
standaloneMuonTracking.MaxHitsPerCompositeID = 3              # 0 = use all hits; >0 caps per layer to top-N by proximity score
# conditions for second track in the events
standaloneMuonTracking.NeighbourTrackMaxDist = 10.0           # cm (= 50 mm)
standaloneMuonTracking.NeighbourTrackMinHits = 5              # must have >=4 hits
standaloneMuonTracking.NeighbourTrackMaxChi2NDF = 3.0         # must have chi2/NDF <= 3
# Guided crowded-layer hit selection
standaloneMuonTracking.doCrowdedLayerHitSelection = False     # enable guided crowded-layer hit selection

# ── Collections ──────────────────────────────────────────────────────────────
standaloneMuonTracking.InputHitCollection = ["MSTrackerHits"]
standaloneMuonTracking.OutputTrackCollection = ["StandaloneMuonTracks"]
standaloneMuonTracking.OutputHitCollection = ["StandaloneMuonTrackHits"]
standaloneMuonTracking.InputRecoSimLinkCollection = ["MSTrackerHitRelations"]
standaloneMuonTracking.OutputRecoSimLinkCollection = ["StandaloneMuonTrackHitSimLinks"]

# Input/Output service
iosvc = IOSvc()
iosvc.Input = "input_digi.root"
iosvc.Output = "output_standalone_muon_tracks.root"

# Application manager configuration
ApplicationMgr(
    TopAlg=[standaloneMuonTracking],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO
)

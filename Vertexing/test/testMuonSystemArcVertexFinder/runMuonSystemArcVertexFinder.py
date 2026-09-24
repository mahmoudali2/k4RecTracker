"""
Run the muon-system arc vertexer on tracks with muon-chamber hits
(StandaloneMuonTracks) to find LLP decays BEHIND the solenoid.

Author: Mahmoud Althakeel

Example
-------
    k4run runMuonSystemArcVertexFinder.py \
        --inputFile  hnl_tracks.root \
        --outputFile muon_arc_vertices.root

    --loose   write every candidate that is upstream and behind the solenoid,
              with no DCA / chi2 / charge cut (to scan cuts offline); the pair
              DCA, s1, s2 are stored in the vertex parameters
"""
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", required=True)
parser.add_argument("--outputFile", required=True)
parser.add_argument("--inputTracks", default="StandaloneMuonTracks")
parser.add_argument("--loose", action="store_true", help="no DCA/chi2 cut, keep all pairs (cut scans)")
parser.add_argument("--noMaterial", action="store_true",
                    help="no multiple scattering in the track errors (hit resolution + field term only)")
parser.add_argument("--fieldUncertainty", type=float, default=0.5,
                    help="field mismatch [T] in the bending-plane error; 0 = off")
parser.add_argument("--noCurvature", action="store_true", help="straight line in xy instead of a circle")
args = parser.parse_args()

svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

try:
    from Configurables import MuonSystemArcVertexFinder
except ImportError:
    from Vertexing.VertexingConf import MuonSystemArcVertexFinder

finder = MuonSystemArcVertexFinder(
    "MuonSystemArcVertexFinder",
    InputTracks=[args.inputTracks],
    OutputVertices=["MuonSystemVertices"],
    # material between solenoid and muon system, and in the muon-system yoke
    CalorimeterX0=90.0,
    YokeX0=40.0,
    MomentumResolution=0.20,
    FitCurvature=not args.noCurvature,
    UseMaterialEffects=not args.noMaterial,
    FieldUncertainty=args.fieldUncertainty,
    OutputLevel=INFO,
)
if args.loose:
    finder.MaxPairDCA = 1.0e9
    finder.MaxChi2 = 1.0e12
    finder.MinUpstreamDistance = -1.0e9
    finder.RequireOppositeCharge = False
    finder.AllowMultipleVertices = True
    finder.KeepOverlappingPairs = True
    finder.RequireBetweenIPAndHits = False

ApplicationMgr(
    TopAlg=[finder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc")],
    OutputLevel=INFO,
)

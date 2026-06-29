# Example steering for the LeastSquaresVertexFitter Gaudi transformer.
#
# It reads a podio/EDM4hep file that already contains a reconstructed
# edm4hep::TrackCollection and fits those tracks to a common vertex, writing an
# edm4hep::VertexCollection to the output file.
#
# Adapt "Input", the track collection name and (optionally) the beam-spot
# settings to your sample.

from k4FWCore import ApplicationMgr, IOSvc
from Gaudi.Configuration import INFO, DEBUG

io_svc = IOSvc("IOSvc")
io_svc.Input = "tracks_input_edm4hep.root"
io_svc.Output = "vertices_output_edm4hep.root"

from Configurables import LeastSquaresVertexFitter

vertexFitter = LeastSquaresVertexFitter(
    "LeastSquaresVertexFitter",
    InputTracks=["TracksFromGenParticles"],
    OutputVertices=["PrimaryVertices"],
    # which track state to use (1 = AtIP, the perigee state)
    TrackStateLocation=1,
    MinimumNumberOfTracks=2,
    # fit control (defaults shown explicitly for clarity)
    MaxIterations=100,
    ConvergenceThreshold=1.0e-12,
    SeedStartRadius=-1.0,  # negative => start the seed at the perigee
    # Beam-spot (prior-vertex) constraint.
    #   - Enable it for PROMPT vertices: for back-to-back / collinear tracks
    #     (e.g. Z->mumu at rest) the vertex is otherwise unconstrained along the
    #     dimuon axis and scatters by tens of mm; the constraint pins it and the
    #     resolution collapses to a few um.
    #   - Keep it OFF for DISPLACED / secondary vertices (e.g. Bs->mumu), where
    #     it would wrongly pull the vertex toward the origin.
    UseBeamSpotConstraint=True,
    BeamSpotPosition=[0.0, 0.0, 0.0],  # mm
    BeamSpotSize=[0.010, 0.010, 4.0],  # mm (sigma_x, sigma_y, sigma_z), illustrative FCC-ee-like
    MarkAsPrimaryVertex=True,
    OutputLevel=DEBUG,
)

from Configurables import EventDataSvc

ApplicationMgr(
    TopAlg=[vertexFitter],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    StopOnSignal=True,
)

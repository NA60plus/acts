#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

srcdir = Path(__file__).resolve().parent.parent.parent.parent
outputDir = Path.cwd()

field = acts.examples.MagneticFieldMapXyz("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/bfield/NewBFieldNA60plus_longsetup.txt")
field_rotated = acts.examples.MagneticFieldMapXyz("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/bfield/NewBFieldNA60plus_longsetupRotated.txt")

matDeco = acts.IMaterialDecorator.fromFile(
    "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geomVTMSLong/material-map_VT_MSlongsetup.json"
)
jsonFile = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geomVTMSLong/tgeo-config_VT_MSlongsetup.json"
tgeo_fileName = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geomVTMSLong/geom_VT_MSlongsetup.root"
jsonDigi = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geomVTMSLong/digismearVTMS2.json"
jsonSeedVT = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geomVTMSLong/seed_configVTMS.json"
jsonSeedMS = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geomVTMSLong/seed_configMS.json"

logLevel = acts.logging.INFO
customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

detector, trackingGeometry, decorators = acts.examples.TGeoDetector.create(
    jsonFile=str(jsonFile),
    fileName=str(tgeo_fileName),
    surfaceLogLevel=customLogLevel(),
    layerLogLevel=customLogLevel(),
    volumeLogLevel=customLogLevel(),
    mdecorator=matDeco,
)


s = runTruthTrackingKalman(
    trackingGeometry,
    field,
    nev=10000,
    digiConfigFile=jsonDigi,
    outputDir=outputDir,
)

kalmanOptions = {
    "multipleScattering": True,
    "energyLoss": True,
    "reverseFilteringMomThreshold": 0,
    "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
    "level": acts.logging.INFO,
}

s.addAlgorithm(
    acts.examples.RefittingAlgorithm(
        acts.logging.INFO,
        inputTracks="kf_tracks",
        outputTracks="gsf_tracks",
        fit=acts.examples.makeKalmanFitterFunction(trackingGeometry, field, **kalmanOptions),
    )
)

matcher = acts.examples.TrackTruthMatcher(
    level=acts.logging.INFO,
    inputTracks="gsf_tracks",
    inputParticles="particles_selected",
    inputMeasurementParticlesMap="measurement_particles_map",
    outputTrackParticleMatching="track_particle_matching",
    outputParticleTrackMatching="ckf_particle_track",
)
s.addAlgorithm(matcher)

s.addWriter(
    acts.examples.TrackFitterPerformanceWriter(
        level=acts.logging.INFO,
        inputTracks="gsf_tracks",
        inputParticles="particles_input",
        inputTrackParticleMatching="track_particle_matching",
        filePath=str(outputDir / "performance_refitter.root"),
    )
)


s.run()

#!/usr/bin/env python3

from pathlib import Path

import acts
import acts.examples

from acts.examples import TGeoDetector



from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

if "__main__" == __name__:

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

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

    detector, trackingGeometry, decorators = TGeoDetector.create(
        jsonFile=str(jsonFile),
        fileName=str(tgeo_fileName),
        surfaceLogLevel=customLogLevel(),
        layerLogLevel=customLogLevel(),
        volumeLogLevel=customLogLevel(),
        mdecorator=matDeco,
    )


    runTruthTrackingKalman(
        trackingGeometry,
        field,
        nev=10000,
        digiConfigFile=jsonDigi,
        outputDir=Path.cwd(),
    ).run()

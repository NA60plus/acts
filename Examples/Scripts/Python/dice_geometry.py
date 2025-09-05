#!/usr/bin/env python3

import os
import json
from pathlib import Path
from acts.examples import TGeoDetector
from acts import UnitConstants as u
import argparse

import acts
from acts import MaterialMapJsonConverter
from acts.examples import (
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)


def runGeometry(
    trackingGeometry,
    decorators,
    outputDir: Path,
    events=1,
    outputObj=True,
    outputCsv=True,
    outputJson=True,
):
    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0
        ithread = 0

        context = AlgorithmContext(ialg, ievt, eventStore, ithread)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            # if not os.path.isdir(outputDir + "/csv"):
            #    os.makedirs(outputDir + "/csv")
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(outputDir / "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=outputDir / "obj"
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            # if not os.path.isdir(outputDir + "/json"):
            #    os.makedirs(outputDir + "/json")
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(outputDir / "json"),
                writePerEvent=True,
                writeSensitive=True,
            )
            writer.write(context)

            jmConverterCfg = MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                processNonMaterial=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=str(outputDir / "geometry-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)

def buildDICEgeometry(geometryFile = '/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/fullgeo/geometry.root', matDeco = None, addVS = True, addMS = True):
    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet
    equidistant = TGeoDetector.Config.BinningType.equidistant

    volumes = []
    if addVS:
        volumes.append(
            Volume(
                name="VS",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(
                    negative=False,
                    central=False,
                    positive=True),
                subVolumeName=LayerTriplet(
                    positive="VTContainer"),
                sensitiveNames=LayerTriplet(
                    positive=["PixelSensor"]),
                sensitiveAxes=LayerTriplet("XYZ"),  # Single value for all
                rRange=LayerTriplet(
                    positive=(0 * u.mm, 5000 * u.mm)),
                zRange=LayerTriplet(
                    positive=(-100 * u.mm, 400 * u.mm)),
                splitTolR=LayerTriplet(
                    positive=-1.),
                splitTolZ=LayerTriplet(
                    positive=1.),
                binning0=LayerTriplet(
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ))
    if addMS:
        volumes.append(
            Volume(
                name="MS",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(
                    negative=False,
                    central=False,
                    positive=True),
                subVolumeName=LayerTriplet(
                    positive="*"),
                sensitiveNames=LayerTriplet(
                    positive=["MSSensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),  # Single value for all
                rRange=LayerTriplet(
                    positive=(0 * u.mm, 10000 * u.mm)),
                zRange=LayerTriplet(
                    positive=(1000 * u.mm, 10000 * u.mm)),
                splitTolR=LayerTriplet(
                    positive=-1.),
                splitTolZ=LayerTriplet(
                    positive=1.),
                binning0=LayerTriplet(
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            )
            )
        
    config = TGeoDetector.Config(
        fileName=geometryFile,
        surfaceLogLevel=acts.logging.VERBOSE,
        layerLogLevel=acts.logging.VERBOSE,
        volumeLogLevel=acts.logging.VERBOSE,
        buildBeamPipe=False,
        unitScalor=10.0,
        volumes=volumes,
        materialDecorator=matDeco,
    )

    return TGeoDetector(config)
        

if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="Run DICE geometry simulation.")
    parser.add_argument("--geometry-file", type=str, default='/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/fullgeo/geometry.root', help="Path to the geometry ROOT file.")
    parser.add_argument("--remove-vs", action="store_true", help="Remove Vertex Spectrometer")
    parser.add_argument("--remove-ms", action="store_true", help="Remove Muon Spectrometer")
    args = parser.parse_args()
    
    detector = buildDICEgeometry(geometryFile=args.geometry_file,
                                 addVS=not args.remove_vs, addMS=not args.remove_ms)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    print("Tracking geometry built successfully.")

    runGeometry(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        outputDir=Path(""),
        events=1,
        outputObj=True,
        outputCsv=True,
        outputJson=True,
    )
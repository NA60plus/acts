#!/usr/bin/env python3

import os
import argparse
from dice_geometry import buildDICEgeometry
import acts
from acts import (
    SurfaceMaterialMapper,
    VolumeMaterialMapper,
    Navigator,
    Propagator,
    StraightLineStepper,
    MaterialMapJsonConverter,
)
from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    MaterialMapping,
    JsonMaterialWriter,
    JsonFormat,
)


def runMaterialMapping(
    trackingGeometry,
    decorators,
    outputDir,
    inputDir,
    mapName="material-map",
    mapFormat=JsonFormat.Json,
    mapSurface=True,
    mapVolume=True,
    readCachedSurfaceInformation=False,
    mappingStep=1,
    s=None,
    nevents=1000,
):
    s = s or Sequencer(numThreads=1,events=nevents)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    wb = WhiteBoard(acts.logging.INFO)

    context = AlgorithmContext(0, 0, wb, 0)

    for decorator in decorators:
        assert decorator.decorate(context) == ProcessCode.SUCCESS

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material-tracks",
            fileList=[
                os.path.join(
                    inputDir,
                    (
                        mapName + "_tracks.root"
                        if readCachedSurfaceInformation
                        else "geant4_material_tracks.root"
                    ),
                )
            ],
            readCachedSurfaceInformation=readCachedSurfaceInformation,
        )
    )

    stepper = StraightLineStepper()

    mmAlgCfg = MaterialMapping.Config(context.geoContext, context.magFieldContext)
    mmAlgCfg.trackingGeometry = trackingGeometry
    mmAlgCfg.inputMaterialTracks = "material-tracks"

    if mapSurface:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
            resolveSensitive=True,
            resolveMaterial=True,
            resolvePassive=True,
        )
        propagator = Propagator(stepper, navigator)
        mapper = SurfaceMaterialMapper(level=acts.logging.INFO, propagator=propagator)
        mmAlgCfg.materialSurfaceMapper = mapper

    if mapVolume:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
        )
        propagator = Propagator(stepper, navigator)
        mapper = VolumeMaterialMapper(
            level=acts.logging.INFO, propagator=propagator, mappingStep=mappingStep
        )
        mmAlgCfg.materialVolumeMapper = mapper

    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=True,
        processApproaches=True,
        processRepresenting=True,
        processBoundaries=True,
        processVolumes=True,
        context=context.geoContext,
    )

    jmw = JsonMaterialWriter(
        level=acts.logging.VERBOSE,
        converterCfg=jmConverterCfg,
        fileName=os.path.join(outputDir, mapName),
        writeFormat=mapFormat,
    )

    mmAlgCfg.materialWriters = [jmw]

    s.addAlgorithm(MaterialMapping(level=acts.logging.INFO, config=mmAlgCfg))

    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=mmAlgCfg.mappingMaterialCollection,
            filePath=os.path.join(
                outputDir,
                mapName + "_tracks.root",
            ),
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s


if "__main__" == __name__:
    parser = argparse.ArgumentParser(description="Script to generate ACTS material map")
    parser.add_argument("--geometry-file", type=str, default='/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/fullgeo/geometry.root', help="Path to the geometry ROOT file.")

    parser.add_argument(
        "-i",
        "--inFile",
        type=str,
        default="geometry-map.json",
        help="Output filename for the generated material map. Supported formats: JSON, CBOR.",
    )
    parser.add_argument(
        "-o",
        "--outFile",
        type=str,
        default="material-map.json",
        help="Output filename for the generated material map. Supported formats: JSON, CBOR.",
    )
    
    parser.add_argument("--remove-vs", action="store_true", help="Remove Vertex Spectrometer")
    parser.add_argument("--remove-ms", action="store_true", help="Remove Muon Spectrometer")
    parser.add_argument("-n", "--nevents", type=int, default=1, help="Number of events to process")
    
    args = parser.parse_args()

    mapName = args.outFile.split(".")[0]
    matDeco = acts.IMaterialDecorator.fromFile(args.inFile)

    detector = buildDICEgeometry(geometryFile=args.geometry_file,
                                matDeco=matDeco, addVS=not args.remove_vs, addMS=not args.remove_ms)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=os.getcwd(),
        inputDir=os.getcwd(),
        readCachedSurfaceInformation=False,
        mapName=mapName,
        mapFormat=JsonFormat.Json,
        nevents=args.nevents,
    ).run()

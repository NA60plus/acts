#!/usr/bin/env python3

import argparse
from dice_geometry import buildDICEgeometry
import os
import acts
from acts import (
    MaterialMapper,
    IntersectionMaterialAssigner,
    BinnedSurfaceMaterialAccumulater,
    MaterialMapJsonConverter,
    logging,
    GeometryContext,
    DetectorBuilder,
    GeometryIdGenerator,
)

from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    CoreMaterialMapping,
    JsonMaterialWriter,
    RootMaterialWriter,
    JsonFormat,
)


def runMaterialMapping(surfaces, inputFile, outputFile, outputMap, loglevel):
    # Create a sequencer
    print("Creating the sequencer with 1 thread (inter event information needed)")

    s = Sequencer(numThreads=1)

    # IO for material tracks reading
    wb = WhiteBoard(acts.logging.INFO)

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material-tracks",
            fileList=[inputFile],
            readCachedSurfaceInformation=False,
        )
    )

    # Assignment setup : Intersection assigner
    materialAssingerConfig = IntersectionMaterialAssigner.Config()
    materialAssingerConfig.surfaces = surfaces
    materialAssinger = IntersectionMaterialAssigner(materialAssingerConfig, loglevel)

    # Accumulation setup : Binned surface material accumulater
    materialAccumulaterConfig = BinnedSurfaceMaterialAccumulater.Config()
    materialAccumulaterConfig.materialSurfaces = surfaces
    materialAccumulater = BinnedSurfaceMaterialAccumulater(
        materialAccumulaterConfig, loglevel
    )

    # Mapper setup
    materialMapperConfig = MaterialMapper.Config()
    materialMapperConfig.assignmentFinder = materialAssinger
    materialMapperConfig.surfaceMaterialAccumulater = materialAccumulater
    materialMapper = MaterialMapper(materialMapperConfig, loglevel)

    # Add the map writer(s)
    mapWriters = []
    # json map writer
    context = AlgorithmContext(0, 0, wb, 0)
    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=True,
        processApproaches=True,
        processRepresenting=True,
        processBoundaries=False,
        processVolumes=True,
        processNonMaterial=True,
        processDenseVolumes=True,
        context=context.geoContext,
    )
    mapWriters.append(
        JsonMaterialWriter(
            level=loglevel,
            converterCfg=jmConverterCfg,
            fileName=outputMap + "",
            writeFormat=JsonFormat.Json,
        )
    )
    mapWriters.append(RootMaterialWriter(level=loglevel, filePath=outputMap + ".root"))

    # Mapping Algorithm
    coreMaterialMappingConfig = CoreMaterialMapping.Config()
    coreMaterialMappingConfig.materialMapper = materialMapper
    coreMaterialMappingConfig.inputMaterialTracks = "material-tracks"
    coreMaterialMappingConfig.mappedMaterialTracks = "mapped-material-tracks"
    coreMaterialMappingConfig.unmappedMaterialTracks = "unmapped-material-tracks"
    coreMaterialMappingConfig.materiaMaplWriters = mapWriters
    coreMaterialMapping = CoreMaterialMapping(coreMaterialMappingConfig, loglevel)
    s.addAlgorithm(coreMaterialMapping)

    # Add the mapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks="mapped-material-tracks",
            filePath=outputFile + "_mapped.root",
            storeSurface=True,
            storeVolume=True,
        )
    )

    # Add the unmapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks="unmapped-material-tracks",
            filePath=outputFile + "_unmapped.root",
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
    materialSurfaces = trackingGeometry.extractMaterialSurfaces()


    outputDir=os.path.join(
                os.getcwd(),
                (
                    "material-map_tracks"
                ),
            )
    inputDir=os.path.join(
                os.getcwd(),
                (
                    "geant4_material_tracks.root"
                ),
            )

    gContext = GeometryContext()
    logLevel = logging.INFO
    runMaterialMapping(
        materialSurfaces, inputDir, outputDir, mapName, logLevel
    ).run()

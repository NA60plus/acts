from pathlib import Path
from typing import Optional, Union, List
from enum import Enum
from collections import namedtuple

import acts
import acts.examples

u = acts.UnitConstants

SeedingAlgorithm = Enum(
    "SeedingAlgorithm",
    "Default TruthSmeared TruthEstimated Orthogonal HoughTransform Gbts",
)

TruthSeedRanges = namedtuple(
    "TruthSeedRanges",
    ["rho", "z", "phi", "eta", "absEta", "pt", "nHits","keep"],
    defaults=[(None, None)] * 8,
)

ParticleSmearingSigmas = namedtuple(
    "ParticleSmearingSigmas",
    ["d0", "d0PtA", "d0PtB", "z0", "z0PtA", "z0PtB", "t0", "phi", "theta", "pRel"],
    defaults=[None] * 10,
)

SeedFinderConfigArg = namedtuple(
    "SeedFinderConfig",
    [
        "maxSeedsPerSpM",
        "cotThetaMax",
        "sigmaScattering",
        "radLengthPerSeed",
        "minPt",
        "impactMax",
        "deltaPhiMax",
        "interactionPointCut",
        "deltaZMax",
        "maxPtScattering",
        "zBinEdges",
        "zBinsCustomLooping",
        "rRangeMiddleSP",
        "useVariableMiddleSPRange",
        "binSizeR",
        "seedConfirmation",
        "centralSeedConfirmationRange",
        "forwardSeedConfirmationRange",
        "verbose",
        "deltaR",  # (min,max)
        "deltaRBottomSP",  # (min,max)
        "deltaRTopSP",  # (min,max)
        "deltaRMiddleSPRange",  # (min,max)
        "collisionRegion",  # (min,max)
        "r",  # (min,max)
        "z",  # (min,max)
        "rMiddle"
    ],
    defaults=[None] * 19 + [(None, None)] * 8,
)
SeedFinderOptionsArg = namedtuple(
    "SeedFinderOptions", ["beamPos", "bFieldInZ"], defaults=[(None, None), None]
)

SeedFilterConfigArg = namedtuple(
    "SeedFilterConfig",
    [
        "impactWeightFactor",
        "zOriginWeightFactor",
        "compatSeedWeight",
        "compatSeedLimit",
        "numSeedIncrement",
        "seedWeightIncrement",
        "seedConfirmation",
        "maxSeedsPerSpMConf",
        "maxQualitySeedsPerSpMConf",
        "useDeltaRorTopRadius",
        "deltaRMin",
        "verbose"
    ],
    defaults=[None] * 12,
)


SeedFinderConfigArgNA60 = namedtuple(
    "SeedFinderConfigNA60",
    [
        "maxSeedsPerSpM",
        "cotThetaMax",
        "sigmaScattering",
        "radLengthPerSeed",
        "minPt",
        "impactMax",
        "deltaPhiMax",
        "interactionPointCut",
        "deltaZMax",
        "maxPtScattering",
        "zBinEdges",
        "zBinsCustomLooping",
        "skipZMiddleBinSearch",
        "rRangeMiddleSP",
        "useVariableMiddleSPRange",
        "binSizeR",
        "seedConfirmation",
        "seedConfirmationRange",
        "verbose", #added by me
        "deltaY",  # (min,max)
        "deltaYBottomSP",  # (min,max)
        "deltaYTopSP",  # (min,max)
        "deltaYMiddleSPRange",  # (min,max)
        "collisionRegion",  # (min,max)
        "y",  # (min,max)
        "z",  # (min,max)
        "zOutermostLayers",  # (min,max)
        "yMiddle" #added by me
    ],
    defaults=[None] * 19 + [(None, None)] * 9,
)

SeedFinderOptionsArgNA60 = namedtuple(
    "SeedFinderOptionsNA60", ["beamPos", "bFieldInZ"], defaults=[(None, None), None]
)

SeedFilterConfigArgNA60 = namedtuple(
    "SeedFilterConfigNA60",
    [
        "impactWeightFactor",
        "zOriginWeightFactor",
        "compatSeedWeight",
        "compatSeedLimit",
        "numSeedIncrement",
        "seedWeightIncrement",
        "seedConfirmation",
        "maxSeedsPerSpMConf",
        "maxQualitySeedsPerSpMConf",
        "verbose", #added by me
        "useDeltaRorTopRadius",
        "deltaYMin", 
    ],
    defaults=[None] * 12,
)




SpacePointGridConfigArg = namedtuple(
    "SeedGridConfig",
    [
        "rMax",
        "zBinEdges",
        "phiBinDeflectionCoverage",
        "impactMax",
        "deltaRMax",
        "maxPhiBins",
        "phi",  # (min,max)
    ],
    defaults=[None] * 6 + [(None, None)] * 1,
)

SeedingAlgorithmConfigArg = namedtuple(
    "SeedingAlgorithmConfig",
    [
        "allowSeparateRMax",
        "zBinNeighborsTop",
        "zBinNeighborsBottom",
        "numPhiNeighbors",
        "useExtraCuts",
    ],
    defaults=[None] * 5,
)

TruthEstimatedSeedingAlgorithmConfigArg = namedtuple(
    "TruthSeederConfig",
    [
        "deltaR",  # (min,max)
    ],
    defaults=[(None, None)],
)

TrackSelectorConfig = namedtuple(
    "TrackSelectorConfig",
    [
        "loc0",
        "loc1",
        "time",
        "eta",
        "absEta",
        "pt",
        "phi",
        "nMeasurementsMin",
        "maxHoles",
        "maxOutliers",
        "maxSharedHits",
        "maxChi2",
        "nMeasurementsGroupMin",
    ],
    defaults=[(None, None)] * 7 + [None] * 6,
)

CkfConfig = namedtuple(
    "CkfConfig",
    [
        "chi2CutOff",
        "numMeasurementsCutOff",
        "maxSteps",
        "seedDeduplication",
        "stayOnSeed",
        "pixelVolumes",
        "stripVolumes",
        "maxPixelHoles",
        "maxStripHoles",
    ],
    defaults=[15.0, 10, None, None, None, None, None, None, None],
)

AmbiguityResolutionConfig = namedtuple(
    "AmbiguityResolutionConfig",
    ["maximumSharedHits", "nMeasurementsMin", "maximumIterations"],
    defaults=[None] * 3,
)

ScoreBasedAmbiguityResolutionConfig = namedtuple(
    "ScoreBasedAmbiguityResolutionConfig",
    [
        "minScore",
        "minScoreSharedTracks",
        "maxShared",
        "maxSharedTracksPerMeasurement",
        "pTMax",
        "pTMin",
        "phiMax",
        "phiMin",
        "etaMax",
        "etaMin",
        "useAmbiguityFunction",
    ],
    defaults=[None] * 11,
)

AmbiguityResolutionMLConfig = namedtuple(
    "AmbiguityResolutionMLConfig",
    ["maximumSharedHits", "nMeasurementsMin", "maximumIterations"],
    defaults=[None] * 3,
)

AmbiguityResolutionMLDBScanConfig = namedtuple(
    "AmbiguityResolutionMLDBScanConfig",
    ["nMeasurementsMin", "epsilonDBScan", "minPointsDBScan"],
    defaults=[None] * 3,
)

SeedFilterMLDBScanConfig = namedtuple(
    "SeedFilterMLDBScanConfig",
    ["epsilonDBScan", "minPointsDBScan", "minSeedScore"],
    defaults=[None] * 3,
)


class VertexFinder(Enum):
    Truth = (1,)
    AMVF = (2,)
    Iterative = (3,)


@acts.examples.NamedTypeArgs(
    seedingAlgorithm=SeedingAlgorithm,
    truthSeedRanges=TruthSeedRanges,
    particleSmearingSigmas=ParticleSmearingSigmas,
    seedFinderConfigArg=SeedFinderConfigArg,
    seedFinderOptionsArg=SeedFinderOptionsArg,
    seedFilterConfigArg=SeedFilterConfigArg,
    spacePointGridConfigArg=SpacePointGridConfigArg,
    seedingAlgorithmConfigArg=SeedingAlgorithmConfigArg,
    truthEstimatedSeedingAlgorithmConfigArg=TruthEstimatedSeedingAlgorithmConfigArg,
    logLevel=acts.logging.Level,
)
def addSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Optional[Union[Path, str]] = None,
    layerMappingConfigFile: Optional[Union[Path, str]] = None,
    connector_inputConfigFile: Optional[Union[Path, str]] = None,
    seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.Default,
    truthSeedRanges: Optional[TruthSeedRanges] = TruthSeedRanges(),
    particleSmearingSigmas: ParticleSmearingSigmas = ParticleSmearingSigmas(),
    initialSigmas: Optional[list] = None,
    initialVarInflation: Optional[list] = None,
    seedFinderConfigArg: SeedFinderConfigArg = SeedFinderConfigArg(),
    seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(),
    seedFilterConfigArg: SeedFilterConfigArg = SeedFilterConfigArg(),
    spacePointGridConfigArg: SpacePointGridConfigArg = SpacePointGridConfigArg(),
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(),
    houghTransformConfig: acts.examples.HoughTransformSeeder.Config = acts.examples.HoughTransformSeeder.Config(),
    truthEstimatedSeedingAlgorithmConfigArg: TruthEstimatedSeedingAlgorithmConfigArg = TruthEstimatedSeedingAlgorithmConfigArg(),
    particleHypothesis: Optional[
        acts.ParticleHypothesis
    ] = acts.ParticleHypothesis.pion,
    inputParticles: str = "particles",
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
    verbose: bool = False,
    inputSourceLinks="sourcelinks",
    inputMeasurements="measurements",
    doTrkVtx=True,
    noGuessing = False,
    suffix = "",
    trkVtxOnly = False,
    addDeltas = True,
    projective = False,
    zPerigee = 0,
    det_suffix = "",
    suffixSpacepoint = ""
) -> None:
    """This function steers the seeding
    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection. Not required for SeedingAlgorithm.TruthSmeared.
    seedingAlgorithm : SeedingAlgorithm, Default
        seeding algorithm to use: one of Default (no truth information used), TruthSmeared, TruthEstimated
    truthSeedRanges : TruthSeedRanges(rho, z, phi, eta, absEta, pt, nHits)
        TruthSeedSelector configuration. Each range is specified as a tuple of (min,max).
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TruthSeedSelector.hpp
        If specified as None, don't run ParticleSmearing at all (and use addCKFTracks(selectedParticles="particles"))
    particleSmearingSigmas : ParticleSmearingSigmas(d0, d0PtA, d0PtB, z0, z0PtA, z0PtB, t0, phi, theta, pRel)
        ParticleSmearing configuration.
        Defaults specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSmearing.hpp
    initialSigmas : list
        Sets the initial covariance matrix diagonal. This is ignored in case of TruthSmearing.
        Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp
    initialVarInflation : list
        List of 6 scale factors to inflate the initial covariance matrix
        Defaults (all 1) specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSmearing.hpp
    seedFinderConfigArg : SeedFinderConfigArg(maxSeedsPerSpM, cotThetaMax, sigmaScattering, radLengthPerSeed, minPt, impactMax, deltaPhiMax, interactionPointCut, deltaZMax, maxPtScattering, zBinEdges, zBinsCustomLooping, rRangeMiddleSP, useVariableMiddleSPRange, binSizeR, seedConfirmation, centralSeedConfirmationRange, forwardSeedConfirmationRange, deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z)
        SeedFinderConfig settings. deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z.
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFinderOptionsArg :  SeedFinderOptionsArg(bFieldInZ, beamPos)
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFilterConfigArg : SeedFilterConfigArg(compatSeedWeight, compatSeedLimit, numSeedIncrement, seedWeightIncrement, seedConfirmation, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf, useDeltaRorTopRadius)
                                Defaults specified in Core/include/Acts/Seeding/SeedFilterConfig.hpp
    spacePointGridConfigArg : SpacePointGridConfigArg(rMax, zBinEdges, phiBinDeflectionCoverage, phi, maxPhiBins, impactMax)
                                SpacePointGridConfigArg settings. phi is specified as a tuple of (min,max).
        Defaults specified in Core/include/Acts/Seeding/SpacePointGrid.hpp
    seedingAlgorithmConfigArg : SeedingAlgorithmConfigArg(allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors, useExtraCuts)
                                Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/SeedingAlgorithm.hpp
    truthEstimatedSeedingAlgorithmConfigArg : TruthEstimatedSeedingAlgorithmConfigArg(deltaR)
        Currently only deltaR=(min,max) range specified here.
    particleHypothesis : Optional[acts.ParticleHypothesis]
        The hypothesis used for track finding. Defaults to pion.
    inputParticles : str, "particles"
        input particles name in the WhiteBoard
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    rnd : RandomNumbers, None
        random number generator. Only used by SeedingAlgorithm.TruthSmeared.
    """

    logLevel = acts.examples.defaultLogging(s, logLevel)()
    logger = acts.logging.getLogger("addSeeding")
    logger.setLevel(logLevel)

    if truthSeedRanges is not None:
        selectedParticles = "truth_seeds_selected"+det_suffix
        addSeedingTruthSelection(
            s,
            inputParticles,
            selectedParticles,
            truthSeedRanges,
            logLevel,
            det_suffix
        )
    else:
        selectedParticles = inputParticles

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        logger.info("Using smeared truth particles for seeding")
        addTruthSmearedSeeding(
            s,
            rnd,
            selectedParticles,
            particleSmearingSigmas,
            initialSigmas,
            initialVarInflation,
            particleHypothesis,
            logLevel,
        )
    else:
        spacePoints = addSpacePointsMaking(
            s,
            trackingGeometry,
            geoSelectionConfigFile,
            logLevel,
            inputSourceLinks=inputSourceLinks,
            inputMeasurements=inputMeasurements,
            suffix=suffixSpacepoint
        )

        if doTrkVtx :
            addTrackletVertexing(
                s,
                inputSpacePoints=spacePoints,
                noGuessing=noGuessing,
                verbose=verbose,
                useFit=True,
                nbins=60,
                addDeltas=addDeltas,
                projective=projective,
                zPerigee=zPerigee
            )

        if trkVtxOnly:
            return s
            
        # Run either: truth track finding or seeding
        if seedingAlgorithm == SeedingAlgorithm.TruthEstimated:
            logger.info("Using truth track finding from space points for seeding")
            seeds = addTruthEstimatedSeeding(
                s,
                spacePoints,
                selectedParticles,
                truthEstimatedSeedingAlgorithmConfigArg,
                logLevel,
                suffixSeed=suffix
            )
        elif seedingAlgorithm == SeedingAlgorithm.Default:
            logger.info("Using default seeding")
            seeds = addStandardSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
                suffixSeed=suffix
            )
        elif seedingAlgorithm == SeedingAlgorithm.Orthogonal:
            logger.info("Using orthogonal seeding")
            seeds = addOrthogonalSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.HoughTransform:
            logger.info("Using Hough Transform seeding")
            houghTransformConfig.inputSpacePoints = [spacePoints]
            houghTransformConfig.inputMeasurements = "measurements"+det_suffix
            houghTransformConfig.inputSourceLinks = "sourcelinks"+det_suffix
            houghTransformConfig.outputProtoTracks = "prototracks"+det_suffix
            houghTransformConfig.outputSeeds = "seeds"
            houghTransformConfig.trackingGeometry = trackingGeometry
            seeds = addHoughTransformSeeding(s, houghTransformConfig, logLevel)
        elif seedingAlgorithm == SeedingAlgorithm.Gbts:
            logger.info("Using Gbts seeding")
            # output of algs changed, only one output now
            seeds = addGbtsSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                trackingGeometry,
                logLevel,
                layerMappingConfigFile,
                geoSelectionConfigFile,
                connector_inputConfigFile,
            )
        else:
            logger.fatal("unknown seedingAlgorithm %s", seedingAlgorithm)

        parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
            level=logLevel,
            inputSeeds=seeds,
            outputTrackParameters=suffix+"estimatedparameters",
            outputSeeds=suffix+"estimatedseeds",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            #verbose=True,
            **acts.examples.defaultKWArgs(
                initialSigmas=initialSigmas,
                initialVarInflation=initialVarInflation,
                particleHypothesis=particleHypothesis,
                #verbose=True
            ),
        )
        s.addAlgorithm(parEstimateAlg)

        prototracks = suffix+"seed-prototracks"
        s.addAlgorithm(
            acts.examples.SeedsToPrototracks(
                level=logLevel,
                inputSeeds=seeds,
                outputProtoTracks=prototracks
            )
        )

        if outputDirRoot is not None:
            addSeedPerformanceWriters(
                s,
                outputDirRoot,
                seeds,
                prototracks,
                selectedParticles,
                inputParticles,
                parEstimateAlg.config.outputTrackParameters,
                logLevel,
                suffix,
                det_suffix,
                noGuessing
            )

        if outputDirCsv is not None:
            outputDirCsv = Path(outputDirCsv)

            if not outputDirCsv.exists():
                outputDirCsv.mkdir()

            csvSeedWriter = acts.examples.CsvSeedWriter(
                level=logLevel,
                inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
                inputSimSeeds=seeds,
                inputSimHits="simhits"+det_suffix,
                inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
                inputMeasurementSimHitsMap="measurement_simhits_map"+det_suffix,
                outputDir=str(outputDirCsv),
                fileName=str(f"seed.csv"),
            )
            s.addWriter(csvSeedWriter)

    return s


@acts.examples.NamedTypeArgs(
    seedingAlgorithm=SeedingAlgorithm,
    truthSeedRanges=TruthSeedRanges,
    particleSmearingSigmas=ParticleSmearingSigmas,
    seedFinderConfigArg=SeedFinderConfigArg,
    seedFinderOptionsArg=SeedFinderOptionsArg,
    seedFilterConfigArg=SeedFilterConfigArg,
    spacePointGridConfigArg=SpacePointGridConfigArg,
    seedingAlgorithmConfigArg=SeedingAlgorithmConfigArg,
    truthEstimatedSeedingAlgorithmConfigArg=TruthEstimatedSeedingAlgorithmConfigArg,
    logLevel=acts.logging.Level,
)
def addNewSeeding(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    geoSelectionConfigFile: Optional[Union[Path, str]] = None,
    truthSeedRanges: Optional[TruthSeedRanges] = None,
    initialSigmas: Optional[list] = None,
    initialVarInflation: Optional[list] = None,
    seedFinderConfigArg: SeedFinderConfigArgNA60 = SeedFinderConfigArgNA60(),
    seedFinderOptionsArg: SeedFinderOptionsArgNA60 = SeedFinderOptionsArgNA60(),
    seedFilterConfigArg: SeedFilterConfigArgNA60 = SeedFilterConfigArgNA60(),
    spacePointGridConfigArg: SpacePointGridConfigArg = SpacePointGridConfigArg(),
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(),
    particleHypothesis: Optional[
        acts.ParticleHypothesis
    ] = acts.ParticleHypothesis.pion,
    verbose: bool = False,
    inputParticles: str = "particles",
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    inputSourceLinks="sourcelinks",
    inputMeasurements="measurements",
    doTrkVtx=True,
    noGuessing = False,
    trkVtxOnly=False,
    addDeltas=True,
    projective=False,
    zPerigee = 0,
    suffix = "",
    det_suffix = ""
) -> None:

    logLevel = acts.examples.defaultLogging(s, logLevel)()
    logger = acts.logging.getLogger("addSeeding")
    logger.setLevel(logLevel)

    if truthSeedRanges is not None:
        selectedParticles = "truth_seeds_selected"+det_suffix
        addSeedingTruthSelection(
            s,
            inputParticles,
            selectedParticles,
            truthSeedRanges,
            logLevel,
            det_suffix
        )
    else:
        selectedParticles = inputParticles


    spacePoints = addSpacePointsMaking(
        s,
        trackingGeometry,
        geoSelectionConfigFile,
        logLevel,
        inputSourceLinks=inputSourceLinks,
        inputMeasurements=inputMeasurements,
        suffix=suffix
    )
    if doTrkVtx :
        addTrackletVertexing(
            s,
            noGuessing=noGuessing,
            addDeltas=addDeltas,
            projective=projective,
            zPerigee = zPerigee,
            inputSpacePoints=spacePoints,#+det_suffix,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
        )

    if trkVtxOnly:
        return s
    
    logger.info("Using default seeding")
    seeds = addStandardSeedingNA60(
        s,
        spacePoints,
        seedingAlgorithmConfigArg,
        seedFinderConfigArg,
        seedFinderOptionsArg,
        seedFilterConfigArg,
        spacePointGridConfigArg,
        logLevel,
        suffixSeed=suffix
    )

    parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
        level=logLevel,
        inputSeeds=seeds,
        outputTrackParameters=suffix+"estimatedparameters",
        outputSeeds=suffix+"estimatedseeds",
        trackingGeometry=trackingGeometry,
        magneticField=field,
        **acts.examples.defaultKWArgs(
            initialSigmas=initialSigmas,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis
        ),
    )
    s.addAlgorithm(parEstimateAlg)

    prototracks = suffix+"seed-prototracks"
    s.addAlgorithm(
        acts.examples.SeedsToPrototracks(
            level=logLevel,
            inputSeeds=seeds,
            outputProtoTracks=prototracks
        )
    )

    if outputDirRoot is not None:
        addSeedPerformanceWriters(
            s,
            outputDirRoot,
            seeds,
            prototracks,
            selectedParticles,
            inputParticles,
            parEstimateAlg.config.outputTrackParameters,
            logLevel,
            suffix,
            det_suffix,
            noGuessing
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)

        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        csvSeedWriter = acts.examples.CsvSeedWriter(
            level=logLevel,
            inputTrackParameters=parEstimateAlg.config.outputTrackParameters,
            inputSimSeeds=seeds,
            inputSimHits="simhits"+det_suffix,
            inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
            inputMeasurementSimHitsMap="measurement_simhits_map"+det_suffix,
            outputDir=str(outputDirCsv),
            fileName=str(f"seed.csv"),
        )
        s.addWriter(csvSeedWriter)


def addSeedingTruthSelection(
    s: acts.examples.Sequencer,
    inputParticles: str,
    outputParticles: str,
    truthSeedRanges: TruthSeedRanges,
    logLevel: acts.logging.Level = None,
    det_suffix = ""
):
    """adds truth particles filtering before filtering
    For parameters description see addSeeding
    """
    selAlg = acts.examples.TruthSeedSelector(
        **acts.examples.defaultKWArgs(
            ptMin=truthSeedRanges.pt[0],
            ptMax=truthSeedRanges.pt[1],
            etaMin=truthSeedRanges.eta[0],
            etaMax=truthSeedRanges.eta[1],
            nHitsMin=truthSeedRanges.nHits[0],
            nHitsMax=truthSeedRanges.nHits[1],
            rhoMin=truthSeedRanges.rho[0],
            rhoMax=truthSeedRanges.rho[1],
            zMin=truthSeedRanges.z[0],
            zMax=truthSeedRanges.z[1],
            phiMin=truthSeedRanges.phi[0],
            phiMax=truthSeedRanges.phi[1],
            absEtaMin=truthSeedRanges.absEta[0],
            absEtaMax=truthSeedRanges.absEta[1],
            keepPrimary=truthSeedRanges.keep[0],
            keepSecondary=truthSeedRanges.keep[1],
        ),
        level=logLevel,
        inputParticles=inputParticles,
        inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
        outputParticles=outputParticles,
    )
    s.addAlgorithm(selAlg)


def addTruthSmearedSeeding(
    s: acts.examples.Sequencer,
    rnd: Optional[acts.examples.RandomNumbers],
    selectedParticles: str,
    particleSmearingSigmas: ParticleSmearingSigmas,
    initialSigmas: Optional[List[float]],
    initialVarInflation: List[float],
    particleHypothesis: Optional[acts.ParticleHypothesis],
    logLevel: acts.logging.Level = None,
):
    """adds algorithm that would mimic detector response uncertainties for truth seeding
    For parameters description see addSeeding
    """

    rnd = rnd or acts.examples.RandomNumbers(seed=42)
    # Run particle smearing
    ptclSmear = acts.examples.ParticleSmearing(
        level=logLevel,
        inputParticles=selectedParticles,
        outputTrackParameters="estimatedparameters",
        randomNumbers=rnd,
        # gaussian sigmas to smear particle parameters
        **acts.examples.defaultKWArgs(
            sigmaD0=particleSmearingSigmas.d0,
            sigmaD0PtA=particleSmearingSigmas.d0PtA,
            sigmaD0PtB=particleSmearingSigmas.d0PtB,
            sigmaZ0=particleSmearingSigmas.z0,
            sigmaZ0PtA=particleSmearingSigmas.z0PtA,
            sigmaZ0PtB=particleSmearingSigmas.z0PtB,
            sigmaT0=particleSmearingSigmas.t0,
            sigmaPhi=particleSmearingSigmas.phi,
            sigmaTheta=particleSmearingSigmas.theta,
            sigmaPRel=particleSmearingSigmas.pRel,
            initialSigmas=initialSigmas,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis,
        ),
    )
    s.addAlgorithm(ptclSmear)

    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=logLevel,
        inputParticles=selectedParticles,
        inputMeasurementParticlesMap="measurement_particles_map",
        outputProtoTracks="truth_particle_tracks",
    )
    s.addAlgorithm(truthTrkFndAlg)


def addTruthEstimatedSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    inputParticles: str,
    TruthEstimatedSeedingAlgorithmConfigArg: TruthEstimatedSeedingAlgorithmConfigArg,
    logLevel: acts.logging.Level = None,
    suffixSeed=""
):
    """adds truth seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    truthSeeding = acts.examples.TruthSeedingAlgorithm(
        level=logLevel,
        inputParticles=inputParticles,
        inputMeasurementParticlesMap="measurement_particles_map",
        inputSpacePoints=[spacePoints],
        outputParticles="truth_seeded_particles",
        outputProtoTracks="truth_particle_tracks",
        outputSeeds=suffixSeed+"seeds",
        **acts.examples.defaultKWArgs(
            deltaRMin=TruthEstimatedSeedingAlgorithmConfigArg.deltaR[0],
            deltaRMax=TruthEstimatedSeedingAlgorithmConfigArg.deltaR[1],
        ),
    )
    sequence.addAlgorithm(truthSeeding)

    return truthSeeding.config.outputSeeds


def addSpacePointsMaking(
    sequence: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geoSelectionConfigFile: Union[Path, str],
    logLevel: acts.logging.Level = None,
    inputSourceLinks="sourcelinks",
    inputMeasurements="measurements",
    suffix = ""
):
    """adds space points making
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputSourceLinks=inputSourceLinks,
        inputMeasurements=inputMeasurements,
        outputSpacePoints=suffix+"spacepoints",
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
    )
    sequence.addAlgorithm(spAlg)
    return spAlg.config.outputSpacePoints


def addStandardSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    logLevel: acts.logging.Level = None,
    outputPrimaryVertex: str="OutputRecPrimaryVertex",
    suffixSeed= ""
):
    """adds standard seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()

    seedFinderConfig = acts.SeedFinderConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMin=seedFinderConfigArg.deltaR[0],
            deltaRMax=seedFinderConfigArg.deltaR[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            deltaRMiddleMinSPRange=seedFinderConfigArg.deltaRMiddleSPRange[0],
            deltaRMiddleMaxSPRange=seedFinderConfigArg.deltaRMiddleSPRange[1],
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            zBinEdges=seedFinderConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            binSizeR=seedFinderConfigArg.binSizeR,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
            rMinMiddle = seedFinderConfigArg.rMiddle[0],
            rMaxMiddle = seedFinderConfigArg.rMiddle[1],
            verbose =  seedFinderConfigArg.verbose            
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=seedFinderConfig.deltaRMin,
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfig.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfig.forwardSeedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
            verbose= seedFilterConfigArg.verbose
        )
    )

    gridConfig = acts.SpacePointGridConfig(
        **acts.examples.defaultKWArgs(
            minPt=seedFinderConfig.minPt,
            rMax=(
                seedFinderConfig.rMax
                if spacePointGridConfigArg.rMax is None
                else spacePointGridConfigArg.rMax
            ),
            zMax=seedFinderConfig.zMax,
            zMin=seedFinderConfig.zMin,
            deltaRMax=(
                seedFinderConfig.deltaRMax
                if spacePointGridConfigArg.deltaRMax is None
                else spacePointGridConfigArg.deltaRMax
            ),
            cotThetaMax=seedFinderConfig.cotThetaMax,
            phiMin=spacePointGridConfigArg.phi[0],
            phiMax=spacePointGridConfigArg.phi[1],
            maxPhiBins=spacePointGridConfigArg.maxPhiBins,
            impactMax=spacePointGridConfigArg.impactMax,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            phiBinDeflectionCoverage=spacePointGridConfigArg.phiBinDeflectionCoverage,
        )
    )

    gridOptions = acts.SpacePointGridOptions(
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptions.bFieldInZ,
        )
    )

    seedingAlg = acts.examples.SeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds=suffixSeed+"seeds",
        **acts.examples.defaultKWArgs(
            allowSeparateRMax=seedingAlgorithmConfigArg.allowSeparateRMax,
            zBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            zBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom,
            numPhiNeighbors=seedingAlgorithmConfigArg.numPhiNeighbors,
            useExtraCuts=seedingAlgorithmConfigArg.useExtraCuts,
        ),
        gridConfig=gridConfig,
        gridOptions=gridOptions,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        inputPrimaryVertex=outputPrimaryVertex
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds

def addStandardSeedingNA60(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg,
    seedFinderConfigArg: SeedFinderConfigArgNA60,
    seedFinderOptionsArg: SeedFinderOptionsArgNA60,
    seedFilterConfigArg: SeedFilterConfigArgNA60,
    spacePointGridConfigArg: SpacePointGridConfigArg,
    logLevel: acts.logging.Level = None,
    outputPrimaryVertex: str="OutputRecPrimaryVertex",
    suffixSeed= ""
):
    """adds standard seeding
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    seedFinderConfig = acts.SeedFinderConfigNA60(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.y[0],
            rMax=seedFinderConfigArg.y[1],
            deltaYMin=seedFinderConfigArg.deltaY[0],
            deltaYMax=seedFinderConfigArg.deltaY[1],
            deltaYMinTopSP=(
                seedFinderConfigArg.deltaY[0]
                if seedFinderConfigArg.deltaYTopSP[0] is None
                else seedFinderConfigArg.deltaYTopSP[0]
            ),
            deltaYMaxTopSP=(
                seedFinderConfigArg.deltaY[1]
                if seedFinderConfigArg.deltaYTopSP[1] is None
                else seedFinderConfigArg.deltaYTopSP[1]
            ),
            deltaYMinBottomSP=(
                seedFinderConfigArg.deltaY[0]
                if seedFinderConfigArg.deltaYBottomSP[0] is None
                else seedFinderConfigArg.deltaYBottomSP[0]
            ),
            deltaYMaxBottomSP=(
                seedFinderConfigArg.deltaY[1]
                if seedFinderConfigArg.deltaYBottomSP[1] is None
                else seedFinderConfigArg.deltaYBottomSP[1]
            ),
            deltaYMiddleMinSPRange=seedFinderConfigArg.deltaYMiddleSPRange[0],
            deltaYMiddleMaxSPRange=seedFinderConfigArg.deltaYMiddleSPRange[1],
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            zOutermostLayers=(
                seedFinderConfigArg.zOutermostLayers[0]
                if seedFinderConfigArg.zOutermostLayers[0] is not None
                else seedFinderConfigArg.z[0],
                seedFinderConfigArg.zOutermostLayers[1]
                if seedFinderConfigArg.zOutermostLayers[1] is not None
                else seedFinderConfigArg.z[1],
            ),
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            zBinEdges=seedFinderConfigArg.zBinEdges,
            zBinsCustomLooping=seedFinderConfigArg.zBinsCustomLooping,
            skipZMiddleBinSearch=seedFinderConfigArg.skipZMiddleBinSearch,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            binSizeR=seedFinderConfigArg.binSizeR,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            seedConfirmationRange=seedFinderConfigArg.seedConfirmationRange,
            yMinMiddle = seedFinderConfigArg.yMiddle[0],
            yMaxMiddle = seedFinderConfigArg.yMiddle[1],
            verbose = seedFinderConfigArg.verbose
            

        ),
    )
    seedFinderOptions = acts.SeedFinderOptionsNA60(
        **acts.examples.defaultKWArgs(
            beamPos=acts.Vector2(0.0, 0.0)
            if seedFinderOptionsArg.beamPos == (None, None)
            else acts.Vector2(
                seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    
    )
    seedFilterConfig = acts.SeedFilterConfigNA60(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaYMin=seedFinderConfig.deltaYMin,
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            seedConfirmationRange=seedFinderConfig.seedConfirmationRange,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
            verbose=seedFilterConfigArg.verbose
        )
    )

    gridConfig = acts.PlanarSpacePointGridConfig(
        **acts.examples.defaultKWArgs(
            zMax=seedFinderConfig.zMax,
            zMin=seedFinderConfig.zMin,
            xMax=seedFinderConfig.zMax,
            xMin=seedFinderConfig.zMin,
            zBinEdges=spacePointGridConfigArg.zBinEdges,
            xBinEdges=spacePointGridConfigArg.zBinEdges,
        )
    )
    gridOptions = acts.PlanarSpacePointGridOptions(
        **acts.examples.defaultKWArgs(
            bFieldInZ=seedFinderOptions.bFieldInZ,
        )
    )

    seedingAlg = acts.examples.SeedingAlgorithmNA60(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds=suffixSeed+"seeds",
        **acts.examples.defaultKWArgs(
            allowSeparateRMax=seedingAlgorithmConfigArg.allowSeparateRMax,
            zBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            zBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom,
            xBinNeighborsTop=seedingAlgorithmConfigArg.zBinNeighborsTop,
            xBinNeighborsBottom=seedingAlgorithmConfigArg.zBinNeighborsBottom
        ),
        gridConfig=gridConfig,
        gridOptions=gridOptions,
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        inputPrimaryVertex=outputPrimaryVertex
    )
    sequence.addAlgorithm(seedingAlg)
    return seedingAlg.config.outputSeeds


def addOrthogonalSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    logLevel: acts.logging.Level = None,
):
    """adds orthogonal seeding algorithm
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    seedFinderConfig = acts.SeedFinderOrthogonalConfig(
        **acts.examples.defaultKWArgs(
            rMin=seedFinderConfigArg.r[0],
            rMax=seedFinderConfigArg.r[1],
            deltaRMinTopSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRTopSP[0] is None
                else seedFinderConfigArg.deltaRTopSP[0]
            ),
            deltaRMaxTopSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRTopSP[1] is None
                else seedFinderConfigArg.deltaRTopSP[1]
            ),
            deltaRMinBottomSP=(
                seedFinderConfigArg.deltaR[0]
                if seedFinderConfigArg.deltaRBottomSP[0] is None
                else seedFinderConfigArg.deltaRBottomSP[0]
            ),
            deltaRMaxBottomSP=(
                seedFinderConfigArg.deltaR[1]
                if seedFinderConfigArg.deltaRBottomSP[1] is None
                else seedFinderConfigArg.deltaRBottomSP[1]
            ),
            collisionRegionMin=seedFinderConfigArg.collisionRegion[0],
            collisionRegionMax=seedFinderConfigArg.collisionRegion[1],
            zMin=seedFinderConfigArg.z[0],
            zMax=seedFinderConfigArg.z[1],
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            cotThetaMax=seedFinderConfigArg.cotThetaMax,
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            radLengthPerSeed=seedFinderConfigArg.radLengthPerSeed,
            minPt=seedFinderConfigArg.minPt,
            impactMax=seedFinderConfigArg.impactMax,
            deltaPhiMax=seedFinderConfigArg.deltaPhiMax,
            interactionPointCut=seedFinderConfigArg.interactionPointCut,
            deltaZMax=seedFinderConfigArg.deltaZMax,
            maxPtScattering=seedFinderConfigArg.maxPtScattering,
            rRangeMiddleSP=seedFinderConfigArg.rRangeMiddleSP,
            useVariableMiddleSPRange=seedFinderConfigArg.useVariableMiddleSPRange,
            seedConfirmation=seedFinderConfigArg.seedConfirmation,
            centralSeedConfirmationRange=seedFinderConfigArg.centralSeedConfirmationRange,
            forwardSeedConfirmationRange=seedFinderConfigArg.forwardSeedConfirmationRange,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfigArg.deltaR[0]
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )
    seedingAlg = acts.examples.SeedingOrthogonalAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
    )
    sequence.addAlgorithm(seedingAlg)

    return seedingAlg.config.outputSeeds


def addHoughTransformSeeding(
    sequence: acts.examples.Sequencer,
    config: acts.examples.HoughTransformSeeder.Config,
    logLevel: acts.logging.Level = None,
):
    """
    Configures HoughTransform (HT) for seeding, instead of extra proxy config objects it takes
    directly the HT example algorithm config.
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    ht = acts.examples.HoughTransformSeeder(config=config, level=logLevel)
    sequence.addAlgorithm(ht)
    # potentially HT can be extended to also produce seeds, but it is not implemented yet
    # configuration option (outputSeeds) exists
    return ht.config.outputSeeds


def addGbtsSeeding(
    sequence: acts.examples.Sequencer,
    spacePoints: str,
    seedFinderConfigArg: SeedFinderConfigArg,
    seedFinderOptionsArg: SeedFinderOptionsArg,
    seedFilterConfigArg: SeedFilterConfigArg,
    trackingGeometry: acts.TrackingGeometry,
    logLevel: acts.logging.Level = None,
    layerMappingConfigFile: Union[Path, str] = None,
    geoSelectionConfigFile: Union[Path, str] = None,
    connector_inputConfigFile: Union[Path, str] = None,
):
    """Gbts seeding"""

    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    layerMappingFile = str(layerMappingConfigFile)  # turn path into string
    connector_inputFile = str(connector_inputConfigFile)
    seedFinderConfig = acts.SeedFinderGbtsConfig(
        **acts.examples.defaultKWArgs(
            sigmaScattering=seedFinderConfigArg.sigmaScattering,
            maxSeedsPerSpM=seedFinderConfigArg.maxSeedsPerSpM,
            minPt=seedFinderConfigArg.minPt,
            connector_input_file=connector_inputFile,
            m_useClusterWidth=False,
        ),
    )
    seedFinderOptions = acts.SeedFinderOptions(
        **acts.examples.defaultKWArgs(
            beamPos=(
                acts.Vector2(0.0, 0.0)
                if seedFinderOptionsArg.beamPos == (None, None)
                else acts.Vector2(
                    seedFinderOptionsArg.beamPos[0], seedFinderOptionsArg.beamPos[1]
                )
            ),
            bFieldInZ=seedFinderOptionsArg.bFieldInZ,
        )
    )
    seedFilterConfig = acts.SeedFilterConfig(
        **acts.examples.defaultKWArgs(
            maxSeedsPerSpM=seedFinderConfig.maxSeedsPerSpM,
            deltaRMin=(
                seedFinderConfigArg.deltaR[0]
                if seedFilterConfigArg.deltaRMin is None
                else seedFilterConfigArg.deltaRMin
            ),
            impactWeightFactor=seedFilterConfigArg.impactWeightFactor,
            zOriginWeightFactor=seedFilterConfigArg.zOriginWeightFactor,
            compatSeedWeight=seedFilterConfigArg.compatSeedWeight,
            compatSeedLimit=seedFilterConfigArg.compatSeedLimit,
            numSeedIncrement=seedFilterConfigArg.numSeedIncrement,
            seedWeightIncrement=seedFilterConfigArg.seedWeightIncrement,
            seedConfirmation=seedFilterConfigArg.seedConfirmation,
            # curvatureSortingInFilter=seedFilterConfigArg.curvatureSortingInFilter,
            maxSeedsPerSpMConf=seedFilterConfigArg.maxSeedsPerSpMConf,
            maxQualitySeedsPerSpMConf=seedFilterConfigArg.maxQualitySeedsPerSpMConf,
            useDeltaRorTopRadius=seedFilterConfigArg.useDeltaRorTopRadius,
        )
    )

    seedingAlg = acts.examples.GbtsSeedingAlgorithm(
        level=logLevel,
        inputSpacePoints=[spacePoints],
        outputSeeds="seeds",
        seedFilterConfig=seedFilterConfig,
        seedFinderConfig=seedFinderConfig,
        seedFinderOptions=seedFinderOptions,
        layerMappingFile=layerMappingFile,
        geometrySelection=acts.examples.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
        inputSourceLinks="sourcelinks",
        trackingGeometry=trackingGeometry,
        fill_module_csv=False,
        inputClusters="clusters",
    )

    sequence.addAlgorithm(seedingAlg)
    return seedingAlg.config.outputSeeds


def addSeedPerformanceWriters(
    sequence: acts.examples.Sequencer,
    outputDirRoot: Union[Path, str],
    seeds: str,
    prototracks: str,
    selectedParticles: str,
    inputParticles: str,
    outputTrackParameters: str,
    logLevel: acts.logging.Level = None,
    suffix ="",
    det_suffix ="",
    noGuessing = False
):
    """Writes seeding related performance output"""
    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir()

    sequence.addWriter(
        acts.examples.SeedingPerformanceWriter(
            level=customLogLevel(minLevel=acts.logging.DEBUG),
            inputSeeds=seeds,
            inputParticles=selectedParticles,
            inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
            filePath=str(outputDirRoot / "performance_seeding")+suffix+".root",
            verbose=False
        )
    )

    sequence.addWriter(
        acts.examples.RootTrackParameterWriter(
            level=customLogLevel(),
            inputTrackParameters=outputTrackParameters,
            inputProtoTracks=prototracks,
            inputParticles=inputParticles,
            inputSimHits="simhits"+det_suffix,
            inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
            inputMeasurementSimHitsMap="measurement_simhits_map"+det_suffix,
            filePath=str(outputDirRoot / "estimatedparams")+suffix+".root",
            treeName="estimatedparams",
        )
    )
    if not noGuessing:
        sequence.addWriter(
            acts.examples.TrackletVertexingPerformanceWriter(
                level=logLevel,
                inputSeeds=seeds,
                inputRecPrimaryVertex="OutputFitPrimaryVertex",
                inputGenPrimaryVertex="OutputGenPrimaryVertex",
                filePath = str(outputDirRoot / "performance_tracklet_vertexing.root"),
                fileMode = "RECREATE",
                verbose=False,
                inputFitFunction="OutputFitFuncVtx",
                inputZTracklets="OutputZTracklets",
                inputZTrackletsPeak="OutputZTrackletsPeak"
            )
        )
acts.examples.NamedTypeArgs(
    config=SeedFilterMLDBScanConfig,
)


def addSeedFilterML(
    s,
    config: SeedFilterMLDBScanConfig = SeedFilterMLDBScanConfig(),
    onnxModelFile: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)()
    from acts.examples.onnx import SeedFilterMLAlgorithm

    inputParticles = "particles"
    selectedParticles = "truth_seeds_selected"
    seeds = "seeds"
    estParams = "estimatedparameters"

    filterML = SeedFilterMLAlgorithm(
        level=customLogLevel,
        inputTrackParameters="estimatedparameters",
        inputSimSeeds="seeds",
        inputSeedFilterNN=onnxModelFile,
        outputTrackParameters="filtered-parameters",
        outputSimSeeds="filtered-seeds",
        **acts.examples.defaultKWArgs(
            epsilonDBScan=config.epsilonDBScan,
            minPointsDBScan=config.minPointsDBScan,
            minSeedScore=config.minSeedScore,
        ),
    )
    s.addAlgorithm(filterML)
    s.addWhiteboardAlias(seeds, "filtered-seeds")
    s.addWhiteboardAlias("estimatedparameters", "filtered-parameters")

    prototracks = "seed-prototracks-ML"
    s.addAlgorithm(
        acts.examples.SeedsToPrototracks(
            level=customLogLevel,
            inputSeeds=seeds,
            outputProtoTracks=prototracks,
        )
    )

    if outputDirRoot is not None:
        addSeedPerformanceWriters(
            s,
            outputDirRoot,
            seeds,
            prototracks,
            selectedParticles,
            inputParticles,
            estParams,
            customLogLevel
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)

        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        csvSeedWriter = acts.examples.CsvSeedWriter(
            level=customLogLevel,
            inputTrackParameters=estParams,
            inputSimSeeds=seeds,
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            outputDir=str(outputDirCsv),
            fileName=str(f"seed.csv"),
        )
        s.addWriter(csvSeedWriter)

    return s


def addKalmanTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    directNavigation: bool = False,
    reverseFilteringMomThreshold: float = 0 * u.GeV,
    inputProtoTracks: str = "truth_particle_tracks",
    multipleScattering: bool = True,
    energyLoss: bool = True,
    clusters: str = None,
    calibrator: acts.examples.MeasurementCalibrator = acts.examples.makePassThroughCalibrator(),
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if directNavigation:
        srfSortAlg = acts.examples.SurfaceSortingAlgorithm(
            level=customLogLevel(),
            inputProtoTracks=inputProtoTracks,
            inputSimHits="simhits",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            outputProtoTracks="sorted_truth_particle_tracks",
        )
        s.addAlgorithm(srfSortAlg)
        inputProtoTracks = srfSortAlg.config.outputProtoTracks

    kalmanOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "reverseFilteringMomThreshold": reverseFilteringMomThreshold,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "level": customLogLevel(),
    }

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        inputClusters=clusters if clusters is not None else "",
        outputTracks="kf_tracks",
        pickTrack=-1,
        fit=acts.examples.makeKalmanFitterFunction(
            trackingGeometry, field, **kalmanOptions
        ),
        calibrator=calibrator,
    )
    s.addAlgorithm(fitAlg)
    s.addWhiteboardAlias("tracks", fitAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=fitAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="kf_track_particle",
        outputParticleTrackMatching="kf_particle_track",
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track", matchAlg.config.outputParticleTrackMatching
    )

    return s


def addTruthTrackingGsf(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputProtoTracks: str = "truth_particle_tracks",
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    gsfOptions = {
        "betheHeitlerApprox": acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
        "maxComponents": 12,
        "componentMergeMethod": acts.examples.ComponentMergeMethod.maxWeight,
        "mixtureReductionAlgorithm": acts.examples.MixtureReductionAlgorithm.KLDistance,
        "weightCutoff": 1.0e-4,
        "level": customLogLevel(),
    }

    gsfAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        outputTracks="gsf_tracks",
        pickTrack=-1,
        fit=acts.examples.makeGsfFitterFunction(trackingGeometry, field, **gsfOptions),
        calibrator=acts.examples.makePassThroughCalibrator(),
    )
    s.addAlgorithm(gsfAlg)
    s.addWhiteboardAlias("tracks", gsfAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=gsfAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="gsf_track_particle",
        outputParticleTrackMatching="gsf_particle_track",
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track", matchAlg.config.outputParticleTrackMatching
    )

    return s

@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
    ckfConfig=CkfConfig,
)
def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    trackSelectorConfig: Optional[
        Union[TrackSelectorConfig, List[TrackSelectorConfig]]
    ] = None,
    ckfConfig: CkfConfig = CkfConfig(),
    twoWay: bool = False,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrajectories: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat=False,
    inputSourceLinks="sourcelinks",
    inputMeasurements="measurements",
    outputPrimaryVertex: str="OutputRecPrimaryVertex",
    suffixOut = "",
    suffixIn = "",
    det_suffix = ""
) -> None:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    trackSelectorConfig : TrackSelectorConfig(loc0, loc1, time, eta, absEta, pt, phi, minMeasurements)
        TrackSelector configuration. Each range is specified as a tuple of (min,max).
        Specify as a list(TrackSelectorConfig) for eta-dependent cuts, with binning specified by absEta[1].
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackSelector.hpp
    writeTrajectories : bool, True
        write trackstates_ckf.root and tracksummary_ckf.root ntuples? These can be quite large.
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    tslist = (
        []
        if trackSelectorConfig is None
        else (
            [trackSelectorConfig]
            if type(trackSelectorConfig) is TrackSelectorConfig
            else trackSelectorConfig
        )
    )
    cutSets = [
        acts.TrackSelector.Config(
            **acts.examples.defaultKWArgs(
                loc0Min=c.loc0[0],
                loc0Max=c.loc0[1],
                loc1Min=c.loc1[0],
                loc1Max=c.loc1[1],
                timeMin=c.time[0],
                timeMax=c.time[1],
                phiMin=c.phi[0],
                phiMax=c.phi[1],
                etaMin=c.eta[0],
                etaMax=c.eta[1],
                absEtaMin=c.absEta[0],
                absEtaMax=c.absEta[1] if len(tslist) == 1 else None,
                ptMin=c.pt[0],
                ptMax=c.pt[1],
                minMeasurements=c.nMeasurementsMin,
                maxHoles=c.maxHoles,
                maxOutliers=c.maxOutliers,
                maxSharedHits=c.maxSharedHits,
                maxChi2=c.maxChi2,
                measurementCounter=c.nMeasurementsGroupMin,
            )
        )
        for c in tslist
    ]
    if len(tslist) == 0:
        trkSelCfg = None
    elif len(tslist) == 1:
        trkSelCfg = cutSets[0]
    else:
        trkSelCfg = acts.TrackSelector.EtaBinnedConfig(
            cutSets=cutSets,
            absEtaEdges=[cutSets[0].absEtaMin] + [c.absEta[1] for c in tslist],
        )

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [
                (
                    acts.GeometryIdentifier(),
                    (
                        [],
                        [ckfConfig.chi2CutOff],
                        [ckfConfig.numMeasurementsCutOff],
                    ),
                )
            ]
        ),
        inputMeasurements=inputMeasurements,
        inputSourceLinks=inputSourceLinks,
        #prima era suffixIn
        inputInitialTrackParameters=suffixIn+"estimatedparameters",
        inputSeeds=(
            suffixIn+"estimatedseeds"
            if ckfConfig.seedDeduplication or ckfConfig.stayOnSeed
            else ""
        ),
        outputTracks=suffixOut+"ckfTracks",
        inputPrimaryVertex=outputPrimaryVertex,
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field, customLogLevel()
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            maxSteps=ckfConfig.maxSteps,
            twoWay=twoWay,
            seedDeduplication=ckfConfig.seedDeduplication,
            stayOnSeed=ckfConfig.stayOnSeed,
            pixelVolumes=ckfConfig.pixelVolumes,
            stripVolumes=ckfConfig.stripVolumes,
            maxPixelHoles=ckfConfig.maxPixelHoles,
            maxStripHoles=ckfConfig.maxStripHoles,
        ),
    )
    s.addAlgorithm(trackFinder)
    s.addWhiteboardAlias(suffixOut+"tracks", trackFinder.config.outputTracks)

    matcher = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=trackFinder.config.outputTracks,
        inputParticles="particles_selected"+det_suffix,
        inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
        outputTrackParticleMatching=suffixOut+"ckf_track_particle",
        outputParticleTrackMatching=suffixOut+"ckf_particle_track",
    )
    s.addAlgorithm(matcher)

    s.addWhiteboardAlias(
        suffixOut+"track_particle", matcher.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        suffixOut+"particle_track", matcher.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name=suffixOut+"ckf",
        tracks=trackFinder.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        logLevel=logLevel,
        writeCovMat=writeCovMat,
        suffix=suffixOut,
        det_suffix=det_suffix
    )

    return s

def addMatching(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    magneticField: acts.MagneticFieldProvider,
    inputTrackParametersMS = "ambitrackspars",
    inputTrackParametersVT = "ambitrackspars",
    inputTrackContainerMS = "ambitracks",
    inputTrackContainerVT = "ambitracks",
    outputTrackParameters = "outputTrackParameters",
    outputTracks = "outputTracks",
    outputMatchedTracks = "outputMatchedTracks",
    inputParticles = "particles",
    inputMeasurementParticlesMapVT = "measurement_particles_map_vt",
    inputMeasurementParticlesMapMS = "measurement_particles_map_ms",
    useRecVtx = False,
    px = 0,
    py = 0,
    pz = 0,
    chi2max = 10000000,
    logLevel: Optional[acts.logging.Level] = None
) -> None:

    customLogLevel = acts.examples.defaultLogging(s, logLevel)
    """
    converter = acts.examples.TracksToParameters(
        level=customLogLevel(),
        inputTracks=inputTrackContainerMS,
        outputTrackParameters="ambitrackspars",
    )
    s.addAlgorithm(converter)
    """
    gsfOptions = {
        "betheHeitlerApprox": acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
        "maxComponents": 12,
        "componentMergeMethod": acts.examples.ComponentMergeMethod.maxWeight,
        "mixtureReductionAlgorithm": acts.examples.MixtureReductionAlgorithm.KLDistance,
        "weightCutoff": 1.0e-4,
        "level": acts.logging.INFO,
    }



    kalmanOptions = {
        "multipleScattering": True,
        "energyLoss": True,
        "reverseFilteringMomThreshold": 0 * u.GeV,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "level": customLogLevel(),
    }

    #trackParameters = converter.config.outputTrackParameters
    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackMatcher = acts.examples.MatchingAlgorithm(
        level=customLogLevel(),
        inputTrackParametersMS=inputTrackParametersMS,
        inputTrackParametersVT=inputTrackParametersVT,
        inputTrackContainerMS=inputTrackContainerMS,
        inputTrackContainerVT=inputTrackContainerVT,
        outputTrackParameters=outputTrackParameters,
        outputTracks=outputTracks,
        outputMatchedTracks = outputMatchedTracks,
        inputParticles = inputParticles,
        inputMeasurementParticlesMapVT = inputMeasurementParticlesMapVT,
        inputMeasurementParticlesMapMS = inputMeasurementParticlesMapMS,
        useRecVtx = useRecVtx,
        px = px,
        py = py,
        pz = pz,
        chi2max = chi2max,
        trackingGeometry = trackingGeometry,
        magneticField = magneticField,
        
        #fit=acts.examples.makeKalmanFitterFunction(
        #    trackingGeometry, magneticField, **kalmanOptions
        #)
        fit=acts.examples.makeGsfFitterFunction(
                trackingGeometry, magneticField, **gsfOptions
            )
    )
    s.addAlgorithm(trackMatcher)
    #s.addWhiteboardAlias("match", trackFinder.config.outputTracks)

    """
    matcher = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=trackFinder.config.outputTracks,
        inputParticles="particles_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching=suffixOut+"ckf_track_particle",
        outputParticleTrackMatching=suffixOut+"ckf_particle_track",
    )
    s.addAlgorithm(matcher)
    s.addWhiteboardAlias(
        suffixOut+"track_particle", matcher.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        suffixOut+"particle_track", matcher.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name=suffixOut+"ckf",
        tracks=trackFinder.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        logLevel=logLevel,
        writeCovMat=writeCovMat,
        suffix=suffixOut
    )
    """
    return s


def addGx2fTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    inputProtoTracks: str = "truth_particle_tracks",
    multipleScattering: bool = False,
    energyLoss: bool = False,
    nUpdateMax: int = 5,
    relChi2changeCutOff: float = 1e-7,
    clusters: str = None,
    calibrator: acts.examples.MeasurementCalibrator = acts.examples.makePassThroughCalibrator(),
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    gx2fOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "nUpdateMax": nUpdateMax,
        "relChi2changeCutOff": relChi2changeCutOff,
        "level": customLogLevel(),
    }

    fitAlg = acts.examples.TrackFittingAlgorithm(
        level=customLogLevel(),
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputProtoTracks=inputProtoTracks,
        inputInitialTrackParameters="estimatedparameters",
        inputClusters=clusters if clusters is not None else "",
        outputTracks="gx2f_tracks",
        pickTrack=-1,
        fit=acts.examples.makeGlobalChiSquareFitterFunction(
            trackingGeometry, field, **gx2fOptions
        ),
        calibrator=calibrator,
    )
    s.addAlgorithm(fitAlg)
    s.addWhiteboardAlias("tracks", fitAlg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=fitAlg.config.outputTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="gx2f_track_particle",
        outputParticleTrackMatching="gx2f_particle_track",
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track", matchAlg.config.outputParticleTrackMatching
    )

    return s


def addTrackWriters(
    s: acts.examples.Sequencer,
    name: str,
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeStates: bool = True,
    writeSummary: bool = True,
    writeCKFperformance: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat=False,
    suffix = "",
    det_suffix = ""
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        if writeStates:
            # write track states from CKF
            trackStatesWriter = acts.examples.RootTrackStatesWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                # @note The full particles collection is used here to avoid lots of warnings
                # since the unselected CKF track might have a majority particle not in the
                # filtered particle collection. This could be avoided when a separate track
                # selection algorithm is used.
                inputParticles="particles_selected"+det_suffix,
                inputTrackParticleMatching=suffix+"track_particle",
                inputSimHits="simhits"+det_suffix,
                inputMeasurementSimHitsMap="measurement_simhits_map"+det_suffix,
                filePath=str(outputDirRoot / f"trackstates_{name}.root"),
                treeName="trackstates",
            )
            s.addWriter(trackStatesWriter)

        if writeSummary:
            # write track summary from CKF
            trackSummaryWriter = acts.examples.RootTrackSummaryWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                # @note The full particles collection is used here to avoid lots of warnings
                # since the unselected CKF track might have a majority particle not in the
                # filtered particle collection. This could be avoided when a separate track
                # selection algorithm is used.
                inputParticles="particles_selected"+det_suffix,
                inputTrackParticleMatching=suffix+"track_particle",
                filePath=str(outputDirRoot / f"tracksummary_{name}.root"),
                treeName="tracksummary",
                writeCovMat=writeCovMat,
            )
            s.addWriter(trackSummaryWriter)

        if writeCKFperformance:
            # Write CKF performance data
            ckfPerfWriter = acts.examples.CKFPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="truth_seeds_selected"+det_suffix,
                inputTrackParticleMatching= suffix+"track_particle",
                inputParticleTrackMatching= suffix+"particle_track",
                filePath=str(outputDirRoot / f"performance_{name}.root"),
            )
            s.addWriter(ckfPerfWriter)

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        if writeSummary:
            csvWriter = acts.examples.CsvTrackWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputMeasurementParticlesMap="measurement_particles_map",
                outputDir=str(outputDirCsv),
                fileName=str(f"tracks_{name}.csv"),
            )
            s.addWriter(csvWriter)


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
)
def addTrackSelection(
    s: acts.examples.Sequencer,
    trackSelectorConfig: TrackSelectorConfig,
    inputTracks: str,
    outputTracks: str,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.TrackSelectorAlgorithm:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # single cut config for implicit single bin eta configuration
    selectorConfig = acts.TrackSelector.Config(
        **acts.examples.defaultKWArgs(
            loc0Min=trackSelectorConfig.loc0[0],
            loc0Max=trackSelectorConfig.loc0[1],
            loc1Min=trackSelectorConfig.loc1[0],
            loc1Max=trackSelectorConfig.loc1[1],
            timeMin=trackSelectorConfig.time[0],
            timeMax=trackSelectorConfig.time[1],
            phiMin=trackSelectorConfig.phi[0],
            phiMax=trackSelectorConfig.phi[1],
            etaMin=trackSelectorConfig.eta[0],
            etaMax=trackSelectorConfig.eta[1],
            absEtaMin=trackSelectorConfig.absEta[0],
            absEtaMax=trackSelectorConfig.absEta[1],
            ptMin=trackSelectorConfig.pt[0],
            ptMax=trackSelectorConfig.pt[1],
            minMeasurements=trackSelectorConfig.nMeasurementsMin,
        )
    )

    trackSelector = acts.examples.TrackSelectorAlgorithm(
        level=customLogLevel(),
        inputTracks=inputTracks,
        outputTracks=outputTracks,
        selectorConfig=selectorConfig,
    )

    s.addAlgorithm(trackSelector)

    return trackSelector


ExaTrkXBackend = Enum("ExaTrkXBackend", "Torch Onnx")


def addExaTrkX(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    geometrySelection: Union[Path, str],
    modelDir: Union[Path, str],
    outputDirRoot: Optional[Union[Path, str]] = None,
    backend: Optional[ExaTrkXBackend] = ExaTrkXBackend.Torch,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Run the particle selection
    # The pre-selection will select truth particles satisfying provided criteria
    # from all particles read in by particle reader for further processing. It
    # has no impact on the truth hits themselves
    s.addAlgorithm(
        acts.examples.TruthSeedSelector(
            level=customLogLevel(),
            ptMin=500 * u.MeV,
            nHitsMin=9,
            inputParticles="particles_initial",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputParticles="particles_seed_selected",
        )
    )

    # Create space points
    s.addAlgorithm(
        acts.examples.SpacePointMaker(
            level=customLogLevel(),
            inputSourceLinks="sourcelinks",
            inputMeasurements="measurements",
            outputSpacePoints="spacepoints",
            trackingGeometry=trackingGeometry,
            geometrySelection=acts.examples.readJsonGeometryList(
                str(geometrySelection)
            ),
        )
    )

    metricLearningConfig = {
        "level": customLogLevel(),
        "embeddingDim": 8,
        "rVal": 1.6,
        "knnVal": 100,
    }

    filterConfig = {
        "level": customLogLevel(),
        "cut": 0.01,
    }

    gnnConfig = {
        "level": customLogLevel(),
        "cut": 0.5,
    }

    if backend == ExaTrkXBackend.Torch:
        metricLearningConfig["modelPath"] = str(modelDir / "embed.pt")
        metricLearningConfig["numFeatures"] = 3
        filterConfig["modelPath"] = str(modelDir / "filter.pt")
        filterConfig["nChunks"] = 10
        filterConfig["numFeatures"] = 3
        gnnConfig["modelPath"] = str(modelDir / "gnn.pt")
        gnnConfig["undirected"] = True
        gnnConfig["numFeatures"] = 3

        graphConstructor = acts.examples.TorchMetricLearning(**metricLearningConfig)
        edgeClassifiers = [
            acts.examples.TorchEdgeClassifier(**filterConfig),
            acts.examples.TorchEdgeClassifier(**gnnConfig),
        ]
        trackBuilder = acts.examples.BoostTrackBuilding(customLogLevel())
    elif backend == ExaTrkXBackend.Onnx:
        metricLearningConfig["modelPath"] = str(modelDir / "embedding.onnx")
        metricLearningConfig["spacepointFeatures"] = 3
        filterConfig["modelPath"] = str(modelDir / "filtering.onnx")
        gnnConfig["modelPath"] = str(modelDir / "gnn.onnx")

        graphConstructor = acts.examples.OnnxMetricLearning(**metricLearningConfig)
        edgeClassifiers = [
            acts.examples.OnnxEdgeClassifier(**filterConfig),
            acts.examples.OnnxEdgeClassifier(**gnnConfig),
        ]
        trackBuilder = acts.examples.CugraphTrackBuilding(customLogLevel())

    findingAlg = acts.examples.TrackFindingAlgorithmExaTrkX(
        level=customLogLevel(),
        inputSpacePoints="spacepoints",
        outputProtoTracks="exatrkx_prototracks",
        graphConstructor=graphConstructor,
        edgeClassifiers=edgeClassifiers,
        trackBuilder=trackBuilder,
    )
    s.addAlgorithm(findingAlg)
    s.addWhiteboardAlias("prototracks", findingAlg.config.outputProtoTracks)

    matchAlg = acts.examples.ProtoTrackTruthMatcher(
        level=customLogLevel(),
        inputProtoTracks=findingAlg.config.outputProtoTracks,
        inputParticles="particles",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputProtoTrackParticleMatching="exatrkx_prototrack_particle",
        outputParticleProtoTrackMatching="exatrkx_particle_prototrack",
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "prototrack_particle", matchAlg.config.outputProtoTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_prototrack", matchAlg.config.outputParticleProtoTrackMatching
    )

    # Write truth track finding / seeding performance
    if outputDirRoot is not None:
        s.addWriter(
            acts.examples.TrackFinderPerformanceWriter(
                level=customLogLevel(),
                inputProtoTracks=findingAlg.config.outputProtoTracks,
                # the original selected particles after digitization
                inputParticles="particles_initial",
                inputMeasurementParticlesMap="measurement_particles_map",
                inputProtoTrackParticleMatching=matchAlg.config.outputProtoTrackParticleMatching,
                filePath=str(Path(outputDirRoot) / "performance_track_finding.root"),
            )
        )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionConfig,
)
def addAmbiguityResolution(
    s,
    config: AmbiguityResolutionConfig = AmbiguityResolutionConfig(),
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrajectories: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat=False,
    suffixIn = "",
    suffixOut = "ambi",
    det_suffix = ""
) -> None:
    from acts.examples import GreedyAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    alg = GreedyAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=suffixIn+"tracks",
        outputTracks=suffixOut+"tracks",
        **acts.examples.defaultKWArgs(
            maximumSharedHits=config.maximumSharedHits,
            nMeasurementsMin=config.nMeasurementsMin,
            maximumIterations=config.maximumIterations,
        ),
    )
    s.addAlgorithm(alg)
    s.addWhiteboardAlias(suffixOut+"tracks", alg.config.outputTracks)

    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=alg.config.outputTracks,
        inputParticles="particles"+det_suffix,
        inputMeasurementParticlesMap="measurement_particles_map"+det_suffix,
        outputTrackParticleMatching=suffixOut+"ambi_track_particle",
        outputParticleTrackMatching=suffixOut+"ambi_particle_track",
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        suffixOut+"track_particle", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        suffixOut+"particle_track", matchAlg.config.outputParticleTrackMatching
    )

    addTrackWriters(
        s,
        name=suffixOut,
        tracks=suffixOut+"tracks",
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        logLevel=logLevel,
        writeCovMat=writeCovMat,
        suffix =  suffixOut,
        det_suffix = det_suffix
    )
    
    converter = acts.examples.TracksToParameters(
        level=customLogLevel(),
        inputTracks=suffixOut+"tracks",
        outputTrackParameters=suffixOut+"trackpars",
    )
    s.addAlgorithm(converter)

    return s


@acts.examples.NamedTypeArgs(
    config=ScoreBasedAmbiguityResolutionConfig,
)
def addScoreBasedAmbiguityResolution(
    s,
    config: ScoreBasedAmbiguityResolutionConfig = ScoreBasedAmbiguityResolutionConfig(),
    tracks: str = "tracks",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    ambiVolumeFile: Optional[Union[Path, str]] = None,
    writeTrajectories: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat=False,
) -> None:
    from acts.examples import ScoreBasedAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, acts.logging.INFO)

    algScoreBased = ScoreBasedAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=tracks,
        configFile=ambiVolumeFile,
        outputTracks="ambiTracksScoreBased",
        **acts.examples.defaultKWArgs(
            minScore=config.minScore,
            minScoreSharedTracks=config.minScoreSharedTracks,
            maxShared=config.maxShared,
            maxSharedTracksPerMeasurement=config.maxSharedTracksPerMeasurement,
            phiMax=config.phiMax,
            phiMin=config.phiMin,
            etaMax=config.etaMax,
            etaMin=config.etaMin,
            useAmbiguityFunction=config.useAmbiguityFunction,
        ),
    )
    s.addAlgorithm(algScoreBased)
    s.addWhiteboardAlias("tracks", algScoreBased.config.outputTracks)

    addTrackWriters(
        s,
        name="ambi_scorebased",
        tracks=algScoreBased.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        logLevel=logLevel,
        writeCovMat=writeCovMat,
    )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionMLConfig,
)
def addAmbiguityResolutionML(
    s,
    config: AmbiguityResolutionMLConfig = AmbiguityResolutionMLConfig(),
    onnxModelFile: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrajectories: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples.onnx import AmbiguityResolutionMLAlgorithm
    from acts.examples import GreedyAmbiguityResolutionAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    algML = AmbiguityResolutionMLAlgorithm(
        level=customLogLevel(),
        inputTracks="tracks",
        inputDuplicateNN=onnxModelFile,
        outputTracks="ambiTracksML",
        **acts.examples.defaultKWArgs(
            nMeasurementsMin=config.nMeasurementsMin,
        ),
    )

    algGreedy = GreedyAmbiguityResolutionAlgorithm(
        level=customLogLevel(),
        inputTracks=algML.config.outputTracks,
        outputTracks="ambiTracksMLGreedy",
        **acts.examples.defaultKWArgs(
            maximumSharedHits=config.maximumSharedHits,
            nMeasurementsMin=config.nMeasurementsMin,
            maximumIterations=config.maximumIterations,
        ),
    )

    s.addAlgorithm(algML)
    s.addAlgorithm(algGreedy)

    addTrackWriters(
        s,
        name="ambiML",
        tracks=algGreedy.config.outputTracks,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    config=AmbiguityResolutionMLDBScanConfig,
)
def addAmbiguityResolutionMLDBScan(
    s,
    config: AmbiguityResolutionMLDBScanConfig = AmbiguityResolutionMLDBScanConfig(),
    onnxModelFile: Optional[Union[Path, str]] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrajectories: bool = True,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    from acts.examples import AmbiguityResolutionMLDBScanAlgorithm

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    alg = AmbiguityResolutionMLDBScanAlgorithm(
        level=customLogLevel(),
        inputTracks="tracks",
        inputDuplicateNN=onnxModelFile,
        outputTracks="ambiTracksMLDBScan",
        **acts.examples.defaultKWArgs(
            nMeasurementsMin=config.nMeasurementsMin,
            epsilonDBScan=config.epsilonDBScan,
            minPointsDBScan=config.minPointsDBScan,
        ),
    )
    s.addAlgorithm(alg)

    addTrackWriters(
        s,
        name="ambiMLDBScan",
        trajectories=alg.config.outputTracks,
        outputDirRoot=outputDirRoot,
        outputDirCsv=outputDirCsv,
        writeStates=writeTrajectories,
        writeSummary=writeTrajectories,
        writeCKFperformance=True,
        logLevel=logLevel,
    )

    return s


@acts.examples.NamedTypeArgs(
    trackSelectorConfig=TrackSelectorConfig,
)
def addVertexFitting(
    s,
    field,
    tracks: Optional[str] = "tracks",
    trackParameters: Optional[str] = None,
    outputProtoVertices: str = "protovertices",
    outputVertices: str = "fittedVertices",
    inputParticles: str = "particles_input",
    selectedParticles: str = "particles_selected",
    inputVertices: str = "vertices_input",
    seeder: Optional[acts.VertexSeedFinder] = acts.VertexSeedFinder.GaussianSeeder,
    vertexFinder: VertexFinder = VertexFinder.Truth,
    maxIterations: Optional[int] = None,
    initialVariances: Optional[List[float]] = None,
    useTime: Optional[bool] = False,
    seeder: Optional[acts.VertexSeedFinder] = acts.VertexSeedFinder.GaussianSeeder,
    spatialBinExtent: Optional[float] = None,
    temporalBinExtent: Optional[float] = None,
    trackSelectorConfig: Optional[TrackSelectorConfig] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    suffixIn: str = "",
    suffixOut: str = "",
    significanceCutSeeding: Optional[float] = 25.158658887529214,
    maximumChi2cutForSeeding: Optional[float] = 40.420277753035,
    maxVertices: Optional[int] = 1,
    createSplitVertices: Optional[bool] = False,
    splitVerticesTrkInvFraction: Optional[int] = 2,
    reassignTracksAfterFirstFit: Optional[bool] = False,
    doMaxTracksCut: Optional[bool] = False,
    maxTracks: Optional[int] = 5000,
    cutOffTrackWeight: Optional[float] = 0.6249340642230282,
    cutOffTrackWeightReassign: Optional[float] = 1,
    rejectedFraction: Optional[float] = 0.5,
) -> None:
    """This function steers the vertex fitting

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from
        addVertexFitting)
    field : magnetic field
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    vertexFinder : VertexFinder, Truth
        vertexFinder algorithm: one of Truth, AMVF, Iterative
    seeder : enum member
        determines vertex seeder for AMVF, can be acts.seeder.GaussianSeeder or
        acts.seeder.AdaptiveGridSeeder
    useTime : bool, False
        determines whether time information is used in vertex seeder, finder,
        and fitter
        only implemented for the AMVF and the AdaptiveGridSeeder
    spatialBinExtent : float, None
        spatial bin extent for the AdaptiveGridSeeder
    temporalBinExtent : float, None
        temporal bin extent for the AdaptiveGridSeeder
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    """
    from acts.examples import (
        TruthVertexFinder,
        VertexFitterAlgorithm,
        IterativeVertexFinderAlgorithm,
        AdaptiveMultiVertexFinderAlgorithm,
        VertexPerformanceWriter,
    )

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if tracks is not None and trackSelectorConfig is not None:
        trackSelector = addTrackSelection(
            s,
            trackSelectorConfig,
            inputTracks=suffixIn+tracks,
            outputTracks="selectedTracksVertexing",
            logLevel=customLogLevel(),
        )
        tracks = trackSelector.config.outputTracks

    if trackParameters is None:
        converter = acts.examples.TracksToParameters(
            level=customLogLevel(),
            inputTracks=suffixIn+tracks,
            outputTrackParameters="selectedTracksParametersVertexing",
        )
        s.addAlgorithm(converter)
        trackParameters = converter.config.outputTrackParameters

    tracks = tracks if tracks is not None else ""

    if vertexFinder == VertexFinder.Truth:
        findVertices = TruthVertexFinder(
            level=customLogLevel(),
            inputTracks=suffixIn+tracks,
            inputParticles=selectedParticles,
            inputMeasurementParticlesMap="measurement_particles_map",
            outputProtoVertices=outputProtoVertices,
            excludeSecondaries=True,
        )
        s.addAlgorithm(findVertices)
        fitVertices = VertexFitterAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            inputProtoVertices=findVertices.config.outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
        )
        s.addAlgorithm(fitVertices)
    elif vertexFinder == VertexFinder.Iterative:
        findVertices = IterativeVertexFinderAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            outputProtoVertices=outputProtoVertices,
            outputVertices=outputVertices,
            significanceCutSeeding=significanceCutSeeding,
            maximumChi2cutForSeeding=maximumChi2cutForSeeding,
            maxVertices=maxVertices,
            createSplitVertices=createSplitVertices,
            splitVerticesTrkInvFraction=splitVerticesTrkInvFraction,
            reassignTracksAfterFirstFit=reassignTracksAfterFirstFit,
            doMaxTracksCut=doMaxTracksCut,
            maxTracks=maxTracks,
            cutOffTrackWeight=cutOffTrackWeight,
            cutOffTrackWeightReassign=cutOffTrackWeightReassign,
            rejectedFraction=rejectedFraction
        )
        s.addAlgorithm(findVertices)
    elif vertexFinder == VertexFinder.AMVF:
        findVertices = AdaptiveMultiVertexFinderAlgorithm(
            level=customLogLevel(),
            inputTrackParameters=trackParameters,
            outputProtoVertices=outputProtoVertices,
            outputVertices=outputVertices,
            bField=field,
            seedFinder=seeder,
            **acts.examples.defaultKWArgs(
                maxIterations=maxIterations,
                initialVariances=initialVariances,
                useTime=useTime,
                spatialBinExtent=spatialBinExtent,
                temporalBinExtent=temporalBinExtent,
            ),
        )
        s.addAlgorithm(findVertices)
    else:
        raise RuntimeError("Invalid finder argument")

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        s.addWriter(
            VertexPerformanceWriter(
                level=customLogLevel(),
                inputVertices=outputVertices,
                inputTracks=suffixIn+tracks,
                inputTruthVertices=inputVertices,
                inputParticles=inputParticles,
                inputSelectedParticles=selectedParticles,
                inputTrackParticleMatching=suffixIn+"track_particle",
                bField=field,
                treeName="vertexing",
                filePath=str(outputDirRoot / "performance_vertexing")+suffixOut+".root",
            )
        )

    return s

def addSecondaryVertexFitting(
    s,
    field,
    tracks: Optional[str] = "tracks",
    trackParameters: Optional[str] = None,
    outputProtoVertices: str = "protovertices",
    outputMasses: str = "masses",
    trackSelectorConfig: Optional[TrackSelectorConfig] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    suffixIn: str = "",
    suffixOut: str = "",
) -> None:
    
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if tracks is not None and trackSelectorConfig is not None:
        trackSelector = addTrackSelection(
            s,
            trackSelectorConfig,
            inputTracks=suffixIn+tracks,
            outputTracks="selectedTracksVertexing",
            logLevel=customLogLevel(),
        )
        tracks = trackSelector.config.outputTracks

    if trackParameters is None:
        converter = acts.examples.TracksToParameters(
            level=customLogLevel(),
            inputTracks=suffixIn+tracks,
            outputTrackParameters="selectedTracksParametersVertexing",
        )
        s.addAlgorithm(converter)
        trackParameters = converter.config.outputTrackParameters


    findVertices =  acts.examples.SecondaryVertexFinderAlgorithm(
        level=customLogLevel(),
        bField=field,
        inputTrackParameters=trackParameters,
        outputProtoVertices=outputProtoVertices,
        outputMasses=outputMasses
    )
    s.addAlgorithm(findVertices)

    return s


def addSingleSeedVertexFinding(
    s,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    inputSpacePoints: Optional[str] = "spacepoints",
    outputVertices: Optional[str] = "fittedSeedVertices",
) -> None:
    from acts.examples import (
        SingleSeedVertexFinderAlgorithm,
        VertexPerformanceWriter,
    )

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    findSingleSeedVertex = SingleSeedVertexFinderAlgorithm(
        level=customLogLevel(),
        inputSpacepoints=inputSpacePoints,
        outputVertices=outputVertices,
    )
    s.addAlgorithm(findSingleSeedVertex)

    inputParticles = "particles_input"
    selectedParticles = "particles_selected"

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            VertexPerformanceWriter(
                level=customLogLevel(),
                inputAllTruthParticles=inputParticles,
                inputSelectedTruthParticles=selectedParticles,
                useTracks=False,
                inputVertices=outputVertices,
                treeName="seedvertexing",
                filePath=str(outputDirRoot / "performance_seedvertexing.root"),
            )
        )

    return s



def addUsedMeasurementsFilter(
    s: acts.examples.Sequencer,
    inputSourceLinks: str="sourcelinks",
    outputSourceLinks: str = "outputsourcelinks",
    inputTracks: str = "tracks",
    logLevel: Optional[acts.logging.Level] = None,
    ) -> None:

    logLevel = acts.examples.defaultLogging(s, logLevel)()

    selAlg = acts.examples.FilterMeasurementsAlgorithm(
        level=logLevel,
        inputSourceLinks=inputSourceLinks,
        inputTracks=inputTracks,
        outputSourceLinks=outputSourceLinks
    )

    s.addAlgorithm(selAlg)

    return s



def addTrackletVertexing(
    s: acts.examples.Sequencer,
    inputSpacePoints: str="spacepoints",
    inputParticles: str="particles",
    inputMeasurementParticlesMap: str="measurement_particles_map",
    outputRecPrimaryVertex: str="OutputRecPrimaryVertex",
    outputGenPrimaryVertex: str="OutputGenPrimaryVertex",
    zmax: float=170,
    zmin: float=0,
    deltaPhi: float=0.1,
    deltaThetaMax: float=0.04,
    deltaThetaMin: float=-0.15,
    verbose: bool=False,
    doMCtruth: bool=True,
    noGuessing: bool=False,
    useFit: bool=True,
    nbins: int=60,
    addDeltas: bool=False,
    projective: bool=False,
    zPerigee: float=0,
    logLevel: Optional[acts.logging.Level] = None,

    ) -> None:

    logLevel = acts.examples.defaultLogging(s, logLevel)()

    selAlg = acts.examples.TrackletVertexingAlgorithm(
        level=logLevel,
        inputSpacePoints=inputSpacePoints,
        inputSpacePointsMC=[inputSpacePoints],
        inputParticles=inputParticles,
        outputRecPrimaryVertex=outputRecPrimaryVertex,
        outputGenPrimaryVertex=outputGenPrimaryVertex,
        outputFitPrimaryVertex="OutputFitPrimaryVertex",
        inputMeasurementParticlesMap=inputMeasurementParticlesMap,
        zmax=zmax,
        zmin=zmin,
        deltaPhi=deltaPhi,
        deltaThetaMin=deltaThetaMin,
        deltaThetaMax=deltaThetaMax,
        verbose=False,#verbose,
        doMCtruth=doMCtruth,
        noGuessing=noGuessing,
        useFit=useFit,
        nbins=nbins,
        addDeltas=addDeltas,
        projective=projective,
        zPerigee=zPerigee,
        outputFitFunction="OutputFitFuncVtx",
        outputZTracklets="OutputZTracklets",
        outputZTrackletsPeak="OutputZTrackletsPeak"
    )

    s.addAlgorithm(selAlg)
    return s


def addContainerMerger(
    s: acts.examples.Sequencer,
    inputTrackParameters=["tracks"],
    outputTrackParameters="mergetracks",
    inputTracks=["tracks"],
    outputTracks="mergetracks",
    logLevel: Optional[acts.logging.Level] = None,

    ) -> None:

    logLevel = acts.examples.defaultLogging(s, logLevel)()

    selAlg = acts.examples.MergeContainersAlgorithm(
        level=logLevel,
        inputTrackParameters=inputTrackParameters,
        outputTrackParameters=outputTrackParameters,
        inputTracks=inputTracks,
        outputTracks=outputTracks
    )

    s.addAlgorithm(selAlg)
    return s
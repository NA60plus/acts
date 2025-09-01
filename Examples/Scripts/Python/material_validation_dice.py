#!/usr/bin/env python3

import os
import argparse

import acts
from acts.examples import Sequencer, RootMaterialTrackWriter
from acts.examples.simulation import addParticleGun, EtaConfig, ParticleConfig, MomentumConfig, PhiConfig
from dice_geometry import buildDICEgeometry
from acts import UnitConstants as u


def runMaterialValidation(
    nevents,
    ntracks,
    trackingGeometry,
    decorators,
    field,
    outputDir,
    outputName="propagation-material",
    s=None,
):
    # Create a sequencer
    s = s or Sequencer(events=nevents, numThreads=-1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    nav = acts.Navigator(trackingGeometry=trackingGeometry)
    stepper = acts.StraightLineStepper()
    # stepper = acts.EigenStepper(field)

    prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        rnd=rnd,
        momentumConfig = MomentumConfig(1 * u.GeV, 10.0 * u.GeV, transverse=False),
        etaConfig = EtaConfig(0, 8, uniform = True),
        phiConfig = PhiConfig(0.0, 360.0 * u.degree),            
        particleConfig = ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
        vtxGen = acts.examples.GaussianVertexGenerator(
                            mean=acts.Vector4(0, 0, 0, 0),
                            stddev=acts.Vector4(0, 0, 0, 0),
                        )
    )

    trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
        level=acts.logging.INFO,
        inputParticles="particles_generated",
        outputTrackParameters="params_particles_generated",
    )
    s.addAlgorithm(trkParamExtractor)

    alg = acts.examples.PropagationAlgorithm(
        propagatorImpl=prop,
        level=acts.logging.INFO,
        sterileLogger=True,
        recordMaterialInteractions=True,
        inputTrackParameters="params_particles_generated",
        outputSummaryCollection="propagation_summary",
        outputMaterialCollection="material_tracks",
    )
    s.addAlgorithm(alg)

    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=alg.config.outputMaterialCollection,
            filePath=os.path.join(outputDir, (outputName + ".root")),
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to process"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=1000, help="Number of tracks per event"
    )
    p.add_argument(
        "-m", "--map", type=str, help="Input file (optional) for the material map"
    )
    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="propagation-material",
        help="Output file name",
    )
    p.add_argument("--remove-vs", action="store_true", help="Remove Vertex Spectrometer")
    p.add_argument("--remove-ms", action="store_true", help="Remove Muon Spectrometer")
    
    args = p.parse_args()

    matDeco = (
        acts.IMaterialDecorator.fromFile(args.map) if args.map != None else None
    )

    detector = buildDICEgeometry(matDeco=matDeco, addVS=not args.remove_vs, addMS=not args.remove_ms)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    field = acts.ConstantBField(acts.Vector3(0, 0, 0 * acts.UnitConstants.T))

    runMaterialValidation(
        args.events,
        args.tracks,
        trackingGeometry,
        decorators,
        field,
        outputDir=os.getcwd(),
        outputName=args.output,
    ).run()

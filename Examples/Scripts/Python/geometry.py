#!/usr/bin/env python3

import os
import json

import acts
from acts import MaterialMapJsonConverter
from acts.examples import (
    WhiteBoard,
    TGeoDetector,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)

import json

def read_json(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
            return data
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except json.JSONDecodeError:
        print(f"Error: The file '{file_path}' is not a valid JSON file.")

def save_json(file_path, data):
    try:
        with open(file_path, 'w', encoding='utf-8') as file:
            json.dump(data, file, indent=4)
    except Exception as e:
        print(f"Error: Unable to save JSON file. {e}")
        
def runGeometry(
    trackingGeometry,
    decorators,
    outputDir,
    events=1,
    outputObj=True,
    outputCsv=True,
    outputJson=True,
):
    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)

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
                outputDir=os.path.join(outputDir, "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=os.path.join(outputDir, "obj")
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            # if not os.path.isdir(outputDir + "/json"):
            #    os.makedirs(outputDir + "/json")
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=os.path.join(outputDir, "json"),
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
                fileName=os.path.join(outputDir, "geometry-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)


if "__main__" == __name__:
    jsonFile="/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geoRuben/tgeoRubenVol.json"
    tgeo_fileName= "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/geoRuben/geometry_Ruben.root"
    logLevel=acts.logging.VERBOSE
    logLevel=acts.logging.INFO
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

    detector, trackingGeometry, decorators = TGeoDetector.create(jsonFile=str(jsonFile),
         fileName=str(tgeo_fileName),
         surfaceLogLevel=customLogLevel(),
         layerLogLevel=customLogLevel(),
         volumeLogLevel=customLogLevel(),
         #mdecorator=matDeco,
     )

    runGeometry(trackingGeometry, decorators, outputDir=os.getcwd())


    file_path = "geometry-map.json"  # Change this to your JSON file path
    json_data = read_json(file_path)
    if json_data:
        surface = json_data["Surfaces"]["entries"]
        for auto in surface:
            if "approach" in auto:
                if auto["approach"] == 2:
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                elif auto["layer"] == 20 and auto["approach"] == 1:
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                elif auto["layer"] == 10 and auto["approach"] == 1:
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                    
    save_json(file_path, json_data)

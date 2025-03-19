// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/MatchingPlotTool.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <TEfficiency.h>


ActsExamples::MatchingPlotTool::MatchingPlotTool(
    const ActsExamples::MatchingPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("MatchingPlotTool", lvl)) {}

void ActsExamples::MatchingPlotTool::book(
    MatchingPlotTool::MatchingPlotCache& matchingPlotCache) const {
  PlotHelpers::Binning bZgen = m_cfg.varBinning.at("Zgen");
  PlotHelpers::Binning bZrec = m_cfg.varBinning.at("Zrec");
  PlotHelpers::Binning bRes = m_cfg.varBinning.at("Res");
  PlotHelpers::Binning bTar = m_cfg.varBinning.at("Tar");
  ACTS_DEBUG("Initialize the histograms for efficiency plots");
  // efficiency vs pT
  matchingPlotCache.resVtxz = PlotHelpers::bookHisto(
      "resVtxz", "Residuals; z_{rec}-z_{gen};Counts", bRes);
  matchingPlotCache.resVtxz_vs_zgen = PlotHelpers::bookHisto(
      "resVtxz_vs_zgen", "Residuals; z_{gen};z_{rec}-z_{gen};Counts", bZgen, bRes);
  matchingPlotCache.zrec_vs_zgen = PlotHelpers::bookHisto(
      "zrec_vs_zgen", "; z_{gen};z_{rec};Counts", bZgen, bZrec);
  matchingPlotCache.eff_vs_zgen = PlotHelpers::bookEff(
      "eff_vs_zgen", "Tracklet vertexing efficiency; z_{gen};Efficiency", bTar);
}

void ActsExamples::MatchingPlotTool::clear(MatchingPlotCache& matchingPlotCache) const {
  delete matchingPlotCache.resVtxz;
  delete matchingPlotCache.resVtxz_vs_zgen;
  delete matchingPlotCache.zrec_vs_zgen;
  delete matchingPlotCache.eff_vs_zgen;
}

void ActsExamples::MatchingPlotTool::write(
    const MatchingPlotTool::MatchingPlotCache& matchingPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  matchingPlotCache.resVtxz->Write();
  matchingPlotCache.resVtxz_vs_zgen->Write();
  matchingPlotCache.zrec_vs_zgen->Write();
  matchingPlotCache.eff_vs_zgen->Write();
}

void ActsExamples::MatchingPlotTool::fill(MatchingPlotTool::MatchingPlotCache& matchingPlotCache, double zrec, double zgen) const {
  PlotHelpers::fillHisto(matchingPlotCache.resVtxz, zrec-zgen, 1);
  PlotHelpers::fillHisto(matchingPlotCache.resVtxz_vs_zgen, zgen, zrec-zgen, 1);
  PlotHelpers::fillHisto(matchingPlotCache.zrec_vs_zgen, zgen, zrec, 1);
  PlotHelpers::fillEff(matchingPlotCache.eff_vs_zgen, zgen, abs(zrec-zgen)<6);
}

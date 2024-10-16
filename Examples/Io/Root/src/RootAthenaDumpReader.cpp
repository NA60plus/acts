// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootAthenaDumpReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <TChain.h>
#include <boost/container/static_vector.hpp>

enum SpacePointType { ePixel = 1, eStrip = 2 };

namespace ActsExamples {

RootAthenaDumpReader::RootAthenaDumpReader(
    const RootAthenaDumpReader::Config& config, Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.inputfile.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.treename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_inputchain = std::make_shared<TChain>(m_cfg.treename.c_str());

  m_outputPixelSpacePoints.initialize(m_cfg.outputPixelSpacePoints);
  m_outputStripSpacePoints.initialize(m_cfg.outputStripSpacePoints);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  m_outputClusters.initialize(m_cfg.outputClusters);

  // Set the branches

  // Set object pointer
  CLhardware = nullptr;
  CLparticleLink_eventIndex = nullptr;
  CLparticleLink_barcode = nullptr;
  CLbarcodesLinked = nullptr;
  CLparticle_charge = nullptr;
  CLphis = nullptr;
  CLetas = nullptr;
  CLtots = nullptr;
  CLlocal_cov = nullptr;
  Part_vParentID = nullptr;
  Part_vParentBarcode = nullptr;
  SPtopStripDirection = nullptr;
  SPbottomStripDirection = nullptr;
  SPstripCenterDistance = nullptr;
  SPtopStripCenterPosition = nullptr;
  TRKproperties = nullptr;
  TRKpattern = nullptr;
  TRKmeasurementsOnTrack_pixcl_sctcl_index = nullptr;
  TRKoutliersOnTrack_pixcl_sctcl_index = nullptr;
  TRKperigee_position = nullptr;
  TRKperigee_momentum = nullptr;
  DTTtrajectory_eventindex = nullptr;
  DTTtrajectory_barcode = nullptr;
  DTTstTruth_subDetType = nullptr;
  DTTstTrack_subDetType = nullptr;
  DTTstCommon_subDetType = nullptr;

  m_inputchain->SetBranchAddress("run_number", &run_number);
  m_inputchain->SetBranchAddress("event_number", &event_number);
  m_inputchain->SetBranchAddress("nSE", &nSE);
  m_inputchain->SetBranchAddress("SEID", SEID);
  m_inputchain->SetBranchAddress("nCL", &nCL);
  m_inputchain->SetBranchAddress("CLindex", CLindex);
  m_inputchain->SetBranchAddress("CLhardware", &CLhardware);
  m_inputchain->SetBranchAddress("CLx", CLx);
  m_inputchain->SetBranchAddress("CLy", CLy);
  m_inputchain->SetBranchAddress("CLz", CLz);
  m_inputchain->SetBranchAddress("CLbarrel_endcap", CLbarrel_endcap);
  m_inputchain->SetBranchAddress("CLlayer_disk", CLlayer_disk);
  m_inputchain->SetBranchAddress("CLeta_module", CLeta_module);
  m_inputchain->SetBranchAddress("CLphi_module", CLphi_module);
  m_inputchain->SetBranchAddress("CLside", CLside);
  m_inputchain->SetBranchAddress("CLmoduleID", CLmoduleID);
  m_inputchain->SetBranchAddress("CLparticleLink_eventIndex",
                                 &CLparticleLink_eventIndex);
  m_inputchain->SetBranchAddress("CLparticleLink_barcode",
                                 &CLparticleLink_barcode);
  m_inputchain->SetBranchAddress("CLbarcodesLinked", &CLbarcodesLinked);
  m_inputchain->SetBranchAddress("CLparticle_charge", &CLparticle_charge);
  m_inputchain->SetBranchAddress("CLphis", &CLphis);
  m_inputchain->SetBranchAddress("CLetas", &CLetas);
  m_inputchain->SetBranchAddress("CLtots", &CLtots);
  m_inputchain->SetBranchAddress("CLloc_direction1", CLloc_direction1);
  m_inputchain->SetBranchAddress("CLloc_direction2", CLloc_direction2);
  m_inputchain->SetBranchAddress("CLloc_direction3", CLloc_direction3);
  m_inputchain->SetBranchAddress("CLJan_loc_direction1", CLJan_loc_direction1);
  m_inputchain->SetBranchAddress("CLJan_loc_direction2", CLJan_loc_direction2);
  m_inputchain->SetBranchAddress("CLJan_loc_direction3", CLJan_loc_direction3);
  m_inputchain->SetBranchAddress("CLpixel_count", CLpixel_count);
  m_inputchain->SetBranchAddress("CLcharge_count", CLcharge_count);
  m_inputchain->SetBranchAddress("CLloc_eta", CLloc_eta);
  m_inputchain->SetBranchAddress("CLloc_phi", CLloc_phi);
  m_inputchain->SetBranchAddress("CLglob_eta", CLglob_eta);
  m_inputchain->SetBranchAddress("CLglob_phi", CLglob_phi);
  m_inputchain->SetBranchAddress("CLeta_angle", CLeta_angle);
  m_inputchain->SetBranchAddress("CLphi_angle", CLphi_angle);
  m_inputchain->SetBranchAddress("CLnorm_x", CLnorm_x);
  m_inputchain->SetBranchAddress("CLnorm_y", CLnorm_y);
  m_inputchain->SetBranchAddress("CLnorm_z", CLnorm_z);
  m_inputchain->SetBranchAddress("CLlocal_cov", &CLlocal_cov);
  m_inputchain->SetBranchAddress("nPartEVT", &nPartEVT);
  m_inputchain->SetBranchAddress("Part_event_number", Part_event_number);
  m_inputchain->SetBranchAddress("Part_barcode", Part_barcode);
  m_inputchain->SetBranchAddress("Part_px", Part_px);
  m_inputchain->SetBranchAddress("Part_py", Part_py);
  m_inputchain->SetBranchAddress("Part_pz", Part_pz);
  m_inputchain->SetBranchAddress("Part_pt", Part_pt);
  m_inputchain->SetBranchAddress("Part_eta", Part_eta);
  m_inputchain->SetBranchAddress("Part_vx", Part_vx);
  m_inputchain->SetBranchAddress("Part_vy", Part_vy);
  m_inputchain->SetBranchAddress("Part_vz", Part_vz);
  m_inputchain->SetBranchAddress("Part_radius", Part_radius);
  m_inputchain->SetBranchAddress("Part_status", Part_status);
  m_inputchain->SetBranchAddress("Part_charge", Part_charge);
  m_inputchain->SetBranchAddress("Part_pdg_id", Part_pdg_id);
  m_inputchain->SetBranchAddress("Part_passed", Part_passed);
  m_inputchain->SetBranchAddress("Part_vProdNin", Part_vProdNin);
  m_inputchain->SetBranchAddress("Part_vProdNout", Part_vProdNout);
  m_inputchain->SetBranchAddress("Part_vProdStatus", Part_vProdStatus);
  m_inputchain->SetBranchAddress("Part_vProdBarcode", Part_vProdBarcode);
  m_inputchain->SetBranchAddress("Part_vParentID", &Part_vParentID);
  m_inputchain->SetBranchAddress("Part_vParentBarcode", &Part_vParentBarcode);
  m_inputchain->SetBranchAddress("nSP", &nSP);
  m_inputchain->SetBranchAddress("SPindex", SPindex);
  m_inputchain->SetBranchAddress("SPx", SPx);
  m_inputchain->SetBranchAddress("SPy", SPy);
  m_inputchain->SetBranchAddress("SPz", SPz);
  m_inputchain->SetBranchAddress("SPCL1_index", SPCL1_index);
  m_inputchain->SetBranchAddress("SPCL2_index", SPCL2_index);
  m_inputchain->SetBranchAddress("SPisOverlap", SPisOverlap);
  m_inputchain->SetBranchAddress("SPradius", SPradius);
  m_inputchain->SetBranchAddress("SPcovr", SPcovr);
  m_inputchain->SetBranchAddress("SPcovz", SPcovz);
  m_inputchain->SetBranchAddress("SPhl_topstrip", SPhl_topstrip);
  m_inputchain->SetBranchAddress("SPhl_botstrip", SPhl_botstrip);
  m_inputchain->SetBranchAddress("SPtopStripDirection", SPtopStripDirection);
  m_inputchain->SetBranchAddress("SPbottomStripDirection",
                                 SPbottomStripDirection);
  m_inputchain->SetBranchAddress("SPstripCenterDistance",
                                 SPstripCenterDistance);
  m_inputchain->SetBranchAddress("SPtopStripCenterPosition",
                                 SPtopStripCenterPosition);

  m_inputchain->SetBranchAddress("nTRK", &nTRK);
  m_inputchain->SetBranchAddress("TRKindex", TRKindex);
  m_inputchain->SetBranchAddress("TRKtrack_fitter", TRKtrack_fitter);
  m_inputchain->SetBranchAddress("TRKparticle_hypothesis",
                                 TRKparticle_hypothesis);
  m_inputchain->SetBranchAddress("TRKproperties", &TRKproperties);
  m_inputchain->SetBranchAddress("TRKpattern", &TRKpattern);
  m_inputchain->SetBranchAddress("TRKndof", TRKndof);
  m_inputchain->SetBranchAddress("TRKmot", TRKmot);
  m_inputchain->SetBranchAddress("TRKoot", TRKoot);
  m_inputchain->SetBranchAddress("TRKchiSq", TRKchiSq);
  m_inputchain->SetBranchAddress("TRKmeasurementsOnTrack_pixcl_sctcl_index",
                                 &TRKmeasurementsOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKoutliersOnTrack_pixcl_sctcl_index",
                                 &TRKoutliersOnTrack_pixcl_sctcl_index);
  m_inputchain->SetBranchAddress("TRKcharge", TRKcharge);
  m_inputchain->SetBranchAddress("TRKperigee_position", &TRKperigee_position);
  m_inputchain->SetBranchAddress("TRKperigee_momentum", &TRKperigee_momentum);
  m_inputchain->SetBranchAddress("TTCindex", TTCindex);
  m_inputchain->SetBranchAddress("TTCevent_index", TTCevent_index);
  m_inputchain->SetBranchAddress("TTCparticle_link", TTCparticle_link);
  m_inputchain->SetBranchAddress("TTCprobability", TTCprobability);
  m_inputchain->SetBranchAddress("nDTT", &nDTT);
  m_inputchain->SetBranchAddress("DTTindex", DTTindex);
  m_inputchain->SetBranchAddress("DTTsize", DTTsize);
  m_inputchain->SetBranchAddress("DTTtrajectory_eventindex",
                                 &DTTtrajectory_eventindex);
  m_inputchain->SetBranchAddress("DTTtrajectory_barcode",
                                 &DTTtrajectory_barcode);
  m_inputchain->SetBranchAddress("DTTstTruth_subDetType",
                                 &DTTstTruth_subDetType);
  m_inputchain->SetBranchAddress("DTTstTrack_subDetType",
                                 &DTTstTrack_subDetType);
  m_inputchain->SetBranchAddress("DTTstCommon_subDetType",
                                 &DTTstCommon_subDetType);

  m_inputchain->Add(m_cfg.inputfile.c_str());
  ACTS_DEBUG("Adding file " << m_cfg.inputfile << " to tree" << m_cfg.treename);

  m_events = m_inputchain->GetEntries();

  ACTS_DEBUG("End of constructor. In total available events=" << m_events);
}  // constructor

ActsExamples::ProcessCode ActsExamples::RootAthenaDumpReader::read(
    const ActsExamples::AlgorithmContext& ctx) {
  ACTS_DEBUG("Reading event " << ctx.eventNumber);
  auto entry = ctx.eventNumber;
  if (entry >= m_events) {
    ACTS_ERROR("event out of bounds");
    return ProcessCode::ABORT;
  }

  std::lock_guard<std::mutex> lock(m_read_mutex);

  m_inputchain->GetEntry(entry);

  // Loop on clusters (measurements)
  ACTS_DEBUG("Found " << nSP << " space points");
  ACTS_DEBUG("Found " << nCL << " clusters / measurements");

  ClusterContainer clusters;
  clusters.resize(nCL);

  // We cannot use im for the index since we might skip measurements
  std::unordered_map<int, std::size_t> imIdxMap;

  for (int im = 0; im < nCL; im++) {
    int bec = CLbarrel_endcap[im];
    int lydisk = CLlayer_disk[im];
    int etamod = CLeta_module[im];
    int phimod = CLphi_module[im];
    int side = CLside[im];
    // ULong64_t moduleID = CLmoduleID     [im];

    ACTS_VERBOSE(bec << " " << lydisk << " " << etamod << " " << phimod << " "
                     << side << " ");

    // Make cluster
    // TODO refactor Cluster class so it is not so tedious
    Cluster cluster;

    const auto& etas = CLetas->at(im);
    const auto& phis = CLetas->at(im);
    const auto& tots = CLtots->at(im);

    const auto totalTot = std::accumulate(tots.begin(), tots.end(), 0);

    const auto [minEta, maxEta] = std::minmax_element(etas.begin(), etas.end());
    const auto [minPhi, maxPhi] = std::minmax_element(phis.begin(), phis.end());

    cluster.sizeLoc0 = *maxEta - *minEta;
    cluster.sizeLoc1 = *maxPhi - *minPhi;

    if (totalTot == 0.0) {
      ACTS_VERBOSE("total time over threshold is 0, set all activations to 0");
      nTotalTotZero++;
    }

    for (const auto& [eta, phi, tot] : Acts::zip(etas, phis, tots)) {
      // Make best out of what we have:
      // Weight the overall collected charge corresponding to the
      // time-over-threshold of each cell Use this as activation (does this make
      // sense?)
      auto activation =
          (totalTot != 0.0) ? CLcharge_count[im] * tot / totalTot : 0.0;

      // This bases every cluster at zero, but shouldn't matter right now
      ActsFatras::Segmentizer::Bin2D bin;
      bin[0] = eta - *minEta;
      bin[1] = phi - *minPhi;

      // Of course we have no Segment2D because this is not Fatras
      cluster.channels.emplace_back(bin, ActsFatras::Segmentizer::Segment2D{},
                                    activation);
    }

    ACTS_VERBOSE("Cluster " << im << ": " << cluster.channels.size()
                            << "cells, dimensions: " << cluster.sizeLoc0 << ", "
                            << cluster.sizeLoc1);

    clusters[im] = cluster;
  }

  if (measurements.size() < static_cast<std::size_t>(nCL)) {
    ACTS_WARNING("Could not convert " << nCL - measurements.size() << " / "
                                      << nCL << " measurements");
  }

  if (particles.size() - prevParticlesSize > 0) {
    ACTS_DEBUG("Created " << particles.size() - prevParticlesSize
                          << " dummy particles");
  }

  if (nTotalTotZero > 0) {
    ACTS_WARNING(nTotalTotZero << " / " << nCL
                               << " clusters have zero time-over-threshold");
  }

  return {clusters, measurements, measPartMap, imIdxMap};
}

std::tuple<SimSpacePointContainer, SimSpacePointContainer,
           SimSpacePointContainer>
RootAthenaDumpReader::readSpacepoints(
    const std::optional<std::unordered_map<int, std::size_t>>& imIdxMap) const {
  SimSpacePointContainer pixelSpacePoints;
  pixelSpacePoints.reserve(nSP);

  SimSpacePointContainer stripSpacePoints;
  stripSpacePoints.reserve(nSP);

  SimSpacePointContainer spacePoints;
  spacePoints.reserve(nSP);

  // Loop on space points
  for (int isp = 0; isp < nSP; isp++) {
    Acts::Vector3 globalPos{SPx[isp], SPy[isp], SPz[isp]};
    double sp_covr = SPcovr[isp];
    double sp_covz = SPcovz[isp];

    // PIX=1  STRIP = 2
    auto type = SPCL2_index[isp] == -1 ? ePixel : eStrip;

    ACTS_VERBOSE("SP:: " << type << " [" << globalPos.transpose() << "] "
                         << spCovr << " " << spCovz);

    boost::container::static_vector<Acts::SourceLink, 2> sLinks;

    const auto cl1Index = SPCL1_index[isp];
    assert(cl1Index >= 0 && cl1Index < nCL);

    auto getGeoId =
        [&](auto athenaId) -> std::optional<Acts::GeometryIdentifier> {
      if (m_cfg.geometryIdMap == nullptr) {
        return Acts::GeometryIdentifier{athenaId};
      }
      if (m_cfg.geometryIdMap->left.find(athenaId) ==
          m_cfg.geometryIdMap->left.end()) {
        return std::nullopt;
      }
      return m_cfg.geometryIdMap->left.at(athenaId);
    };

    auto cl1GeoId = getGeoId(CLmoduleID[cl1Index]);
    if (!cl1GeoId) {
      ACTS_WARNING("Could not find geoId for spacepoint cluster 1");
      continue;
    }

    if (imIdxMap && !imIdxMap->contains(cl1Index)) {
      ACTS_WARNING("Measurement 1 for spacepoint " << isp << " not created");
      continue;
    }

    IndexSourceLink first(*cl1GeoId,
                          imIdxMap ? imIdxMap->at(cl1Index) : cl1Index);
    sLinks.emplace_back(first);

    // First create pixel spacepoint here, later maybe overwrite with strip
    // spacepoint
    SimSpacePoint sp(globalPos, std::nullopt, spCovr, spCovz, std::nullopt,
                     sLinks);

    if (type == ePixel) {
      pixelSpacePoints.push_back(sp);
    } else {
      const auto cl2Index = SPCL2_index[isp];
      assert(cl2Index >= 0 && cl2Index < nCL);

      auto cl2GeoId = getGeoId(CLmoduleID[cl1Index]);
      if (!cl2GeoId) {
        ACTS_WARNING("Could not find geoId for spacepoint cluster 2");
        continue;
      }

      if (imIdxMap && !imIdxMap->contains(cl2Index)) {
        ACTS_WARNING("Measurement 2 for spacepoint " << isp << " not created");
        continue;
      }

      IndexSourceLink second(*cl2GeoId,
                             imIdxMap ? imIdxMap->at(cl2Index) : cl2Index);
      sLinks.emplace_back(second);

      using Vector3f = Eigen::Matrix<float, 3, 1>;
      const Vector3f topStripDirection{SPtopStripDirection->at(isp).at(0),
                                       SPtopStripDirection->at(isp).at(1),
                                       SPtopStripDirection->at(isp).at(2)};
      const Vector3f bottomStripDirection{
          SPbottomStripDirection->at(isp).at(0),
          SPbottomStripDirection->at(isp).at(1),
          SPbottomStripDirection->at(isp).at(2)};
      const Vector3f stripCenterDistance{SPstripCenterDistance->at(isp).at(0),
                                         SPstripCenterDistance->at(isp).at(1),
                                         SPstripCenterDistance->at(isp).at(2)};
      const Vector3f topStripCenterPosition{
          SPtopStripCenterPosition->at(isp).at(0),
          SPtopStripCenterPosition->at(isp).at(1),
          SPtopStripCenterPosition->at(isp).at(2)};

      sp = SimSpacePoint(globalPos, std::nullopt, spCovr, spCovz, std::nullopt,
                         sLinks, SPhl_topstrip[isp], SPhl_botstrip[isp],
                         topStripDirection.cast<double>(),
                         bottomStripDirection.cast<double>(),
                         stripCenterDistance.cast<double>(),
                         topStripCenterPosition.cast<double>());

      stripSpacePoints.push_back(sp);
    }

    spacePoints.push_back(sp);
  }

  ACTS_DEBUG("Created " << pixelSpacePoints.size() << " "
                        << " pixel space points");

  ACTS_DEBUG("Created " << spacePoints.size() << " overall space points");

  m_outputPixelSpacePoints(ctx, std::move(pixelSpacePoints));
  m_outputStripSpacePoints(ctx, std::move(stripSpacePoints));
  m_outputSpacePoints(ctx, std::move(spacePoints));
  m_outputClusters(ctx, std::move(clusters));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples

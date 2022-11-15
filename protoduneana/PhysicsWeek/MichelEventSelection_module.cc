///////////////////////////////////////////////////////////////////////
// Class:       MichelReco
// Module Type: analyzer
// File:        MichelReco_test.cc
// 
// Michel electron event selection based on finding a cluster of 
// Michel-like hits near the end of reconstructed tracks 
// clusters are required in all three planes in order to select
// an event as from a Michel electron
//
// Generated at Tue Mar 21 06:06:49 2017 by Aiden Reynolds using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////

#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/Simulation/LArVoxelData.h"
#include "larsim/Simulation/LArVoxelList.h"
#include "larsim/Simulation/SimListUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardata/Utilities/DatabaseUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/ArtDataHelper/MVAReader.h"

#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <fstream>
#include <tuple>

// Maximum number of tracks
//constexpr int kMaxTrack = 1000; // unused

namespace MichelReco {

class MichelReco;

class MichelReco : public art::EDAnalyzer {
public:
  explicit MichelReco(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelReco(MichelReco const &) = delete;
  MichelReco(MichelReco &&) = delete;
  MichelReco & operator = (MichelReco const &) = delete;
  MichelReco & operator = (MichelReco &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions
  void beginJob() override;
  void endJob() override;
  void beginRun(const art::Run& run) override;
  void reconfigure(fhicl::ParameterSet const& p) ;

  // My functions

  // check if a position is inside the fiducial volume
  bool insideFidVol(geo::Point_t const&);
  
  // check if a hit is close to the projected 2d end of a track
  bool hitCloseToTrackEnd(detinfo::DetectorPropertiesData const& detProp,
                          double radius, geo::Point_t const& end,
                          double end2D[2], recob::Hit hit, 
                          geo::GeometryCore const & geom);  

  // match a reco track to the true particles track id
  int trackMatching (int trackIndex, 
                     art::FindManyP<recob::Hit> hitsFromTracks);

  // check if a true particle is a decaying muon
  bool isMuonDecaying (simb::MCParticle particle, 
                       std::vector<simb::MCParticle> particles);
  
  //  takes vector of hits and checks that the majority
  //  of the charge in these hits is from a michel 
  bool areHitsMichel(detinfo::DetectorClocksData const& clockData,
                     const std::vector<recob::Hit> & hits );
  
  calo::CalorimetryAlg fCalorimetryAlg;

private:
  size_t fEvNumber;

  TCanvas * c;
  
  // Reconstruction parameters
  double fRadiusThreshold;       // Radius for hits in event selection
  int fNumberThreshold;            // Selection threshold for number of close hits 
  double fCNNThreshold;            // CNN output threshold for selection
  
  int fNumberCloseHitsStart[3];               // Number of hits within selection radius
  int fNumberCloseHitsEnd[3];               // Number of hits within selection radius
  
  // Purity and Efficiency Numbers
  double fYesSelected;       // Number of correctly selected events for each set of reco params
  double fNoSelected;        // Number of wrongly selected events "      "
  int fNMichel;

  // fhicl parameters
  int fBestview;                        // Best plane, default collection
  art::InputTag fTrackModuleLabel;      // Module label from track recnstrucion
  art::InputTag fNNetModuleLabel;       // Module label for CNN
  art::InputTag fParticleModuleLabel;   // Module label for particle simulation
  double fFidVolCut;                    // Size of cut to select a fiducial volume
};


MichelReco::MichelReco(fhicl::ParameterSet const & p)
  : 
  EDAnalyzer(p),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
  reconfigure(p);
}

void MichelReco::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
}

void MichelReco::endJob()
{
  mf::LogVerbatim("MichelReco") << "MichelReco finished job";
  
  std::ofstream michelRecoOut;
  michelRecoOut.open("test.txt");
  michelRecoOut << "nMichel Purity Efficiency" << std::endl;

  double purity =  fYesSelected / ( fYesSelected + fNoSelected );
  double efficiency = fYesSelected/ fNMichel;

  michelRecoOut << fNMichel << " " << purity << " " << efficiency << std::endl;
  michelRecoOut.close();
}

void MichelReco::beginRun(const art::Run&)
{
  art::ServiceHandle<sim::LArG4Parameters> larParameters;
}

void MichelReco::analyze(art::Event const & evt)
{

  fEvNumber = evt.id().event();
  std::cout << "Michel event selection is on event " << fEvNumber << std::endl;

  // find all decaying muons in truth
  auto const & particles = *evt.getValidHandle<std::vector<simb::MCParticle>>( fParticleModuleLabel );
  for ( auto const & particle : particles ) {
    if ( abs(particle.PdgCode()) == 13 ) {
      geo::Point_t const mcEnd{ particle.EndX(), particle.EndY(), particle.EndZ() };
      if ( isMuonDecaying( particle, particles ) && insideFidVol( mcEnd ) ) {
        fNMichel += 1;
      }
    }
  }
  
  // Get the reconstructed objects, tracks and hits
  auto const & trackHandle = evt.getValidHandle<std::vector<recob::Track>>( fTrackModuleLabel );
  auto const & tracks = * trackHandle;
  art::FindManyP<recob::Hit> hitsFromTracks( trackHandle, evt, fTrackModuleLabel );

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  // loop over tracks
  for (auto const & track : tracks) {
    
    // find best true trask ID for this track
    //int bestTrackId = trackMatching( &track - &tracks[0], hitsFromTracks );
    
    // find track start and end
    auto const [start_pt, start] = std::make_tuple(track.Vertex(), track.Vertex<TVector3>());
    auto const [end_pt, end] = std::make_tuple(track.End(), track.End<TVector3>());

    // Check if the track starts or ends in the fiducial volume
    bool startInFidVol = insideFidVol( start_pt );
    bool endInFidVol = insideFidVol( end_pt );
    if ( startInFidVol || endInFidVol ) {

      // reset hit counts
      for (int plane = 0; plane < 3; plane++) {
        fNumberCloseHitsStart[ plane ] = 0; fNumberCloseHitsEnd[ plane ] = 0;
      }

      // get the hit results helper
      anab::MVAReader<recob::Hit,4> hitResults( evt, fNNetModuleLabel );

      // loop over hits store the selected hits
      std::vector< recob::Hit > taggedHitsStart;
      std::vector< recob::Hit > taggedHitsEnd;
      for ( size_t h = 0; h < hitResults.size(); ++h ) {
        
        // cnn output vector: [em,trk,none,michel]
        std::array<float,4> cnn_out = hitResults.getOutput( h );
        
        // keep hits with large CNN values
        if ( cnn_out[hitResults.getIndex("michel")] > fCNNThreshold ) {
          
          // need to check at both ends of the track
          if ( startInFidVol ) {
            
            // In each plane check if the michel tagged hits are close to the end of the track
            for ( int plane = 0; plane < 3; plane++ ) {
              
              // project the track onto the plane
              geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
              auto const & trackStart2D = pma::GetProjectionToPlane(start,
                                                                    plane,
                                                                    geom.FindTPCAtPosition(start_pt).TPC,
                                                                    geom.PositionToCryostatID(start_pt).Cryostat);
              double start2D[2] = {trackStart2D.X(), trackStart2D.Y()};

              // get the hit itself
              const recob::Hit & hit = hitResults.item(h);
              if (hit.View() == plane) {

                // check that the hit is close to the track endpoint
                // and add to tagged hits if so
                if ( hitCloseToTrackEnd( detProp, fRadiusThreshold, start_pt, start2D, hit, geom ) ) {
                  fNumberCloseHitsStart[ plane ] += 1;
                  if ( plane == 2 ) { taggedHitsStart.push_back( hit ); };
                }

              }

            }

          }
          if ( endInFidVol ) {

            // In each plane check if the michel tagged hits are close to the end of the track
            for ( int plane = 0; plane < 3; plane++ ) {
             
              // project the track onto the plane 
              geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
              auto const & trackEnd2D = pma::GetProjectionToPlane(end,
                                                                  plane,
                                                                  geom.FindTPCAtPosition(end_pt).TPC,
                                                                  geom.PositionToCryostatID(end_pt).Cryostat);
              double end2D[2] = {trackEnd2D.X(), trackEnd2D.Y()};

              // get the hit itself
              const recob::Hit & hit = hitResults.item(h);
              if (hit.View() == plane) {

                // check that the hit is close to the track endpoint
                // and add to tagged hits if so
                if ( hitCloseToTrackEnd( detProp, fRadiusThreshold, end_pt, end2D, hit, geom ) ) {
                  fNumberCloseHitsEnd[ plane ] += 1;
                  if ( plane == 2 ) { taggedHitsEnd.push_back( hit ); }
                }

              }

            }

          }

        }

      } // end of loop over hits

      // event selection decision: clusters of hits in all planes at one end of the track
      bool startSelected = ( fNumberCloseHitsStart[0] > fNumberThreshold && fNumberCloseHitsStart[1] > fNumberThreshold && fNumberCloseHitsStart[2] > fNumberThreshold );
      bool endSelected = ( fNumberCloseHitsEnd[0] > fNumberThreshold && fNumberCloseHitsEnd[1] > fNumberThreshold && fNumberCloseHitsEnd[2] > fNumberThreshold );
      bool eventSelected = ( startSelected || endSelected );
      if ( startSelected && endSelected ) { eventSelected = false; }
      
      // check if the event was correclty tagged
      //  i.e. do the tagged hits correspond to a michel 
      bool particleIsMichel(false);
      if (eventSelected) {
        if ( startSelected ) { particleIsMichel = areHitsMichel(clockData, taggedHitsStart); }
        if ( endSelected ) { particleIsMichel = areHitsMichel(clockData, taggedHitsEnd); }
      }
      
      // update purity and efficiecny numbers and draw example events
      if ( particleIsMichel && eventSelected ) {
        
        fYesSelected += 1; 
        
        // plot the region around track end to see what the event is
        // if ( startSelected ) {
        //   
        //   for (int plane = 0; plane < 3; plane++) { 
        //   
        //     geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
        //     auto const & trackStart2D = pma::GetProjectionToPlane(trackStart, plane, geom.FindTPCAtPosition(start).TPC, geom.FindCryostatAtPosition(start));
        //     double start2D[2] = {trackStart2D.X(), trackStart2D.Y()};

        //     TH2D * histCNN = new TH2D("cnn", "Correctly Tagged", 200, start2D[0] - 20.0, start2D[0] +20.0, 200, start2D[1] - 20.0, start2D[1] + 20.0); 
        //     TH2D * histCharge = new TH2D("charge", "Correctly Tagged", 200, start2D[0] - 20.0, start2D[0] + 20.0, 200, start2D[1] - 20.0, start2D[1] + 20.0); 
        //     histCNN -> SetStats(0);
        //     histCharge -> SetStats(0);

        //     for ( size_t h = 0; h< hitResults.size(); ++h ) {
        //       const recob::Hit & hit = hitResults.item(h);
        //       auto const & hitLocation = pma::WireDriftToCm(hit.WireID().Wire, hit.PeakTime(), plane, geom.FindTPCAtPosition(start).TPC, geom.FindCryostatAtPosition(start)); 
        //     
        //       if (hit.View() == plane) {
        //         histCNN -> Fill(hitLocation.X(), hitLocation.Y(), hitResults.getOutput(h)[hitResults.getIndex("michel")]);
        //         histCharge -> Fill(hitLocation.X(), hitLocation.Y(), hit.Integral());
        //       }
        //     } 
        //   
        //     histCNN -> Draw("colz");
        //     c->SaveAs(Form("imgs/correctCNN%d_plane%d.png",bestTrackId, plane));
        //     histCharge-> Draw("colz");
        //     c->SaveAs(Form("imgs/correctCharge%d_plane%d.png",bestTrackId, plane));
        // 
        //   }
        // }
        // if ( endSelected ) {
        //   
        //   for (int plane = 0; plane < 3; plane++) { 
        //   
        //     geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
        //     auto const & trackEnd2D = pma::GetProjectionToPlane(trackEnd, plane, geom.FindTPCAtPosition(end).TPC, geom.FindCryostatAtPosition(end));
        //     double end2D[2] = {trackEnd2D.X(), trackEnd2D.Y()};

        //     TH2D * histCNN = new TH2D("cnn", "Correctly Tagged", 200, end2D[0] - 20.0, end2D[0] +20.0, 200, end2D[1] - 20.0, end2D[1] + 20.0); 
        //     TH2D * histCharge = new TH2D("charge", "Correctly Tagged", 200, end2D[0] - 20.0, end2D[0] + 20.0, 200, end2D[1] - 20.0, end2D[1] + 20.0); 
        //     histCNN -> SetStats(0);
        //     histCharge -> SetStats(0);

        //     for ( size_t h = 0; h< hitResults.size(); ++h ) {
        //       const recob::Hit & hit = hitResults.item(h);
        //       auto const & hitLocation = pma::WireDriftToCm(hit.WireID().Wire, hit.PeakTime(), plane, geom.FindTPCAtPosition(end).TPC, geom.FindCryostatAtPosition(end)); 
        //     
        //       if (hit.View() == plane) {
        //         histCNN -> Fill(hitLocation.X(), hitLocation.Y(), hitResults.getOutput(h)[hitResults.getIndex("michel")]);
        //         histCharge -> Fill(hitLocation.X(), hitLocation.Y(), hit.Integral());
        //       }
        //     } 
        //   
        //     histCNN -> Draw("colz");
        //     c->SaveAs(Form("imgs/correctCNN%d_plane%d.png",bestTrackId, plane));
        //     histCharge-> Draw("colz");
        //     c->SaveAs(Form("imgs/correctCharge%d_plane%d.png",bestTrackId, plane));
        // 
        //   }

        // }

      }
      if ( !particleIsMichel && eventSelected ) {
        
        fNoSelected += 1; 
        
        // plot the region around track end to see what the event is
        // if ( startSelected ) {
        //   
        //   for (int plane = 0; plane < 3; plane++) { 
        //   
        //     geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
        //     auto const & trackStart2D = pma::GetProjectionToPlane(trackStart, plane, geom.FindTPCAtPosition(start).TPC, geom.FindCryostatAtPosition(start));
        //     double start2D[2] = {trackStart2D.X(), trackStart2D.Y()};

        //     TH2D * histCNN = new TH2D("cnn", "Incorrectly Tagged", 200, start2D[0] - 20.0, start2D[0] +20.0, 200, start2D[1] - 20.0, start2D[1] + 20.0); 
        //     TH2D * histCharge = new TH2D("charge", "Incorrectly Tagged", 200, start2D[0] - 20.0, start2D[0] + 20.0, 200, start2D[1] - 20.0, start2D[1] + 20.0); 
        //     histCNN -> SetStats(0);
        //     histCharge -> SetStats(0);

        //     for ( size_t h = 0; h< hitResults.size(); ++h ) {
        //       const recob::Hit & hit = hitResults.item(h);
        //       auto const & hitLocation = pma::WireDriftToCm(hit.WireID().Wire, hit.PeakTime(), plane, geom.FindTPCAtPosition(start).TPC, geom.FindCryostatAtPosition(start)); 
        //     
        //       if (hit.View() == plane) {
        //         histCNN -> Fill(hitLocation.X(), hitLocation.Y(), hitResults.getOutput(h)[hitResults.getIndex("michel")]);
        //         histCharge -> Fill(hitLocation.X(), hitLocation.Y(), hit.Integral());
        //       }
        //     } 
        //   
        //     histCNN -> Draw("colz");
        //     c->SaveAs(Form("imgs/incorrectCNN%d_plane%d.png",bestTrackId, plane));
        //     histCharge-> Draw("colz");
        //     c->SaveAs(Form("imgs/incorrectCharge%d_plane%d.png",bestTrackId, plane));
        // 
        //   }
        // }
        // if ( endSelected ) {
        //   
        //   for (int plane = 0; plane < 3; plane++) { 
        //   
        //     geo::GeometryCore const & geom = *art::ServiceHandle<geo::Geometry>();
        //     auto const & trackEnd2D = pma::GetProjectionToPlane(trackEnd, plane, geom.FindTPCAtPosition(end).TPC, geom.FindCryostatAtPosition(end));
        //     double end2D[2] = {trackEnd2D.X(), trackEnd2D.Y()};

        //     TH2D * histCNN = new TH2D("cnn", "Incorrectly Tagged", 200, end2D[0] - 20.0, end2D[0] +20.0, 200, end2D[1] - 20.0, end2D[1] + 20.0); 
        //     TH2D * histCharge = new TH2D("charge", "Incorrectly Tagged", 200, end2D[0] - 20.0, end2D[0] + 20.0, 200, end2D[1] - 20.0, end2D[1] + 20.0); 
        //     histCNN -> SetStats(0);
        //     histCharge -> SetStats(0);

        //     for ( size_t h = 0; h< hitResults.size(); ++h ) {
        //       const recob::Hit & hit = hitResults.item(h);
        //       auto const & hitLocation = pma::WireDriftToCm(hit.WireID().Wire, hit.PeakTime(), plane, geom.FindTPCAtPosition(end).TPC, geom.FindCryostatAtPosition(end)); 
        //     
        //       if (hit.View() == plane) {
        //         histCNN -> Fill(hitLocation.X(), hitLocation.Y(), hitResults.getOutput(h)[hitResults.getIndex("michel")]);
        //         histCharge -> Fill(hitLocation.X(), hitLocation.Y(), hit.Integral());
        //       }
        //     } 
        //   
        //     histCNN -> Draw("colz");
        //     c->SaveAs(Form("imgs/incorrectCNN%d_plane%d.png",bestTrackId, plane));
        //     histCharge-> Draw("colz");
        //     c->SaveAs(Form("imgs/incorrectCharge%d_plane%d.png",bestTrackId, plane));
        // 
        //   }
        // 
        // }

      }

    }

  } // end of loop over tracks

}
// // Checks if the hits from a given track match the track with a given index
// int MichelReco::trackMatching( int trackIndex, art::FindManyP<recob::Hit> hitsFromTracks )
// {
//   std::map<int,double> trackID_E;
//   art::ServiceHandle<cheat::BackTrackerService> bt_serv;

//   for (size_t h = 0; h < hitsFromTracks.at(trackIndex).size(); ++h)
//   {
//     for (auto const & id : bt_serv->HitToTrackIDEs(hitsFromTracks.at(trackIndex)[h]))
//     {
//       trackID_E[id.trackID] += id.energy;
//     }
//   }

//   double max_e = 0.0; double tot_e = 0.0;
//   int best_id = 0;
//   for (std::map<int,double>::iterator it = trackID_E.begin(); it != trackID_E.end(); ++it)
//   {
//     tot_e += it->second;
//     if (it->second > max_e)
//     {
//       max_e = it->second;
//       best_id = it->first;
//     }
//   }

//   if ((max_e > 0.0) && (true/*max_e > 0.5*trackID_sumen[best_id]*/) && (tot_e > 0.0))
//   {
//     return best_id;
//   }
//   else
//   {
//     return -999;
//   }

// }

// checks if the particle responsible for a track has an associated muon decay
bool MichelReco::isMuonDecaying ( simb::MCParticle particle, std::vector<simb::MCParticle> particles ) { 

  bool hasElectron = false, hasNuMu = false, hasNuE = false;

  unsigned int nSec = particle.NumberDaughters();
  for (size_t d = 0; d < nSec; ++d)
  {
    for (auto const & daughter : particles)
    {
      if (daughter.Mother() == particle.TrackId()) 
      {
        int d_pdg = abs(daughter.PdgCode());
        if (d_pdg == 11) hasElectron = true;
        else if (d_pdg == 14) hasNuMu = true;
        else if (d_pdg == 12) hasNuE = true;
      }
    }
  }
  return (hasElectron && hasNuMu && hasNuE);

}

// checks if a position is inside the fiducial volume
bool MichelReco::insideFidVol( geo::Point_t const& pos ) {

  geo::GeometryCore const& geom = *art::ServiceHandle<geo::Geometry>();
  bool inside = false;

  geo::TPCID idtpc = geom.FindTPCAtPosition(pos);

  if (geom.HasTPC(idtpc)) {

    const geo::TPCGeo& tpcgeo = geom.GetElement(idtpc);
    double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
    double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
    double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();
    
    for (auto const& tpcg : geom.Iterate<geo::TPCGeo>()) {
        if (tpcg.MinX() < minx) minx = tpcg.MinX();
        if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
        if (tpcg.MinY() < miny) miny = tpcg.MinY();
        if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
        if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
        if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
    }

    double dista = fabs(minx - pos.X());
    double distb = fabs(pos.X() - maxx);
    if ((pos.X() > minx) && (pos.X() < maxx) &&
        (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
      
    dista = fabs(miny - pos.Y());
    distb = fabs(pos.Y()-maxy);
    if (inside && (pos.Y() > miny) && (pos.Y() < maxy) &&
        (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
    else inside = false;
    
    dista = fabs(minz - pos.Z());
    distb = fabs(pos.Z() - maxz);
    if (inside && (pos.Z() > minz) && (pos.Z() < maxz) &&
        (dista > fFidVolCut) && (distb > fFidVolCut)) inside = true;
    else inside = false;

  }
  return inside;
}

// checks if a 2d hit is close to a given position
bool MichelReco::hitCloseToTrackEnd(detinfo::DetectorPropertiesData const& detProp,
                                    double radius, geo::Point_t const& end, double end2D[2], recob::Hit hit, geo::GeometryCore const & geom ) {

  bool close = false;
 
  auto const & hitLocation = pma::WireDriftToCm(detProp,
                                                hit.WireID().Wire, hit.PeakTime(), hit.View(), geom.FindTPCAtPosition(end).TPC, geom.PositionToCryostatID(end).Cryostat);

  double deltaXToHit = hitLocation.X() - end2D[0];
  double deltaYToHit = hitLocation.Y() - end2D[1];

  double displacementToHit = pow(pow(deltaXToHit,2)+pow(deltaYToHit,2),0.5);

  if (displacementToHit < radius) close = true; 

  return close;

}

// checks if vector of hits come mostly from michel charge
bool MichelReco::areHitsMichel(detinfo::DetectorClocksData const& clockData,
                               const std::vector< recob::Hit > & hits ) {

  const simb::MCParticle *  mcParticle = 0;

  art::ServiceHandle< cheat::BackTrackerService > bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  std::unordered_map< int, double > trkIDE;
  for ( auto const & hit : hits ) {
    for ( auto const & ide : bt_serv->HitToTrackIDEs(clockData, hit) ) { trkIDE[ide.trackID] += ide.energy; }
  }
  
  int best_id(0);
  double tot_e(0), max_e(0);
  for ( auto const & contrib : trkIDE ) {

    tot_e += contrib.second;
    if ( contrib.second > max_e ) {
      max_e = contrib.second;
      best_id = contrib.first;
    }
  }

  if ( (max_e > 0) && (tot_e > 0) ) {
    if ( best_id < 0 ) {
      best_id = -best_id;
    }
    mcParticle = pi_serv->TrackIdToParticle_P( best_id );
  }
 
  if (mcParticle != 0) { 
    return ( abs(mcParticle->PdgCode()) && ( mcParticle->Process() == "Decay" || mcParticle->Process() == "muMinusCaptureAtRest") ) ;
  }
  else {
    std::cout << "No match found" << hits.size() << std::endl;
    return false; 
  }

}

void MichelReco::reconfigure( fhicl::ParameterSet const& p ) {

  c = new TCanvas( "canv", "canv" );

  fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
  fNNetModuleLabel = p.get<std::string>("NNetModuleLabel");
  fParticleModuleLabel = p.get<std::string>("ParticleModuleLabel");
  
  fBestview = p.get<int>("Bestview");

  fRadiusThreshold = p.get<double>("MichelCloseHitRadius");
  fNumberThreshold = p.get<int>("CloseHitsThreshold"); 
  fCNNThreshold = p.get<double>("CNNThreshSelect");
  fFidVolCut = p.get<double>("FidVolCut"); 
  
  for ( int plane = 0; plane < 3; plane++ ) { 
    fNumberCloseHitsStart[ plane ] = 0; fNumberCloseHitsEnd[ plane ] = 0;
  }

  fYesSelected = 0;
  fNoSelected = 0;
  fNMichel = 0;

  return;
}

} // MichelReco namespace

DEFINE_ART_MODULE(MichelReco::MichelReco)

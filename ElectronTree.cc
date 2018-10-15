
// -*- C++ -*-
//
// Package:    
// Class:      ElectronTree
// 
/**\class ElectronTree ElectronTree.cc Demo/ElectronTree/src/ElectronTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Steenu Johnson
//         Created:  Mon May 28 20:56:30 CEST 2018
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//
// class decleration
//

struct Lepton{
  TLorentzVector v;
  int Charge;
  //  float trkIso;
  int LT; //Loose=1 Tight=2
};
std::vector<Lepton> Muons;
std::vector<Lepton> Electrons;

class ElectronTree : public edm::EDAnalyzer {
public:
  explicit ElectronTree(const edm::ParameterSet&);
  ~ElectronTree();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void ClearVariables();
  void Sort(int opt);
  double electronEA(reco::GsfElectron const& it);
  double electronPFIsolation(reco::GsfElectron const& it, double rhoValue_);
  bool electronMediumCuts(reco::GsfElectron const& it, const reco::Vertex PV_);
  bool electronTriggerCuts(reco::GsfElectron const& it);
  const reco::Candidate* getmother(const reco::Candidate* particle);
  void FillMothers(const reco::Candidate* particle);
  // ----------member data ---------------------------
  unsigned int minTracks_;
  bool isMC;
  TH1D *demohisto;
  TTree *OutputTree;
  int NTowers;
  double TowerEta[5000];
  double TowerPhi[5000];
  double TowerEmEt[5000];
  double TowerEmEnergy[5000];
  double TowerEnergy[5000];
  double TowerHadEt[5000];
  double TowerHadEnergy[5000];
  int NElectrons;
  double ElectronPt[1000];
  double ElectronEta[1000];
  double ElectronPhi[1000];
  double ElectronTrackIso[1000];
  double ElectronEcalIso[1000];
  double ElectronHcalIso[1000];
  double ElectrondeltaEtaSuperClusterTrackAtVtx[1000];
  double ElectrondeltaPhiSuperClusterTrackAtVtx[1000];
  double ElectronsigmaIetaIeta[1000];
  double ElectronhadronicOverEm[1000];
  double ElectronFbrem[1000];
  double ElectronEoverP[1000];
  double ElectronecalEnergy[1000];
  double Electronp_in[1000];
  double ElectronDxy[1000];
  double ElectronDz[1000];
  double ElectronnumberOfLostHits[1000];
  double ElectronsuperClustereta[1000];
  double ElectronPFIsolation[1000];
  bool ElectronIsMedium[1000];
  bool ElectronIsTrigger[1000];
  int NJets;
  double JetPt[1000];
  double JetEta[1000];
  double JetPhi[1000];
  double JetPx[1000];
  double JetPy[1000];
  double JetPz[1000];
  double JetEnergy[1000];
  double JetEmEnergy[1000];
  double JetHadronEnergy[1000];
  int NMuons;
  double MuonPt[1000];
  double MuonEta[1000];
  double MuonPhi[1000];
  double MuonEnergy[1000];
  double METPt;
  double METPhi;
  bool isZfromMu;
  bool isZfromE;
  int NMC;
  double MCPt[5000];
  double MCEta[5000];
  double MCPhi[5000];
  double MCEnergy[5000];
  double MCCharge[5000];
  double MCMass[5000];
  int MCIndex[5000];
  int MCId[5000];
  int MCStatus[5000];
  int NMCmother[5000];
  int MCMotherIndex[5000];
  int NMCdaughter[5000];
  //math::PtEtaPhiMLorentzVector TowerVector[5000];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

ElectronTree::ElectronTree(const edm::ParameterSet& iConfig) :
  minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0)),
  isMC(iConfig.getParameter<bool>("MC"))
{

  //now do what ever initialization is needed

     
  edm::Service<TFileService> fs;
  demohisto = fs->make<TH1D>("towers" , "Towers" , 100 , 0 , 5000 );
  OutputTree = fs->make<TTree>("RecoTree", "Reco Tree version 0");

  //initialize Tree Branches

  OutputTree->Branch("NTowers",&NTowers,"NTowers/I");
  OutputTree->Branch("TowerEta",&TowerEta,"TowerEta[NTowers]/D");
  OutputTree->Branch("TowerPhi",&TowerPhi,"TowerPhi[NTowers]/D");
  OutputTree->Branch("TowerEmEt",&TowerEmEt,"TowerEmEt[NTowers]/D");
  OutputTree->Branch("TowerEmEnergy",&TowerEmEnergy,"TowerEmEnergy[NTowers]/D");
  OutputTree->Branch("TowerHadEt",&TowerHadEt,"TowerHadEt[NTowers]/D");
  OutputTree->Branch("TowerHadEnergy",&TowerHadEnergy,"TowerHadEnergy[NTowers]/D");
  OutputTree->Branch("TowerEnergy",&TowerEnergy,"TowerEnergy[NTowers]/D");
  OutputTree->Branch("NElectrons",&NElectrons,"NElectrons/I");
  OutputTree->Branch("ElectronPt",&ElectronPt,"ElectronPt[NElectrons]/D");
  OutputTree->Branch("ElectronEta",&ElectronEta,"ElectronEta[NElectrons]/D");
  OutputTree->Branch("ElectronPhi",&ElectronPhi,"ElectronPhi[NElectrons]/D");
  OutputTree->Branch("ElectronTrackIso", &ElectronTrackIso, "ElectronTrackIso[NElectrons]/D");
  OutputTree->Branch("ElectronEcalIso", &ElectronEcalIso, "ElectronEcalIso[NElectrons]/D");
  OutputTree->Branch("ElectronHcalIso", &ElectronHcalIso, "ElectronHcalIso[NElectrons]/D");
  OutputTree->Branch("ElectrondeltaEtaSuperClusterTrackAtVtx",&ElectrondeltaEtaSuperClusterTrackAtVtx,"ElectrondeltaEtaSuperClusterTrackAtVtx[NElectrons]/D");
  OutputTree->Branch("ElectrondeltaPhiSuperClusterTrackAtVtx",&ElectrondeltaEtaSuperClusterTrackAtVtx,"ElectrondeltaEtaSuperClusterTrackAtVtx[NElectrons]/D");
  OutputTree->Branch("ElectronsigmaIetaIeta",&ElectronsigmaIetaIeta,"ElectronsigmaIetaIeta[NElectrons]/D");
  OutputTree->Branch("ElectronhadronicOverEm",&ElectronhadronicOverEm,"ElectronhadronicOverEm[NElectrons]/D");
  OutputTree->Branch("ElectronFbrem",&ElectronFbrem,"ElectronFbrem[NElectrons]/D");
  OutputTree->Branch("ElectronEoverP",&ElectronEoverP,"ElectronEoverP[NElectrons]/D");
  OutputTree->Branch("ElectronecalEnergy",&ElectronecalEnergy,"ElectronecalEnergy[NElectrons]/D");
  OutputTree->Branch("Electronp_in",&Electronp_in,"Electronp_in[NElectrons]/D");
  OutputTree->Branch("ElectronDxy",&ElectronDxy,"ElectronDxy[NElectrons]/D");
  OutputTree->Branch("ElectronDz",&ElectronDz,"ElectronDz[NElectrons]/D");
  OutputTree->Branch("ElectronnumberOfLostHits",&ElectronnumberOfLostHits,"ElectronnumberOfLostHits[NElectrons]/D");
  OutputTree->Branch("ElectronsuperClustereta",&ElectronsuperClustereta,"ElectronsuperClustereta[NElectrons]/D");
  OutputTree->Branch("ElectronPFIsolation",&ElectronPFIsolation,"ElectronPFIsolation[NElectrons]/D");
  OutputTree->Branch("ElectronIsMedium", &ElectronIsMedium, "ElectronIsMedium[NElectrons]/D");
  OutputTree->Branch("ElectronIsTrigger", &ElectronIsTrigger, "ElectronIsTrigger[NElectrons]/D");
  OutputTree->Branch("NJets", &NJets, "NJets/I");
  OutputTree->Branch("JetPx", &JetPx, "JetPx[NJets]/D");
  OutputTree->Branch("JetPy", &JetPy, "JetPy[NJets]/D");
  OutputTree->Branch("JetPz", &JetPz, "JetPz[NJets]/D");  
  OutputTree->Branch("JetPt", &JetPt, "JetPt[NJets]/D");
  OutputTree->Branch("JetEta", &JetEta, "JetEta[NJets]/D");
  OutputTree->Branch("JetPhi", &JetPhi, "JetPhi[NJets]/D");
  OutputTree->Branch("JetEnergy", &JetEnergy, "JetEnergy[NJets]/D");
  OutputTree->Branch("JetEmEnergy", &JetEmEnergy, "JetEmEnergy[NJets]/D");
  OutputTree->Branch("JetHadronEnergy", &JetHadronEnergy, "JetHadronEnergy[NJets]/D");
  OutputTree->Branch("NMuons", &NMuons, "NMuons/I");
  OutputTree->Branch("MuonPt", &MuonPt, "MuonPt[NMuons]/D");  
  OutputTree->Branch("MuonEta", &MuonEta, "MuonEta[NMuons]/D");
  OutputTree->Branch("MuonPhi", &MuonPhi, "MuonPhi[NMuons]/D");
  OutputTree->Branch("MuonEnergy", &MuonEnergy, "MuonEnergy[NMuons]/D");
  OutputTree->Branch("METPt", &METPt, "METPt/D");
  OutputTree->Branch("METPhi", &METPhi, "METPhi/D");
  OutputTree->Branch("NMC", &NMC, "NMC/I");
  OutputTree->Branch("MCPt", &MCPt, "MCPt[NMC]/D");
  OutputTree->Branch("MCId", &MCId, "MCId[NMC]/I");
  OutputTree->Branch("MCEta", &MCEta, "MCEta[NMC]/D");
  OutputTree->Branch("MCPhi", &MCPhi, "MCPhi[NMC]/D");
  OutputTree->Branch("MCEnergy", &MCEnergy, "MCEnergy[NMC]/D");
  OutputTree->Branch("MCMass", &MCMass, "MCMass[NMC]/D");
  OutputTree->Branch("MCCharge", &MCCharge, "MCCharge[NMC]/D");
  OutputTree->Branch("MCIndex", &MCIndex, "MCIndex[NMC]/I");
  OutputTree->Branch("MCStatus", &MCStatus, "MCStatus[NMC]/I");
  OutputTree->Branch("NMCmother", &NMCmother, "NMCmother[NMC]/I");
  OutputTree->Branch("MCMotherIndex",&MCMotherIndex,"MCMotherIndex[NMC]/I");
  OutputTree->Branch("NMCdaughter", &NMCdaughter, "NMCdaughter[NMC]/I");
  //  OutputTree->Branch("TowerVector",&TowerVector,"TowerVector[NTowers]/D");
}


ElectronTree::~ElectronTree()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ElectronTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  ClearVariables();

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);
  Handle<CaloTowerCollection> ct;
  iEvent.getByLabel("towerMaker",ct); 
  Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel("gsfElectrons",electrons);
  Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel("ak5PFJets",jets);
  Handle<reco::MuonCollection> muon;
  iEvent.getByLabel("muons",muon);
  Handle<reco::METCollection> mets;
  iEvent.getByLabel("htMetAK5",mets);
  Handle<double> rho_;
  iEvent.getByLabel("fixedGridRhoFastjetAll",rho_);
  double rhoValue = rho_.isValid() ? (double)*(rho_.product()) : 0.0;
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("offlinePrimaryVertices",vertices);
  //     LogInfo("Demo") << "number of tracks "<<tracks->size();
  demohisto->Fill(ct->size());
  if (vertices->empty()) return;
  const reco::Vertex &PV = vertices->front();

  for(const reco::Muon &mu : *muon){
    Lepton m;
    m.v.SetPxPyPzE(mu.px(),mu.py(),mu.pz(),mu.energy());
    m.Charge = mu.charge();
    if(m.v.Pt()>10 && fabs(m.v.Eta())<2.4) //contains both loose or tight electron
      Muons.push_back(m);
  }
  for(reco::GsfElectron const& ele : *electrons){
    Lepton e;
    e.v.SetPxPyPzE(ele.px(),ele.py(),ele.pz(),ele.energy());
    e.Charge = ele.charge();
    if(e.v.Pt()>10 && fabs(e.v.Eta())<2.4) //contains both loose or tight electron                                                               
      Electrons.push_back(e);
  }
  
  Sort(1);
  Sort(2);
  if((int)Muons.size()>1){
    double mass=(Muons.at(0).v + Muons.at(1).v).M();
    if(mass>81 && mass<101) isZfromMu=true;
  }

  if((int)Electrons.size()>1){
    double mass=(Electrons.at(0).v + Electrons.at(1).v).M();
    if(mass>81 && mass<101) isZfromE=true;
  }

  //Filling the tree
  //  if(isZfromMu && (int)electrons->size()>0){
  // if(isZfromE){
  if((int)electrons->size()>0){
    for(CaloTower const& cal : *ct){
      TowerEta[NTowers]=cal.eta();
      TowerPhi[NTowers]=cal.phi();
      TowerEmEt[NTowers]=cal.emEt();
      TowerEmEnergy[NTowers]=cal.emEnergy();
      TowerHadEt[NTowers]=cal.hadEt();
      TowerHadEnergy[NTowers]=cal.hadEnergy();
      TowerEnergy[NTowers]=cal.energy();
    //    TowerVector[NTowers]=cal.p4();
      if(NTowers>4999) cout<<"Change array size of towers"<<endl;
      NTowers++;
    }
    for(reco::GsfElectron const& ele : *electrons){
      ElectronPt[NElectrons]=ele.pt();
      ElectronEta[NElectrons]=ele.eta();
      ElectronPhi[NElectrons]=ele.phi();
      ElectronTrackIso[NElectrons] = ele.dr03TkSumPt();
      ElectronEcalIso[NElectrons]=ele.dr03EcalRecHitSumEt();
      ElectronHcalIso[NElectrons]=ele.dr03HcalTowerSumEt();
      ElectrondeltaEtaSuperClusterTrackAtVtx[NElectrons]=ele.deltaEtaSuperClusterTrackAtVtx();
      ElectrondeltaPhiSuperClusterTrackAtVtx[NElectrons]=ele.deltaPhiSuperClusterTrackAtVtx();
      ElectronsigmaIetaIeta[NElectrons]=ele.sigmaIetaIeta();
      ElectronhadronicOverEm[NElectrons]=ele.hadronicOverEm();
      ElectronFbrem[NElectrons]=ele.fbrem();
      ElectronEoverP[NElectrons]=ele.eSuperClusterOverP();
      ElectronecalEnergy[NElectrons]=ele.ecalEnergy();
      Electronp_in[NElectrons]=ele.trackMomentumAtVtx().R();
      ElectronDxy[NElectrons]=ele.gsfTrack()->dxy(PV.position());
      ElectronDz[NElectrons]=ele.gsfTrack()->dz(PV.position());
      ElectronnumberOfLostHits[NElectrons]=ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      ElectronsuperClustereta[NElectrons]=ele.superCluster()->eta();
      ElectronPFIsolation[NElectrons]=electronPFIsolation(ele,rhoValue);
      ElectronIsMedium[NElectrons]=electronMediumCuts(ele,PV);
      ElectronIsTrigger[NElectrons]=electronTriggerCuts(ele);
      if(NElectrons>999) cout<<"Change array size of electrons"<<endl;
      NElectrons++;
    }
    for(reco::PFJet const& j : *jets){
      if(j.pt()>20){
	JetPx[NJets] = j.px();
	JetPy[NJets] = j.py();
	JetPz[NJets] = j.pz();
	JetPt[NJets] = j.pt();
	JetEta[NJets] = j.eta();
	JetPhi[NJets] = j.phi();
	JetEnergy[NJets] = j.energy();
	JetEmEnergy[NJets] = j.chargedEmEnergy()+j.neutralEmEnergy();
	JetHadronEnergy[NJets] = j.chargedHadronEnergy()+ j.neutralHadronEnergy();
	if(NJets>999) cout<<"Change array size of electrons"<<endl;
	NJets++;
      }
    }
    for(reco::Muon const& mu : *muon){
      MuonPt[NMuons] = mu.pt(); 
      MuonEta[NMuons] = mu.eta();                                            
      MuonPhi[NMuons] = mu.phi();                                            
      MuonEnergy[NMuons] = mu.energy();
      if(NMuons>999) cout<<"more than 1000 muons, please change array size"<<endl;
      NMuons++;

    }
    reco::MET const& met=mets->front();
    METPt = met.pt();
    METPhi = met.phi();
    if( minTracks_ <= tracks->size() ) {
      LogInfo("Demo") << "number of tracks "<<tracks->size();
    }
    if(isMC){
      Handle<reco::GenParticleCollection> genp;
      iEvent.getByLabel("genParticles",genp);
      for(reco::GenParticle const& genparticle : *genp){
	if(genparticle.pt()>1 && genparticle.status()==1){
	  MCId[NMC] = genparticle.pdgId();
	  MCPt[NMC] = genparticle.pt();
	  MCEta[NMC] = genparticle.eta();
	  MCPhi[NMC] = genparticle.phi();
	  MCEnergy[NMC] = genparticle.energy();
	  MCMass[NMC] = genparticle.mass();
	  MCCharge[NMC] = genparticle.charge();
	  MCIndex[NMC] = NMC;
	  MCStatus[NMC] = genparticle.status();
	  NMCmother[NMC] = genparticle.numberOfMothers();
	  MCMotherIndex[NMC]=54321;
	  NMCdaughter[NMC] = genparticle.numberOfDaughters();
	  if(genparticle.numberOfMothers() > 0){
	    FillMothers(&genparticle);
	  }
	  if(NMC>4999) cout<<"more than 5000 truth particles, please change array size"<<endl;    
	  NMC++;
	}
      }
    }
    OutputTree->Fill();
  }
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElectronTree::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronTree::endJob() {
}
void ElectronTree::ClearVariables(){
  NTowers=0;
  isZfromMu=false;
  isZfromE=false;
  double d1=54321;
  Muons.clear();
  Electrons.clear();
  for(int i=0;i<2000;i++){
    TowerEta[i]=54321;
    TowerPhi[i]=54321;
    TowerEmEt[i]=54321;
    TowerEmEnergy[i]=54321;
    TowerEnergy[i]=54321;
    TowerHadEt[i]=54321;
    TowerHadEnergy[i]=54321;
  }
  NElectrons=0;
  for(int i=0;i<1000;i++){
    ElectronEta[i]=54321;
    ElectronPhi[i]=54321;
    ElectronPt[i]=54321;
    ElectronTrackIso[i] = 54321;
    ElectronEcalIso[i]=54321;
    ElectronHcalIso[i]=54321;
    ElectrondeltaEtaSuperClusterTrackAtVtx[i]=d1;
    ElectrondeltaPhiSuperClusterTrackAtVtx[i]=d1;
    ElectronsigmaIetaIeta[i]=d1;
    ElectronhadronicOverEm[i]=d1;
    ElectronFbrem[i]=d1;
    ElectronEoverP[i]=d1;
    ElectronecalEnergy[i]=d1;
    Electronp_in[i]=d1;
    ElectronDxy[i]=d1;
    ElectronDz[i]=d1;
    ElectronnumberOfLostHits[i]=d1;
    ElectronsuperClustereta[i]=d1;
    ElectronPFIsolation[i]=d1;
    ElectronIsMedium[i]=false;
    ElectronIsTrigger[i]=false;
  }
  NJets = 0;
  for(int i =0;i<1000;i++){
    JetPx[i]=d1;
    JetPy[i]=d1;
    JetPz[i]=d1;
    JetPt[i]=d1;
    JetEta[i]=d1;
    JetPhi[i]=d1;
    JetEnergy[i]= d1;
    JetEmEnergy[i] =d1;
    JetHadronEnergy[i] =d1;
  }
  NMuons = 0;
  for(int i =0;i<1000;i++){
    MuonPt[i]=d1;
    MuonEta[i]=d1;
    MuonPhi[i]=d1;
    MuonEnergy[i]= d1;
  }
  NMC=0;
  for(int i=0; i<5000; i++){
    MCPt[i]= d1;
    MCEta[i]= d1;
    MCPhi[i]= d1;
    MCEnergy[i]= d1;
    MCMass[i]= d1;
    MCCharge[i]= d1;
    MCId[i]= 0;
    MCIndex[i] = d1;
    MCStatus[i] = 0;
    NMCmother[i]=d1;
    MCMotherIndex[i]=d1;
    NMCdaughter[i]=d1;
  }

  METPt=d1;
  METPhi=d1;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronTree);
void ElectronTree::Sort(int opt)
{
  if(opt==1){
    for(int i=0; i<(int)Muons.size()-1; i++){
      for(int j=i+1; j<(int)Muons.size(); j++){
	if(Muons.at(i).v.Pt()<Muons.at(j).v.Pt())
	  std::swap(Muons.at(i),Muons.at(j));
      }
    }
  }
  if(opt==2){
    for(int i=0;i<(int)Electrons.size()-1;i++){
      for(int j=i+1; j<(int)Electrons.size(); j++){
        if(Electrons.at(i).v.Pt()<Electrons.at(j).v.Pt())
	  std::swap(Electrons.at(i),Electrons.at(j));
      }
    }
  }
}
double ElectronTree::electronEA(reco::GsfElectron const& it){
  double eta= fabs(it.superCluster()->eta());
  if(eta<1.0000) return 0.13;
  else if(eta<1.4790) return 0.14;
  else if(eta<2.0000) return 0.07;
  else if(eta<2.2000) return 0.09;
  else if(eta<2.3000) return 0.11;
  else if(eta<2.4000) return 0.11;
  else return 0.14; 
}
double ElectronTree::electronPFIsolation(reco::GsfElectron const& it,double rhoValue_){
  reco::GsfElectron::PflowIsolationVariables pfIso = it.pfIsolationVariables();
  double isoChargedHadrons = pfIso.chargedHadronIso;

  double isoNeutralHadrons = pfIso.neutralHadronIso;
  double isoPhotons        = pfIso.photonIso;
  double isoNeutrals       = isoNeutralHadrons + isoPhotons;
  double EA  = electronEA(it);
  double iso = isoChargedHadrons + TMath::Max( 0. , isoNeutrals - rhoValue_*EA );
  return iso;
}

bool ElectronTree::electronMediumCuts(reco::GsfElectron const& it, const reco::Vertex PV_){
  bool isMediumBarrel = ( fabs(it.superCluster()->position().eta()) <= 1.479     &&
			  fabs(it.deltaEtaSuperClusterTrackAtVtx()) < 0.004 &&
			  fabs(it.deltaPhiSuperClusterTrackAtVtx()) < 0.06 &&
			  it.sigmaIetaIeta() < 0.01 &&
			  it.hadronicOverEm() < 0.12 &&
			  fabs(it.gsfTrack()->dxy(PV_.position())) < 0.02 &&
			  fabs(it.gsfTrack()->dz(PV_.position())) < 0.1 &&
			  fabs( 1.0/it.ecalEnergy() - it.eSuperClusterOverP()/it.ecalEnergy() ) < 0.05 &&
			  it.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 1 );
  bool isMediumEndcap = ( fabs(it.superCluster()->position().eta()) > 1.479     &&
			  fabs(it.superCluster()->position().eta()) <  2.500   &&
                          fabs(it.deltaEtaSuperClusterTrackAtVtx()) < 0.007 &&
                          fabs(it.deltaPhiSuperClusterTrackAtVtx()) < 0.03 &&
                          it.sigmaIetaIeta() < 0.03 &&
                          it.hadronicOverEm() < 0.1 &&
                          fabs(it.gsfTrack()->dxy(PV_.position())) < 0.02 &&
                          fabs(it.gsfTrack()->dz(PV_.position())) < 0.1 &&
                          fabs( 1.0/it.ecalEnergy() - it.eSuperClusterOverP()/it.ecalEnergy() ) < 0.05 &&
                          it.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 1 );
  return isMediumBarrel || isMediumEndcap;

}

bool ElectronTree::electronTriggerCuts(reco::GsfElectron const& it){
  bool isTriggerBarrel = ( fabs(it.superCluster()->position().eta()) <= 1.479     &&
			   fabs(it.deltaEtaSuperClusterTrackAtVtx()) < 0.007 &&
			   fabs(it.deltaPhiSuperClusterTrackAtVtx()) < 0.15 &&
			   it.sigmaIetaIeta() < 0.01 &&
			   it.hadronicOverEm() < 0.12 &&
			   it.dr03EcalRecHitSumEt()/it.pt() < 0.2 &&
			   it.dr03HcalTowerSumEt()/it.pt() < 0.2 &&
			   it.dr03TkSumPt()/it.pt()  < 0.2 );
  bool isTriggerEndcap = ( fabs(it.superCluster()->position().eta()) > 1.479     &&
			   fabs(it.superCluster()->position().eta()) <  2.500   &&
			   fabs(it.deltaEtaSuperClusterTrackAtVtx()) < 0.009 &&
			   fabs(it.deltaPhiSuperClusterTrackAtVtx()) < 0.10 &&
			   it.sigmaIetaIeta() < 0.03 &&
			   it.hadronicOverEm() < 0.1 &&
			   it.dr03EcalRecHitSumEt()/it.pt() < 0.2 &&
                           it.dr03HcalTowerSumEt()/it.pt() < 0.2 &&
                           it.dr03TkSumPt()/it.pt()  < 0.2 );

  return isTriggerBarrel || isTriggerEndcap;

}
const reco::Candidate* ElectronTree::getmother(const reco::Candidate* particle)//Finds mother of the particle
{
  const reco::Candidate* mo;
  mo = particle;
  if(particle->numberOfMothers()==0) return mo ;
  if(particle->mother(0)->pdgId() !=particle->pdgId()){
    mo = particle->mother(0);
    return mo;
  }
  if(particle->mother(0)->pdgId() == particle->pdgId()){
    mo = getmother(particle->mother(0));
  }
  return mo;
}
void ElectronTree::FillMothers(const reco::Candidate* particle){ //Fill properties of mother
  if(particle->numberOfMothers()==0) return ;  
  
  for(int m = 1;m<NMC;m++){
    if(MCId[m]==getmother(particle)->pdgId()&&MCId[m]==2212&&(TMath::Abs(MCEta[m]-getmother(particle)->eta())<0.001)&&(TMath::Abs(MCPhi[m]-getmother(particle)->phi())<0.001)){ //)
      MCMotherIndex[NMC] = m;
      return;
    }
    if(MCId[m]==getmother(particle)->pdgId()&&(TMath::Abs(MCPt[m]-getmother(particle)->pt())/getmother(particle)->pt()<0.00001)&&(TMath::Abs(MCEta[m]-getmother(particle)->eta())<0.001)&&(TMath::Abs(MCPhi[m]-getmother(particle)->phi())<0.001)){
      MCMotherIndex[NMC] = m;
      //cout<<"Found already stored mother of "<<MCId[NMC]<<" with PID "<<MCId[m]<<" at position "<<m<<" pt diff is "<<(TMath::Abs(MCPt[m]-getmother(particle)->pt())/getmother(particle)->pt())<<endl;
      return;
    }
  }

  MCMotherIndex[NMC] = NMC+1;
  NMC++;
  MCId[NMC] = getmother(particle)->pdgId();
  MCPt[NMC] = getmother(particle)->pt();
  MCEta[NMC] = getmother(particle)->eta();
  MCPhi[NMC] = getmother(particle)->phi();
  MCEnergy[NMC] = getmother(particle)->energy();
  MCMass[NMC] = getmother(particle)->mass();
  MCCharge[NMC] = getmother(particle)->charge();
  MCStatus[NMC] = getmother(particle)->status();
  NMCmother[NMC] = getmother(particle)->numberOfMothers();
  NMCdaughter[NMC] = getmother(particle)->numberOfDaughters();
  MCIndex[NMC] = NMC;
  MCMotherIndex[NMC] = 54321;
  if(getmother(particle)->numberOfMothers()>0){
    particle = getmother(particle);
    FillMothers(particle);
  }

}







 

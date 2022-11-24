#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <bitset>
#include <algorithm>

unsigned const NMAX(256);

struct LHEEvent {
  unsigned iProc;
  unsigned nP{};
  int id[NMAX]{};
  int status[NMAX]{};
  double px[NMAX]{};
  double py[NMAX]{};
  double pz[NMAX]{};
  double e[NMAX]{};
  double m[NMAX]{};
};

class LHAupLHEFCustom : public Pythia8::LHAupLHEF {
public:
  LHAupLHEFCustom(Pythia8::Info* info, const char* filename) :
    LHAupLHEF(info, filename)
  {}

  bool useExternal() override { return true; }
};

class Channel {
public:
  Channel(std::string const& configFile, unsigned order, unsigned nPartons, std::string const& lheFile);
  ~Channel();

  unsigned order() const { return order_; }
  unsigned nPartons() const { return nPartons_; }
  void setParams(double meanAbsW, double wnorm);
  void printParams() const;
  double weight() const;
  void init();
  void setUNLOPSOpt(bool counterTerm = false);
  HepMC::GenEvent const* getEvent(bool counterTerm, bool& atEnd);
  LHEEvent getLHEEvent() const;
  unsigned numSkipped() const { return numSkipped_; }

private:
  Pythia8::Pythia pythia_;
  Pythia8::LHAupLHEF* lheInput_;
  HepMC::Pythia8ToHepMC toHepMC_;
  HepMC::GenEvent genEvent_;
  std::string lheFile_;
  unsigned order_;
  unsigned nPartons_;
  double meanAbsW_;
  double wnorm_;
  unsigned numSkipped_;
  unsigned nEvents_;
};

Channel::Channel(std::string const& _configFile, unsigned _order, unsigned _nPartons, std::string const& _lheFile) :
  pythia_("", false),
  lheInput_(nullptr),
  lheFile_(_lheFile),
  order_(_order),
  nPartons_(_nPartons),
  meanAbsW_(0.),
  wnorm_(0.),
  numSkipped_(0),
  nEvents_(0)
{
  pythia_.settings.mode("Beams:frameType", 4);
  pythia_.settings.mode("Merging:nRequested", _nPartons);

  // Switch off warnings for parton-level events.
  toHepMC_.set_print_inconsistency(false);
  toHepMC_.set_free_parton_exception(false);
  // Do not store cross section information, as this will be done manually.
  toHepMC_.set_store_pdf(false);
  toHepMC_.set_store_proc(false);
  toHepMC_.set_store_xsec(false);

  std::ifstream configStream(_configFile);

  std::string line;
  bool found(false);
  while (std::getline(configStream, line)) {
    if (line.find("[PYTHIA]") != std::string::npos) {
      pythia_.readFile(configStream);
      found = true;
      break;
    }
  }

  configStream.close();

  if (!found)
    pythia_.readFile(_configFile);
}

Channel::~Channel()
{
  pythia_.setLHAupPtr(nullptr);
  delete lheInput_;
}

void
Channel::setParams(double _meanAbsW, double _wnorm)
{
  meanAbsW_ = _meanAbsW;
  wnorm_ = _wnorm;
}

void
Channel::printParams() const
{
  std::cout.precision(10);
  std::cout << "Channel mean|W|: " << meanAbsW_ << std::endl;
  std::cout << "Channel wnorm: " << wnorm_ << std::endl;
}

double
Channel::weight() const
{
  if (nPartons_ == 0) // no counter term -> reduce occurrence by half
    return 0.5 * meanAbsW_;
  else
    return meanAbsW_;
}

void
Channel::init()
{
  delete lheInput_;

  lheInput_ = new LHAupLHEFCustom(&pythia_.info, lheFile_.c_str());

  pythia_.setLHAupPtr(lheInput_);

  // need to have at least one merging-related flag to be on to turn on merging
  setUNLOPSOpt();

  pythia_.init();

  if (std::abs(lheInput_->strategy()) == 4)
    setParams(lheInput_->xMax(0), lheInput_->xMax(0));
  else
    setParams(lheInput_->xSec(0), 1.);

  numSkipped_ = 0;
}

void
Channel::setUNLOPSOpt(bool _counterTerm/* = false*/)
{
  pythia_.settings.flag("Merging:doUNLOPSTree", false);
  pythia_.settings.flag("Merging:doUNLOPSLoop", false);
  pythia_.settings.flag("Merging:doUNLOPSSubt", false);
  pythia_.settings.flag("Merging:doUNLOPSSubtNLO", false);

  if (_counterTerm) {
    switch (order_) {
    case 0:
      pythia_.settings.flag("Merging:doUNLOPSSubt", true);
      break;
    case 1:
      pythia_.settings.flag("Merging:doUNLOPSSubtNLO", true);
      break;
    default:
      break;
    }
    pythia_.settings.mode("Merging:nRecluster", 1);
  }
  else {
    switch (order_) {
    case 0:
      pythia_.settings.flag("Merging:doUNLOPSTree", true);
      break;
    case 1:
      pythia_.settings.flag("Merging:doUNLOPSLoop", true);
      break;
    default:
      break;
    }
    pythia_.settings.mode("Merging:nRecluster", 0);
  }
}

HepMC::GenEvent const*
Channel::getEvent(bool _counterTerm, bool& _atEnd)
{
  setUNLOPSOpt(_counterTerm);

  bool success(false);
  try {
    success = pythia_.next();
  }
  catch (std::invalid_argument& ex) {
    std::cout << lheFile_ << " " << nEvents_ << std::endl;
    std::cout << "iProc " << lheInput_->idProcess() << std::endl;
    std::cout << "nP " << lheInput_->sizePart() << std::endl;
    for (unsigned iP(0); iP != lheInput_->sizePart(); ++iP) {
      std::cout << "id[" << iP << "]" << lheInput_->id(iP) << std::endl;
      std::cout << "status[" << iP << "]" << lheInput_->status(iP) << std::endl;
      std::cout << "px[" << iP << "]" << lheInput_->px(iP) << std::endl;
      std::cout << "py[" << iP << "]" << lheInput_->py(iP) << std::endl;
      std::cout << "pz[" << iP << "]" << lheInput_->pz(iP) << std::endl;
      std::cout << "e[" << iP << "]" << lheInput_->e(iP) << std::endl;
      std::cout << "m[" << iP << "]" << lheInput_->m(iP) << std::endl;
    }
    throw;
  }
    
  if (!success) {
    if (pythia_.info.atEndOfFile())
      _atEnd = true;

    return nullptr;
  }

  double ckkwWeight(pythia_.info.mergingWeightNLO());
  if (!(ckkwWeight == ckkwWeight)) { // NaN
    ++numSkipped_;
    return nullptr;
  }

  double weight(ckkwWeight);

  if (std::abs(pythia_.info.lhaStrategy()) == 4) {
    // we need to make all events unweighted
    weight *= pythia_.info.weight() / wnorm_;
  }
  else {
    weight *= pythia_.info.weight();
  }

  if (_counterTerm)
    weight *= -1.;

  genEvent_.clear();

  // Set hepmc event weight.
  genEvent_.weights().push_back(weight);

  // Fill HepMC event.
  toHepMC_.fill_next_event(pythia_, &genEvent_);

  return &genEvent_;
}

LHEEvent
Channel::getLHEEvent() const
{
  LHEEvent outEvent{};

  outEvent.iProc = lheInput_->idProcess();
  outEvent.nP = lheInput_->sizePart();
  for (unsigned iP(0); iP != outEvent.nP; ++iP) {
    outEvent.id[iP] = lheInput_->id(iP);
    outEvent.status[iP] = lheInput_->status(iP);
    outEvent.px[iP] = lheInput_->px(iP);
    outEvent.py[iP] = lheInput_->py(iP);
    outEvent.pz[iP] = lheInput_->pz(iP);
    outEvent.e[iP] = lheInput_->e(iP);
    outEvent.m[iP] = lheInput_->m(iP);
  }

  return outEvent;
}

int
main(int argc, char* argv[])
{
  std::string configFile(argv[1]);
  std::ifstream configStream(configFile);
  std::string line;

  std::vector<Channel*> channels;

  while (std::getline(configStream, line)) {
    if (line.find("[CHANNELS]") != std::string::npos)
      break;
  }

  while (std::getline(configStream, line)) {
    if (line.size() == 0)
      continue;

    if (line.find("[PYTHIA]") != std::string::npos)
      break;

    std::stringstream sstream;
    sstream.str(line);

    unsigned order;
    unsigned npart;
    std::string lhe;
    sstream >> order >> npart >> lhe;

    std::cout << line << std::endl;

    channels.push_back(new Channel(configFile, order, npart, lhe));
  }

  configStream.close();

  double wtotal(0.);
  for (auto* ch : channels) {
    ch->init();
    wtotal += ch->weight();
  }

  for (auto* ch : channels)
    ch->printParams();

  TFile* outputFile(TFile::Open(argv[2], "recreate"));
  auto* output(new TTree("events", "Events"));

  double weight;
  unsigned short channel;
  bool counterTerm;
  unsigned short nPart;
  int pPid[NMAX];
  float pPt[NMAX];
  float pEta[NMAX];
  float pPhi[NMAX];
  float pM[NMAX];
  bool hasBoson;
  float vPid;
  float vPt;
  float vEta;
  float vPhi;
  float vM;
  int lPid[2];
  float lPt[2];
  float lEta[2];
  float lPhi[2];
  unsigned short nNeutrino;
  int nPid[NMAX];
  float nPt[NMAX];
  float nEta[NMAX];
  float nPhi[NMAX];
  unsigned short nJ;
  float jPt[NMAX];
  float jEta[NMAX];
  float jPhi[NMAX];
  float jM[NMAX];

  output->Branch("weight", &weight, "weight/D");
  output->Branch("channel", &channel, "channel/s");
  output->Branch("counterTerm", &counterTerm, "counterTerm/O");
  output->Branch("parton.size", &nPart, "size/s");
  output->Branch("parton.pid", pPid, "pid[parton.size]/I");
  output->Branch("parton.pt", pPt, "pt[parton.size]/F");
  output->Branch("parton.eta", pEta, "eta[parton.size]/F");
  output->Branch("parton.phi", pPhi, "phi[parton.size]/F");
  output->Branch("parton.m", pM, "m[parton.size]/F");
  output->Branch("hasBoson", &hasBoson, "hasBoson/O");
  output->Branch("boson.pt", &vPt, "pt/F");
  output->Branch("boson.eta", &vEta, "eta/F");
  output->Branch("boson.phi", &vPhi, "phi/F");
  output->Branch("boson.m", &vM, "m/F");
  output->Branch("lepton.pid", lPid, "pid[2]/I");
  output->Branch("lepton.pt", lPt, "pt[2]/F");
  output->Branch("lepton.eta", lEta, "eta[2]/F");
  output->Branch("lepton.phi", lPhi, "phi[2]/F");
  output->Branch("neutrino.size", &nNeutrino, "size/s");
  output->Branch("neutrino.pid", nPid, "pid[neutrino.size]/I");
  output->Branch("neutrino.pt", nPt, "pt[neutrino.size]/F");
  output->Branch("neutrino.eta", nEta, "eta[neutrino.size]/F");
  output->Branch("neutrino.phi", nPhi, "phi[neutrino.size]/F");
  output->Branch("jet.size", &nJ, "size/s");
  output->Branch("jet.pt", jPt, "pt[jet.size]/F");
  output->Branch("jet.eta", jEta, "eta[jet.size]/F");
  output->Branch("jet.phi", jPhi, "phi[jet.size]/F");
  output->Branch("jet.m", jM, "m[jet.size]/F");

  auto chEnd(channels.end());

  auto* jetDef(new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));

  std::vector<fastjet::PseudoJet> fjinputs;
  fjinputs.reserve(4096);

  std::function<HepMC::GenParticle const*(HepMC::GenParticle const&)> findBoson;
  findBoson = [&findBoson](HepMC::GenParticle const& _part)->HepMC::GenParticle const* {
    auto* vtx(_part.end_vertex());
    if (vtx == nullptr)
      return nullptr;

    for (auto mItr(vtx->particles_out_const_begin()); mItr != vtx->particles_out_const_end(); ++mItr) {
      if ((*mItr)->pdg_id() == _part.pdg_id())
        return findBoson(**mItr);
    }

    return &_part;
  };

  std::function<HepMC::GenParticle const&(HepMC::GenParticle const&)> findAncestor;
  findAncestor = [&findAncestor](HepMC::GenParticle const& _part)->HepMC::GenParticle const& {
    auto* vtx(_part.production_vertex());
    assert(vtx != nullptr);

    auto mItr(vtx->particles_in_const_begin());
    auto mEnd(vtx->particles_in_const_end());
    for (; mItr != mEnd; ++mItr) {
      if ((*mItr)->pdg_id() == _part.pdg_id()) {
        return findAncestor(**mItr);
      }
    }

    return _part;
  };

  std::function<HepMC::GenParticle const&(HepMC::GenParticle const&)> findDecendant;
  findDecendant = [&findDecendant](HepMC::GenParticle const& _part)->HepMC::GenParticle const& {
    auto* vtx(_part.end_vertex());
    if (vtx == nullptr)
      return _part;

    auto mItr(vtx->particles_out_const_begin());
    auto mEnd(vtx->particles_out_const_end());
    for (; mItr != mEnd; ++mItr) {
      if ((*mItr)->pdg_id() == _part.pdg_id())
        return findDecendant(**mItr);
    }

    throw std::runtime_error("lepton disappeared");
  };

  std::function<bool(HepMC::GenParticle const&, HepMC::GenParticle const&)> matchParent;
  matchParent = [&matchParent](HepMC::GenParticle const& _part, HepMC::GenParticle const& _boson)->bool {
    auto* vtx(_part.production_vertex());
    if (vtx == nullptr)
      return false;

    auto mItr(vtx->particles_in_const_begin());
    auto mEnd(vtx->particles_in_const_end());
    for (; mItr != mEnd; ++mItr) {
      if (*mItr == &_boson)
        break;
      else if ((*mItr)->pdg_id() == _part.pdg_id())
        return matchParent(**mItr, _boson);
    }

    return mItr != mEnd;
  };

  long long nEvents(-1);
  if (argc > 3)
    nEvents = std::atoi(argv[3]);

  TRandom3 rand;
  rand.SetSeed();

  long long iEvent(0);
  bool atEnd(false);

  while (iEvent++ != nEvents) {
    Channel* ch(nullptr);
    HepMC::GenEvent const* evt(nullptr);

    while (true) {
      double r(rand.Uniform(wtotal));
  
      channel = 0;
      double t(0.);
      while (channel != channels.size()) {
        ch = channels[channel];
        t += ch->weight();
        if (r < t)
          break;
  
        ++channel;
      }
  
      if (ch->nPartons() == 0)
        counterTerm = false;
      else
        counterTerm = rand.Uniform() < 0.5;
  
      evt = ch->getEvent(counterTerm, atEnd);
      if (evt != nullptr || atEnd)
        break;
    }

    if (atEnd)
      break;

    continue;

    weight = evt->weights()[0];

    auto lheEvent(ch->getLHEEvent());

    nPart = 0;
    for (unsigned iP(0); iP != lheEvent.nP; ++iP) {
      if (lheEvent.status[iP] != 1)
        continue;

      pPid[nPart] = lheEvent.id[iP];

      double px(lheEvent.px[iP]);
      double py(lheEvent.py[iP]);
      double pz(lheEvent.pz[iP]);
      double pt(std::sqrt(px * px + py * py));
      double p(std::sqrt(pt * pt + pz * pz));

      pPt[nPart] = pt;
      pEta[nPart] = 0.5 * std::log((p + pz) / (p - pz));
      pPhi[nPart] = std::atan2(py, px);
      pM[nPart] = lheEvent.m[iP];

      ++nPart;
    }

    typedef std::pair<HepMC::GenParticle const*, HepMC::GenParticle const*> LeptonPair;

    std::vector<LeptonPair> leptonPairs;

    for (auto&& pItr(evt->particles_begin()); pItr != evt->particles_end(); ++pItr) {
      auto& part(**pItr);

      if (part.status() != 1)
        continue;

      int pdgId(part.pdg_id());

      if (pdgId != -11 && pdgId != -13)
        continue;

      auto& leptonAncestor(findAncestor(part));
      auto* vtx(leptonAncestor.production_vertex());
      assert(vtx != nullptr);

      for (auto dItr(vtx->particles_out_const_begin()); dItr != vtx->particles_out_const_end(); ++dItr) {
        if ((*dItr)->pdg_id() == -pdgId) {
          auto& antilepton(findDecendant(**dItr));
          leptonPairs.emplace_back(&part, &antilepton);
          break;
        }
      }
    }

    hasBoson = !leptonPairs.empty();

    LeptonPair const* bosonDecay(nullptr);

    if (hasBoson) {
      double mMax(0.);
      for (auto& leptonPair : leptonPairs) {
        auto&& lP(leptonPair.first->momentum());
        auto&& aP(leptonPair.second->momentum());

        TLorentzVector vP(lP.x() + aP.x(), lP.y() + aP.y(), lP.z() + aP.z(), lP.e() + aP.e());
        double m(vP.M());

        if (m > mMax) {
          mMax = m;
          vPt = vP.Pt();
          vEta = vP.Eta();
          vPhi = TVector2::Phi_mpi_pi(vP.Phi());
          vM = m;

          lPid[0] = leptonPair.first->pdg_id();
          lPt[0] = lP.perp();
          lEta[0] = lP.eta();
          lPhi[0] = TVector2::Phi_mpi_pi(lP.phi());

          lPid[1] = leptonPair.second->pdg_id();
          lPt[1] = aP.perp();
          lEta[1] = aP.eta();
          lPhi[1] = TVector2::Phi_mpi_pi(aP.phi());

          bosonDecay = &leptonPair;
        }
      }
    }
    else {
      vPt = 0.;
      vEta = 0.;
      vPhi = 0.;
      vM = 0.;

      std::fill_n(lPid, 2, 0);
      std::fill_n(lPt, 2, 0.);
      std::fill_n(lEta, 2, 0.);
      std::fill_n(lPhi, 2, 0.);
    }

    fjinputs.clear();

    nNeutrino = 0;

    for (auto&& pItr(evt->particles_begin()); pItr != evt->particles_end(); ++pItr) {
      auto& part(**pItr);
    
      if (part.status() != 1)
        continue;

      auto&& p4(part.momentum());      
      unsigned absId(std::abs(part.pdg_id()));

      if (absId == 12 || absId == 14 || absId == 16) {
        nPid[nNeutrino] = part.pdg_id();
        nPt[nNeutrino] = p4.perp();
        nEta[nNeutrino] = p4.eta();
        nPhi[nNeutrino] = TVector2::Phi_mpi_pi(p4.phi());
        ++nNeutrino;
        continue;
      }

      if (hasBoson && (absId == 11 || absId == 13) && (&part == bosonDecay->first || &part == bosonDecay->second))
        continue;

      fjinputs.emplace_back(p4.px(), p4.py(), p4.pz(), p4.e());
    }

    fastjet::ClusterSequence seq(fjinputs, *jetDef);
    auto jets(fastjet::sorted_by_pt(seq.inclusive_jets(15.)));

    nJ = jets.size();
    for (unsigned iJ(0); iJ != nJ; ++iJ) {
      fastjet::PseudoJet const& jet(jets[iJ]);

      jPt[iJ] = jet.pt();
      jEta[iJ]= jet.eta();
      jPhi[iJ] = TVector2::Phi_mpi_pi(jet.phi());
      jM[iJ] = jet.m();
    }

    output->Fill();
  }

  outputFile->cd();
  output->Write();
  delete outputFile;

  for (unsigned iC(0); iC != channels.size(); ++iC) {
    std::cout << "Channel " << iC << ": " << channels[iC]->numSkipped() << " events with not-a-number weight" << std::endl;
    delete channels[iC];
  }

  delete jetDef;

  return 0;
}

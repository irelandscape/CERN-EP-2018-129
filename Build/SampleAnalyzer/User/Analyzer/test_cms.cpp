#include "SampleAnalyzer/User/Analyzer/test_cms.h"
using namespace MA5;
using namespace std;

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
  bool test_cms::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;

  Manager()->AddRegionSelection("diphoton");
  Manager()->AddRegionSelection("diptau");
  Manager()->AddHisto("low_pTmiss", 20, 0, 1000);
  Manager()->AddHisto("high_pTmiss", 20, 0, 1000);
  Manager()->AddHisto("eth", 20, 0, 1000);
  Manager()->AddHisto("uth", 20, 0, 1000);
  Manager()->AddHisto("thth", 20, 0, 1000);

  Manager()->AddCut("leading_photons_pT");
  Manager()->AddCut("diphoton_M");
  Manager()->AddCut("azimuthal_sep");
  Manager()->AddCut("jets_azimutal_sep");
  Manager()->AddCut("multijet_backgrounds");
  Manager()->AddCut("kinetic_requirements");
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void test_cms::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  cout << "END   Finalization" << endl;
}

void test_cms::addIfSignalPhoton (const EventFormat& event,
                                  vector<const RecPhotonFormat*>& signal_photons,
                                  const RecPhotonFormat photon)
{
  /* Insert photon in a vector sorted by pT.
   * pT must be at least 20GeV
   */

  auto eta = photon.abseta();
  auto photon_pt = photon.pt();

  if (eta > 1.44 && (eta < 1.57 || eta > 2.5))
    return;

  // Must be at least 20 GeV
  if (photon_pt < 20)
    return;

  double iPi = PHYSICS->Isol->eflow->sumIsolation(photon,
                                                  event.rec(),
                                                  0.3,
                                                  0.,
                                                  IsolationEFlow::TRACK_COMPONENT);
  double In = PHYSICS->Isol->eflow->sumIsolation(photon,
                                                 event.rec(),
                                                 0.3,
                                                 0.,
                                                 IsolationEFlow::NEUTRAL_COMPONENT);
  double Igam = PHYSICS->Isol->eflow->sumIsolation(photon,
                                                   event.rec(),
                                                   0.3,
                                                   0.,
                                                   IsolationEFlow::PHOTON_COMPONENT);

  bool signal_photon = false;

  if (eta < 1.44)
  {
    // barrel region
    if ((Igam - photon_pt) < (0.7 + (0.005 * photon_pt)) && 
        In < (1.0  + (0.04 * photon_pt)) && 
        iPi < 1.5)
    {
      signal_photon = true;
    }
  }
  else
  {
    // end-cap region
    if ((Igam - photon_pt) < (1.0 + (0.005 * photon_pt)) && 
        In < (1.5  + (0.04 * photon_pt)) && 
        iPi < 1.2)
    {
      signal_photon = true;
    }
  }

  if (!signal_photon)
    return;

  // This is a signal photon. Insert in vector ordered by decreasing pT
  bool inserted = false;
  for (auto it = signal_photons.begin();
       it != signal_photons.end();
       it++)
  {
    if ((*it)->pt() < photon_pt)
    {
      signal_photons.insert(it, &photon);
      inserted = true;
      break;
    }
  }

  if (!inserted)
    signal_photons.push_back(&photon);
}

void test_cms::addIfSignalLepton(const EventFormat& event,
                                 vector<const RecLeptonFormat*>& signal_leptons,
                                 const RecLeptonFormat& lepton,
                                 double min_pt,
                                 double min_eta,
                                 double DR)
{
  auto pt = lepton.pt();
  auto eta = lepton.abseta();

  if (pt < min_pt)
    return;

  if (eta >= min_eta)
    return;

  // Based on other code found in madanalysis code
  double all_E = PHYSICS->Isol->eflow->sumIsolation(lepton,
                                                    event.rec(),
                                                    DR,
                                                    0.,
                                                    IsolationEFlow::ALL_COMPONENTS);

  if (all_E < 0.2 * pt)
    signal_leptons.push_back(&lepton);
}

const RecTauFormat* test_cms::getPairTau (const EventFormat& event,
                                          const RecTauFormat& tau)
{
  const RecTauFormat* best_tau = NULL;

  for (auto &tau2 : event.rec()->taus())
  {
    if (&tau2 == &tau)
      continue;

    if (tau2.pt() < 40 || (best_tau && tau2.pt() < best_tau->pt()))
      continue;

    if (tau.dr(tau2) >= 0.5)
      continue;

    if (tau.abseta() >= 2.1)
      continue;

    best_tau = &tau2;
  }

  return best_tau;
}

bool test_cms::vetoSignalPhotons(MALorentzVector& pTmiss,
                                  vector<const RecPhotonFormat*>& signal_photons,
                                  double& M, // Return the diphoton invariant mass
                                  double& pT) // Return the diphoton transverse momentum
{
  // Must have at least 2 signal photons with enough momentum
  if(!Manager()->ApplyCut(signal_photons.size() >= 2 &&
                          signal_photons[0]->pt() >= 30 &&
                          signal_photons[1]->pt() >= 20,
                          "leading_photons_pT"))
  {
    return true;
  }

  // Check diphoton invariant mass is < 65GeV
  MALorentzVector diPhoton;
  diPhoton.SetPxPyPzE(signal_photons[0]->px()+signal_photons[1]->px(),
                      signal_photons[0]->py()+signal_photons[1]->py(),
                      signal_photons[0]->pz()+signal_photons[1]->pz(),
                      signal_photons[0]->e()+signal_photons[1]->e());

  M = diPhoton.M();
  if(!Manager()->ApplyCut(M < 95, "diphoton_M"))
  {

    return true;
  }

  pT = diPhoton.Pt();

  // diphoton azimuthal separation requirement
  if (!Manager()->ApplyCut(fabs(diPhoton.DeltaPhi(pTmiss)) > 2.1,
                           "azimuthal_sep"))
  {
    return true;
  }

  return false;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool test_cms::Execute(SampleFormat& sample, const EventFormat& event)
{
  vector<const RecLeptonFormat*> electrons;
  vector<const RecPhotonFormat*> signal_photons;
  vector<const RecLeptonFormat*> signal_electrons;
  vector<const RecLeptonFormat*> signal_muons;
  // Each vector index correspond to one of the selected tauh pair.
  vector<const RecLeptonFormat*> vector_tauh1;
  vector<const RecLeptonFormat*> vector_tauh2;
  MALorentzVector pTmiss = event.rec()->MET().momentum();
  double pTmiss_amplitude = pTmiss.Pt();
  double diphoton_M;
  double diphoton_pt;


  if (event.rec() == NULL || event.mc() == NULL)
      return true;
  //
  // Required for histograms
  Manager()->InitializeForNewEvent(event.mc()->weight());

  // Find leading photons
  for (auto &photon : event.rec()->photons())
  {
    this->addIfSignalPhoton(event, signal_photons, photon);
  }

  // Signal electrons
  for (auto &electron : event.rec()->electrons())
  {
    this->addIfSignalLepton(event, signal_electrons, electron, 26, 2.1, 0.3);
  }
  
  // Signal muons
  for (auto &muon : event.rec()->muons())
  {
    this->addIfSignalLepton(event, signal_electrons, muon, 26, 2.4, 0.3);
  }
  
  if (pTmiss_amplitude > 105)
  {
    // etauh
    for (auto &electron : signal_electrons)
    {
      for (auto& tau : event.rec()->taus())
      {
        if (electron->dr(tau) > 0.5 || tau.pt() <= 20 || tau.abseta() >= 2.3)
          continue;

        Manager()->FillHisto("eth", pTmiss_amplitude);
      }
    }
    
    // utauh
    for (auto &muon : signal_muons)
    {
      for (auto& tau : event.rec()->taus())
      {
        if (muon->dr(tau) > 0.5 || tau.pt() <= 20 || tau.abseta() >= 2.3)
          continue;

        Manager()->FillHisto("uth", pTmiss_amplitude);
      }
    }

    // tauhtauh
    for (auto& tau : event.rec()->taus())
    {
      const RecTauFormat* tau2;
      MALorentzVector diTauE;

      if ((tau2 = getPairTau(event, tau)) != NULL)
      {
        if ((tau.pt() + tau2->pt()) < 65)
          continue;

        diTauE.SetPxPyPzE(tau.px()+tau2->px(),
                          tau.py()+tau2->py(),
                          tau.pz()+tau2->pz(),
                          tau.e()+tau2->e());

        if (diTauE.M() >= 125)
          continue;

        Manager()->FillHisto("thth", pTmiss_amplitude);
      }
    }
  }

  // continue diphoton channel
  if (this->vetoSignalPhotons(pTmiss, 
                               signal_photons, 
                               diphoton_M,
                               diphoton_pt))
  {
    return true;
  }

  // Check azimuthal separation between pTmiss and jets with pT > 50Gev
  vector<const RecJetFormat*> signal_jets;
  for (auto &jet : event.rec()->jets())
  {
    if (jet.pt() > 50)
      signal_jets.push_back(&jet);
  }

  unsigned jet_count = 0;
  for (auto jet: signal_jets)
  {
    if (!Manager()->ApplyCut(jet->dphi_0_pi(pTmiss) > 0.5,
                             "jets_azimutal_sep"))
    {
      return true;
    }

    if (jet->pt() >= 30)
      ++jet_count;
  }

  if (!Manager()->ApplyCut(jet_count < 3, "multijet_backgrounds"))
  {
    return true;
  }

  // Kinematic requirements
  if (pTmiss_amplitude > 50 && pTmiss_amplitude < 130 &&
      (signal_photons[0]->pt() / diphoton_M) > 0.45 &&
      (signal_photons[1]->pt() / diphoton_M) > 0.25 &&
      diphoton_pt > 75)
  {
    Manager()->FillHisto("low_pTmiss", pTmiss_amplitude);
  }
  else if (pTmiss_amplitude > 130 &&
           (signal_photons[0]->pt() / diphoton_M) > 0.5 &&
           (signal_photons[1]->pt() / diphoton_M) > 0.25 &&
           diphoton_pt > 90)
  {
    Manager()->FillHisto("high_pTmiss", pTmiss_amplitude);
  }
  else
  {
    Manager()->ApplyCut(false, "kinetic_requirements");
  }

  return true;
}


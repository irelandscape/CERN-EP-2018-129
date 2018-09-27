#ifndef analysis_test_cms_h
#define analysis_test_cms_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class test_cms : public AnalyzerBase
{
  INIT_ANALYSIS(test_cms,"test_cms")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
  void addIfSignalPhoton (const EventFormat& event,
                          std::vector<const RecPhotonFormat*>& signal_photons,
                          const RecPhotonFormat photon);
  bool vetoSignalPhotons(MALorentzVector& pTmiss,
                         std::vector<const RecPhotonFormat*>& signal_photons,
                         double& M,
                         double& pT);
  void addIfSignalLepton(const EventFormat& event,
                         std::vector<const RecLeptonFormat*>& signal_leptons,
                         const RecLeptonFormat& lepton,
                         double min_pt,
                         double min_eta,
                         double DR);
  const RecTauFormat* getPairTau (const EventFormat& event,
                                  const RecTauFormat& tau);
                                        
};
}

#endif

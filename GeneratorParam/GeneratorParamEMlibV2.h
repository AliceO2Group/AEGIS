#ifndef GENERATORPARAMEMLIBV2_H
#define GENERATORPARAMEMLIBV2_H
#include "TRandom.h"
#include "TObject.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "GeneratorParamLibBase.h"

class iostream;
class TRandom;
class TF1;

using namespace std;

class GeneratorParamEMlibV2 : public GeneratorParamLibBase {
  
public:
  
  enum Particle_t{kPizero=0, kEta, kRho0, kOmega, kEtaprime, kPhi, kJpsi, kPsi2S, kUpsilon,
    kSigma0, kK0s, kDeltaPlPl, kDeltaPl, kDeltaMi, kDeltaZero,
    kRhoPl, kRhoMi, kK0star, kK0l, kLambda, kKPl, kKMi,
    kOmegaPl, kOmegaMi, kXiPl, kXiMi, kSigmaPl, kSigmaMi,
    kDirectRealGamma, kDirectVirtGamma, kNParticles};

   static const int kNHadrons = kNParticles-2; // total number of particles minus DirectRealGamma and DirectVirtGamma

  enum CollisionSystem_t {kpp900GeV=0x000, kpp2760GeV=0x64, kpp7TeV=0xC8, kpPb=0x12C, kPbPb=0x190};
  
  enum Centrality_t{ kpp = 0x0, k0005=0x1, k0510=0x2, k1020=0x3, k2030=0x4, k3040=0x5, k4050=0x6, k5060=0x7,
    k0010=0x8, k2040=0x9, k4060=0xA, k6080=0xB, k0020=0xC, k0040=0xD, k2080=0xE, k4080=0xF, k2050=0x10, kCentralities=0x11};
  
  enum v2Sys_t{kLoV2Sys=-1, kNoV2Sys=0, kUpV2Sys=+1};

  GeneratorParamEMlibV2() { };

  static void SelectParams( Int_t collisionSystem,
                            Int_t centSelect      = kpp,
                            Int_t v2sys           = kNoV2Sys)
  {
    fgSelectedCollisionsSystem  = collisionSystem;
    fgSelectedCentrality        = centSelect;
    fgSelectedV2Systematic      = v2sys;
  }
  
  GenFunc   GetPt(Int_t param, const char * tname=0) const;
  GenFunc   GetY(Int_t param, const char * tname=0) const;
  GenFuncIp GetIp(Int_t param, const char * tname=0) const;
  GenFunc   GetV2(Int_t param, const char * tname=0) const;
  
  // General functions
  static Bool_t SetPtParametrizations(TString fileName, TString dirName);
  static Bool_t SetFlowParametrizations(TString fileName, TString dirName) ;
  static void   SetMtScalingFactors(TString fileName, TString dirName);
  static Bool_t SetPtYDistributions(TString fileName, TString dirName);
  static TF1*   GetPtParametrization(Int_t np);
  static TH1D*  GetMtScalingFactors();
  static TH2F*  GetPtYDistribution(Int_t np);

  static Int_t fgSelectedCollisionsSystem;                                                      // selected pT parameter
  static Int_t fgSelectedCentrality;                                                            // selected Centrality
  static Int_t fgSelectedV2Systematic;                                                          // selected v2 systematics, usefully values: -1,0,1

  static Double_t PtExponential(const Double_t *pt, const Double_t *param);
  static Double_t PtModifiedHagedornPowerlaw(const Double_t *pt, const Double_t *param);
  static Double_t IntegratedKrollWada(const Double_t *mh, const Double_t *);
  
  static Double_t YFlat(Double_t y);
  static TF1*     MtScal(Int_t np, TString name, Bool_t isMeson);
  static Double_t V2Param(const Double_t *px, const Double_t *param);
  static Double_t EtScalingV2(Double_t ptPart, Int_t iPart, Int_t iRef) ;
  static Double_t V2Flat(const Double_t *px, const Double_t *param);
  static Double_t KEtScal(Double_t pt, Int_t np, Int_t nq=2);
  static Double_t GetTAA(Int_t cent);
  
  static Double_t CrossOverLc(double a, double b, double x);
  static Double_t CrossOverRc(double a, double b, double x);

  static const Double_t fgkV2param[kCentralities][16];                     // parameters of pi v2
  static const Double_t fgkRawPtOfV2Param[kCentralities][10];              // parameters of the raw pt spectrum of v2 analysys
  static const Double_t fgkThermPtParam[kCentralities][2];                 // parameters of thermal gamma pt
  static const Double_t fgkHM[kNHadrons];                                         // particle masses
  static const Double_t fgkMtFactor[3][kNHadrons];                                // mt scaling factor

  // direct gamma
  static Double_t PtPromptRealGamma(const Double_t *px, const Double_t *dummy);
  static Double_t PtPromptVirtGamma(const Double_t *px, const Double_t *dummy);
  static Double_t PtThermalRealGamma(const Double_t *px, const Double_t *dummy);
  static Double_t PtThermalVirtGamma(const Double_t *px, const Double_t *dummy);
  
  static Int_t    IpDirectRealGamma(TRandom *ran);
  static Double_t PtDirectRealGamma(const Double_t *px, const Double_t *dummy);
  static Double_t YDirectRealGamma(const Double_t *py, const Double_t *dummy);
  static Double_t V2DirectRealGamma(const Double_t *px, const Double_t *dummy);
  
  static Int_t    IpDirectVirtGamma(TRandom *ran);
  static Double_t PtDirectVirtGamma(const Double_t *px, const Double_t *dummy);
  static Double_t YDirectVirtGamma(const Double_t *py, const Double_t *dummy);
  static Double_t V2DirectVirtGamma(const Double_t *px, const Double_t *dummy);
  
  // Pizero
  static Int_t    IpPizero(TRandom *ran);
  static Double_t PtPizero(const Double_t *px, const Double_t *dummy);
  static Double_t YPizero(const Double_t *py, const Double_t *dummy);
  static Double_t V2Pizero(const Double_t *px, const Double_t *dummy);
  
  // Eta
  static Int_t    IpEta(TRandom *ran);
  static Double_t PtEta(const Double_t *px, const Double_t *dummy);
  static Double_t YEta(const Double_t *py, const Double_t *dummy);
  static Double_t V2Eta(const Double_t *px, const Double_t *dummy);
  
  // Rho
  static Int_t    IpRho0(TRandom *ran);
  static Double_t PtRho0(const Double_t *px, const Double_t *dummy);
  static Double_t YRho0(const Double_t *py, const Double_t *dummy);
  static Double_t V2Rho0(const Double_t *py, const Double_t *dummy);
  
  // Omega
  static Int_t    IpOmega(TRandom *ran);
  static Double_t PtOmega(const Double_t *px, const Double_t *dummy);
  static Double_t YOmega(const Double_t *py, const Double_t *dummy);
  static Double_t V2Omega(const Double_t *py, const Double_t *dummy);
  
  // Etaprime
  static Int_t    IpEtaprime(TRandom *ran);
  static Double_t PtEtaprime(const Double_t *px, const Double_t *dummy);
  static Double_t YEtaprime(const Double_t *py, const Double_t *dummy);
  static Double_t V2Etaprime(const Double_t *py, const Double_t *dummy);
  
  // Phi
  static Int_t    IpPhi(TRandom *ran);
  static Double_t PtPhi(const Double_t *px, const Double_t *dummy);
  static Double_t YPhi(const Double_t *py, const Double_t *dummy);
  static Double_t V2Phi(const Double_t *py, const Double_t *dummy);
  
  // Jpsi
  static Int_t    IpJpsi(TRandom *ran);
  static Double_t PtJpsi(const Double_t *px, const Double_t *dummy);
  static Double_t YJpsi(const Double_t *py, const Double_t *dummy);
  static Double_t V2Jpsi(const Double_t *py, const Double_t *dummy);

  // psi(2S)
  static Int_t    IpPsi2S(TRandom *ran);
  static Double_t PtPsi2S(const Double_t *px, const Double_t *dummy);
  static Double_t YPsi2S(const Double_t *py, const Double_t *dummy);
  static Double_t V2Psi2S(const Double_t *py, const Double_t *dummy);

  // Upsilon
  static Int_t    IpUpsilon(TRandom *ran);
  static Double_t PtUpsilon(const Double_t *px, const Double_t *dummy);
  static Double_t YUpsilon(const Double_t *py, const Double_t *dummy);
  static Double_t V2Upsilon(const Double_t *py, const Double_t *dummy);
  
  // Sigma
  static Int_t    IpSigma0(TRandom *ran);
  static Double_t PtSigma0(const Double_t *px, const Double_t *dummy);
  static Double_t YSigma0(const Double_t *py, const Double_t *dummy);
  static Double_t V2Sigma0(const Double_t *px, const Double_t *dummy);
  
  // K0short
  static Int_t    IpK0short(TRandom *ran);
  static Double_t PtK0short(const Double_t *px, const Double_t *dummy);
  static Double_t YK0short(const Double_t *py, const Double_t *dummy);
  static Double_t V2K0short(const Double_t *px, const Double_t *dummy);

  // K0long
  static Int_t    IpK0long(TRandom *ran);
  static Double_t PtK0long(const Double_t *px, const Double_t *dummy);
  static Double_t YK0long(const Double_t *py, const Double_t *dummy);
  static Double_t V2K0long(const Double_t *px, const Double_t *dummy);

  // Lambda
  static Int_t    IpLambda(TRandom *ran);
  static Double_t PtLambda(const Double_t *px, const Double_t *dummy);
  static Double_t YLambda(const Double_t *py, const Double_t *dummy);
  static Double_t V2Lambda(const Double_t *px, const Double_t *dummy);
  
  // Delta++
  static Int_t    IpDeltaPlPl(TRandom *ran);
  static Double_t PtDeltaPlPl(const Double_t *px, const Double_t *dummy);
  static Double_t YDeltaPlPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaPlPl(const Double_t *px, const Double_t *dummy);
  
  // Delta+
  static Int_t    IpDeltaPl(TRandom *ran);
  static Double_t PtDeltaPl(const Double_t *px, const Double_t *dummy);
  static Double_t YDeltaPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaPl(const Double_t *px, const Double_t *dummy);
  
  // Delta-
  static Int_t    IpDeltaMi(TRandom *ran);
  static Double_t PtDeltaMi(const Double_t *px, const Double_t *dummy);
  static Double_t YDeltaMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaMi(const Double_t *px, const Double_t *dummy);
  
  // Delta0
  static Int_t    IpDeltaZero(TRandom *ran);
  static Double_t PtDeltaZero(const Double_t *px, const Double_t *dummy);
  static Double_t YDeltaZero(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaZero(const Double_t *px, const Double_t *dummy);
  
  // Rho+
  static Int_t    IpRhoPl(TRandom *ran);
  static Double_t PtRhoPl(const Double_t *px, const Double_t *dummy);
  static Double_t YRhoPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2RhoPl(const Double_t *px, const Double_t *dummy);
  
  // Rho-
  static Int_t    IpRhoMi(TRandom *ran);
  static Double_t PtRhoMi(const Double_t *px, const Double_t *dummy);
  static Double_t YRhoMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2RhoMi(const Double_t *px, const Double_t *dummy);
  
  // K0*
  static Int_t    IpK0star(TRandom *ran);
  static Double_t PtK0star(const Double_t *px, const Double_t *dummy);
  static Double_t YK0star(const Double_t *py, const Double_t *dummy);
  static Double_t V2K0star(const Double_t *px, const Double_t *dummy);

  // K+
  static Int_t    IpKPl(TRandom *ran);
  static Double_t PtKPl(const Double_t *px, const Double_t *dummy);
  static Double_t YKPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2KPl(const Double_t *px, const Double_t *dummy);

  // K-
  static Int_t    IpKMi(TRandom *ran);
  static Double_t PtKMi(const Double_t *px, const Double_t *dummy);
  static Double_t YKMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2KMi(const Double_t *px, const Double_t *dummy);

  // Omega+
  static Int_t    IpOmegaPl(TRandom *ran);
  static Double_t PtOmegaPl(const Double_t *px, const Double_t *dummy);
  static Double_t YOmegaPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2OmegaPl(const Double_t *px, const Double_t *dummy);

  // Omega-
  static Int_t    IpOmegaMi(TRandom *ran);
  static Double_t PtOmegaMi(const Double_t *px, const Double_t *dummy);
  static Double_t YOmegaMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2OmegaMi(const Double_t *px, const Double_t *dummy);

  // Xi+
  static Int_t    IpXiPl(TRandom *ran);
  static Double_t PtXiPl(const Double_t *px, const Double_t *dummy);
  static Double_t YXiPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2XiPl(const Double_t *px, const Double_t *dummy);

  // Xi-
  static Int_t    IpXiMi(TRandom *ran);
  static Double_t PtXiMi(const Double_t *px, const Double_t *dummy);
  static Double_t YXiMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2XiMi(const Double_t *px, const Double_t *dummy);

  // Sigma+
  static Int_t    IpSigmaPl(TRandom *ran);
  static Double_t PtSigmaPl(const Double_t *px, const Double_t *dummy);
  static Double_t YSigmaPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2SigmaPl(const Double_t *px, const Double_t *dummy);

  // Sigma-
  static Int_t    IpSigmaMi(TRandom *ran);
  static Double_t PtSigmaMi(const Double_t *px, const Double_t *dummy);
  static Double_t YSigmaMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2SigmaMi(const Double_t *px, const Double_t *dummy);

private:
  static TF1*     fPtParametrization[kNHadrons];     // pt paramtrizations
  static TF1*     fPtParametrizationProton;   // pt paramtrization
  static TH1D*    fMtFactorHisto;             // mt scaling factors
  static TH2F*    fPtYDistribution[kNHadrons];       // pt-y distributions
  static TF1*     fV2Parametrization[kNHadrons+1];     // pt paramtrizations
  static Int_t    fV2RefParameterization[kNHadrons+1]; // ID of a hadron used for parameterization of V2 for Et scaling

  ClassDef(GeneratorParamEMlibV2,7);
};

#endif

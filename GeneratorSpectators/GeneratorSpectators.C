#include "FairGenerator.h"
#include "GeneratorSpectators.h"
#include "Generators/GeneratorTGenerator.h"
#include <fairlogger/Logger.h>

class GeneratorSpectatorsO2 : public o2::eventgen::GeneratorTGenerator {

public:
  GeneratorSpectatorsO2() : GeneratorTGenerator("GeneratorSpectators", "GeneratorSpectators") {
    spec = new GeneratorSpectators();
    setTGenerator(spec);
  };

  ~GeneratorSpectatorsO2() { delete spec; };

  Bool_t Init() override {
    spec->Init();
    return true;
  }

  void updateHeader(o2::dataformats::MCEventHeader *eventHeader) final {
    // not updating the header, just printing particles for the time being
    LOG(info) << Form("--- Number of particles: %zu ---", mParticles.size());
    for (int i = 0; i < mParticles.size(); ++i) {
      auto &p = mParticles[i];
      LOG(info) << Form(
          "PdgCode: %d, FirstMother: %d, FirstDaughter: %d, LastDaughter: %d",
          p.GetPdgCode(), p.GetFirstMother(), p.GetFirstDaughter(),
          p.GetLastDaughter());
    }
  }

  GeneratorSpectators* getGenerator() const { return spec; }

private:
  GeneratorSpectators *spec = nullptr;
};

FairGenerator *GeneratorSingleNeutron() {
  auto wrap = new GeneratorSpectatorsO2();
  auto spec = wrap->getGenerator();
  spec->SetNpart(1);
  spec->SetMomentum(5360. / 2.);
  spec->SetDirection(0, 0., 0., 1.);
  spec->SetDivergence(0.000032);  // beam divergence in murad
  return wrap;
}

FairGenerator *GeneratorNeutrons(int nNeutrons = -1, float b = -1., bool useFluctuation = false) {
  auto wrap = new GeneratorSpectatorsO2();
  auto spec = wrap->getGenerator();
  spec->SetImpactParameter(b);
  spec->SetNpartFluctuation(useFluctuation);
  spec->SetNpart(nNeutrons);
  spec->SetMomentum(5360. / 2.);
  spec->SetDirection(0, 0., 0., 1.);
  spec->SetDivergence(0.000032);  // beam divergence in murad
  return wrap;
}


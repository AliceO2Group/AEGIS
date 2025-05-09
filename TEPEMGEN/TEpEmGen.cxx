/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 *                                                                        *
 * Copyright(c) 1997, 1998, 2002, Adrian Alscher and Kai Hencken          *
 * See $ALICE_ROOT/EpEmGen/diffcross.f for full Copyright notice          *
 *                                                                        *
 *                                                                        *
 * Copyright(c) 2002 Kai Hencken, Yuri Kharlov, Serguei Sadovsky          *
 * See $ALICE_ROOT/EpEmGen/epemgen.f for full Copyright notice            *
 *                                                                        *
 **************************************************************************/
//------------------------------------------------------------------------
// TEpEmGen is an interface class to fortran event generator of
// single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 GeV/c
//%
// Yuri.Kharlov@cern.ch
// 9 October 2002
//
// Revised on September 2018 for ALICEo2: Roberto Preghenella (preghenella@bo.infn.it)
//------------------------------------------------------------------------

#include "TRandom.h"

#include "TEpEmGen.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TEcommon.h"

#ifndef WIN32
# define ee_init  ee_init_
# define ee_event ee_event_
# define eernd    eernd_
#else
# define ee_init  EE_INIT
# define ee_event EE_EVENT
# define eernd    EERND
#endif

extern "C" {
  void ee_init  (Double_t &ymin, Double_t &ymax, Double_t &xmin, Double_t &xmax, Double_t &cm_energy, Double_t &Z);
  void ee_event (Double_t &ymin, Double_t &ymax, Double_t &xmin, Double_t &xmax,
	         Double_t &yE,   Double_t &yP,   Double_t &xE,   Double_t &xP, 
		 Double_t &phi,  Double_t &w);
  Double_t eernd(Int_t*) {
    Double_t r;
    do r=gRandom->Rndm(); while(0 >= r || r >= 1);
    return r;
  }
}

ClassImp(TEpEmGen)

//------------------------------------------------------------------------------
TEpEmGen::TEpEmGen() : TGenerator("TEpEmGen","TEpEmGen")
{
// TEpEmGen constructor: creates a TClonesArray in which it will store all
// particles. Note that there may be only one functional TEpEmGen object
// at a time, so it's not use to create more than one instance of it.

}

//------------------------------------------------------------------------------
TEpEmGen::~TEpEmGen()
{
  // Destroys the object, deletes and disposes all TParticles currently on list.
  fParticles->Delete();
  delete fParticles;
  fParticles = nullptr; // to prevent double deletion with TGenerator destructor
}

//______________________________________________________________________________
void TEpEmGen::GenerateEvent(Double_t ymin, Double_t ymax, Double_t ptmin, Double_t ptmax,
			     Double_t &yElectron, Double_t &yPositron,
			     Double_t &xElectron, Double_t &xPositron,
			     Double_t &phi12,     Double_t &weight)
{
  //produce one event
  ee_event(ymin,ymax,ptmin,ptmax,
	   yElectron,yPositron,xElectron,xPositron,
	   phi12,weight);
}

//______________________________________________________________________________
void TEpEmGen::Initialize(Double_t ymin, Double_t ymax, Double_t ptmin, Double_t ptmax, Double_t cm_energy, Double_t Z)
{
  // Initialize EpEmGen
  Double_t ptminMeV = ptmin*1000;
  Double_t ptmaxMeV = ptmax*1000;
  ee_init(ymin,ymax,ptminMeV,ptmaxMeV,cm_energy,Z);
}

//______________________________________________________________________________
Int_t TEpEmGen::ImportParticles(TClonesArray *particles, Option_t *option)
{
  if (particles == 0) return 0;
  TClonesArray &clonesParticles = *particles;
  clonesParticles.Clear();
  Int_t numpart = fParticles->GetEntries();
  for (Int_t i = 0; i<numpart; i++) {
    TParticle *particle = (TParticle *)fParticles->At(i);
    new(clonesParticles[i]) TParticle(*particle);
  }
  return numpart;
}


//______________________________________________________________________________
Double_t TEpEmGen::GetXsection()
{
  // Return cross section accumulated so far
  return EEVENT.Xsecttot;
}

//______________________________________________________________________________
Double_t TEpEmGen::GetDsection()
{
  // Return cross section error accumulated so far
  return EEVENT.Dsecttot;
}

// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP
// Collaboration

#include "DPPythia8Generator.h"

#include <cmath>
#include <set>
#include <vector>

#include "BeamSmearingUtils.h"
#include "FairPrimaryGenerator.h"
#include "Pythia8/Pythia.h"
#include "ShipUnit.h"
#include "TDatabasePDG.h"  // for TDatabasePDG
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"

using ShipUnit::c_light;
using ShipUnit::cm;
using ShipUnit::mm;
const Int_t debug = 1;

// -----   Default constructor   -------------------------------------------
DPPythia8Generator::DPPythia8Generator() {
  // fHadDecay = false;
  fId = 2212;           // proton
  fMom = 400;           // proton
  fDP = 9900015;        // DP  pdg code
  fLmin = 5000. * cm;   // mm minimum  decay position z  ROOT units !
  fLmax = 12000. * cm;  // mm maximum decay position z
  fFDs = 7.7 / 10.4;  // correction for Pythia6 to match measured Ds production
  fpbrem = kFALSE;
  fpbremPDF = 0;
  fsmearBeam = 8 * mm;  // default value for smearing beam (8 mm)
  fPaintBeam = 5 * cm;  // default value for painting beam (5 cm)
  fdy = kFALSE;
  fDPminM = 0.5;
  fInputFile = nullptr;
  fnRetries = 0;
  fnDPtot = 0;
  fShipEventNr = 0;
  fPythia = new Pythia8::Pythia();
  // fPythiaHadDecay =  new Pythia8::Pythia();
}
// -------------------------------------------------------------------------

void DPPythia8Generator::Print() { fPythia->settings.listAll(); };

// -----   Default constructor   -------------------------------------------
Bool_t DPPythia8Generator::Init() {
  if (fUseRandom1) fRandomEngine = std::make_shared<PyTr1Rng>();
  if (fUseRandom3) fRandomEngine = std::make_shared<PyTr3Rng>();
  fPythia->setRndmEnginePtr(fRandomEngine);
  // fPythiaHadDecay->setRndmEnginePtr(fRandomEngine);
  fn = 0;

  if (!fpbrem) {
    if (debug) {
      std::cout << "Beam Momentum " << fMom << std::endl;
    }
    fPythia->settings.mode("Beams:idA", fId);
    fPythia->settings.mode("Beams:idB", 2212);
    fPythia->settings.mode("Beams:frameType", 2);
    fPythia->settings.parm("Beams:eA", fMom);  // codespell:ignore parm
    fPythia->settings.parm("Beams:eB", 0.);    // codespell:ignore parm

    if (fdy)
      fPythia->settings.parm("PhaseSpace:mHatMin",  // codespell:ignore parm
                             fDPminM);

  } else {
    if (!fpbremPDF) {
      LOG(fatal) << "Failed in retrieving dark photon PDF for production by "
                    "proton bremstrahlung!";
      return kFALSE;
    }
  }

  TDatabasePDG* pdgBase = TDatabasePDG::Instance();
  Double_t root_ctau = pdgBase->GetParticle(fDP)->Lifetime();

  if (debug) {
    std::cout << "Final particle parameters for PDGID " << fDP << ":"
              << std::endl;
    List(fDP);
  }
  if (debug) {
    std::cout << "tau root PDG database " << root_ctau
              << "[s] ctau root = " << root_ctau * 3e10 << "[cm]" << std::endl;
  }
  fctau = fPythia->particleData.tau0(fDP);  //* 3.3333e-12
  if (debug) {
    std::cout << "ctau pythia " << fctau << "[mm]" << std::endl;
  }
  int initPass = fPythia->init();
  if (debug) {
    std::cout << "Pythia initialisation bool: " << initPass << std::endl;
  }

  if (!initPass) {
    LOG(fatal) << "Pythia initialisation failed";
    return kFALSE;
  }

  return kTRUE;
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
DPPythia8Generator::~DPPythia8Generator() {}
// -------------------------------------------------------------------------

// -----   Passing the event   ---------------------------------------------
Bool_t DPPythia8Generator::ReadEvent(FairPrimaryGenerator* cpg) {
  auto [dx, dy] = CalculateBeamOffset(fsmearBeam, fPaintBeam);
  Double_t tp, tS, zp, xp, yp, zS, xS, yS, pz, px, py, e, w;
  Double_t tm, zm, xm, ym, pmz, pmx, pmy, em;
  Int_t im;

  Int_t imN, imtmp, idtmp;
  Double_t tN, zN, xN, yN, pzN, pxN, pyN, eN;
  Int_t imout, idpout;

  int iDP =
      0;  // index of the chosen DP, also ensures that at least 1 DP is produced
  std::vector<int>
      dec_chain;  // pythia indices of the particles to be stored on the stack
  std::vector<int> dpvec;  // pythia indices of DP particles

  do {
    // bit for proton brem
    if (fpbrem) {
      fPythia->event.reset();
      double dpmom = 0;
      double thetain = 0;
      fpbremPDF->GetRandom2(dpmom, thetain);
      double dpm = fPythia->particleData.m0(fDP);
      double dpe = sqrt(dpmom * dpmom + dpm * dpm);
      double phiin = 2. * M_PI * gRandom->Rndm();

      if (debug > 1) {
        std::cout << " Adding DP gun with p " << dpmom << " m " << dpm << " e "
                  << dpe << " theta,phi " << thetain << "," << phiin
                  << std::endl
                  << std::flush;
      }
      fPythia->event.append(fDP, 1, 0, 0, dpmom * sin(thetain) * cos(phiin),
                            dpmom * sin(thetain) * sin(phiin),
                            dpmom * cos(thetain), dpe, dpm);
    }

    if (!fPythia->next()) LOG(fatal) << "fPythia->next() failed";

    // fPythia->event.list();
    for (int i = 0; i < fPythia->event.size(); i++) {
      // find all DP
      if (abs(fPythia->event[i].id()) == fDP) {
        dpvec.push_back(i);
      }
    }
    iDP = dpvec.size();
    fnDPtot += iDP;
    if (iDP == 0) {
      // fLogger->Info(MESSAGE_ORIGIN,"Event without DP. Retry.");
      // fPythia->event.list();
      fnRetries +=
          1;  // can happen if phasespace does not allow meson to decay to DP
    } else {
      // for mesons, could have more than one ... but for DY prod, need to take
      // the last one... int r =  int( gRandom->Uniform(0,iDP) );
      int r = iDP - 1;

      int i = dpvec[r];
      // production vertex
      zp = fPythia->event[i].zProd();
      xp = fPythia->event[i].xProd();
      yp = fPythia->event[i].yProd();
      tp = fPythia->event[i].tProd();
      // momentum
      pz = fPythia->event[i].pz();
      px = fPythia->event[i].px();
      py = fPythia->event[i].py();
      e = fPythia->event[i].e();
      // old decay vertex
      if (debug > 1) {
        std::cout << " Debug: decay product of A: "
                  << fPythia->event[fPythia->event[i].daughter1()].id() << " "
                  << fPythia->event[fPythia->event[i].daughter2()].id()
                  << std::endl;
        //  new decay vertex
        Double_t LS = gRandom->Uniform(fLmin / mm, fLmax / mm);  // in mm
        Double_t p = TMath::Sqrt(px * px + py * py + pz * pz);
        Double_t lam = LS / p;
        xS = xp + lam * px;
        yS = yp + lam * py;
        zS = zp + lam * pz;
        Double_t gam = e / TMath::Sqrt(e * e - p * p);
        Double_t beta = p / e;
        tS = tp + LS / beta;  // all in Pythia units
        w = TMath::Exp(-LS / (beta * gam * fctau)) *
            ((fLmax / mm - fLmin / mm) / (beta * gam * fctau));

        // direct mother of the DP
        im = (Int_t)fPythia->event[i].mother1();

        // look upstream in the mother1 chain for the first proton or neutron
        imN = -1;
        imtmp = im;
        while (imtmp > 0) {
          idtmp = TMath::Abs(fPythia->event[imtmp].id());
          if (idtmp == 2212 || idtmp == 2112) {
            imN = imtmp;
            break;
          }
          imtmp = (Int_t)fPythia->event[imtmp].mother1();
        }

        // if found, store proton/neutron ancestor first
        if (imN > 0 && imN != im) {
          zN = fPythia->event[imN].zProd();
          xN = fPythia->event[imN].xProd();
          yN = fPythia->event[imN].yProd();
          pzN = fPythia->event[imN].pz();
          pxN = fPythia->event[imN].px();
          pyN = fPythia->event[imN].py();
          eN = fPythia->event[imN].e();
          tN = fPythia->event[imN].tProd();

          cpg->AddTrack(
              (Int_t)fPythia->event[imN].id(), pxN, pyN, pzN, xN * mm + dx,
              yN * mm + dy, zN * mm, -1, false, eN, tN * mm / c_light,
              w);  // proton/neutron ancestor is the root of the exported chain

          dec_chain.push_back(imN);

          if (debug > 1)
            std::cout << std::endl
                      << " insert nucleon ancestor id " << imN
                      << " pdg=" << fPythia->event[imN].id() << " pz = " << pzN
                      << " [GeV], z = " << zN << " [mm] t = " << tN << " [mm/c]"
                      << std::endl;
        }
        // direct mother of DP
        zm = fPythia->event[im].zProd();
        xm = fPythia->event[im].xProd();
        ym = fPythia->event[im].yProd();
        pmz = fPythia->event[im].pz();
        pmx = fPythia->event[im].px();
        pmy = fPythia->event[im].py();
        em = fPythia->event[im].e();
        tm = fPythia->event[im].tProd();

        // if nucleon ancestor exists, mother points to track 0
        imout = -1;
        if (imN > 0 && imN != im) imout = 0;

        cpg->AddTrack(
            (Int_t)fPythia->event[im].id(), pmx, pmy, pmz, xm * mm + dx,
            ym * mm + dy, zm * mm, imout, false, em, tm * mm / c_light,
            w);  // convert pythia's (x,y,z[mm], t[mm/c]) to ([cm], [s])

        // if nucleon ancestor exists, DP points to track 1, otherwise to track
        // 0
        idpout = 0;
        if (imN > 0 && imN != im) idpout = 1;

        cpg->AddTrack(fDP, px, py, pz, xp * mm + dx, yp * mm + dy, zp * mm,
                      idpout, false, e, tp * mm / c_light, w);

        // bookkeep the indices of stored particles
        dec_chain.push_back(im);
        dec_chain.push_back(i);

        if (debug > 1)
          std::cout << std::endl
                    << " insert mother id " << im
                    << " pdg=" << fPythia->event[im].id() << " pmz = " << pmz
                    << " [GeV],  zm = " << zm << " [mm] tm = " << tm
                    << " [mm/c]" << std::endl;
        if (debug > 1)
          std::cout << " ----> insert DP id " << i << " pdg=" << fDP
                    << " pz = " << pz << " [GeV] zp = " << zp
                    << " [mm] tp = " << tp << " [mm/c]" << std::endl;
        iDP = i;
      }
    }
    while (iDP ==
           0);  // ----------- avoid rare empty events w/o any DP's produced

    if (fShipEventNr % 100 == 0) {
      LOGF(info, "ship event %i / pythia event-nr %i", fShipEventNr, fn);
    }
    fShipEventNr += 1;
    // fill a container with pythia indices of the DP decay chain
    if (debug > 1) std::cout << "Filling daughter particles" << std::endl;
    // if (!hadDecay){
    for (int k = 0; k < fPythia->event.size(); k++) {
      // if daughter of DP, copy
      if (debug > 1)
        std::cout << k << " pdg =" << fPythia->event[k].id() << " mum "
                  << fPythia->event[k].mother1() << std::endl;
      im = fPythia->event[k].mother1();
      while (im > 0) {
        if (im == iDP) {
          break;
        }  // pick the decay products of only 1 chosen DP
        // if ( abs(fPythia->event[im].id())==fDP && im == iDP ){break;}
        else {
          im = fPythia->event[im].mother1();
        }
      }
      if (im < 1) {
        if (debug > 1) std::cout << "reject" << std::endl;
        continue;
      }
      if (debug > 1) std::cout << "accept" << std::endl;
      dec_chain.push_back(k);
    }

    // go over daughters and store them on the stack
    // original code started from +2 because dec_chain contained:
    //   [mother, DP]
    // now it can contain:
    //   [nucleon ancestor, mother, DP]
    // therefore, if the first stored particle is a proton/neutron,
    // daughters start from +3; otherwise they start from +2.
    std::vector<int>::iterator itbeg = dec_chain.begin() + 2;
    if (dec_chain.size() > 2 &&
        (TMath::Abs(fPythia->event[dec_chain[0]].id()) == 2212 ||
         TMath::Abs(fPythia->event[dec_chain[0]].id()) == 2112)) {
      itbeg = dec_chain.begin() + 3;
    }

    for (std::vector<int>::iterator it = itbeg; it != dec_chain.end(); ++it) {
      // pythia index of the particle to store
      int k = *it;
      // find mother position on the output stack: impy -> im
      int impy = fPythia->event[k].mother1();
      std::vector<int>::iterator itm =
          std::find(dec_chain.begin(), dec_chain.end(), impy);
      im = -1;  // safety
      if (itm != dec_chain.end())
        im = itm - dec_chain.begin();  // convert iterator into sequence number

      Bool_t wanttracking = false;
      if (fPythia->event[k].isFinal()) {
        wanttracking = true;
      }
      pz = fPythia->event[k].pz();
      px = fPythia->event[k].px();
      py = fPythia->event[k].py();
      e = fPythia->event[k].e();
      if (fextFile) {
        im += 1;
      };
      cpg->AddTrack((Int_t)fPythia->event[k].id(), px, py, pz, xS * mm, yS * mm,
                    zS * mm, im, wanttracking, e, tS * mm / c_light, w);
    }
    return kTRUE;
  }
  // -------------------------------------------------------------------------
  void DPPythia8Generator::SetParameters(char* par) {
    // Set Parameters
    fPythia->readString(par);
    if (debug) {
      std::cout << "fPythia->readString(\"" << par << "\")" << std::endl;
    }
  }

  // -------------------------------------------------------------------------

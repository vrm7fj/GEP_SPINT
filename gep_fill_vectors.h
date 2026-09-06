#ifndef GEP_FILL_VECTORS_H
#define GEP_FILL_VECTORS_H

#include "TChain.h"
#include "TTreeFormula.h"

#include <vector>
#include <iostream>

#include "gep_config.h"

// ============================================================
// Structure containing quantities needed after event selection
// ============================================================

struct GEPData {

  std::vector<double> ft_x;
  std::vector<double> ft_y;
  std::vector<double> ft_xp;
  std::vector<double> ft_yp;

  std::vector<double> fpp_x;
  std::vector<double> fpp_y;
  std::vector<double> fpp_xp;
  std::vector<double> fpp_yp;

  std::vector<double> t_vz;
  std::vector<double> t_vx;
  std::vector<double> t_vy;
  
  std::vector<double> er_ft;

  std::vector<double> eresidu_ft;
  std::vector<double> eresidv_ft;
  std::vector<double> eresidu_fpp;
  std::vector<double> eresidv_fpp;

  std::vector<double> mod_ft;
  std::vector<double> mod_fpp;
  
};


// ============================================================
// Fill vectors
// ============================================================

void FillVectors(TChain *C, GEPData &data) {

  // ----------------------------------------------------------
  // Branch variables
  // ----------------------------------------------------------

  Double_t gemFT_x[MAXHIT];
  Double_t gemFT_y[MAXHIT];
  Double_t gemFT_xp[MAXHIT];
  Double_t gemFT_yp[MAXHIT];

  Double_t gemFPP_x[MAXHIT];
  Double_t gemFPP_y[MAXHIT];
  Double_t gemFPP_xp[MAXHIT];
  Double_t gemFPP_yp[MAXHIT];

  Double_t vz[MAXHIT];
  Double_t vx[MAXHIT];
  Double_t vy[MAXHIT];

  Double_t eresidu_ft[MAXHIT];
  Double_t eresidv_ft[MAXHIT];
  Double_t eresidu_fpp[MAXHIT];
  Double_t eresidv_fpp[MAXHIT];


  // ----------------------------------------------------------
  // Enable branches
  // ----------------------------------------------------------

  C->SetBranchStatus("*", 0);

  C->SetBranchStatus("sbs.tr.*", 1);
  C->SetBranchStatus("sbs.gemFT.track.*", 1);
  C->SetBranchStatus("sbs.gemFPP.track.*", 1);
  C->SetBranchStatus("sbs.hcal.*", 1);
  C->SetBranchStatus("earm.ecal.*", 1);
  C->SetBranchStatus("heep.*", 1);
  C->SetBranchStatus("scalhel.*", 1);


  // ----------------------------------------------------------
  // Set branch addresses
  // ----------------------------------------------------------

  C->SetBranchAddress("sbs.gemFT.track.x",  gemFT_x);
  C->SetBranchAddress("sbs.gemFT.track.y",  gemFT_y);
  C->SetBranchAddress("sbs.gemFT.track.xp", gemFT_xp);
  C->SetBranchAddress("sbs.gemFT.track.yp", gemFT_yp);
  C->SetBranchAddress("sbs.gemFT.hit.eresidu", eresidu_ft);
  C->SetBranchAddress("sbs.gemFT.hit.eresidv", eresidv_ft);

  C->SetBranchAddress("sbs.gemFPP.track.x",  gemFPP_x);
  C->SetBranchAddress("sbs.gemFPP.track.y",  gemFPP_y);
  C->SetBranchAddress("sbs.gemFPP.track.xp", gemFPP_xp);
  C->SetBranchAddress("sbs.gemFPP.track.yp", gemFPP_yp);
  C->SetBranchAddress("sbs.gemFPP.hit.eresidu", eresidu_fpp);
  C->SetBranchAddress("sbs.gemFPP.hit.eresidv", eresidv_fpp);

  C->SetBranchAddress("sbs.tr.vz", vz);
  C->SetBranchAddress("sbs.tr.vx", vx);
  C->SetBranchAddress("sbs.tr.vy", vy);

  // ----------------------------------------------------------
  // Global cut
  // ----------------------------------------------------------

  TTreeFormula *GlobalCut =
    new TTreeFormula("GlobalCut", globalcut, C);


  // ----------------------------------------------------------
  // Event loop
  // ----------------------------------------------------------

  Long64_t nevent = 0;

  int treenum = -1;
  int oldtreenum = -1;

  while (C->GetEntry(nevent)) {

    treenum = C->GetTreeNumber();

    if (treenum != oldtreenum) {

      oldtreenum = treenum;
      GlobalCut->UpdateFormulaLeaves();

    }


    if (nevent % 1000 == 0) {

      std::cout
        << "Event " << nevent
        << ", file = "
        << C->GetFile()->GetName()
        << std::endl;

    }


    bool passedcut = GlobalCut->EvalInstance(0) != 0;


    if (passedcut) {

      data.t_vz.push_back(vz[0]);
      data.t_vx.push_back(vx[0]);
      data.t_vy.push_back(vy[0]);

      data.ft_x.push_back(gemFT_x[0]);
      data.ft_y.push_back(gemFT_y[0]);
      data.ft_xp.push_back(gemFT_xp[0]);
      data.ft_yp.push_back(gemFT_yp[0]);

      data.fpp_x.push_back(gemFPP_x[0]);
      data.fpp_y.push_back(gemFPP_y[0]);
      data.fpp_xp.push_back(gemFPP_xp[0]);
      data.fpp_yp.push_back(gemFPP_yp[0]);

    }

    nevent++;
  }


  std::cout << std::endl;
  std::cout << "Total selected events: "
            << data.ft_x.size()
            << std::endl;


  delete GlobalCut;
}

#endif

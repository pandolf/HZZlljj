{
//=========Macro generated from canvas: c_fitfunc_btag0mmc_fitfunc_btag0mm-ext/CANFITFUNC
//=========  (Mon Sep 12 11:12:14 2011) by ROOT version5.27/06b
   TCanvas *c_fitfunc_btag0mmc_fitfunc_btag0mm-ext = new TCanvas("c_fitfunc_btag0mmc_fitfunc_btag0mm-ext", "CANFITFUNC",0,0,900,900);
   c_fitfunc_btag0mmc_fitfunc_btag0mm-ext->Range(81.24999,-0.04636785,768.75,0.1309846);
   c_fitfunc_btag0mmc_fitfunc_btag0mm-ext->SetFillColor(0);
   c_fitfunc_btag0mmc_fitfunc_btag0mm-ext->SetBorderMode(0);
   c_fitfunc_btag0mmc_fitfunc_btag0mm-ext->SetBorderSize(2);
   
   TF1 *fit_poly3_btag0mm = new TF1("fit_poly3_btag0mm","pol3",150,700);
   fit_poly3_btag0mm->SetFillColor(19);
   fit_poly3_btag0mm->SetFillStyle(0);
   fit_poly3_btag0mm->SetLineWidth(3);
   fit_poly3_btag0mm->SetChisquare(9.673839e-06);
   fit_poly3_btag0mm->SetNDF(11);
   fit_poly3_btag0mm->GetXaxis()->SetTitle("M_{H}");
   fit_poly3_btag0mm->GetYaxis()->SetTitle("#varepsilon");
   fit_poly3_btag0mm->SetParameter(0,-0.2493327);
   fit_poly3_btag0mm->SetParError(0,0.008406717);
   fit_poly3_btag0mm->SetParLimits(0,0,0);
   fit_poly3_btag0mm->SetParameter(1,0.001984463);
   fit_poly3_btag0mm->SetParError(1,7.327784e-05);
   fit_poly3_btag0mm->SetParLimits(1,0,0);
   fit_poly3_btag0mm->SetParameter(2,-3.567324e-06);
   fit_poly3_btag0mm->SetParError(2,1.964261e-07);
   fit_poly3_btag0mm->SetParLimits(2,0,0);
   fit_poly3_btag0mm->SetParameter(3,2.032139e-09);
   fit_poly3_btag0mm->SetParError(3,1.652929e-10);
   fit_poly3_btag0mm->SetParLimits(3,0,0);
   fit_poly3_btag0mm->Draw("L");
   
   TGraph *graph = new TGraph(15);
   graph->SetName("effgr_btag0mm");
   graph->SetTitle("Efficiency vs Mass (mm , btag0)");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(21);
   graph->SetPoint(0,190,0.01303832);
   graph->SetPoint(1,200,0.02062848);
   graph->SetPoint(2,210,0.02866847);
   graph->SetPoint(3,230,0.04367485);
   graph->SetPoint(4,250,0.05558928);
   graph->SetPoint(5,300,0.08004429);
   graph->SetPoint(6,350,0.09530594);
   graph->SetPoint(7,400,0.1048654);
   graph->SetPoint(8,425,0.1046195);
   graph->SetPoint(9,450,0.1054648);
   graph->SetPoint(10,500,0.1051818);
   graph->SetPoint(11,525,0.1027745);
   graph->SetPoint(12,550,0.1026393);
   graph->SetPoint(13,575,0.09964557);
   graph->SetPoint(14,600,0.09472636);
   
   TH1F *effgr_btag0mm4__4 = new TH1F("effgr_btag0mm4__4","Efficiency vs Mass (mm , btag0)",100,149,641);
   effgr_btag0mm4__4->SetMinimum(0);
   effgr_btag0mm4__4->SetMaximum(0.1147074);
   effgr_btag0mm4__4->SetDirectory(0);
   effgr_btag0mm4__4->SetStats(0);
   effgr_btag0mm4__4->GetXaxis()->SetTitle("M_{H}");
   effgr_btag0mm4__4->GetYaxis()->SetTitle("#varepsilon");
   graph->SetHistogram(effgr_btag0mm4);
   
   
   TF1 *fit_poly3_btag0mm = new TF1("fit_poly3_btag0mm","pol3",149,641);
   fit_poly3_btag0mm->SetFillColor(19);
   fit_poly3_btag0mm->SetFillStyle(0);
   fit_poly3_btag0mm->SetLineWidth(3);
   fit_poly3_btag0mm->SetChisquare(9.673839e-06);
   fit_poly3_btag0mm->SetNDF(11);
   fit_poly3_btag0mm->SetParameter(0,-0.2493327);
   fit_poly3_btag0mm->SetParError(0,0.008406717);
   fit_poly3_btag0mm->SetParLimits(0,0,0);
   fit_poly3_btag0mm->SetParameter(1,0.001984463);
   fit_poly3_btag0mm->SetParError(1,7.327784e-05);
   fit_poly3_btag0mm->SetParLimits(1,0,0);
   fit_poly3_btag0mm->SetParameter(2,-3.567324e-06);
   fit_poly3_btag0mm->SetParError(2,1.964261e-07);
   fit_poly3_btag0mm->SetParLimits(2,0,0);
   fit_poly3_btag0mm->SetParameter(3,2.032139e-09);
   fit_poly3_btag0mm->SetParError(3,1.652929e-10);
   fit_poly3_btag0mm->SetParLimits(3,0,0);
   graph->GetListOfFunctions()->Add(fit_poly3_btag0mm);
   graph->Draw("p");
   c_fitfunc_btag0mmc_fitfunc_btag0mm-ext->Modified();
   c_fitfunc_btag0mmc_fitfunc_btag0mm-ext->cd();
   c_fitfunc_btag0mmc_fitfunc_btag0mm-ext->SetSelected(c_fitfunc_btag0mmc_fitfunc_btag0mm-ext);
}

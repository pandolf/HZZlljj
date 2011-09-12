{
//=========Macro generated from canvas: c_fitfunc_btag1mmc_fitfunc_btag1mm-ext/CANFITFUNC
//=========  (Mon Sep 12 11:12:14 2011) by ROOT version5.27/06b
   TCanvas *c_fitfunc_btag1mmc_fitfunc_btag1mm-ext = new TCanvas("c_fitfunc_btag1mmc_fitfunc_btag1mm-ext", "CANFITFUNC",0,0,900,900);
   c_fitfunc_btag1mmc_fitfunc_btag1mm-ext->Range(81.24999,-0.01390607,768.75,0.1251547);
   c_fitfunc_btag1mmc_fitfunc_btag1mm-ext->SetFillColor(0);
   c_fitfunc_btag1mmc_fitfunc_btag1mm-ext->SetBorderMode(0);
   c_fitfunc_btag1mmc_fitfunc_btag1mm-ext->SetBorderSize(2);
   
   TF1 *fit_poly3_btag1mm = new TF1("fit_poly3_btag1mm","pol3",150,700);
   fit_poly3_btag1mm->SetFillColor(19);
   fit_poly3_btag1mm->SetFillStyle(0);
   fit_poly3_btag1mm->SetLineWidth(3);
   fit_poly3_btag1mm->SetChisquare(0.0001539626);
   fit_poly3_btag1mm->SetNDF(11);
   fit_poly3_btag1mm->GetXaxis()->SetTitle("M_{H}");
   fit_poly3_btag1mm->GetYaxis()->SetTitle("#varepsilon");
   fit_poly3_btag1mm->SetParameter(0,-0.09034264);
   fit_poly3_btag1mm->SetParError(0,0.0335378);
   fit_poly3_btag1mm->SetParLimits(0,0,0);
   fit_poly3_btag1mm->SetParameter(1,0.0006907409);
   fit_poly3_btag1mm->SetParError(1,0.000292335);
   fit_poly3_btag1mm->SetParLimits(1,0,0);
   fit_poly3_btag1mm->SetParameter(2,-4.590385e-07);
   fit_poly3_btag1mm->SetParError(2,7.836234e-07);
   fit_poly3_btag1mm->SetParLimits(2,0,0);
   fit_poly3_btag1mm->SetParameter(3,-2.759169e-10);
   fit_poly3_btag1mm->SetParError(3,6.594202e-10);
   fit_poly3_btag1mm->SetParLimits(3,0,0);
   fit_poly3_btag1mm->Draw("L");
   
   TGraph *graph = new TGraph(15);
   graph->SetName("effgr_btag1mm");
   graph->SetTitle("Efficiency vs Mass (mm , btag1)");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(21);
   graph->SetPoint(0,190,0.02297734);
   graph->SetPoint(1,200,0.0277686);
   graph->SetPoint(2,210,0.0319374);
   graph->SetPoint(3,230,0.0407416);
   graph->SetPoint(4,250,0.04794084);
   graph->SetPoint(5,300,0.06679446);
   graph->SetPoint(6,350,0.08351642);
   graph->SetPoint(7,400,0.09586941);
   graph->SetPoint(8,425,0.1058206);
   graph->SetPoint(9,450,0.09967036);
   graph->SetPoint(10,500,0.09712458);
   graph->SetPoint(11,525,0.1092896);
   graph->SetPoint(12,550,0.1053757);
   graph->SetPoint(13,575,0.1053006);
   graph->SetPoint(14,600,0.0978125);
   
   TH1F *effgr_btag1mm5__5 = new TH1F("effgr_btag1mm5__5","Efficiency vs Mass (mm , btag1)",100,149,641);
   effgr_btag1mm5__5->SetMinimum(0);
   effgr_btag1mm5__5->SetMaximum(0.1179209);
   effgr_btag1mm5__5->SetDirectory(0);
   effgr_btag1mm5__5->SetStats(0);
   effgr_btag1mm5__5->GetXaxis()->SetTitle("M_{H}");
   effgr_btag1mm5__5->GetYaxis()->SetTitle("#varepsilon");
   graph->SetHistogram(effgr_btag1mm5);
   
   
   TF1 *fit_poly3_btag1mm = new TF1("fit_poly3_btag1mm","pol3",149,641);
   fit_poly3_btag1mm->SetFillColor(19);
   fit_poly3_btag1mm->SetFillStyle(0);
   fit_poly3_btag1mm->SetLineWidth(3);
   fit_poly3_btag1mm->SetChisquare(0.0001539626);
   fit_poly3_btag1mm->SetNDF(11);
   fit_poly3_btag1mm->SetParameter(0,-0.09034264);
   fit_poly3_btag1mm->SetParError(0,0.0335378);
   fit_poly3_btag1mm->SetParLimits(0,0,0);
   fit_poly3_btag1mm->SetParameter(1,0.0006907409);
   fit_poly3_btag1mm->SetParError(1,0.000292335);
   fit_poly3_btag1mm->SetParLimits(1,0,0);
   fit_poly3_btag1mm->SetParameter(2,-4.590385e-07);
   fit_poly3_btag1mm->SetParError(2,7.836234e-07);
   fit_poly3_btag1mm->SetParLimits(2,0,0);
   fit_poly3_btag1mm->SetParameter(3,-2.759169e-10);
   fit_poly3_btag1mm->SetParError(3,6.594202e-10);
   fit_poly3_btag1mm->SetParLimits(3,0,0);
   graph->GetListOfFunctions()->Add(fit_poly3_btag1mm);
   graph->Draw("p");
   c_fitfunc_btag1mmc_fitfunc_btag1mm-ext->Modified();
   c_fitfunc_btag1mmc_fitfunc_btag1mm-ext->cd();
   c_fitfunc_btag1mmc_fitfunc_btag1mm-ext->SetSelected(c_fitfunc_btag1mmc_fitfunc_btag1mm-ext);
}

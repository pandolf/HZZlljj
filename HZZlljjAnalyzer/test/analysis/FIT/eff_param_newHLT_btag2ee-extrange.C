{
//=========Macro generated from canvas: c_fitfunc_btag2eec_fitfunc_btag2ee-ext/CANFITFUNC
//=========  (Mon Sep 12 11:12:14 2011) by ROOT version5.27/06b
   TCanvas *c_fitfunc_btag2eec_fitfunc_btag2ee-ext = new TCanvas("c_fitfunc_btag2eec_fitfunc_btag2ee-ext", "CANFITFUNC",0,0,900,900);
   c_fitfunc_btag2eec_fitfunc_btag2ee-ext->Range(81.24999,-0.014756,768.75,0.05287109);
   c_fitfunc_btag2eec_fitfunc_btag2ee-ext->SetFillColor(0);
   c_fitfunc_btag2eec_fitfunc_btag2ee-ext->SetBorderMode(0);
   c_fitfunc_btag2eec_fitfunc_btag2ee-ext->SetBorderSize(2);
   
   TF1 *fit_poly3_btag2ee = new TF1("fit_poly3_btag2ee","pol3",150,700);
   fit_poly3_btag2ee->SetFillColor(19);
   fit_poly3_btag2ee->SetFillStyle(0);
   fit_poly3_btag2ee->SetLineWidth(3);
   fit_poly3_btag2ee->SetChisquare(3.200862e-05);
   fit_poly3_btag2ee->SetNDF(11);
   fit_poly3_btag2ee->GetXaxis()->SetTitle("M_{H}");
   fit_poly3_btag2ee->GetYaxis()->SetTitle("#varepsilon");
   fit_poly3_btag2ee->SetParameter(0,-0.04310413);
   fit_poly3_btag2ee->SetParError(0,0.01529186);
   fit_poly3_btag2ee->SetParLimits(0,0,0);
   fit_poly3_btag2ee->SetParameter(1,0.0002520757);
   fit_poly3_btag2ee->SetParError(1,0.0001332928);
   fit_poly3_btag2ee->SetParLimits(1,0,0);
   fit_poly3_btag2ee->SetParameter(2,1.139391e-08);
   fit_poly3_btag2ee->SetParError(2,3.573002e-07);
   fit_poly3_btag2ee->SetParLimits(2,0,0);
   fit_poly3_btag2ee->SetParameter(3,-3.383753e-10);
   fit_poly3_btag2ee->SetParError(3,3.006686e-10);
   fit_poly3_btag2ee->SetParLimits(3,0,0);
   fit_poly3_btag2ee->Draw("L");
   
   TGraph *graph = new TGraph(15);
   graph->SetName("effgr_btag2ee");
   graph->SetTitle("Efficiency vs Mass (ee , btag2)");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(21);
   graph->SetPoint(0,190,0.003592738);
   graph->SetPoint(1,200,0.005120317);
   graph->SetPoint(2,210,0.006829595);
   graph->SetPoint(3,230,0.009961336);
   graph->SetPoint(4,250,0.0157118);
   graph->SetPoint(5,300,0.02589709);
   graph->SetPoint(6,350,0.0312076);
   graph->SetPoint(7,400,0.03762599);
   graph->SetPoint(8,425,0.03783302);
   graph->SetPoint(9,450,0.04465634);
   graph->SetPoint(10,500,0.04565562);
   graph->SetPoint(11,525,0.04058707);
   graph->SetPoint(12,550,0.04241438);
   graph->SetPoint(13,575,0.04134965);
   graph->SetPoint(14,600,0.03963849);
   
   TH1F *effgr_btag2ee3__3 = new TH1F("effgr_btag2ee3__3","Efficiency vs Mass (ee , btag2)",100,149,641);
   effgr_btag2ee3__3->SetMinimum(0);
   effgr_btag2ee3__3->SetMaximum(0.04986191);
   effgr_btag2ee3__3->SetDirectory(0);
   effgr_btag2ee3__3->SetStats(0);
   effgr_btag2ee3__3->GetXaxis()->SetTitle("M_{H}");
   effgr_btag2ee3__3->GetYaxis()->SetTitle("#varepsilon");
   graph->SetHistogram(effgr_btag2ee3);
   
   
   TF1 *fit_poly3_btag2ee = new TF1("fit_poly3_btag2ee","pol3",149,641);
   fit_poly3_btag2ee->SetFillColor(19);
   fit_poly3_btag2ee->SetFillStyle(0);
   fit_poly3_btag2ee->SetLineWidth(3);
   fit_poly3_btag2ee->SetChisquare(3.200862e-05);
   fit_poly3_btag2ee->SetNDF(11);
   fit_poly3_btag2ee->SetParameter(0,-0.04310413);
   fit_poly3_btag2ee->SetParError(0,0.01529186);
   fit_poly3_btag2ee->SetParLimits(0,0,0);
   fit_poly3_btag2ee->SetParameter(1,0.0002520757);
   fit_poly3_btag2ee->SetParError(1,0.0001332928);
   fit_poly3_btag2ee->SetParLimits(1,0,0);
   fit_poly3_btag2ee->SetParameter(2,1.139391e-08);
   fit_poly3_btag2ee->SetParError(2,3.573002e-07);
   fit_poly3_btag2ee->SetParLimits(2,0,0);
   fit_poly3_btag2ee->SetParameter(3,-3.383753e-10);
   fit_poly3_btag2ee->SetParError(3,3.006686e-10);
   fit_poly3_btag2ee->SetParLimits(3,0,0);
   graph->GetListOfFunctions()->Add(fit_poly3_btag2ee);
   graph->Draw("p");
   c_fitfunc_btag2eec_fitfunc_btag2ee-ext->Modified();
   c_fitfunc_btag2eec_fitfunc_btag2ee-ext->cd();
   c_fitfunc_btag2eec_fitfunc_btag2ee-ext->SetSelected(c_fitfunc_btag2eec_fitfunc_btag2ee-ext);
}

{
//=========Macro generated from canvas: c_fitfunc_btag2mmc_fitfunc_btag2mm-ext/CANFITFUNC
//=========  (Mon Sep 12 11:12:14 2011) by ROOT version5.27/06b
   TCanvas *c_fitfunc_btag2mmc_fitfunc_btag2mm-ext = new TCanvas("c_fitfunc_btag2mmc_fitfunc_btag2mm-ext", "CANFITFUNC",0,0,900,900);
   c_fitfunc_btag2mmc_fitfunc_btag2mm-ext->Range(81.24999,-0.01603383,768.75,0.05813761);
   c_fitfunc_btag2mmc_fitfunc_btag2mm-ext->SetFillColor(0);
   c_fitfunc_btag2mmc_fitfunc_btag2mm-ext->SetBorderMode(0);
   c_fitfunc_btag2mmc_fitfunc_btag2mm-ext->SetBorderSize(2);
   
   TF1 *fit_poly3_btag2mm = new TF1("fit_poly3_btag2mm","pol3",150,700);
   fit_poly3_btag2mm->SetFillColor(19);
   fit_poly3_btag2mm->SetFillStyle(0);
   fit_poly3_btag2mm->SetLineWidth(3);
   fit_poly3_btag2mm->SetChisquare(9.375312e-06);
   fit_poly3_btag2mm->SetNDF(11);
   fit_poly3_btag2mm->GetXaxis()->SetTitle("M_{H}");
   fit_poly3_btag2mm->GetYaxis()->SetTitle("#varepsilon");
   fit_poly3_btag2mm->SetParameter(0,-0.05164656);
   fit_poly3_btag2mm->SetParError(0,0.008275988);
   fit_poly3_btag2mm->SetParLimits(0,0,0);
   fit_poly3_btag2mm->SetParameter(1,0.0003299199);
   fit_poly3_btag2mm->SetParError(1,7.213833e-05);
   fit_poly3_btag2mm->SetParLimits(1,0,0);
   fit_poly3_btag2mm->SetParameter(2,-1.72873e-07);
   fit_poly3_btag2mm->SetParError(2,1.933716e-07);
   fit_poly3_btag2mm->SetParLimits(2,0,0);
   fit_poly3_btag2mm->SetParameter(3,-1.797923e-10);
   fit_poly3_btag2mm->SetParError(3,1.627225e-10);
   fit_poly3_btag2mm->SetParLimits(3,0,0);
   fit_poly3_btag2mm->Draw("L");
   
   TGraph *graph = new TGraph(15);
   graph->SetName("effgr_btag2mm");
   graph->SetTitle("Efficiency vs Mass (mm , btag2)");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(21);
   graph->SetPoint(0,190,0.004266586);
   graph->SetPoint(1,200,0.005902092);
   graph->SetPoint(2,210,0.007932441);
   graph->SetPoint(3,230,0.01171251);
   graph->SetPoint(4,250,0.01808697);
   graph->SetPoint(5,300,0.02649332);
   graph->SetPoint(6,350,0.03636048);
   graph->SetPoint(7,400,0.04092116);
   graph->SetPoint(8,425,0.04289337);
   graph->SetPoint(9,450,0.0454602);
   graph->SetPoint(10,500,0.0462595);
   graph->SetPoint(11,525,0.0485992);
   graph->SetPoint(12,550,0.04840192);
   graph->SetPoint(13,575,0.04729116);
   graph->SetPoint(14,600,0.0444953);
   
   TH1F *effgr_btag2mm6__6 = new TH1F("effgr_btag2mm6__6","Efficiency vs Mass (mm , btag2)",100,149,641);
   effgr_btag2mm6__6->SetMinimum(0);
   effgr_btag2mm6__6->SetMaximum(0.05303246);
   effgr_btag2mm6__6->SetDirectory(0);
   effgr_btag2mm6__6->SetStats(0);
   effgr_btag2mm6__6->GetXaxis()->SetTitle("M_{H}");
   effgr_btag2mm6__6->GetYaxis()->SetTitle("#varepsilon");
   graph->SetHistogram(effgr_btag2mm6);
   
   
   TF1 *fit_poly3_btag2mm = new TF1("fit_poly3_btag2mm","pol3",149,641);
   fit_poly3_btag2mm->SetFillColor(19);
   fit_poly3_btag2mm->SetFillStyle(0);
   fit_poly3_btag2mm->SetLineWidth(3);
   fit_poly3_btag2mm->SetChisquare(9.375312e-06);
   fit_poly3_btag2mm->SetNDF(11);
   fit_poly3_btag2mm->SetParameter(0,-0.05164656);
   fit_poly3_btag2mm->SetParError(0,0.008275988);
   fit_poly3_btag2mm->SetParLimits(0,0,0);
   fit_poly3_btag2mm->SetParameter(1,0.0003299199);
   fit_poly3_btag2mm->SetParError(1,7.213833e-05);
   fit_poly3_btag2mm->SetParLimits(1,0,0);
   fit_poly3_btag2mm->SetParameter(2,-1.72873e-07);
   fit_poly3_btag2mm->SetParError(2,1.933716e-07);
   fit_poly3_btag2mm->SetParLimits(2,0,0);
   fit_poly3_btag2mm->SetParameter(3,-1.797923e-10);
   fit_poly3_btag2mm->SetParError(3,1.627225e-10);
   fit_poly3_btag2mm->SetParLimits(3,0,0);
   graph->GetListOfFunctions()->Add(fit_poly3_btag2mm);
   graph->Draw("p");
   c_fitfunc_btag2mmc_fitfunc_btag2mm-ext->Modified();
   c_fitfunc_btag2mmc_fitfunc_btag2mm-ext->cd();
   c_fitfunc_btag2mmc_fitfunc_btag2mm-ext->SetSelected(c_fitfunc_btag2mmc_fitfunc_btag2mm-ext);
}

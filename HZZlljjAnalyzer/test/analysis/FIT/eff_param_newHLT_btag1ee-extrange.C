{
//=========Macro generated from canvas: c_fitfunc_btag1eec_fitfunc_btag1ee-ext/CANFITFUNC
//=========  (Mon Sep 12 11:12:14 2011) by ROOT version5.27/06b
   TCanvas *c_fitfunc_btag1eec_fitfunc_btag1ee-ext = new TCanvas("c_fitfunc_btag1eec_fitfunc_btag1ee-ext", "CANFITFUNC",0,0,900,900);
   c_fitfunc_btag1eec_fitfunc_btag1ee-ext->Range(81.24999,-0.0125434,768.75,0.1128906);
   c_fitfunc_btag1eec_fitfunc_btag1ee-ext->SetFillColor(0);
   c_fitfunc_btag1eec_fitfunc_btag1ee-ext->SetBorderMode(0);
   c_fitfunc_btag1eec_fitfunc_btag1ee-ext->SetBorderSize(2);
   
   TF1 *fit_poly3_btag1ee = new TF1("fit_poly3_btag1ee","pol3",150,700);
   fit_poly3_btag1ee->SetFillColor(19);
   fit_poly3_btag1ee->SetFillStyle(0);
   fit_poly3_btag1ee->SetLineWidth(3);
   fit_poly3_btag1ee->SetChisquare(3.753885e-05);
   fit_poly3_btag1ee->SetNDF(11);
   fit_poly3_btag1ee->GetXaxis()->SetTitle("M_{H}");
   fit_poly3_btag1ee->GetYaxis()->SetTitle("#varepsilon");
   fit_poly3_btag1ee->SetParameter(0,-0.08651818);
   fit_poly3_btag1ee->SetParError(0,0.01656027);
   fit_poly3_btag1ee->SetParLimits(0,0,0);
   fit_poly3_btag1ee->SetParameter(1,0.000639916);
   fit_poly3_btag1ee->SetParError(1,0.000144349);
   fit_poly3_btag1ee->SetParLimits(1,0,0);
   fit_poly3_btag1ee->SetParameter(2,-3.873668e-07);
   fit_poly3_btag1ee->SetParError(2,3.86937e-07);
   fit_poly3_btag1ee->SetParLimits(2,0,0);
   fit_poly3_btag1ee->SetParameter(3,-3.283154e-10);
   fit_poly3_btag1ee->SetParError(3,3.25608e-10);
   fit_poly3_btag1ee->SetParLimits(3,0,0);
   fit_poly3_btag1ee->Draw("L");
   
   TGraph *graph = new TGraph(15);
   graph->SetName("effgr_btag1ee");
   graph->SetTitle("Efficiency vs Mass (ee , btag1)");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(21);
   graph->SetPoint(0,190,0.0196389);
   graph->SetPoint(1,200,0.02317839);
   graph->SetPoint(2,210,0.02743467);
   graph->SetPoint(3,230,0.03633261);
   graph->SetPoint(4,250,0.04280739);
   graph->SetPoint(5,300,0.06202897);
   graph->SetPoint(6,350,0.07660714);
   graph->SetPoint(7,400,0.08571377);
   graph->SetPoint(8,425,0.09338933);
   graph->SetPoint(9,450,0.09142358);
   graph->SetPoint(10,500,0.09189457);
   graph->SetPoint(11,525,0.09780505);
   graph->SetPoint(12,550,0.09353309);
   graph->SetPoint(13,575,0.09177141);
   graph->SetPoint(14,600,0.08648302);
   
   TH1F *effgr_btag1ee2__2 = new TH1F("effgr_btag1ee2__2","Efficiency vs Mass (ee , btag1)",100,149,641);
   effgr_btag1ee2__2->SetMinimum(0);
   effgr_btag1ee2__2->SetMaximum(0.1056217);
   effgr_btag1ee2__2->SetDirectory(0);
   effgr_btag1ee2__2->SetStats(0);
   effgr_btag1ee2__2->GetXaxis()->SetTitle("M_{H}");
   effgr_btag1ee2__2->GetYaxis()->SetTitle("#varepsilon");
   graph->SetHistogram(effgr_btag1ee2);
   
   
   TF1 *fit_poly3_btag1ee = new TF1("fit_poly3_btag1ee","pol3",149,641);
   fit_poly3_btag1ee->SetFillColor(19);
   fit_poly3_btag1ee->SetFillStyle(0);
   fit_poly3_btag1ee->SetLineWidth(3);
   fit_poly3_btag1ee->SetChisquare(3.753885e-05);
   fit_poly3_btag1ee->SetNDF(11);
   fit_poly3_btag1ee->SetParameter(0,-0.08651818);
   fit_poly3_btag1ee->SetParError(0,0.01656027);
   fit_poly3_btag1ee->SetParLimits(0,0,0);
   fit_poly3_btag1ee->SetParameter(1,0.000639916);
   fit_poly3_btag1ee->SetParError(1,0.000144349);
   fit_poly3_btag1ee->SetParLimits(1,0,0);
   fit_poly3_btag1ee->SetParameter(2,-3.873668e-07);
   fit_poly3_btag1ee->SetParError(2,3.86937e-07);
   fit_poly3_btag1ee->SetParLimits(2,0,0);
   fit_poly3_btag1ee->SetParameter(3,-3.283154e-10);
   fit_poly3_btag1ee->SetParError(3,3.25608e-10);
   fit_poly3_btag1ee->SetParLimits(3,0,0);
   graph->GetListOfFunctions()->Add(fit_poly3_btag1ee);
   graph->Draw("p");
   c_fitfunc_btag1eec_fitfunc_btag1ee-ext->Modified();
   c_fitfunc_btag1eec_fitfunc_btag1ee-ext->cd();
   c_fitfunc_btag1eec_fitfunc_btag1ee-ext->SetSelected(c_fitfunc_btag1eec_fitfunc_btag1ee-ext);
}

{
//=========Macro generated from canvas: c_fitfunc_btag0eec_fitfunc_btag0ee-ext/CANFITFUNC
//=========  (Mon Sep 12 11:12:14 2011) by ROOT version5.27/06b
   TCanvas *c_fitfunc_btag0eec_fitfunc_btag0ee-ext = new TCanvas("c_fitfunc_btag0eec_fitfunc_btag0ee-ext", "CANFITFUNC",0,0,900,900);
   c_fitfunc_btag0eec_fitfunc_btag0ee-ext->Range(81.24999,-0.04997628,768.75,0.1222031);
   c_fitfunc_btag0eec_fitfunc_btag0ee-ext->SetFillColor(0);
   c_fitfunc_btag0eec_fitfunc_btag0ee-ext->SetBorderMode(0);
   c_fitfunc_btag0eec_fitfunc_btag0ee-ext->SetBorderSize(2);
   
   TF1 *fit_poly3_btag0ee = new TF1("fit_poly3_btag0ee","pol3",150,700);
   fit_poly3_btag0ee->SetFillColor(19);
   fit_poly3_btag0ee->SetFillStyle(0);
   fit_poly3_btag0ee->SetLineWidth(3);
   fit_poly3_btag0ee->SetChisquare(6.126377e-05);
   fit_poly3_btag0ee->SetNDF(11);
   fit_poly3_btag0ee->GetXaxis()->SetTitle("M_{H}");
   fit_poly3_btag0ee->GetYaxis()->SetTitle("#varepsilon");
   fit_poly3_btag0ee->SetParameter(0,-0.2613789);
   fit_poly3_btag0ee->SetParError(0,0.02115577);
   fit_poly3_btag0ee->SetParLimits(0,0,0);
   fit_poly3_btag0ee->SetParameter(1,0.002071747);
   fit_poly3_btag0ee->SetParError(1,0.000184406);
   fit_poly3_btag0ee->SetParLimits(1,0,0);
   fit_poly3_btag0ee->SetParameter(2,-3.840909e-06);
   fit_poly3_btag0ee->SetParError(2,4.943126e-07);
   fit_poly3_btag0ee->SetParLimits(2,0,0);
   fit_poly3_btag0ee->SetParameter(3,2.252034e-09);
   fit_poly3_btag0ee->SetParError(3,4.159647e-10);
   fit_poly3_btag0ee->SetParLimits(3,0,0);
   fit_poly3_btag0ee->Draw("L");
   
   TGraph *graph = new TGraph(15);
   graph->SetName("effgr_btag0ee");
   graph->SetTitle("Efficiency vs Mass (ee , btag0)");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(21);
   graph->SetPoint(0,190,0.01092572);
   graph->SetPoint(1,200,0.01849003);
   graph->SetPoint(2,210,0.023671);
   graph->SetPoint(3,230,0.03649743);
   graph->SetPoint(4,250,0.05117248);
   graph->SetPoint(5,300,0.07607853);
   graph->SetPoint(6,350,0.09178921);
   graph->SetPoint(7,400,0.09928678);
   graph->SetPoint(8,425,0.09322475);
   graph->SetPoint(9,450,0.1001037);
   graph->SetPoint(10,500,0.0968918);
   graph->SetPoint(11,525,0.09116618);
   graph->SetPoint(12,550,0.09133907);
   graph->SetPoint(13,575,0.0887107);
   graph->SetPoint(14,600,0.08541508);
   
   TH1F *effgr_btag0ee1__1 = new TH1F("effgr_btag0ee1__1","Efficiency vs Mass (ee , btag0)",100,149,641);
   effgr_btag0ee1__1->SetMinimum(0);
   effgr_btag0ee1__1->SetMaximum(0.1090215);
   effgr_btag0ee1__1->SetDirectory(0);
   effgr_btag0ee1__1->SetStats(0);
   effgr_btag0ee1__1->GetXaxis()->SetTitle("M_{H}");
   effgr_btag0ee1__1->GetYaxis()->SetTitle("#varepsilon");
   graph->SetHistogram(effgr_btag0ee1);
   
   
   TF1 *fit_poly3_btag0ee = new TF1("fit_poly3_btag0ee","pol3",149,641);
   fit_poly3_btag0ee->SetFillColor(19);
   fit_poly3_btag0ee->SetFillStyle(0);
   fit_poly3_btag0ee->SetLineWidth(3);
   fit_poly3_btag0ee->SetChisquare(6.126377e-05);
   fit_poly3_btag0ee->SetNDF(11);
   fit_poly3_btag0ee->SetParameter(0,-0.2613789);
   fit_poly3_btag0ee->SetParError(0,0.02115577);
   fit_poly3_btag0ee->SetParLimits(0,0,0);
   fit_poly3_btag0ee->SetParameter(1,0.002071747);
   fit_poly3_btag0ee->SetParError(1,0.000184406);
   fit_poly3_btag0ee->SetParLimits(1,0,0);
   fit_poly3_btag0ee->SetParameter(2,-3.840909e-06);
   fit_poly3_btag0ee->SetParError(2,4.943126e-07);
   fit_poly3_btag0ee->SetParLimits(2,0,0);
   fit_poly3_btag0ee->SetParameter(3,2.252034e-09);
   fit_poly3_btag0ee->SetParError(3,4.159647e-10);
   fit_poly3_btag0ee->SetParLimits(3,0,0);
   graph->GetListOfFunctions()->Add(fit_poly3_btag0ee);
   graph->Draw("p");
   c_fitfunc_btag0eec_fitfunc_btag0ee-ext->Modified();
   c_fitfunc_btag0eec_fitfunc_btag0ee-ext->cd();
   c_fitfunc_btag0eec_fitfunc_btag0ee-ext->SetSelected(c_fitfunc_btag0eec_fitfunc_btag0ee-ext);
}

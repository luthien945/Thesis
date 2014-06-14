{
	gStyle->SetOptStat(0);

	TFile *f11 = new TFile("./acc_QA.root");
	f11->cd();
	TH1D *acc11 = (TH1D *) f11->Get("AcceptanceC");
	TH1D *Rback11 = (TH1D *) f11->Get("back_Ratio");

	TFile *f10 = new TFile("../Run10/acc_QA.root");
	f10->cd();
	TH1D *acc10 = (TH1D *) f10->Get("AcceptanceC");
	TH1D *Rback10 = (TH1D *) f10->Get("back_Ratio");

	double LX1=0.0;
	double HX1=5.7;
	double LY1=0.92;
	double HY1=1.08;	

	double LX2=0.0;
	double HX2=5.7;
	double LY2=0.9;
	double HY2=1.05;

	double mX1=LX1+0.1*(HX1-LX1); 
	double mY1=LY1+0.8*(HY1-LY1); 

	double mX2=LX2+0.1*(HX2-LX2); 
	double mY2=LY2+0.8*(HY2-LY2); 

	TLine *Line0 = new TLine(LX1,1.,HX1,1.);
	Line0->SetLineStyle(7);
	Line0->SetLineWidth(2);

	TLine *LineX1 = new TLine(LX1,LY2,LX1,HY2);
	LineX1->SetLineStyle(1);
	LineX1->SetLineWidth(3);

	TLine *LineX2 = new TLine(HX1,LY2,HX1,HY2);
	LineX2->SetLineStyle(1);
	LineX2->SetLineWidth(3);

	TLine *LineY1 = new TLine(LX1,LY2,HX1,LY2);
	LineY1->SetLineStyle(1);
	LineY1->SetLineWidth(3);

	TLine *LineY2 = new TLine(LX1,HY2,HX1,HY2);
	LineY2->SetLineStyle(1);
	LineY2->SetLineWidth(3);

	TCanvas *c1 = new TCanvas("c1","acc",800,600);
	c1->cd();
	TH2F *frame1 = new TH2F("frame1","",1,LX1,HX1,1,LY2,HY2);
	frame1->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frame1->GetYaxis()->SetTitle("Acceptence Factor");
	frame1->GetYaxis()->SetNdivisions(503);
	frame1->GetYaxis()->SetTitleOffset(1);
	frame1->SetLineWidth(2);
	frame1->Draw();
	Line0->Draw("same");
	acc10->SetLineColor(2);
	acc10->Draw("epsame");
	acc11->SetMarkerColor(1);
	acc11->SetMarkerStyle(4);
	acc11->SetLineColor(1);
	acc11->SetLineWidth(2);
	acc11_1 = (TH1D *)acc11->Clone("acc11_1");
	acc11_1->SetMarkerSize(1.9);
	acc11->Draw("epsame");
	acc11_1->Draw("epsame");
	LineX1->Draw("same");
	LineX2->Draw("same");
	LineY1->Draw("same");
	LineY2->Draw("same");
	TLegend *t1 = new TLegend(0.2,0.3,0.4,0.5);
	t1->SetFillColor(10);
	t1->SetBorderSize(0);
	t1->SetTextSize(0.035);
	t1->AddEntry(acc10,"Run10","epl");	
	t1->AddEntry(acc11,"Run11","epl");	
	t1->Draw("same");
	c1->SaveAs("acc_compare.gif");


	TCanvas *c2 = new TCanvas("c2","Rback",800,600);
	c2->cd();
	TH2F *frame2 = new TH2F("frame2","",1,LX1,HX1,1,LY2,HY2);
	frame2->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frame2->GetYaxis()->SetTitle("LikeSign/MixedEvent_unLikeSign");
	frame2->SetLineWidth(2);
	frame2->Draw();
	Line0->Draw("same");
	Rback10->SetLineColor(2);
	Rback10->Draw("epsame");
	Rback11->SetMarkerColor(1);
	Rback11->SetMarkerStyle(4);
	Rback11->SetLineColor(1);
	Rback11->SetLineWidth(2);
	Rback11_1 = (TH1D *)Rback11->Clone("Rback11_1");
	Rback11_1->SetMarkerSize(1.6);
	Rback11->Draw("epsame");
	Rback11_1->Draw("epsame");
	LineX1->Draw("same");
	LineX2->Draw("same");
	LineY1->Draw("same");
	LineY2->Draw("same");
	TLegend *t2 = new TLegend(0.2,0.3,0.4,0.5);
	t2->SetFillColor(10);
	t2->SetBorderSize(0);
	t2->SetTextSize(0.035);
	t2->AddEntry(Rback10,"Run10","epl");	
	t2->AddEntry(Rback11,"Run11","epl");	
	t2->Draw("same");
	c2->SaveAs("Rback_compare.gif");
}

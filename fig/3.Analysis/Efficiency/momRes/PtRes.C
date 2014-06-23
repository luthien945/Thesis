#include "TObjArray.h"
//#define QA
//#define EFFICIENCY
Double_t Res(Double_t *x, Double_t *par){
	return TMath::Sqrt(par[0]*par[0]*x[0]*x[0]+par[1]*par[1]);
}

///////////////////////////////////////////////////////////////////////////////
// Crystal-Ball function on both sides
// par[0] - N, par[1] - mean , par[2] - sigma, par[3] - n, par[4] - alpha
// par[5] - m, par[6] - beta
///////////////////////////////////////////////////////////////////////////////
Double_t CrystalBall2(Double_t *x, Double_t *par)
{
  Double_t N = par[0];
  Double_t mu = par[1];
  Double_t s = par[2];
  Double_t n = par[3];
  Double_t alpha = par[4];
  Double_t m = par[5];
  Double_t beta = par[6];
    
  Double_t A = TMath::Power(n/fabs(alpha), n) * TMath::Exp(-alpha*alpha/2.);
  Double_t B = n/fabs(alpha) - fabs(alpha);
  
  Double_t C = TMath::Power(m/fabs(beta), m) * TMath::Exp(-beta*beta/2.);
  Double_t D = m/fabs(beta) - fabs(beta);

  Double_t norm = (x[0]-mu)/s;

  if(norm < -alpha) {
    return N * A * TMath::Power(B-norm, -n);
  } else if(norm < beta) {
    return N * TMath::Exp(-0.5*norm*norm);
  } else { 
    return N * C * TMath::Power(D+norm, -m);
  }
}

void PtRes(){
	gROOT->Reset();
	gROOT->LoadMacro("rootinit.C");
	rootinit();
	gROOT->LoadMacro("/Users/luthien/worksdir/Run11AuAu200/Effciency/Ratio_Err.C");
	// open embedding file
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
  char *infileP = "./Positron_total.root";
  char *infileN = "./Electron_total.root";
	TFile *f1 = new TFile(infileN);
	f1->cd();
	TH2F *hPtDiff       = (TH2F *) f1->Get("PtDiff");
	TH2F *hPtDiffScale  = (TH2F *) f1->Get("PtDiffScale");
	//TFile *f2 = new TFile(infileP);
	//f2->cd();
	//hPtDiff->Add((TH2F *) f2->Get("PtDiff"));
	//hPtDiffScale->Add((TH1F *)f2->Get("PtDiffScale"));

	TObjArray PtDiff_Fit;
	TF1 *g1 = new TF1("g1","gaus",-0.5,0.5);
	hPtDiff->FitSlicesY(g1,21,401,0,"QNR",&PtDiff_Fit);
	
	TH1D *fit[3];
	double max[3] = {1800.,0.0035,0.025};
	double min[3] = {0.,-0.001,0.006};
	TCanvas *c1 = new TCanvas("c1","",800,600);
	//c1->Divide(3,1);
	TF1 *myresfit = new TF1("myresfit",Res,0.2,4,2);
  myresfit->SetParName(0, "a");
  myresfit->SetParName(1, "b");
	myresfit->SetParameter(0,3.74e-2);
	myresfit->SetParameter(1,8.2e-2);
	fit[2]->Fit("myresfit","er","",0.2,4);
  c1->cd();
  fit[2]->SetMaximum(0.025);
  fit[2]->GetYaxis()->SetTitleOffset(1.1);
  fit[2]->GetYaxis()->SetTitle("#sigma_{p_{T}}/p_{T}");
	fit[2]->GetXaxis()->SetTitle("p_{T}(GeV/c)");
  fit[2]->Draw("ep");
  myresfit->Draw("same");
	c1->SaveAs("PtRes.gif");
	c1->SaveAs("PtRes.eps");

	TCanvas *c2 = new TCanvas("c2","",800,600);
	c2->SetLogy();
	c2->SetGridx();
	c2->SetGridy();
	c2->cd();
	TH1F *hPtDiffScale_PY = (TH1F *) hPtDiffScale->ProjectionY();
	TF1 *momshape = new TF1("momshape",CrystalBall2,-0.2,0.2,7);
  momshape->SetParName(0,"Counts");
  momshape->SetParName(1,"#mu");
  momshape->SetParName(2,"#sigma_{p_{T}}/p_{T}");
  momshape->SetParName(3,"n");
  momshape->SetParName(4,"#alpha");
  momshape->SetParName(5,"m");
  momshape->SetParName(6,"#beta");
	momshape->SetParameters(10000., -1e-4, 0.01, 1.35, 1.8, 3.3, 1.8);
	hPtDiffScale_PY->Fit(momshape,"NOR","",-0.2,0.2);
	double par[7];
	momshape->GetParameters(par);
	momshape->SetParameters(par);
	hPtDiffScale_PY->Fit(momshape,"NOR","",-0.2,0.2);
	momshape->GetParameters(par);
	momshape->SetParameters(par);
	hPtDiffScale_PY->GetXaxis()->SetTitle("(p_{T}(RC)-p_{T}(MC))/p_{T}(MC)");
	hPtDiffScale_PY->GetXaxis()->SetRangeUser(-0.2,0.2);
	hPtDiffScale_PY->Fit(momshape,"MER","",-0.2,0.2);
	hPtDiffScale_PY->Draw("ep");
	c2->SaveAs("momShape.gif");
	c2->SaveAs("momShape.eps");

}

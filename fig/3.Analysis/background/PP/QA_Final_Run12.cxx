{
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
//#include "../sb.h"

	////////////////////
	//  backGround
	const Int_t mback=2;    // 1: LikeSign, 2: MixEvent, 3: LikeSign Normalized MixEvent
	const Int_t mAcc=1;     // 1: 1D afig0eptance correction  2: 2D 

	const double pt1=0.;
	const double pt2=5;
	const double m1=0.;
	const double m2=3;

	const Int_t cent1=2;
	const Int_t cent2=10;
	///////////////////
	// normalization range
	//double lowM=1.0;  
	//double hiM=2.00;
	double lowM=0.4;  
	double hiM=1.;
	double lowMls=0.4;  // Low limit for likeS background  

	char TitleC[256]="MinBias";

	////////////////////////////////////////////////////////////////////////////////////////////

	gROOT->Macro("rootinit.C");
	//gROOT->Macro("/Users/jiezhao/run10/Final/eeStyle.C");
	//gROOT->Macro("/home/lenovo/Desktop/Windows/jiezhao/jiezhao/run10/Final/eeStyle.C");

	gStyle->SetOptTitle(0);
	gStyle->SetPadTickX(0);
	gStyle->SetPadTickY(0);
	//gStyle->SetErrorX(0); 
	//gStyle->SetEndErrorSize(10); 

	//TFile *f=new TFile("../../rootfile/200mb_star.root");
	TFile *f = new TFile("../rootfile/RawHist/Run12pp200_diE_MB.root");

	TH1D *Refmult;
	//Refmult =(TH1D *) f->Get("RefMult_c");

	TH1D *Hcentrality;
	Hcentrality =(TH1D *) f->Get("Central");
	double Nevents=Hcentrality->Integral(cent1,cent2);	
	cout<<Hcentrality->GetEntries()<<endl;	
	cout<<Hcentrality->Integral()<<endl;	

	TH2D *RepMass = (TH2D *)f->Get("MassVsPtVsCen_npe_unLikeSign");               //  x,mass;  y,pt;   z,cent
	TH2D *ReeMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignNN");      
	TH2D *RppMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignPP");     	
	TH2D *MepMass = (TH2D *)f->Get("MassVsPtVsCen_npe_unLikesign_Mix"); 
	TH2D *MeeMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignNN_Mix"); 
	TH2D *MppMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignPP_Mix"); 		

	//RepMass->GetYaxis()->SetRangeUser(pt1,pt2);
	//ReeMass->GetYaxis()->SetRangeUser(pt1,pt2);
	//RppMass->GetYaxis()->SetRangeUser(pt1,pt2);
	//MepMass->GetYaxis()->SetRangeUser(pt1,pt2);
	//MeeMass->GetYaxis()->SetRangeUser(pt1,pt2);
	//MppMass->GetYaxis()->SetRangeUser(pt1,pt2);

	// for 2D 
	TH2D *h1_mp=(TH2D *) RepMass->Clone("h1mp");      //x mass;  y  pt
	TH2D *h2_mp=(TH2D *) ReeMass->Clone("h2mp");
	TH2D *h3_mp=(TH2D *) RppMass->Clone("h3mp");
	TH2D *h4_mp=(TH2D *) MepMass->Clone("h4mp");
	TH2D *h5_mp=(TH2D *) MeeMass->Clone("h5mp");
	TH2D *h6_mp=(TH2D *) MppMass->Clone("h6mp");

	// for 1D
	int ptb1=RepMass->GetYaxis()->FindBin(pt1*1.00001);
	int ptb2=RepMass->GetYaxis()->FindBin(pt2*1.00001)-1;
	TH1D *new11=(TH1D *) RepMass->ProjectionX("new11",ptb1,ptb2,"e");
	TH1D *new22=(TH1D *) ReeMass->ProjectionX("new22",ptb1,ptb2,"e");
	TH1D *new33=(TH1D *) RppMass->ProjectionX("new33",ptb1,ptb2,"e");
	TH1D *new44=(TH1D *) MepMass->ProjectionX("new44",ptb1,ptb2,"e");
	TH1D *new55=(TH1D *) MeeMass->ProjectionX("new55",ptb1,ptb2,"e");
	TH1D *new66=(TH1D *) MppMass->ProjectionX("new66",ptb1,ptb2,"e");

	TH1D *n11=(TH1D *) RepMass->ProjectionX("n11",ptb1,ptb2,"e");
	TH1D *n22=(TH1D *) ReeMass->ProjectionX("n22",ptb1,ptb2,"e");
	TH1D *n33=(TH1D *) RppMass->ProjectionX("n33",ptb1,ptb2,"e");
	TH1D *n44=(TH1D *) MepMass->ProjectionX("n44",ptb1,ptb2,"e");
	TH1D *n55=(TH1D *) MeeMass->ProjectionX("n55",ptb1,ptb2,"e");
	TH1D *n66=(TH1D *) MppMass->ProjectionX("n66",ptb1,ptb2,"e");


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//0-5, 1000 , 0.005
	//0.180  0.782 ,1.019, 3.096, 3.683
	//30     158     203   619     737 
  
	//const int Mpots=20;
	//int Mbin[Mpots+1]={0,10,20,30,50,70,90,120,160,200,240,280,330,390,460,540,620,700,780,860,1200};
	const int Mpots=34;
	double mass[Mpots];
	int Mbin[Mpots+1]={0,5,10,15,20,25,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240,260,280,300,330,360,390,440,500,580,660,740,860,1000,1200};
	double mass_Width[Mpots];

	const int Ppots=26;
	int Pbin[Ppots+1]={0,6,12,16,18,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,85,90,100,120,160,200};
	double Pt[Ppots];
	double Pt_Width[Ppots];

	// new 2D binning
	double massb[Mpots+1];
	double pb[Ppots+1];

	/////////////////////////////////////////
	const Int_t nBinM= MepMass->GetNbinsX();
	const Int_t nBinP= MepMass->GetNbinsY();
	cout<<"defult massbin: "<<nBinM<<" ptbin: "<<nBinP<<endl;
	/////////////////////////////////////////

	//bin center and bin width
	double binW=MepMass->GetXaxis()->GetBinWidth(18);    //mass bin
	double binWP=MepMass->GetYaxis()->GetBinWidth(18);   //pt bin 

	for(int i=0;i<Mpots;i++){	
		mass[i]=(Mbin[i]+Mbin[i+1])*binW/2.0;
		mass_Width[i]=(Mbin[i+1]-Mbin[i])*binW/2.0;   // afig0tully it's half binWidth
	}	

	for(int i=0;i<Ppots;i++){	
		Pt[i]=(Pbin[i]+Pbin[i+1])*binWP/2.0;
		Pt_Width[i]=(Pbin[i+1]-Pbin[i])*binWP/2.0;    // afig0tully it's half binWidth
	}


	//new bin 
	for(int i=0;i<Mpots+1;i++){	
		massb[i]=binW*Mbin[i];
			cout<<"mass :"<<massb[i]<<endl;
	}	

	for(int i=0;i<Ppots+1;i++){	
		pb[i]=binWP*Pbin[i];
			cout<<"pt: "<<pb[i]<<" :"<<i<<endl;
	}
	////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	TH1D *LikeBack=new TH1D("LikeBack","LikeBack",Mpots,massb);
	TH1D *MixBack1=new TH1D("MixBack1","MixBack1",Mpots,massb);
	TH1D *MixBack2=new TH1D("MixBack2","MixBack2",Mpots,massb);

	TH2D *LikeBack2D=new TH2D("LikeBack2D","LikeBack2D",Mpots,massb,Ppots,pb);
	TH2D *MixBack12D=new TH2D("MixBack12D","MixBack12D",Mpots,massb,Ppots,pb);
	TH2D *MixBack22D=new TH2D("MixBack22D","MixBack22D",Mpots,massb,Ppots,pb);


	double LikeV[Mpots],LikeV_err[Mpots];
	double Afig0[Mpots]; 
	double LikeV2[Mpots][Ppots],LikeV2_err[Mpots][Ppots];	

	LikeBack->Sumw2();
	MixBack1->Sumw2();
	MixBack2->Sumw2();

	LikeBack2D->Sumw2();
	MixBack12D->Sumw2();
	MixBack22D->Sumw2();

	//  1D
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	new11->Rebin(Mpots,"ah11",massb);
	new22->Rebin(Mpots,"ah22",massb);
	new33->Rebin(Mpots,"ah33",massb);
	new44->Rebin(Mpots,"ah44",massb);
	new55->Rebin(Mpots,"ah55",massb);
	new66->Rebin(Mpots,"ah66",massb);

	TH1D *h11=ah11->Clone();
	TH1D *h22=ah22->Clone();
	TH1D *h33=ah33->Clone();
	TH1D *h44=ah44->Clone();
	TH1D *h55=ah55->Clone();
	TH1D *h66=ah66->Clone();	

	if(mAcc==1){	
		for(int i=1;i<Mpots+1;i++){
			double mixpp=h55->GetBinContent(i);
			double mixee=h66->GetBinContent(i);	
			//if(mixpp==0.000){mixpp=mixee*0.5,mixee=mixee*0.5;}
			//if(mixee==0.000){mixee=mixpp*0.5,mixpp=mixpp*0.5;}			

			double afig0 = 1.;
			//if(mixpp==0.00&& mixee==0.00)afig0 = 0.00;
			if(mixpp==0.00 || mixee==0.00)afig0 = 0.00;
			else afig0 = h44->GetBinContent(i)/(pow(mixee*mixpp,0.5));

			//Afig0[i]=afig0;		

			double smpp=h22->GetBinContent(i);
			double smee=h33->GetBinContent(i);

			double verse=0.00;
			if(smee==0.00 || smpp==0.00)verse = 0.00;
			else verse=pow(smpp*smee,-0.5);

			LikeV[i-1] = pow(h22->GetBinContent(i)*h33->GetBinContent(i),0.5)*afig0;
			LikeV_err[i-1]= 0.5*verse*sqrt(pow(h22->GetBinContent(i)*h33->GetBinError(i),2.) + pow(h33->GetBinContent(i)*h22->GetBinError(i),2.))*afig0;

		}

		// set histo vaule
		for(int i=1;i<Mpots+1;i++){
			LikeBack->SetBinContent(i,LikeV[i-1]);
			LikeBack->SetBinError(i,LikeV_err[i-1]);		
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	


	// 2D
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TH2D *MassP1=new TH2D("MassP1","MassP1",Mpots,massb,Ppots,pb);
	TH2D *MassP2=new TH2D("MassP2","MassP2",Mpots,massb,Ppots,pb);
	TH2D *MassP3=new TH2D("MassP3","MassP3",Mpots,massb,Ppots,pb);
	TH2D *MassP4=new TH2D("MassP4","MassP4",Mpots,massb,Ppots,pb);
	TH2D *MassP5=new TH2D("MassP5","MassP5",Mpots,massb,Ppots,pb);
	TH2D *MassP6=new TH2D("MassP6","MassP6",Mpots,massb,Ppots,pb);

	MassP1->Sumw2();
	MassP2->Sumw2();
	MassP3->Sumw2();
	MassP4->Sumw2();
	MassP5->Sumw2();
	MassP6->Sumw2();

	if(mAcc==2){		
		for(int i=0;i<Mpots;i++){
			for(int j=0;j<Ppots;j++){
				double ber=0.,bct=0.;
				bct=h1_mp->IntegralAndError(Mbin[i]+1,Mbin[i+1],Pbin[j]+1,Pbin[j+1],ber);
				MassP1->SetBinContent(i+1,j+1,bct);
				MassP1->SetBinError(i+1,j+1,ber);

				bct=h2_mp->IntegralAndError(Mbin[i]+1,Mbin[i+1],Pbin[j]+1,Pbin[j+1],ber);
				MassP2->SetBinContent(i+1,j+1,bct);
				MassP2->SetBinError(i+1,j+1,ber);

				bct=h3_mp->IntegralAndError(Mbin[i]+1,Mbin[i+1],Pbin[j]+1,Pbin[j+1],ber);
				MassP3->SetBinContent(i+1,j+1,bct);
				MassP3->SetBinError(i+1,j+1,ber);

				bct=h4_mp->IntegralAndError(Mbin[i]+1,Mbin[i+1],Pbin[j]+1,Pbin[j+1],ber);
				MassP4->SetBinContent(i+1,j+1,bct);
				MassP4->SetBinError(i+1,j+1,ber);

				bct=h5_mp->IntegralAndError(Mbin[i]+1,Mbin[i+1],Pbin[j]+1,Pbin[j+1],ber);
				MassP5->SetBinContent(i+1,j+1,bct);
				MassP5->SetBinError(i+1,j+1,ber);

				bct=h6_mp->IntegralAndError(Mbin[i]+1,Mbin[i+1],Pbin[j]+1,Pbin[j+1],ber);
				MassP6->SetBinContent(i+1,j+1,bct);
				MassP6->SetBinError(i+1,j+1,ber);

			}
		}
		///////////////////////////////////////////
		for(int i=1;i<Mpots+1;i++){
			for(int j=1;j<Ppots+1;j++){
				double mixpp=MassP5->GetBinContent(i,j);
				double mixee=MassP6->GetBinContent(i,j);	
				//if(mixpp==0.000){mixpp=mixee*0.5,mixee=mixee*0.5;}
				//if(mixee==0.000){mixee=mixpp*0.5,mixpp=mixpp*0.5;}			

				double afig0 = 1.;
				//if(mixpp==0.00&& mixee==0.00)afig0 = 0.00;
				if(mixpp==0.00 || mixee==0.00)afig0 = 0.00;
				else afig0 = MassP4->GetBinContent(i,j)/(pow(mixee*mixpp,0.5));

				// if(afig0<0.000001)cout<<"ssss: "<<afig0<<" : "<<massb[i]<<" "<<pb[j]<<endl;

				double smpp=MassP2->GetBinContent(i,j);
				double smee=MassP3->GetBinContent(i,j);

				double verse=0.00;
				if(smee==0.00 || smpp==0.00)verse = 0.00;
				else verse=pow(smpp*smee,-0.5);

				LikeV2[i-1][j-1] = pow(MassP2->GetBinContent(i,j)*MassP3->GetBinContent(i,j),0.5)*afig0;
				LikeV2_err[i-1][j-1]= 0.5*verse*sqrt(pow(MassP2->GetBinContent(i,j)*MassP3->GetBinError(i,j),2.) + pow(MassP3->GetBinContent(i,j)*MassP2->GetBinError(i,j),2.))*afig0;


			}
		}

		for(int i=1;i<Mpots+1;i++){
			for(int j=1;j<Ppots+1;j++){		
				LikeBack2D->SetBinContent(i,j,LikeV2[i-1][j-1]);
				LikeBack2D->SetBinError(i,j,LikeV2_err[i-1][j-1]);	

			}  
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

	//  MixEvent back ground  
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	TH1D *H11 = (TH1D*) h11->Clone();
	TH1D *H22 = (TH1D*) h22->Clone();
	TH1D *H33 = (TH1D*) h33->Clone();
	TH1D *H44 = (TH1D*) h44->Clone();
	TH1D *H55 = (TH1D*) h55->Clone();
	TH1D *H66 = (TH1D*) h66->Clone();

	int LM  =new11->FindBin(lowM);
	int HM  =new11->FindBin(hiM);

	int LMLS=H11->FindBin(lowMls);

	double N11=n11->Integral(LM,HM);
	double N22=n22->Integral(LM,HM);
	double N33=n33->Integral(LM,HM);
	double N44=n44->Integral(LM,HM);
	double N55=n55->Integral(LM,HM);
	double N66=n66->Integral(LM,HM);	

	double T1 =n55->Integral();
	double T2 =n66->Integral();
	double T12=n44->Integral();
	double T1err =sqrt(T1);
	double T2err =sqrt(T2);
	double T12err=sqrt(T12);

	double defultRerr=2.*sqrt(pow(T1*T2err/(2.*sqrt(T1*T2)*T12),2.)+pow(T2*T1err/(2.*sqrt(T1*T2)*T12),2.)+pow(sqrt(T1*T2)*T12err/(T12*T12),2.) );
	double defultR=2.*sqrt(T1*T2)/T12;

	cout<<defultR<<endl;
	//defultR=1.;

	double A1=N22/N55;
	double A2=N33/N66;	

	double dA1=sqrt(N22/(N55*N55)+N22*N22/(N55*N55*N55));
	double dA2=sqrt(N33/(N66*N66)+N33*N33/(N66*N66*N66));

	double B1 =n55->Integral(); 
	double B2 =n66->Integral(); 
	double B12=n44->Integral();

	B1=B1*A1;
	B2=B2*A2;

	double norerr=0.5*sqrt(pow(dA1*A2,2.0)+pow(dA2*A1,2.0))/sqrt(A1*A2);

	double nor=sqrt(A1*A2);
	H44->Scale(nor);

	cout<<"(Nor)Likepp: "<<N22<<" Likemm: "<<N33<<" unLike: "<<N11<<" : "<<nor<<" : "<<norerr<<" : "<<norerr/nor<<endl;

	cout<<" MixEvent normalization Factor: "<<nor<<"     "<<lowM<<" - "<<hiM<<endl;

	H55->Scale(A1);
	H66->Scale(A2);

	///////////////////////////////////////////
	//double B1err =sqrt(B1);
	//double B2err =sqrt(B2);
	//double B12err=sqrt(B12);
	////double norerr=sqrt( pow(B1err*B2/B12,2.)+pow(B1*B2err/B12,2.)+pow(B1*B2*B12err/(B12*B12),2.) );  
	//double norerr=2.*sqrt(pow(B1*B2err/(2.*sqrt(B1*B2)*B12),2.)+pow(B2*B1err/(2.*sqrt(B1*B2)*B12),2.)+pow(sqrt(B1*B2)*B12err/(B12*B12),2.) );


	//double nor=2.0*sqrt(B1*B2)/(B12*defultR);
	//double norb=2.0*sqrt(B1b*B2b)/(B12*defultR);
	//H44->Scale(nor);

	//cout<<"Likepp: "<<B1<<" Likemm: "<<B2<<" unLike: "<<B12<<" : "<<nor<<" : "<<norerr<<" : "<<norerr/nor<<endl;

	//cout<<" MixEvent normalization Factor(1): "<<nor<<"     "<<lowM<<" - "<<hiM<<endl;
	//cout<<" MixEvent normalization Factor(2): "<<norb<<"     "<<lowMb<<" - "<<hiMb<<endl;

	//H55->Scale(A1);    
	//H66->Scale(A2);
	////////////////////////////////////////////// 

	//1D
	if(mAcc==1){
		MixBack1=(TH1D*)LikeBack->Clone();

	}


	//2D
	if(mAcc==2){
		MixBack12D=(TH2D*) MassP4->Clone();
		MixBack12D->Scale(nor);

		for(int i=1;i<Mpots+1;i++){
			for(int j=1;j<Ppots+1;j++){
				MixBack12D->SetBinContent(i,j,LikeV2[i-1][j-1]);
				MixBack12D->SetBinError(i,j,LikeV2_err[i-1][j-1]);		
			}	
		}	
	}	


	// Like + Mix 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(mAcc==1)TH1D *mLikeSign = (TH1D*) LikeBack->Clone();
	if(mAcc==2)TH1D *mLikeSign = (TH1D*) LikeBack2D->ProjectionX("mLikeSign");
	TH1D *mixBack = (TH1D*) H44->Clone();

	//LM =mixBack->FindBin(lowM);
	//HM =mixBack->FindBin(hiM);

	//double mAA=mLikeSign->Integral(LM,HM);
	//double mBB=mixBack->Integral(LM,HM);
	//mixBack->Scale(mAA/mBB);

	if(mAcc==1){
		MixBack2=(TH1D*)LikeBack->Clone();

		for(int i=LMLS;i<Mpots;i++){
			MixBack2->SetBinContent(i,mixBack->GetBinContent(i));   	
			MixBack2->SetBinError(i,mixBack->GetBinError(i)); 
		}
	}	

	if(mAcc==2){
		MixBack22D=(TH2D*) MassP1->Clone();
		MixBack22D->Scale(nor);

		for(int i=1;i<LMLS;i++){
			for(int j=1;j<Ppots+1;j++){
				MixBack22D->SetBinContent(i,j,LikeV2[i-1][j-1]);
				MixBack22D->SetBinError(i,j,LikeV2_err[i-1][j-1]);		
			}	
		}	
	}	

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	TH1D *signal1D=new TH1D("signal1D","signal1D",Mpots,massb);
	TH1D *BackGround=new TH1D("BackGround","BackGround",Mpots,massb);
	TH2D *signal2D=new TH2D("signal2D","signal2D",Mpots,massb,Ppots,pb);

	signal1D->Sumw2();	
	BackGround->Sumw2();	
	signal2D->Sumw2();	

	if(mback==1){	
		if(mAcc==1){
			signal1D->Add(h11,+1.);
			signal1D->Add(LikeBack,-1.);	
			BackGround=(TH1D*) LikeBack->Clone();
		}

		if(mAcc==2){
			signal2D->Add(MassP1,+1.);
			signal2D->Add(LikeBack2D,-1.);
			int ptbL=signal2D->GetYaxis().FindBin(pt1*1.00001);
			int ptbH=signal2D->GetYaxis().FindBin(pt2*1.00001)-1;
			signal1D=(TH1D*) signal2D->ProjectionX("signal1D_x",ptbL,ptbH);	
			BackGround=(TH1D*) MixBack12D->ProjectionX("BackGround_x",ptbL,ptbH);
		}

	}

	if(mback==2){
		if(mAcc==1){
			signal1D->Add(h11,+1.);
			signal1D->Add(MixBack1,-1.);	
			BackGround=(TH1D*) MixBack1->Clone();
		}

		if(mAcc==2){
			signal2D->Add(MassP1,+1.);
			signal2D->Add(MixBack12D,-1.);
			int ptbL=signal2D->GetYaxis().FindBin(pt1*1.00001);
			int ptbH=signal2D->GetYaxis().FindBin(pt2*1.00001)-1;
			signal1D=(TH1D*) signal2D->ProjectionX("signal1D_x",ptbL,ptbH);	
			BackGround=(TH1D*) MixBack12D->ProjectionX("BackGround_x",ptbL,ptbH);	
		}
	}	

	if(mback==3){
		if(mAcc==1){
			signal1D->Add(h11,+1.);
			signal1D->Add(MixBack2,-1.);	
			BackGround=(TH1D*) MixBack2->Clone();
		}

		if(mAcc==2){
			signal2D->Add(MassP1,+1.);
			signal2D->Add(MixBack22D,-1.);
			int ptbL=signal2D->GetYaxis().FindBin(pt1*1.00001);
			int ptbH=signal2D->GetYaxis().FindBin(pt2*1.00001)-1;
			signal1D=(TH1D*) signal2D->ProjectionX("signal1D_x",ptbL,ptbH);	
			BackGround=(TH1D*) MixBack22D->ProjectionX("BackGround_x",ptbL,ptbH);	
		}
	}




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//   Plots
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TCanvas *fig0 = new TCanvas("fig0", "fig0",0,0,800,600);
	fig0.cd();
	fig0->Divide(2,2,0.02,0,0);


	//TH1D HLikeM=(*H22)/(*H55);
	TH1D *HLikeM=H22->Clone();
	HLikeM->Divide(H22,H55,1.,1.,"");
	//TH1D HLikeP=(*H33)/(*H66);
	TH1D *HLikeP=H33->Clone();
	HLikeP->Divide(H33,H66,1.,1.,"");
	TH1D *HLikeS=H22->Clone();
	HLikeS->Add(H33);
	TH1D *HLikeSM=H55->Clone();
	HLikeSM->Add(H66);

	//TH1D HLikeT=(*HLikeS)/(*HLikeSM);
	TH1D *HLikeT=HLikeS->Clone();
	HLikeT->Divide(HLikeS,HLikeSM,1.,1.,"");
	//TH1D HMixRatio=(*H44)/(*HLikeSM);
	TH1D *HMixRatio=H44->Clone();
	HMixRatio->Divide(H44,HLikeSM,1.,1.,"");
	//TH1D HLikevsMix=(*BackGround)/(*H44);
	TH1D *HLikevsMix=BackGround->Clone();
	HLikevsMix->Divide(BackGround,H44,1.,1.,"");

	HLikeM.SetMarkerStyle(24);
	HLikeP.SetMarkerStyle(20);
	HMixRatio.SetMarkerStyle(20);	
	HLikeT.SetMarkerStyle(20);
	HLikevsMix.SetMarkerStyle(20);

	HLikeM.SetMarkerColor(2);
	HLikeP.SetMarkerColor(2);
	HMixRatio.SetMarkerColor(1);	
	HLikeT.SetMarkerColor(kViolet+5);
	HLikevsMix.SetMarkerColor(2);

	HLikeM.SetLineColor(2);
	HLikeP.SetLineColor(2);
	HMixRatio.SetLineColor(1);	
	HLikeT.SetLineColor(kViolet+5);
	HLikevsMix.SetLineColor(1);

	HLikeM.SetLineWidth(2);
	HLikeP.SetLineWidth(2);
	HMixRatio.SetLineWidth(2);	
	HLikeT.SetLineWidth(2);
	HLikevsMix.SetLineWidth(2);

	HLikeM.SetMarkerSize(1.0);
	HLikeP.SetMarkerSize(1.2);
	HMixRatio.SetMarkerSize(1.2);	
	HLikeT.SetMarkerSize(1.2);
	HLikevsMix.SetMarkerSize(1.);	

	double LX1=0.0;
	double HX1=3.0;
	double LY1=0.5;
	double HY1=2.5;	

	double LX2=0.0;
	double HX2=3.0;
	double LY2=0.5;
	double HY2=2.5;

	double mX1=LX1+0.1*(HX1-LX1); 
	double mY1=LY1+0.8*(HY1-LY1); 

	double mX2=LX2+0.1*(HX2-LX2); 
	double mY2=LY2+0.8*(HY2-LY2); 

	TLine *Line0 = new TLine(LX1,1.,HX1,1.);
	Line0->SetLineStyle(7);
	Line0->SetLineWidth(2);

	TLine *LineX1 = new TLine(LX1,LY1,LX1,HY1);
	LineX1->SetLineStyle(1);
	LineX1->SetLineWidth(3);

	TLine *LineX2 = new TLine(HX1,LY1,HX1,HY1);
	LineX2->SetLineStyle(1);
	LineX2->SetLineWidth(3);

	TLine *LineY1 = new TLine(LX1,LY1,HX1,LY1);
	LineY1->SetLineStyle(1);
	LineY1->SetLineWidth(3);

	TLine *LineY2 = new TLine(LX1,HY1,HX1,HY1);
	LineY2->SetLineStyle(1);
	LineY2->SetLineWidth(3);

	//1
	/////////////////////////////////////////
	fig0_1->cd();
	fig0_1->SetFrameLineWidth(2);
	TH2F *frame1 = new TH2F("frame1","",1,LX1,HX1,1,LY1,HY1);
	frame1->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frame1->GetYaxis()->SetTitle("Like/Mix_Like (--)");
	frame1->Draw();

	HLikeM.Draw("zpsame");
	Line0->Draw("same");

	TLatex *texx = new TLatex(mX1,mY1,"(a)  #frac{LikeS^{- -}(same)}{LikeS^{- -}(mix)} ");
	texx->SetTextSize(0.06);
	texx->SetLineWidth(2);
	texx->SetTextFont(42);
	texx->SetTextColor(kGray+2);
	texx->Draw();

	//2
	/////////////////////////////////////////
	fig0_2->cd();
	fig0_2->SetFrameLineWidth(2);
	TH2F *frame2 = new TH2F("frame2","",1,LX1,HX1,1,LY1,HY1);
	frame2->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frame2->GetYaxis()->SetTitle("Like/Mix_Like (++)");
	frame2->Draw();

	HLikeP.Draw("zpsame");
	Line0->Draw("same");

	TLatex *texx = new TLatex(mX1,mY1,"(b)  #frac{LikeS^{+ +}(same)}{LikeS^{+ +}(mix)} ");
	texx->SetTextSize(0.06);
	texx->SetLineWidth(2);
	texx->SetTextFont(42);
	texx->SetTextColor(kGray+2);
	texx->Draw();

	//3
	/////////////////////////////////////////
	fig0_3->cd();
	fig0_3->SetFrameLineWidth(2);
	TH2F *frame3 = new TH2F("frame3","",1,LX1,HX1,1,LY1,HY1);
	frame3->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frame3->GetYaxis()->SetTitle("Like/Mix_Like");
	frame3->Draw();

	HLikeT.Draw("zpsame");
	Line0->Draw("same");

	TLatex *texx = new TLatex(mX1,mY1,"(c)  #frac{LikeS (same)}{LikeS (mix)}  (++and--)");
	texx->SetTextSize(0.06);
	texx->SetLineWidth(2);
	texx->SetTextFont(42);
	texx->SetTextColor(kGray+2);
	texx->Draw();

	//4
	/////////////////////////////////////////
	//LY2=0.1;
	//HY2=1.9;
	//mY2=LY2+0.8*(HY2-LY2); 

	fig0_4->cd();
	fig0_4->SetFrameLineWidth(2);



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////
	//for(int i=0;i<nF;i++){
	//	Ftmass[i]=Ftmass[i];
	//	Ftmass_Width[i]=Ftmass_Width[i];
	//	FtSB_err[i]=FtSB[i];
	//	FtSB[i]=1.;
	//}  

	//FtHSB=new TGraphErrors(nF,Ftmass,FtSB,Ftmass_Width,FtSB_err);
	//FtHSB->SetMarkerStyle(20);
	////FtHSB->SetMarkerColor(kGray+2);
	//FtHSB->SetMarkerColor(kGreen-10);
	//FtHSB->SetFillStyle(1001);
	//FtHSB->SetFillColor(kGreen-10);
	//FtHSB->SetMarkerSize(1.2);
	//FtHSB->SetLineColor(0);	
	//FtHSB->SetLineWidth(2);	


	TCanvas *fig1 = new TCanvas("fig1", "fig1",0,0,800,600);	
	fig1->cd();

	TH2F *frame = new TH2F("frame","",1,LX1,HX1,1,LY1,HY1);
	frame->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frame->GetYaxis()->SetTitle("Like/Mix_unLike");
	frame->SetLineWidth(2);
	frame->Draw();

	TF1 *Fit1=new TF1("Fit1","pol0",0.4,2);
	Fit1->SetLineWidth(5);
	Fit1->SetLineColor(4);
	Fit1->SetRange(lowM,hiM);

	HLikevsMix.Fit(Fit1,"NOR");
	//Fit1->Draw("same");
	double nw=0.05;
	TBox *norRe=new TBox(lowM,Fit1->Eval(lowM)-nw,hiM,Fit1->Eval(lowM)+nw);
	norRe->SetFillColor(kGray+1);
	norRe->SetFillStyle(1001);
	norRe->SetLineColor(0);

	//FtHSB->Draw("e3same");
	norRe->Draw("same");

	//HLikevsMix.Draw("e1psame");
	//Line0->Draw("same");

	//double flowM=1.,fhiM=5.;
	double flowM=lowM,fhiM=3;
	TF1 *Fit2=new TF1("Fit2","1+exp((x-[0])/[1])",lowM,3);
	Fit2->SetParameter(0,100);
	Fit2->SetParameter(1,10);
	Fit2->SetRange(flowM,fhiM);

	HLikevsMix.Fit(Fit2,"NER");

	//Create a histogram to hold the confidence intervals
	TH1D *hint=(TH1D*) HLikevsMix.Clone();
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
	//Now the "hint" histogram has the fitted function values as the 
	//bin contents and the confidence intervals as bin errors
	//hint->SetStats(kFALSE);
	hint->SetFillColor(5);
	//hint->SetFillStyle(3002);
	//hint->Draw("e3 same");
	TH1D *hintL=hint.Clone();
	TH1D *hintH=hint.Clone();
	int nFbs=hintL->GetNbinsX();
	for(int i=0;i<nFbs;i++){
		double bc=hint.GetBinContent(i);
		double be=hint.GetBinError(i);
		hintL->SetBinContent(i,bc-be);
		hintH->SetBinContent(i,bc+be);
		hintL->SetBinError(i,be/10.);
		hintH->SetBinError(i,be/10.);
	}

	//hintL->Draw("same");
	//hintH->Draw("same");

	HLikevsMix.Draw("e1psame");

	Line0->Draw("same");
	Fit2->Draw("same");

	double fpar0 =Fit2->GetParameter(0);
	double fpar1 =Fit2->GetParameter(1);
	double fperr0=Fit2->GetParError(0);
	double fperr1=Fit2->GetParError(1);
	//TF1 *Fit2L=new TF1("Fit2L","1+exp((x-[0])/[1])",flowM,fhiM);
	//TF1 *Fit2H=new TF1("Fit2H","1+exp((x-[0])/[1])",flowM,fhiM);
	//Fit2L->SetParameter(0,fpar0);
	//Fit2L->SetParameter(1,fpar1-fperr1);
	//Fit2H->SetParameter(0,fpar0);
	//Fit2H->SetParameter(1,fpar1+fperr1);

	TF1 *Fit2L=new TF1("Fit2L","pol5",flowM,fhiM);
	TF1 *Fit2H=new TF1("Fit2H","pol5",flowM,fhiM);

	hintL->Fit(Fit2L,"NOR");
	hintH->Fit(Fit2H,"NOR");
	const int npars=6;
	double parL[npars];
	double parH[npars];
	Fit2L->GetParameters(parL);
	Fit2H->GetParameters(parH);

	Fit2 ->SetLineStyle(1);
	Fit2L->SetLineStyle(9);
	Fit2H->SetLineStyle(9);
	Fit2 ->SetLineWidth(2);
	Fit2L->SetLineWidth(2);
	Fit2H->SetLineWidth(2);
	Fit2 ->SetLineColor(4);
	Fit2L->SetLineColor(6);
	Fit2H->SetLineColor(6);

	Fit2L->Draw("same");
	Fit2H->Draw("same");

	ofstream outData;
	outData.open("LikeS_star_Nor.h");

	outData<<"const int npars="<<npars<<";"<<endl;;
	for(int i=0;i<npars;i++){
		if(i==0)outData<<"double parL[npars]={";
		outData<<parL[i]<<", ";
		if(i==(npars-1))outData<<" }; "<<endl;
	}
	for(int i=0;i<npars;i++){
		if(i==0)outData<<"double parH[npars]={";
		outData<<parH[i]<<", ";
		if(i==(npars-1))outData<<" }; "<<endl;
	}

	outData<<"double flowM="<<flowM<<";"<<endl;;
	outData<<"double fhiM= "<<fhiM <<";"<<endl;;
	outData<<"double fpar0="<<fpar0<<";"<<endl;;
	outData<<"double fpar1="<<fpar1<<";"<<endl;;
	outData<<"TF1 *fFit2 =new TF1(\"fFit2\", \"1+exp((x-[0])/[1])\",flowM,fhiM);"<<endl;
	outData<<"TF1 *fFit2L=new TF1(\"fFit2L\",\"pol5\",flowM,fhiM);"<<endl;
	outData<<"TF1 *fFit2H=new TF1(\"fFit2H\",\"pol5\",flowM,fhiM);"<<endl;
	outData<<"fFit2 ->SetParameter(0,fpar0);"<<endl;
	outData<<"fFit2 ->SetParameter(1,fpar1);"<<endl;
	outData<<"fFit2L->SetParameters("<<parL[0]<<","<<parL[1]<<","<<parL[2]<<","<<parL[3]<<","<<parL[4]<<","<<parL[5]<<");"<<endl;
	outData<<"fFit2H->SetParameters("<<parH[0]<<","<<parH[1]<<","<<parH[2]<<","<<parH[3]<<","<<parH[4]<<","<<parH[5]<<");"<<endl;

	//outData<<"double fpar0="<<fpar0<<";"<<endl;;
	//outData<<"double fpar1="<<fpar1<<";"<<endl;;
	//outData<<"double fperr0="<<fperr0<<";"<<endl;;
	//outData<<"double fperr1="<<fperr1<<";"<<endl;;
	//outData<<"double flowM="<<flowM<<";"<<endl;;
	//outData<<"double fhiM= "<<fhiM <<";"<<endl;;
	//outData<<"TF1 *fFit2 =new TF1(\"fFit2\", \"1+exp((x-[0])/[1])\",flowM,fhiM);"<<endl;
	//outData<<"TF1 *fFit2L=new TF1(\"fFit2L\",\"1+exp((x-[0])/[1])\",flowM,fhiM);"<<endl;
	//outData<<"TF1 *fFit2H=new TF1(\"fFit2H\",\"1+exp((x-[0])/[1])\",flowM,fhiM);"<<endl;
	//outData<<"fFit2 ->SetParameter(0,fpar0);"<<endl;
	//outData<<"fFit2 ->SetParameter(1,fpar1);"<<endl;
	//outData<<"fFit2L->SetParameter(0,fpar0);"<<endl;
	//outData<<"fFit2L->SetParameter(1,fpar1-fperr1);"<<endl;
	//outData<<"fFit2H->SetParameter(0,fpar0);"<<endl;
	//outData<<"fFit2H->SetParameter(1,fpar1+fperr1);"<<endl;

	outData.close();

	//return;

	LineX1->Draw();
	LineX2->Draw();
	LineY1->Draw();
	LineY2->Draw();

	TLatex *tex100 = new TLatex(2.5,1.04," 100% error ");
	tex100->SetTextSize(0.04);
	tex100->SetLineWidth(2);
	tex100->SetTextFont(42);
	tex100->SetTextColor(kGreen-8);
	tex100->Draw();


	TLatex *texx = new TLatex(mX1,mY1," #frac{LikeS (same)}{unLikeS (mix)} ");
	texx->SetTextSize(0.05);
	texx->SetLineWidth(2);
	texx->SetTextFont(42);
	texx->SetTextColor(kGray+2);
	texx->Draw();


	char Tname[256];
	sprintf(Tname,"p + p  #sqrt{s_{NN}} = 200 GeV (%s)",TitleC);

	Title = new TLegend(0.38,0.2,0.9,0.3);
	Title->SetFillColor(10);
	Title->SetFillStyle(0);
	Title->SetBorderSize(0);
	Title->SetTextSize(0.05);
	Title->SetTextFont(22);
	Title->SetTextColor(kGray+2);
	Title->SetHeader(Tname); 
	//Title->Draw();

	fig1->SaveAs("ratio_background.gif");
	fig1->SaveAs("ratio_background.eps");

  fig0->cd();
  fig0_4->cd();
	TH2F *frame4 = new TH2F("frame4","",1,LX2,HX2,1,LY2,HY2+2);
	//TH2F *frame4 = new TH2F("frame4","",1,LX2,HX2,1,0.1,1.9);
	frame4->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frame4->GetYaxis()->SetTitle("Mix_unLike/Mix_Like");
	frame4->Draw();

	//HMixRatio.Draw("zpsame");

	//cout<<HMixRatio.GetBinCenter(20)<<" "<<HMixRatio.GetBinWidth(20)<<"  "<<HMixRatio.GetBinContent(20)<<" "<<HMixRatio.GetBinError(20)<<endl;
	//cout<<HMixRatio.GetBinCenter(28)<<" "<<HMixRatio.GetBinWidth(28)<<"  "<<HMixRatio.GetBinContent(28)<<" "<<HMixRatio.GetBinError(28)<<endl;

	Line0->Draw("same");

	TLatex *texx = new TLatex(mX2,mY2,"(d)  #frac{LikeS (same)}{unLikeS (mix)} ");
	texx->SetTextSize(0.06);
	texx->SetLineWidth(2);
	texx->SetTextFont(42);
	texx->SetTextColor(kGray+2);
	texx->Draw();

	norRe->Draw("same");
	HLikevsMix.Draw("e1psame");
	Line0->Draw("same");
	Fit2->Draw("same");
	Fit2L->Draw("same");
	Fit2H->Draw("same");

	fig0->SaveAs("ratio.gif");
	fig0->SaveAs("ratio.eps");


	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////	
	fig2=new TCanvas("fig2","fig2",850,600);
	fig2->cd();


	TH2F *framed=new TH2F("framed","",1,LX1,HX1,1,0.85,1.15);
	framed->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	framed->GetYaxis()->SetTitle("++/--");
	framed->Draw();

	//TH1D rPPMM=(*H33)/(*H22);
	TH1F *rPPMM=H33->Clone();
	rPPMM->Divide(H33,H22,1.,1.,"");
	//TH1D mPPMM=(*H66)/(*H55);
	TH1F *mPPMM=H66->Clone();
	mPPMM->Divide(H66,H55,1.,1.,"");

	rPPMM->SetMarkerSize(1.5);
	rPPMM->SetMarkerStyle(20);
	rPPMM->SetMarkerColor(1);
	mPPMM->SetMarkerSize(1.5);
	mPPMM->SetMarkerStyle(20);
	mPPMM->SetMarkerColor(2);
	rPPMM->SetLineColor(1);
	mPPMM->SetLineColor(2);
	rPPMM->SetLineWidth(2);
	mPPMM->SetLineWidth(2);

	rPPMM->Draw("zpsame");
	mPPMM->Draw("zpsame");
	Line0->Draw("same");

	TH1F *LA=(TH1F*) rPPMM.Clone();
	TH1F *LB=(TH1F*) mPPMM.Clone();

	legd = new TLegend(0.28,0.65,0.7,0.85);
	legd->SetFillColor(10);
	legd->SetFillStyle(0);
	legd->SetBorderSize(0);
	legd->SetTextSize(0.05);
	legd->SetTextFont(22);
	legd->SetTextColor(kGray+2);
	legd->SetHeader(Tname); 
	legd->AddEntry(LA,"sameEvent","lp");
	legd->AddEntry(LB,"mixEvent","lp");
	legd->Draw();

	fig2->SaveAs("ppmm.gif");
	fig2->SaveAs("ppmm.eps");



	////////////////////////////////////////////////////////////////////////////////	
	////////////////////////////////////////////////////////////////////////////////
	fig3=new TCanvas("fig3","fig3",850,600);
	fig3->cd();

	TH2F *framed6=new TH2F("framed6","",1,LX1,HX1,1,0.8,1.1);
	framed6->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	framed6->GetYaxis()->SetTitle("sum/geometry");
	framed6->Draw();

	TH1F *AAA=(TH1F *) H22->Clone();
	TH1F *BBB=(TH1F *) H33->Clone();
	TH1F *CCC=(TH1F *) H55->Clone();
	TH1F *DDD=(TH1F *) H66->Clone();

	TH1F *Rsunppmm=AAA.Clone();
	Rsunppmm->Add(AAA,BBB,1.,1.);

	TH1F *Msunppmm=CCC.Clone();
	Msunppmm->Add(CCC,DDD,1.,1.);

	TH1F *Rgmppmm=(TH1F *) H22->Clone();
	TH1F *Mgmppmm=(TH1F *) H22->Clone();


	for(int i=0;i<Mpots;i++){

		float vb=2.*sqrt(AAA->GetBinContent(i)*BBB->GetBinContent(i));
		float vb_err;
		if(vb!=0.)vb_err=(1./(vb))*(sqrt( pow(AAA->GetBinContent(i)*BBB->GetBinError(i),2.)+ pow(BBB->GetBinContent(i)*AAA->GetBinError(i),2.) ));
		else vb_err=0.;

		Rgmppmm->SetBinContent(i,vb);
		Rgmppmm->SetBinError(i,vb_err);


		vb=2.*sqrt(CCC->GetBinContent(i)*DDD->GetBinContent(i));
		if(vb!=0.)vb_err=(1./(vb))*(sqrt( pow(CCC->GetBinContent(i)*DDD->GetBinError(i),2.)+ pow(DDD->GetBinContent(i)*CCC->GetBinError(i),2.) ));
		else vb_err=0.;

		Mgmppmm->SetBinContent(i,vb);
		Mgmppmm->SetBinError(i,vb_err);

	}	

	Rsunppmm.Scale(Rgmppmm->Integral()/Rsunppmm.Integral());
	Msunppmm.Scale(Mgmppmm->Integral()/Msunppmm.Integral());


	//TH1F RLike=(*Rsunppmm)/(*Rgmppmm);
	TH1F *RLike=Rsunppmm->Clone();
	RLike->Divide(Rsunppmm,Rgmppmm,1.,1.,"");


	TH1F *MLike=Msunppmm->Clone();
	MLike->Divide(Msunppmm,Mgmppmm,1.,1.,"");

	RLike.SetMarkerStyle(20);
	RLike.SetMarkerSize(1.2);
	RLike.SetMarkerColor(1);

	MLike.SetMarkerStyle(24);
	MLike.SetMarkerSize(1.5);
	MLike.SetMarkerColor(2);
	MLike.SetLineWidth(3);
	MLike.SetLineColor(2);

	RLike.Draw("epsame");
	MLike.Draw("epsame");
	RLike.SetName("RLike");
	MLike.SetName("MLike");

	legd6 = new TLegend(0.27,0.65,0.7,0.85);
	legd6->SetFillColor(10);
	legd6->SetFillStyle(0);
	legd6->SetBorderSize(0);
	legd6->SetTextSize(0.05);
	legd6->SetTextFont(22);
	legd6->SetTextColor(kGray+2);
	legd6->SetHeader(Tname); 
	legd6->AddEntry("RLike","sameEvent","lp");
	legd6->AddEntry("MLike","mixEvent","lp");
	legd6->Draw();

	Line0->Draw("same");

	fig3->SaveAs("geometry_ratio.gif");
	fig3->SaveAs("geometry_ratio.eps");


	//////////////////////////////////////////////////////////////////	
	//////////////////////////////////////////////////////////////////	
	fig0->cd();	
	fig0_4->cd();

	TH1F *Afig0ept=(TH1F *) H44->Clone();

	//TH1F gmAcc=(*Afig0ept)/(*Mgmppmm);
	TH1F *gmAcc=Afig0ept->Clone();
	gmAcc->Divide(Afig0ept,Mgmppmm,1.,1.,"");

	gmAcc.SetMarkerStyle(24);
	gmAcc.SetMarkerSize(1.2);
	gmAcc.SetMarkerColor(2);
	gmAcc.SetLineWidth(3);
	gmAcc.SetLineColor(2);

	//gmAcc.Draw("epsame");
	gmAcc.SetName("gmAcc");

	TH1D *bHMixRatio=HMixRatio.Clone();

	//fig0->Update();	
	//legd7 = new TLegend(0.4,0.15,0.85,0.45);
	//legd7->SetFillColor(10);
	//legd7->SetFillStyle(0);
	//legd7->SetBorderSize(0);
	//legd7->SetTextSize(0.05);
	//legd7->SetTextFont(22);
	//legd7->SetTextColor(kGray+2);
	//legd7->SetHeader(Tname); 
	//legd7->AddEntry(bHMixRatio,"sum ++ and --","p");
	//legd7->AddEntry("gmAcc","geometry mean","p");
	//legd7->Draw();

	//fig0->SaveAs("ratio.gif");
	//fig0->SaveAs("ratio.eps");
	//fig0_4->SaveAs("Acc.gif");
	//fig0_4->SaveAs("Acc.eps");


	////////////////////////////////////////////

	TH2D *Rep2D=h1_mp.Clone();
	TH2D *Rpp2D=h2_mp.Clone();
	TH2D *Rmm2D=h3_mp.Clone();
	TH2D *Mpp2D=h5_mp.Clone();
	TH2D *Mmm2D=h6_mp.Clone();

	Rep2D->Sumw2();
	Rpp2D->Sumw2();
	Rmm2D->Sumw2();
	Mpp2D->Sumw2();
	Mmm2D->Sumw2();

	Rep2D->Rebin2D(20,5);
	Rpp2D->Rebin2D(20,5);
	Rmm2D->Rebin2D(20,5);
	Mpp2D->Rebin2D(20,5);
	Mmm2D->Rebin2D(20,5);

	Mpp2D->Scale(A1);	
	Mmm2D->Scale(A2);
	Rpp2D->Add(Rmm2D);
	Mpp2D->Add(Mmm2D);

	TH2D *BACK=(TH2D*)Rpp2D.Clone();
	BACK->Add(Mpp2D,-1.);
	TH2D *R2d=BACK->Clone();	

  TCanvas *ctmp = new TCanvas("ctmp","",800,600);
  ctmp->SetFrameLineWidth(3);
  ctmp->SetLogy();
  ctmp->cd();
  int pimass = Rep2D->GetXaxis()->FindBin(0.2);
  TH1D *Rep1D_pT = (TH1D *)Rep2D->ProjectionY("Rep1D_pT",pimass,-1);
  Rep1D_pT->SetMinimum(1e-1);
  Rep1D_pT->GetXaxis()->SetTitle("p_{T}^{ee} GeV/c");
  Rep1D_pT->SetLineColor(4);
  Rep1D_pT->SetLineWidth(3);
  Rep1D_pT->SetLineStyle(7);
  Rep1D_pT->Draw("HIST");
  TH1D *Rpp1D_pT = (TH1D *)Rpp2D->ProjectionY("Rpp1D_pT",pimass,-1);
  Rpp1D_pT->SetLineWidth(3);
  Rpp1D_pT->Draw("same");
  TH1D *Mpp1D_pT = (TH1D *)Mpp2D->ProjectionY("Mpp1D_pT",pimass,-1);
  Mpp1D_pT->SetLineWidth(3);
  Mpp1D_pT->SetLineColor(2);
  Mpp1D_pT->Draw("same");
  

  TLegend *tctmp = new TLegend(0.6,0.55,0.8,0.8);
  tctmp->SetHeader("M_{ee}>0.2GeV/c^{2}");
  tctmp->SetFillColor(10);
  tctmp->SetBorderSize(0);
  tctmp->SetTextSize(0.045);
  tctmp->AddEntry(Rpp1D_pT,"Real event likesign");
  tctmp->AddEntry(Mpp1D_pT,"Mixed event likesign");
  tctmp->AddEntry(Rep1D_pT,"Real event unlikesign");
  tctmp->Draw("same");

  ctmp->SaveAs("nor_pT.gif");

  TCanvas *ctmp2 = new TCanvas("ctmp2","",800,600);
  ctmp2->SetFrameLineWidth(3);
  ctmp2->SetLogy();
  ctmp2->cd();

  TH1D *Rep1D_m = (TH1D *)Rep2D->ProjectionX("Rep1D_m");
  Rep1D_m->SetMinimum(1e-1);
  Rep1D_m->GetXaxis()->SetTitle("M_{ee} GeV/c");
  Rep1D_m->SetLineColor(4);
  Rep1D_m->SetLineWidth(3);
  Rep1D_m->SetLineStyle(7);
  Rep1D_m->Draw("HIST");
  TH1D *Rpp1D_m = (TH1D *)Rpp2D->ProjectionX("Rpp1D_m");
  Rpp1D_m->SetLineWidth(3);
  Rpp1D_m->Draw("same");
  TH1D *Mpp1D_m = (TH1D *)Mpp2D->ProjectionX("Mpp1D_m");
  Mpp1D_m->SetLineWidth(3);
  Mpp1D_m->SetLineColor(2);
  Mpp1D_m->Draw("same");
  

  TLegend *tctmp2 = new TLegend(0.6,0.7,0.8,0.8);
  tctmp2->SetFillColor(10);
  tctmp2->SetBorderSize(0);
  tctmp2->SetTextSize(0.045);
  tctmp2->AddEntry(Rpp1D_m,"Real event likesign");
  tctmp2->AddEntry(Mpp1D_m,"Mixed event likesign");
  tctmp2->AddEntry(Rep1D_m,"Real event unlikesign");
  tctmp2->Draw("same");

  ctmp2->SaveAs("nor_m.gif");


	const int nX=R2d->GetNbinsX();
	const int nY=R2d->GetYaxis()->FindBin(pt2-1e-8);;
	double sigma=1.;
	double mean=R2d->Integral(1,-1,1,nY)/(nX*nY);
	double sum2=0.;
	TH1D *SigmaD0=new TH1D("SigmaD0","#sigma(N_{#pm#pm}-B_{#pm#pm})0",150,-20,20);
	TH1D *SigmaD1=new TH1D("SigmaD1","#sigma(N_{#pm#pm}-B_{#pm#pm})1",150,-20,20);
	TH1D *SigmaD2=new TH1D("SigmaD2","#sigma(N_{#pm#pm}-B_{#pm#pm})2",150,-20,20);

	for(int i=0;i<nX;i++){
		for(int j=0;j<nY;j++){
			if(R2d->GetBinError(i+1,j+1)==0.)continue;
			double xmass=R2d->GetXaxis()->GetBinCenter(i+1);
			if(xmass<m2)SigmaD0->Fill(R2d->GetBinContent(i+1,j+1)/R2d->GetBinError(i+1,j+1));
			if(xmass<1.0)SigmaD1->Fill(R2d->GetBinContent(i+1,j+1)/R2d->GetBinError(i+1,j+1));
			if(xmass>=1.&&xmass<2.0)SigmaD2->Fill(R2d->GetBinContent(i+1,j+1)/R2d->GetBinError(i+1,j+1));
		}
	}


	for(int i=0;i<nX;i++){
		for(int j=0;j<nY;j++){
      double Rpp = Rpp2D->GetBinContent(i+1,j+1);
      double Mpp = Mpp2D->GetBinContent(i+1,j+1);
      double dRpp = Rpp2D->GetBinError(i+1,j+1);
      double dMpp = Mpp2D->GetBinError(i+1,j+1);
			if(Rpp==0) R2d->SetBinContent(i+1,j+1,-111);
      else R2d->SetBinContent(i+1,j+1,(Rpp-Mpp)/sqrt(dRpp*dRpp+dMpp*dMpp));
      //R2d->SetBinContent(i+1,j+1,(Rpp-Mpp)/sqrt(dRpp*dRpp+dMpp*dMpp));
    }
  }

  sigma=sqrt(sum2/(nX*nY));
	//R2d->Scale(1./sigma);
	cout<<sigma<<endl;


	//TH2D R2d=(*Rpp2D)/(*Mpp2D);	
	TCanvas fign("fign","fign",800,600);
	fign.cd();

	TPad *pl2=new TPad("pl2","pl2",0.0,0.0,1.0,1.0,0,0,0);
	pl2->SetLeftMargin(0.13);
	pl2->SetRightMargin(0.1);
	pl2->SetTopMargin(0.06);
	pl2->SetBottomMargin(.12);
	//pl2->SetTickx(1);
	//pl2->SetTicky(1);
	//pl2->SetGrid(1,1);
	pl2->Draw();
	pl2->cd();
	//pl2->SetLogy();

	TH2F *framen=new TH2F("framen","",1,0.,3,1,pt1,pt2);
	framen->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	framen->GetYaxis()->SetTitle("p_{T}^{ee} (GeV/c) ");
	framen->GetZaxis()->SetTitle("N_{#pm#pm}-B_{#pm#pm}");
	framen->Draw();

	gStyle->SetNumberContours(20);
	gStyle->SetPalette(1);
	R2d.GetZaxis()->SetRangeUser(-10,10);
	R2d.Draw("colzsame");
	TBox NR(lowM,0.,hiM,2.1);
	NR.SetFillColor(0);
	NR.SetFillStyle(0);
	NR.SetLineColor(1);
	NR.SetLineWidth(3);
	NR.Draw();


	legn = new TLegend(0.33,0.7,0.9,0.88);
	legn->SetFillColor(10);
	legn->SetFillStyle(0);
	legn->SetBorderSize(0);
	legn->SetTextSize(0.05);
	legn->SetTextFont(62);
	legn->SetTextColor(2);
	//legn->SetTextColor(kGray+2);
	legn->SetHeader(Tname); 
	legn->AddEntry("","        (N_{#pm#pm}-B_{#pm#pm})/#sigma(N_{#pm#pm}-B_{#pm#pm})","");
	legn->Draw();

	fign.SaveAs("norRegion.gif");
	fign.SaveAs("paper_nRegion.gif");
	fign.SaveAs("paper_nRegion.eps");

	TCanvas figS("figS","figS",1200,450);
	figS.cd();
	figS.Divide(3,1,0,0,0);

	figS.cd(1);
	SigmaD0->Draw();	
	SigmaD0->SetLineColor(2);	
	SigmaD0->Fit("gaus","","",-10,10);	
	SigmaD0->GetYaxis()->SetTitle("counts");	

	TLegend Ld1(0.18,0.3,0.58,0.4);
	Ld1->SetBorderSize(0);
	Ld1->SetFillStyle(0);
	Ld1->SetTextColor(2);
	Ld1->SetHeader("0-3 (GeV/c^{2})");
	Ld1->Draw();

	figS.cd(2);
	SigmaD1->Draw();	
	SigmaD1->SetLineColor(2);	
	SigmaD1->Fit("gaus","","",-10,10);	
	SigmaD1->GetXaxis()->SetTitle("(N_{#pm#pm}-B_{#pm#pm})/#sigma(N_{#pm#pm}-B_{#pm#pm})");	
	SigmaD1->GetXaxis()->CenterTitle();	

	TLegend Ld2(0.05,0.3,0.55,0.4);
	Ld2->SetBorderSize(0);
	Ld2->SetFillStyle(0);
	Ld2->SetTextColor(2);
	Ld2->SetHeader("0-1 (GeV/c^{2})");
	Ld2->Draw();

	figS.cd(3);
	SigmaD2->SetLineColor(2);	
	SigmaD2->Fit("gaus","","",-10,10);	
	SigmaD2->Draw();	

	TLegend Ld3(0.05,0.3,0.55,0.4);
	Ld3->SetBorderSize(0);
	Ld3->SetFillStyle(0);
	Ld3->SetTextColor(2);
	Ld3->SetHeader("1-2 (GeV/c^{2})");
	Ld3->Draw();

	figS->SaveAs("sigma.gif");
	figS->SaveAs("sigma.eps");

	//paper
	///////////////////////////////////////////////////
	TCanvas *figF = new TCanvas("figF", "figF",0,0,800,600);	
  figF->SetFrameLineWidth(3);
	figF->cd();

	//LY2=0.1;
	//HY2=1.9;
	//TH2F *frameF = new TH2F("frameF","",1,LX2,HX2,1,0.1,1.9);
	TH2F *frameF = new TH2F("frameF","",1,LX2,HX2,1,0.9,1.1);
	frameF->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frameF->GetYaxis()->SetTitle("acceptence factor");
  frameF->GetYaxis()->SetTitleOffset(1.5);
	frameF->Draw();

	TH1D *hACC=(TH1D*)HMixRatio.Clone();
	hACC.SetMarkerColor(2);
	hACC.SetMarkerSize(1.8);
	hACC.Draw("zpsame");
	Line0.Draw("same");

	mY2=LY2+0.8*(HY2-LY2); 

	TLatex *texx = new TLatex((LX2+HX2)*0.35,1.08," #frac{unLikeS (mix)}{LikeS (mix)} ");
	texx->SetTextSize(0.05);
	texx->SetLineWidth(2);
	texx->SetTextFont(42);
	texx->SetTextColor(kGray+2);
	texx->Draw();

	Title = new TLegend(0.38,0.2,0.9,0.3);
	Title->SetFillColor(10);
	Title->SetFillStyle(0);
	Title->SetBorderSize(0);
	Title->SetTextSize(0.05);
	Title->SetTextFont(62);
	Title->SetTextColor(1);
	Title->SetHeader(Tname); 
	Title->Draw();

	figF->SaveAs("paper_Acc.gif");
	figF->SaveAs("paper_Acc.eps");

	//
	TCanvas *figF1 = new TCanvas("figF1", "figF1",0,0,800,600);	
	figF1->cd();
	figF1->SetTicks(1,1);
	//figF1->SetGrid(1,1);

	LY1=0.94;
	HY1=1.06;
	//LY1=0.96;
	//HY1=1.04;
	TH2F *frameF1 = new TH2F("frameF1","",1,LX1,HX1,1,LY1,HY1);
	frameF1->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
	frameF1->GetYaxis()->SetTitle("Like/Mix_unLike");
	frameF1->GetYaxis()->SetNdivisions(510);
	frameF1->GetXaxis()->SetNdivisions(510);
	frameF1->SetLineWidth(2);
	frameF1->Draw();

	//FtHSB->Draw("e3same");
	norRe->Draw("same");

	HLikevsMix.SetMarkerSize(1.2);
	HLikevsMix.Draw("e1psame");
	Line0->Draw("same");

	HLikevsMix.Fit(Fit2,"NOR");
	Fit2->Draw("same");
	Fit2L->Draw("same");
	Fit2H->Draw("same");

	//FtHSB->Draw("e3same");

	//HLikevsMix.SetMarkerSize(1.2);
	//HLikevsMix.Draw("ze1psame");
	//Line0->Draw("same");

	//TF1 *Fit1=new TF1("Fit1","pol0",0.75,3.0);
	//Fit1->SetLineWidth(3);
	//Fit1->SetLineColor(4);
	//Fit1->SetRange(lowM,hiM);

	//HLikevsMix.Fit(Fit1,"NOR");
	//Fit1->Draw("same");

	TLine *Line0 = new TLine(LX1,1.,HX1,1.);
	Line0->SetLineStyle(7);
	Line0->SetLineWidth(2);

	TLine *LineX1 = new TLine(LX1,LY1,LX1,HY1);
	LineX1->SetLineStyle(1);
	LineX1->SetLineWidth(3);

	TLine *LineX2 = new TLine(HX1,LY1,HX1,HY1);
	LineX2->SetLineStyle(1);
	LineX2->SetLineWidth(3);

	TLine *LineY1 = new TLine(LX1,LY1,HX1,LY1);
	LineY1->SetLineStyle(1);
	LineY1->SetLineWidth(3);

	TLine *LineY2 = new TLine(LX1,HY1,HX1,HY1);
	LineY2->SetLineStyle(1);
	LineY2->SetLineWidth(3);

	LineX1->Draw();
	LineX2->Draw();
	LineY1->Draw();
	LineY2->Draw();

	mX1=LX1+0.1*(HX1-LX1); 
	mY1=LY1+0.77*(HY1-LY1); 

	TLatex *texx = new TLatex(mX1,mY1," #frac{LikeS (same)}{unLikeS (mix)} ");
	texx->SetTextSize(0.05);
	texx->SetLineWidth(2);
	texx->SetTextFont(42);
	texx->SetTextColor(kGray+2);
	texx->Draw();

	mT = new TLegend(0.5,0.16,0.8,0.21);
	mT->SetFillColor(10);
	mT->SetFillStyle(0);
	mT->SetBorderSize(0);
	mT->SetTextSize(0.04);
	mT->SetTextFont(42);
	mT->SetTextColor(kGreen-8);
	//mT->AddEntry(FtHSB,"100% error estimated by S/B","F");
	mT->Draw();

	//	Title = new TLegend(0.38,0.15,0.9,0.25);
	Title = new TLegend(0.15,0.8,0.9,0.92);
	Title->SetFillColor(10);
	Title->SetFillStyle(0);
	Title->SetBorderSize(0);
	Title->SetTextSize(0.05);
	Title->SetTextFont(62);
	Title->SetTextColor(1);
	Title->SetHeader(Tname); 
	Title->Draw();

	Title2 = new TLegend(0.14,0.16,0.45,0.31);
	Title2->SetFillColor(10);
	Title2->SetFillStyle(0);
	Title2->SetBorderSize(0);
	Title2->SetTextSize(0.035);
	Title2->SetTextFont(42);
	Title2->SetTextColor(kGray+1);
	Title2->AddEntry(Fit2,"Fitting","l" ); 
	Title2->AddEntry(Fit2H,"95\% confidence level","l"); 
	Title2->AddEntry(norRe,"normalization region","F"); 
	Title2->Draw();

	figF1->SaveAs("paper_ratio.gif");
	figF1->SaveAs("paper_ratio.eps");


	TFile fout("acc_QA.root","recreate");
	fout.cd();
	hACC->SetName("AcceptanceC");
	hACC->Write();
	HLikeT->SetName("LikeS_Ratio");
	HLikeT->Write();
	HLikevsMix->SetName("back_Ratio");
	HLikevsMix->Write();
  signal1D->SetName("signal");
  signal1D->Write();
  BackGround->SetName("background");
  BackGround->Write();

	fout.Write();

	return;
}







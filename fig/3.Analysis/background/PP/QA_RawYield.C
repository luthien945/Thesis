{
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "LikeS_star_Nor.h"
//#include "../sb.h"

	////////////////////
	//  backGround
	const Int_t mback=3;    // 1: LikeSign, 2: MixEvent, 3: LikeSign Normalized MixEvent
	const Int_t mAcc=2;     // 1: 1D afig0eptance correction  2: 2D 

	const double pt1=0.;
	const double pt2=2.;

	const Int_t cent1=2;
	const Int_t cent2=10;
	///////////////////
	// normalization range
	//double lowM=1.0;  
	//double hiM=2.00;
	double lowM=0.4;  
	double hiM=1.5;
	double lowpT=0;  
	double hipT=2.;
	double lowMls=0.4;  // Low limit for likeS background  

	char TitleC[256]="MinBias";

  fFit2 ->SetParameter(0,fpar0);
  fFit2 ->SetParameter(1,fpar1);
  fFit2L->SetParameters(1.01111,-0.0695908,0.148037,-0.135587,0.0561521,-0.00257007);
  fFit2H->SetParameters(0.888892,0.735515,-1.6236,1.76809,-0.884687,0.17969);
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
	double nevent=Hcentrality->Integral();	
	cout<<Hcentrality->GetEntries()<<endl;	
	cout<<Hcentrality->Integral()<<endl;	

	TH2D *RepMass = (TH2D *)f->Get("MassVsPtVsCen_npe_unLikeSign");               //  x,mass;  y,pt;   z,cent
	TH2D *ReeMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignNN");      
	TH2D *RppMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignPP");     	
	TH2D *MepMass = (TH2D *)f->Get("MassVsPtVsCen_npe_unLikesign_Mix"); 
	TH2D *MeeMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignNN_Mix"); 
	TH2D *MppMass = (TH2D *)f->Get("MassVsPtVsCen_npe_LikeSignPP_Mix"); 		

	//RepMass->GetZaxis()->SetRange(cent1,cent2);
	//ReeMass->GetZaxis()->SetRange(cent1,cent2);
	//RppMass->GetZaxis()->SetRange(cent1,cent2);
	//MepMass->GetZaxis()->SetRange(cent1,cent2);
	//MeeMass->GetZaxis()->SetRange(cent1,cent2);
	//MppMass->GetZaxis()->SetRange(cent1,cent2);

	// for 2D 
	TH2D *h1_mp=(TH2D *) RepMass->Clone("h1mp");      //x mass;  y  pt
	TH2D *h2_mp=(TH2D *) ReeMass->Clone("h2mp");
	TH2D *h3_mp=(TH2D *) RppMass->Clone("h3mp");
	TH2D *h4_mp=(TH2D *) MepMass->Clone("h4mp");
	TH2D *h5_mp=(TH2D *) MeeMass->Clone("h5mp");
	TH2D *h6_mp=(TH2D *) MppMass->Clone("h6mp");

	// for 1D
	int ptb1=RepMass->GetYaxis()->FindBin(pt1+1e-8);
	int ptb2=RepMass->GetYaxis()->FindBin(pt2-1e-8);
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
	//const int Mpots=40;
	//int Mbin[Mpots+1]={0,5,10,15,20,25,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240,260,280,300,330,360,390,420,460,500,540,580,620,660,700,740,780,820,860,1000,1200};
  const int Mpots = 43;
  int Mbin[Mpots+1]={0,2,4,6,8,10,12,14,16,18,20,22,25,35,40,62,
    80,102,126,150,156,157,158,160,178,193,200,202,203,204,207,243,
    295,370,440,520,570,
    600,610,615,626,645,695,
    800};
  double mass[Mpots];
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
		//	cout<<"mass :"<<massb[i]<<endl;
	}	

	for(int i=0;i<Ppots+1;i++){	
		pb[i]=binWP*Pbin[i];
		//	cout<<"pt: "<<pb[i]<<" :"<<i<<endl;
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

		for(int i=1;i<Mpots;i++){
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

	LM =mixBack->FindBin(lowM);
	HM =mixBack->FindBin(hiM);

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
		MixBack22D=(TH2D*) MassP4->Clone();
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
      for(int i = 1; i<Mpots+1;i++) {
        double bct = signal1D->GetBinContent(i);
        double ber = signal1D->GetBinError(i);
        if(mass[i-1]<=lowMls || mass[i-1]>3.) continue;
        bct -= BackGround->GetBinContent(i)*(fFit2->Eval(mass[i-1])-1);
        ber = sqrt(ber*ber+pow(BackGround->GetBinError(i)*(fFit2->Eval(mass[i-1])-1),2.0));
        signal1D->SetBinContent(i,bct);
        signal1D->SetBinError(i,ber);
        BackGround->SetBinContent(i,BackGround->GetBinContent(i)*fFit2->Eval(mass[i-1]));
        BackGround->SetBinError(i,BackGround->GetBinError(i)*fFit2->Eval(mass[i-1]));
      } 
		}

		if(mAcc==2){
			signal2D->Add(MassP1,+1.);
			signal2D->Add(MixBack22D,-1.);
			int ptbL=signal2D->GetYaxis().FindBin(pt1*1.00001);
			int ptbH=signal2D->GetYaxis().FindBin(pt2*1.00001)-1;
			signal1D=(TH1D*) signal2D->ProjectionX("signal1D_x",ptbL,ptbH);	
			BackGround=(TH1D*) MixBack22D->ProjectionX("BackGround_x",ptbL,ptbH);	
      for(int i = 1; i<Mpots+1;i++) {
        double bct = signal1D->GetBinContent(i);
        double ber = signal1D->GetBinError(i);
        if(mass[i-1]<=lowMls || mass[i-1]>3.) continue;
        bct -= BackGround->GetBinContent(i)*(fFit2->Eval(mass[i-1])-1);
        ber = sqrt(ber*ber+pow(BackGround->GetBinError(i)*(fFit2->Eval(mass[i-1])-1),2.0));
        signal1D->SetBinContent(i,bct);
        signal1D->SetBinError(i,ber);
        BackGround->SetBinContent(i,BackGround->GetBinContent(i)*fFit2->Eval(mass[i-1]));
        BackGround->SetBinError(i,BackGround->GetBinError(i)*fFit2->Eval(mass[i-1]));
      } 
		}
	}

  for(int i = 1; i<Mpots+1;i++) {
    
    h11->SetBinError(i,h11->GetBinError(i)/nevent/mass_Width[i-1]/2.0);
    h11->SetBinContent(i,h11->GetBinContent(i)/nevent/mass_Width[i-1]/2.0);

    mixBack->SetBinError(i,mixBack->GetBinError(i)/nevent/mass_Width[i-1]/2.0);
    mixBack->SetBinContent(i,mixBack->GetBinContent(i)/nevent/mass_Width[i-1]/2.0);
    mLikeSign->SetBinError(i,mLikeSign->GetBinError(i)/nevent/mass_Width[i-1]/2.0);
    mLikeSign->SetBinContent(i,mLikeSign->GetBinContent(i)/nevent/mass_Width[i-1]/2.0);
    BackGround->SetBinError(i,BackGround->GetBinError(i)/nevent/mass_Width[i-1]/2.0);
    BackGround->SetBinContent(i,BackGround->GetBinContent(i)/nevent/mass_Width[i-1]/2.0);
    signal1D->SetBinError(i,signal1D->GetBinError(i)/nevent/mass_Width[i-1]/2.0);
    signal1D->SetBinContent(i,signal1D->GetBinContent(i)/nevent/mass_Width[i-1]/2.0);
  }




	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//   Plots
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas *fig0 = new TCanvas("fig0","",800,600);
  fig0->SetFrameLineWidth(3);
  fig0->cd();
  fig0->SetLogy();
  TH1F *frame0 = new TH1F("frame0","",1,0,4);
  frame0->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
  frame0->GetYaxis()->SetTitle("dN/dM_{ee} (GeV/c^{2})^{-1}");
  frame0->SetMaximum(1);
  frame0->SetMinimum(2e-9);
  frame0->Draw();
  h11->SetLineWidth(2);
  h11->Draw("HISTsame");
  mLikeSign->SetLineWidth(2);
  mLikeSign->SetLineColor(6);
  mLikeSign->SetLineStyle(4);
  mLikeSign->Draw("HISTsame");
  mixBack->SetLineWidth(2);
  mixBack->SetLineColor(8);
  mixBack->SetLineStyle(7);
  mixBack->Draw("HISTsame");
  BackGround->SetLineWidth(2);
  BackGround->SetLineColor(4);
  BackGround->SetLineStyle(9);
  BackGround->Draw("HISTsame");
  signal1D->SetMarkerStyle(20);
  signal1D->SetMarkerSize(1.2);
  signal1D->SetMarkerColor(2);
  signal1D->SetLineColor(2);
  signal1D->Draw("epsame");

  TLegend *t0 = new TLegend(0.5,0.45,0.8,0.8);
  t0->SetFillColor(10);
  t0->SetBorderSize(0);
  t0->SetTextSize(0.045);
  t0->SetHeader("Run12 p+p #sqrt{s}=200GeV");
  t0->AddEntry(h11,"foreground","l");
  t0->AddEntry(mLikeSign,"likesign","l");
  t0->AddEntry(mixBack,"mixevent","l");
  t0->AddEntry(BackGround,"likeS + Mix(likeS residual)","l");
  t0->AddEntry(signal1D,"Raw signal","epl");
  t0->Draw("same");

  fig0->SaveAs("RawSignal.gif");
  fig0->SaveAs("RawSignal.pdf");

  char name[50];
  sprintf(name,"RawYield_%iD.root",mAcc);
  TFile *f1 = new TFile(name,"recreate");
  f1->cd();
  signal1D->SetName("Signal");
  mLikeSign->SetName("LikeSign");
  mixBack->SetName("mixBack");
  BackGround->SetName("mixLSRD");

  signal1D->Write();
  mLikeSign->Write();
  mixBack->Write();
  BackGround->Write();
	return;
}







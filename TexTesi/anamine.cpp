#include <TGeoManager.h>
#include <TString.h>
#include <TGeoNode.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TApplication.h> 
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TDirectoryFile.h>

#include "/mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu/include/EDepSim/TG4Event.h"
#include "/mnt/c/Linux/Dune/edep-sim/edep-gcc-7-x86_64-linux-gnu/include/EDepSim/TG4HitSegment.h"

#include "/mnt/c/Linux/Dune/kloe-simu/include/struct.h"
#include "/mnt/c/Linux/Dune/kloe-simu/include/utils.h"

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

TString path = "/mnt/c/Linux/Dune/kloe-simu/files/ana/";
TFile fout(path + "result.root","RECREATE");

double kloe_center[] = {0.0000000 * 10, -238.47300 * 10, 2391.0000 * 10};
double emcalo_barrel_mod_halfw =  29.285 * 10;
double emcalo_barrel_mod_halfl = 215. * 10;
double emcalo_barrel_mod_h = 11.5 * 2 * 10;
double emcalo_int_radius = 2000;
double emcalo_ext_radius = TMath::Sqrt((emcalo_int_radius+emcalo_barrel_mod_h)*(emcalo_int_radius+emcalo_barrel_mod_h)+emcalo_barrel_mod_halfw*emcalo_barrel_mod_halfw);
//double dr = sqrt( yv_true*yv_true + zv_true*zv_true ) - emcalo_int_radius;
double threshold = 35.;

TCut fiducial = TString::Format("maxde_frontoutlayer < %f && abs(xv_reco) < 1500 ",threshold).Data();
TCut fiducial_truth = "dr < 176 && abs(xv_true) < 1500";
TCut CC = "isCC == 1";
TCut muonOK = "isPmuOK == 1 && !(pxmu_reco == 0 && pymu_reco == 0 && pzmu_reco == 0)";
TCut QE = "NHitLayer[0] == 1 ||  NHitLayer[0] == 2";
TCut QualityCut = "chi2_cr < 1E5";
TCut lowE = "Enu_true/1E3<20";

TCanvas c1("c1","c1",1000,1000,1000,1000);

void EvalChi2(TH1D& h1, TH1D& h2, double& chi2, int& ndof)
{
  double n1, n2;
  double vchi2 = 0.;
  
  chi2 = 0.;
  
  ndof =  h1.GetNbinsX();

  for(int i = 0; i < ndof; i++)
  {
    n1 = h1.GetBinContent(i+1);
    n2 = h2.GetBinContent(i+1);
    
    if(n1 != 0 || n2 != 0)
      vchi2 = (n1 - n2) * (n1 - n2) / (n1 + n2);
    else
      vchi2 = 0.;
    
    chi2 += vchi2;
  }
}

void FillHisto(TChain& t1, TChain& t2, TH1D& h1, TH1D& h2, int nEv)
{
  TString name1 = TString::Format("h1_n%d",nEv);
  TString name2 = TString::Format("h2_n%d",nEv); 

  h1.SetName(name1.Data());
  h2.SetName(name2.Data());

  int long n1 = t1.Draw(TString::Format("pmu_reco/1E3>>%s",name1.Data()).Data(),fiducial + CC + muonOK + QualityCut+lowE,"E0", nEv);
  int long n2 = t2.Draw(TString::Format("pmu_reco/1E3>>%s",name2.Data()).Data(),fiducial + CC + muonOK + QualityCut,"E0", nEv);
  
  std::cout << "nEv1: " << nEv << " -> " << n1 << " " << std::endl;
  std::cout << "nEv2: " << nEv << " -> " << n2 << " " << std::endl;
}

void processChain(TChain& t)
{
  
  std::cout << fiducial.GetName() << " " << fiducial.GetTitle() << std::endl;

  /////////////////////////////////////////////// all yx vertex position
  
  c1.DrawFrame(kloe_center[0]-3000,kloe_center[1]-3000,kloe_center[0]+3000,kloe_center[1]+3000);
  t.Draw("yv_true:xv_true","","same",10000,0);
  TBox b(kloe_center[0]-emcalo_barrel_mod_halfl,kloe_center[1]-emcalo_ext_radius,kloe_center[0]+emcalo_barrel_mod_halfl,kloe_center[1]+emcalo_ext_radius);
  b.SetFillStyle(0);
  b.Draw();
   
  c1.SaveAs(path + "result.pdf(","pdf"); 
  c1.SaveAs(path+"xy_true.png");
 
  

  ////////////////////////////////////////////////// all yz vertex position
  
 
  c1.DrawFrame(kloe_center[2]-3000,kloe_center[1]-3000,kloe_center[2]+3000,kloe_center[1]+3000);
  t.Draw("yv_true:zv_true","","same",10000,0);
  TEllipse elint(kloe_center[2],kloe_center[1],emcalo_int_radius,emcalo_int_radius);
  TEllipse elext(kloe_center[2],kloe_center[1],emcalo_ext_radius,emcalo_ext_radius);
  elint.SetFillStyle(0);
  elext.SetFillStyle(0);
  elint.Draw();
  elext.Draw();
  

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "yz_true.png");
  
  
  ////////////////////////////////////////////////////////////// int and ext yx vertex position
  
  c1.DrawFrame(kloe_center[0]/1000-3,kloe_center[1]/1000-3,kloe_center[0]/1000+3,kloe_center[1]/1000+3);
  TH1F *htemp = (TH1F*)gPad->GetPrimitive("hframe");
  htemp->GetYaxis()->SetTitle("y[m]");
  htemp->GetXaxis()->SetTitle("x[m]");

  t.SetMarkerColor(kBlue);
  t.Draw("yv_true/1000:xv_true/1000",!fiducial,"same",10000,0);
  
  t.SetMarkerColor(kBlack);
  t.Draw("yv_true/1000:xv_true/1000",fiducial,"same",10000,0);
  TBox b2(kloe_center[0]*0.001-emcalo_barrel_mod_halfl*0.001,kloe_center[1]*0.001-emcalo_ext_radius*0.001,kloe_center[0]*0.001+emcalo_barrel_mod_halfl*0.001,kloe_center[1]*0.001+emcalo_ext_radius*0.001);
  b2.SetFillStyle(0);
  b2.Draw();
  
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "xy_fiducial.png");
  

  ///////////////////////////////////////////////////////////// int and ext yz vertex position
  
  c1.DrawFrame(kloe_center[2]/1000-3,kloe_center[1]/1000-3,kloe_center[2]/1000+3,kloe_center[1]/1000+3);
  TH1F *htemp2 = (TH1F*)gPad->GetPrimitive("hframe");
  htemp2->GetYaxis()->SetTitle("z[m]");
  htemp2->GetXaxis()->SetTitle("y[m]");
  t.SetMarkerColor(kBlue);
  t.Draw("yv_true/1000:zv_true/1000",!fiducial,"same",10000,0);
  t.SetMarkerColor(kBlack);
  t.Draw("yv_true/1000:zv_true/1000",fiducial,"same",10000,0);
  TEllipse elint2(kloe_center[2]*0.001,kloe_center[1]*0.001,emcalo_int_radius*0.001,emcalo_int_radius*0.001);
  TEllipse elext2(kloe_center[2]*0.001,kloe_center[1]*0.001,emcalo_ext_radius*0.001,emcalo_ext_radius*0.001);
  elint2.SetFillStyle(0);
  elext2.SetFillStyle(0);
  elint2.Draw();
  elext2.Draw();

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "yz_fiducial.png");
  
 
  /////////////////////////////////////////////////////////// neutrino energy before and after cut
  
  c1.SetRightMargin(0.13);
  c1.SetLeftMargin(0.13);
  
  TH1D h_Enu("h_Enu","",80,0,40);
  long int nCC = t.Draw("Enu_true/1E3>>h_Enu",CC,"");
  TH1D h_Enu_int("h_Enu_int","",80,0,40);
  long int nFiducial = t.Draw("Enu_true/1E3>>h_Enu_int",fiducial + CC,"");
  
  gStyle->SetOptStat(0);
  h_Enu.SetTitle("");//Neutrino energy distribution 
  h_Enu.GetXaxis()->SetTitle("E_{#nu}[GeV]");
  h_Enu.GetYaxis()->SetTitle("events");
  h_Enu.SetFillColor(kGray);
  h_Enu.SetLineColor(kBlack);
  h_Enu.SetFillStyle(3001);
  h_Enu_int.SetFillColor(kBlue-10);
  h_Enu_int.SetLineColor(kBlue);
  h_Enu_int.SetFillStyle(3001);
  h_Enu.Draw();
  h_Enu_int.Draw("same");
  
  TLegend* legend1 = new TLegend(0.6,0.75,0.9,0.9);
  legend1->AddEntry("h_Enu", "All neutrinos", "f");
  legend1->AddEntry((TObject*)0, TString::Format("entries: %.0f",h_Enu.GetEntries()).Data(), "");  
  legend1->AddEntry("h_Enu_int", "Neutrinos surviving the cut", "f");
  legend1->AddEntry((TObject*)0, TString::Format("entries: %.0f",h_Enu_int.GetEntries()).Data(), ""); 
  legend1->SetFillColor(0);
  legend1->SetBorderSize(1);
  legend1->Draw("");
  
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "Enu_cut.png");

  ///////////////////////////////////////////////////////// fiducial volume selection efficiency
  
  TEfficiency eff_fid(h_Enu_int,h_Enu);
  eff_fid.SetTitle(";E [GeV];#varepsilon_{sel}(E)");
  eff_fid.Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "fiducial_eff.png");
  std::cout << "Fiducial Volume selection efficiency: " << double(nFiducial)/double(nCC) << std::endl;

  ///////////////////////////////////////////////////////// fraction of CC events
  
  long int nEntries = t.GetEntries();
  std::cout << "CC fraction: " << double(nCC)/double(nEntries) << "; " << nEntries <<"; " << nCC <<std::endl;

  ///////////////////////////////////////////////////////// pmutrue of all muons and correctly reconstructed muons
  
  TH1D h_pmu("h_pmu","",80,0,40); 
  TH1D h_pmu_ok("h_pmu_ok","",80,0,40);
  long int nintmu = t.Draw("pmu_true/1E3>>h_pmu",fiducial + CC,"");
  long int nintmu_ok = t.Draw("pmu_true/1E3>>h_pmu_ok",fiducial + CC + muonOK,"");
  h_pmu.SetTitle("");//Neutrino energy distribution 
  h_pmu.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  h_pmu.GetYaxis()->SetTitle("events");
  h_pmu.SetFillColor(kBlue-10);
  h_pmu.SetLineColor(kBlue);
  h_pmu.SetFillStyle(3001);
  h_pmu_ok.SetFillColor(kGreen-10);
  h_pmu_ok.SetLineColor(kGreen);
  h_pmu_ok.SetFillStyle(3001);
  h_pmu.Draw();
  h_pmu_ok.Draw("same");
  
  TLegend* legend2 = new TLegend(0.6,0.75,0.9,0.9);
  legend2->AddEntry("h_pmu", "All muons", "f");
  legend2->AddEntry((TObject*)0, TString::Format("entries: %.0f",h_pmu.GetEntries()).Data(), ""); 
  legend2->AddEntry("h_pmu_ok", "Correctly reconstructed muons", "f");
  legend2->AddEntry((TObject*)0, TString::Format("entries: %.0f",h_pmu_ok.GetEntries()).Data(), ""); 
  legend2->SetFillColor(0);
  legend2->SetBorderSize(1);
  legend2->Draw("");
  
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "pmutrue_reco.png");

  ///////////////////////////////////////////////////////// reconstruction efficiency
  
  TEfficiency eff_murec(h_pmu_ok,h_pmu);
  eff_murec.Draw("");
  eff_murec.SetTitle(";p [GeV/c];#varepsilon_{reco}(p)");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "eff_reco.png");
  std::cout << "Reconstucted muon efficiency: " << double(nintmu_ok)/double(nintmu) << " "<< nintmu_ok << " " << nintmu<<std::endl;


  ///////////////////////////////////////////////////////// pmutrue vs pmureco of correctly reconstructed muons after cut
  
  TH2D h_p("h_p","",200,0,40,200,0,40);
  t.Draw("pmu_reco/1E3:pmu_true/1E3>>h_p",fiducial + CC + muonOK,"");
  h_p.SetTitle("");
  h_p.GetXaxis()->SetTitle("p_{#mu}^{true}[GeV/c]");
  h_p.GetYaxis()->SetTitle("p_{#mu}^{reco}[GeV/c]");
  h_p.Draw("colz");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "pmu_recoVSpm_true.png");

  ///////////////////////////////////////////////////////// pmutrue + pmureco of correctly reconstructed muons after cut
  
  TH1D h_pmureco_cut("h_pmureco_cut","",80,0,40);
  TH1D h_pmutrue_cut("h_pmutrue_cut","",80,0,40);
  long int pr= t.Draw("pmu_reco/1E3>>h_pmureco_cut",fiducial + CC + muonOK,""); 
  long int pt= t.Draw("pmu_true/1E3>>h_pmutrue_cut",fiducial + CC + muonOK,"");
  h_pmureco_cut.SetTitle("");//Neutrino energy distribution 
  h_pmureco_cut.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  h_pmureco_cut.GetYaxis()->SetTitle("events");
  h_pmureco_cut.SetFillColor(kBlue-10);
  h_pmureco_cut.SetLineColor(kBlue);
  h_pmureco_cut.SetFillStyle(3001);
  h_pmutrue_cut.SetFillColor(kGreen-10);
  h_pmutrue_cut.SetLineColor(kGreen);
  h_pmutrue_cut.SetFillStyle(3001);

  h_pmureco_cut.Draw();
  h_pmutrue_cut.Draw("same");
  
  TLegend* legend3 = new TLegend(0.6,0.75,0.9,0.9);
  legend3->AddEntry("h_pmutrue_cut", "p_{#mu}^{true}", "f");
  legend3->AddEntry("h_pmureco_cut", "p_{#mu}^{reco}", "f");
  legend3->SetFillColor(0);
  legend3->SetBorderSize(1);
  legend3->Draw("");


  c1.SaveAs(path + "result.pdf","pdf"); 
  c1.SaveAs(path + "pmu_recoANDpm_true.png");  

}

void anamine()
{
  gROOT->SetBatch(1);
  TChain tNom("tout","tNom");
  TChain tH1Y("tout","tH1Y");
  tNom.Add(path + "nominal/*");
  tH1Y.Add(path + "shift/*");
  std::cout << "H1Y: " << tH1Y.GetEntries() << "; " <<std::endl;
  
  //////////////////////////////////////////////SetALIASES
  
  tNom.SetAlias("pmu_true","sqrt(pxmu_true*pxmu_true+pymu_true*pymu_true+pzmu_true*pzmu_true)");
  tH1Y.SetAlias("pmu_true","sqrt(pxmu_true*pxmu_true+pymu_true*pymu_true+pzmu_true*pzmu_true)");
  
  tNom.SetAlias("pmu_reco","sqrt(pxmu_reco*pxmu_reco+pymu_reco*pymu_reco+pzmu_reco*pzmu_reco)");
  tH1Y.SetAlias("pmu_reco","sqrt(pxmu_reco*pxmu_reco+pymu_reco*pymu_reco+pzmu_reco*pzmu_reco)");
  
  tNom.SetAlias("ptmu_true","sqrt(pymu_true*pymu_true+pzmu_true*pzmu_true)");
  tH1Y.SetAlias("ptmu_true","sqrt(pymu_true*pymu_true+pzmu_true*pzmu_true)");

  tNom.SetAlias("ptmu_reco","sqrt(pymu_reco*pymu_reco+pzmu_reco*pzmu_reco)");
  tH1Y.SetAlias("ptmu_reco","sqrt(pymu_reco*pymu_reco+pzmu_reco*pzmu_reco)");
  
  tNom.SetAlias("reldptmu","1-ptmu_true/ptmu_reco");
  tH1Y.SetAlias("reldptmu","1-ptmu_true/ptmu_reco");

  tNom.SetAlias("reldplmu","1-abs(pxmu_true)/abs(pxmu_reco)");
  tH1Y.SetAlias("reldplmu","1-abs(pxmu_true)/abs(pxmu_reco)");
  

  tNom.SetAlias("red_chi2_cr","2*chi2_cr/nHit");
  tH1Y.SetAlias("red_chi2_cr","2*chi2_cr/nHit");
  
  tNom.SetAlias("dr","sqrt( (yv_true+2384.73)*(yv_true+2384.73) + (zv_true-23910)*(zv_true-23910) ) - 2000");
  tH1Y.SetAlias("dr","sqrt( (yv_true+2384.73)*(yv_true+2384.73) + (zv_true-23910)*(zv_true-23910) ) - 2000");

  tNom.SetAlias("red_chi2_ln","2*chi2_ln/nHit");
  tH1Y.SetAlias("red_chi2_ln","2*chi2_ln/nHit");
  ////////////////////////////////////////////////////////////////////

  processChain(tNom);


  int long NlowEok = tNom.GetEntries(lowE);
  
  
  ///////////////////////////////////////////////////////////////////Enu NOM vs SHIFT

  TH1D hEnuNom("hEnuNom","",100,0,20);
  TH1D hEnuH1Y("hEnuH1Y","",100,0,20);
  int long n1 = tNom.Draw("Enu_true/1E3>>hEnuNom",lowE,"E0");
  int long n2 = tH1Y.Draw("Enu_true/1E3>>hEnuH1Y",lowE,"E0",NlowEok,0);

  if(n1 != n2)
  {
    std::cout << "error: samples are different" << n1 << " != " << n2 << std::endl;
    return;
  }

  hEnuNom.SetTitle("");//Neutrino energy distribution 
  hEnuNom.GetXaxis()->SetTitle("E_{#nu}[GeV]");
  hEnuNom.GetYaxis()->SetTitle("events");
  //hEnuNom.SetFillColor(kBlue-10);
  hEnuNom.SetLineColor(kBlue);
  hEnuNom.SetFillStyle(3001);
  //hEnuH1Y.SetFillColor(kRed-10);
  hEnuH1Y.SetLineColor(kRed);
  hEnuH1Y.SetFillStyle(3001);
  c1.cd();
  hEnuNom.Draw();
  hEnuH1Y.Draw("same");
  
  TLegend* legend4 = new TLegend(0.6,0.75,0.9,0.9);
  legend4->AddEntry("hEnuNom", "E_{#nu} from nominal flux", "f");
  legend4->AddEntry("hEnuH1Y", "E_{#nu} from shifted flux", "f");
  legend4->SetFillColor(0);
  legend4->SetBorderSize(1);
  legend4->Draw("");

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "EnuNOMvsSHIFT.png");


  /////////////////////////////////////////////////////////////////MOMENTUM HISTOGRAMS FROM TWO FLUXES
  
  TH1D hNom("hNom","",40,0,16);
  TH1D hH1Y("hH1Y","",40,0,16);
  
  TH1D hNomHalf1("hNomHalf1","",40,0,16);
  TH1D hNomHalf2("hNomHalf2","",40,0,16);

  TH1D hNom_true("hNom_true","",40,0,16);
  TH1D hH1Y_true("hH1Y_true","",40,0,16);
  
  TH1D hNom_true_all("hNom_true_all","",40,0,16);
  TH1D hH1Y_true_all("hH1Y_true_all","",40,0,16);

  TH1D hNom_true_fid("hNom_true_fid","",40,0,16);
  TH1D hH1Y_true_fid("hH1Y_true_fid","",40,0,16);
  
  TH1D hNom_true_fid_QE("hNom_true_fid_QE","",40,0,16);
  TH1D hH1Y_true_fid_QE("hH1Y_true_fid_QE","",40,0,16);
  
  TH1D hNom_true_fid_Q("hNom_true_fid_Q","",40,0,16);
  TH1D hH1Y_true_fid_Q("hH1Y_true_fid_Q","",40,0,16);
  
  TH1D hNom_reco_fid_Q("hNom_reco_fid_Q","",40,0,16);
  TH1D hH1Y_reco_fid_Q("hH1Y_reco_fid_Q","",40,0,16);
  
  ////////////Control half statistic histograms from same flux (preco)
  
  tNom.Draw("pmu_reco/1E3>>hNomHalf1",fiducial + CC + muonOK + lowE,"E0",tNom.GetEntries()*0.5,0);
  tNom.Draw("pmu_reco/1E3>>hNomHalf2",fiducial + CC + muonOK + lowE,"E0",tNom.GetEntries()*0.5,tNom.GetEntries()*0.5);
  
  hNomHalf1.GetXaxis()->SetTitle("p_{#mu}^{reco}[GeV/c]");
  hNomHalf1.GetYaxis()->SetTitle("events");
  //hNomHalf1.SetFillColor(kBlue-10);
  hNomHalf1.SetLineColor(kBlue);
  hNomHalf1.SetFillStyle(3001);
  //hNomHalf2.SetFillColor(kGreen-10);
  hNomHalf2.SetLineColor(kGreen);
  hNomHalf2.SetFillStyle(3001);

  c1.SetLogy(true);
  hNomHalf1.Draw();
  hNomHalf2.Draw("same");
  TLegend* legend = new TLegend(0.6,0.75,0.9,0.9);
  legend->AddEntry("hNomHalf1", TString::Format("entries: %.0f",hNomHalf1.GetEntries()).Data(), "f");
  legend->AddEntry("hNomHalf2", TString::Format("entries: %.0f",hNomHalf2.GetEntries()).Data(), "f");
  legend->SetFillColor(0);
  legend->SetBorderSize(1);
  legend->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "precoHalf.png");
  c1.SetLogy(false);

  ////////////Two fluxes only correctly reconstructed muons (preco)
  
  tNom.Draw("pmu_reco/1E3>>hNom",fiducial + CC + muonOK + lowE,"E0");
  tH1Y.Draw("pmu_reco/1E3>>hH1Y",fiducial + CC + muonOK + lowE,"E0",NlowEok,0);

  hNom.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom.GetYaxis()->SetTitle("events");
  //hNom.SetFillColor(kBlue-10);
  hNom.SetLineColor(kBlue);
  //hNom.SetFillStyle(3001);
  //hH1Y.SetFillColor(kRed-10);
  hH1Y.SetLineColor(kRed);
  //hH1Y.SetFillStyle(3001);

  c1.SetLogy(true);
  hNom.Draw();
  hH1Y.Draw("same");
  TLegend* legend4a = new TLegend(0.6,0.75,0.9,0.9);
  legend4a->AddEntry("hNom", "p_{#mu}^{reco} nominal", "f");
  legend4a->AddEntry((TObject*)0, TString::Format("entries: %.0f",hNom.GetEntries()).Data(), ""); 
  legend4a->AddEntry("hH1Y", "p_{#mu}^{reco} shifted", "f");
  legend4a->AddEntry((TObject*)0, TString::Format("entries: %.0f",hH1Y.GetEntries()).Data(), "");
  legend4a->SetFillColor(0);
  legend4a->SetBorderSize(1);
  legend4a->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "precoFluxesMuonOK.png");
  c1.SetLogy(false);


  ////////////Two fluxes only correctly reconstructed muons (ptrue)
  
  tNom.Draw("pmu_true/1E3>>hNom_true",fiducial + CC + muonOK + lowE,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true",fiducial + CC + muonOK + lowE,"E0",NlowEok,0);
  
  hNom_true.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom_true.GetYaxis()->SetTitle("events");
  //hNom_true.SetFillColor(kBlue-10);
  hNom_true.SetLineColor(kBlue);
  hNom_true.SetFillStyle(3001);
  //hH1Y_true.SetFillColor(kRed-10);
  hH1Y_true.SetLineColor(kRed);
  hH1Y_true.SetFillStyle(3001);


  c1.SetLogy(true);
  hNom_true.Draw();
  hH1Y_true.Draw("same");
  TLegend* legend4b = new TLegend(0.6,0.75,0.9,0.9);
  legend4b->AddEntry("hNom_true", "p_{#mu}^{true} nominal", "f");
  legend4b->AddEntry((TObject*)0, TString::Format("entries: %.0f",hNom_true.GetEntries()).Data(), ""); 
  legend4b->AddEntry("hH1Y_true", "p_{#mu}^{true} shifted", "f");
  legend4b->AddEntry((TObject*)0, TString::Format("entries: %.0f",hH1Y_true.GetEntries()).Data(), ""); 
  legend4b->SetFillColor(0);
  legend4b->SetBorderSize(1);
  legend4b->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "ptrueFluxesMuonOK.png");
  c1.SetLogy(false);
  
  ////////////Two fluxes all muons (ptrue)
  
  tNom.Draw("pmu_true/1E3>>hNom_true_all", fiducial_truth + CC + lowE,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_all", fiducial_truth + CC + lowE,"E0",NlowEok,0);
  
  hNom_true_all.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom_true_all.GetYaxis()->SetTitle("events");
  //hNom_true_fid.SetFillColor(kBlue-10);
  hNom_true_all.SetLineColor(kBlue);
  hNom_true_all.SetFillStyle(3001);
  //hH1Y_true_fid.SetFillColor(kRed-10);
  hH1Y_true_all.SetLineColor(kRed);
  hH1Y_true_all.SetFillStyle(3001);
  
  c1.SetLogy(true);
  hNom_true_all.Draw();
  hH1Y_true_all.Draw("same");
  TLegend* legend5all = new TLegend(0.6,0.75,0.9,0.9);
  legend5all->AddEntry("hNom_true_all", "p_{#mu}^{true} nominal", "f");
  legend5all->AddEntry((TObject*)0, TString::Format("entries: %.0f",hNom_true_all.GetEntries()).Data(), ""); 
  legend5all->AddEntry("hH1Y_true_all", "p_{#mu}^{true} shifted", "f");
  legend5all->AddEntry((TObject*)0, TString::Format("entries: %.0f",hH1Y_true_all.GetEntries()).Data(), ""); 
  legend5all->SetFillColor(0);
  legend5all->SetBorderSize(1);
  legend5all->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "ptrueFluxeseverything.png");
  c1.SetLogy(false);

  
  ////////////Two fluxes all muons in the fiducial (ptrue)

  tNom.Draw("pmu_true/1E3>>hNom_true_fid",fiducial + CC + lowE,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_fid",fiducial + CC + lowE,"E0",NlowEok,0);
  
  hNom_true_fid.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom_true_fid.GetYaxis()->SetTitle("events");
  //hNom_true_fid.SetFillColor(kBlue-10);
  hNom_true_fid.SetLineColor(kBlue);
  hNom_true_fid.SetFillStyle(3001);
  //hH1Y_true_fid.SetFillColor(kRed-10);
  hH1Y_true_fid.SetLineColor(kRed);
  hH1Y_true_fid.SetFillStyle(3001);


  c1.SetLogy(true);
  hNom_true_fid.Draw();
  hH1Y_true_fid.Draw("same");
  TLegend* legend5 = new TLegend(0.6,0.75,0.9,0.9);
  legend5->AddEntry("hNom_true_fid", "p_{#mu}^{true} nominal", "f");
  legend5->AddEntry((TObject*)0, TString::Format("entries: %.0f",hNom_true_fid.GetEntries()).Data(), ""); 
  legend5->AddEntry("hH1Y_true_fid", "p_{#mu}^{true} shifted", "f");
  legend5->AddEntry((TObject*)0, TString::Format("entries: %.0f",hH1Y_true_fid.GetEntries()).Data(), ""); 
  legend5->SetFillColor(0);
  legend5->SetBorderSize(1);
  legend5->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "ptrueFluxesAll.png");
  c1.SetLogy(false);
 
  
  ///////////////////////////////////QEcut (p_reco) 
  
  tNom.Draw("pmu_reco/1E3>>hNom_true_fid_QE",fiducial + CC + muonOK + lowE + QE,"E0");
  tH1Y.Draw("pmu_reco/1E3>>hH1Y_true_fid_QE",fiducial + CC + muonOK + lowE + QE,"E0",NlowEok,0);
  
  hNom_true_fid_QE.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom_true_fid_QE.GetYaxis()->SetTitle("events");
  hNom_true_fid_QE.SetFillColor(kBlue-10);
  hNom_true_fid_QE.SetLineColor(kBlue);
  hNom_true_fid_QE.SetFillStyle(3001);
  hH1Y_true_fid_QE.SetFillColor(kRed-10);
  hH1Y_true_fid_QE.SetLineColor(kRed);
  hH1Y_true_fid_QE.SetFillStyle(3001);
  
  c1.SetLogy(true);
  hNom_true_fid_QE.Draw();
  hH1Y_true_fid_QE.Draw("same");
  TLegend* legend6 = new TLegend(0.6,0.75,0.9,0.9);
  legend6->AddEntry("hNom_true_fid_QE", "p_{#mu}^{reco} nominal", "f");
  legend6->AddEntry("hH1Y_true_fid_QE", "p_{#mu}^{reco} shifted", "f");
  legend6->SetFillColor(0);
  legend6->SetBorderSize(1);
  legend6->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "precoFluxesQE.png");
  c1.SetLogy(false);
  
  /////////////////////////////////ptrueVSpreco QE
  
  TH2D h_pQE("h_pQE","",200,0,40,200,0,20);
  h_pQE.GetXaxis()->SetTitle("p_{#mu}^{true} [GeV/c]");
  h_pQE.GetYaxis()->SetTitle("p_{#mu}^{reco} [GeV/c]");
  tNom.Draw("pmu_reco/1E3:pmu_true/1E3>>h_pQE",fiducial + CC + muonOK + lowE + QE,"");
  h_pQE.Draw("colz");

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "ptrueVSprecoQE.png");
  
  
  ///////////////////////////////////Quality cut (p_true) 
  
  tNom.Draw("pmu_true/1E3>>hNom_true_fid_Q",fiducial + CC + muonOK + lowE + QualityCut,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_fid_Q",fiducial + CC + muonOK + lowE + QualityCut,"E0",NlowEok,0);
  
  hNom_true_fid_Q.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom_true_fid_Q.GetYaxis()->SetTitle("events");
  //hNom_true_fid_Q.SetFillColor(kBlue-10);
  hNom_true_fid_Q.SetLineColor(kBlue);
  //hNom_true_fid_Q.SetFillStyle(3001);
  //hH1Y_true_fid_Q.SetFillColor(kRed-10);
  hH1Y_true_fid_Q.SetLineColor(kRed);
  //hH1Y_true_fid_Q.SetFillStyle(3001);

  c1.SetLogy(true);
  hNom_true_fid_Q.Draw();
  hH1Y_true_fid_Q.Draw("same");
  TLegend* legendQ = new TLegend(0.6,0.75,0.9,0.9);
  legendQ->AddEntry("hNom_true_fid_Q", "p_{#mu}^{true} nominal", "f");
  legendQ->AddEntry((TObject*)0, TString::Format("entries: %.0f",hNom_true_fid_Q.GetEntries()).Data(), "");  
  legendQ->AddEntry("hH1Y_true_fid_Q", "p_{#mu}^{true} shifted", "f");
  legendQ->AddEntry((TObject*)0, TString::Format("entries: %.0f",hH1Y_true_fid_Q.GetEntries()).Data(), "");
  
  legendQ->SetFillColor(0);
  legendQ->SetBorderSize(1);
  legendQ->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "ptrueFluxesQ.png");
  c1.SetLogy(false);
  
  ///////////////////////////////////Quality cut (p_reco) 
  
  tNom.Draw("pmu_reco/1E3>>hNom_reco_fid_Q",fiducial + CC + muonOK + lowE + QualityCut,"E0");
  tH1Y.Draw("pmu_reco/1E3>>hH1Y_reco_fid_Q",fiducial + CC + muonOK + lowE + QualityCut,"E0",NlowEok,0);
  
  hNom_reco_fid_Q.GetXaxis()->SetTitle("p_{#mu}[GeV/c]");
  hNom_reco_fid_Q.GetYaxis()->SetTitle("events");
  //hNom_true_fid_Q.SetFillColor(kBlue-10);
  hNom_reco_fid_Q.SetLineColor(kBlue);
  //hNom_true_fid_Q.SetFillStyle(3001);
  //hH1Y_true_fid_Q.SetFillColor(kRed-10);
  hH1Y_reco_fid_Q.SetLineColor(kRed);
  //hH1Y_true_fid_Q.SetFillStyle(3001);

  c1.SetLogy(true);
  hNom_reco_fid_Q.Draw();
  hH1Y_reco_fid_Q.Draw("same");
  TLegend* legendQr = new TLegend(0.6,0.75,0.9,0.9);
  legendQr->AddEntry("hNom_reco_fid_Q", "p_{#mu}^{reco} nominal", "f");
  legendQr->AddEntry((TObject*)0, TString::Format("entries: %.0f",hNom_reco_fid_Q.GetEntries()).Data(), "");  
  legendQr->AddEntry("hH1Y_reco_fid_Q", "p_{#mu}^{reco} shifted", "f");
  legendQr->AddEntry((TObject*)0, TString::Format("entries: %.0f",hH1Y_reco_fid_Q.GetEntries()).Data(), "");
  
  legendQr->SetFillColor(0);
  legendQr->SetBorderSize(1);
  legendQr->Draw("");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "precoFluxesQ.png");
  c1.SetLogy(false);
  
  /////////////////////////////////ptrueVSpreco Quality Cut
  
  TH2D h_pQ("h_pQ","",200,0,20,200,0,20);
  h_pQ.GetXaxis()->SetTitle("p_{#mu}^{true} [GeV/c]");
  h_pQ.GetYaxis()->SetTitle("p_{#mu}^{reco} [GeV/c]");
  tNom.Draw("pmu_reco/1E3:pmu_true/1E3>>h_pQ",fiducial + CC + muonOK + lowE + QualityCut,"");
  h_pQ.Draw("colz");

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "ptrueVSprecoQ.png");
  
  ////////////////////////////////Resolutions
  
  TH2D h_reldpt_chi2("h_reldpt_chi2","",100,0,1E6,100,-1,1);
  tNom.Draw("reldptmu:chi2_cr>>h_reldpt_chi2",fiducial + CC + muonOK,"");
  h_reldpt_chi2.SetTitle("");
  h_reldpt_chi2.GetXaxis()->SetTitle("#chi^{2}_{cr}");
  h_reldpt_chi2.GetYaxis()->SetTitle("res");
  h_reldpt_chi2.Draw("colz");

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "pt_chi2_cr.png");
  
  TH2D h_reldpt_redchi2("h_reldpt_redchi2","",100,0,10000,100,-1,1);
  tNom.Draw("reldptmu:red_chi2_cr>>h_reldpt_redchi2",fiducial + CC + muonOK,"");
  h_reldpt_redchi2.Draw("colz");

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "pt_redchi2_cr.png");
  
  TH2D h_reldpl_chi2("h_reldpl_chi2","",100,0,100,100,-1,1);
  tNom.Draw("reldplmu:chi2_ln>>h_reldpl_chi2",fiducial + CC + muonOK,"");
  h_reldpl_chi2.SetTitle("");
  h_reldpl_chi2.GetXaxis()->SetTitle("#chi^{2}_{ln}");
  h_reldpl_chi2.GetYaxis()->SetTitle("res");
  h_reldpl_chi2.Draw("colz");
  

  c1.SaveAs(path + "result.pdf","pdf");
  c1.SaveAs(path + "pl_chi2_ln.png");
  
  TH2D h_reldpl_redchi2("h_reldpl_redchi2","",100,0,10,100,-1,1);
  tNom.Draw("reldplmu:red_chi2_ln>>h_reldpl_redchi2",fiducial + CC + muonOK,"");
  h_reldpl_redchi2.Draw("colz");

  c1.SaveAs(path + "result.pdf","pdf");  
  c1.SaveAs(path + "pl_redchi2_ln.png");
  
  /////////////////////////////////////////////////////CHI2TEST//////////////////////////////////////////////////////////////////
  /*
  TH1D hDiff("hDiff","",40,0,16);
  TH1D hChi2("hChi2","",40,0,16);
  TH1D hChi2_half("hChi2_half","",40,0,16);
  TH1D hChi2_true("hChi2_true","",40,0,16);
  TH1D hChi2_true_fid("hChi2_true_fid","",40,0,16);
  TH1D hChi2_true_fid_QE("hChi2_true_fid_QE","",40,0,16);
  TH1D hChi2_true_fid_Q("hChi2_true_fid_Q","",40,0,16);
  TH1D hChi2_reco_fid_Q("hChi2_reco_fid_Q","",40,0,16);
  */
  int ndof;
  double chi2tot = 0.;
  double chi2 = 0.;
  
  double chi2tot_half = 0.;
  double chi2_half = 0.;

  double chi2tot_true = 0.;
  double chi2_true = 0.;
  
  double chi2tot_true_all = 0.;
  double chi2_true_all = 0.;

  double chi2tot_true_fid = 0.;
  double chi2_true_fid = 0.;
  
  double chi2tot_true_fid_QE = 0.;
  double chi2_true_fid_QE = 0.;
  
  double chi2tot_true_fid_Q = 0.;
  double chi2_true_fid_Q = 0.;
  
  double chi2tot_reco_fid_Q = 0.;
  double chi2_reco_fid_Q = 0.;
  
  EvalChi2(hNom, hH1Y, chi2, ndof);
  EvalChi2(hNomHalf1, hNomHalf2, chi2_half, ndof);
  EvalChi2(hNom_true, hH1Y_true, chi2_true, ndof);
  EvalChi2(hNom_true_all, hH1Y_true_all, chi2_true_all, ndof);
  EvalChi2(hNom_true_fid, hH1Y_true_fid, chi2_true_fid, ndof);
  EvalChi2(hNom_true_fid_Q, hH1Y_true_fid_Q, chi2_true_fid_Q, ndof);
  EvalChi2(hNom_reco_fid_Q, hH1Y_reco_fid_Q, chi2_reco_fid_Q, ndof);
	
  /*
  for(int i = 0; i < hDiff.GetNbinsX(); i++)
  {
    hDiff.SetBinContent(i+1,abs(hNom.GetBinContent(i+1)-hH1Y.GetBinContent(i+1))/hNom.GetBinContent(i+1));

    chi2 = (hNom.GetBinContent(i+1)-hH1Y.GetBinContent(i+1))*(hNom.GetBinContent(i+1)-hH1Y.GetBinContent(i+1))/(hNom.GetBinContent(i+1)+hH1Y.GetBinContent(i+1));
	chi2_half = (hNomHalf1.GetBinContent(i+1)-hNomHalf2.GetBinContent(i+1))*(hNomHalf1.GetBinContent(i+1)-hNomHalf2.GetBinContent(i+1))/(hNomHalf1.GetBinContent(i+1)+hNomHalf2.GetBinContent(i+1));
    chi2_true = (hNom_true.GetBinContent(i+1)-hH1Y_true.GetBinContent(i+1))*(hNom_true.GetBinContent(i+1)-hH1Y_true.GetBinContent(i+1))/(hNom_true.GetBinContent(i+1)+hH1Y_true.GetBinContent(i+1));
    chi2_true_fid = (hNom_true_fid.GetBinContent(i+1)-hH1Y_true_fid.GetBinContent(i+1))*(hNom_true_fid.GetBinContent(i+1)-hH1Y_true_fid.GetBinContent(i+1))/(hNom_true_fid.GetBinContent(i+1)+hH1Y_true_fid.GetBinContent(i+1));
    chi2_true_fid_QE = (hNom_true_fid_QE.GetBinContent(i+1)-hH1Y_true_fid_QE.GetBinContent(i+1))*(hNom_true_fid_QE.GetBinContent(i+1)-hH1Y_true_fid_QE.GetBinContent(i+1))/(hNom_true_fid_QE.GetBinContent(i+1)+hH1Y_true_fid_QE.GetBinContent(i+1));
	chi2_true_fid_Q = (hNom_true_fid_Q.GetBinContent(i+1)-hH1Y_true_fid_Q.GetBinContent(i+1))*(hNom_true_fid_Q.GetBinContent(i+1)-hH1Y_true_fid_Q.GetBinContent(i+1))/(hNom_true_fid_Q.GetBinContent(i+1)+hH1Y_true_fid_Q.GetBinContent(i+1));
	chi2_reco_fid_Q = (hNom_reco_fid_Q.GetBinContent(i+1)-hH1Y_reco_fid_Q.GetBinContent(i+1))*(hNom_reco_fid_Q.GetBinContent(i+1)-hH1Y_reco_fid_Q.GetBinContent(i+1))/(hNom_reco_fid_Q.GetBinContent(i+1)+hH1Y_reco_fid_Q.GetBinContent(i+1));
	
    chi2tot += chi2;
	chi2tot_half += chi2_half;
    chi2tot_true += chi2_true;
    chi2tot_true_fid += chi2_true_fid;
	chi2tot_true_fid_QE += chi2_true_fid_QE;
	chi2tot_true_fid_Q += chi2_true_fid_Q;
	chi2tot_reco_fid_Q += chi2_reco_fid_Q;

    hChi2.SetBinContent(i+1,chi2);
	hChi2_half.SetBinContent(i+1,chi2_half);
    hChi2_true.SetBinContent(i+1,chi2_true);
    hChi2_true_fid.SetBinContent(i+1,chi2_true_fid);
	hChi2_true_fid_QE.SetBinContent(i+1,chi2_true_fid_QE);
	hChi2_true_fid_Q.SetBinContent(i+1,chi2_true_fid_Q);
	hChi2_reco_fid_Q.SetBinContent(i+1,chi2_reco_fid_Q);
  }
  */
  //hDiff.Draw(); 
  //c1.SaveAs(path + "result.pdf","pdf");

  //hChi2.Draw("");
  std::cout << "chi2: " << chi2 << std::endl;
  std::cout << "prob: " << TMath::Prob(chi2,hNom.GetNbinsX()) << std::endl;
  //c1.SaveAs(path + "result.pdf","pdf");
  
  //hChi2_half.Draw("");
  std::cout << "chi2_half: " << chi2_half  << std::endl;
  std::cout << "prob_half: " << TMath::Prob(chi2_half,hNomHalf1.GetNbinsX()) << std::endl;
  //c1.SaveAs(path + "result.pdf","pdf");

  //hChi2_true.Draw("");
  std::cout << "chi2_true: " << chi2_true << std::endl;
  std::cout << "prob_true: " << TMath::Prob(chi2_true,hNom_true.GetNbinsX()) << std::endl;
  
  //c1.SaveAs(path + "result.pdf","pdf");
  
  std::cout << "chi2_true_all: " << chi2_true_all << std::endl;
  std::cout << "prob_true_all: " << TMath::Prob(chi2_true_all,hNom_true_all.GetNbinsX()) << std::endl;

  //hChi2_true_fid.Draw("");
  std::cout << "chi2_true_fid: " << chi2_true_fid << std::endl;
  std::cout << "prob_true_fid: " << TMath::Prob(chi2_true_fid,hNom_true_fid.GetNbinsX()) << std::endl;
  //c1.SaveAs(path + "result.pdf","pdf");
  
  //hChi2_true_fid_QE.Draw("");
  //std::cout << "chi2_true_fid_QE: " << chi2tot_true_fid_QE << std::endl;
  //std::cout << "prob_true_fid_QE: " << TMath::Prob(chi2tot_true_fid_QE,hChi2_true_fid_QE.GetNbinsX()) << std::endl;
  //c1.SaveAs(path + "result.pdf)","pdf");
  
  //hChi2_true_fid_Q.Draw("");
  std::cout << "chi2_true_fid_Q: " << chi2_true_fid_Q << std::endl;
  std::cout << "prob_true_fid_Q: " << TMath::Prob(chi2_true_fid_Q,hNom_true_fid_Q.GetNbinsX()) << std::endl;
  //c1.SaveAs(path + "result.pdf)","pdf");
  
  //hChi2_reco_fid_Q.Draw("");
  std::cout << "chi2_reco_fid_Q: " << chi2_reco_fid_Q << std::endl;
  std::cout << "prob_reco_fid_Q: " << TMath::Prob(chi2_reco_fid_Q,hNom_reco_fid_Q.GetNbinsX()) << std::endl;
  //c1.SaveAs(path + "result.pdf)","pdf");
  
 
  /////////////////////////////////////////////////////KOLMOGOROV AND ANDERSON//////////////////////////////////////////////////////////////////
  TH1D* hNomb = new TH1D("hNomb","",1000,0,16);
  TH1D* hH1Yb = new TH1D("hH1Yb","",1000,0,16);  
  TH1D* hNomHalf1b = new TH1D("hNomHalf1b","",1000,0,16);
  TH1D* hNomHalf2b = new TH1D("hNomHalf2b","",1000,0,16);
  TH1D* hNom_trueb = new TH1D("hNom_trueb","",1000,0,16);
  TH1D* hH1Y_trueb = new TH1D("hH1Y_trueb","",1000,0,16);
  TH1D* hNom_true_fidb = new TH1D("hNom_true_fidb","",1000,0,16);
  TH1D* hH1Y_true_fidb = new TH1D("hH1Y_true_fidb","",1000,0,16);
  TH1D* hNom_true_allb = new TH1D("hNom_true_allb","",1000,0,16);
  TH1D* hH1Y_true_allb = new TH1D("hH1Y_true_allb","",1000,0,16);
  TH1D* hNom_true_fid_QEb  = new TH1D("hNom_true_fid_QEb","",1000,0,16);
  TH1D* hH1Y_true_fid_QEb  = new TH1D("hH1Y_true_fid_QEb","",1000,0,16);
  TH1D* hNom_true_fid_Qb  = new TH1D("hNom_true_fid_Qb","",1000,0,16);
  TH1D* hH1Y_true_fid_Qb  = new TH1D("hH1Y_true_fid_Qb","",1000,0,16);
  TH1D* hNom_reco_fid_Qb  = new TH1D("hNom_reco_fid_Qb","",1000,0,16);
  TH1D* hH1Y_reco_fid_Qb  = new TH1D("hH1Y_reco_fid_Qb","",1000,0,16);
  
  tNom.Draw("pmu_reco/1E3>>hNomb",fiducial + CC + muonOK + lowE,"E0");
  tH1Y.Draw("pmu_reco/1E3>>hH1Yb",fiducial + CC + muonOK + lowE,"E0",NlowEok,0);
  tNom.Draw("pmu_reco/1E3>>hNomHalf1b",fiducial + CC + muonOK + lowE,"E0",tNom.GetEntries()*0.5,0);
  tNom.Draw("pmu_reco/1E3>>hNomHalf2b",fiducial + CC + muonOK + lowE,"E0",tNom.GetEntries()*0.5,tNom.GetEntries()*0.5);
  tNom.Draw("pmu_true/1E3>>hNom_trueb",fiducial + CC + muonOK + lowE,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_trueb",fiducial + CC + muonOK + lowE,"E0",NlowEok,0);
  tNom.Draw("pmu_true/1E3>>hNom_true_fidb",fiducial + CC + lowE,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_fidb",fiducial + CC + lowE,"E0",NlowEok,0);
  tNom.Draw("pmu_true/1E3>>hNom_true_allb",fiducial_truth + CC + lowE,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_allb",fiducial_truth + CC + lowE,"E0",NlowEok,0);
  tNom.Draw("pmu_reco/1E3>>hNom_true_fid_QEb",fiducial + CC + lowE + muonOK + QE,"E0");
  tH1Y.Draw("pmu_reco/1E3>>hH1Y_true_fid_QEb",fiducial + CC + lowE + muonOK + QE,"E0",NlowEok,0);
  tNom.Draw("pmu_true/1E3>>hNom_true_fid_Qb",fiducial + CC + lowE + muonOK + QualityCut,"E0");
  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_fid_Qb",fiducial + CC + lowE + muonOK + QualityCut,"E0",NlowEok,0);
  tNom.Draw("pmu_reco/1E3>>hNom_reco_fid_Qb",fiducial + CC + lowE + muonOK + QualityCut,"E0");
  tH1Y.Draw("pmu_reco/1E3>>hH1Y_reco_fid_Qb",fiducial + CC + lowE + muonOK + QualityCut,"E0",NlowEok,0);
/*  
  TH1* hNomc = hNomb.GetCumulative();
  TH1* hH1Yc = hH1Yb.GetCumulative();
  TH1* hNomHalf1c = hNomHalf1b.GetCumulative();
  TH1* hNomHalf2c = hNomHalf2b.GetCumulative();
  TH1* hNom_truec = hNom_trueb.GetCumulative();
  TH1* hH1Y_truec = hH1Y_trueb.GetCumulative();
  TH1* hNom_true_fidc = hNom_true_fidb.GetCumulative();
  TH1* hH1Y_true_fidc = hH1Y_true_fidb.GetCumulative();
  
  hNomc->Scale(1/(hNomc->Integral()));
  hH1Yc->Scale(1/(hH1Yc->Integral()));
  hNomHalf1c->Scale(1/(hNomHalf1c->Integral()));
  hNomHalf2c->Scale(1/(hNomHalf2c->Integral()));
  hNom_truec->Scale(1/(hNom_truec->Integral()));
  hH1Y_truec->Scale(1/(hH1Y_truec->Integral()));
  hNom_true_fidc->Scale(1/(hNom_true_fidc->Integral()));
  hH1Y_true_fidc->Scale(1/(hH1Y_true_fidc->Integral()));
  

  double k = 0.;
  double k_half = 0.;
  double k_true = 0.;
  double k_true_fid = 0.;

  for(int i = 0; i < hNomc->GetNbinsX(); i++)
  {
	double diff = std::abs(hNomc->GetBinContent(i+1)-hH1Yc->GetBinContent(i+1));
	double diff_half = std::abs(hNomHalf1c->GetBinContent(i+1)-hNomHalf2c->GetBinContent(i+1));
	double diff_true = std::abs(hNom_truec->GetBinContent(i+1)-hH1Y_truec->GetBinContent(i+1));
	double diff_true_fid = std::abs(hNom_true_fidc->GetBinContent(i+1)-hH1Y_true_fidc->GetBinContent(i+1));
    
	if(diff > k) {k = diff;}
	if(diff_half > k_half) {k_half = diff_half;}
	if(diff_true > k_true) {k_true = diff_true;}
	if(diff_true_fid > k_true_fid) {k_true_fid = diff_true_fid;}
  }
  
  double kexp = -2*k*k * ((hNomb.GetEntries()*hH1Yb.GetEntries())/(hNomb.GetEntries()+hH1Yb.GetEntries()));
  double k_halfexp = -2*k_half*k_half * ((hNomHalf1b.GetEntries()*hNomHalf2b.GetEntries())/(hNomHalf1b.GetEntries()+hNomHalf2b.GetEntries()));
  double k_trueexp = -2*k_true*k_true * ((hNom_trueb.GetEntries()*hH1Y_trueb.GetEntries())/(hNom_trueb.GetEntries()+hH1Y_trueb.GetEntries()));
  double k_true_fidexp = - 2 * k_true_fid * k_true_fid * ((hNom_true_fidb.GetEntries()*hH1Y_true_fidb.GetEntries())/(hNom_true_fidb.GetEntries()+hH1Y_true_fidb.GetEntries()));
  
  double pk = 2* std::exp(kexp);
  double pk_half = 2* std::exp(k_halfexp);
  double pk_true = 2* std::exp(k_trueexp);
  double pk_true_fid = 2* std::exp(k_true_fidexp);

  hNomc->SetMarkerColor(kBlack);
  hNomc->SetLineColor(kBlack);
  hH1Yc->SetMarkerColor(kRed);
  hH1Yc->SetLineColor(kRed);
  hNomc->Draw("");
  hH1Yc->Draw("same");
  c1.SaveAs(path + "result.pdf)","pdf");
 */ 
  double k = hNomb->KolmogorovTest(hH1Yb,"M");
  double pk = hNomb->KolmogorovTest(hH1Yb,"");
  double k_half = hNomHalf1b->KolmogorovTest(hNomHalf2b,"M");
  double pk_half = hNomHalf1b->KolmogorovTest(hNomHalf2b,"");
  double k_true = hNom_trueb->KolmogorovTest(hH1Y_trueb,"M");
  double pk_true = hNom_trueb->KolmogorovTest(hH1Y_trueb,"");
  double k_true_fid = hNom_true_fidb->KolmogorovTest(hH1Y_true_fidb,"M");
  double pk_true_fid = hNom_true_fidb->KolmogorovTest(hH1Y_true_fidb,"");
  double k_true_all = hNom_true_allb->KolmogorovTest(hH1Y_true_allb,"M");
  double pk_true_all = hNom_true_allb->KolmogorovTest(hH1Y_true_allb,"");
  double k_true_fid_QE = hNom_true_fid_QEb->KolmogorovTest(hH1Y_true_fid_QEb,"M");
  double pk_true_fid_QE = hNom_true_fid_QEb->KolmogorovTest(hH1Y_true_fid_QEb,"");
  double k_true_fid_Q = hNom_true_fid_Qb->KolmogorovTest(hH1Y_true_fid_Qb,"M");
  double pk_true_fid_Q = hNom_true_fid_Qb->KolmogorovTest(hH1Y_true_fid_Qb,"");
  double k_reco_fid_Q = hNom_reco_fid_Qb->KolmogorovTest(hH1Y_reco_fid_Qb,"M");
  double pk_reco_fid_Q = hNom_reco_fid_Qb->KolmogorovTest(hH1Y_reco_fid_Qb,"");

  
  std::cout << "kolmogorov test: " << k << "; n:" << hNom.GetEntries() << "; m:" << hH1Y.GetEntries() <<"; prob:" << pk << "; sigma:" << TMath::ErfInverse(1-pk)*sqrt(2) << std::endl;
  std::cout << "kolmogorov test (half): " << k_half << "; n:" << hNomHalf1.GetEntries() << "; m:" << hNomHalf2.GetEntries() << "; prob:" << pk_half << "; sigma:" << TMath::ErfInverse(1-pk_half)*sqrt(2) <<std::endl;
  std::cout << "kolmogorov test (true): " << k_true << "; n:" << hNom_true.GetEntries() << "; m:" << hH1Y_true.GetEntries() << "; prob:" << pk_true <<"; sigma:" << TMath::ErfInverse(1-pk_true)*sqrt(2) << std::endl;
  std::cout << "kolmogorov test (true_fid): " << k_true_fid << "; n:" << hNom_true_fid.GetEntries() << "; m:" << hH1Y_true_fid.GetEntries() << "; prob:" << pk_true_fid <<"; sigma:" << TMath::ErfInverse(1-pk_true_fid)*sqrt(2) << std::endl;
  std::cout << "kolmogorov test (true_all): " << k_true_all << "; n:" << hNom_true_all.GetEntries() << "; m:" << hH1Y_true_all.GetEntries() << "; prob:" << pk_true_all <<"; sigma:" << TMath::ErfInverse(1-pk_true_all)*sqrt(2) << std::endl;
  std::cout << "kolmogorov test (true_fid_QE): " << k_true_fid_QE << "; n:" << hNom_true_fid_QE.GetEntries() << "; m:" << hH1Y_true_fid_QE.GetEntries() << "; prob:" << pk_true_fid_QE << "; sigma:" << TMath::ErfInverse(1-pk_true_fid_QE)*sqrt(2) <<std::endl;
  std::cout << "kolmogorov test (true_fid_Q): " << k_true_fid_Q << "; n:" << hNom_true_fid_Q.GetEntries() << "; m:" << hH1Y_true_fid_Q.GetEntries() << "; prob:" << pk_true_fid_Q << "; sigma:" << TMath::ErfInverse(1-pk_true_fid_Q)*sqrt(2) <<std::endl;
  std::cout << "kolmogorov test (reco_fid_Q): " << k_reco_fid_Q << "; n:" << hNom_reco_fid_Q.GetEntries() << "; m:" << hH1Y_reco_fid_Q.GetEntries() << "; prob:" << pk_reco_fid_Q << "; sigma:" << TMath::ErfInverse(1-pk_reco_fid_Q)*sqrt(2) <<std::endl;
  
  double a = hNomb->AndersonDarlingTest(hH1Yb,"T");
  double pa = hNomb->AndersonDarlingTest(hH1Yb,"");
  double a_half = hNomHalf1b->AndersonDarlingTest(hNomHalf2b,"T");
  double pa_half = hNomHalf1b->AndersonDarlingTest(hNomHalf2b,"");
  double a_true = hNom_trueb->AndersonDarlingTest(hH1Y_trueb,"T");
  double pa_true = hNom_trueb->AndersonDarlingTest(hH1Y_trueb,"");
  double a_true_fid = hNom_true_fidb->AndersonDarlingTest(hH1Y_true_fidb,"T");
  double pa_true_fid = hNom_true_fidb->AndersonDarlingTest(hH1Y_true_fidb,"");
  double a_true_all = hNom_true_allb->AndersonDarlingTest(hH1Y_true_allb,"T");
  double pa_true_all = hNom_true_allb->AndersonDarlingTest(hH1Y_true_allb,"");
  double a_true_fid_QE = hNom_true_fid_QEb->AndersonDarlingTest(hH1Y_true_fid_QEb,"T");
  double pa_true_fid_QE = hNom_true_fid_QEb->AndersonDarlingTest(hH1Y_true_fid_QEb,"");
  double a_true_fid_Q = hNom_true_fid_Qb->AndersonDarlingTest(hH1Y_true_fid_Qb,"T");
  double pa_true_fid_Q = hNom_true_fid_Qb->AndersonDarlingTest(hH1Y_true_fid_Qb,"");
  double a_reco_fid_Q = hNom_reco_fid_Qb->AndersonDarlingTest(hH1Y_true_fid_Qb,"T");
  double pa_reco_fid_Q = hNom_reco_fid_Qb->AndersonDarlingTest(hH1Y_reco_fid_Qb,"");

  
  std::cout << "Anderson test: " << a << "; n:" << hNom.GetEntries() << "; m:" << hH1Y.GetEntries() <<"; prob:" << pa << "; sigma:" << TMath::ErfInverse(1-pa)*sqrt(2) <<std::endl;
  std::cout << "Anderson test (half): " << a_half << "; n:" << hNomHalf1.GetEntries() << "; m:" << hNomHalf2.GetEntries() << "; prob:" << pa_half <<"; sigma:" << TMath::ErfInverse(1-pa_half)*sqrt(2) << std::endl;
  std::cout << "Anderson test (true): " << a_true << "; n:" << hNom_true.GetEntries() << "; m:" << hH1Y_true.GetEntries() << "; prob:" << pa_true <<"; sigma:" << TMath::ErfInverse(1-pa_true)*sqrt(2) << std::endl;
  std::cout << "Anderson test (true_fid): " << a_true_fid << "; n:" << hNom_true_fid.GetEntries() << "; m:" << hH1Y_true_fid.GetEntries() << "; prob:" << pa_true_fid <<"; sigma:" << TMath::ErfInverse(1-pa_true_fid)*sqrt(2) << std::endl;
  std::cout << "Anderson test (true_all): " << a_true_all << "; n:" << hNom_true_all.GetEntries() << "; m:" << hH1Y_true_all.GetEntries() << "; prob:" << pa_true_all <<"; sigma:" << TMath::ErfInverse(1-pa_true_all)*sqrt(2) << std::endl;
  std::cout << "Anderson test (true_fid_QE): " << a_true_fid_QE << "; n:" << hNom_true_fid_QE.GetEntries() << "; m:" << hH1Y_true_fid_QE.GetEntries() << "; prob:" << pa_true_fid_QE <<"; sigma:" << TMath::ErfInverse(1-pa_true_fid_QE)*sqrt(2) << std::endl;
  std::cout << "Anderson test (true_fid_Q): " << a_true_fid_Q << "; n:" << hNom_true_fid_Q.GetEntries() << "; m:" << hH1Y_true_fid_Q.GetEntries() << "; prob:" << pa_true_fid_Q <<"; sigma:" << TMath::ErfInverse(1-pa_true_fid_Q)*sqrt(2) << std::endl;
  std::cout << "Anderson test (reco_fid_Q): " << a_reco_fid_Q << "; n:" << hNom_reco_fid_Q.GetEntries() << "; m:" << hH1Y_reco_fid_Q.GetEntries() << "; prob:" << pa_reco_fid_Q <<"; sigma:" << TMath::ErfInverse(1-pa_reco_fid_Q)*sqrt(2) <<std::endl;
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*
  int nEvents[] = {1000, 10000, 50000, 100000, 500000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1424112};
  int nn = sizeof(nEvents)/sizeof(int);
  
  TH1D* h1[nn];
  TH1D* h2[nn]; 
  
  TGraph grChi2(nn);
  grChi2.SetMarkerStyle(kPlus);
  grChi2.SetMarkerSize(1.5);
  grChi2.SetMarkerColor(kRed);
  
  TGraph grPVal(nn);
  grPVal.SetMarkerStyle(kPlus);
  grPVal.SetMarkerSize(1.5);
  grPVal.SetMarkerColor(kRed);
  
  TGraph grSens(nn);
  grSens.SetMarkerStyle(kPlus);
  grSens.SetMarkerSize(1.5);
  grSens.SetMarkerColor(kRed);
  
  double chi;
  int n;
  
  for(int i = 0; i < nn; i++)
  {
    h1[i] = new TH1D("","",40,0,16);
    h2[i] = new TH1D("","",40,0,16);
    
    FillHisto(tNom, tH1Y, *(h1[i]), *(h2[i]), nEvents[i]);
    
    EvalChi2(*(h1[i]), *(h2[i]), chi, n);
    
    double p_value = TMath::Prob(chi,n);
    double nsigmas = TMath::ErfInverse(1-p_value)*sqrt(2);

    std::cout << "Nevents: " << nEvents[i] << std::endl;
    std::cout << "chi2   : " << chi << std::endl;
    std::cout << "prob   : " << p_value << std::endl;
    std::cout << "Nsigma : " << nsigmas << std::endl;
    
    grChi2.SetPoint(i, nEvents[i], chi);
    grPVal.SetPoint(i, nEvents[i], p_value);
    grSens.SetPoint(i, nEvents[i], nsigmas);
  }
  
  c1.cd();
  grChi2.Draw("APL");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SetLogy(true);
  grPVal.Draw("APL");
  c1.SaveAs(path + "result.pdf","pdf");
  c1.SetLogy(false);
  grSens.Draw("APL");
  c1.SaveAs(path + "result.pdf)","pdf");
  */ 
}

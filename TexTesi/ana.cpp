
TString path = "/wd/dune-it/beam_monitoring/script";

TFile fout(path + "/../result.root","RECREATE");

double kloe_center[] = {0.0000000 * 10, -238.47300 * 10, 2391.0000 * 10};

double emcalo_barrel_mod_halfw =  29.285 * 10;
double emcalo_barrel_mod_halfl = 215. * 10;
double emcalo_barrel_mod_h = 11.5 * 2 * 10;
double emcalo_int_radius = 2000;
double emcalo_ext_radius = TMath::Sqrt((emcalo_int_radius+emcalo_barrel_mod_h)*(emcalo_int_radius+emcalo_barrel_mod_h)+emcalo_barrel_mod_halfw*emcalo_barrel_mod_halfw);

double threshold = 40.;

TCut fiducial = TString::Format("maxde_frontoutlayer < %f && abs(xv_reco) < 1500 ",threshold).Data();
TCut CC = "isCC == 1";
TCut muonOK = "isPmuOK == 1 && !(pxmu_reco == 0 && pymu_reco == 0 && pzmu_reco == 0)";
TCut QE = "NHitLayer[0] == 1 ||  NHitLayer[0] == 2";
TCut lowE = "Enu_true/1E3<20";
TCut QualityCut = "red_chi2_ln < 3 && red_chi2_cr < 1000";

TCanvas c1("","",1000,1000);

void processChain(TChain& t)
{
  // all yx vertex position
  c1.DrawFrame(kloe_center[0]-3000,kloe_center[1]-3000,kloe_center[0]+3000,kloe_center[1]+3000);

  t.Draw("yv_true:xv_true","","same",10000,0);

  TBox b(kloe_center[0]-emcalo_barrel_mod_halfl,kloe_center[1]-emcalo_ext_radius,kloe_center[0]+emcalo_barrel_mod_halfl,kloe_center[1]+emcalo_ext_radius);
  b.SetFillStyle(0);
  b.Draw();

  c1.SaveAs(path + "/../result.pdf(","pdf");

  // all yz vertex position
  c1.DrawFrame(kloe_center[2]-3000,kloe_center[1]-3000,kloe_center[2]+3000,kloe_center[1]+3000);

  t.Draw("yv_true:zv_true","","same",10000,0);

  TEllipse elint(kloe_center[2],kloe_center[1],emcalo_int_radius,emcalo_int_radius);
  TEllipse elext(kloe_center[2],kloe_center[1],emcalo_ext_radius,emcalo_ext_radius);
  elint.SetFillStyle(0);
  elext.SetFillStyle(0);
  elint.Draw();
  elext.Draw();

  c1.SaveAs(path + "/../result.pdf","pdf");
  
  // int and ext yx vertex position
  c1.DrawFrame(kloe_center[0]-3000,kloe_center[1]-3000,kloe_center[0]+3000,kloe_center[1]+3000);
  t.SetMarkerColor(kRed);
  t.Draw("yv_true:xv_true",!fiducial,"same",10000,0);
  t.SetMarkerColor(kBlack);
  t.Draw("yv_true:xv_true",fiducial,"same",10000,0);

  c1.SaveAs(path + "/../result.pdf","pdf");

  // int and ext yz vertex position
  c1.DrawFrame(kloe_center[2]-3000,kloe_center[1]-3000,kloe_center[2]+3000,kloe_center[1]+3000);
  t.SetMarkerColor(kRed);
  t.Draw("yv_true:zv_true",!fiducial,"same",10000,0);
  t.SetMarkerColor(kBlack);
  t.Draw("yv_true:zv_true",fiducial,"same",10000,0);

  c1.SaveAs(path + "/../result.pdf","pdf");
 
  // neutrino energy
  TH1D h_Enu("h_Enu","",80,0,40);
  long int nAll = t.Draw("Enu_true/1E3>>h_Enu","","");
  h_Enu.Draw();

  c1.SaveAs(path + "/../result.pdf","pdf");

  // fiducial vollume selection efficiency
  TH1D h_Enu_int("h_Enu_int","",80,0,40);
  long int nFiducial = t.Draw("Enu_true/1E3>>h_Enu_int",fiducial,"");

  TEfficiency eff_fid(h_Enu_int,h_Enu);
  eff_fid.Draw("");

  c1.SaveAs(path + "/../result.pdf","pdf");

  std::cout << "Fiducial Volume selection efficiency: " << double(nFiducial)/double(nAll) << std::endl;

  long int nCC = t.GetEntries(fiducial + CC);

  std::cout << "CC fraction: " << double(nCC)/double(nFiducial) << std::endl;

  TH1D h_pmu("h_pmu","",80,0,40);
  long int nintmu = t.Draw("pmu_true/1E3>>h_pmu",fiducial + CC,"");
  h_pmu.Draw();

  c1.SaveAs(path + "/../result.pdf","pdf");

  TH1D h_pmu_ok("h_pmu_ok","",80,0,40);
  long int nintmu_ok = t.Draw("pmu_true/1E3>>h_pmu_ok",fiducial + CC + muonOK,"");

  TEfficiency eff_murec(h_pmu_ok,h_pmu);
  eff_murec.Draw("");

  c1.SaveAs(path + "/../result.pdf","pdf");

  std::cout << "Reconstucted muon: " << double(nintmu_ok)/double(nintmu) << std::endl;

  TH2D h_p("h_p","",200,0,40,200,0,40);
  t.Draw("pmu_reco/1E3:pmu_true/1E3>>h_p",fiducial + CC + muonOK,"");
  h_p.Draw("colz");

  c1.SaveAs(path + "/../result.pdf","pdf");
  
  TH2D h_reldpt_chi2("h_reldpt_chi2","",100,0,1E6,100,-1,1);
  t.Draw("reldptmu:chi2_cr>>h_reldpt_chi2",fiducial + CC + muonOK,"");
  h_reldpt_chi2.Draw("colz");

  c1.SaveAs(path + "/../result.pdf","pdf");
  
  TH2D h_reldpt_redchi2("h_reldpt_redchi2","",100,0,10000,100,-1,1);
  t.Draw("reldptmu:red_chi2_cr>>h_reldpt_redchi2",fiducial + CC + muonOK,"");
  h_reldpt_redchi2.Draw("colz");

  c1.SaveAs(path + "/../result.pdf","pdf");
  
  TH2D h_reldpl_chi2("h_reldpl_chi2","",100,0,100,100,-1,1);
  t.Draw("reldplmu:chi2_ln>>h_reldpl_chi2",fiducial + CC + muonOK,"");
  h_reldpl_chi2.Draw("colz");

  c1.SaveAs(path + "/../result.pdf","pdf");
  
  TH2D h_reldpl_redchi2("h_reldpl_redchi2","",100,0,10,100,-1,1);
  t.Draw("reldplmu:red_chi2_ln>>h_reldpl_redchi2",fiducial + CC + muonOK,"");
  h_reldpl_redchi2.Draw("colz");

  c1.SaveAs(path + "/../result.pdf","pdf");  
  
  long int nintmu_ok_qc = t.GetEntries(fiducial + CC + muonOK + QualityCut);
  
  std::cout << "Quality cut: " << double(nintmu_ok_qc)/double(nintmu_ok) << std::endl;

  TH1D h_pmu_1D("h_pmu_1D","",80,0,40);
  t.Draw("pmu_reco/1E3>>h_pmu_1D",fiducial + CC + muonOK,"E0");
  t.SetLineColor(kRed);
  t.SetMarkerColor(kRed);
  t.Draw("pmu_true/1E3",fiducial + CC + muonOK,"same");
  t.SetLineColor(kBlack);
  t.SetMarkerColor(kBlack);

  c1.SaveAs(path + "/../result.pdf","pdf");  

}

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

  int long n1 = t1.Draw(TString::Format("pmu_reco/1E3>>%s",name1.Data()).Data(),fiducial + CC + muonOK + QualityCut,"E0", nEv);
  int long n2 = t2.Draw(TString::Format("pmu_reco/1E3>>%s",name2.Data()).Data(),fiducial + CC + muonOK + QualityCut,"E0", nEv);
  
  std::cout << "nEv1: " << nEv << " -> " << n1 << " " << std::endl;
  std::cout << "nEv2: " << nEv << " -> " << n2 << " " << std::endl;
}

void ana()
{
  gROOT->SetBatch();

  ///////////////////////////////////////////////////////////////////////////
  // Read chain
  TChain tNom("tout","tNom");
  TChain tH1Y("tout","tH1Y");

  tNom.Add(path + "/../nominal/*");
  tH1Y.Add(path + "/../H1ShiftY/*");
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // Define alias
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

  tNom.SetAlias("red_chi2_ln","2*chi2_ln/nHit");
  tH1Y.SetAlias("red_chi2_ln","2*chi2_ln/nHit");
  ///////////////////////////////////////////////////////////////////////////
  
  
  ///////////////////////////////////////////////////////////////////////////
  // Some plots
  processChain(tNom);
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // Ensuring TChains have the same event number
  int long NlowEok = tNom.GetEntries(lowE);
  
  tNom.Draw(">>entrylistNom",lowE,"entrylist");
  
  TEntryList *entrylistNom = (TEntryList*)gDirectory->Get("entrylistNom");
  tNom.SetEntryList(entrylistNom);
  
  tH1Y.Draw(">>entrylistH1Y",lowE,"entrylist",NlowEok,0);
  
  TEntryList *entrylistH1Y = (TEntryList*)gDirectory->Get("entrylistH1Y");
  tH1Y.SetEntryList(entrylistH1Y);
  
  TH1D hEnuNom("hEnuNom","",100,0,20);
  TH1D hEnuH1Y("hEnuH1Y","",100,0,20);

  int long n1 = tNom.Draw("Enu_true/1E3>>hEnuNom","","E0");
  int long n2 = tH1Y.Draw("Enu_true/1E3>>hEnuH1Y","","E0");

  if(n1 != n2)
  {
    std::cout << "error: samples are different" << n1 << " != " << n2 << std::endl;
    return;
  }
  
  std::cout << "Total number of events (CC+NC): " << n1 << std::endl;
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Comparing neutrino energy
  hEnuNom.SetLineColor(kBlack);
  hEnuH1Y.SetLineColor(kRed);

  c1.cd();
  hEnuNom.Draw();
  hEnuH1Y.Draw("same");
  c1.SaveAs(path + "/../result.pdf","pdf");
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Reconstructed Pmu
  TH1D hNom("hNom","",40,0,16);
  TH1D hH1Y("hH1Y","",40,0,16);

  tNom.Draw("pmu_reco/1E3>>hNom",fiducial + CC + muonOK + QualityCut,"E0");
  hNom.SetMarkerColor(kBlack);
  hNom.SetLineColor(kBlack);

  tH1Y.Draw("pmu_reco/1E3>>hH1Y",fiducial + CC + muonOK + QualityCut,"E0");
  hH1Y.SetMarkerColor(kRed);
  hH1Y.SetLineColor(kRed);

  c1.SetLogy(true);
  hNom.Draw();
  hH1Y.Draw("same");
  c1.SaveAs(path + "/../result.pdf","pdf");

  c1.SetLogy(false);
  
  double chi2;
  int ndof;
  
  EvalChi2(hNom, hH1Y, chi2, ndof);

  std::cout << "Reco Pmu" << std::endl;
  std::cout << "chi2: " << chi2 << std::endl;
  std::cout << "prob: " << TMath::Prob(chi2,ndof) << std::endl;
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // True Pmu
  TH1D hNom_true("hNom_true","",40,0,16);
  TH1D hH1Y_true("hH1Y_true","",40,0,16);

  tNom.Draw("pmu_true/1E3>>hNom_true",fiducial + CC + muonOK + QualityCut,"E0");
  hNom_true.SetMarkerColor(kBlack);
  hNom_true.SetLineColor(kBlack);

  tH1Y.Draw("pmu_true/1E3>>hH1Y_true",fiducial + CC + muonOK + QualityCut,"E0");
  hH1Y_true.SetMarkerColor(kRed);
  hH1Y_true.SetLineColor(kRed);

  c1.SetLogy(true);
  hNom_true.Draw();
  hH1Y_true.Draw("same");
  c1.SaveAs(path + "/../result.pdf","pdf");

  c1.SetLogy(false);
  
  EvalChi2(hNom_true, hH1Y_true, chi2, ndof);

  std::cout << "True Pmu" << std::endl;
  std::cout << "chi2: " << chi2 << std::endl;
  std::cout << "prob: " << TMath::Prob(chi2,ndof) << std::endl;
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // True Pmu fiducial
  TH1D hNom_true_fid("hNom_true_fid","",40,0,16);
  TH1D hH1Y_true_fid("hH1Y_true_fid","",40,0,16);

  tNom.Draw("pmu_true/1E3>>hNom_true_fid",fiducial + CC,"E0");
  hNom_true_fid.SetMarkerColor(kBlack);
  hNom_true_fid.SetLineColor(kBlack);

  tH1Y.Draw("pmu_true/1E3>>hH1Y_true_fid",fiducial + CC,"E0");
  hH1Y_true_fid.SetMarkerColor(kRed);
  hH1Y_true_fid.SetLineColor(kRed);

  c1.SetLogy(true);
  hNom_true_fid.Draw();
  hH1Y_true_fid.Draw("same");
  c1.SaveAs(path + "/../result.pdf","pdf");

  c1.SetLogy(false);
  
  EvalChi2(hNom_true_fid, hH1Y_true_fid, chi2, ndof);

  std::cout << "All fiducial volume true Pmu" << std::endl;
  std::cout << "chi2: " << chi2 << std::endl;
  std::cout << "prob: " << TMath::Prob(chi2,ndof) << std::endl;
  ///////////////////////////////////////////////////////////////////////////
  
  
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
  
  for(int i = 0; i < nn; i++)
  {
    h1[i] = new TH1D("","",40,0,16);
    h2[i] = new TH1D("","",40,0,16);
    
    FillHisto(tNom, tH1Y, *(h1[i]), *(h2[i]), nEvents[i]);
    
    EvalChi2(*(h1[i]), *(h2[i]), chi2, ndof);
    
    double p_value = TMath::Prob(chi2,ndof);
    double nsigmas = TMath::ErfInverse(1-p_value)*sqrt(2);

    std::cout << "Nevents: " << nEvents[i] << std::endl;
    std::cout << "chi2   : " << chi2 << std::endl;
    std::cout << "prob   : " << p_value << std::endl;
    std::cout << "Nsigma : " << nsigmas << std::endl;
    
    grChi2.SetPoint(i, nEvents[i], chi2);
    grPVal.SetPoint(i, nEvents[i], p_value);
    grSens.SetPoint(i, nEvents[i], nsigmas);
  }
  
  c1.cd();
  grChi2.Draw("APL");
  c1.SaveAs(path + "/../result.pdf","pdf");
  c1.SetLogy(true);
  grPVal.Draw("APL");
  c1.SaveAs(path + "/../result.pdf","pdf");
  c1.SetLogy(false);
  grSens.Draw("APL");
  c1.SaveAs(path + "/../result.pdf)","pdf");
}


void PhotonTimePlotter(bool PrintPlot=true)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  vector<TString> label = {
    "1_LASC_10p5x10p5_mm2",
    "4_LASC_4x4_mm2",
    "4_LASC_4x4_mm2_spec_refl_in_gaps"};

  //1 filename per chain!!!!
  vector<TString> filename = {
    "/eos/user/f/fmonti/SiPM_Optical_Simulation/LASC_10p5x10p5_TILE_10p5x10p5/surface5_refl97_Tile11_t*.5x11.5x4x_planar_geometry_0.24844x0.24844x0.8_tilt0_PDE_sourceuniform_seed*.root",
    "/eos/user/f/fmonti/SiPM_Optical_Simulation/4LASC_4x4_TILE_11p5x11p5/surface5_refl97_Tile11_t*.5x11.5x4x_planar_geometry_0.20567x0.20567x0.8_tilt0_PDE_sourceuniform_seed*.root",
    "/eos/user/f/fmonti/SiPM_Optical_Simulation/4LASC_4x4_TILE_11p5x11p5_allbackaluminum/surface5_refl97_Tile11_t*.5x11.5x4x_planar_geometry_0.20567x0.20567x0.8_tilt0_PDE_sourceuniform_seed*.root"};
  assert(label.size()==filename.size());

  map<TString,TChain*> ch;
  for(int i=0;i<label.size();++i)
  {
    ch[label[i]] = new TChain("ntu","ntu");
    ch[label[i]] -> Add(filename[i].Data());
  }

  /*
  TChain *ch_LASC = new TChain("ntu","ntu");
  ch_LASC->Add("/eos/user/f/fmonti/SiPM_Optical_Simulation/LASC/surface5_refl97_Tile11.5x11.5x4x_planar_geometry_0.264x0.264x0.8_tilt0_PDE_sourceuniform_seed56.root");

  TChain *ch_LASC_nowrap = new TChain("ntu","ntu");
  ch_LASC_nowrap->Add("/eos/user/f/fmonti/SiPM_Optical_Simulation/LASC/surface0_Tile11.5x11.5x4x_planar_geometry_0.264x0.264x0.8_tilt0_PDE_sourceuniform_seed56.root");

  TChain *ch_LASC_norefl = new TChain("ntu","ntu");
  ch_LASC_norefl->Add("/eos/user/f/fmonti/SiPM_Optical_Simulation/LASC/surface5_refl97_Tile11.5x11.5x4x_planar_geometry_10.5x10.5x0.8_tilt0_PDE_sourceuniform_seed56.root");

  TChain *ch_STD = new TChain("ntu","ntu");
  ch_STD->Add("/eos/user/f/fmonti/SiPM_Optical_Simulation/generic_sipm_4x4/surface5_refl97_Tile11.5x11.5x4x_planar_geometry_4x4x0.8_tilt0_PDE_sourceuniform_seed542.root");
  */

  //draw Ncollected photons VS impact point
  TCanvas *c_amp_max = new TCanvas();
  map<TString,TProfile*> p_amp_x; 
  for(int i=0;i<label.size();++i)
  {
    ch[label[i]] -> Draw(Form("@Ph_time.size(Ph_time<=10):mu_x_hit>>p_%s(22,-5.5,5.5)",label[i].Data()),"abs(mu_y_hit)>1.05 && abs(mu_y_hit)<4.45","proff");
    p_amp_x[label[i]] = (TProfile*)gDirectory->Get(Form("p_%s",label[i].Data()));
    p_amp_x[label[i]]->SetMarkerStyle(20);
    p_amp_x[label[i]]->SetMarkerColor(i+1);
    p_amp_x[label[i]]->SetTitle(label[i]);
    if(PrintPlot)
      c_amp_max->Print(Form("p_amp_x_%s.png",label[i].Data()));
  }
  p_amp_x[label[0]]->Draw();
  p_amp_x[label[0]]->GetYaxis()->SetRangeUser(2000,3000);
  p_amp_x[label[0]]->GetXaxis()->SetTitle("impact position (mm)");
  p_amp_x[label[0]]->GetYaxis()->SetTitle("<#detected photons> within 10 ns");
  for(i=1;i<label.size();++i)
    p_amp_x[label[i]]->Draw("same");
  c_amp_max->BuildLegend(0.15,0.75,0.44,0.87);
  if(PrintPlot)
    c_amp_max->Print("p_amp_x_all.png");

  //draw Ncollected photons VS impact point in 2D
  map<TString,TProfile2D*> p2_amp_x; 
  for(int i=0;i<label.size();++i)
  {
    ch[label[i]] -> Draw(Form("@Ph_time.size(Ph_time<=10):mu_x_hit:mu_y_hit>>p2_%s(22,-5.5,5.5)",label[i].Data()),"abs(mu_y_hit)>1.05 && abs(mu_y_hit)<4.45","proffCOLZ");
    p2_amp_x[label[i]] = (TProfile2D*)gDirectory->Get(Form("p2_%s",label[i].Data()));
    if(PrintPlot)
      c_amp_max->Print(Form("p2_amp_x_%s.png",label[i].Data()));
  }




  /*
  ch_STD->Draw("@Ph_time.size():mu_x_hit>>p_STD(22,-5.5,5.5)","abs(mu_y_hit)<2.","proff");
  TProfile* p_STD = (TProfile*)gDirectory->Get("p_STD");
  p_STD->SetMarkerStyle(20);
  p_STD->SetMarkerColor(kRed);
  p_STD->SetTitle("STD SiPM 4x4");

  ch_LASC->Draw("@Ph_time.size():mu_x_hit>>p_LASC(22,-5.5,5.5)","abs(mu_y_hit)<2.","proff");
  TProfile* p_LASC = (TProfile*)gDirectory->Get("p_LASC");
  p_LASC->SetMarkerStyle(20);
  p_LASC->SetMarkerColor(kBlue);
  p_LASC->SetTitle("LASC SiPM 10.5x10.5");

  ch_LASC_nowrap->Draw("@Ph_time.size():mu_x_hit>>p_LASC_nowrap(22,-5.5,5.5)","abs(mu_y_hit)<2.","proff");
  TProfile* p_LASC_nowrap = (TProfile*)gDirectory->Get("p_LASC_nowrap");
  p_LASC_nowrap->SetMarkerStyle(20);
  p_LASC_nowrap->SetMarkerColor(kGreen);
  p_LASC_nowrap->SetTitle("LASC SiPM 10.5x10.5 w\o wrap");

  ch_LASC_norefl->Draw("@Ph_time.size():mu_x_hit>>p_LASC_norefl(22,-5.5,5.5)","abs(mu_y_hit)<2.","proff");
  TProfile* p_LASC_norefl = (TProfile*)gDirectory->Get("p_LASC_norefl");
  p_LASC_norefl->SetMarkerStyle(20);
  p_LASC_norefl->SetMarkerColor(kMagenta);
  p_LASC_norefl->SetTitle("LASC SiPM 10.5x10.5 w\o wrap");
  

  p_STD->Draw();
  p_STD->GetYaxis()->SetRangeUser(1200,2300);
  p_STD->GetXaxis()->SetTitle("impact position (mm)");
  p_STD->GetYaxis()->SetTitle("<#detected photons> within 10 ns");
  p_LASC->Draw("SAME");
  p_LASC_nowrap->Draw("SAME");
  p_LASC_norefl->Draw("SAME");
  c_amp_max->BuildLegend(0.15,0.75,0.44,0.87);
  c_amp_max->Print("amp_max_LASC_STD_comparison.png");
  c_amp_max->Print("amp_max_LASC_STD_comparison.pdf");
  */
  /*
  Long64_t n_STD_centr, n_STD_edge, n_LASC_centr, n_LASC_edge;
  TCanvas *c_wf = new TCanvas();
  n_STD_centr = ch_STD->Draw("eventNb","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  n_STD_edge = ch_STD->Draw("eventNb","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");
  n_LASC_centr = ch_LASC->Draw("eventNb","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  n_LASC_edge = ch_LASC->Draw("eventNb","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");
  n_LASC_nowrap_centr = ch_LASC_nowrap->Draw("eventNb","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  n_LASC_nowrap_edge = ch_LASC_nowrap->Draw("eventNb","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");
  n_LASC_norefl_centr = ch_LASC_norefl->Draw("eventNb","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  n_LASC_norefl_edge = ch_LASC_norefl->Draw("eventNb","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");


  
  ch_STD->Draw("Ph_time>>h_time_STD_centr(1000,0,10)","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  TH1F* h_time_STD_centr = (TH1F*)gDirectory->Get("h_time_STD_centr");
  h_time_STD_centr->Scale(1./n_STD_centr);
  h_time_STD_centr->SetTitle("sipm 4x4, central events");
  h_time_STD_centr->SetLineColor(1);
  h_time_STD_centr->SetLineWidth(2);

  ch_STD->Draw("Ph_time>>h_time_STD_edge(1000,0,10)","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");
  TH1F* h_time_STD_edge = (TH1F*)gDirectory->Get("h_time_STD_edge");
  h_time_STD_edge->Scale(1./n_STD_edge);
  h_time_STD_edge->SetTitle("sipm 4x4, edge events");
  h_time_STD_edge->SetLineColor(2);
  h_time_STD_edge->SetLineWidth(2);

  ch_LASC->Draw("Ph_time>>h_time_LASC_centr(1000,0,10)","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  TH1F* h_time_LASC_centr = (TH1F*)gDirectory->Get("h_time_LASC_centr");
  h_time_LASC_centr->Scale(1./n_LASC_centr);
  h_time_LASC_centr->SetTitle("LASC sipm, central events");
  h_time_LASC_centr->SetLineColor(3);
  h_time_LASC_centr->SetLineWidth(2);

  ch_LASC->Draw("Ph_time>>h_time_LASC_edge(1000,0,10)","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");
  TH1F* h_time_LASC_edge = (TH1F*)gDirectory->Get("h_time_LASC_edge");
  h_time_LASC_edge->Scale(1./n_LASC_edge);
  h_time_LASC_edge->SetTitle("LASC sipm, edge events");
  h_time_LASC_edge->SetLineColor(4);
  h_time_LASC_edge->SetLineWidth(2);


  ch_LASC_nowrap->Draw("Ph_time>>h_time_LASC_nowrap_centr(1000,0,10)","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  TH1F* h_time_LASC_nowrap_centr = (TH1F*)gDirectory->Get("h_time_LASC_nowrap_centr");
  h_time_LASC_nowrap_centr->Scale(1./n_LASC_nowrap_centr);
  h_time_LASC_nowrap_centr->SetTitle("LASC_nowrap sipm, central events");
  h_time_LASC_nowrap_centr->SetLineColor(5);
  h_time_LASC_nowrap_centr->SetLineWidth(2);

  ch_LASC_nowrap->Draw("Ph_time>>h_time_LASC_nowrap_edge(1000,0,10)","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");
  TH1F* h_time_LASC_nowrap_edge = (TH1F*)gDirectory->Get("h_time_LASC_nowrap_edge");
  h_time_LASC_nowrap_edge->Scale(1./n_LASC_nowrap_edge);
  h_time_LASC_nowrap_edge->SetTitle("LASC_nowrap sipm, edge events");
  h_time_LASC_nowrap_edge->SetLineColor(6);
  h_time_LASC_nowrap_edge->SetLineWidth(2);

  ch_LASC_norefl->Draw("Ph_time>>h_time_LASC_norefl_centr(1000,0,10)","abs(mu_x_hit)<2. && abs(mu_y_hit)<2.");
  TH1F* h_time_LASC_norefl_centr = (TH1F*)gDirectory->Get("h_time_LASC_norefl_centr");
  h_time_LASC_norefl_centr->Scale(1./n_LASC_norefl_centr);
  h_time_LASC_norefl_centr->SetTitle("LASC_norefl sipm, central events");
  h_time_LASC_norefl_centr->SetLineColor(7);
  h_time_LASC_norefl_centr->SetLineWidth(2);

  ch_LASC_norefl->Draw("Ph_time>>h_time_LASC_norefl_edge(1000,0,10)","!(abs(mu_x_hit)<2. && abs(mu_y_hit)<2.)");
  TH1F* h_time_LASC_norefl_edge = (TH1F*)gDirectory->Get("h_time_LASC_norefl_edge");
  h_time_LASC_norefl_edge->Scale(1./n_LASC_norefl_edge);
  h_time_LASC_norefl_edge->SetTitle("LASC_norefl sipm, edge events");
  h_time_LASC_norefl_edge->SetLineColor(8);
  h_time_LASC_norefl_edge->SetLineWidth(2);


  h_time_STD_centr->Draw();
  h_time_STD_centr->GetXaxis()->SetRangeUser(0.1,10);
  h_time_STD_centr->GetXaxis()->SetTitle("time (ns)");
  h_time_STD_centr->GetYaxis()->SetTitle("<#detected photons>");
  h_time_STD_edge->Draw("SAME");
  h_time_LASC_centr->Draw("SAME");
  h_time_LASC_edge->Draw("SAME");
  h_time_LASC_nowrap_centr->Draw("SAME");
  h_time_LASC_nowrap_edge->Draw("SAME");
  h_time_LASC_norefl_centr->Draw("SAME");
  h_time_LASC_norefl_edge->Draw("SAME");
  c_wf->BuildLegend(0.15,0.67,0.53,0.88);
  c_wf->Print("arrival_time_LASC_vs_STD.png");
  c_wf->Print("arrival_time_LASC_vs_STD.pdf");

  h_time_STD_centr->GetCumulative()->Draw();
  h_time_STD_edge->GetCumulative()->Draw("SAME");
  h_time_LASC_centr->GetCumulative()->Draw("SAME");
  h_time_LASC_edge->GetCumulative()->Draw("SAME");
  h_time_LASC_nowrap_centr->GetCumulative()->Draw("SAME");
  h_time_LASC_nowrap_edge->GetCumulative()->Draw("SAME");
  h_time_LASC_norefl_centr->GetCumulative()->Draw("SAME");
  h_time_LASC_norefl_edge->GetCumulative()->Draw("SAME");
  c_wf->BuildLegend(0.15,0.67,0.53,0.88);
  c_wf->Print("cumul_arrival_time_LASC_vs_STD_10ns.png");
  c_wf->Print("cumul_arrival_time_LASC_vs_STD_10ns.pdf");

  TCanvas *c_wf_zoom = new TCanvas();
  h_time_STD_centr->GetXaxis()->SetRangeUser(0.1,1.);
  h_time_STD_centr->GetCumulative()->Draw();
  h_time_STD_edge->GetCumulative()->Draw("SAME");
  h_time_LASC_centr->GetCumulative()->Draw("SAME");
  h_time_LASC_edge->GetCumulative()->Draw("SAME");
  h_time_LASC_nowrap_centr->GetCumulative()->Draw("SAME");
  h_time_LASC_nowrap_edge->GetCumulative()->Draw("SAME");
  h_time_LASC_norefl_centr->GetCumulative()->Draw("SAME");
  h_time_LASC_norefl_edge->GetCumulative()->Draw("SAME");
  c_wf_zoom->BuildLegend(0.15,0.67,0.53,0.88);
  c_wf_zoom->Print("cumul_arrival_time_LASC_vs_STD_1ns.png");
  c_wf_zoom->Print("cumul_arrival_time_LASC_vs_STD_1ns.pdf");
  */  
}

void digiPlotter(bool PrintPlot = true)
{
  const double time_offset = 10;
  
  cout<<"In function digiPlotter.C"<<endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  vector<TString> label = {
    "1_LASC_10p5x10p5_mm2",
    "4_LASC_4x4_mm2",
    "4_LASC_4x4_mm2_spec_refl_in_gaps"};

  //1 filename per chain!!!!
  vector<TString> filename = {
    "/eos/user/f/fmonti/SiPM_Optical_Simulation/LASC_10p5x10p5_TILE_10p5x10p5/digi/digi_DCR0.000000MHz_surface5_refl97_Tile11.5x11.5x4x_planar_geometry_0.24844x0.24844x0.8_tilt0_PDE_sourceuniform_seed*.root",
    "/eos/user/f/fmonti/SiPM_Optical_Simulation/4LASC_4x4_TILE_11p5x11p5/digi/digi_DCR0.000000MHz_surface5_refl97_Tile11.5x11.5x4x_planar_geometry_0.20567x0.20567x0.8_tilt0_PDE_sourceuniform_seed*.root",
    "/eos/user/f/fmonti/SiPM_Optical_Simulation/4LASC_4x4_TILE_11p5x11p5_allbackaluminum/digi/digi_DCR0.000000MHz_surface5_refl97_Tile11.5x11.5x4x_planar_geometry_0.20567x0.20567x0.8_tilt0_PDE_sourceuniform_seed*.root"};
  assert(label.size()==filename.size());

  map<TString,TChain*> ch;
  for(int i=0;i<label.size();++i)
  {
    ch[label[i]] = new TChain("digi","digi");
    ch[label[i]] -> Add(filename[i].Data());
  }

  //draw LDE20 for central and for border
  TCanvas *c_time = new TCanvas();
  map<TString,TH1F*> h_time; 
  for(int i=0;i<label.size();++i)
  {
    ch[label[i]] -> Draw(Form("LDE20-10>>h_all_%s(60,0.455,1.055)",label[i].Data()),"","");
    ch[label[i]] -> Draw(Form("LDE20-10>>h_sipm_%s(60,0.455,1.055)",label[i].Data()),"fabs(mu_x_hit)>1 && fabs(mu_x_hit)<3.75 && fabs(mu_y_hit)>1 && fabs(mu_y_hit)<3.75","same");
    ch[label[i]] -> Draw(Form("LDE20-10>>h_gap_%s(60,0.455,1.055)",label[i].Data()),"fabs(mu_x_hit)<0.7 || fabs(mu_y_hit)<0.7","same");
    h_time[label[i]+"all"] = (TH1F*)gDirectory->Get(Form("h_all_%s",label[i].Data()));
    h_time[label[i]+"sipm"] = (TH1F*)gDirectory->Get(Form("h_sipm_%s",label[i].Data()));
    h_time[label[i]+"gap"] = (TH1F*)gDirectory->Get(Form("h_gap_%s",label[i].Data()));
    h_time[label[i]+"all"] -> SetLineColor(4);
    h_time[label[i]+"sipm"] -> SetLineColor(7);
    h_time[label[i]+"gap"] -> SetLineColor(8);
    h_time[label[i]+"all"] -> SetTitle("all");
    h_time[label[i]+"sipm"] -> SetTitle("underneath sipm");
    h_time[label[i]+"gap"] -> SetTitle("central gap");
    h_time[label[i]+"all"] -> Draw();
    h_time[label[i]+"sipm"] -> Draw("same");
    h_time[label[i]+"gap"] -> Draw("same");
    h_time[label[i]+"all"]->GetXaxis()->SetTitle("time(thr 20ph) (ns)");
    c_time->BuildLegend(0.15,0.75,0.44,0.87);
    if(PrintPlot)
      c_time->Print(Form("h_time_%s.png",label[i].Data()));
  }

  //draw LDE20 VS impact point
  TCanvas *c_time_pos = new TCanvas();
  map<TString,TProfile*> p_time_pos; 
  for(int i=0;i<label.size();++i)
  {
    if(label[i]!="1_LASC_10p5x10p5_mm2")
      ch[label[i]] -> Draw(Form("LDE20-10:mu_x_hit>>p_%s(15,-5.5,5.5)",label[i].Data()),"fabs(mu_y_hit)>1 && fabs(mu_y_hit)<3.75","proff");
    else
      ch[label[i]] -> Draw(Form("LDE20-10:mu_x_hit>>p_%s(15,-5.5,5.5)",label[i].Data()),"fabs(mu_y_hit)<4.","proff");
    p_time_pos[label[i]] = (TProfile*)gDirectory->Get(Form("p_%s",label[i].Data()));
    p_time_pos[label[i]]->SetMarkerStyle(20);
    p_time_pos[label[i]]->SetMarkerColor(i+1);
    p_time_pos[label[i]]->SetTitle(label[i]);
    if(PrintPlot)
      c_time_pos->Print(Form("p_time_pos_%s.png",label[i].Data()));
  }
  p_time_pos[label[0]]->Draw();
  p_time_pos[label[0]]->GetYaxis()->SetRangeUser(0.70,0.80);
  p_time_pos[label[0]]->GetXaxis()->SetTitle("impact position (mm)");
  p_time_pos[label[0]]->GetYaxis()->SetTitle("time(thr 20ph) (ns)");
  for(i=1;i<label.size();++i)
    p_time_pos[label[i]]->Draw("same");
  c_time_pos->BuildLegend(0.15,0.75,0.44,0.87);
  if(PrintPlot)
    c_time_pos->Print("p_time_pos_all.png");

  //draw Ncollected photons VS impact point in 2D
  map<TString,TProfile2D*> p2_time_pos; 
  for(int i=0;i<label.size();++i)
  {
    ch[label[i]] -> Draw(Form("LDE20-10:mu_x_hit:mu_y_hit>>p2_%s(15,-5.5,5.5,15,-5.5,5.5)",label[i].Data()),"","proffCOLZ");
    p2_time_pos[label[i]] = (TProfile2D*)gDirectory->Get(Form("p2_%s",label[i].Data()));
    p2_time_pos[label[i]] -> GetXaxis() -> SetTitle("x impact point (mm)");
    p2_time_pos[label[i]] -> GetYaxis() -> SetTitle("y impact point (mm)");
    //p2_time_pos[label[i]] -> GetZaxis() -> SetRangeUser(0.7,0.8);
    if(PrintPlot)
      c_time_pos->Print(Form("p2_time_pos_%s.png",label[i].Data()));
  }

}

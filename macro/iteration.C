

void iteration()
{
  /*
  //creates x position of N cells uniformly distributed on a length l
  const int Ncellx=20;
  const float l=10.5;
  const int Nrepit=20;

  
  cout<<endl;
  cout<<"x = |";
  for(int irepit=0;irepit<Nrepit;++irepit)
    for(int i=0;i<Ncellx;++i)
      cout<< (l/Ncellx)*(i+0.5) - l*0.5 <<"|";
  cout<<endl;
  
  cout<<endl;
  cout<<"y = |";
  for(int i=0;i<Ncellx;++i)
    for(int irepit=0;irepit<Nrepit;++irepit)
      cout<< (l/Ncellx)*(i+0.5) - l*0.5 <<"|";
  cout<<endl;
  */

  //distribute ncells on 4 different sectors 
  const int Ncellx=10;
  const int Ncelly=10;
  const double xcenter=2.75;//offset of single SiPM center 
  const double ycenter=2.75;//offset of single SiPM center 
  const float l=4.;
  const int Nrepit=10;

  cout<<"|";
  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncellx)*(ix+0.5) - l*0.5 - xcenter<<"|";

  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncellx)*(ix+0.5) - l*0.5 + xcenter<<"|";

  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncellx)*(ix+0.5) - l*0.5 - xcenter<<"|";

  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncellx)*(ix+0.5) - l*0.5 + xcenter<<"|";

  cout<<"\n\n\n";
  cout<<"|";
  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncelly)*(iy+0.5) - l*0.5 + xcenter<<"|";

  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncelly)*(iy+0.5) - l*0.5 + xcenter<<"|";

  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncelly)*(iy+0.5) - l*0.5 - xcenter<<"|";

  for(int iy=0;iy<Ncelly;++iy)
    for(int ix=0;ix<Ncellx;++ix)
      cout<< (l/Ncelly)*(iy+0.5) - l*0.5 - xcenter<<"|";

  cout<<"\n\n\n";  
}

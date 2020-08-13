inline double PDK2(double a, double b, double c)
{
  return (a-b-c)*(a+b+c)*(a-b+c)*(a+b-c)/(4*a*a);
}

inline double PDK(double a, double b, double c){
  return TMath::Sqrt( PDK2(a,b,c) );
}
void AssignRandom(double W, const vector<double>& masses, vector<double>& invMass){
  
  double TCM= std::accumulate(masses.begin(),masses.end(), W,  std::minus<double>());

  auto const Nt= invMass.size();
  
  vector<Double_t> rno(Nt);
  rno[0] = 0;
 
  for (int n=1; n<Nt-1; n++)  rno[n]=gRandom->Rndm();   // fNt-2 random numbers
  std::sort(rno.begin(),rno.end()-1);  // sort them
  //  std::sort(rno.begin()+1,rno.end(),std::greater <double>());  // sort them
  rno[Nt-1] = 1;
  double sum=0;
  for (int n=0; n<Nt; n++) {
    sum      += masses[n];
    invMass[n] = rno[n]*TCM + sum;
  }

}
void AssignRandomEquDist(double W, const vector<double>& masses, vector<double>& invMass){
  
  double TCM= std::accumulate(masses.begin(),masses.end(), W,  std::minus<double>());

  auto const Nt= invMass.size();
  
  vector<Double_t> rno(Nt);
  rno[0] = 0;
 
  for (int n=1; n<Nt-1; n++)  rno[n]=(n)*(1./(Nt-1));   // fNt-2 random numbers
  // std::sort(rno.begin(),rno.end()-1);  // sort them
  //  std::sort(rno.begin()+1,rno.end(),std::greater <double>());  // sort them
  rno[Nt-1] = 1;
  double sum=0;
  for (int n=0; n<Nt; n++) {
    sum      += masses[n];
    invMass[n] = rno[n]*TCM + sum;
  }
  cout<<" AssignRandomEquDist "<<endl;
  for(auto& ma:invMass)
    cout<<ma<<" ";
  cout<<endl;
}
/*
void AssignRandomEquVel(double W, const vector<double>& masses, vector<double>& invMass){
  //in rest frame sum of W = (m1^2+(m1*v)^2)^0.5 +  (m2^2+(m2*v)^2)^0.5 +...
  // W = m1(1+(v)^2)^0.5 +  m2(1+(v)^2)^0.5 +...
  // W = (m1 + m2 +..)(1+(v)^2)^0.5
  // = > v^2 = (W^2/(Msum^2) - 1)  and v = sqrt(W^2/(Msum^2) - 1)
  double Msum = std::accumulate(masses.begin(),masses.end(), 0,  std::plus<double>());
  double v = sqrt(W*W/(Msum*Msum) - 1 );

  //calculate momentum of particles
  auto const Nt= invMass.size();
  
  vector<Double_t> rno(Nt);
  rno[0] = 0;
 
  for (int n=1; n<Nt-1; n++)  rno[n]=(n)*(1./(Nt-1));   // fNt-2 random numbers
  std::sort(rno.begin(),rno.end()-1);  // sort them
  //  std::sort(rno.begin()+1,rno.end(),std::greater <double>());  // sort them
  rno[Nt-1] = 1;
  double sum=0;
  for (int n=0; n<Nt; n++) {
    sum      += masses[n];
    invMass[n] = rno[n]*TCM + sum;
  }

}
*/
double CalcWeight(const vector<double>& masses,const vector<double>& invMass){
   
  //
  //-----> compute the weight of the current event
  //
  double wt=1;
  auto const Nt= invMass.size();
  for (int n=0; n<Nt-1; n++) {
    wt*= PDK(invMass[n+1],invMass[n],masses[n+1]);
  }
  
  return wt;
}


//////////////////////////////////////////////////////////////////////////////////////
void MaxOfPhaseSpace(int Nparticles,double W, long nsamples=1E3){

  TLorentzVector Mother(0,0,0,W);
  Double_t masses[Nparticles];

  masses[0]=5;
  for(int i=1;i<Nparticles;i++){
    masses[i]=0.5;
  }

  vector<double> vmasses(masses,masses+Nparticles);

  
  TGenPhaseSpace event;
  event.SetDecay(Mother, Nparticles, masses);

  double maxWeight=0;
  for(int ievent=0;ievent<nsamples;ievent++){
    Double_t weight = event.Generate()/event.GetWtMax();
    if(weight==0) exit(0);
    //cout<<weight<<" "<<maxWeight<<endl;
    if(weight>maxWeight) maxWeight=weight;
  }
  std::cout<<"Max sampled "<<maxWeight<<"  max TGEn  "<<1./event.GetWtMax()<<endl;

  double maxWeight2=0;
  vector<double> invariantM(Nparticles);
  for(int ievent=0;ievent<nsamples;ievent++){
    AssignRandom(W,vmasses,invariantM);
    Double_t weight = CalcWeight(vmasses,invariantM);
    //cout<<weight<<" "<<maxWeight2<<endl;
    if(weight>maxWeight2){
      maxWeight2=weight;
      cout<<ievent<<" "<<weight<<" masses ";
      for(int im=1;im<invariantM.size();im++)
	cout<<(invariantM[im])<<" ";
      cout<<endl;
    }
   }
  
  std::cout<<"Max sampled "<<maxWeight2<<"  max TGEn  "<<1./event.GetWtMax()<<endl;
  AssignRandomEquDist(W,vmasses,invariantM);
  std::cout<<"Max EquDist "<< CalcWeight(vmasses,invariantM) <<"  max TGEn  "<<1./event.GetWtMax()<<endl;
  
}



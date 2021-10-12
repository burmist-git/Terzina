//my
#include "PrimaryGeneratorAction.hh"
#include "VolumeStructures.hh"

//G4
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

//root
#include "TMath.h"

PrimaryGeneratorAction::PrimaryGeneratorAction() :
  _particleGun(0),
  _particleName("pi+"),
  _particleMomentum(3.0*GeV),
  _PhiAngle(0.0*deg),
  _ThetaAngle(0.0*deg),
  _singlePhoton(false)
{
  _particleGun = new G4ParticleGun(1);  
  _BunchXID = 0;

  //backGen = new backgroundGen();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete _particleGun;
  //delete backGen;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle;
  // Correct for center of bar
  G4double xInit, yInit, zInit;
  G4double dX, dY, dZ;
  G4double Ekin, m;
  //G4int pdgID;
  //G4int i;
  //i = 0;

  //_particleMomentum = 7000.0*GeV;
  //_particleMomentum = 20.0*GeV;

  xInit = 15.0*cm;
  yInit = 15.0*cm;
  zInit = -19.0*cm;
  
  if(_particleName != "opticalphoton"){
    ///////////////////////
    _BunchXID++;
    particle = particleTable->FindParticle(_particleName);
    m = particle->GetPDGMass();
    Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum + m*m) - m);
    
    //G4cout<<"_particleMomentum = "<<_particleMomentum<<G4endl
    //	<<"_ThetaAngle       = "<<_ThetaAngle*180.0/TMath::Pi()<<G4endl
    //	<<"_PhiAngle         = "<<_PhiAngle*180.0/TMath::Pi()<<G4endl;
    
    //_ThetaAngle = _ThetaAngle + (-1 + 2*G4UniformRand())*2*TMath::Pi()/360;//mearing of one degree
    //_PhiAngle   = _PhiAngle   + (-1 + 2*G4UniformRand())*2*TMath::Pi()/360;//mearing of one degree
    
    dX =  TMath::Sin(_ThetaAngle)*TMath::Cos(_PhiAngle);
    dY =  TMath::Sin(_ThetaAngle)*TMath::Sin(_PhiAngle);
    dZ =  TMath::Cos(_ThetaAngle);
    
    G4ThreeVector dir(dX, dY, dZ);
    _particleGun->SetParticleDefinition(particle);
    _particleGun->SetParticleMomentumDirection(dir);
    _particleGun->SetParticleEnergy(Ekin);  
    _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
    //for(i = 0;i<10;i++)
    _particleGun->GeneratePrimaryVertex(anEvent);
  }
  else {
    _BunchXID++;

      //World
    G4double world_sizeX = 60.0*cm;
    G4double world_sizeY = 60.0*cm;
    G4double world_sizeZ = 40.0*cm;
    G4double d_l = 1.0*mm;

    G4double eMin_ev = 1.0*eV;
    G4double eMax_ev = 5.0*eV;

    Int_t nn_part = 1;

    for(int i = 0;i<nn_part;i++){
    
      //xInit = (-1 + 2*G4UniformRand())*(world_sizeY/2.0 - d_l);
      xInit = (-1 + 2*G4UniformRand())*200.0;
      //yInit = (-1 + 2*G4UniformRand())*(world_sizeY/2.0 - d_l);
      //yInit = sqrt(30*30-xInit*xInit);
      //xInit = 150;
      //yInit = 150;
      yInit = 0.0;
      zInit = -world_sizeZ/2.0 + d_l;

      
      particle = particleTable->FindParticle(_particleName);
      m = particle->GetPDGMass();
      Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum + m*m) - m);

      //_ThetaAngle = (-1 + 2*G4UniformRand())*TMath::Pi();//mearing of one degree
      //_ThetaAngle = (-1 + 2*G4UniformRand())*TMath::Pi()/180.0*10.0;//mearing of one degree
      _ThetaAngle = (-1 + 2*G4UniformRand())*TMath::Pi()/180.0*0.0;//mearing of one degree
      _PhiAngle   = 0;
    
      dX =  TMath::Sin(_ThetaAngle)*TMath::Cos(_PhiAngle);
      dY =  TMath::Sin(_ThetaAngle)*TMath::Sin(_PhiAngle);
      dZ =  TMath::Cos(_ThetaAngle);
      
      //dX = 0.0;
      //dY = 0.0;
      //dZ = 1.0;
      
      G4ThreeVector dir(dX, dY, dZ);
      _particleGun->SetParticleDefinition(particle);
      _particleGun->SetParticleMomentumDirection(dir);
      _particleGun->SetParticleEnergy(G4UniformRand()*(eMax_ev - eMin_ev) + eMin_ev);  
      _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
      
      _particleGun->GeneratePrimaryVertex(anEvent);
    }
  }
  
}

G4int PrimaryGeneratorAction::GenFlatInt(G4int iMin,G4int iMax){
  G4int val;
  val = (G4int)floor((iMax - iMin + 1)*G4UniformRand() + iMin);
  return val;
}

void PrimaryGeneratorAction::generateThetaAndPhi(){
  _PhiAngle = G4UniformRand()*2*TMath::Pi();
  _ThetaAngle = TMath::Pi() - genCos2dist();
}

G4double PrimaryGeneratorAction::genCos2dist(){
  G4double theta = -999.0;//deg 
  G4double x = -999.0;
  G4double y = -999.0;
  while(theta==-999.0){
    x = G4UniformRand()*(70.0*TMath::Pi()/180.0); //rad
    y = G4UniformRand();
    if(TMath::Power(TMath::Cos(x),1.85)>y){
      theta = x;
    }
  }  
  return theta;
}

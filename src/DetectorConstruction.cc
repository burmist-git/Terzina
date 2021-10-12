//My
#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"

//G4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Color.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
//#include <G4GDMLParser.hh>

//root 
#include "TMath.h"

DetectorConstruction::DetectorConstruction()
{
  //magField = new MagneticField();
  worldVisAtt = new G4VisAttributes();
  quartzVisAtt = new G4VisAttributes();
  aerogelVisAtt = new G4VisAttributes();
  sensitiveVisAtt = new G4VisAttributes();
  pmtboxVisAtt = new G4VisAttributes();
  absVisAtt = new G4VisAttributes();
  // Define Materials to be used
  DefineMaterials();
}

DetectorConstruction::~DetectorConstruction()
{
  //delete magField;
  delete worldVisAtt;
  delete quartzVisAtt;
  delete sensitiveVisAtt;
  delete pmtboxVisAtt;
  delete absVisAtt;
  delete stepLimit;
}

void DetectorConstruction::DefineMaterials()
{
  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  // Define elements
  //G4Element* H = 
  //new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01*g/mole);
  G4Element* C = 
    new G4Element("Carbon",   symbol = "C", z = 6., a = 12.01*g/mole);
  G4Element* N = 
    new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01*g/mole);
  G4Element* O =
    new G4Element("Oxygen",   symbol = "O", z = 8., a = 16.00*g/mole);
  G4Element* Si = 
    new G4Element("Silicon",  symbol = "Si", z = 14., a = 28.09*g/mole);
  G4Element* Al = 
    new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98*g/mole);
  G4Element* F = 
    new G4Element("Fluorine", symbol = "F", z = 9., a = 19.0*g/mole);

  // Quartz Material (SiO2_cladd)
  SiO2_cladd = new G4Material("quartzCladd", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_cladd->AddElement(Si, natoms = 1);
  SiO2_cladd->AddElement(O , natoms = 2);

  // Quartz Material (SiO2_coat)
  SiO2_coat = new G4Material("quartzCoat", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_coat->AddElement(Si, natoms = 1);
  SiO2_coat->AddElement(O , natoms = 2);

  // Quartz Material (SiO2)
  SiO2 = new G4Material("quartz", density = 2.200*g/cm3, ncomponents = 2);
  SiO2->AddElement(Si, natoms = 1);
  SiO2->AddElement(O , natoms = 2);

  // Aerogel Material (SiO2)
  Aerogel = new G4Material("Aerogel", density = 2000*g/m3, ncomponents = 2);
  Aerogel->AddElement(Si, natoms = 1);
  Aerogel->AddElement(O , natoms = 2);

  // C4F10
  C4F10 = new G4Material("fluorocarbon", density = 1.8/1000*g/cm3, ncomponents = 2);
  C4F10->AddElement(C, natoms = 4);
  C4F10->AddElement(F, natoms = 10);

  // Air
  Air = new G4Material("Air", density = 1.290*mg/cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  // Aluminum
  Aluminum = new G4Material("Aluminum", density = 2.7*g/cm3, ncomponents = 1);
  Aluminum->AddElement(Al, fractionmass = 1.0);

  // Aluminum of the mirror
  AluminumMirr = new G4Material("mirrAluminum", density = 2.7*g/cm3, ncomponents = 1);
  AluminumMirr->AddElement(Al, fractionmass = 1.0);

  /*
  // Assign Materials
  world.material = Air;
  secA.material = SiO2;
  secB.material = SiO2;
  secC.material = SiO2;
  secWin.material = SiO2;
  fiberCorr.material = SiO2;
  fiberClad.material = SiO2_cladd;
  //fiberCoat.material = SiO2_coat;
  fiberCoat.material = Aluminum;
  sensitive.material = Aluminum;
  abs1.material = Aluminum;
  */

  //
  // Generate and Add Material Properties Table
  //						
  const G4int num = 36;
  G4double WaveLength[num];
  G4double Absorption[num];      // Default value for absorption
  G4double AirAbsorption[num];
  G4double AirRefractiveIndex[num];
  G4double PhotonEnergy[num];

  // Absorption of quartz per 1m
  G4double QuartzAbsorption[num] =
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,
     0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,
     0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,
     0.990610945};

  G4double C4F10RefractiveIndex[num];
  G4double AerogelRefractiveIndex[num];

  for (int i=0; i<num; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    // Aerogel and C4F10
    Absorption[i] = 100*m;
    // Air
    AirAbsorption[i] = 4.*cm;
    AirRefractiveIndex[i] = 1.;
    // C4F10
    C4F10RefractiveIndex[i] = 1.0014;
    // Aerogel
    AerogelRefractiveIndex[i] = 1.03;
    PhotonEnergy[num - (i+1)] = twopi*hbarc/WaveLength[i];
    // Absorption is given per length and G4 needs mean free path
    // length, calculate it here
    // mean free path length - taken as probablility equal 1/e
    // that the photon will be absorbed
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    //EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*
    //epotekBarJoint.thickness;
  }

  G4double QuartzRefractiveIndex[num] =
    {1.456535,1.456812,1.4571  ,1.457399,1.457712,1.458038,
     1.458378,1.458735,1.459108,1.4595  ,1.459911,1.460344,
     1.460799,1.46128 ,1.461789,1.462326,1.462897,1.463502,
     1.464146,1.464833,1.465566,1.46635 ,1.46719 ,1.468094,
     1.469066,1.470116,1.471252,1.472485,1.473826,1.475289,
     1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

  G4double CladdingRefractiveIndex[num];

  for(int i=0; i<num; i++){
    CladdingRefractiveIndex[i] = TMath::Sqrt(QuartzRefractiveIndex[i]*QuartzRefractiveIndex[i]-0.22*0.22); 
  }

  // Assign absorption and refraction to materials

  // Quartz
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
  QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

  // C4F10
  G4MaterialPropertiesTable* C4F10MPT = new G4MaterialPropertiesTable();
  C4F10MPT->AddProperty("RINDEX", PhotonEnergy, C4F10RefractiveIndex, num);
  C4F10MPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption, num);

  // Aerogel
  G4MaterialPropertiesTable* AerogelMPT = new G4MaterialPropertiesTable();
  AerogelMPT->AddProperty("RINDEX", PhotonEnergy, AerogelRefractiveIndex, num);
  AerogelMPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption, num);
  
  // Cladding (of the fiber) only for the fiber aplication
  G4MaterialPropertiesTable* CladdingMPT = new G4MaterialPropertiesTable();
  CladdingMPT->AddProperty("RINDEX", PhotonEnergy, CladdingRefractiveIndex, num);
  CladdingMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);

  // Assign this material to the bars
  SiO2->SetMaterialPropertiesTable(QuartzMPT);
  SiO2_cladd->SetMaterialPropertiesTable(CladdingMPT);
  C4F10->SetMaterialPropertiesTable(C4F10MPT);
  Aerogel->SetMaterialPropertiesTable(AerogelMPT);
  Air->SetMaterialPropertiesTable(AirMPT);

}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //World
  G4double world_sizeX = 60.0*cm;
  G4double world_sizeY = 60.0*cm;
  G4double world_sizeZ = 40.0*cm;

  //
  G4double global_shift_X0 = 0.0*cm;
  G4double global_shift_Y0 = 0.0*cm;
  G4double global_shift_Z0 = -world_sizeZ/2.0 + 3.0*cm;  

  //Projection sphere
  G4double proj_sphere_thickness = 1.0*mm;
  G4double proj_sphere_R = 20.8*cm;
  G4double proj_sphere_x0 = 0.0*mm;
  G4double proj_sphere_y0 = 0.0*mm;
  G4double proj_sphere_z0 = 315.0*mm - proj_sphere_R + global_shift_Z0;
  proj_sphere_z0=proj_sphere_z0-0;
  G4double proj_sphere_alpha = 16.7*deg;
  G4double proj_sphere_angle = 0.0*deg;

  //Aluminum absorber of the projection sphere
  G4double proj_sphere_al_absorber_thickness = 0.1*mm;
  G4double proj_sphere_al_absorber_R = proj_sphere_R;
  G4double proj_sphere_al_absorber_x0 = 0.0*mm;
  G4double proj_sphere_al_absorber_y0 = 0.0*mm;
  G4double proj_sphere_al_absorber_z0 = proj_sphere_z0 + proj_sphere_thickness;
  G4double proj_sphere_al_absorber_alpha = proj_sphere_alpha;
  
  //Small mirror
  G4double mirror_small_thickness = 1.0*mm;
  G4double mirror_small_R = 30.6*cm;
  //G4double mirror_small_R = 48.0*cm;
  G4double mirror_small_x0 = 0.0*mm;
  G4double mirror_small_y0 = 0.0*mm;
  G4double mirror_small_z0 = 79.0*mm - mirror_small_R + global_shift_Z0;
  mirror_small_z0=mirror_small_z0-4;
  G4double mirror_small_alpha = 17.8*deg;
  //G4double mirror_small_alpha = 10.0*deg;
  G4double mirror_small_angle = 0.0*deg;

  //Aluminum absorber of the small mirror
  G4double mirror_small_al_absorber_thickness = 0.1*mm;
  G4double mirror_small_al_absorber_R = mirror_small_R;
  G4double mirror_small_al_absorber_x0 = 0.0*mm;
  G4double mirror_small_al_absorber_y0 = 0.0*mm;
  G4double mirror_small_al_absorber_z0 = mirror_small_z0-mirror_small_thickness;
  G4double mirror_small_al_absorber_alpha = mirror_small_alpha;
  
  //Large mirror
  G4double mirror_large_thickness = 1.0*mm;
  G4double mirror_large_R = 74.2*cm;
  G4double mirror_large_x0 = 0.0*mm;
  G4double mirror_large_y0 = 0.0*mm;
  G4double mirror_large_z0 = 350.0*mm - mirror_large_R + global_shift_Z0;
  G4double mirror_large_alpha = 23.0*deg;
  G4double mirror_large_angle = 0.0*deg;
  
  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;

  // 
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  //G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,C4F10,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);

  //
  //Projection sphere
  //
  /*
  G4VSolid *projection_sphere_solid = new G4Sphere("Sensitive",                           // Name
						   proj_sphere_R,                         // Inner radius
						   proj_sphere_R + proj_sphere_thickness, // Outer radius
						   0.*deg,                                // Starting phi
						   360.*deg,                              // Delta phi
						   0.*deg,                                // Starting theta
						   proj_sphere_alpha);                    // Delta theta
  */
  G4VSolid *projection_sphere_solid = new G4Sphere("projection_sphere_solid",             // Name
						   proj_sphere_R,                         // Inner radius
						   proj_sphere_R + proj_sphere_thickness, // Outer radius
						   0.*deg,                                // Starting phi
						   360.*deg,                              // Delta phi
						   0.*deg,                                // Starting theta
						   proj_sphere_alpha);                    // Delta theta
  //G4LogicalVolume *sensitive_logical = new G4LogicalVolume(projection_sphere_solid,AluminumMirr,"Sensitive");
  G4LogicalVolume *sensitive_logical = new G4LogicalVolume(projection_sphere_solid,AluminumMirr,"sensitive_logical");
  sensitive_logical->GetName(); 
  Ta.setX(proj_sphere_x0);
  Ta.setY(proj_sphere_y0);
  Ta.setZ(proj_sphere_z0);
  Ra.rotateX(proj_sphere_angle);
  Tr = G4Transform3D(Ra, Ta);
  /*
  G4VPhysicalVolume *projection_sphere_physical = new G4PVPlacement(Tr,                //Transformation
								    sensitive_logical, //its logical volume				 
								    "Sensitive",       //its name
								    world_logical,     //its mother  volume
								    false,             //no boolean operation
								    0);	               //copy number
  */
  G4VPhysicalVolume *projection_sphere_physical = new G4PVPlacement(Tr,                //Transformation
								    sensitive_logical, //its logical volume				 
								    "projection_sphere_physical",       //its name
								    world_logical,     //its mother  volume
								    false,             //no boolean operation
								    0);	               //copy number
  Ra.rotateX(-proj_sphere_angle);

  //
  //Projection sphere absorber
  //
  G4VSolid *projection_sphere_absorber_solid = new G4Sphere("projection_sphere_absorber_solid",    // Name
							    proj_sphere_al_absorber_R,                         // Inner radius
							    proj_sphere_al_absorber_R + proj_sphere_al_absorber_thickness, // Outer radius
							    0.*deg,                                // Starting phi
							    360.*deg,                              // Delta phi
							    0.*deg,                                // Starting theta
							    proj_sphere_al_absorber_alpha);                    // Delta theta
  G4LogicalVolume *projection_sphere_absorber_logical = new G4LogicalVolume(projection_sphere_absorber_solid,Aluminum,"projection_sphere_absorber_logical");
  projection_sphere_absorber_logical->GetName(); 
  Ta.setX(proj_sphere_al_absorber_x0);
  Ta.setY(proj_sphere_al_absorber_y0);
  Ta.setZ(proj_sphere_al_absorber_z0);
  Ra.rotateX(0.0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *projection_sphere_absorber_physical = new G4PVPlacement(Tr,                //Transformation
									     projection_sphere_absorber_logical, //its logical volume				 
									     "projection_sphere_absorber_physical",       //its name
									     world_logical,     //its mother  volume
									     false,             //no boolean operation
									     0);	               //copy number
  Ra.rotateX(-0.0);
  
  //
  //Small mirror
  //
  G4VSolid *small_mirror_solid = new G4Sphere("small_mirror_solid",                    // Name
					      mirror_small_R,                          // Inner radius
					      mirror_small_R + mirror_small_thickness, // Outer radius
					      0.*deg,                                  // Starting phi
					      360.*deg,                                // Delta phi
					      0.*deg,                                  // Starting theta
					      mirror_small_alpha);                     // Delta theta
  G4LogicalVolume *small_mirror_logical = new G4LogicalVolume(small_mirror_solid,AluminumMirr,"small_mirror_logical");
  small_mirror_logical->GetName(); 
  Ta.setX(mirror_small_x0);
  Ta.setY(mirror_small_y0);
  Ta.setZ(mirror_small_z0);
  Ra.rotateX(mirror_small_angle);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *small_mirror_physical = new G4PVPlacement(Tr,                      //Transformation
							       small_mirror_logical,    //its logical volume				 
							       "small_mirror_physical", //its name
							       world_logical,           //its mother  volume
							       false,                   //no boolean operation
							       0);	                //copy number
  Ra.rotateX(-mirror_small_angle);

  //
  //Small mirror absorber
  //
  G4VSolid *small_mirror_absorber_solid = new G4Sphere("small_mirror_absorber_solid",   // Name
						       mirror_small_al_absorber_R,      // Inner radius
						       mirror_small_al_absorber_R + mirror_small_al_absorber_thickness, // Outer radius
						       0.*deg,                          // Starting phi
						       360.*deg,                        // Delta phi
						       0.*deg,                          // Starting theta
						       mirror_small_al_absorber_alpha); // Delta theta
  G4LogicalVolume *small_mirror_absorber_logical = new G4LogicalVolume(small_mirror_absorber_solid,Aluminum,"small_mirror_absorber_logical");
  small_mirror_absorber_logical->GetName(); 
  Ta.setX(mirror_small_al_absorber_x0);
  Ta.setY(mirror_small_al_absorber_y0);
  Ta.setZ(mirror_small_al_absorber_z0);
  Ra.rotateX(0.0);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *small_mirror_absorber_physical = new G4PVPlacement(Tr,                      //Transformation
									small_mirror_absorber_logical,    //its logical volume				 
									"small_mirror_absorber_physical", //its name
									world_logical,           //its mother  volume
									false,                   //no boolean operation
									0);	                //copy number
  Ra.rotateX(-0.0);



  
  //
  //Large mirror
  //
  G4VSolid *large_mirror_solid = new G4Sphere("large_mirror_solid",                    // Name
					      mirror_large_R,                          // Inner radius
					      mirror_large_R + mirror_large_thickness, // Outer radius
					      0.*deg,                                  // Starting phi
					      360.*deg,                                // Delta phi
					      0.*deg,                                  // Starting theta
					      mirror_large_alpha);                     // Delta theta
  G4LogicalVolume *large_mirror_logical = new G4LogicalVolume(large_mirror_solid,AluminumMirr,"large_mirror_logical");
  large_mirror_logical->GetName(); 
  Ta.setX(mirror_large_x0);
  Ta.setY(mirror_large_y0);
  Ta.setZ(mirror_large_z0);
  Ra.rotateX(mirror_large_angle);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *large_mirror_physical = new G4PVPlacement(Tr,                      //Transformation
							       large_mirror_logical,    //its logical volume				 
							       "large_mirror_physical", //its name
							       world_logical,           //its mother  volume
							       false,                   //no boolean operation
							       0);	                //copy number
  Ra.rotateX(-mirror_large_angle);
  //Tr = G4Transform3D(Ra, Ta);


  /*
  //
  // Sensitive volume
  //
  G4VSolid *sensitive_solid = new G4Box("Sensitive", sensitive_sizeX/2.0, sensitive_sizeY/2.0, sensitive_sizeZ/2.0);
  G4LogicalVolume *sensitive_logical = new G4LogicalVolume(sensitive_solid, Aluminum,"Sensitive");
  Ta.setX(sensitive_X0);
  Ta.setY(sensitive_Y0);
  Ta.setZ(sensitive_Z0);
  Ra.rotateX(-sensitive_angle);
  Tr = G4Transform3D(Ra, Ta);
  G4VPhysicalVolume *sensitive_physical_top = new G4PVPlacement(Tr,                //Transformation
								sensitive_logical, //its logical volume				 
								"Sensitive",       //its name
								c4f10_body_logical,     //its mother  volume
								false,	           //no boolean operation
								0);	           //copy number
  Ra.rotateX(sensitive_angle);
  sensitive_physical_top->GetName();
  */
  

  //
  // Set Visualization Attributes
  //
  //G4Color blue        = G4Color(0., 0., 1.);
  //G4Color green       = G4Color(0., 1., 0.);
  G4Color red         = G4Color(1., 0., 0.);
  G4Color white       = G4Color(1., 1., 1.);
  //G4Color cyan        = G4Color(0., 1., 1.);
  G4Color DircColor   = G4Color(0.0, 0.0, 1.0, 0.2);
  G4Color SensColor   = G4Color(0.0, 1.0, 1.0, 0.1);

  /*
  worldVisAtt->SetColor(white);
  worldVisAtt->SetVisibility(true);
  quartzVisAtt->SetColor(DircColor);
  quartzVisAtt->SetVisibility(true);
  
  sensitiveVisAtt->SetColor(red);
  sensitiveVisAtt->SetVisibility(true);
  absVisAtt->SetColor(SensColor);
  absVisAtt->SetVisibility(true);

  sphericalMirror_body_logical->SetVisAttributes(quartzVisAtt);
  flatMirror_logical->SetVisAttributes(quartzVisAtt);

  sensitive_logical->SetVisAttributes(sensitiveVisAtt);
  al_body_logical->SetVisAttributes(absVisAtt);
  c4f10_body_logical->SetVisAttributes(absVisAtt);

  aerogel_body_logical->SetVisAttributes(quartzVisAtt);


  */

  //world.logical->SetVisAttributes(worldVisAtt);
  
  //
  // Define Optical Borders
  //

  // Surface for killing photons at borders
  const G4int num1 = 2;
  G4double Ephoton[num1] = {1.5*eV, 5.8*eV};

  G4OpticalSurface* OpVolumeKillSurface =
    new G4OpticalSurface("VolumeKillSurface");
  OpVolumeKillSurface->SetType(dielectric_metal);
  OpVolumeKillSurface->SetFinish(polished);
  OpVolumeKillSurface->SetModel(glisur);


  G4double ReflectivityKill[num1] = {0., 0.};
  G4double EfficiencyKill[num1] = {1., 1.};
  G4MaterialPropertiesTable* VolumeKill = new G4MaterialPropertiesTable();
  VolumeKill->AddProperty("REFLECTIVITY", Ephoton, ReflectivityKill, num1);
  VolumeKill->AddProperty("EFFICIENCY",   Ephoton, EfficiencyKill,   num1);
  OpVolumeKillSurface->SetMaterialPropertiesTable(VolumeKill);
  new G4LogicalSkinSurface("SensitiveSurface", 
  			   sensitive_logical, OpVolumeKillSurface);
  



  // Define mirror surface
  const G4int num2 = 36;
  G4double EfficiencyMirrors[num2];
  G4double WaveLength[num2];
  G4double PhotonEnergy[num2];
  G4double MirrorReflectivity[num2];
  for (G4int i=0; i<num2; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    PhotonEnergy[num2 - (i+1)] = twopi*hbarc/WaveLength[i];
    EfficiencyMirrors[i] = 0.0;
    MirrorReflectivity[i]=0.85;
  }
  /*
  G4double MirrorReflectivity[num2]=
    {0.87,0.88,0.885,0.89,0.895,0.9,0.905,0.91,0.915,0.92,0.923,0.9245,
     0.926,0.928,0.93,0.935,0.936,0.937,0.938,0.94,0.94,0.939,0.9382,
     0.938,0.937,0.937,0.936,0.935,0.934,0.932,0.93,0.928,0.926,0.924,
     0.922,0.92};
  */
  G4MaterialPropertiesTable* MirrorMPT = new G4MaterialPropertiesTable();
  MirrorMPT->AddProperty("RELECTIVITY", PhotonEnergy, MirrorReflectivity, num2);
  MirrorMPT->AddProperty("EFFICIENCY" , PhotonEnergy, EfficiencyMirrors,  num2);

  G4OpticalSurface* OpMirrorSurface = new G4OpticalSurface("MirrorSurface");
  OpMirrorSurface->SetType(dielectric_metal);
  OpMirrorSurface->SetFinish(polished);
  OpMirrorSurface->SetModel(glisur);

  new G4LogicalSkinSurface("MirrorSurfT",
			   small_mirror_logical, OpMirrorSurface);
  new G4LogicalSkinSurface("MirrorSurfT",
			   large_mirror_logical, OpMirrorSurface);

  OpMirrorSurface->SetMaterialPropertiesTable(MirrorMPT);
  ///////////////////////////////////////


  
  // 
  // Sensitive detector definition
  //

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SensitiveDetector* aSD = new SensitiveDetector("fTOF");
  SDman->AddNewDetector(aSD);
  sensitive_logical->SetSensitiveDetector(aSD);
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4double maxStep   = 0.1*mm;
  G4double maxLength = 2.0*m;
  G4double maxTime   = 20.0*ns; 
  G4double minEkin   = 1.0/100*MeV;
  G4double mionRang  = 0.01*mm;
  stepLimit = new G4UserLimits(maxStep,maxLength,maxTime,minEkin,mionRang);
  /*
  secA.logical->SetUserLimits(stepLimit);
  secB.logical->SetUserLimits(stepLimit);
  secC.logical->SetUserLimits(stepLimit);
  secWin.logical->SetUserLimits(stepLimit);
  */

  //G4GDMLParser parser;
  //parser.Write("CpFM.gdml", world.physical);

  return world_physical;
}

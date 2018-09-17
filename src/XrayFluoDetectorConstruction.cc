//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: XrayFluoDetectorConstruction.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//    Nov 2002 Alfonso Mantero materials added,
//             Material selection implementation
// 16 Jul 2003 Alfonso Mantero Detector type selection added + minor fixes
// -------------------------------------------------------------------

#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoSD.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"


#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "XrayFluoNistMaterials.hh"
#include "G4SDManager.hh"


// #include "G4Region.hh"
// #include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoDetectorConstruction::XrayFluoDetectorConstruction()
  : aNavigator(0), detectorType(0),sampleGranularity(false), phaseSpaceFlag(false),
    DeviceSizeX(0), DeviceSizeY(0),DeviceThickness(0),
    solidWorld(0),logicWorld(0),physiWorld(0),
    solidHPGe(0),logicHPGe(0),physiHPGe(0),
    solidSample (0),logicSample(0),physiSample (0),
    solidDia1(0),logicDia1(0),physiDia1(0),
    solidDia3(0),logicDia3(0),physiDia3(0),
    solidOhmicPos(0),logicOhmicPos(0), physiOhmicPos(0),
    solidWindow(0), logicWindow(0), physiWindow(0),
    solidOhmicNeg(0),logicOhmicNeg(0), physiOhmicNeg(0),
    solidPixel(0),logicPixel(0), physiPixel(0),
    OhmicPosMaterial(0), OhmicNegMaterial(0),
    pixelMaterial(0),sampleMaterial(0),
    Dia1Material(0),Dia3Material(0),
    defaultMaterial(0), windowMaterial (0), Teflon(0), SampleRegion(0)
{
  materials = XrayFluoNistMaterials::GetInstance();

  HPGeSD.Put(0);

  aNavigator = new G4Navigator();

  DefineDefaultMaterials();

  NbOfPixelRows     =  1; // should be 1
  NbOfPixelColumns  =  1; // should be 1
  NbOfPixels        =  NbOfPixelRows*NbOfPixelColumns;
  PixelSizeXY       =  std::sqrt(40.) * mm;
  PixelThickness = 2.7 * mm; //should be 3.5 mm

  G4cout << "PixelThickness(mm): "<< PixelThickness/mm << G4endl;
  G4cout << "PixelSizeXY(cm): "<< PixelSizeXY/cm << G4endl;

  ContactSizeXY     = PixelSizeXY; //std::sqrt(40) * mm; //should be the same as PixelSizeXY
  SampleThickness = 0 * mm;
  SampleSizeXY = 0. * mm; //need to change to 5 mm it's the diameter of the target.

  SampleRadius=0.25*mm;
  Dia1Thickness = 1. *mm;
  Dia3Thickness = 1. *mm;
  Dia1SizeXY = 3. *cm;
  Dia3SizeXY = 3. *cm;

  Ablatorthickness = 175. * um;
  Culayerthickness = 25. * um;
  Stepthickness = 310. * um;
  Allayerthickness = 50. * um;

  InnerRadius=105*um; //MixCuCH
  SecondRadius=125*um;// DenseMixCuCH
  ThirdRadius=145*um;//DenseCH 7.25 g/cc,
  OuterRadius=145*um; //No this layer
//20*(r2^3-r1^3)+r1^3*1.1+(r3^3-r2^3)*10=212.5^3*1.1 Cu mass conservation
//(r3^3-r2^3)*10=(212.5^3-182.5^3)*0.87 CH mass convervation
//when r1=130 um, r2=((212.5^3*1.1-(212.5^3-182.5^3)*0.87)/20+18.9/20*130^3)^(1/3)
//r3=(((212.5^3*1.1-(212.5^3-182.5^3)*0.87)/20+18.9/20*130^3)+(212.5^3-182.5^3)*0.087)^(1/3)
  filtersizeXY=2.5*cm;
  Teflonthickness=10*mm;



  G4double filterthickness[15];
  for(G4int t=0; t<5; t++)
    filterthickness[t]=0.1*mm;

  filterthickness[5]=0.15*mm;
  filterthickness[6]=0.5*mm;
  filterthickness[7]=0.5*mm;
  filterthickness[8]=1.56*mm;
  filterthickness[9]=1*mm;
  filterthickness[10]=2*mm;
  filterthickness[11]=3*mm;
  filterthickness[12]=4*mm;
  filterthickness[13]=6.4*mm;
  filterthickness[14]=6.4*mm;




  DiaInnerSize = 2.9 * cm; //(Hole in the detector's diaphragm) it was 1 mm


  OhmicNegThickness = 1e-6*cm;// 0.005
  OhmicPosThickness = 1e-6*cm;// 0.005
  windowThickness = 3 * mm; //changed it for OMEGAEP window
  ThetaHPGe = 135. * deg;
  PhiHPGe = 225. * deg;

  ThetaDia1 = 135. * deg;
  PhiDia1 = 90. * deg;
  AlphaDia1 = 225. * deg;

  AlphaDia3 = 180. * deg;
  Dia3Dist =  66.5 * mm;
  Dia3InnerSize = 1. * mm;
  ThetaDia3 = 180. * deg;
  PhiDia3 = 90. * deg;

  DistDia = 66.5 * mm;
  DistDe =DistDia+ (Dia1Thickness
		    +PixelThickness)/2+OhmicPosThickness+windowThickness ;

  grainDia = 1 * mm;
  PixelCopyNb=0;
  grainCopyNb=0;
  G4String defaultDetectorType = "sili";

  targetCuts = new G4ProductionCuts();//added for cuts
  SampleRegion=0;
  ComputeApparateParameters();

//   G4String regName = "SampleRegion";
//   sampleRegion = new G4Region(regName);

  if (!phaseSpaceFlag) SetDetectorType(defaultDetectorType);

  // create commands for interactive definition of the apparate

  detectorMessenger = new XrayFluoDetectorMessenger(this);

  G4cout << "XrayFluoDetectorConstruction created" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoDetectorConstruction* XrayFluoDetectorConstruction::instance = 0;

XrayFluoDetectorConstruction* XrayFluoDetectorConstruction::GetInstance()
{
  if (instance == 0)
    {
      instance = new XrayFluoDetectorConstruction;

    }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetDetectorType(G4String type)
{
  if (type=="sili")
    {
      detectorType = XrayFluoSiLiDetectorType::GetInstance();
    }
   else if (type=="hpge")
     {
       detectorType = XrayFluoHPGeDetectorType::GetInstance();
    }/*
   else if (type=="aifira")
     {
       detectorType = XrayFluoAifSiLi::GetInstance();
       }*/
  else
    {
      G4ExceptionDescription execp;
      execp << type + "detector type unknown";
      G4Exception("XrayFluoDataSet::LoadData()","example-xray_fluorescence06",
	  FatalException, execp);
    }
  //GeometryHasBeenModified invoked by the messenger

}

XrayFluoVDetectorType* XrayFluoDetectorConstruction::GetDetectorType() const
{
  return detectorType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorConstruction::~XrayFluoDetectorConstruction()

{
  delete detectorMessenger;
  delete detectorType;
  G4cout << "XrayFluoDetectorConstruction deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XrayFluoDetectorConstruction::Construct()
{
  return ConstructApparate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::DefineDefaultMaterials()
{


  //define materials of the apparate

  sampleMaterial = materials->GetMaterial("Dolorite");
  Dia1Material = materials->GetMaterial("G4_Pb");
  Dia3Material = materials->GetMaterial("G4_Galactic");
  pixelMaterial = materials->GetMaterial("SiLi");
  //OhmicPosMaterial = materials->GetMaterial("G4_Cu");
  OhmicPosMaterial = materials->GetMaterial("G4_Ni");
  OhmicNegMaterial = materials->GetMaterial("G4_Pb");
  defaultMaterial = materials->GetMaterial("G4_Galactic");
  windowMaterial = materials->GetMaterial("G4_Pyrex_Glass");//changed for OMEGAEP
  CH = materials->GetMaterial("G4_POLYSTYRENE");
  Cu = materials->GetMaterial("G4_Cu");
  Al = materials->GetMaterial("G4_Al");
  Pb = materials->GetMaterial("G4_Pb");
  Teflon = materials-> GetMaterial("G4_TEFLON");
  Mylar = materials -> GetMaterial("G4_MYLAR");
  Ti = materials -> GetMaterial("G4_Ti");
  Fe = materials -> GetMaterial("G4_Fe");
  Mo = materials -> GetMaterial("G4_Mo");
  Ag = materials -> GetMaterial("G4_Ag");
  Sn = materials -> GetMaterial("G4_Sn");
  Ta = materials -> GetMaterial("G4_Ta");
  Au = materials -> GetMaterial("G4_Au");

  MixCuCH = materials -> GetMaterial("MixCuCH");
  DenseMixCuCH = materials -> GetMaterial("DenseMixCuCH");
  CuFoam  = materials -> GetMaterial("CuFoam");
  DenseCH = materials -> GetMaterial("DenseCH");
  DenseCu = materials -> GetMaterial("DenseCu");

    G4Material*  filtermateriallist[15];
    filtermateriallist[0]=materials->GetMaterial("G4_Al");
    filtermateriallist[1]=materials -> GetMaterial("G4_Ti");
    filtermateriallist[2]= materials -> GetMaterial("G4_Fe");
    filtermateriallist[3]=Cu;
    filtermateriallist[4]=Mo;
    filtermateriallist[5]=Ag;
    filtermateriallist[6]=Sn;
    filtermateriallist[7]=Ta;
    filtermateriallist[8]=Au;
    for(G4int f=9; f<15; f++)
      filtermateriallist[f]=Pb;









  //I don't understand why need the variable "materials", why not use G4Material directly.
  //but it does not work Cu = G4Material::GetMaterial("G4_Cu");
  //I don't know why.
}

void XrayFluoDetectorConstruction::SetOhmicPosThickness(G4double val)
{

  if (!phaseSpaceFlag) {


    if (val == 0.0) {
      OhmicPosMaterial = materials->GetMaterial("G4_Galactic");
    }
    else {
      OhmicPosThickness = val;
      //OhmicPosMaterial = materials->GetMaterial("G4_Cu");
      OhmicPosMaterial = materials->GetMaterial("G4_Ni");
    }

  }
  else{
    G4cout << "Not available in this configuration" << G4endl;
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XrayFluoDetectorConstruction::ConstructApparate()
{
  // complete the apparate parameters definition

  //ComputeApparateParameters();

  //world and associated navigator

  //solidWorld = new G4Sphere("World",	      		        //its name
  solidWorld = new G4Sphere("World",0.,
   WorldSizeXY/2,0., twopi, 0., pi);;	//its size

  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
  physiWorld = new G4PVPlacement(0,			//no rotation
				 G4ThreeVector(),	//at (0,0,0)
				 "World",		//its name
				 logicWorld,		//its logical volume
				 0,			//its mother  volume
				 false,			//no boolean operation
				 0);			//copy number

  aNavigator->SetWorldVolume(physiWorld);

  G4cout << "World constructed" << G4endl;
  G4bool ConstructDetector=0;
if (ConstructDetector){
  solidBMXS= new G4Box("BMXS",
              10.*cm/2, 10.*cm/2, 20.*cm/2);
  solidBMXSCollimator= new G4Tubs("Collimator", 0., 0.64*cm, //1/2 inch diameter collimator
                                10.*cm/2, 0, 360*degree );//name, r inner, r Outer, length, start degree, end degree

  logicBMXS = new G4LogicalVolume(solidBMXS,
                                    Pb,
                                    "BMXS");
  logicBMXSCollimator = new G4LogicalVolume(solidBMXSCollimator,
                                          defaultMaterial,
                                          "Collimator");

  physiBMXS = new G4PVPlacement(0, G4ThreeVector(0, 0, 15*cm),
                                  "BMXS",
                                logicBMXS,
                               physiWorld,
                              false,
                            0);// when place new BMXS connon, this need to be 1
  physiBMXSCollimator = new G4PVPlacement(0, G4ThreeVector(0,0,-5.*cm),
                                        "Collimator",
                                        logicBMXSCollimator,
                                       physiBMXS,
                                      false,
                                      0);


  physiBMXS55 = new G4PVPlacement(0, G4ThreeVector(15*cm, 0, 15*cm),
                                  "BMXS55",
                                logicBMXS,
                               physiWorld,
                              false,
                            0);// when place new BMXS connon, this need to be 1
    G4cout << "BMXS constructed" << G4endl;

    G4cout << "begin BMXS filters" << G4endl;

    solidCartridge = new G4Box("Cartridge", 3*cm, 3*cm, 5*cm);
    logicCartridge = new G4LogicalVolume(solidCartridge,
                                        defaultMaterial,
                                        "Cartridge");
    physiCartridge = new G4PVPlacement(0, G4ThreeVector(0,0,5*cm),
                                      "Cartridge",
                                       logicCartridge,
                                     physiBMXS,
                                   false,
                                 0);

     solidWindow=0; logicWindow=0; physiWindow=0;

 	  if (windowThickness > 0.)
 	    { solidWindow = new G4Box("Window",		//its name
 					filtersizeXY/2, filtersizeXY/2,windowThickness/2);

 	    logicWindow = new G4LogicalVolume(solidWindow,    //its solid
 						windowMaterial, //its material
 						"Window");      //its name
 	    physiWindow = new G4PVPlacement(0,
 					      G4ThreeVector(0.,
 							    0.,
 							    -5*cm+windowThickness/2),
 					      "Window",
 					      logicWindow,
 					      physiCartridge,
 					      false,
 					      0);

 	    }
    // I would need solidTeflonBlock, solid100umfilter, solidMylar, and some other thickness filter solid
    solidTeflonBlock= new G4Box("TeflonBlock",
                   filtersizeXY/2, filtersizeXY/2, Teflonthickness/2);
    logicTeflonBlock= new G4LogicalVolume(solidTeflonBlock,
                                          Teflon,
                                          "TeflonBlock");
    if(logicTeflonBlock->GetMaterial())
      G4cout<<"Teflon "<<G4endl;

    physiTeflonBlock = new G4PVPlacement(0,G4ThreeVector(0,0, -4.5*cm+windowThickness),
                                        "TeflonBlock",
                                       logicTeflonBlock,
                                     physiCartridge,
                                   false,
                                 0);


}


/*
    //solid100umfilter= new G4Box("100umfilter",
    //            filtersizeXY/2, filtersizeXY/2, (G4doulbe thinfilter=100*um)/2);
    solidMylar = new G4Box("Mylar",
                          filtersizeXY/2,filtersizeXY/2, (G4double mylarthickness=250*um)/2);
    logicMylar = new G4LogicalVolume(solidMylar,
                                     Mylar,
                                    "Mylar");
    solidIP=new G4Box("IP",
                     filtersizeXY/2,filtersizeXY/2, (G4double IPthickness=500*um)/2));
    logicIP=new G4LogicalVolume(solidIP,
                              defaultMaterial,
                            "IP");


    //solidAgfilter= new G4Box("Agfilter",
    //                        filtersizeXY/2, filtersizeXY/2, 0.15*mm);
    //solid500umfilter= new G4Box("500umfilter",
    //                filtersizeXY/2, filtersizeXY/2, 250*um);

    //solidAufilter= new G4Box("Aufilter",
    //                                filtersizeXY/2, filtersizeXY/2, 1.56*mm);
    //solidPbfilter1= new G4Box("Pbfilter1",
    //                                filtersizeXY/2, filtersizeXY/2, 1.*mm);
    //solidPbfilter2= new G4Box("Pbfilter1",
    //                                filtersizeXY/2, filtersizeXY/2, 2.*mm);

    G4double Mylar1Position[15];
    Mylar1Position[0]=TeflonPosition+Teflonthickness/2+mylarthickness/2;
    for (G4int mp=1;mp<15;mp++)
    {
      Mylar1Position[mp]=Mylar1Position[mp-1]+filterthickness[mp-1]+2*mylarthickness+IPthickness;
    }


    G4Box * solidfilter[15];
    G4LogicalVolume * logicfilter[15];
    G4PVPlacement * physifilter[15];

    G4String name="filter";
    for (G4int sf=0; sf<15, sf++)
    {
      solidfilter[sf]= new G4Box(name+G4String::to_string(sf+1),
      filtersizeXY/2, filtersizeXY/2, filterthickness[sf]/2);

      logicfilter[sf]= new G4LogicalVolume(solidfilter[sf],
                                          filtermateriallist[sf],
                                          name+G4String::to_string(sf+1));
      physifilter[sf]= new G4PVPlacement(0, G4ThreeVector(0,0,
                      Mylar1Position[sf]+(filterthickness[sf]+mylarthickness)/2),
                                        name+G4String::to_string(sf+1),
                                      logicfilter[sf],
                                    physiCartridge,
                                  false,
                                0);
      physiMylar1 = new G4PVPlacement(0,G4ThreeVector(0,0,Mylar1Position[sf]),
                                    "Mylar1",
                                  logicMylar,
                                physiCartridge,
                              false,
                            sf);
      physiMylar2 =new G4PVPlacement(0,G4ThreeVector(0,0,
        Mylar1Position[sf]+filterthickness[sf]+IPthickness+mylarthickness),
                                    "Mylar2",
                                  logicMylar,
                                physiCartridge,
                              false,
                            sf);
      physiIP =new  G4PVPlacement(0,G4ThreeVector(0,0,Mylar1Position[sf]+filterthickness[sf]+(IPthickness+mylarthickness)/2),
                                    "IP",
                                  logicIP,
                                physiCartridge,
                              false,
                            sf);
    }
*/



/*
  //HPGeDetector

  if (!phaseSpaceFlag) {

    solidHPGe = 0;  physiHPGe = 0;  logicHPGe=0;
    solidPixel=0; logicPixel=0; physiPixel=0;

    if (DeviceThickness > 0.)
      {
	solidHPGe = new G4Box("HPGeDetector",		//its name
			      DeviceSizeX/2,DeviceSizeY/2,DeviceThickness/2);//size


	logicHPGe = new G4LogicalVolume(solidHPGe,	//its solid
					defaultMaterial,	//its material
					"HPGeDetector");	//its name

	zRotPhiHPGe.rotateX(PhiHPGe);
	G4double x,y,z;
	z = DistDe * std::cos(ThetaHPGe);
	y =DistDe * std::sin(ThetaHPGe);
	x = 0.*cm;
	physiHPGe = new G4PVPlacement(G4Transform3D(zRotPhiHPGe,G4ThreeVector(x,y,z)),
				      "HPGeDetector",	//its name
				      logicHPGe,	//its logical volume
				      physiWorld,	//its mother  volume
				      false,		//no boolean operation
				      0);		//copy number
      }
    // Pixel




    for ( G4int j=0; j < NbOfPixelColumns ; j++ )
      { for ( G4int i=0; i < NbOfPixelRows ; i++ )
	{
	  solidPixel=0; logicPixel=0;   physiPixel=0;
	  if (PixelThickness > 0.)
	    solidPixel = new G4Box("Pixel",
				   PixelSizeXY/2,PixelSizeXY/2, PixelThickness/2);

	  logicPixel = new G4LogicalVolume(solidPixel,
					   pixelMaterial,	//its material
					   "Pixel");	        //its name

	  /*
	    zRotPhiHPGe.rotateX(PhiHPGe);
	    G4double x,y,z;
	    z = DistDe * std::cos(ThetaHPGe);
	    y =DistDe * std::sin(ThetaHPGe);
	    x = 0.*cm;*/


/*

	  physiPixel = new G4PVPlacement(0,
					 G4ThreeVector(0,
						       i*PixelSizeXY,
						       j*PixelSizeXY ),
					 "Pixel",
					 logicPixel,	 //its logical volume
					 physiHPGe, //its mother  volume
					 false,	 //no boolean operation
					 PixelCopyNb);//copy number






	  // OhmicNeg

	  solidOhmicNeg=0; logicOhmicNeg=0; physiOhmicNeg=0;

	  if (OhmicNegThickness > 0.)
	    { solidOhmicNeg = new G4Box("OhmicNeg",		//its name
					PixelSizeXY/2,PixelSizeXY/2,OhmicNegThickness/2);

	    logicOhmicNeg = new G4LogicalVolume(solidOhmicNeg,    //its solid
						OhmicNegMaterial, //its material
						"OhmicNeg");      //its name

	    physiOhmicNeg = new G4PVPlacement(0,
					      G4ThreeVector
					      (0.,
					       0.,
					       (PixelThickness+OhmicNegThickness)/2),
					      "OhmicNeg",        //its name
					      logicOhmicNeg,     //its logical volume
					      physiHPGe,        //its mother
					      false,             //no boulean operat
					      PixelCopyNb);                //copy number

	    }
	  // OhmicPos
	  solidOhmicPos=0; logicOhmicPos=0; physiOhmicPos=0;

	  if (OhmicPosThickness > 0.)
	    { solidOhmicPos = new G4Box("OhmicPos",		//its name
					PixelSizeXY/2,PixelSizeXY/2,OhmicPosThickness/2);

	    logicOhmicPos = new G4LogicalVolume(solidOhmicPos,    //its solid
						OhmicPosMaterial, //its material
						"OhmicPos");      //its name

	    physiOhmicPos = new G4PVPlacement(0,
					      G4ThreeVector(0.,
							    0.,
							    (-PixelThickness-OhmicPosThickness)/2),
					      "OhmicPos",
					      logicOhmicPos,
					      physiHPGe,
					      false,
					      PixelCopyNb);

	    }

	  /////////// widow place here! ////////////////
	  // OhmicPos
	  solidWindow=0; logicWindow=0; physiWindow=0;

	  if (windowThickness > 0.)
	    { solidWindow = new G4Box("Window",		//its name
					PixelSizeXY/2,PixelSizeXY/2,windowThickness/2);

	    logicWindow = new G4LogicalVolume(solidWindow,    //its solid
						windowMaterial, //its material
						"Window");      //its name

	    physiWindow = new G4PVPlacement(0,
					      G4ThreeVector(0.,
							    0.,
							    ((-PixelThickness-windowThickness)/2)
							    -OhmicPosThickness),
					      "OhmicWindow",
					      logicWindow,
					      physiHPGe,
					      false,
					      PixelCopyNb);

	    }



	  PixelCopyNb += PixelCopyNb;
	  G4cout << "PixelCopyNb: " << PixelCopyNb << G4endl;
	}

      }

  }
  */

  //Sample

  if (sampleGranularity) {

    solidSample=0;  logicSample=0;  physiSample=0;
    if (SampleThickness > 0.)
      {
	solidSample = new G4Box("Sample",		//its name
				SampleSizeXY/2,SampleSizeXY/2,SampleThickness);//size

	logicSample= new G4LogicalVolume(solidSample,	//its solid
					 defaultMaterial,	//its material
					 "Sample");	//its name

	physiSample = new G4PVPlacement(0,			//no rotation
					G4ThreeVector(),	//at (0,0,0)
					"Sample",	//its name
					logicSample,	//its logical volume
					physiWorld,	//its mother  volume
					false,		//no boolean operation
					0);		//copy number

      }




    G4int nbOfGrainsX = ((G4int)(SampleSizeXY/grainDia)) -1 ;

    // y dim of a max density plane is 2rn-(n-1)ar, wehere a = (1-(std::sqrt(3)/2)), n is
    // number of rows and r the radius of the grain. so the Y-dim of the sample must
    // be greater or equal to this. It results that nmust be <= (SampleY-a)/(1-a).
    // Max Y shift of the planes superimposing along Z axis is minor (2/std::sqrt(3)r)

    G4double a = (1.-(std::sqrt(3.)/2.));
    G4int nbOfGrainsY =  (G4int) ( ((SampleSizeXY/(grainDia/2.)) -a)/(2.-a) ) -1;

    // same for the z axis, but a = 2 * (std::sqrt(3) - std::sqrt(2))/std::sqrt(3)

    G4double b = 2. * (std::sqrt(3.) - std::sqrt(2.))/std::sqrt(3.);
    G4int nbOfGrainsZ =  (G4int) ( ((SampleThickness/(grainDia/2.)) -b)/(2.-b) )-1;

    if (SampleThickness > 0.){

      solidGrain=0; logicGrain=0; physiGrain=0;
      solidGrain = new G4Sphere("Grain",0.,
				grainDia/2,0., twopi, 0., pi);

      logicGrain = new G4LogicalVolume(solidGrain,
				       sampleMaterial,	//its material
				       "Grain");	        //its name
      G4ThreeVector grainPosition;
      G4double grainInitPositionX = 0;
      G4double grainInitPositionY = 0;
      G4double grainInitPositionZ = (-1.*SampleThickness/2.+grainDia/2.);
      G4double grainStepX = grainDia = 0;
      G4double grainStepY = grainDia*(1.-(0.5-(std::sqrt(3.)/4.)));
      G4double grainStepZ = grainDia*std::sqrt(2./3.);

      for ( G4int k=0; k < nbOfGrainsZ ; k++ ) {
	for ( G4int j=0; j < nbOfGrainsY ; j++ ) {
	  for ( G4int i=0; i < nbOfGrainsX ; i++ ) {

	    // Now we identify the layer and the row where the grain is , to place it in the right position



	    if (k%3 == 0) { // first or (4-multiple)th layer: structure is ABCABC
	      grainInitPositionY = (-1.*SampleSizeXY/2.+grainDia/2.);
	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia/2.);
	      }

	      else if ( ((j+1) % 2)  == 0 ) {
		grainInitPositionX = (-1.*SampleSizeXY/2.+ grainDia);
	      }

	    }
	    else if ( ((k+2) % 3) == 0 ) { // B-layer

	      grainInitPositionY = ( (-1.*SampleSizeXY/2.) + (grainDia/2.)*(1. + (1./std::sqrt(3.)) ) );

	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia);
	      }

	      else if ( (j+1)%2  == 0 ) {
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia/2);
	      }

	    }

	    else if ( (k+1)%3 == 0 ) { // B-layer

	      grainInitPositionY = (-1.*SampleSizeXY/2.+(grainDia/2.)*(1.+2./std::sqrt(3.)) );

	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia/2.);
	      }

	      else if ( (j+1)%2  == 0 ) {
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia);
	      }

	    }

	    physiGrain = new G4PVPlacement(0,
					   G4ThreeVector( grainInitPositionX + i*grainStepX,
							  grainInitPositionY + j*grainStepY,
							  grainInitPositionZ + k*grainStepZ),
					   "Grain",
					   logicGrain,	 //its logical volume
					   physiSample, //its mother  volume
					   false,	 //no boolean operation
					   grainCopyNb);//copy number

	    grainCopyNb = grainCopyNb +1;
	  }
	}
      }
    }
    G4cout << "I hope nothing here" << G4endl;

  }
  else {
    G4cout << "begin target constrction" << G4endl;


    solidSample=0;  logicSample=0;  physiSample=0;
    if (SampleThickness > 0.)
      {
	solidSample = new G4Box("Sample",		//its name
				SampleSizeXY/2,SampleSizeXY/2,SampleThickness/2);//size

	logicSample= new G4LogicalVolume(solidSample,	//its solid
					 defaultMaterial,	//its material
					 "Sample");	//its name

	physiSample = new G4PVPlacement(0,			//no rotation
					G4ThreeVector(0,0,0),	//at (0,0,0)
					"Sample",	//its name
					logicSample,	//its logical volume
					physiWorld,	//its mother  volume
					false,		//no boolean operation
					0);		//copy number

    solidAblator= new G4Tubs("Ablator",
                   0, SampleSizeXY/2, Ablatorthickness/2, 0, 360.*degree ); // name Ablator,
                   //inner radus 0, Outer radius half of SampleSizeXY,
                   //25 um thickness, start from 0, segment angle 360.
    logicAblator= new G4LogicalVolume(solidAblator, CH, "Ablator");
    //need to set CH material
    physiAblator= new G4PVPlacement(0, G4ThreeVector(0,0,Ablatorthickness/2),"Ablator",
                  logicAblator,
                  physiSample,
                  false,
                  0);

   solidCulayer= new G4Tubs("Culayer",
                  0, SampleSizeXY/2, Culayerthickness/2, 0, 360.*degree ); // name Culayer,
                  //inner radus 0, Outer radius half of SampleSizeXY,
                  //25 um thickness, start from 0, segment angle 360.
   logicCulayer= new G4LogicalVolume(solidCulayer, Cu, "Culayer");
   //need to set Cu material
   physiCulayer= new G4PVPlacement(0, G4ThreeVector(0,0,Ablatorthickness+Culayerthickness/2),"Culayer",
                  logicCulayer,
                 physiSample,
               false,
             0);

   solidAllayer= new G4Tubs("Steplayer",
                  0, SampleSizeXY/2, Stepthickness/2, 0, 360.*degree ); // name Allayer,
                  //inner radus 0, Outer radius half of SampleSizeXY,
                  //25 um thickness, start from 0, segment angle 360.
   logicAllayer= new G4LogicalVolume(solidAllayer, CH, "Steplayer");
   physiAllayer= new G4PVPlacement(0, G4ThreeVector(0,0,Ablatorthickness+Culayerthickness+Stepthickness/2),"Steplayer",
                  logicAllayer,
                 physiSample,
               false,
             0);

     }else if(SampleRadius>0)
     {
         solidSample = new G4Box("Sample",		//its name
               SampleRadius,SampleRadius,SampleRadius);//size

         logicSample= new G4LogicalVolume(solidSample,	//its solid
                  defaultMaterial,	//its material
                  "Sample");	//its name

         physiSample = new G4PVPlacement(0,			//no rotation
                 G4ThreeVector(0,0,0),	//at (0,0,0)
                 "Sample",	//its name
                 logicSample,	//its logical volume
                 physiWorld,	//its mother  volume
                 false,		//no boolean operation
                 0);		//copy number

         solidInnerlayer= new G4Sphere("Inner",
                        0, InnerRadius , 0,360*degree, 0, 180.*degree );
                        //inner radus 0, Outer radius ,
                        //25 um thickness, start from 0, segment angle 360.
         logicInnerlayer= new G4LogicalVolume(solidInnerlayer, MixCuCH, "Inner");
         physiInnerlayer= new G4PVPlacement(0, G4ThreeVector(0,0,0),"Inner",
                        logicInnerlayer,
                       physiSample,
                     false,
                   0);

         solid2ndlayer= new G4Sphere("Secondlayer",
                        InnerRadius, SecondRadius,  0, 360.*degree, 0, 180*degree);
         logic2ndlayer= new G4LogicalVolume(solid2ndlayer, DenseMixCuCH, "Secondlayer");
         physi2ndlayer= new G4PVPlacement(0, G4ThreeVector(0,0,0),"Secondlayer",
                        logic2ndlayer,
                       physiSample,
                     false,
                   0);
        if (ThirdRadius>SecondRadius)
        {
           solid3rdlayer= new G4Sphere("ThirdLayer",
                          SecondRadius, ThirdRadius,  0, 360.*degree, 0, 180*degree);
           logic3rdlayer= new G4LogicalVolume(solid3rdlayer, DenseCH, "ThirdLayer");

           physi3rdlayer= new G4PVPlacement(0, G4ThreeVector(0,0,0),"ThirdLayer",
                          logic3rdlayer,
                         physiSample,
                       false,
                     0);
           if (OuterRadius>ThirdRadius)
           {
               solidOuterlayer= new G4Sphere("Outerlayer",
                              ThirdRadius, OuterRadius,  0, 360.*degree, 0, 180*degree );
               logicOuterlayer= new G4LogicalVolume(solidOuterlayer, DenseMixCuCH, "Outerlayer");
               physiOuterlayer= new G4PVPlacement(0, G4ThreeVector(0,0,0),"Outerlayer",
                              logicOuterlayer,
                             physiSample,
                           false,
                         0);
           }
        }
     }
     G4cout<<"begin region cut"<<G4endl;

     if(SampleRegion)
       delete SampleRegion;


     SampleRegion = new G4Region("target");

     G4cout<<"region"<<G4endl;

     SampleRegion->SetProductionCuts(targetCuts);
     G4cout<<"cut"<<G4endl;

     SampleRegion->AddRootLogicalVolume(logicSample);

     G4cout<<"end region cut"<<G4endl;



  }
/*
  if (!phaseSpaceFlag) {
    //Diaphragm1

    solidDia1 = 0;  physiDia1 = 0;  logicDia1=0;

    if (Dia1Thickness > 0.)
      {
	solidDia1 = new G4Tubs("Diaphragm1",		//its name
			       DiaInnerSize/2,
			       Dia1SizeXY/2,
			       Dia1Thickness/2,
			       0,
			       360);//size


	logicDia1 = new G4LogicalVolume(solidDia1,	//its solid
					Dia1Material,	//its material
					"Diaphragm1");	//its name

	zRotPhiDia1.rotateX(AlphaDia1);
	G4double x,y,z;
	z = DistDia * std::cos(ThetaDia1);
	y =DistDia * std::sin(ThetaDia1);
	x = 0.*cm;
	physiDia1 = new G4PVPlacement(G4Transform3D(zRotPhiDia1,G4ThreeVector(x,y,z)),
				      "Diaphragm1",	//its name
				      logicDia1,	//its logical volume
				      physiWorld,	//its mother  volume
				      false,		//no boolean operation
				      0);		//copy number
      }

    //Diaphragm3

    solidDia3 = 0;  physiDia3 = 0;  logicDia3 =0;

    if (Dia3Thickness > 0.)
      {
	solidDia3 = new G4Tubs("Diaphragm3",
			       Dia3InnerSize/2,
			       Dia3SizeXY/2,
			       Dia3Thickness/2,
			       0,
			       360);


      logicDia3 = new G4LogicalVolume(solidDia3,	//its solid
				      Dia3Material,	//its material
				      "Diaphragm3");	//its name

      zRotPhiDia3.rotateX(AlphaDia3);
      G4double x,y,z;
      z = Dia3Dist * std::cos(ThetaDia3);
      y =Dia3Dist * std::sin(ThetaDia3);
      x = 0.*cm;
      physiDia3 = new G4PVPlacement(G4Transform3D(zRotPhiDia3,G4ThreeVector(x,y,z)),                                           "Diaphragm3",	//its name
				    logicDia3,	//its logical volume
				    physiWorld,	//its mother  volume
				    false,		//no boolean operation
				    0);		//copy number
      }
  }
  */
  G4cout << "finished target constrction" << G4endl;



  // Visualization attributes

  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes * yellow= new G4VisAttributes( G4Colour(255/255. ,255/255. ,51/255. ));
  G4VisAttributes * red= new G4VisAttributes( G4Colour(255/255. , 0/255. , 0/255. ));
  G4VisAttributes * blue= new G4VisAttributes( G4Colour(0/255. , 0/255. ,  255/255. ));
  G4VisAttributes * grayc= new G4VisAttributes( G4Colour(128/255. , 128/255. ,  128/255. ));
  G4VisAttributes * lightGray= new G4VisAttributes( G4Colour(178/255. , 178/255. ,  178/255. ));
  G4VisAttributes * green= new G4VisAttributes( G4Colour(0/255. , 255/255. ,  0/255. ));

  yellow->SetVisibility(true);
  yellow->SetForceSolid(true);
  red->SetVisibility(true);
//  red->SetForceSolid(true);
  blue->SetVisibility(true);
//  blue->SetForceSolid(true);
  green->SetVisibility(true);
  green->SetForceSolid(true);
  grayc->SetVisibility(true);
  grayc->SetForceSolid(true);
  lightGray->SetVisibility(true);
//  lightGray->SetForceSolid(true);
  simpleBoxVisAtt->SetVisibility(true); //so that it's not solid
  if (!phaseSpaceFlag && logicPixel) {
    logicPixel->SetVisAttributes(red); //modified!!!
    logicHPGe->SetVisAttributes(blue);

    logicDia1->SetVisAttributes(lightGray);
    logicDia3->SetVisAttributes(lightGray);

    logicOhmicNeg->SetVisAttributes(yellow);
    logicOhmicPos->SetVisAttributes(yellow);

    logicWindow->SetVisAttributes(green);

  }

  //logicAblator->SetVisAttributes(blue);
  //logicCulayer->SetVisAttributes(red);
  //logicAllayer->SetVisAttributes(lightGray);
  logicInnerlayer->SetVisAttributes(blue);
  logic2ndlayer->SetVisAttributes(red);
  logic3rdlayer->SetVisAttributes(lightGray);
  //logicOuterlayer->SetVisAttributes(lightGray);

  G4cout << "Finished Visualization" << G4endl;
  if(ConstructDetector)
  {
    logicBMXS->SetVisAttributes(simpleBoxVisAtt);
    logicBMXSCollimator->SetVisAttributes(blue);
    logicCartridge->SetVisAttributes(simpleBoxVisAtt);
    logicTeflonBlock->SetVisAttributes(green);
  }
  /*
  for(G4int i=0;i<15;i++ )
  {
    logicfilter[i]->SetVisAttributes(yellow);
    logicfilter[i]->SetForceSolid(true);
    logicIP[i]->SetVisAttributes(blue);
    logicIP[i]->SetForceSolid(true);

  }
  */



  logicSample->SetVisAttributes(simpleBoxVisAtt);

  if (sampleGranularity) logicSample->SetVisAttributes(simpleBoxVisAtt); // mandatory



  if (sampleGranularity)  logicGrain->SetVisAttributes(grayc);

  //always return the physical World

  PrintApparateParameters();

  G4cout << "finished construction" << G4endl;


  return physiWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::ConstructSDandField()
{
  if (!phaseSpaceFlag)
    {

      //
      // Sensitive Detectors
      //
      if (HPGeSD.Get() == 0)
	{
	  XrayFluoSD* SD = new XrayFluoSD ("HPGeSD",this);
	  HPGeSD.Put( SD );
	}
      G4SDManager::GetSDMpointer()->AddNewDetector(HPGeSD.Get());
      if (logicPixel)
	SetSensitiveDetector(logicPixel,HPGeSD.Get());
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::PrintApparateParameters()
{
  G4cout << "-----------------------------------------------------------------------"
	 << G4endl
	 << "The sample is a box whose size is: "
	 << G4endl
	 << SampleThickness/cm
	 << " cm * "
	 << SampleSizeXY/cm
	 << " cm * "
	 << SampleSizeXY/cm
	 << " cm"
	 << G4endl
	 <<" Material: " << logicSample->GetMaterial()->GetName()
	 <<G4endl;
  if (!phaseSpaceFlag) {
    G4cout <<"The Detector is a slice  " << DeviceThickness/(1.e-6*m) <<  " micron thick of " << pixelMaterial->GetName()
	   <<G4endl
	   << "The Anode is a slice " << OhmicPosThickness/mm << "mm thick of "<< OhmicPosMaterial->GetName()
	   <<G4endl;
      }
  G4cout <<"-------------------------------------------------------------------------"
	 << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::UpdateGeometry()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::Clean();
  G4LogicalVolumeStore::Clean();
  G4SolidStore::Clean();

//  if (sampleRegion)
  //  sampleRegion->RemoveRootLogicalVolume(logicSample);

  zRotPhiHPGe.rotateX(-1.*PhiHPGe);
  zRotPhiDia1.rotateX(-1.*AlphaDia1);
  zRotPhiDia3.rotateX(-1.*AlphaDia3);

  //Triggers a new call of Construct() and of all the geometry resets.
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::DeleteGrainObjects()
{
  if (sampleGranularity) {
    delete solidGrain;
    delete logicGrain;
    delete physiGrain;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XrayFluoDetectorConstruction::GetDetectorPosition() const
{
  G4double 	z = DistDe * std::cos(ThetaHPGe);
  G4double	y = DistDe * std::sin(ThetaHPGe);
  G4double	x = 0.*cm;
//G4ThreeVector position(x,y,z);
  G4ThreeVector position(0,0,0);

  return position;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSampleMaterial(G4String newMaterial)
{
    G4cout << "Material Change in Progress " << newMaterial << G4endl;
    sampleMaterial = materials->GetMaterial(newMaterial);
    logicSample->SetMaterial(sampleMaterial);
    PrintApparateParameters();
    //GeometryHasBeenModified is called by the messenger
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliMagF.h"
#include "AliESDfriend.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGeoGlobalMagField.h>
#endif

const int kNSlots = 2;
const float kHugeChi=1e9;
const int kAccBit=0x1<<15;
const int kTestBit=0x1<<14;

//-------------- cuts ------------->>>
const int kMinCls = 70;
const float kMaxDCAY = 3;
const float kMaxDCAZ = 3;
const float kMaxEta = 0.8;
const float kMinPt = 1.;
//-------------- cuts -------------<<<

TTree *dbgTreeCl=0,*dbgTreeTr=0;
TFile *DebugFile=0;
AliTPCclusterMI *clsODbg=0,*clsHDbg=0;
AliExternalTrackParam *trODbg=0,*trHDbg=0;
typedef struct {
  Float_t q2pt;
  Float_t chi2;
  Bool_t mtcO;
  Bool_t mtcH;
} dbgCl_t;
dbgCl_t dbgCl;

typedef struct {
  ULong64_t flgO;
  ULong64_t flgH;
  Float_t chi2mtc;
  Float_t chi2O;
  Float_t chi2H;
  Int_t   nclO;
  Int_t   nclH;
  Bool_t mtcO;
  Bool_t mtcH;
} dbgTr_t;
dbgTr_t dbgTr;
  
int currEvent=-1;
int currSlot=0;
TFile *flIn[kNSlots]={0};
AliESDEvent *esdEv[kNSlots]={0};
AliESDfriend *esdFr[kNSlots]={0};
TTree *esdTree[kNSlots]={0};

float GetRowX(int row);
void SetSlot(int s=0);
void PrintTrack(Int_t i);
void PrintTracks();
void ConnectFriends();
Int_t LoadESD(const char* path, Bool_t friends=kFALSE);
Int_t LoadEvent(Int_t iev);
Bool_t AcceptTrack(AliESDtrack* trc);
Float_t CompareTracks(AliExternalTrackParam &par0,AliExternalTrackParam &par1);
void ProcessEvent();
void FillSeedsInfo(int itr0,int itr1);
void CloseDbgTree();
void BookDbgTree(const char* dbgName);

void CompOfflHLT(const char* pathOffl, const char* pathHLT, const char* dbgName="dbgTree.root")
{
  //
  Bool_t friends=kTRUE;
  //
  BookDbgTree(dbgName);
  //
  SetSlot(0);
  int nevOffl = LoadESD(pathOffl,friends);
  SetSlot(1);
  int nevHLT  = LoadESD(pathHLT, friends);
  //
  if (nevOffl<0 || nevOffl!=nevHLT) {
    printf("No or mismatch in N events: %d %d\n",nevOffl,nevHLT);
    exit(1);
  }
  //
  for (currEvent=0;currEvent<nevOffl;currEvent++) {
    SetSlot(0);
    LoadEvent(currEvent);
    SetSlot(1);
    LoadEvent(currEvent);
    if (currEvent==0) esdEv[0]->InitMagneticField();
    ProcessEvent();
  }
  //
  CloseDbgTree();
}

//===========================================================
void ProcessEvent()
{
  int ntr0 = esdEv[0]->GetNumberOfTracks();
  int ntr1 = esdEv[1]->GetNumberOfTracks();
  if (!ntr0 || !ntr1) return;
  int bestMatch[ntr0]={0};
  float chiMatch[ntr0]={0};
  for (int itr0=0;itr0<ntr0;itr0++) {
    bestMatch[itr0] = -1;
    chiMatch[itr0] = kHugeChi;
    AliESDtrack* trc0 = esdEv[0]->GetTrack(itr0);    
    if (!AcceptTrack(trc0)) continue;
    *trODbg = *trc0->GetInnerParam();
    //
    for (int itr1=0;itr1<ntr1;itr1++) {
      AliESDtrack* trc1 = esdEv[1]->GetTrack(itr1);
      if (trc0->Charge()!=trc1->Charge()) continue;
      if (!AcceptTrack(trc1)) continue;
      //
      *trHDbg = *trc1->GetInnerParam();
      //
      float chi2 = CompareTracks(*trODbg,*trHDbg);
      if (chi2<chiMatch[itr0]) {
	chiMatch[itr0] = chi2;
	bestMatch[itr0] = itr1;
      }
      //
    }
    if (bestMatch[itr0]<0) continue;
    //
    //    double qy2x = par0.GetY()/par0.GetX()*par0.Charge();
    //    printf("Ev.%4d | tr %4d -> %4d | %f %+.4f (%+f %+f)\n",currEvent,itr0,bestMatch[itr0],chiMatch[itr0],qy2x,par0.GetY(),par0.GetX());
    AliESDtrack* trc1 = esdEv[1]->GetTrack(bestMatch[itr0]);
    //
    dbgTr.chi2mtc = dbgCl.chi2 = chiMatch[itr0];
    //
    dbgTr.chi2O = trc0->GetTPCchi2();
    dbgTr.flgO = trc0->GetStatus();
    dbgTr.nclO = trc0->GetTPCncls();
    //
    dbgTr.chi2H = trc1->GetTPCchi2();
    dbgTr.flgH = trc1->GetStatus();
    dbgTr.nclH = trc1->GetTPCncls();
    //
    dbgTreeTr->Fill();
    FillSeedsInfo(itr0,bestMatch[itr0]);
    //
  }  
}

//___________________________________________________________
void FillSeedsInfo(int itr0,int itr1)
{
  const AliESDtrack* trc0 = esdEv[0]->GetTrack(itr0);
  const AliESDtrack* trc1 = esdEv[1]->GetTrack(itr1);
  const AliESDfriendTrack* trf0 = trc0->GetFriendTrack();
  const AliESDfriendTrack* trf1 = trc1->GetFriendTrack();
  if (!trf0||!trf1) return;
  const AliTPCseed* sd0 = (AliTPCseed*)trf0->GetTPCseed();
  const AliTPCseed* sd1 = (AliTPCseed*)trf1->GetTPCseed();
  if (!sd0||!sd1) return;
  //
  dbgCl.mtcO = trc0->IsOn(AliESDtrack::kITSrefit);
  dbgCl.mtcH = trc1->IsOn(AliESDtrack::kITSrefit);  
  dbgCl.q2pt = trc0->GetParameter()[4];

  // fill clusters difference
  for (int ir=0;ir<159;ir++) {
    const AliTPCclusterMI* cls0 = sd0->GetClusterFast(ir);
    const AliTPCclusterMI* cls1 = sd1->GetClusterFast(ir);
    if (!cls0 || !cls1 || cls0->GetDetector()!=cls1->GetDetector()) continue;
    //
    *clsODbg = *cls0;
    *clsHDbg = *cls1;    
    dbgTreeCl->Fill();
    //
  }
}

//___________________________________________________________
Float_t CompareTracks(AliExternalTrackParam &par0,AliExternalTrackParam &par1)
{
  float chi2 = 2*kHugeChi;
  double alp0 = par0.GetAlpha(), alp1 = par1.GetAlpha();
  double dalpha = TMath::Abs(alp0-alp1);
  if (dalpha>25*TMath::DegToRad()) return chi2; // at least neighbouring sectors
  //
  if (dalpha>FLT_EPSILON && !par1.Rotate(alp0)) return chi2;
  if (TMath::Abs(par0.GetX()-par1.GetX())>FLT_EPSILON) {
    AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    Double_t xyz[3],b[3];
    par1.GetXYZ(xyz);
    fld->Field(xyz,b);
    if (!par1.PropagateToBxByBz(par0.GetX(),b)) return chi2;
  }
  //
  chi2 = par0.GetPredictedChi2(&par1);
  //
  return chi2;
}


//___________________________________________________________
Bool_t AcceptTrack(AliESDtrack* trc)
{
  if (trc->TestBit(kTestBit)) return trc->TestBit(kAccBit); // already tested
  trc->SetBit(kTestBit);
  trc->ResetBit(kAccBit);
  if (!trc->IsOn(AliESDtrack::kTPCin)) return kFALSE;
  if (trc->GetNcls(1)<kMinCls) return kFALSE;
  if (trc->Pt()<kMinPt) return kFALSE;
  double x = trc->GetInnerParam()->GetX();
  if (x<80 || x>85) return kFALSE;
  float dy,dz;
  if (TMath::Abs(trc->GetTPCInnerParam()->GetX())>2) return kFALSE;
  trc->GetImpactParametersTPC(dy,dz);
  if (TMath::Abs(dy)>kMaxDCAY || TMath::Abs(dz)>kMaxDCAZ) return kFALSE;
  double eta = trc->Eta();
  if (TMath::Abs(eta)>kMaxEta) return kFALSE;
  trc->SetBit(kAccBit);
  return kTRUE;
}

//===========================================================
Int_t LoadESD(const char* path, Bool_t friends)
{
  flIn[currSlot] = TFile::Open(path);
  if (!flIn[currSlot]) {
    printf("Failed to open %s\n",path);
    return -1;
  }
  //
  esdTree[currSlot] = (TTree*) flIn[currSlot]->Get("esdTree");
  if (!esdTree[currSlot]) {
    printf("No ESDtree found in %s\n",path);
    return -1;
  }
  //
  int nev = (int)esdTree[currSlot]->GetEntries();
  printf("Loaded esdTree with %d entries for slot %d from %s\n",nev,currSlot,path);
  esdEv[currSlot] = new AliESDEvent();
  esdEv[currSlot]->ReadFromTree(esdTree[currSlot]);
  if (friends) ConnectFriends();
  //
  return nev;
}

Int_t LoadEvent(Int_t iev)
{
  if (!esdEv[currSlot]) return -1; 
  esdEv[currSlot]->Reset();
  if (esdFr[currSlot]) esdFr[currSlot]->Reset();
  esdTree[currSlot]->GetEntry(iev);
  if (esdFr[currSlot]) esdEv[currSlot]->SetESDfriend(esdFr[currSlot]);
  esdEv[currSlot]->ConnectTracks();
  return 0;
}


void PrintTracks()
{
  for (int i=0;i<esdEv[currSlot]->GetNumberOfTracks();i++) PrintTrack(i);
}

void PrintTrack(Int_t i)
{
  AliESDtrack* tr = esdEv[currSlot]->GetTrack(i);
  if (!tr) return;
  double p[3];
  tr->GetPxPyPz(p);
  printf("%3d: itsRF%d tpcRF%d trdRF%d tofRF%d | P:%6.2f Px:%+6.2f Py:%+6.2f Pz:%+6.2f\n",
	 i,
	 tr->IsOn(AliESDtrack::kITSrefit),
	 tr->IsOn(AliESDtrack::kTPCrefit),
	 tr->IsOn(AliESDtrack::kTRDrefit),
	 tr->IsOn(AliESDtrack::kTOFrefit),
	 tr->GetP(),p[0],p[1],p[2]);
}


void ConnectFriends()
{
  // Connect the friends tree as soon as available.
  //
  // Handle the friends first
  //
  if (!esdTree[currSlot]->FindBranch("ESDfriend.")) {
    // Try to add ESDfriend. branch as friend
      TString esdFriendTreeFName;
      esdFriendTreeFName = (esdTree[currSlot]->GetCurrentFile())->GetName();    
      TString basename = gSystem->BaseName(esdFriendTreeFName);
      Int_t index = basename.Index("#")+1;
      basename.Remove(index);
      basename += "AliESDfriends.root";
      TString dirname = gSystem->DirName(esdFriendTreeFName);
      dirname += "/";
      esdFriendTreeFName = dirname + basename;
      //
      TTree* cTree = esdTree[currSlot]->GetTree();
      if (!cTree) cTree = esdTree[currSlot];      
      cTree->AddFriend("esdFriendTree", esdFriendTreeFName.Data());
      cTree->SetBranchStatus("ESDfriend.", 1);
      esdFr[currSlot] = (AliESDfriend*)(esdEv[currSlot]->FindListObject("AliESDfriend"));
      if (esdFr[currSlot]) cTree->SetBranchAddress("ESDfriend.", &esdFr[currSlot]);
  }
}


//__________________________________________________
void SetSlot(int s)
{
  if (s>=kNSlots) {
    printf("max %d slots supported\n",kNSlots);
    exit(1);
  }
  currSlot = s;
}

float GetRowX(int row)
{
  const float kPadRowX[159] = {
    85.225, 85.975, 86.725, 87.475, 88.225, 88.975, 89.725, 90.475, 91.225, 91.975, 92.725, 93.475, 94.225, 94.975, 95.725,
    96.475, 97.225, 97.975, 98.725, 99.475,100.225,100.975,101.725,102.475,103.225,103.975,104.725,105.475,106.225,106.975,
    107.725,108.475,109.225,109.975,110.725,111.475,112.225,112.975,113.725,114.475,115.225,115.975,116.725,117.475,118.225,
    118.975,119.725,120.475,121.225,121.975,122.725,123.475,124.225,124.975,125.725,126.475,127.225,127.975,128.725,129.475,
    130.225,130.975,131.725,135.100,136.100,137.100,138.100,139.100,140.100,141.100,142.100,143.100,144.100,145.100,146.100,
    147.100,148.100,149.100,150.100,151.100,152.100,153.100,154.100,155.100,156.100,157.100,158.100,159.100,160.100,161.100,
    162.100,163.100,164.100,165.100,166.100,167.100,168.100,169.100,170.100,171.100,172.100,173.100,174.100,175.100,176.100,
    177.100,178.100,179.100,180.100,181.100,182.100,183.100,184.100,185.100,186.100,187.100,188.100,189.100,190.100,191.100,
    192.100,193.100,194.100,195.100,196.100,197.100,198.100,199.350,200.850,202.350,203.850,205.350,206.850,208.350,209.850,
    211.350,212.850,214.350,215.850,217.350,218.850,220.350,221.850,223.350,224.850,226.350,227.850,229.350,230.850,232.350,
    233.850,235.350,236.850,238.350,239.850,241.350,242.850,244.350,245.850
  };
  return kPadRowX[row];
}

//______________________________________
void BookDbgTree(const char* dbgName)
{
  //
  DebugFile = new TFile(dbgName,"recreate");
  //
  dbgTreeCl = new TTree("clsTr","clsTree");
  clsODbg = new AliTPCclusterMI();
  clsHDbg = new AliTPCclusterMI();  
  dbgTreeCl->Branch("clsO","AliTPCclusterMI",&clsODbg);
  dbgTreeCl->Branch("clsH","AliTPCclusterMI",&clsHDbg);
  dbgTreeCl->Branch("dbgCl",&dbgCl,"q2pt/F:chi2/F:mtcO/O:mtcH/O");
  //
  dbgTreeCl->SetAlias("dy","clsO.fY-clsH.fY");
  dbgTreeCl->SetAlias("dz","clsO.fZ-clsH.fZ");
  dbgTreeCl->SetAlias("qy2x","clsO.fY/clsO.fX*sign(q2pt)");
  dbgTreeCl->SetAlias("sect","(clsO.fDetector%18)");
  //
  dbgTreeTr = new TTree("trTr","TracksTree");
  trODbg = new AliExternalTrackParam();
  trHDbg = new AliExternalTrackParam();
  dbgTreeTr->Branch("trO","AliExternalTrackParam",&trODbg);
  dbgTreeTr->Branch("trH","AliExternalTrackParam",&trHDbg);
  dbgTreeTr->Branch("dbgTr",&dbgTr,"flgO/l:flgH/l:chi2mtc/F:chi2O/F:chi2H/F:nclO/I:nclH/I:mtcO/O:mtcH/O");
  //
}

//______________________________________
void CloseDbgTree()
{
  DebugFile->cd();
  dbgTreeCl->Write();
  delete dbgTreeCl;
  dbgTreeTr->Write();
  delete dbgTreeTr;
  //
  DebugFile->Close();
  delete DebugFile;
}

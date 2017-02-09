
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliMagF.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGeoGlobalMagField.h>
#endif

const int kNSlots = 2;
const float kHugeChi=1e9;
const int kAccBit=0x1<<15;
const int kTestBit=0x1<<14;

int currSlot=0;
TFile *flIn[kNSlots]={0};
AliESDEvent *esdEv[kNSlots]={0};
AliESDfriend *esdFr[kNSlots]={0};
TTree *esdTree[kNSlots]={0};

void SetSlot(int s=0);
void PrintTrack(Int_t i);
void PrintTracks();
void ConnectFriends();
Int_t LoadESD(const char* path, Bool_t friends=kFALSE);
Int_t LoadEvent(Int_t iev);
Bool_t AcceptTrack(AliESDtrack* trc);
Float_t CompareTracks(AliESDtrack* trc0, AliESDtrack* trc1);
void ProcessEvent();


void CompOfflHLT(const char* pathOffl, const char* pathHLT, Bool_t friends=kTRUE)
{
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
  for (int iev=0;iev<nevOffl;iev++) {
    SetSlot(0);
    LoadEvent(iev);
    SetSlot(1);
    LoadEvent(iev);
    if (iev==0) esdEv[0]->InitMagneticField();
    ProcessEvent();
  }
  //
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
    for (int itr1=0;itr1<ntr1;itr1++) {
      AliESDtrack* trc1 = esdEv[1]->GetTrack(itr1);
      if (trc0->Charge()!=trc1->Charge()) continue;
      if (!AcceptTrack(trc1)) continue;
      //
      float chi2 = CompareTracks(trc0,trc1);
      if (chi2<chiMatch[itr0]) {
	chiMatch[itr0] = chi2;
	bestMatch[itr0] = itr1;
      }
      //
    }
    printf("%d -> %d | %f\n",itr0,bestMatch[itr0],chiMatch[itr0]);
  }
  
}

//___________________________________________________________
Float_t CompareTracks(AliESDtrack* trc0, AliESDtrack* trc1)
{
  float chi2 = 2*kHugeChi;
  AliExternalTrackParam* par0 = (AliExternalTrackParam*)trc0->GetInnerParam();
  AliExternalTrackParam* par1 = (AliExternalTrackParam*)trc1->GetInnerParam();
  double alp0 = par0->GetAlpha(), alp1 = par1->GetAlpha();
  double dalpha = TMath::Abs(alp0-alp1);
  if (dalpha>25*TMath::DegToRad()) return chi2; // at least neighbouring sectors
  //
  if (dalpha>1e-3 && !par1->Rotate(alp0)) return chi2;
  if (TMath::Abs(par0->GetX()-par1->GetX())>1e-3) {
    AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
    Double_t xyz[3],b[3];
    par1->GetXYZ(xyz);
    fld->Field(xyz,b);
    if (!par1->PropagateToBxByBz(par0->GetX(),b)) return chi2;
  }
  //
  chi2 = par0->GetPredictedChi2(par1);
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
  if (trc->GetNcls(1)<70) return kFALSE;
  float dy,dz;
  if (TMath::Abs(trc->GetTPCInnerParam()->GetX())<2) return kFALSE;
  trc->GetImpactParametersTPC(dy,dz);
  if (TMath::Abs(dy)>3 || TMath::Abs(dz)>3) return kFALSE;
  double eta = trc->Eta();
  if (TMath::Abs(eta)>0.8) return kFALSE;
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

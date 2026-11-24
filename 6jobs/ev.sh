export LD_LIBRARY_PATH=/lcrc/project/eic/tgibbons/xxhash-install/lib:/lcrc/project/eic/tgibbons/gif-lib:$LD_LIBRARY_PATH
cd /lcrc/project/eic/tgibbons/G4ChargeExchange/build/root
cat > /tmp/ev.C <<'EOF'
void ev(const char* d){
  TString a=TString(d)+"/cuts/cut_angles_all.root", b=TString(d)+"/cuts/cut_angles_all_1.root";
  Long64_t n0=-1,n1=-1;
  { TFile* f=TFile::Open(a); if(f && !f->IsZombie()){ TTree* t=(TTree*)f->Get("AnglesRaw"); n0=t?t->GetEntries():-2; f->Close(); } }
  { TFile* f=TFile::Open(b); if(f && !f->IsZombie()){ TTree* t=(TTree*)f->Get("AnglesRaw"); n1=t?t->GetEntries():-2; f->Close(); } }
  printf("%s base=%lld cont=%lld\n", d, n0, n1);
}
EOF
for d in newPhysics_801_100m_*/*GeV/*block; do
  [ -e $d/cuts/cut_angles_all_1.root ] && root -l -b -q "/tmp/ev.C(\"$d\")" 2>&1 | tail -4
done | sed 's/newPhysics_801_100m_//'

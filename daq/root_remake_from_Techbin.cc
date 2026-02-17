#include <iostream>
#include <fstream>
#include <ctime>
#include <filesystem>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "/home/kinoshita/hm_includes/KR_math.cc"

#define RECORD_PULSE_SHAPE 1
#define BIN_NANOSEC 2
#define ABS_INTEGRAL 0

using namespace std;

Short_t reversendian(Short_t x){
    Short_t ux = x;
    Short_t dx = x;
    ux <<= 8;
    ux &= 0xff00;
    dx >>= 8;
    dx &= 0x00ff;
    Short_t y = ux | dx;
    return y;
}

class Techheader {
    private:
        bool waveFlag;
        Long64_t TDCbit;
        Short_t TDCFPbit;
        Double_t TDC;
        Short_t CH;
        //Short_t QDC;
        Short_t points;
    public:
        Techheader(Short_t a, Short_t b, Short_t c, Short_t d, Short_t e, Short_t f){
            Short_t templ1 = reversendian(a);
            Short_t templ2 = reversendian(b);
            Short_t templ3 = reversendian(c);
            Short_t templ4 = reversendian(d);
            Short_t templ5 = reversendian(e);
            Short_t templ6 = reversendian(f);
            Short_t headerShort[6] = {templ1, templ2, templ3, templ4, templ5, templ6};

            waveFlag = (headerShort[0]>>15)&0x0001;

            TDCbit = 0;
            for (int i=0;i<=3;i++){
                TDCbit <<= 16;
                TDCbit &= 0xffffffffffff0000;
                TDCbit |= (headerShort[i] & 0x000000000000ffff);
            }
            TDCbit <<= 1;
            TDCbit >>= 10;
            TDCbit &= 0x003fffffffffffff;

            TDCFPbit = headerShort[3];
            TDCFPbit <<= 7;
            TDCFPbit >>= 8;
            TDCFPbit &= 0x00ff;

            TDC = ((Double_t)TDCbit) * 1.0 * 1e-9 + ((Double_t)TDCFPbit) * 3.90625 * 1e-12;

            CH = headerShort[3] & 0x0001;
            CH <<= 2;
            CH |= (headerShort[4]>>13)&0x0007;
            CH++;

            points = headerShort[5];
        }
        ~Techheader(){}
        bool GetWaveFlag(){
            return waveFlag;
        }
        Long64_t GetTDCbit(){
            return TDCbit;
        }
        Short_t GetTDCFPbit(){
            return TDCFPbit;
        }
        Double_t GetTDC(){
            return TDC;
        }
        Short_t GetCH(){
            return CH;
        }
        Short_t GetPoints(){
            return points;
        }
        bool CheckSameTDC(Techheader target){
            if (TDCbit == target.GetTDCbit() && TDCFPbit == target.GetTDCFPbit()){
                return true;
            }
            else{
                return false;
            }
        }
        void output_allinfo(){
            cout << "Wave Flag : " << waveFlag << endl;
            cout << "TDC : " << TDC << endl;
            cout << "CH : " << CH << endl;
            cout << "Points : " << points << endl;
        }
};

Int_t Convert_TechbinROOT(TString infilename, TString outfilename = "",Int_t ADCmin = 0, Int_t ADCmax = 230, Int_t ADCtailmin = 50, Int_t Pedestalpoints = 50){
    //Set up

    gStyle -> SetOptStat(0);
    gStyle -> SetLabelSize(0.04, "XYZ");
    gStyle -> SetTitleSize(0.05, "T");
    gStyle -> SetTitleSize(0.04, "XYZ");
    gStyle -> SetTitleOffset(1.1, "X");
    gStyle -> SetTitleOffset(1.2, "Y");
    gStyle -> SetPadLeftMargin(0.14);
    gStyle -> SetPadBottomMargin(0.14);

    FILE *fip;
    fip = fopen(infilename.Data(),"rb");

    if (fip == NULL){
        cout << "Error : Could not open file" << endl;
        cout << infilename << endl;
        return -1;
    }
    TString ofile_name;
    if (outfilename == TString("")){
        ofile_name = Form("%s.root",infilename.Remove(infilename.Length()-4).Data());
    }
    else{
        ofile_name = outfilename;
    }

    TFile* fop = new TFile(ofile_name.Data(),"RECREATE");

    Int_t    EventID;
    UShort_t Ch;
    Double_t Pedestal;
    Double_t PedestalSigma;
    Double_t ADC;
    Double_t ADCtail;
    Double_t PeakHeight;
    Double_t PeakTime;
    Double_t StartTime;
    Int_t    StartBin;
    Double_t LiveT;
    Double_t TDC;

    TTree *tree_temp = new TTree("tree_temp", "");
    tree_temp -> Branch("EventID",       &EventID,       "EventID/I");
    tree_temp -> Branch("Ch",            &Ch,            "Ch/s");
    tree_temp -> Branch("Pedestal",      &Pedestal,      "Pedestal/D");
    tree_temp -> Branch("PedestalSigma", &PedestalSigma, "PedestalSigma/D");
    tree_temp -> Branch("ADC",           &ADC,           "ADC/D");
    tree_temp -> Branch("ADCtail",       &ADCtail,       "ADCtail/D");
    tree_temp -> Branch("PeakHeight",    &PeakHeight,    "PeakHeight/D");
    tree_temp -> Branch("PeakTime",      &PeakTime,      "PeakTime/D");
    tree_temp -> Branch("StartTime",     &StartTime,     "StartTime/D");
    tree_temp -> Branch("StartBin",      &StartBin,      "StartBin/I");
    tree_temp -> Branch("LiveT",         &LiveT,         "LiveT/D");
    tree_temp -> Branch("TDC",           &TDC,           "TDC/D");

    TTree *tree_out = new TTree("tree", "");
    tree_out -> Branch("EventID",       &EventID,       "EventID/I");
    tree_out -> Branch("Ch",            &Ch,            "Ch/s");
    tree_out -> Branch("Pedestal",      &Pedestal,      "Pedestal/D");
    tree_out -> Branch("PedestalSigma", &PedestalSigma, "PedestalSigma/D");
    tree_out -> Branch("ADC",           &ADC,           "ADC/D");
    tree_out -> Branch("ADCtail",       &ADCtail,       "ADCtail/D");
    tree_out -> Branch("PeakHeight",    &PeakHeight,    "PeakHeight/D");
    tree_out -> Branch("PeakTime",      &PeakTime,      "PeakTime/D");
    tree_out -> Branch("StartTime",     &StartTime,     "StartTime/D");
    tree_out -> Branch("StartBin",      &StartBin,      "StartBin/I");
    tree_out -> Branch("LiveT",         &LiveT,         "LiveT/D");
    tree_out -> Branch("TDC",           &TDC,           "TDC/D");

    fseek(fip,0,SEEK_END);
    Long_t totalBytes = ftell(fip);
    cout << "Total Binary File Bytes:" << totalBytes << " Byte" << endl;
    fseek(fip,0,SEEK_SET);

    Long_t readBytes = 0;

    Double_t pretime[8] = {0.0};
    Long64_t TDCbit = 0;
    Short_t TDCFPbit = 0;
    Int_t preID = -1;
    Int_t Eventcounter = 0;
    Int_t datapoint = 0;

    vector<Double_t> VCh;
    vector<Double_t> VLiveT;

    while (readBytes < totalBytes){
        Short_t Bheader[6];

        if (fread(&Bheader,sizeof(Short_t),6,fip) < 6){
            cout << "Error : Failed to read file" << endl;
            return -1;
        }
        readBytes += sizeof(Short_t)*6;

        Techheader head = Techheader(Bheader[0],Bheader[1],Bheader[2],Bheader[3],Bheader[4],Bheader[5]);

        //head.output_allinfo();
        Ch = head.GetCH();
        TDCbit = head.GetTDCbit();
        TDCFPbit = head.GetTDCFPbit();
        TDC = head.GetTDC();

        datapoint = head.GetPoints();

        EventID = Eventcounter;
        LiveT = TDC - pretime[Ch-1];

        pretime[Ch-1] = TDC;

        VCh.push_back(Ch);
        VLiveT.push_back(TDC);

        Char_t hCharWave[4];
        if (fread(&hCharWave,sizeof(Char_t),4,fip) < 1){
            cout << "Error : Failed to read file" << endl;
            return -1;
        }
        readBytes += sizeof(Char_t)*4;

        Short_t wave_data[datapoint];
        Short_t wave_data_bc[datapoint];
        if (fread(&wave_data_bc,sizeof(Short_t),datapoint,fip) < 1){
            cout << "Error : Failed to read file" << endl;
            return -1;
        }
        readBytes += sizeof(Short_t)*datapoint;

        TH1D *h1 = new TH1D(Form("h%d_Ch%d",EventID,Ch),"",datapoint,0,datapoint*BIN_NANOSEC);
        for (int k=0;k<datapoint;k++){
            wave_data[k] = reversendian(wave_data_bc[k]);
            h1 -> SetBinContent(k+1,wave_data[k]);
        }

        Double_t tPedestal = 0;
        for (Int_t k = 1; k <= Pedestalpoints; k++) {
            tPedestal += h1 -> GetBinContent(k);
        }
        tPedestal /= (Double_t)Pedestalpoints;

        // subtract pedestal
        for (Int_t k = 1; k <= datapoint; k++) {
            h1 -> SetBinContent(k, h1 -> GetBinContent(k) - tPedestal);
        }

        // get pedestal fluctuation
        Double_t tPedestalSigma = 0;
        for (Int_t k = 1; k <= Pedestalpoints; k++) {
            tPedestalSigma += TMath::Power(h1 -> GetBinContent(k), 2);
        }
        tPedestalSigma = TMath::Sqrt(tPedestalSigma / Pedestalpoints);

        Double_t tStartTime = 0;
        Int_t tStartBin = 0;
        
        // get start time
        for (Int_t k = h1 -> GetMaximumBin(); 0 < k; k--) {
            if ((h1 -> GetBinContent(k-1) < 2.0 * tPedestalSigma && h1 -> GetBinContent(k) < 2.0 * tPedestalSigma) && (h1 -> GetBinContent(k) < h1 -> GetBinContent(k+1) && h1 -> GetBinContent(k) < h1 -> GetBinContent(k+2))) {
                tStartTime = h1 -> GetBinCenter(k);
                tStartBin  = k;
                break;
            }
        }

        Pedestal   = tPedestal;
        PedestalSigma = tPedestalSigma;
        PeakHeight = h1 -> GetMaximum();
        PeakTime   = h1 -> GetBinCenter(h1 -> GetMaximumBin());
        StartTime  = tStartTime;
        StartBin   = tStartBin;

        #if ABS_INTEGRAL
        for (int i=h1 -> FindBin(tStartTime + ADCmin); i<=h1 -> FindBin(tStartTime + ADCmax)-1; i++)
        {
            ADC += h1 -> GetBinContent(i);
            if (i >= h1 -> FindBin(tStartTime + ADCtailmin))
            {
                ADCtail += fabs(h1 -> GetBinContent(i));
            }
        }
        #else
        ADC        = h1 -> Integral(h1 -> FindBin(tStartTime + ADCmin), h1 -> FindBin(tStartTime + ADCmax) - 1);
        ADCtail    = h1 -> Integral(h1 -> FindBin(tStartTime + ADCtailmin), h1 -> FindBin(tStartTime + ADCmax) - 1);
        #endif

        #if RECORD_PULSE_SHAPE
        h1 -> Write();
        #endif

        delete h1;
        Eventcounter++;
        tree_temp -> Fill();
        if (Eventcounter % 10000 == 0){
            cout << "Converted " << Eventcounter << " Event" << endl;
        }
    }
    cout << "Converted " << Eventcounter << " Event" << endl;

    vector<Int_t> indexV = KR::HeapSort(VLiveT,VCh);
    for (int i=0; i<static_cast<int>(indexV.size()); i++)
    {
        tree_temp -> GetEntry(indexV[i]);
        tree_out -> Fill();
        if ((i+1) % 10000 == 0)
        {
            cout << "Sorted & Filled " << i+1 << " Events" << endl;
        }
    }
    cout << "Sorted & Filled " << static_cast<int>(indexV.size()) << " Event" << endl;

    delete tree_temp;

    tree_out -> Write();
    fop -> Close();

    return 0;
}

Int_t extract_channel(TString infilename, Int_t targetChNo, TString outfilename = "")
{
    TString ifile_name = infilename;
    TString ofile_name;
    if (outfilename == TString("")){
        ofile_name = Form("%s.root",(infilename.Remove(infilename.Length()-5)+Form("_Ch%d",targetChNo)).Data());
    }
    else{
        ofile_name = outfilename;
    }
    TFile* fip = new TFile(ifile_name.Data());
    TFile* fop = new TFile(ofile_name.Data(),"RECREATE");

    Int_t    inEventID;
    UShort_t inCh;
    Double_t inPedestal;
    Double_t inPedestalSigma;
    Double_t inADC;
    Double_t inADCtail;
    Double_t inPeakHeight;
    Double_t inPeakTime;
    Double_t inStartTime;
    Int_t    inStartBin;
    Double_t inLiveT;
    Double_t inTDC;

    fip -> cd();

    TTree *tree_in = (TTree*)fip -> Get("tree");
    tree_in -> SetBranchAddress("EventID",       &inEventID);
    tree_in -> SetBranchAddress("Ch",            &inCh);
    tree_in -> SetBranchAddress("Pedestal",      &inPedestal);
    tree_in -> SetBranchAddress("PedestalSigma", &inPedestalSigma);
    tree_in -> SetBranchAddress("ADC",           &inADC);
    tree_in -> SetBranchAddress("ADCtail",       &inADCtail);
    tree_in -> SetBranchAddress("PeakHeight",    &inPeakHeight);
    tree_in -> SetBranchAddress("PeakTime",      &inPeakTime);
    tree_in -> SetBranchAddress("StartTime",     &inStartTime);
    tree_in -> SetBranchAddress("StartBin",      &inStartBin);
    tree_in -> SetBranchAddress("LiveT",         &inLiveT);
    tree_in -> SetBranchAddress("TDC",           &inTDC);

    Int_t    EventID;
    Double_t Pedestal;
    Double_t PedestalSigma;
    Double_t ADC;
    Double_t ADCtail;
    Double_t PeakHeight;
    Double_t PeakTime;
    Double_t StartTime;
    Int_t    StartBin;
    Double_t LiveT;
    Double_t TDC;

    fop -> cd();

    TTree *tree_out = new TTree("tree", "");
    tree_out -> Branch("EventID",       &EventID,       "EventID/I");
    tree_out -> Branch("Pedestal",      &Pedestal,      "Pedestal/D");
    tree_out -> Branch("PedestalSigma", &PedestalSigma, "PedestalSigma/D");
    tree_out -> Branch("ADC",           &ADC,           "ADC/D");
    tree_out -> Branch("ADCtail",       &ADCtail,       "ADCtail/D");
    tree_out -> Branch("PeakHeight",    &PeakHeight,    "PeakHeight/D");
    tree_out -> Branch("PeakTime",      &PeakTime,      "PeakTime/D");
    tree_out -> Branch("StartTime",     &StartTime,     "StartTime/D");
    tree_out -> Branch("StartBin",      &StartBin,      "StartBin/I");
    tree_out -> Branch("LiveT",         &LiveT,         "LiveT/D");
    tree_out -> Branch("TDC",           &TDC,         "TDC/D");

    Int_t counter = 0;
    for (int i=0; i<tree_in -> GetEntries(); i++)
    {
        fip -> cd();
        tree_in -> GetEntry(i);
        if (inCh != targetChNo)
        {
            //cout << "Passed EventID:" << EventID << endl;
            continue;
        }

        EventID = inEventID;
        Pedestal = inPedestal;
        PedestalSigma = inPedestalSigma;
        ADC = inADC;
        ADCtail = inADCtail;
        PeakHeight = inPeakHeight;
        PeakTime = inPeakTime;
        StartTime = inStartTime;
        StartBin = inStartBin;
        LiveT = inLiveT;
        TDC = inTDC;

        fop -> cd();
        tree_out -> Fill();

        #if RECORD_PULSE_SHAPE
        TH1D *h1 = (TH1D*)fip -> Get(Form("h%d_Ch%d",EventID,targetChNo));
        if (h1 == NULL){
            cout << "Not Found histgram:" << Form("h%d_Ch%d",EventID,targetChNo) << endl;
        }
        else{
            h1 -> SetName(Form("h%d",EventID));
            h1 -> Write();
        }
        delete h1;
        #endif

        counter++;
        if(counter%10000 == 0){
            cout << "Extracted " << counter << " events" << endl;
        }
    }
    cout << "Extracted " << counter << " Events" << endl;
    fop -> cd();
    tree_out -> Write();
    fip -> Close();
    fop -> Close();
    return 0;
}
void ForExtract_Channel(TString infilename, Int_t max, Int_t min = 1)
{
    for (int i=min; i<=max; i++)
    {
        extract_channel(infilename, i);
        cout << "Outputed : Channel " << i << endl;
    }
}

Int_t main(Int_t argc,Char_t **argv){
    // integral range of pulse shape
    const Double_t ADC_MIN  = 0;   // ns, from StartTime
    const Double_t ADC_MAX  = 100; // ns, from StartTime
    const Double_t ADCtail_MIN  = 60; // ns, from StartTime

    const Int_t date = 20260107;
    const TString foldername = Form("../kaiseki/data/%d/",date);
    //const TString foldername = "";
    const Int_t testNo = 1;
    const TString chara = "Cs";
    TString fin_name  = Form("%s%d_test%d_%s.bin",foldername.Data(),date,testNo,chara.Data());
    TString fout_name = Form("%s/root_file/%.0lf_root_%d_test%d_%s_%.0lf.root",foldername.Data(),ADCtail_MIN,date,testNo,chara.Data(),ADC_MAX);
    TString info_name = Form("%s/info_file/%.0lf_root_%d_test%d_%s_%.0lf.txt",foldername.Data(),ADCtail_MIN,date,testNo,chara.Data(),ADC_MAX);

    Convert_TechbinROOT(fin_name,fout_name,ADC_MIN,ADC_MAX,ADCtail_MIN,30);
    cout << endl << "Finished convert binary to ROOT" << endl << endl;
    
    TString exfin_name = fout_name;
    //TString exfout_name = Form("%s%.0lf_root_%d_test%d_Cs_%.0lf_xxx.root",foldername.Data(),ADCtail_MIN,date,testNo,ADC_MAX);
    ForExtract_Channel(exfin_name, 8);
    cout << endl << "Finished extract ch" << endl;
}
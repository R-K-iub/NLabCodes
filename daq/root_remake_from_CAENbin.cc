#include <iostream>
#include <fstream>
#include <ctime>
#include <filesystem>
#include <stdio.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TDirectory.h"

#define RECORD_PULSE_SHAPE 1
#define TerminalOutStep 10000

using namespace std;

class CAENheader {
    private:
        uint32_t eventsize;
        uint32_t boardID;
        uint32_t options;
        uint32_t channelmask;
        uint32_t eventcount;
        double timetag;
        uint32_t channelnum;
    public:
        CAENheader(uint32_t a, uint32_t b, uint32_t c, uint32_t d){
            eventsize = a&0x0FFFFFFF;
            boardID = (b>>27)&0x0000001F;
            options = (b>>8)&0x0000FFFF;
            channelmask = (b&0x000000FF)|((c>>16)&0x0000FF00);
            eventcount = c&0x00FFFFFF;
            timetag = d*8e-9;
            channelnum = 0;
            for (int i=0;i<channelmask;i++){
                if ((channelmask>>i)&0x1){
                    channelnum++;
                }
            }
        }
        ~CAENheader(){}
        uint32_t GetEventSize(){
            return eventsize;
        }
        uint32_t GetBoardID(){
            return boardID;
        }
        uint32_t GetOptions(){
            return options;
        }
        uint32_t GetChannelMask(){
            return channelmask;
        }
        uint32_t GetChannelNum(){
            return channelnum;
        }
        uint32_t GetEventCount(){
            return eventcount;
        }
        double GetTimeTag(){
            return timetag;
        }
        bool CheckChannelUse(int t){
            bool targetFlag = false;
            if (((channelmask>>t)&0x1) == 1)
                targetFlag = true;
            return targetFlag;
        }
        void OutputInfo(){
            cout << "Event Size : " << eventsize << endl;
            cout << "Borad ID : " << boardID << endl;
            cout << "Options : " << options << endl;
            cout << "Channel Mask : " << channelmask << endl;
            cout << "Channel Number : " << channelnum << endl;
            cout << "Event Count : " << eventcount << endl;
            cout << "Trigger Time Tag : " << timetag << endl;
        }
};

Int_t Convert_CAENbinROOT(TString infilename, TString outfilename = "", TString Infofilename = "", Int_t ADCmin = 0, Int_t ADCmax = 230, Int_t ADCtailmin = 50, Int_t Pedestalpoints = 50){
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

    ofstream infoos;
    bool InfoFlag = false;
    if (Infofilename != TString(""))
    {
        infoos.open(Infofilename.Data());
        InfoFlag = true;
    }

    Int_t    EventID;
    UShort_t Channels;
    UShort_t UseChannel;
    Double_t Pedestal[16];
    Double_t PedestalSigma[16];
    Double_t ADC[16];
    Double_t ADCtail[16];
    Double_t PeakHeight[16];
    Double_t PeakTime[16];
    Double_t StartTime[16];
    Int_t    StartBin[16];
    Double_t LiveT = 0.0;

    TTree *tree_out = new TTree("tree", "");
    tree_out -> Branch("EventID",       &EventID,       "EventID/I");
    tree_out -> Branch("Channels",      &Channels,      "Channels/s");
    tree_out -> Branch("UseChannel",    &UseChannel,    "UseChannel/s");
    tree_out -> Branch("Pedestal",      Pedestal,      "Pedestal[Channels]/D");
    tree_out -> Branch("PedestalSigma", PedestalSigma, "PedestalSigma[Channels]/D");
    tree_out -> Branch("ADC",           ADC,           "ADC[Channels]/D");
    tree_out -> Branch("ADCtail",       ADCtail,       "ADCtail[Channels]/D");
    tree_out -> Branch("PeakHeight",    PeakHeight,    "PeakHeight[Channels]/D");
    tree_out -> Branch("PeakTime",      PeakTime,      "PeakTime[Channels]/D");
    tree_out -> Branch("StartTime",     StartTime,     "StartTime[Channels]/D");
    tree_out -> Branch("StartBin",      StartBin,      "StartBin[Channels]/I");
    tree_out -> Branch("LiveT",         &LiveT,         "LiveT/D");

    fseeko64(fip,0,SEEK_END);
    loff_t totalBytes8 = ftello64(fip);
    Long64_t totalBytes = static_cast<Long64_t>(totalBytes8);
    cout << "Total Binary File Bytes:" << totalBytes << " Byte" << endl;
    fseeko64(fip,0,SEEK_SET);

    Long64_t readBytes = 0;

    UInt_t datapoint;
    UShort_t channelmask;
    UShort_t channelnum;
    UInt_t eventID;

    Double_t pretime = 0.0;
    Int_t preID = -1;
    Int_t Eventcounter = 0;

    while (readBytes < totalBytes){
        uint32_t Bheader[4];

        if (fread(&Bheader,sizeof(uint32_t),4,fip) < 4){
            cout << Eventcounter+1 << endl;
            cout << "Error : Failed to read file" << endl;
            Int_t errorc = ferror(fip);
            cout << "Error code : " << errorc << endl;
            cout << "Byte Progress : " << readBytes << "/" << totalBytes << endl;
            return -1;
        }
        readBytes += sizeof(uint32_t)*4;

        CAENheader head = CAENheader(Bheader[0],Bheader[1],Bheader[2],Bheader[3]);

        //head.OutputInfo();
        channelmask = head.GetChannelMask();
        channelnum = head.GetChannelNum();
        eventID = head.GetEventCount()-1;

        Double_t timetag_d = head.GetTimeTag();// double型のtimetag_dにdouble型のtimetagと8e-9の積を代入
        Double_t livetime = timetag_d - pretime;// double型のlivetimeにtimetag_dとpretimeの差を代入
        if(livetime<0)livetime+=4294967295* 8e-9 *0.5;// livetimeが0未満なら右辺の計算結果をlivetimeに足す
        pretime = timetag_d;// pretimeにtimetag_dを代入
        

        UInt_t wave_data_size = (head.GetEventSize()-4);

        UInt_t wave_data_per_ch = wave_data_size/channelnum;
        datapoint = wave_data_per_ch*2;

        EventID = eventID;
        Channels = channelnum;
        LiveT = livetime;
        UseChannel = channelmask;
        if (EventID-preID != 1){
            cout << "Skipped Events: Between " << EventID << " & " << preID << endl;
            if (InfoFlag)
            {
                infoos << "Skipped Events: Between " << EventID << " & " << preID << endl;
            }
        }
        preID = EventID;

        Int_t last_channel = -1;
        for (int j=0;j<channelnum;j++){
            UShort_t wave_data[wave_data_per_ch*2];
            uint32_t wave_data_bc[wave_data_per_ch];
            if (fread(&wave_data_bc,sizeof(uint32_t),wave_data_per_ch,fip) < 1){
                cout << "Error : Failed to read file" << endl;
                return -1;
            }
            readBytes += sizeof(uint32_t)*wave_data_per_ch;
            for (int k=last_channel+1;k<channelnum;k++){
                if(head.CheckChannelUse(k)){
                    last_channel = k;
                    break;
                }
            }

            TH1D *h1 = new TH1D(Form("h%d_Ch%d",eventID,last_channel),"",datapoint,0,datapoint*2);
            for (int k=0;k<wave_data_per_ch;k++){
                wave_data[k*2] = wave_data_bc[k]&0x3FFF;
                wave_data[k*2+1] = (wave_data_bc[k]>>16)&0x3FFF;
                h1 -> SetBinContent(k*2+1,wave_data[k*2]);
                h1 -> SetBinContent(k*2+2,wave_data[k*2+1]);
            }

            Double_t tPedestal = 0;
            for (Int_t k = 1; k <= Pedestalpoints; k++) {
                tPedestal += h1 -> GetBinContent(k);
            }
            tPedestal /= (Double_t)Pedestalpoints;

            // subtract pedestal
            for (Int_t k = 1; k <= datapoint; k++) {
                h1 -> SetBinContent(k, tPedestal - h1 -> GetBinContent(k));
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

            Pedestal[j]   = tPedestal;
            PedestalSigma[j] = tPedestalSigma;
            ADC[j]        = h1 -> Integral(h1 -> FindBin(tStartTime + ADCmin), h1 -> FindBin(tStartTime + ADCmax) - 1);
            PeakHeight[j] = h1 -> GetMaximum();
            PeakTime[j]   = h1 -> GetBinCenter(h1 -> GetMaximumBin());
            StartTime[j]  = tStartTime;
            StartBin[j]   = tStartBin;

            ADCtail[j]    = h1 -> Integral(h1 -> FindBin(tStartTime + ADCtailmin), h1 -> FindBin(tStartTime + ADCmax) - 1);

            #if RECORD_PULSE_SHAPE
            h1 -> Write();
            #endif

            delete h1;
        }
        Eventcounter++;
        tree_out -> Fill();
        if (Eventcounter%TerminalOutStep == 0){
            cout << "Converted " << Eventcounter << " Event" << endl;
            cout << "Byte Progress : " << readBytes << "/" << totalBytes << " = " << static_cast<Double_t>(readBytes)/static_cast<Double_t>(totalBytes)*100.0 << "%" << endl;
            if (InfoFlag)
            {
                infoos << "Converted " << Eventcounter << " Event" << endl;
            }
        }
    }
    cout << "Converted " << Eventcounter << " Event" << endl;
    if (InfoFlag)
    {
        infoos << "Converted " << Eventcounter << " Event" << endl;
    }
    if (EventID + 1 - Eventcounter > 0){
        cout << "Unrecording Events: " << EventID + 1 - Eventcounter << endl;
        if (InfoFlag)
        {
            infoos << "Unrecording Events: " << EventID + 1 - Eventcounter << endl;
        }
    }
    tree_out -> Write();
    fop -> Close();
    if (infoos.is_open())
    {
        infoos.close();
    }

    return 0;
}

Int_t extract_channel(TString infilename, Int_t targetChNo, TString outfilename = ""){
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
    UShort_t inChannels;
    UShort_t inUseChannel;
    Double_t inPedestal[16];
    Double_t inPedestalSigma[16];
    Double_t inADC[16];
    Double_t inADCtail[16];
    Double_t inPeakHeight[16];
    Double_t inPeakTime[16];
    Double_t inStartTime[16];
    Int_t    inStartBin[16];
    Double_t inLiveT;

    fip -> cd();

    TTree *tree_in = (TTree*)fip -> Get("tree");
    tree_in -> SetBranchAddress("EventID",       &inEventID);
    tree_in -> SetBranchAddress("Channels",      &inChannels);
    tree_in -> SetBranchAddress("UseChannel",    &inUseChannel);
    tree_in -> SetBranchAddress("Pedestal",      inPedestal);
    tree_in -> SetBranchAddress("PedestalSigma", inPedestalSigma);
    tree_in -> SetBranchAddress("ADC",           inADC);
    tree_in -> SetBranchAddress("ADCtail",       inADCtail);
    tree_in -> SetBranchAddress("PeakHeight",    inPeakHeight);
    tree_in -> SetBranchAddress("PeakTime",      inPeakTime);
    tree_in -> SetBranchAddress("StartTime",     inStartTime);
    tree_in -> SetBranchAddress("StartBin",      inStartBin);
    tree_in -> SetBranchAddress("LiveT",         &inLiveT);

    Int_t    EventID;
    UShort_t Channels;
    UShort_t UseChannel;
    Double_t Pedestal;
    Double_t PedestalSigma;
    Double_t ADC;
    Double_t ADCtail;
    Double_t PeakHeight;
    Double_t PeakTime;
    Double_t StartTime;
    Int_t    StartBin;
    Double_t LiveT;

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

    Int_t counter = 0;

    for (int i=0;i<tree_in->GetEntries();i++){
        fip -> cd();
        tree_in -> GetEntry(i);

        bool targetFlag = false;
        Int_t chid = -1;
        for (int j=0;j<16;j++){
            if ((inUseChannel>>j)&0x1 == 1){
                chid++;
                if (j == targetChNo){
                    targetFlag = true;
                    break;
                }
            }
        }
        if (!targetFlag){
            cout << "Passed EventID:" << inEventID << endl;
            continue;
        }

        EventID = inEventID;
        Pedestal = inPedestal[chid];
        PedestalSigma = inPedestalSigma[chid];
        ADC = inADC[chid];
        ADCtail = inADCtail[chid];
        PeakHeight = inPeakHeight[chid];
        PeakTime = inPeakTime[chid];
        StartTime = inStartTime[chid];
        StartBin = inStartBin[chid];
        LiveT = inLiveT;

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
        if(counter%TerminalOutStep == 0){
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
    const Double_t ADC_MAX  = 230; // ns, from StartTime
    const Double_t ADCtail_MIN  = 30; // ns, from StartTime

    const Int_t date = 20240720;
    const TString foldername = Form("../kaiseki/data/%d/",date);
    //const TString foldername = "../kaiseki/data/20240912/";
    const Int_t testNo = 0;
    const TString chara = "Cs";
    TString fin_name  = Form("%s%d_test%d_%s",foldername.Data(),date,testNo,chara.Data());
    TString fout_name = Form("%s/root_file/%.0lf_root_%d_test%d_%s_%.0lf.root",foldername.Data(),ADCtail_MIN,date,testNo,chara.Data(),ADC_MAX);
    TString info_name = Form("%s/info_file/%.0lf_root_%d_test%d_%s_%.0lf.txt",foldername.Data(),ADCtail_MIN,date,testNo,chara.Data(),ADC_MAX);
    cout << fout_name << endl;

    Convert_CAENbinROOT(fin_name,fout_name,info_name,ADC_MIN,ADC_MAX,ADCtail_MIN);
    cout << endl << "Finished convert binary to ROOT" << endl << endl;
    TString exfin_name = fout_name;
    //TString exfout_name = Form("%s%.0lf_root_%d_test%d_Cs_%.0lf_xxx.root",foldername.Data(),ADCtail_MIN,date,testNo,ADC_MAX);
    const Int_t targetchannel = 0;
    extract_channel(exfin_name,targetchannel);
    cout << endl << "Finished extract ch" << targetchannel << endl;
}


#ifndef AMC_RFM
#define AMC_RFM
#include "TF1.h"
#include "TFile.h"
#include "TString.h"

#include <iostream>

void WriteTF1(TFile *fp, TF1 *f)
{
    if (fp -> GetOption() == "READ")
    {
        std::cerr << "Error : File Write Error (@ROOT_File_Manipulation.cc)" << std::endl;
        std::cerr << "File( " << fp -> GetName() << " )'s Open Mode is Read Only." << std::endl;
    }
    else
    {
        fp -> cd();
        f -> Write();
        if (f -> Write() != 0)
        {
            std::cout << "Log : Function was Written @ " << fp -> GetName() << std::endl;
        }
    }
}
void WriteTF1(TString fpn, TF1 *f)
{
    TFile *fp = new TFile(fpn.Data(), "UPDATE");
    if (fp -> IsZombie())
    {
        std::cerr << "Error : File Open Error (@ROOT_File_Manipulation.cc)" << std::endl;
        std::cerr << "File( " << fpn << " ) open was failed." << std::endl;
    }
    else
    {
        f -> SetName("Calibration_Function");
        WriteTF1(fp, f);
        fp -> Close();
    }
}
void WriteCalibFunction(TFile *fp, TF1 *f)
{
    auto befName = f -> GetName();
    f -> SetName("Calibration_Function");
    WriteTF1(fp, f);
    f -> SetName(befName);
}
void WriteCalibFunction(TString fpn, TF1 *f)
{
    auto befName = f -> GetName();
    f -> SetName("Calibration_Function");
    WriteTF1(fpn, f);
    f -> SetName(befName);
}
TF1* ReadTF1(TFile *fp, TString fname)
{
    TF1 *f = (TF1*)fp->Get(fname.Data());
    if (f == NULL)
    {
        std::cerr << "Error : Read Function Error (@ROOT_File_Manipulation.cc)" << std::endl;
        std::cerr << "Get Function( " << fname << " ) was failed." << std::endl;
    }
    return f;
}
TString ReadCalibFunction(TFile *fp)
{
    TF1 *cf = ReadTF1(fp, TString("Calibration_Function"));
    TString out = "";
    if (cf == NULL)
    {
        std::cerr << "Error : Read Calibration Function Error (@ROOT_File_Manipulation.cc)" << std::endl;
        std::cerr << "Calibration Function Not Found (@ " << fp -> GetName() << " )" << std::endl;
        return out;
    }
    out = Form("%lf*ADC+%lf", cf->GetParameter(0), cf->GetParameter(1));
    std::cout << "Got Function : " << out << std::endl;
    return out;
}
TF1* ReadCalibFunctionTF1(TFile *fp)
{
    TF1 *cf = ReadTF1(fp, TString("Calibration_Function"));
    if (cf == NULL)
    {
        std::cerr << "Error : Read Calibration Function Error (@ROOT_File_Manipulation.cc)" << std::endl;
        std::cerr << "Calibration Function Not Found (@ " << fp -> GetName() << " )" << std::endl;
    }
    return cf;
}

#endif
// auto_fitter_summed.cpp
#include "TFile.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"
#include "TDirectory.h"
#include "TApplication.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>

using namespace std;

const int nCrystals = 16;

int silence_stdout() {
    fflush(stdout);
    int original_stdout = dup(fileno(stdout));
    int null_fd = open("/dev/null", O_WRONLY);
    dup2(null_fd, fileno(stdout));
    close(null_fd);
    return original_stdout;
}

void restore_stdout(int original_stdout) {
    fflush(stdout);
    dup2(original_stdout, fileno(stdout));
    close(original_stdout);
}

double get_max_y(const std::vector<TH1D*>& hists) {
    double maxy = 0;
    for (const auto& h : hists) {
        if (h) maxy = std::max(maxy, h->GetMaximum());
    }
    return maxy;
}

void calibrate_run(int run, vector<TH1D*>& raw, vector<TH1D*>& calibrated,
                   double calibration_params[nCrystals][2], TCanvas* cFit, TCanvas* cCalGraph) {
    TVirtualFitter::SetDefaultFitter("Minuit");
    if (gMinuit) gMinuit->SetPrintLevel(-1);
    double realE[2] = {1460, 2614};

    for (int c = 0; c < nCrystals; ++c) {
        cFit->cd(c+1);
        TH1D* h = raw[c];
        h->GetXaxis()->SetRangeUser(400, 3000);
        h->Draw();

        TSpectrum s1460(1), s2614(1);
        h->GetXaxis()->SetRangeUser(1360, 1560);
        s1460.Search(h, 1, "", 0.1);
        h->GetXaxis()->SetRangeUser(2500, 2800);
        s2614.Search(h, 1, "", 0.1);
        h->GetXaxis()->UnZoom();

        double *x1460 = s1460.GetPositionX();
        double *x2614 = s2614.GetPositionX();
        double peakADC[2] = {x1460[0], x2614[0]};

        for (int i = 0; i < 2; ++i) {
            TF1* fit = new TF1(Form("fit%d_c%d", i, c), "gaus(0)+pol1(3)", peakADC[i] - 10, peakADC[i] + 10);
            fit->FixParameter(0, h->GetBinContent(h->FindBin(peakADC[i])));
            fit->FixParameter(1, peakADC[i]);
            fit->FixParameter(2, 2);
            int fd = silence_stdout();
            h->Fit(fit, "MLISRNQ");
            restore_stdout(fd);
            peakADC[i] = fit->GetParameter(1);
            fit->SetLineColor(i+2);
            fit->Draw("same");
        }

        cCalGraph->cd(c+1);
        TGraph* gr = new TGraph(2, peakADC, realE);
        TF1* lin = new TF1(Form("lin%d", c), "pol1", peakADC[0], peakADC[1]);
        gr->Fit(lin, "MLISR");
        gr->SetMarkerStyle(20);
        gr->Draw("AP");

        calibration_params[c][0] = lin->GetParameter(0);
        calibration_params[c][1] = lin->GetParameter(1);
    }
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <start_run> <end_run>" << endl;
        return 1;
    }

    int runStart = atoi(argv[1]);
    int runEnd = atoi(argv[2]);
TFile* outFile = new TFile("../clover_analysis/cal_sum/summed_calibration_output.root", "RECREATE");
    TDirectory* dir_raw = outFile->mkdir("raw");
    TDirectory* dir_cal = outFile->mkdir("cal");

    vector<TH1D*> summed_raw(nCrystals);
    vector<TH1D*> summed_cal(nCrystals);
    double calib_params[nCrystals][2] = {};

    for (int i = 0; i < nCrystals; ++i) {
        summed_raw[i] = new TH1D(Form("sum_raw_%d", i), Form("Summed Raw %d", i), 6000, 0, 3000);
        summed_cal[i] = new TH1D(Form("sum_cal_%d", i), Form("Summed Calibrated %d", i), 6000, 0, 3000);
        summed_raw[i]->SetDirectory(0);
        summed_cal[i]->SetDirectory(0);
    }

    map<int, vector<TH1D*>> rawRuns;
    map<int, vector<TH1D*>> calRuns;

    TCanvas* cFits = new TCanvas("Fits", "Gaussian Fits", 1200, 900);
    cFits->Divide(4, 4);
    TCanvas* cGraphs = new TCanvas("Graphs", "Calibration Graphs", 1200, 900);
    cGraphs->Divide(4, 4);

    for (int run = runStart; run <= runEnd; ++run) {
      string fileName = Form("../../sorted_files/root_data_mvmelstrun%03d.bin_tree.root", run);
        TFile* inFile = new TFile(fileName.c_str(), "READ");
        if (!inFile || inFile->IsZombie()) continue;

        TTree* tree = (TTree*) inFile->Get("clover");
        if (!tree) continue;

        vector<TH1D*> raw(nCrystals), cal(nCrystals);
        for (int i = 0; i < nCrystals; ++i) {
            raw[i] = new TH1D(Form("raw_r%d_c%d", run, i), "", 6000, 0, 3000);
            cal[i] = new TH1D(Form("cal_r%d_c%d", run, i), "", 6000, 0, 3000);
            raw[i]->SetDirectory(0);
            cal[i]->SetDirectory(0);
        }

        TTreeReader reader(tree);
        TTreeReaderArray<Double_t> amp(reader, "clover_cross.amplitude");
        while (reader.Next()) {
            for (int i = 0; i < nCrystals; ++i) {
                double val = amp.At(i);
                if (val > 10) raw[i]->Fill(val);
            }
        }

        calibrate_run(run, raw, cal, calib_params, cFits, cGraphs);

        // Create canvas for raw per run
        TCanvas* cRawRun = new TCanvas(Form("cRaw_r%d", run), Form("Raw Run %d", run), 1200, 900);
        cRawRun->Divide(4, 4);
        for (int i = 0; i < nCrystals; ++i) {
            cRawRun->cd(i+1);
            raw[i]->Draw("hist");
        }

        // Create canvas for calibrated per run
        TCanvas* cCalRun = new TCanvas(Form("cCal_r%d", run), Form("Calibrated Run %d", run), 1200, 900);
        cCalRun->Divide(4, 4);
        for (int i = 0; i < nCrystals; ++i) {
            cCalRun->cd(i+1);
            cal[i]->Draw("hist");
        }

        for (int i = 0; i < nCrystals; ++i) {
            int nBins = raw[i]->GetNbinsX();
            for (int b = 1; b <= nBins; ++b) {
                double x = raw[i]->GetBinCenter(b);
                double y = raw[i]->GetBinContent(b);
                double E = calib_params[i][0] + calib_params[i][1] * x;
                cal[i]->Fill(E, y);
            }
            summed_raw[i]->Add(raw[i]);
            summed_cal[i]->Add(cal[i]);
        }

        rawRuns[run] = raw;
        calRuns[run] = cal;
    }

    outFile->cd("raw");
    for (auto& run_hvec : rawRuns) {
        for (int i = 0; i < nCrystals; ++i)
            run_hvec.second[i]->Write();
    }

    outFile->cd("cal");
    for (auto& run_hvec : calRuns) {
        for (int i = 0; i < nCrystals; ++i)
            run_hvec.second[i]->Write();
    }

    for (int i = 0; i < nCrystals; ++i) {
        summed_raw[i]->Write();
        summed_cal[i]->Write();
    }

    TCanvas* cSumRaw = new TCanvas("cSumRaw", "Summed Raw", 1200, 900);
    cSumRaw->Divide(4,4);
    TCanvas* cSumCal = new TCanvas("cSumCal", "Summed Calibrated", 1200, 900);
    cSumCal->Divide(4,4);
    TCanvas* cOverlay = new TCanvas("cOverlay", "Raw vs Calibrated", 1200, 900);
    cOverlay->Divide(4,4);

    double maxY = max(get_max_y(summed_raw), get_max_y(summed_cal));

    for (int i = 0; i < nCrystals; ++i) {
        cSumRaw->cd(i+1); summed_raw[i]->SetMaximum(maxY*1.1); summed_raw[i]->Draw("hist");
        cSumCal->cd(i+1); summed_cal[i]->SetMaximum(maxY*1.1); summed_cal[i]->Draw("hist");
        cOverlay->cd(i+1);
        summed_raw[i]->SetLineColor(kBlue);
        summed_cal[i]->SetLineColor(kRed);
        summed_raw[i]->Draw("hist");
        summed_cal[i]->Draw("hist same");
    }

    cFits->Draw();
    cGraphs->Draw();
    cSumRaw->Draw();
    cSumCal->Draw();
    cOverlay->Draw();

    outFile->Write();
    outFile->Close();

// Save calibration parameters to ../clover_analysis
ofstream paramOut(Form("../clover_analysis/cal_params_run%d_%d.txt", runStart, runEnd));
if (paramOut.is_open()) {
    paramOut << "# Crystal  Intercept  Slope\n";
    for (int i = 0; i < nCrystals; ++i) {
        paramOut << i << "  " << calib_params[i][0] << "  " << calib_params[i][1] << "\n";
    }
    paramOut.close();
    cout << "Calibration parameters saved to ../clover_analysis/calibration_parameters_run"
         << runStart << "_" << runEnd << ".txt" << endl;
} else {
    cerr << "Error: Could not write calibration parameters file!" << endl;
}

    cFits->Update();
    cGraphs->Update();
    cSumRaw->Update();
    cSumCal->Update();
    cOverlay->Update();

    app.Run();
    return 0;
}

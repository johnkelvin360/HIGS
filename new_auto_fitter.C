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
#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
using namespace std;

const int nPeaks     = 3;
const int nCrystals = 16;

// --------------------------------------------------
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
// --------------------------------------------------

void auto_fitter() {
    TVirtualFitter::SetDefaultFitter("Minuit");
    if (gMinuit) gMinuit->SetPrintLevel(-1);

    Double_t calibration_parameters[nCrystals][2];

    TCanvas *canvas  = new TCanvas("canvas",  "Fit Histogram with Gaussians", 800,600);
    canvas->Divide(4,4);

    TCanvas *canvas1 = new TCanvas("canvas1", "Calibration_Fit",              800,600);
    canvas1->Divide(4,4);

    TFile* output = new TFile("../../../jka197/clover_analysis/calibration_test_run026_linear.root","recreate");

    // Create ROOT directories
    TDirectory* dir_raw = output->mkdir("raw");
    TDirectory* dir_cal = output->mkdir("cal");

    TH1D *hclover_cross_amplitude[nCrystals];
    TH1D *hclover_cross_amplitude_calibrated[nCrystals];
    TH1D *hclover_cross_amplitude_full            = new TH1D("hclover_cross_amplitude_test_full",            "hclover_cross_amplitude_test_full",            6000,0,3000);
    TH1D *hclover_cross_amplitude_full_calibrated = new TH1D("hclover_cross_amplitude_test_full_calibrated", "hclover_cross_amplitude_test_full_calibrated",6000,0,3000);

    // Prevent deletion of full histograms after file is closed
    hclover_cross_amplitude_full->SetDirectory(0);
    hclover_cross_amplitude_full_calibrated->SetDirectory(0);

    for (int i = 0; i < nCrystals; ++i) {
        hclover_cross_amplitude[i] = new TH1D(
            Form("hclover_cross_amplitude_test_%d",i),
            Form("hclover_cross_amplitude_test_%d",i),
            6000,0,3000
        );
        hclover_cross_amplitude_calibrated[i] = new TH1D(
            Form("hclover_cross_amplitude_test_calibrated_%d",i),
            Form("hclover_cross_amplitude_test_calibrated_%d",i),
            6000,0,3000
        );

        // Prevent deletion of individual histograms after file is closed
        hclover_cross_amplitude[i]->SetDirectory(0);
        hclover_cross_amplitude_calibrated[i]->SetDirectory(0);
    }

    TFile* file1 = new TFile("../../root_data_mvmelstrun026.bin_tree.root","READ");
    TTree *data = dynamic_cast<TTree*>(file1->Get("clover"));
    if (!data) {
        cerr << "Tree 'clover' not found in the file!" << endl;
        file1->Close();
        return;
    }

    TTreeReader reader(data);
    TTreeReaderArray<Double_t> clover_cross_amplitude(reader,"clover_cross.amplitude");

    long int nEntries = data->GetEntries();
    int count = 0;
    while (reader.Next()) {
        for (int j = 0; j < nCrystals; ++j) {
            double Ej = clover_cross_amplitude.At(j);
            if (Ej > 10) {
                hclover_cross_amplitude[j]->Fill(Ej);
                hclover_cross_amplitude_full->Fill(Ej);
            }
        }
        ++count;
        int oneper = (int)(nEntries * 0.01);
        if (oneper > 0 && (count % oneper) == 0) {
            cout << "Events processed " << count
                 << " - percent done " << (count/oneper) << "%" << endl;
        }
    }

    // Write raw histograms
    output->cd("raw");
    for (int j = 0; j < nCrystals; ++j) {
        if (hclover_cross_amplitude[j]) hclover_cross_amplitude[j]->Write();
    }
    if (hclover_cross_amplitude_full) hclover_cross_amplitude_full->Write();

    // Calibration loop
    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        canvas->cd(crystal+1);
        TH1D* h = hclover_cross_amplitude[crystal];
        h->Draw();
        h->GetXaxis()->SetRangeUser(400,3000);

        TSpectrum *spectrum_511  = new TSpectrum(1);
        TSpectrum *spectrum_1460 = new TSpectrum(1);
        TSpectrum *spectrum_2640 = new TSpectrum(1);

        spectrum_511->Search(h,1,"",.1);
        h->GetXaxis()->SetRangeUser(1360,1560);
        spectrum_1460->Search(h,1,"",.1);
        h->GetXaxis()->SetRangeUser(2500,2800);
        spectrum_2640->Search(h,1,"",.1);

        Double_t *xPeaks_511  = spectrum_511->GetPositionX();
        Double_t *xPeaks_1460 = spectrum_1460->GetPositionX();
        Double_t *xPeaks_2640 = spectrum_2640->GetPositionX();

        h->GetXaxis()->UnZoom();
        std::vector<Double_t> xPeaks = { xPeaks_511[0], xPeaks_1460[0], xPeaks_2640[0] };

        double x[2];
        double y[2] = { 1460, 2614 };

        for (int i = 1; i < nPeaks; ++i) {
            TF1 *gauss  = new TF1(Form("gauss%d", i),  "gaus(0)+pol1(3)", xPeaks[i] - 10, xPeaks[i] + 10);
            gauss->FixParameter(0, h->GetBinContent(h->FindBin(xPeaks[i])));
            gauss->FixParameter(1, xPeaks[i]);
            gauss->FixParameter(2, 2);

            int stdout_fd = silence_stdout();
            h->Fit(gauss, "MLISRNQ");
            restore_stdout(stdout_fd);

            TF1 *gauss_final = new TF1(Form("gauss_final%d", i), "gaus(0)+pol1(3)", xPeaks[i] - 10, xPeaks[i] + 10);
            gauss_final->SetParameters(gauss->GetParameters());

            stdout_fd = silence_stdout();
            h->Fit(gauss_final, "MLISRNQ");
            restore_stdout(stdout_fd);

            x[i-1] = gauss_final->GetParameter(1);
            gauss_final->SetLineColor(i + 2);
            gauss_final->Draw("same");

            delete gauss;
        }

        cout << "Crystal #" << crystal
             << "  1460 keV centroid @ channel = " << x[0]
             << " ; 2614 keV centroid @ channel = " << x[1] << endl;

        canvas->Update();
        gPad->SetLogy();

        canvas1->cd(crystal+1);
        TGraph* graph = new TGraph(2, x, y);
        graph->SetName(Form("Calibration_of_Crystal_#_%d",crystal));
        graph->SetTitle(Form("Calibration_of_Crystal_#_%d",crystal));

        TF1 *linear = new TF1(Form("linear%d", crystal), "pol1", x[0], x[1]);
        graph->Fit(linear, "MLISR");

        calibration_parameters[crystal][0] = linear->GetParameter(0);
        calibration_parameters[crystal][1] = linear->GetParameter(1);

        graph->SetMarkerSize(1.2);
        graph->SetMarkerStyle(20);
        graph->Draw("AP");
        canvas1->Update();

        delete linear;
    }

    // Write calibrated histograms
    output->cd("cal");
    for (int j = 0; j < nCrystals; ++j) {
        if (hclover_cross_amplitude_calibrated[j]) hclover_cross_amplitude_calibrated[j]->Write();
    }
    if (hclover_cross_amplitude_full_calibrated) hclover_cross_amplitude_full_calibrated->Write();

    // Save parameters
    ofstream calib_file("../../../jka197/clover_analysis/clover_calibration_parameters.txt");
    if (calib_file.is_open()) {
        calib_file << "# CrystalIndex  Slope(a)  Intercept(b)\n";
        for (int i = 0; i < nCrystals; ++i) {
            calib_file << i << " "
                       << calibration_parameters[i][1] << " "
                       << calibration_parameters[i][0] << "\n";
        }
        calib_file.close();
        cout << "Calibration parameters saved.\n";
    } else {
        cerr << "Could not write calibration parameters.\n";
    }

    output->cd();
    output->Write();
    output->Close();

    // Display canvases
    canvas->Draw();
    canvas1->Draw();
}

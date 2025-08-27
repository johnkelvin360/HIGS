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
#include "TLegend.h"
#include "TColor.h"
#include "TApplication.h"
#include "TPaveText.h"
#include "TMarker.h"
#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <fstream>
using namespace std;

const int nPeaks = 3;
const int nCrystals = 16;

// --------------------------------------------------
int silence_stdout()
{
    fflush(stdout);
    int original_stdout = dup(fileno(stdout));
    int null_fd = open("/dev/null", O_WRONLY);
    dup2(null_fd, fileno(stdout));
    close(null_fd);
    return original_stdout;
}

void restore_stdout(int original_stdout)
{
    fflush(stdout);
    dup2(original_stdout, fileno(stdout));
    close(original_stdout);
}

// --------------------------------------------------
void run_by_run_calibration(int start_run, int end_run)
{
    TVirtualFitter::SetDefaultFitter("Minuit");
    if (gMinuit)
        gMinuit->SetPrintLevel(-1);

    int n_runs = end_run - start_run + 1;

    // Color palette for different runs
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow, kViolet,
                    kSpring, kTeal, kAzure, kPink, kGray, kBlack, kOrange + 7, kSpring + 5};

    // Storage for calibration parameters and histograms
    vector<map<int, pair<double, double>>> calibration_params(n_runs); // [run_index][crystal] -> (slope, intercept)
    vector<vector<TH1D *>> raw_histograms(n_runs);                     // [run_index][crystal]
    vector<vector<TH1D *>> cal_histograms(n_runs);                     // [run_index][crystal]
    vector<TH1D *> sum_raw_histograms(nCrystals);                      // Sum across runs per crystal
    vector<TH1D *> sum_cal_histograms(nCrystals);                      // Sum across runs per crystal

    // Storage for raw amplitudes (alternative approach to avoid re-reading)
    vector<vector<vector<double>>> raw_amplitudes(n_runs); // [run_index][crystal][event]

    // Canvases to keep open
    vector<TCanvas *> individual_raw_canvases;
    vector<TCanvas *> individual_cal_canvases;
    vector<TCanvas *> overlay_raw_canvases;
    vector<TCanvas *> overlay_cal_canvases;
    vector<TCanvas *> sum_canvases;
    vector<TCanvas *> linear_fit_canvases;

    // Initialize sum histograms
    for (int crystal = 0; crystal < nCrystals; ++crystal)
    {
        sum_raw_histograms[crystal] = new TH1D(
            Form("sum_raw_crystal_%d", crystal),
            Form("Sum Raw Energy - Crystal %d;Channel;Counts", crystal),
            6000, 0, 3000);
        sum_cal_histograms[crystal] = new TH1D(
            Form("sum_cal_crystal_%d", crystal),
            Form("Sum Calibrated Energy - Crystal %d;Energy (keV);Counts", crystal),
            6000, 0, 3000);
    }

    // Initialize raw_amplitudes storage
    for (int run_idx = 0; run_idx < n_runs; ++run_idx)
    {
        raw_amplitudes[run_idx].resize(nCrystals);
    }

    // Process each run
    for (int run_idx = 0; run_idx < n_runs; ++run_idx)
    {
        int run_number = start_run + run_idx;
        cout << "==========================================" << endl;
        cout << "Processing run " << run_number << "..." << endl;
        cout << "==========================================" << endl;

        // Initialize histograms for this run
        raw_histograms[run_idx].resize(nCrystals);
        cal_histograms[run_idx].resize(nCrystals);

        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            raw_histograms[run_idx][crystal] = new TH1D(
                Form("raw_run%d_crystal_%d", run_number, crystal),
                Form("Raw Energy - Run %d Crystal %d;Channel;Counts", run_number, crystal),
                6000, 0, 3000);
            cal_histograms[run_idx][crystal] = new TH1D(
                Form("cal_run%d_crystal_%d", run_number, crystal),
                Form("Calibrated Energy - Run %d Crystal %d;Energy (keV);Counts", run_number, crystal),
                6000, 0, 3000);
        }

        // Read data for this run - STEP 1: Fill raw histograms and store amplitudes
        TString input_filename = Form("../../sorted_files/root_data_mvmelstrun%03d.bin_tree.root", run_number);
        cout << "Opening file: " << input_filename << endl;

        TFile *file1 = new TFile(input_filename, "READ");
        if (!file1 || file1->IsZombie())
        {
            cerr << "ERROR: Could not open file for run " << run_number << endl;
            continue;
        }

        TTree *data = dynamic_cast<TTree *>(file1->Get("clover"));
        if (!data)
        {
            cerr << "ERROR: Tree 'clover' not found in run " << run_number << "!" << endl;
            file1->Close();
            continue;
        }

        TTreeReader reader(data);
        TTreeReaderArray<Double_t> clover_cross_amplitude(reader, "clover_cross.amplitude");

        long int nEntries = data->GetEntries();
        cout << "Found " << nEntries << " entries in run " << run_number << endl;

        int count = 0;
        int oneper = (int)(nEntries * 0.01);
        if (oneper == 0)
            oneper = 1;

        cout << "Reading data and storing amplitudes..." << endl;
        while (reader.Next())
        {
            for (int crystal = 0; crystal < nCrystals; ++crystal)
            {
                double Ej = clover_cross_amplitude.At(crystal);
                if (Ej > 10)
                {
                    // STEP 1: Fill raw histograms and store amplitudes for later calibration
                    raw_histograms[run_idx][crystal]->Fill(Ej);
                    raw_amplitudes[run_idx][crystal].push_back(Ej);
                }
            }

            ++count;

            if (count % oneper == 0)
            {
                // Use floating point to avoid integer overflow
                double percent_done = (static_cast<double>(count) * 100.0) / static_cast<double>(nEntries);
                cout << "Processed " << count << " events (" << percent_done << "%)" << endl;
            }
        }

        cout << "Finished reading run " << run_number << ". Processed " << count << " events." << endl;
        file1->Close();

        // STEP 2: Fit peaks in raw histograms to get calibration parameters
        cout << "Starting calibration for run " << run_number << "..." << endl;
        TCanvas *calib_canvas = new TCanvas(
            Form("calib_run%d", run_number),
            Form("Calibration - Run %d", run_number),
            800, 600);
        calib_canvas->Divide(4, 4);

        // Create 4x4 canvas for linear fit plots
        TCanvas *linear_fit_canvas = new TCanvas(
            Form("linear_fit_run%d", run_number),
            Form("Linear Fits - Run %d", run_number),
            800, 600);
        linear_fit_canvas->Divide(4, 4);

        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            cout << "Calibrating crystal " << crystal << "..." << endl;
            calib_canvas->cd(crystal + 1);
            TH1D *h = raw_histograms[run_idx][crystal];
            h->Draw();
            h->GetXaxis()->SetRangeUser(400, 3000);

            TSpectrum *spectrum_511 = new TSpectrum(1);
            TSpectrum *spectrum_1460 = new TSpectrum(1);
            TSpectrum *spectrum_2640 = new TSpectrum(1);

            spectrum_511->Search(h, 1, "", .1);
            h->GetXaxis()->SetRangeUser(1360, 1560);
            spectrum_1460->Search(h, 1, "", .1);
            h->GetXaxis()->SetRangeUser(2500, 2800);
            spectrum_2640->Search(h, 1, "", .1);

            Double_t *xPeaks_511 = spectrum_511->GetPositionX();
            Double_t *xPeaks_1460 = spectrum_1460->GetPositionX();
            Double_t *xPeaks_2640 = spectrum_2640->GetPositionX();

            h->GetXaxis()->UnZoom();
            std::vector<Double_t> xPeaks = {xPeaks_511[0], xPeaks_1460[0], xPeaks_2640[0]};

            double x[2];
            double y[2] = {1460, 2614};

            for (int i = 1; i < nPeaks; ++i)
            {
                TF1 *gauss = new TF1(Form("gauss%d_%d", crystal, i), "gaus(0)+pol1(3)", xPeaks[i] - 10, xPeaks[i] + 10);
                gauss->FixParameter(0, h->GetBinContent(h->FindBin(xPeaks[i])));
                gauss->FixParameter(1, xPeaks[i]);
                gauss->FixParameter(2, 2);

                int stdout_fd = silence_stdout();
                h->Fit(gauss, "MLISRNQ");
                restore_stdout(stdout_fd);

                TF1 *gauss_final = new TF1(Form("gauss_final%d_%d", crystal, i), "gaus(0)+pol1(3)", xPeaks[i] - 10, xPeaks[i] + 10);
                gauss_final->SetParameters(gauss->GetParameters());

                stdout_fd = silence_stdout();
                h->Fit(gauss_final, "MLISRNQ");
                restore_stdout(stdout_fd);

                x[i - 1] = gauss_final->GetParameter(1);
                gauss_final->SetLineColor(i + 2);
                gauss_final->Draw("same");

                delete gauss;
            }

            cout << "Crystal #" << crystal
                 << "  1460 keV centroid @ channel = " << x[0]
                 << " ; 2614 keV centroid @ channel = " << x[1] << endl;

            // Create and fit the linear calibration graph
            TGraph *graph = new TGraph(2, x, y);
            graph->SetName(Form("Calibration_Crystal_%d", crystal));
            graph->SetTitle(Form("Calibration Crystal %d;Channel;Energy (keV)", crystal));
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(1.2);
            graph->SetMarkerColor(kRed);

            TF1 *linear = new TF1(Form("linear_%d_%d", run_idx, crystal), "pol1", x[0] - 100, x[1] + 100);

            int stdout_fd = silence_stdout();
            graph->Fit(linear, "MLISR");
            restore_stdout(stdout_fd);

            calibration_params[run_idx][crystal] = make_pair(
                linear->GetParameter(1), // slope
                linear->GetParameter(0)  // intercept
            );

            cout << "Crystal " << crystal << ": slope = " << linear->GetParameter(1)
                 << ", intercept = " << linear->GetParameter(0) << endl;

            // Draw linear fit on the dedicated 4x4 canvas
            linear_fit_canvas->cd(crystal + 1);
            graph->Draw("AP");
            linear->SetLineColor(kBlue);
            linear->Draw("same");

            // Add text with calibration parameters
            TPaveText *pt = new TPaveText(0.15, 0.7, 0.55, 0.85, "NDC");
            pt->SetFillColor(0);
            pt->SetBorderSize(1);
            pt->SetTextSize(0.06);
            pt->AddText(Form("Slope: %.4f", linear->GetParameter(1)));
            pt->AddText(Form("Intercept: %.1f", linear->GetParameter(0)));
            pt->Draw();

            // Add marker for 511 keV peak if found
            if (xPeaks_511[0] > 0)
            {
                double energy_511 = linear->GetParameter(1) * xPeaks_511[0] + linear->GetParameter(0);
                TMarker *m511 = new TMarker(xPeaks_511[0], energy_511, 29);
                m511->SetMarkerColor(kGreen);
                m511->SetMarkerSize(1.5);
                m511->Draw();
            }

            gPad->SetGrid(1, 1);
        }

        // Save linear fit canvas
        linear_fit_canvas->Draw();
        linear_fit_canvases.push_back(linear_fit_canvas);

        // STEP 3: Apply calibration using stored raw amplitudes (avoiding bin artifacts)
        cout << "Applying calibration to stored amplitudes for run " << run_number << "..." << endl;

        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            double slope = calibration_params[run_idx][crystal].first;
            double intercept = calibration_params[run_idx][crystal].second;

            int crystal_events = raw_amplitudes[run_idx][crystal].size();
            cout << "Calibrating crystal " << crystal << " with " << crystal_events << " events..." << endl;

            // Apply calibration to each stored amplitude
            for (double amplitude : raw_amplitudes[run_idx][crystal])
            {
                double energy = slope * amplitude + intercept;
                cal_histograms[run_idx][crystal]->Fill(energy);

                // Add to sum histograms
                sum_raw_histograms[crystal]->Fill(amplitude);
                sum_cal_histograms[crystal]->Fill(energy);
            }
        }

        cout << "Finished calibrating run " << run_number << endl;

        // Clear amplitudes to save memory for next run
        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            raw_amplitudes[run_idx][crystal].clear();
            raw_amplitudes[run_idx][crystal].shrink_to_fit();
        }

        // Save individual run to ROOT file (including linear fit canvas)
        TString outname = Form("../clover_analysis/cal_sum/calibration_run%03d_linear.root", run_number);
        cout << "Saving results to: " << outname << endl;

        TFile *output = new TFile(outname, "RECREATE");
        TDirectory *dirRaw = output->mkdir("raw");
        TDirectory *dirCal = output->mkdir("cal");
        TDirectory *dirFits = output->mkdir("linear_fits");

        // Write raw histograms
        dirRaw->cd();
        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            raw_histograms[run_idx][crystal]->Write();
        }

        // Write calibrated histograms
        dirCal->cd();
        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            cal_histograms[run_idx][crystal]->Write();
        }

        // Write linear fit canvas
        dirFits->cd();
        linear_fit_canvas->Write();

        output->Close();
        cout << "Saved results for run " << run_number << endl;

        std::ofstream param_file(Form("../clover_analysis/calibration_run%d.txt", run_number));
        for (int c = 0; c < nCrystals; c++)
        {
            param_file << "Crystal " << c
                       << " slope = " << calibration_params[run_idx][c].first
                       << " intercept = " << calibration_params[run_idx][c].second
                       << std::endl;
        }
        param_file.close();

        std::cout << "Finished processing run " << run_number << std::endl;
        delete calib_canvas;
    }

    // Save linear fit canvases to separate files
    for (int run_idx = 0; run_idx < n_runs; ++run_idx)
    {
        int run_number = start_run + run_idx;
        TString linear_fit_outname = Form("../clover_analysis/cal_sum/linear_fits_run%03d.root", run_number);
        TFile *linear_fit_output = new TFile(linear_fit_outname, "RECREATE");
        linear_fit_canvases[run_idx]->Write();
        linear_fit_output->Close();
        cout << "Saved linear fit plots to: " << linear_fit_outname << endl;
    }

    // Create all requested plots in 4x4 grid format with 800x600 canvases and legends

    // 1. Individual raw energy plots per run per crystal (4x4 grid)
    cout << "Creating individual raw energy plots..." << endl;
    for (int run_idx = 0; run_idx < n_runs; ++run_idx)
    {
        int run_number = start_run + run_idx;

        TCanvas *cRaw = new TCanvas(Form("Run%d_Raw", run_number), Form("Run %d Raw", run_number), 800, 600);
        cRaw->Divide(4, 4);

        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            cRaw->cd(crystal + 1);
            TH1D *hraw = raw_histograms[run_idx][crystal];
            hraw->SetLineColor(kRed);
            hraw->SetLineWidth(1);
            hraw->SetFillStyle(0);
            hraw->SetMarkerStyle(0);
            hraw->Draw("HIST");
            gPad->SetLogy();
        }

        // Save to ROOT file
        TString outname = Form("../clover_analysis/cal_sum/calibration_run%03d_linear.root", run_number);
        TFile *output = new TFile(outname, "UPDATE");
        cRaw->Write();
        output->Close();

        individual_raw_canvases.push_back(cRaw);
        cRaw->Draw();
    }

    // 2. Individual calibrated energy plots per run per crystal (4x4 grid)
    cout << "Creating individual calibrated energy plots..." << endl;
    for (int run_idx = 0; run_idx < n_runs; ++run_idx)
    {
        int run_number = start_run + run_idx;

        TCanvas *cCal = new TCanvas(Form("Run%d_Cal", run_number), Form("Run %d Calibrated", run_number), 800, 600);
        cCal->Divide(4, 4);

        for (int crystal = 0; crystal < nCrystals; ++crystal)
        {
            cCal->cd(crystal + 1);
            TH1D *hcal = cal_histograms[run_idx][crystal];
            hcal->SetLineColor(kBlue);
            hcal->SetLineWidth(1);
            hcal->SetFillStyle(0);
            hcal->SetMarkerStyle(0);
            hcal->Draw("HIST");
            gPad->SetLogy();
        }

        // Save to ROOT file
        TString outname = Form("../clover_analysis/cal_sum/calibration_run%03d_linear.root", run_number);
        TFile *output = new TFile(outname, "UPDATE");
        cCal->Write();
        output->Close();

        individual_cal_canvases.push_back(cCal);
        cCal->Draw();
    }

    // 3. Overlay raw energy across runs (4x4 grid)
    cout << "Creating overlay raw plots..." << endl;
    TCanvas *cOR = new TCanvas("OverlayRawRuns", "Overlay Raw Runs", 800, 600);
    cOR->Divide(4, 4);

    for (int crystal = 0; crystal < nCrystals; ++crystal)
    {
        cOR->cd(crystal + 1);

        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        bool first = true;

        for (int run_idx = 0; run_idx < n_runs; ++run_idx)
        {
            int run_number = start_run + run_idx;
            TH1D *hraw = raw_histograms[run_idx][crystal];
            hraw->SetFillStyle(0);
            hraw->SetLineWidth(2);
            hraw->SetLineColor(colors[run_idx % 16]);

            if (first)
            {
                hraw->Draw("HIST");
                first = false;
            }
            else
            {
                hraw->Draw("HIST SAME");
            }

            legend->AddEntry(hraw, Form("Run %d", run_number), "l");
        }

        legend->Draw();
        gPad->SetLogy();
    }

    // Save overlay to sum file
    TString sum_outname = Form("../clover_analysis/cal_sum/sum_calibration_runs_%03d_to_%03d.root", start_run, end_run);
    TFile *sum_output = new TFile(sum_outname, "UPDATE");
    cOR->Write();
    sum_output->Close();

    overlay_raw_canvases.push_back(cOR);
    cOR->Draw();

    // 4. Overlay calibrated energy across runs (4x4 grid)
    cout << "Creating overlay calibrated plots..." << endl;
    TCanvas *cOC = new TCanvas("OverlayCalRuns", "Overlay Calibrated Runs", 800, 600);
    cOC->Divide(4, 4);

    for (int crystal = 0; crystal < nCrystals; ++crystal)
    {
        cOC->cd(crystal + 1);

        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        bool first = true;

        for (int run_idx = 0; run_idx < n_runs; ++run_idx)
        {
            int run_number = start_run + run_idx;
            TH1D *hcal = cal_histograms[run_idx][crystal];
            hcal->SetFillStyle(0);
            hcal->SetLineWidth(2);
            hcal->SetLineColor(colors[run_idx % 16]);

            if (first)
            {
                hcal->Draw("HIST");
                first = false;
            }
            else
            {
                hcal->Draw("HIST SAME");
            }

            legend->AddEntry(hcal, Form("Run %d", run_number), "l");
        }

        legend->Draw();
        gPad->SetLogy();
    }

    // Save overlay to sum file
    sum_output = new TFile(sum_outname, "UPDATE");
    cOC->Write();
    sum_output->Close();

    overlay_cal_canvases.push_back(cOC);
    cOC->Draw();

    // 5. Summed raw energy per crystal (4x4 grid)
    cout << "Creating summed raw plots..." << endl;
    TCanvas *cSR = new TCanvas("SumRaw", "Summed Raw Energies", 800, 600);
    cSR->Divide(4, 4);

    for (int crystal = 0; crystal < nCrystals; ++crystal)
    {
        cSR->cd(crystal + 1);
        sum_raw_histograms[crystal]->SetLineColor(kRed);
        sum_raw_histograms[crystal]->SetLineWidth(2);
        sum_raw_histograms[crystal]->SetFillStyle(0);
        sum_raw_histograms[crystal]->Draw("HIST");
        gPad->SetLogy();
    }

    // Save to sum file
    sum_output = new TFile(sum_outname, "UPDATE");
    cSR->Write();
    sum_output->Close();

    sum_canvases.push_back(cSR);
    cSR->Draw();

    // 6. Summed calibrated energy per crystal (4x4 grid)
    cout << "Creating summed calibrated plots..." << endl;
    TCanvas *cSC = new TCanvas("SumCal", "Summed Calibrated Energies", 800, 600);
    cSC->Divide(4, 4);

    for (int crystal = 0; crystal < nCrystals; ++crystal)
    {
        cSC->cd(crystal + 1);
        sum_cal_histograms[crystal]->SetLineColor(kBlue);
        sum_cal_histograms[crystal]->SetLineWidth(2);
        sum_cal_histograms[crystal]->SetFillStyle(0);
        sum_cal_histograms[crystal]->Draw("HIST");
        gPad->SetLogy();
    }

    // Save to sum file
    sum_output = new TFile(sum_outname, "UPDATE");
    cSC->Write();
    sum_output->Close();

    sum_canvases.push_back(cSC);
    cSC->Draw();

    // 7. Overlay summed raw vs calibrated per crystal (4x4 grid)
    cout << "Creating summed overlay plots..." << endl;
    TCanvas *cSO = new TCanvas("SumOverlay", "Overlay Summed Raw vs Calibrated", 800, 600);
    cSO->Divide(4, 4);

    for (int crystal = 0; crystal < nCrystals; ++crystal)
    {
        cSO->cd(crystal + 1);

        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

        // Set X-axis title to Energy (keV) for both histograms
        sum_raw_histograms[crystal]->SetTitle(Form("Sum Raw vs Calibrated - Crystal %d;Energy (keV);Counts", crystal));
        sum_cal_histograms[crystal]->SetTitle(Form("Sum Raw vs Calibrated - Crystal %d;Energy (keV);Counts", crystal));

        // Summed raw in red
        sum_raw_histograms[crystal]->SetFillStyle(0);
        sum_raw_histograms[crystal]->SetLineWidth(2);
        sum_raw_histograms[crystal]->SetLineColor(kRed);
        sum_raw_histograms[crystal]->Draw("HIST");
        legend->AddEntry(sum_raw_histograms[crystal], "Sum Raw", "l");

        // Summed calibrated in blue
        sum_cal_histograms[crystal]->SetFillStyle(0);
        sum_cal_histograms[crystal]->SetLineWidth(2);
        sum_cal_histograms[crystal]->SetLineColor(kBlue);
        sum_cal_histograms[crystal]->Draw("HIST SAME");
        legend->AddEntry(sum_cal_histograms[crystal], "Sum Calibrated", "l");

        legend->Draw();
        gPad->SetLogy();
    }

    // Save to sum file
    sum_output = new TFile(sum_outname, "UPDATE");
    cSO->Write();
    sum_output->Close();

    sum_canvases.push_back(cSO);
    cSO->Draw();

    // Save calibration parameters
    ofstream calib_file("../clover_analysis/cal_sum/clover_calibration_parameters.txt");
    if (calib_file.is_open())
    {
        calib_file << "# Run CrystalIndex Slope(a) Intercept(b)\n";
        for (int run_idx = 0; run_idx < n_runs; ++run_idx)
        {
            int run_number = start_run + run_idx;
            for (int crystal = 0; crystal < nCrystals; ++crystal)
            {
                calib_file << run_number << " " << crystal << " "
                           << calibration_params[run_idx][crystal].first << " "
                           << calibration_params[run_idx][crystal].second << "\n";
            }
        }
        calib_file.close();
        cout << "Calibration parameters saved to clover_calibration_parameters_multirun.txt" << endl;
    }
    else
    {
        cerr << "Could not write calibration parameters.\n";
    }

    cout << "==========================================" << endl;
    cout << "Processing complete! All canvases are displayed." << endl;
    cout << "Individual run files saved in: ../clover_analysis/cal_sum/" << endl;
    cout << "Sum results saved in: " << sum_outname << endl;
    cout << "Press Ctrl+C to exit and close all canvases." << endl;
    cout << "==========================================" << endl;
}

// Main function
int main(int argc, char *argv[])
{
    // Parse command line arguments
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " <start_run> <end_run>" << endl;
        cout << "Example: " << argv[0] << " 98 100" << endl;
        return 1;
    }

    int start_run = atoi(argv[1]);
    int end_run = atoi(argv[2]);

    if (start_run > end_run)
    {
        cout << "Error: start_run must be <= end_run" << endl;
        return 1;
    }

    // keep canvases open
    TApplication app("App", &argc, argv);

    run_by_run_calibration(start_run, end_run);

    app.Run();

    return 0;
}

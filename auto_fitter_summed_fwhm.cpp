// auto_fitter_summed.cpp
#include "TFile.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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
#include "TLine.h"
#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <fcntl.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;

const int nPeaks = 3;
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

// Function to calculate FWHM of a peak in a histogram (fallback method)
double calculateFWHM(TH1D* hist, double peak_center, double search_range = 50.0) {
    // Find the maximum bin around the peak
    int bin_center = hist->FindBin(peak_center);
    int bin_low = hist->FindBin(peak_center - search_range);
    int bin_high = hist->FindBin(peak_center + search_range);
    
    double max_value = hist->GetBinContent(bin_center);
    if (max_value <= 0) return 0.0;
    double half_max = max_value / 2.0;
    
    // Find left half-max point
    double left_edge = 0.0;
    for (int bin = bin_center; bin >= bin_low; bin--) {
        if (hist->GetBinContent(bin) <= half_max) {
            // Linear interpolation for better precision
            double x1 = hist->GetBinCenter(bin);
            double y1 = hist->GetBinContent(bin);
            double x2 = hist->GetBinCenter(bin + 1);
            double y2 = hist->GetBinContent(bin + 1);
            
            if (y2 != y1) {
                left_edge = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1);
            } else {
                left_edge = x1;
            }
            break;
        }
    }
    
    // Find right half-max point
    double right_edge = 0.0;
    for (int bin = bin_center; bin <= bin_high; bin++) {
        if (hist->GetBinContent(bin) <= half_max) {
            // Linear interpolation for better precision
            double x1 = hist->GetBinCenter(bin - 1);
            double y1 = hist->GetBinContent(bin - 1);
            double x2 = hist->GetBinCenter(bin);
            double y2 = hist->GetBinContent(bin);
            
            if (y2 != y1) {
                right_edge = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1);
            } else {
                right_edge = x2;
            }
            break;
        }
    }
    
    if (left_edge > 0 && right_edge > 0) {
        return right_edge - left_edge;
    }
    
    return 0.0; // Return 0 if FWHM couldn't be calculated
}

// --------------------------------------------------
void run_by_run_calibration(int start_run, int end_run) {
    TVirtualFitter::SetDefaultFitter("Minuit");
    if (gMinuit) gMinuit->SetPrintLevel(-1);
    
    int n_runs = end_run - start_run + 1;
    
    // Color palette for different runs / overlays
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow, kViolet, 
                   kSpring, kTeal, kAzure, kPink, kGray, kBlack, kOrange+7, kSpring+5};
    
    // Storage for calibration parameters and histograms
    vector<map<int, pair<double, double>>> calibration_params(n_runs); // [run_index][crystal] -> (slope, intercept)
    vector<vector<TH1D*>> raw_histograms(n_runs); // [run_index][crystal]
    vector<vector<TH1D*>> cal_histograms(n_runs); // [run_index][crystal]
    vector<TH1D*> sum_raw_histograms(nCrystals); // Sum across runs per crystal
    vector<TH1D*> sum_cal_histograms(nCrystals); // Sum across runs per crystal
    
    // Storage for raw amplitudes (alternative approach to avoid re-reading)
    vector<vector<vector<double>>> raw_amplitudes(n_runs); // [run_index][crystal][event]
    
    // Storage for FWHM values and errors for individual runs
    vector<vector<double>> fwhm_511_individual(n_runs, vector<double>(nCrystals, 0.0));
    vector<vector<double>> fwhm_511_err_individual(n_runs, vector<double>(nCrystals, 0.0));
    vector<vector<double>> fwhm_1460_individual(n_runs, vector<double>(nCrystals, 0.0));
    vector<vector<double>> fwhm_1460_err_individual(n_runs, vector<double>(nCrystals, 0.0));
    vector<vector<double>> fwhm_2614_individual(n_runs, vector<double>(nCrystals, 0.0));
    vector<vector<double>> fwhm_2614_err_individual(n_runs, vector<double>(nCrystals, 0.0));
    
    // Storage for summed FWHM and errors
    vector<double> fwhm_511_sum(nCrystals, 0.0);
    vector<double> fwhm_511_sum_err(nCrystals, 0.0);
    vector<double> fwhm_1460_sum(nCrystals, 0.0);
    vector<double> fwhm_1460_sum_err(nCrystals, 0.0);
    vector<double> fwhm_2614_sum(nCrystals, 0.0);
    vector<double> fwhm_2614_sum_err(nCrystals, 0.0);
    
    // Canvases to keep open
    vector<TCanvas*> individual_raw_canvases;
    vector<TCanvas*> individual_cal_canvases;
    vector<TCanvas*> overlay_raw_canvases;
    vector<TCanvas*> overlay_cal_canvases;
    vector<TCanvas*> sum_canvases;
    vector<TCanvas*> linear_fit_canvases;
    
    // Initialize sum histograms
    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        sum_raw_histograms[crystal] = new TH1D(
            Form("sum_raw_crystal_%d", crystal),
            Form("Sum Raw Energy - Crystal %d;Channel;Counts", crystal),
            6000, 0, 3000
        );
        sum_cal_histograms[crystal] = new TH1D(
            Form("sum_cal_crystal_%d", crystal),
            Form("Sum Calibrated Energy - Crystal %d;Energy (keV);Counts", crystal),
            6000, 0, 3000
        );
    }
    
    // Initialize raw_amplitudes storage
    for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
        raw_amplitudes[run_idx].resize(nCrystals);
    }
    
    // Process each run
    for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
        int run_number = start_run + run_idx;
        cout << "==========================================" << endl;
        cout << "Processing run " << run_number << "..." << endl;
        cout << "==========================================" << endl;
        
        // Initialize histograms for this run
        raw_histograms[run_idx].resize(nCrystals);
        cal_histograms[run_idx].resize(nCrystals);
        
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            raw_histograms[run_idx][crystal] = new TH1D(
                Form("raw_run%d_crystal_%d", run_number, crystal),
                Form("Raw Energy - Run %d Crystal %d;Channel;Counts", run_number, crystal),
                6000, 0, 3000
            );
            cal_histograms[run_idx][crystal] = new TH1D(
                Form("cal_run%d_crystal_%d", run_number, crystal),
                Form("Calibrated Energy - Run %d Crystal %d;Energy (keV);Counts", run_number, crystal),
                6000, 0, 3000
            );
        }
        
        // Read data for this run - STEP 1: Fill raw histograms and store amplitudes
        TString input_filename = Form("../../sorted_files/root_data_mvmelstrun%03d.bin_tree.root", run_number);
        cout << "Opening file: " << input_filename << endl;

        TFile* file1 = new TFile(input_filename, "READ");
        if (!file1 || file1->IsZombie()) {
            cerr << "ERROR: Could not open file for run " << run_number << endl;
            continue;
        }

        TTree *data = dynamic_cast<TTree*>(file1->Get("clover"));
        if (!data) {
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
        if (oneper == 0) oneper = 1;

        cout << "Reading data and storing amplitudes..." << endl;
        while (reader.Next()) {
            for (int crystal = 0; crystal < nCrystals; ++crystal) {
                double Ej = clover_cross_amplitude.At(crystal);
                if (Ej > 10) {
                    // STEP 1: Fill raw histograms and store amplitudes for later calibration
                    raw_histograms[run_idx][crystal]->Fill(Ej);
                    raw_amplitudes[run_idx][crystal].push_back(Ej);
                }
            }

            ++count;
            
            if (count % oneper == 0) {
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
            800, 600
        );
        calib_canvas->Divide(4, 4);

        // Create 4x4 canvas for linear fit plots
        TCanvas *linear_fit_canvas = new TCanvas(
            Form("linear_fit_run%d", run_number),
            Form("Linear Fits - Run %d", run_number),
            800, 600
        );
        linear_fit_canvas->Divide(4, 4);
        
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            cout << "Calibrating crystal " << crystal << "..." << endl;
            calib_canvas->cd(crystal + 1);
            TH1D* h = raw_histograms[run_idx][crystal];
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
            std::vector<Double_t> xPeaks = {0.0, 0.0, 0.0};
            // Fill vector carefully with checks
            xPeaks[0] = (xPeaks_511 && xPeaks_511[0] > 0) ? xPeaks_511[0] : 0.0;
            xPeaks[1] = (xPeaks_1460 && xPeaks_1460[0] > 0) ? xPeaks_1460[0] : 0.0;
            xPeaks[2] = (xPeaks_2640 && xPeaks_2640[0] > 0) ? xPeaks_2640[0] : 0.0;
            
            double x[2];
            double y[2] = {1460, 2614};
            
            // Fit 1460 and 2614 for calibration AND extract gaussian sigma for FWHM
            for (int i = 1; i < nPeaks; ++i) {
                if (xPeaks[i] <= 0) {
                    cout << "  Warning: peak " << i << " not found for crystal " << crystal << " in run " << run_number << endl;
                    continue;
                }
                
                // Fit gaus + pol1 around peak
                TF1 *gauss = new TF1(Form("gauss%d_%d", crystal, i), "gaus(0)+pol1(3)", xPeaks[i] - 10, xPeaks[i] + 10);
                // set initial guesses
                gauss->SetParameter(0, h->GetBinContent(h->FindBin(xPeaks[i]))); // amplitude
                gauss->SetParameter(1, xPeaks[i]); // mean
                gauss->SetParameter(2, 2.0); // sigma
                gauss->SetParameter(3, 0.0); // pol1 const
                gauss->SetParameter(4, 0.0); // pol1 slope

                int stdout_fd = silence_stdout();
                h->Fit(gauss, "MLISRNQ"); // quiet fit
                restore_stdout(stdout_fd);

                // Use final refined fit
                TF1 *gauss_final = new TF1(Form("gauss_final%d_%d", crystal, i), "gaus(0)+pol1(3)", xPeaks[i] - 10, xPeaks[i] + 10);
                gauss_final->SetParameters(gauss->GetParameters());

                stdout_fd = silence_stdout();
                h->Fit(gauss_final, "MLISRNQ");
                restore_stdout(stdout_fd);

                x[i-1] = gauss_final->GetParameter(1);
                double sigma = gauss_final->GetParameter(2);
                double sigma_err = gauss_final->GetParError(2);
                double fwhm = 2.35482 * sigma;
                double fwhm_err = 2.35482 * sigma_err;

                // Save FWHM per peak
                if (i == 1) { // 1460
                    fwhm_1460_individual[run_idx][crystal] = fwhm;
                    fwhm_1460_err_individual[run_idx][crystal] = fwhm_err;
                } else if (i == 2) { // 2614
                    fwhm_2614_individual[run_idx][crystal] = fwhm;
                    fwhm_2614_err_individual[run_idx][crystal] = fwhm_err;
                }

                gauss_final->SetLineColor(i + 2);
                gauss_final->Draw("same");
                
                delete gauss;
            }
            
            // Additionally fit 511 peak (for FWHM only) if found
            if (xPeaks[0] > 0) {
                TF1 *gauss511 = new TF1(Form("gauss511_%d_%d", run_idx, crystal), "gaus(0)+pol1(3)", xPeaks[0] - 10, xPeaks[0] + 10);
                gauss511->SetParameter(0, h->GetBinContent(h->FindBin(xPeaks[0])));
                gauss511->SetParameter(1, xPeaks[0]);
                gauss511->SetParameter(2, 2.0);
                gauss511->SetParameter(3, 0.0);
                gauss511->SetParameter(4, 0.0);

                int stdout_fd = silence_stdout();
                h->Fit(gauss511, "MLISRNQ");
                restore_stdout(stdout_fd);

                TF1 *gauss511_final = new TF1(Form("gauss511_final%d_%d", run_idx, crystal), "gaus(0)+pol1(3)", xPeaks[0] - 10, xPeaks[0] + 10);
                gauss511_final->SetParameters(gauss511->GetParameters());
                
                stdout_fd = silence_stdout();
                h->Fit(gauss511_final, "MLISRNQ");
                restore_stdout(stdout_fd);
                
                double sigma511 = gauss511_final->GetParameter(2);
                double sigma511_err = gauss511_final->GetParError(2);
                fwhm_511_individual[run_idx][crystal] = 2.35482 * sigma511;
                fwhm_511_err_individual[run_idx][crystal] = 2.35482 * sigma511_err;
                
                gauss511_final->SetLineColor(kGreen);
                gauss511_final->Draw("same");
                delete gauss511;
            } else {
                // Not found, leave zero
            }
            
            cout << "Crystal #" << crystal;
            if (xPeaks[1] > 0 && xPeaks[2] > 0) {
                cout << "  1460 keV centroid @ channel = " << x[0]
                     << " ; 2614 keV centroid @ channel = " << x[1] << endl;
            } else {
                cout << "  (warning: centroids missing for this crystal/run)" << endl;
            }
            
            // Create and fit the linear calibration graph using the centroids for 1460 and 2614
            // Only perform calibration if both centroids available
            if (xPeaks[1] > 0 && xPeaks[2] > 0) {
                double xvals[2] = { x[0], x[1] };
                TGraph* graph = new TGraph(2, xvals, y);
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
                
                // Add marker for 511 keV peak in raw channel coordinates (optional)
                if (xPeaks_511 && xPeaks_511[0] > 0) {
                    double energy_511 = linear->GetParameter(1) * xPeaks_511[0] + linear->GetParameter(0);
                    TMarker *m511 = new TMarker(xPeaks_511[0], energy_511, 29);
                    m511->SetMarkerColor(kGreen);
                    m511->SetMarkerSize(1.5);
                    m511->Draw();
                }
                
                // Turn GRID OFF for linear fit panels (user requested no grid lines)
                gPad->SetGrid(0, 0);
            } else {
                // If calibration can't be done, put zeros in calibration params so later calibration step doesn't crash
                calibration_params[run_idx][crystal] = make_pair(0.0, 0.0);
            }
        }
        
        // Save linear fit canvas
        linear_fit_canvas->Draw();
        linear_fit_canvases.push_back(linear_fit_canvas);
        
        // STEP 3: Apply calibration using stored raw amplitudes (avoiding bin artifacts)
        cout << "Applying calibration to stored amplitudes for run " << run_number << "..." << endl;

        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            double slope = calibration_params[run_idx][crystal].first;
            double intercept = calibration_params[run_idx][crystal].second;
            
            int crystal_events = raw_amplitudes[run_idx][crystal].size();
            cout << "Calibrating crystal " << crystal << " with " << crystal_events << " events..." << endl;
            
            // Apply calibration to each stored amplitude
            for (double amplitude : raw_amplitudes[run_idx][crystal]) {
                double energy = slope * amplitude + intercept;
                cal_histograms[run_idx][crystal]->Fill(energy);
                
                // Add to sum histograms
                sum_raw_histograms[crystal]->Fill(amplitude);
                sum_cal_histograms[crystal]->Fill(energy);
            }
        }

        cout << "Finished calibrating run " << run_number << endl;
        
        // If any FWHM values haven't been captured for 1460/2614 via the fit (e.g. fit failed),
        // fall back to the histogram-based half-max method for robustness.
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            if (fwhm_1460_individual[run_idx][crystal] <= 0.0) {
                double fallback = calculateFWHM(cal_histograms[run_idx][crystal], 1460.0, 50.0);
                fwhm_1460_individual[run_idx][crystal] = fallback;
                fwhm_1460_err_individual[run_idx][crystal] = 0.0;
            }
            if (fwhm_2614_individual[run_idx][crystal] <= 0.0) {
                double fallback = calculateFWHM(cal_histograms[run_idx][crystal], 2614.0, 50.0);
                fwhm_2614_individual[run_idx][crystal] = fallback;
                fwhm_2614_err_individual[run_idx][crystal] = 0.0;
            }
            if (fwhm_511_individual[run_idx][crystal] <= 0.0) {
                double fallback = calculateFWHM(cal_histograms[run_idx][crystal], 511.0, 30.0);
                fwhm_511_individual[run_idx][crystal] = fallback;
                fwhm_511_err_individual[run_idx][crystal] = 0.0;
            }
            cout << "Crystal " << crystal << ": 511 keV FWHM = " << fwhm_511_individual[run_idx][crystal]
                 << " ± " << fwhm_511_err_individual[run_idx][crystal]
                 << " ; 1460 keV FWHM = " << fwhm_1460_individual[run_idx][crystal]
                 << " ± " << fwhm_1460_err_individual[run_idx][crystal]
                 << " ; 2614 keV FWHM = " << fwhm_2614_individual[run_idx][crystal]
                 << " ± " << fwhm_2614_err_individual[run_idx][crystal] << " keV" << endl;
        }
        
        // Clear amplitudes to save memory for next run
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            raw_amplitudes[run_idx][crystal].clear();
            raw_amplitudes[run_idx][crystal].shrink_to_fit();
        }
        
        // Save individual run to ROOT file (including linear fit canvas)
        TString outname = Form("../clover_analysis/cal_sum/calibration_run%03d_linear.root", run_number);
        cout << "Saving results to: " << outname << endl;

        TFile *output = new TFile(outname, "RECREATE");
        TDirectory* dirRaw = output->mkdir("raw");
        TDirectory* dirCal = output->mkdir("cal");
        TDirectory* dirFits = output->mkdir("linear_fits");

        // Write raw histograms
        dirRaw->cd();
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            raw_histograms[run_idx][crystal]->Write();
        }

        // Write calibrated histograms
        dirCal->cd();
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            cal_histograms[run_idx][crystal]->Write();
        }

        // Write linear fit canvas
        dirFits->cd();
        linear_fit_canvas->Write();

        output->Close();
        cout << "Saved results for run " << run_number << endl;
        
        std::ofstream param_file(Form("../clover_analysis/calibration_run%d.txt", run_number));
        for (int c = 0; c < nCrystals; c++) {
            param_file << "Crystal " << c
                       << " slope = " << calibration_params[run_idx][c].first
                       << " intercept = " << calibration_params[run_idx][c].second
                       << std::endl;
        }
        param_file.close();

        std::cout << "Finished processing run " << run_number << std::endl;
        delete calib_canvas;
    }
    
    // Calculate FWHM for summed histograms by fitting gaus+pol1 around each peak
    cout << "Calculating FWHM for summed histograms (fits)..." << endl;
    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        TH1D* hs = sum_cal_histograms[crystal];
        if (!hs) continue;
        
        // Fit 511
        {
            // pick a local search region
            double low = 480.0, high = 540.0;
            TF1 *f511 = new TF1(Form("sum_gauss511_%d", crystal), "gaus(0)+pol1(3)", low, high);
            // initial guesses (use bin max near expected)
            int bin511 = hs->FindBin(511.0);
            f511->SetParameter(0, hs->GetBinContent(bin511));
            f511->SetParameter(1, 511.0);
            f511->SetParameter(2, 2.0);
            f511->SetParameter(3, 0.0);
            f511->SetParameter(4, 0.0);

            int stdout_fd = silence_stdout();
            hs->Fit(f511, "RMLISNQ");
            restore_stdout(stdout_fd);

            double sigma = f511->GetParameter(2);
            double sigma_err = f511->GetParError(2);
            fwhm_511_sum[crystal] = 2.35482 * sigma;
            fwhm_511_sum_err[crystal] = 2.35482 * sigma_err;
            delete f511;
        }

        // Fit 1460
        {
            double low = 1400.0, high = 1520.0;
            TF1 *f1460 = new TF1(Form("sum_gauss1460_%d", crystal), "gaus(0)+pol1(3)", low, high);
            int bin1460 = hs->FindBin(1460.0);
            f1460->SetParameter(0, hs->GetBinContent(bin1460));
            f1460->SetParameter(1, 1460.0);
            f1460->SetParameter(2, 2.0);
            f1460->SetParameter(3, 0.0);
            f1460->SetParameter(4, 0.0);

            int stdout_fd = silence_stdout();
            hs->Fit(f1460, "RMLISNQ");
            restore_stdout(stdout_fd);

            double sigma = f1460->GetParameter(2);
            double sigma_err = f1460->GetParError(2);
            fwhm_1460_sum[crystal] = 2.35482 * sigma;
            fwhm_1460_sum_err[crystal] = 2.35482 * sigma_err;
            delete f1460;
        }

        // Fit 2614
        {
            double low = 2560.0, high = 2660.0;
            TF1 *f2614 = new TF1(Form("sum_gauss2614_%d", crystal), "gaus(0)+pol1(3)", low, high);
            int bin2614 = hs->FindBin(2614.0);
            f2614->SetParameter(0, hs->GetBinContent(bin2614));
            f2614->SetParameter(1, 2614.0);
            f2614->SetParameter(2, 2.0);
            f2614->SetParameter(3, 0.0);
            f2614->SetParameter(4, 0.0);

            int stdout_fd = silence_stdout();
            hs->Fit(f2614, "RMLISNQ");
            restore_stdout(stdout_fd);

            double sigma = f2614->GetParameter(2);
            double sigma_err = f2614->GetParError(2);
            fwhm_2614_sum[crystal] = 2.35482 * sigma;
            fwhm_2614_sum_err[crystal] = 2.35482 * sigma_err;
            delete f2614;
        }

        cout << "Summed Crystal " << crystal << ": 511 FWHM = " << fwhm_511_sum[crystal] << " ± " << fwhm_511_sum_err[crystal]
             << " ; 1460 FWHM = " << fwhm_1460_sum[crystal] << " ± " << fwhm_1460_sum_err[crystal]
             << " ; 2614 FWHM = " << fwhm_2614_sum[crystal] << " ± " << fwhm_2614_sum_err[crystal] << endl;
    }
    
    // Save linear fit canvases to separate files
    for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
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
    for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
        int run_number = start_run + run_idx;
        
        TCanvas* cRaw = new TCanvas(Form("Run%d_Raw", run_number), Form("Run %d Raw", run_number), 800, 600);
        cRaw->Divide(4,4);
        
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            cRaw->cd(crystal+1);
            TH1D* hraw = raw_histograms[run_idx][crystal];
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
    for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
        int run_number = start_run + run_idx;
        
        TCanvas* cCal = new TCanvas(Form("Run%d_Cal", run_number), Form("Run %d Calibrated", run_number), 800, 600);
        cCal->Divide(4,4);
        
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            cCal->cd(crystal+1);
            TH1D* hcal = cal_histograms[run_idx][crystal];
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
    TCanvas* cOR = new TCanvas("OverlayRawRuns", "Overlay Raw Runs", 800, 600);
    cOR->Divide(4,4);

    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        cOR->cd(crystal+1);
        
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        bool first = true;
        
        for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
            int run_number = start_run + run_idx;
            TH1D* hraw = raw_histograms[run_idx][crystal];
            hraw->SetFillStyle(0);
            hraw->SetLineWidth(2);
            hraw->SetLineColor(colors[run_idx % 16]);
            
            if (first) {
                hraw->Draw("HIST");
                first = false;
            } else {
                hraw->Draw("HIST SAME");
            }
            
            legend->AddEntry(hraw, Form("Run %d", run_number), "l");
        }
        
        legend->Draw();
        gPad->SetLogy();
    }

    // Save overlay to sum file
    TString sum_outname = Form("../clover_analysis/cal_sum/sum_calibration_runs_%03d_to_%03d.root", start_run, end_run);
    TFile *sum_output = new TFile(sum_outname, "RECREATE");
    cOR->Write();
    sum_output->Close();

    overlay_raw_canvases.push_back(cOR);
    cOR->Draw();

    // 4. Overlay calibrated energy across runs (4x4 grid)
    cout << "Creating overlay calibrated plots..." << endl;
    TCanvas* cOC = new TCanvas("OverlayCalRuns", "Overlay Calibrated Runs", 800, 600);
    cOC->Divide(4,4);

    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        cOC->cd(crystal+1);
        
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        bool first = true;
        
        for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
            int run_number = start_run + run_idx;
            TH1D* hcal = cal_histograms[run_idx][crystal];
            hcal->SetFillStyle(0);
            hcal->SetLineWidth(2);
            hcal->SetLineColor(colors[run_idx % 16]);
            
            if (first) {
                hcal->Draw("HIST");
                first = false;
            } else {
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
    TCanvas* cSR = new TCanvas("SumRaw", "Summed Raw Energies", 800, 600);
    cSR->Divide(4,4);

    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        cSR->cd(crystal+1);
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
    TCanvas* cSC = new TCanvas("SumCal", "Summed Calibrated Energies", 800, 600);
    cSC->Divide(4,4);

    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        cSC->cd(crystal+1);
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
    TCanvas* cSO = new TCanvas("SumOverlay", "Overlay Summed Raw vs Calibrated", 800, 600);
    cSO->Divide(4,4);

    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        cSO->cd(crystal+1);
        
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
    
    // 7b. NEW: Overlay summed calibrated spectra vs one individual calibrated run (first run = start_run)
    cout << "Creating overlay: summed calibrated vs individual run " << start_run << "..." << endl;
    TCanvas* cSumVsRun = new TCanvas("SumCal_vs_Run", Form("Summed Calibrated vs Run %d", start_run), 800, 600);
    cSumVsRun->Divide(4,4);

    int overlay_run_idx = 0; // use the first run processed (start_run)
    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        cSumVsRun->cd(crystal+1);

        TH1D* hSum = sum_cal_histograms[crystal];
        TH1D* hRun = cal_histograms[overlay_run_idx][crystal];

        // If the run calibration is missing (all zeros), still draw sum
        bool hasRunHist = (hRun != nullptr);

        // Draw summed calibrated first
        if (hSum) {
            hSum->SetLineColor(kBlue);
            hSum->SetLineWidth(2);
            hSum->SetFillStyle(0);
            hSum->SetTitle(Form("Sum Cal vs Run %d - Crystal %d;Energy (keV);Counts", start_run, crystal));
            hSum->Draw("HIST");
        }

        // Draw individual run on top with distinct color
        if (hasRunHist) {
            hRun->SetLineColor(colors[overlay_run_idx % 16]);
            hRun->SetLineWidth(2);
            hRun->SetFillStyle(0);
            if (hSum) hRun->Draw("HIST SAME");
            else hRun->Draw("HIST");
        }

        // Legend
        TLegend* leg = new TLegend(0.65, 0.65, 0.92, 0.9);
        if (hSum) leg->AddEntry(hSum, "Summed Calibrated", "l");
        if (hasRunHist) leg->AddEntry(hRun, Form("Run %d Calibrated", start_run), "l");
        leg->Draw();
        gPad->SetLogy();
    }

    // Save this canvas to the sum file
    sum_output = new TFile(sum_outname, "UPDATE");
    cSumVsRun->Write();
    sum_output->Close();
    sum_canvases.push_back(cSumVsRun);
    cSumVsRun->Draw();

    // 8. FWHM comparison plots: Individual runs vs Summed (with error bars)
    cout << "Creating FWHM comparison plots with error bars..." << endl;

    // Create canvases for FWHM comparison (one per peak)
    TCanvas* cFWHM_511 = new TCanvas("FWHM_511_Comparison", "511 keV FWHM Comparison", 1200, 800);
    TCanvas* cFWHM_1460 = new TCanvas("FWHM_1460_Comparison", "1460 keV FWHM Comparison", 1200, 800);
    TCanvas* cFWHM_2614 = new TCanvas("FWHM_2614_Comparison", "2614 keV FWHM Comparison", 1200, 800);

    cFWHM_511->Divide(4, 4);
    cFWHM_1460->Divide(4, 4);
    cFWHM_2614->Divide(4, 4);

    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        // Build TGraphErrors for each peak
        TGraphErrors* g511 = new TGraphErrors();
        TGraphErrors* g1460 = new TGraphErrors();
        TGraphErrors* g2614 = new TGraphErrors();
        
        int idx511 = 0, idx1460 = 0, idx2614 = 0;
        for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
            int run_number = start_run + run_idx;
            double x = run_number;
            // 511
            if (fwhm_511_individual[run_idx][crystal] > 0) {
                g511->SetPoint(idx511, x, fwhm_511_individual[run_idx][crystal]);
                g511->SetPointError(idx511, 0.0, fwhm_511_err_individual[run_idx][crystal]);
                idx511++;
            }
            // 1460
            if (fwhm_1460_individual[run_idx][crystal] > 0) {
                g1460->SetPoint(idx1460, x, fwhm_1460_individual[run_idx][crystal]);
                g1460->SetPointError(idx1460, 0.0, fwhm_1460_err_individual[run_idx][crystal]);
                idx1460++;
            }
            // 2614
            if (fwhm_2614_individual[run_idx][crystal] > 0) {
                g2614->SetPoint(idx2614, x, fwhm_2614_individual[run_idx][crystal]);
                g2614->SetPointError(idx2614, 0.0, fwhm_2614_err_individual[run_idx][crystal]);
                idx2614++;
            }
        }
        
        // 511 canvas
        cFWHM_511->cd(crystal + 1);
        double min_fwhm = 1e6, max_fwhm = 0.0;
        for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
            if (fwhm_511_individual[run_idx][crystal] > 0) {
                min_fwhm = min(min_fwhm, fwhm_511_individual[run_idx][crystal]);
                max_fwhm = max(max_fwhm, fwhm_511_individual[run_idx][crystal]);
            }
        }
        if (fwhm_511_sum[crystal] > 0) {
            min_fwhm = min(min_fwhm, fwhm_511_sum[crystal] * 0.9);
            max_fwhm = max(max_fwhm, fwhm_511_sum[crystal] * 1.1);
        }
        if (min_fwhm >= max_fwhm || min_fwhm == 1e6) {
            min_fwhm = 0.0;
            max_fwhm = 100.0;
        }
        TH1F* frame_511 = new TH1F(Form("frame_511_crystal_%d", crystal), 
                                   Form("511 keV FWHM - Crystal %d;Run Number;FWHM (keV)", crystal),
                                   1, start_run - 0.5, end_run + 0.5);
        frame_511->SetMinimum(min_fwhm);
        frame_511->SetMaximum(max_fwhm);
        frame_511->Draw();
        g511->SetMarkerStyle(20);
        g511->SetMarkerSize(1.0);
        g511->SetMarkerColor(kBlue);
        g511->Draw("P SAME");
        // Summed FWHM line
        if (fwhm_511_sum[crystal] > 0) {
            TLine* l511 = new TLine(start_run - 0.5, fwhm_511_sum[crystal], end_run + 0.5, fwhm_511_sum[crystal]);
            l511->SetLineColor(kRed);
            l511->SetLineWidth(2);
            l511->SetLineStyle(2);
            l511->Draw();
        }
        TLegend* leg511 = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg511->AddEntry(g511, "Individual Runs", "p");
        leg511->AddEntry((TObject*)0, "Summed (red dashed)", "");
        leg511->Draw();
        // Turn GRID OFF for FWHM panels
        gPad->SetGrid(0, 0);
        
        // 1460 canvas
        cFWHM_1460->cd(crystal + 1);
        min_fwhm = 1e6; max_fwhm = 0.0;
        for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
            if (fwhm_1460_individual[run_idx][crystal] > 0) {
                min_fwhm = min(min_fwhm, fwhm_1460_individual[run_idx][crystal]);
                max_fwhm = max(max_fwhm, fwhm_1460_individual[run_idx][crystal]);
            }
        }
        if (fwhm_1460_sum[crystal] > 0) {
            min_fwhm = min(min_fwhm, fwhm_1460_sum[crystal] * 0.9);
            max_fwhm = max(max_fwhm, fwhm_1460_sum[crystal] * 1.1);
        }
        if (min_fwhm >= max_fwhm || min_fwhm == 1e6) {
            min_fwhm = 0.0;
            max_fwhm = 100.0;
        }
        TH1F* frame_1460 = new TH1F(Form("frame_1460_crystal_%d", crystal), 
                                   Form("1460 keV FWHM - Crystal %d;Run Number;FWHM (keV)", crystal),
                                   1, start_run - 0.5, end_run + 0.5);
        frame_1460->SetMinimum(min_fwhm);
        frame_1460->SetMaximum(max_fwhm);
        frame_1460->Draw();
        g1460->SetMarkerStyle(20);
        g1460->SetMarkerSize(1.0);
        g1460->SetMarkerColor(kBlue);
        g1460->Draw("P SAME");
        if (fwhm_1460_sum[crystal] > 0) {
            TLine* l1460 = new TLine(start_run - 0.5, fwhm_1460_sum[crystal], end_run + 0.5, fwhm_1460_sum[crystal]);
            l1460->SetLineColor(kRed);
            l1460->SetLineWidth(2);
            l1460->SetLineStyle(2);
            l1460->Draw();
        }
        TLegend* leg1460 = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg1460->AddEntry(g1460, "Individual Runs", "p");
        leg1460->Draw();
        // Turn GRID OFF for FWHM panels
        gPad->SetGrid(0, 0);
        
        // 2614 canvas
        cFWHM_2614->cd(crystal + 1);
        min_fwhm = 1e6; max_fwhm = 0.0;
        for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
            if (fwhm_2614_individual[run_idx][crystal] > 0) {
                min_fwhm = min(min_fwhm, fwhm_2614_individual[run_idx][crystal]);
                max_fwhm = max(max_fwhm, fwhm_2614_individual[run_idx][crystal]);
            }
        }
        if (fwhm_2614_sum[crystal] > 0) {
            min_fwhm = min(min_fwhm, fwhm_2614_sum[crystal] * 0.9);
            max_fwhm = max(max_fwhm, fwhm_2614_sum[crystal] * 1.1);
        }
        if (min_fwhm >= max_fwhm || min_fwhm == 1e6) {
            min_fwhm = 0.0;
            max_fwhm = 100.0;
        }
        TH1F* frame_2614 = new TH1F(Form("frame_2614_crystal_%d", crystal), 
                                   Form("2614 keV FWHM - Crystal %d;Run Number;FWHM (keV)", crystal),
                                   1, start_run - 0.5, end_run + 0.5);
        frame_2614->SetMinimum(min_fwhm);
        frame_2614->SetMaximum(max_fwhm);
        frame_2614->Draw();
        g2614->SetMarkerStyle(20);
        g2614->SetMarkerSize(1.0);
        g2614->SetMarkerColor(kBlue);
        g2614->Draw("P SAME");
        if (fwhm_2614_sum[crystal] > 0) {
            TLine* l2614 = new TLine(start_run - 0.5, fwhm_2614_sum[crystal], end_run + 0.5, fwhm_2614_sum[crystal]);
            l2614->SetLineColor(kRed);
            l2614->SetLineWidth(2);
            l2614->SetLineStyle(2);
            l2614->Draw();
        }
        TLegend* leg2614 = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg2614->AddEntry(g2614, "Individual Runs", "p");
        leg2614->Draw();
        // Turn GRID OFF for FWHM panels
        gPad->SetGrid(0, 0);
        
        // Save graphs for later writing to ROOT
        g511->SetName(Form("g511_crystal_%d", crystal));
        g1460->SetName(Form("g1460_crystal_%d", crystal));
        g2614->SetName(Form("g2614_crystal_%d", crystal));
        
        // Write the graphs into the sum output file immediately (so they are stored)
        sum_output = new TFile(sum_outname, "UPDATE");
        g511->Write();
        g1460->Write();
        g2614->Write();
        sum_output->Close();
    }

    // Save FWHM comparison canvases
    sum_output = new TFile(sum_outname, "UPDATE");
    cFWHM_511->Write();
    cFWHM_1460->Write();
    cFWHM_2614->Write();
    sum_output->Close();

    // Draw the canvases
    cFWHM_511->Draw();
    cFWHM_1460->Draw();
    cFWHM_2614->Draw();

    // Add to canvas vectors for memory management
    sum_canvases.push_back(cFWHM_511);
    sum_canvases.push_back(cFWHM_1460);
    sum_canvases.push_back(cFWHM_2614);
    
    // Save calibration parameters
    ofstream calib_file("../clover_analysis/cal_sum/clover_calibration_parameters.txt");
    if (calib_file.is_open()) {
        calib_file << "# Run CrystalIndex Slope(a) Intercept(b)\n";
        for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
            int run_number = start_run + run_idx;
            for (int crystal = 0; crystal < nCrystals; ++crystal) {
                calib_file << run_number << " " << crystal << " " 
                          << calibration_params[run_idx][crystal].first << " " 
                          << calibration_params[run_idx][crystal].second << "\n";
            }
        }
        calib_file.close();
        cout << "Calibration parameters saved to clover_calibration_parameters_multirun.txt" << endl;
    } else {
        cerr << "Could not write calibration parameters.\n";
    }
    
    // Save FWHM results to file (averages ± std as before)
    ofstream fwhm_file("../clover_analysis/cal_sum/fwhm_comparison_results.txt");
    if (fwhm_file.is_open()) {
        fwhm_file << "# FWHM Comparison Results\n";
        fwhm_file << "# Crystal Individual_Runs_FWHM_511(keV) Summed_FWHM_511(keV) Individual_Runs_FWHM_1460(keV) Summed_FWHM_1460(keV) Individual_Runs_FWHM_2614(keV) Summed_FWHM_2614(keV)\n";
        
        for (int crystal = 0; crystal < nCrystals; ++crystal) {
            fwhm_file << crystal << " ";
            
            // Calculate average and std dev of individual runs for each peak
            double sum_511 = 0.0, sum_1460 = 0.0, sum_2614 = 0.0;
            double sum_sq_511 = 0.0, sum_sq_1460 = 0.0, sum_sq_2614 = 0.0;
            int count_511 = 0, count_1460 = 0, count_2614 = 0;
            
            for (int run_idx = 0; run_idx < n_runs; ++run_idx) {
                if (fwhm_511_individual[run_idx][crystal] > 0) {
                    sum_511 += fwhm_511_individual[run_idx][crystal];
                    sum_sq_511 += fwhm_511_individual[run_idx][crystal] * fwhm_511_individual[run_idx][crystal];
                    count_511++;
                }
                if (fwhm_1460_individual[run_idx][crystal] > 0) {
                    sum_1460 += fwhm_1460_individual[run_idx][crystal];
                    sum_sq_1460 += fwhm_1460_individual[run_idx][crystal] * fwhm_1460_individual[run_idx][crystal];
                    count_1460++;
                }
                if (fwhm_2614_individual[run_idx][crystal] > 0) {
                    sum_2614 += fwhm_2614_individual[run_idx][crystal];
                    sum_sq_2614 += fwhm_2614_individual[run_idx][crystal] * fwhm_2614_individual[run_idx][crystal];
                    count_2614++;
                }
            }
            
            double avg_511 = (count_511 > 0) ? sum_511 / count_511 : 0.0;
            double avg_1460 = (count_1460 > 0) ? sum_1460 / count_1460 : 0.0;
            double avg_2614 = (count_2614 > 0) ? sum_2614 / count_2614 : 0.0;
            double std_511 = (count_511 > 1) ? sqrt((sum_sq_511 / count_511) - (avg_511 * avg_511)) : 0.0;
            double std_1460 = (count_1460 > 1) ? sqrt((sum_sq_1460 / count_1460) - (avg_1460 * avg_1460)) : 0.0;
            double std_2614 = (count_2614 > 1) ? sqrt((sum_sq_2614 / count_2614) - (avg_2614 * avg_2614)) : 0.0;
            
            fwhm_file << avg_511 << "±" << std_511 << " " << fwhm_511_sum[crystal] << " "
                      << avg_1460 << "±" << std_1460 << " " << fwhm_1460_sum[crystal] << " "
                      << avg_2614 << "±" << std_2614 << " " << fwhm_2614_sum[crystal] << "\n";
        }
        
        fwhm_file.close();
        cout << "FWHM comparison results saved to fwhm_comparison_results.txt" << endl;
    } else {
        cerr << "Could not write FWHM comparison results.\n";
    }
    
    cout << "==========================================" << endl;
    cout << "Processing complete! All canvases are displayed." << endl;
    cout << "Individual run files saved in: ../clover_analysis/cal_sum/" << endl;
    cout << "Sum results saved in: " << sum_outname << endl;
    cout << "Press Ctrl+C to exit and close all canvases." << endl;
    cout << "==========================================" << endl;
}

// Main function
int main(int argc, char* argv[]) {
    // Parse command line arguments
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <start_run> <end_run>" << endl;
        cout << "Example: " << argv[0] << " 98 100" << endl;
        return 1;
    }
    
    int start_run = atoi(argv[1]);
    int end_run = atoi(argv[2]);
    
    if (start_run > end_run) {
        cout << "Error: start_run must be <= end_run" << endl;
        return 1;
    }
    
    // Create ROOT application to keep canvases open
    TApplication app("App", &argc, argv);
    
    // Run the calibration
    run_by_run_calibration(start_run, end_run);
    
    // Keep the application running
    app.Run();
    
    return 0;
}

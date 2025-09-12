// run_by_run_cal.C  (fixed and extended)
// NOTE: function name run_by_run_cal(int runNumber = 26) preserved exactly.

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>

using namespace std;

const int nCrystals = 16;

// Helper: load reference gamma energies from CSV (ignores header + label column)
vector<double> loadReferenceEnergies(const char* filename) {
    vector<double> refEnergies;
    ifstream infile(filename);
    if (!infile.is_open()) {
        // don't print here ? caller will decide
        return refEnergies;
    }

    string line;
    // optional skip header if present (safe)
    if (getline(infile, line)) {
        // if first line contains any letters assume it's header; otherwise rewind to parse it
        bool hasLetter = false;
        for (char c : line) if (isalpha((unsigned char)c)) { hasLetter = true; break; }
        if (hasLetter) {
            // header skipped, continue
        } else {
            // first line actually numeric ? parse it
            stringstream ss(line);
            string energyStr, label;
            if (getline(ss, energyStr, ',')) {
                try {
                    double E = stod(energyStr);
                    refEnergies.push_back(E);
                } catch (...) {}
            }
        }
    }

    while (getline(infile, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string energyStr;
        if (!getline(ss, energyStr, ',')) continue;
        // trim
        auto trim = [](string &s) {
            while (!s.empty() && isspace((unsigned char)s.front())) s.erase(0,1);
            while (!s.empty() && isspace((unsigned char)s.back())) s.pop_back();
        };
        trim(energyStr);
        if (energyStr.empty()) continue;
        try {
            double E = stod(energyStr);
            refEnergies.push_back(E);
        } catch (...) {
            continue;
        }
    }
    infile.close();
    sort(refEnergies.begin(), refEnergies.end());
    refEnergies.erase(unique(refEnergies.begin(), refEnergies.end()), refEnergies.end());
    return refEnergies;
}

// Find closest reference energy for a detected peak (keeps your original logic)
double findClosestReference(double rawCentroid, const vector<double>& refEnergies) {
    double bestE = -1;
    double minDiff = 1e9;
    for (auto E : refEnergies) {
        double diff = fabs(E - rawCentroid);
        if (diff < minDiff) {
            minDiff = diff;
            bestE = E;
        }
    }
    return bestE;
}

void run_by_run_cal(int runNumber = 26) {
    // --- load reference CSV (your original path) ---
    const char* refFile = "../clover_analysis/ref_energies.csv";
    vector<double> refEnergies = loadReferenceEnergies(refFile);
    if (refEnergies.empty()) {
        cerr << "Warning: No reference energies loaded from " << refFile
             << " ? falling back to defaults (511,1460,2614 keV)." << endl;
        refEnergies = {511.0, 1460.0, 2614.0};
    } else {
        cout << "Loaded " << refEnergies.size() << " reference energies from " << refFile << endl;
    }

    // --- open input file & tree ---
    TString inFile = Form("../../sorted_files/root_data_mvmelstrun%03d.bin_tree.root", runNumber);
    TFile* file1 = new TFile(inFile, "READ");
    if (!file1 || file1->IsZombie()) {
        cerr << "Error: could not open input file: " << inFile << endl;
        return;
    }

    TTree* data = dynamic_cast<TTree*>(file1->Get("clover"));
    if (!data) {
        cerr << "Tree 'clover' not found in input file!" << endl;
        file1->Close();
        return;
    }

    // --- configure reading of the original amplitude branch ---
    // NOTE: this assumes the branch is stored in a way compatible with a Double_t array per event.
    // If the branch is a std::vector<Double_t>, you'd need to change to vector<Double_t>* and SetBranchAddress accordingly.
    Double_t clover_cross[nCrystals];
    TBranch* br = data->GetBranch("clover_cross.amplitude");
    if (!br) {
        cerr << "ERROR: branch 'clover_cross.amplitude' not found in the tree." << endl;
        file1->Close();
        return;
    }
    data->SetBranchAddress("clover_cross.amplitude", clover_cross);

    // --- prepare raw histograms ---
    TH1D* hclover[nCrystals];
    for (int i = 0; i < nCrystals; ++i) {
        hclover[i] = new TH1D(Form("hclover_%d", i),
                              Form("Crystal %d (raw);Channel/Measured;Counts", i),
                              6000, 0, 3000);
        hclover[i]->SetDirectory(0);
    }

    // --- prepare calibrated histograms (keV) ---
    TH1D* hclover_cal[nCrystals];
    const int calBins = 3000;   // 1 keV/bin from 0..3000 keV
    const double calMin = 0.0;
    const double calMax = 3000.0;
    for (int i = 0; i < nCrystals; ++i) {
        hclover_cal[i] = new TH1D(Form("hclover_cal_%d", i),
                                  Form("Crystal %d (calibrated);Energy (keV);Counts", i),
                                  calBins, calMin, calMax);
        hclover_cal[i]->SetDirectory(0);
    }

    // --- create summed calibrated histogram (all crystals) ---
    TH1D* hsum = new TH1D("hsum", "Summed calibrated spectrum;Energy (keV);Counts", calBins, calMin, calMax);
    hsum->SetDirectory(0);

    // --- fill raw histograms with progress ---
    Long64_t nEntries = data->GetEntries();
    if (nEntries == 0) {
        cerr << "Tree has 0 entries ? nothing to do." << endl;
        file1->Close();
        return;
    }
    cout << "Filling raw histograms (" << nEntries << " entries)..." << endl;
    Long64_t onePercent = nEntries / 100;
    if (onePercent == 0) onePercent = 1;
    for (Long64_t i = 0; i < nEntries; ++i) {
        data->GetEntry(i);
        for (int j = 0; j < nCrystals; ++j) {
            if (clover_cross[j] > 10) hclover[j]->Fill(clover_cross[j]);
        }
        if (i % onePercent == 0) {
            int percent = static_cast<int>( (double)i / nEntries * 100.0 );
            cout << "\rFilling raw histograms: " << percent << "% completed" << flush;
        }
    }
    cout << "\rFilling raw histograms: 100% completed" << endl;

    // --- create canvases (raw, fit, calibrated) ---
    TCanvas* cHist = new TCanvas("cHist", "Raw Histograms", 1200, 800);
    cHist->Divide(4,4);
    TCanvas* cFit = new TCanvas("cFit", "Reference vs Raw Fits", 1200, 800);
    cFit->Divide(4,4);
    TCanvas* cCal  = new TCanvas("cCal",  "Calibrated Histograms (keV)", 1200, 800);
    cCal->Divide(4,4);

    // --- prepare CSV output for peaks & targets ---
    TString csvName = Form("../clover_analysis/peaks_run%03d.csv", runNumber);
    ofstream peaksCsv(csvName.Data());
    if (!peaksCsv.is_open()) {
        cerr << "WARNING: could not open " << csvName << " for writing peak list." << endl;
    } else {
        peaksCsv << "Run,Crystal,RawCentroid,MatchedRef_keV,Calibrated_keV,Slope,Intercept,IsTarget\n";
    }

    TString tCsvName = Form("../clover_analysis/targets_run%03d.csv", runNumber);
    ofstream targetsCsv(tCsvName.Data());
    if (!targetsCsv.is_open()) {
        cerr << "WARNING: could not open " << tCsvName << " for writing target matches." << endl;
    } else {
        targetsCsv << "Run,Crystal,TargetKEV,RawExpected,RawFound,Calibrated_keV,Distance,Status\n";
    }

    // --- calibration params ---
    Double_t slope[nCrystals], intercept[nCrystals];
    for (int i = 0; i < nCrystals; ++i) { slope[i] = 1.0; intercept[i] = 0.0; }

    // Keep graphs and fits alive so they remain visible on the pad
    vector<TGraph*> savedGraphs; savedGraphs.reserve(nCrystals);
    vector<TF1*>   savedFits;   savedFits.reserve(nCrystals);

    // target energies for vertical lines (keV) and colors
    const double targetE[3] = {511.0, 1460.0, 2614.0};
    const int targetColor[3] = {kRed, kBlue, kMagenta};

    // tolerance (in measured units, keV) to accept a match when checking targets
    const double tolKEV = 12.0;

    // --- loop crystals: find peaks, match to reference list and fit (reference vs raw) ---
    for (int crystal = 0; crystal < nCrystals; ++crystal) {
        // Draw raw histogram first on cHist pad
        cHist->cd(crystal+1);
        hclover[crystal]->Draw();
        gPad->SetLogy();

        TSpectrum spectrum(40);
        int nfound = spectrum.Search(hclover[crystal], 2, "nobackground", 0.05);
        if (nfound <= 0) {
            cerr << "Crystal " << crystal << ": TSpectrum found no peaks. Using identity calibration." << endl;
            slope[crystal] = 1.0;
            intercept[crystal] = 0.0;
            continue;
        }

        Double_t* xPeaks = spectrum.GetPositionX();
        vector<double> rawPeaks;
        for (int p = 0; p < nfound; ++p) rawPeaks.push_back(xPeaks[p]);
        sort(rawPeaks.begin(), rawPeaks.end());

        // --- REPLACEMENT: windowed TSpectrum + Gaussian-fitting for canonical targets (511,1460,2614)
        // We assume the input "measured" values in clover_cross.amplitude are already in keV
        // but may be slightly shifted. For each target energy we search a narrow window,
        // run TSpectrum inside that window, pick the best peak, then fit a Gaussian+linear
        // baseline to get a precise centroid. Those centroids (measured, reference_keV)
        // are then used for the final 3-point pol1 (or fallback) calibration.

        vector<double> matchedRaw, matchedRef; // here raw==measured (keV) values

        // window half-widths (keV) around the target to search for peaks
        const double winHalf[3] = {25.0, 40.0, 60.0};
        // fit half-widths (keV) around found peak for Gaussian+linear fit
        const double fitHalf[3] = {12.0, 20.0, 30.0};

        // store found centroids
        vector<double> foundMeas; foundMeas.reserve(3);
        vector<double> foundRef;  foundRef.reserve(3);

        for (int t = 0; t < 3; ++t) {
            double Etarget = targetE[t];
            double rmin = Etarget - winHalf[t];
            double rmax = Etarget + winHalf[t];
            // clone histogram and restrict range for localized search
            TH1D* hwin = (TH1D*)hclover[crystal]->Clone(Form("hwin_%d_%d", crystal, t));
            hwin->SetDirectory(0);
            hwin->GetXaxis()->SetRangeUser(rmin, rmax);

            // quick check: if region empty, skip
            double integral = hwin->Integral();
            if (integral <= 0) { delete hwin; continue; }

            // TSpectrum search inside the window
            TSpectrum spec(5);
            int nfound_win = spec.Search(hwin, 2, "nobackground", 0.05);
            double bestPos = -1; double bestAmp = -1;
            if (nfound_win > 0) {
                Double_t *px = spec.GetPositionX();
                for (int p = 0; p < nfound_win; ++p) {
                    double pos = px[p];
                    if (pos < rmin || pos > rmax) continue;
                    int bin = hwin->GetXaxis()->FindBin(pos);
                    double amp = hwin->GetBinContent(bin);
                    if (amp > bestAmp) { bestAmp = amp; bestPos = pos; }
                }
            }

            // If TSpectrum didn't find anything, try the local max in the window as fallback
            if (bestPos < 0) {
                int binmin = hwin->GetXaxis()->FindBin(rmin);
                int binmax = hwin->GetXaxis()->FindBin(rmax);
                if (binmax <= binmin) { delete hwin; continue; }
                for (int b = binmin; b <= binmax; ++b) {
                    double amp = hwin->GetBinContent(b);
                    if (amp > bestAmp) {
                        bestAmp = amp;
                        bestPos = hwin->GetXaxis()->GetBinCenter(b);
                    }
                }
            }

            if (bestPos > 0 && bestAmp > 0) {
                // prepare fit window around bestPos
                double fmin = std::max(bestPos - fitHalf[t], rmin);
                double fmax = std::min(bestPos + fitHalf[t], rmax);
                if (fmax <= fmin) { delete hwin; continue; }

                // Fit Gaussian + linear baseline directly on the original histogram hclover[crystal]
                TF1* gauss_lin = new TF1(Form("gauss_lin_%d_%d", crystal, t),
                                         "[0] + [1]*x + [2]*exp(-0.5*((x-[3])/[4])*((x-[3])/[4]))",
                                         fmin, fmax);

                // initial parameters: background ~ average of edges, slope 0, amplitude ~ peak-height, centroid bestPos, sigma guess
                double bkgL = hclover[crystal]->GetBinContent(hclover[crystal]->GetXaxis()->FindBin(fmin));
                double bkgR = hclover[crystal]->GetBinContent(hclover[crystal]->GetXaxis()->FindBin(fmax));
                double bkg = 0.5*(bkgL + bkgR);
                double ampGuess = bestAmp - bkg;
                if (ampGuess < 0) ampGuess = bestAmp;
                double sigmaGuess = max(1.0, fitHalf[t]/6.0);

                gauss_lin->SetParameter(0, bkg);
                gauss_lin->SetParameter(1, 0.0);
                gauss_lin->SetParameter(2, ampGuess);
                gauss_lin->SetParameter(3, bestPos);
                gauss_lin->SetParameter(4, sigmaGuess);
                // constrain sigma to avoid crazy fits
                gauss_lin->SetParLimits(4, 0.5, fitHalf[t]);

                // perform fit quietly; allow it to fail gracefully
                int fitStatus = hclover[crystal]->Fit(gauss_lin, "RQ", "", fmin, fmax);

                bool ok = (fitStatus == 0);
                double centroid = gauss_lin->GetParameter(3);
                double centroidErr = gauss_lin->GetParError(3);
                double amplitude = gauss_lin->GetParameter(2);

                // accept only reasonable fits: amplitude positive and centroid within the window
                if (ok && amplitude > 0 && centroid >= fmin && centroid <= fmax) {
                    foundMeas.push_back(centroid);
                    foundRef.push_back(Etarget);

                    // write this target detection to peaks CSV with IsTarget=1
                    if (peaksCsv.is_open()) {
                        double calibrated = centroid; // before calibration, measured==centroid
                        peaksCsv << runNumber << "," << crystal << "," << centroid << "," << Etarget << "," << calibrated << "," << slope[crystal] << "," << intercept[crystal] << ",1\n";
                    }

                    // draw marker on the raw hist pad
                    cHist->cd(crystal+1);
                    double ytop_local = hclover[crystal]->GetMaximum(); if (ytop_local <= 0) ytop_local = 1.0;
                    TLine *lt_found = new TLine(centroid, 0.5, centroid, ytop_local*1.05);
                    lt_found->SetLineColor(targetColor[t]); lt_found->SetLineWidth(2); lt_found->Draw();
                    TLatex *labf = new TLatex(centroid, ytop_local*1.05, Form("%.0f keV", Etarget));
                    labf->SetTextAlign(22); labf->SetTextSize(0.03); labf->Draw();

                    // annotate on fit pad later (after computing fit)
                }

                delete gauss_lin;
            }
            delete hwin;
        }

        // After searching for canonical targets, decide calibration points
        if (foundMeas.size() < 3) {
            // Not all three canonical targets found; try to supplement using whole-spectrum matches to build >=3 points
            vector<double> rawPeaksCopy = rawPeaks;
            // add any already found target centroids first
            for (size_t k = 0; k < foundMeas.size(); ++k) {
                matchedRaw.push_back(foundMeas[k]);
                matchedRef.push_back(foundRef[k]);
            }
            // now supplement with closest ref-energy matches from whole spectrum until we have 3 unique pairs
            for (double rp : rawPeaksCopy) {
                // avoid duplicates near already used centroids
                bool nearExisting = false;
                for (double used : matchedRaw) if (fabs(used - rp) < 1e-6) { nearExisting = true; break; }
                if (nearExisting) continue;
                double refE = findClosestReference(rp, refEnergies);
                if (refE > 0) {
                    matchedRaw.push_back(rp);
                    matchedRef.push_back(refE);
                    double calibrated = slope[crystal]*rp + intercept[crystal];
                    if (peaksCsv.is_open()) {
                        peaksCsv << runNumber << "," << crystal << "," << rp << "," << refE << "," << calibrated << "," << slope[crystal] << "," << intercept[crystal] << ",0\n";
                    }
                }
                if (matchedRaw.size() >= 3) break;
            }
        } else {
            // we already have three canonical target centroids; use them
            for (size_t k = 0; k < foundMeas.size(); ++k) {
                matchedRaw.push_back(foundMeas[k]);
                matchedRef.push_back(foundRef[k]);
            }
        }

        // --- If still not enough matches, give up and use identity calibration ---
        if (matchedRaw.size() < 3) {
            cerr << "Crystal " << crystal << ": insufficient matches for 3-point calibration (" << matchedRaw.size() << "). Using identity." << endl;
            slope[crystal] = 1.0;
            intercept[crystal] = 0.0;
            // write any matched pairs (if present) to CSV
            if (peaksCsv.is_open()) {
                for (size_t k = 0; k < matchedRaw.size(); ++k) {
                    double rawc = matchedRaw[k];
                    double refk = matchedRef[k];
                    double calibrated = slope[crystal]*rawc + intercept[crystal];
                    peaksCsv << runNumber << "," << crystal << "," << rawc << "," << refk << "," << calibrated << "," << slope[crystal] << "," << intercept[crystal] << ",0\n";
                }
            }
            continue;
        }

        // Fit: X = raw centroid (measured, keV), Y = reference energy (keV) -> ref = a + b*measured
        // --- enforce a 3-point linear regression (pol1) using exactly 3 matched points ---
        // Prefer using the three canonical targets if they are present in matchedRef
        vector<pair<double,double>> pairs;
        for (size_t k = 0; k < matchedRaw.size(); ++k) pairs.emplace_back(matchedRaw[k], matchedRef[k]);
        sort(pairs.begin(), pairs.end(), [](const pair<double,double>& a, const pair<double,double>& b){ return a.first < b.first; });

        // detect if canonical targets are present in the matched list (within tolKEV)
        bool hasTarget[3] = {false,false,false};
        pair<double,double> targetPairs[3];
        for (auto &pr : pairs) {
            for (int tt=0; tt<3; ++tt) {
                if (!hasTarget[tt] && fabs(pr.second - targetE[tt]) <= tolKEV) {
                    hasTarget[tt] = true;
                    targetPairs[tt] = pr;
                }
            }
        }

        vector<pair<double,double>> finalTriplet;
        if (hasTarget[0] && hasTarget[1] && hasTarget[2]) {
            finalTriplet.push_back(targetPairs[0]);
            finalTriplet.push_back(targetPairs[1]);
            finalTriplet.push_back(targetPairs[2]);
            cout << "Crystal " << crystal << ": using canonical targets for 3-point fit (511,1460,2614)\n";
        } else {
            // fallback: pick min, median, max in measured (raw) space
            size_t n = pairs.size();
            finalTriplet.push_back(pairs.front());
            finalTriplet.push_back(pairs[n/2]);
            finalTriplet.push_back(pairs.back());
            cout << "Crystal " << crystal << ": using min/median/max triplet for final 3-point fit\n";
        }

        // Build arrays for the 3 points and do the fit
        double xvals[3], yvals[3];
        for (int k = 0; k < 3; ++k) { xvals[k] = finalTriplet[k].first; yvals[k] = finalTriplet[k].second; }

        cFit->cd(crystal+1);
        TGraph* gr = new TGraph(3, xvals, yvals);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.2);
        gr->SetMarkerColor(kBlack);

        // compute axis padding from full matched pairs for context
        double xmin = pairs.front().first;
        double xmax = pairs.back().first;
        double ymin = pairs.front().second;
        double ymax = pairs.front().second;
        for (auto &pr : pairs) { ymin = min(ymin, pr.second); ymax = max(ymax, pr.second); }
        double xpad = max(5.0, 0.05*(xmax-xmin));
        double ypad = max(5.0, 0.05*(ymax-ymin));
        if (xpad <= 0) xpad = 10;
        if (ypad <= 0) ypad = 10;

        gr->GetXaxis()->SetLimits(xmin - xpad, xmax + xpad);
        gr->SetMinimum(ymin - ypad);
        gr->SetMaximum(ymax + ypad);

        gr->SetTitle(Form("Crystal %d: ref vs measured (3-point fit)", crystal));
        gr->GetXaxis()->SetTitle("Measured (keV)");
        gr->GetYaxis()->SetTitle("Reference Energy (keV)");
        gr->Draw("AP");

        // Fit only the 3-point graph (pol1)
        TF1* lin = new TF1(Form("lin_crystal_%d", crystal), "pol1");
        gr->Fit(lin, "Q"); // quiet fit

        intercept[crystal] = lin->GetParameter(0);
        slope[crystal]     = lin->GetParameter(1);

        // draw fit visibly and keep it alive
        lin->SetLineColor(kRed);
        lin->SetLineWidth(2);
        lin->Draw("same");

        // indicate which 3 points were used in the fit
        TPaveText *pt = new TPaveText(0.12, 0.68, 0.42, 0.92, "NDC");
        pt->SetBorderSize(0); pt->SetFillColor(0);
        pt->AddText(Form("3-point fit used (meas->keV):"));
        pt->AddText(Form("(%.2f -> %.2f), (%.2f -> %.2f), (%.2f -> %.2f)",
                        xvals[0], yvals[0], xvals[1], yvals[1], xvals[2], yvals[2]));
        pt->AddText(Form("slope = %.6f", slope[crystal]));
        pt->AddText(Form("intercept = %.2f", intercept[crystal]));
        pt->Draw();

        // keep graph & fit in vectors so they are not deleted and remain visible
        savedGraphs.push_back(gr);
        savedFits.push_back(lin);

        // Save matched peaks (non-target if any) to CSV (IsTarget=0 for these rows)
        if (peaksCsv.is_open()) {
            for (size_t k = 0; k < pairs.size(); ++k) {
                double rawc = pairs[k].first;
                double refk = pairs[k].second;
                double calibrated = slope[crystal]*rawc + intercept[crystal];
                peaksCsv << runNumber << "," << crystal << "," << rawc << "," << refk << "," << calibrated << "," << slope[crystal] << "," << intercept[crystal] << ",0\n";
            }
        }

        // --- redraw histogram pad so vertical lines are drawn on the raw histogram pad ---
        cHist->cd(crystal+1);
        hclover[crystal]->Draw();
        gPad->SetLogy();

        // Draw thin gray vertical lines at detected raw peaks (optional visual aid)
        double ytop = hclover[crystal]->GetMaximum();
        if (ytop <= 0) ytop = 1.0;
        for (double rp : rawPeaks) {
            TLine *lp = new TLine(rp, 0.5, rp, ytop * 1.05);
            lp->SetLineColor(kGray+2);
            lp->SetLineStyle(7);
            lp->SetLineWidth(1);
            lp->Draw();
        }

        // For each target energy, compute expected measured position (pre-correction) and try to find/mark it
        for (int t = 0; t < 3; ++t) {
            double E = targetE[t];
            double s = slope[crystal];
            double b = intercept[crystal];
            if (fabs(s) < 1e-9) {
                cout << "Run " << runNumber << " Crystal " << crystal << ": slope ~0, cannot project " << E << " keV." << endl;
                if (targetsCsv.is_open()) targetsCsv << runNumber << "," << crystal << "," << E << ",,,," << "NO_SLOPE\n";
                continue;
            }
            // project reference E back into measured space (measured = (E - intercept)/slope)
            double measured_expected = (E - b) / s;

            // find closest detected measured peak (rawPeaks are in measured units)
            if (!rawPeaks.empty()) {
                auto it = std::min_element(rawPeaks.begin(), rawPeaks.end(),
                      [measured_expected](double a, double bb){ return fabs(a - measured_expected) < fabs(bb - measured_expected); });
                double measured_found = *it;
                double dist = fabs(measured_found - measured_expected);

                if (dist <= (double) ( (t==0) ? 12.0 : (t==1 ? 25.0 : 40.0) ) ) {
                    double calib_found = slope[crystal]*measured_found + intercept[crystal];
                    cout << "Run " << runNumber << " Crystal " << crystal 
                         << " : FOUND target " << E << " keV near measured " << measured_found 
                         << " (expected " << measured_expected << ", dist=" << dist << "). Calibrated=" 
                         << calib_found << " keV\n";

                    // record target CSV
                    if (targetsCsv.is_open()) {
                        targetsCsv << runNumber << "," << crystal << "," << E << "," 
                                   << measured_expected << "," << measured_found << "," << calib_found << "," 
                                   << dist << ",FOUND\n";
                    }

                    // draw a thicker colored line at the found measured peak ON THE RAW HIST
                    TLine *lt_found = new TLine(measured_found, 0.5, measured_found, ytop * 1.05);
                    lt_found->SetLineColor(targetColor[t]);
                    lt_found->SetLineWidth(3);
                    lt_found->SetLineStyle(1);
                    lt_found->Draw();

                    TLatex *labf = new TLatex(measured_found, ytop * 1.05, Form("%.0f keV", E));
                    labf->SetTextAlign(22);
                    labf->SetTextSize(0.03);
                    labf->Draw();

                    // append a row marking this as a target in peaksCsv too
                    if (peaksCsv.is_open()) {
                        peaksCsv << runNumber << "," << crystal << "," << measured_found << "," << E << "," << calib_found << "," << slope[crystal] << "," << intercept[crystal] << ",1\n";
                    }

                    // --- add marker and dashed guide on the cFit pad to indicate this target ---
                    cFit->cd(crystal+1);
                    double fit_ymin = ymin - ypad;
                    double fit_ymax = ymax + ypad;
                    if (gr->GetHistogram()) {
                        fit_ymin = gr->GetHistogram()->GetMinimum();
                        fit_ymax = gr->GetHistogram()->GetMaximum();
                    }
                    TLine *guide = new TLine(measured_found, fit_ymin, measured_found, fit_ymax);
                    guide->SetLineStyle(2);
                    guide->SetLineColor(targetColor[t]);
                    guide->SetLineWidth(1);
                    guide->Draw();

                    // marker at (measured_found, E)
                    TMarker *m = new TMarker(measured_found, E, 20);
                    m->SetMarkerSize(1.2);
                    m->SetMarkerColor(targetColor[t]);
                    m->Draw();

                    // label near the marker
                    TLatex *labg = new TLatex(measured_found, E, Form("  %.0f keV", E));
                    labg->SetTextSize(0.03);
                    labg->SetTextColor(targetColor[t]);
                    labg->SetTextAlign(12); // left-bottom
                    labg->Draw();

                    // go back to histogram pad
                    cHist->cd(crystal+1);
                } else {
                    cout << "Run " << runNumber << " Crystal " << crystal 
                         << " : NO match for " << E << " keV (expected measured " << measured_expected 
                         << ", closest detected measured " << measured_found << " dist=" << dist << " > tol)." << endl;
                    if (targetsCsv.is_open()) targetsCsv << runNumber << "," << crystal << "," << E << "," 
                                                         << measured_expected << ",,," << dist << ",NO_MATCH\n";
                }
            } else {
                cout << "Run " << runNumber << " Crystal " << crystal << " : no detected peaks to match " << E << " keV." << endl;
                if (targetsCsv.is_open()) targetsCsv << runNumber << "," << crystal << "," << E << ",,," << ",,NO_PEAKS\n";
            }
        } // end targets loop
    } // end crystals loop

    // --- create output file and clone original tree structure (empty) ---
    TString outFile = Form("../sorted_files/root_data_corrected_run%03d.root", runNumber);
    TFile* fout = new TFile(outFile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        cerr << "ERROR: could not create output file: " << outFile << endl;
        file1->Close();
        if (peaksCsv.is_open()) peaksCsv.close();
        if (targetsCsv.is_open()) targetsCsv.close();
        return;
    }

    TTree* newTree = data->CloneTree(0); // empty clone keeps same branch structure

    // add requested new branch name: clover_cross_c.amplitude
    Double_t clover_cross_c[nCrystals];
    TBranch* newBranch = newTree->Branch("clover_cross_c.amplitude", clover_cross_c,
                                         Form("clover_cross_c.amplitude[%d]/D", nCrystals));
    if (!newBranch) {
        cerr << "ERROR: could not create branch 'clover_cross_c.amplitude'." << endl;
        fout->Close();
        file1->Close();
        if (peaksCsv.is_open()) peaksCsv.close();
        if (targetsCsv.is_open()) targetsCsv.close();
        return;
    }

    // --- Fill new tree AND calibrated histograms: read original tree, compute corrected array, fill clone ---
    cout << "Applying calibration and filling new tree & calibrated histograms (" << nEntries << " entries)..." << endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
        data->GetEntry(i);
        for (int j = 0; j < nCrystals; ++j) {
            clover_cross_c[j] = slope[j] * clover_cross[j] + intercept[j];
            // fill calibrated hist
            if (clover_cross_c[j] >= calMin && clover_cross_c[j] <= calMax) {
                hclover_cal[j]->Fill(clover_cross_c[j]);
                // fill summed calibrated histogram
                hsum->Fill(clover_cross_c[j]);
            } else {
                // optional: clamp to under/overflow bin (here we skip)
            }
        }
        newTree->Fill();

        if (i % onePercent == 0) {
            int percent = static_cast<int>( (double)i / nEntries * 100.0 );
            cout << "\rProcessing: " << percent << "% completed" << flush;
        }
    }
    cout << "\rProcessing: 100% completed" << endl;

    // Write everything
    fout->cd();
    // save histograms (raw)
    TDirectory* dRaw = fout->mkdir("raw");
    dRaw->cd();
    for (int j=0; j<nCrystals; ++j) hclover[j]->Write();
    // save calibrated histograms
    fout->cd();
    TDirectory* dCal = fout->mkdir("calibrated");
    dCal->cd();
    for (int j=0; j<nCrystals; ++j) hclover_cal[j]->Write();
    // save summed calibrated histogram
    hsum->Write();
    // save tree & new branch
    fout->cd();
    newTree->Write();
    fout->Close();

    // close input
    file1->Close();

    // close CSVs
    if (peaksCsv.is_open()) {
        peaksCsv.close();
        cout << "Saved peaks CSV: " << (string)csvName << endl;
    }
    if (targetsCsv.is_open()) {
        targetsCsv.close();
        cout << "Saved targets CSV: " << (string)tCsvName << endl;
    }

    // Draw calibrated histograms with target lines on cCal
    for (int i = 0; i < nCrystals; ++i) {
        cCal->cd(i+1);
        hclover_cal[i]->Draw();
        gPad->SetLogy();
        double ytop = hclover_cal[i]->GetMaximum();
        if (ytop <= 0) ytop = 1.0;
        for (int t = 0; t < 3; ++t) {
            double E = targetE[t];
            if (E >= calMin && E <= calMax) {
                TLine *lt = new TLine(E, 0.5, E, ytop * 1.05);
                lt->SetLineColor(targetColor[t]);
                lt->SetLineWidth(2);
                lt->Draw();
                TLatex *lab = new TLatex(E, ytop * 1.05, Form("%.0f keV", E));
                lab->SetTextAlign(22);
                lab->SetTextSize(0.03);
                lab->Draw();
            }
        }
    }
    cCal->Update();

    // Also draw summed calibrated histogram on its own canvas for quick inspection
    TCanvas* cSum = new TCanvas("cSum", "Summed calibrated spectrum", 1000, 700);
    hsum->Draw();
    gPad->SetLogy();
    double ytop_sum = hsum->GetMaximum();
    if (ytop_sum <= 0) ytop_sum = 1.0;
    for (int t = 0; t < 3; ++t) {
        double E = targetE[t];
        if (E >= calMin && E <= calMax) {
            TLine *lt = new TLine(E, 0.5, E, ytop_sum * 1.05);
            lt->SetLineColor(targetColor[t]);
            lt->SetLineWidth(2);
            lt->Draw();
            TLatex *lab = new TLatex(E, ytop_sum * 1.05, Form("%.0f keV", E));
            lab->SetTextAlign(22);
            lab->SetTextSize(0.03);
            lab->Draw();
        }
    }
    cSum->Update();

    // update visuals
    cHist->Update();
    cFit->Update();

    cout << "Saved corrected file: " << (string)outFile << endl;
    cout << "New branch: clover_cross_c.amplitude (calibrated energies in keV) written." << endl;
}

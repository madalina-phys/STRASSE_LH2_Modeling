// Standard and ROOT libraries
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TText.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm> // For std::min_element and std::max_element

// Including the .C file makes the CryoSolver function available.
// It also transitively includes Definitions.h and Functions.h.
#include "CryoSolver.C"

void Plot_Gradients()
{
    // --- Define a single operating point ---
    double p_cond_run = 101325;//21553;//101325; // ~1 atm
    double q_target_run = 0.005*138;  // W/m^2
    double q_th_run = 0.00*138;      // W/m^2

    std::cout << "\n======================================================\n";
    std::cout << "Running CryoSolver for a single point:\n";
    std::cout << "P_cond = " << p_cond_run << " Pa, q_target = " << q_target_run << " W/m^2, q_th = " << q_th_run << " W/m^2" << std::endl;
    std::cout << "======================================================\n";

    // Call the solver. The debug flag is set to false to keep the output clean.
    // The final state is printed by CryoSolver at the end of the run.
    CryoSolver(p_cond_run, q_target_run, q_th_run, false);
    // After CryoSolver runs, the global variables (P1, T1, etc.) are updated.
    
    // --- Define the path length for each point in the loop ---
    std::vector<double> s_path;
    double current_path = 0.0;

    s_path.push_back(current_path); // 0: Condenser Out
    current_path += y_supply;
    s_path.push_back(current_path); // 1: P1, T1 (Bottom of supply)
    current_path += x_supply;
    s_path.push_back(current_path); // 2: P2, T2 (End of horiz. supply)
    current_path += L_th;
    s_path.push_back(current_path); // 3: P2_th, T2_th (Target inlet)
    current_path += x; // x is the filling level
    s_path.push_back(current_path); // 4: P_x, T_x (L/V interface)
    current_path = s_path[3] + D_target; // Path to top of target
    s_path.push_back(current_path); // 5: P3_th, T3_th (Target outlet)
    current_path += L_th;
    s_path.push_back(current_path); // 6: P3, T3 (Start of return line)
    current_path += x_return;
    s_path.push_back(current_path); // 7: P4, T4 (Bottom of return)
    current_path += y_return;
    s_path.push_back(current_path); // 8: P5, T5 (Top of return)
    current_path += dx_supply_return;
    s_path.push_back(current_path); // 9: P6, T6 (Condenser In)

    // --- Collect Pressure and Temperature data ---
    std::vector<double> pressures = {
        ::P_cond, P1, P2, P2_th, P_x, P3_th, P3, P4, P5, P6
    };
    std::vector<double> temperatures = {
        ::T_cond, T1, T2, T2_th, T_x, T3_th, T3, T4, T5, T6
    };
    std::vector<const char*> point_labels = {
        "Cond. Out", "Supply Bot.", "Supply End", "Targ. In", "L/V Iface", "Targ. Out", "Return Start", "Return Bot.", "Return Top", "Cond. In"
    };

    // --- Create TGraph objects ---
    TGraph *gr_pressure = new TGraph(s_path.size(), &s_path[0], &pressures[0]);
    gr_pressure->SetName("gr_pressure");
    gr_pressure->SetTitle("Pressure Profile along the Loop;Path Length (m);Pressure (Pa)");
    gr_pressure->SetMarkerStyle(20);
    gr_pressure->SetMarkerSize(1.0);
    gr_pressure->SetLineWidth(2);
    gr_pressure->SetLineColor(kBlue);
    gr_pressure->SetMarkerColor(kBlue);

    TGraph *gr_temperature = new TGraph(s_path.size(), &s_path[0], &temperatures[0]);
    gr_temperature->SetName("gr_temperature");
    gr_temperature->SetTitle("Temperature Profile along the Loop;Path Length (m);Temperature (K)");
    gr_temperature->SetMarkerStyle(21);
    gr_temperature->SetMarkerSize(1.0);
    gr_temperature->SetLineWidth(2);
    gr_temperature->SetLineColor(kRed);
    gr_temperature->SetMarkerColor(kRed);

    // --- Create canvas and draw ---
    gStyle->SetOptTitle(1);
    TCanvas *c1 = new TCanvas("c_gradients", "Pressure and Temperature Gradients", 1400, 700);
    c1->Divide(2, 1);

    // Draw Pressure Plot
    c1->cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);

    // --- Manually set Y-axis range for pressure plot to add padding ---
    auto p_min_it = std::min_element(pressures.begin(), pressures.end());
    auto p_max_it = std::max_element(pressures.begin(), pressures.end());
    double p_min_val = *p_min_it;
    double p_max_val = *p_max_it;
    double p_range = p_max_val - p_min_val;
    // Add 10% padding, handle case where range is zero
    double p_padding = (p_range > 1e-9) ? 0.1 * p_range : 0.1 * p_max_val;
    gr_pressure->GetYaxis()->SetRangeUser(p_min_val - p_padding, p_max_val + p_padding);

    gr_pressure->Draw("ALP");
    gr_pressure->GetYaxis()->SetTitleOffset(1.5); // Give more space for axis title
    
    // Add labels to pressure plot
    // Get plot range to calculate a dynamic offset for labels
    double y_min_p = gr_pressure->GetYaxis()->GetXmin();
    double y_max_p = gr_pressure->GetYaxis()->GetXmax();
    double offset_p = (y_max_p - y_min_p) * 0.04; // 4% of the y-axis range

    TText *text_p = new TText();
    text_p->SetTextSize(0.025);
    for(size_t i = 0; i < s_path.size(); ++i) {
        // Alternate labels above and below the point to avoid overlap
        if (i % 2 == 0) {
            // Place label below the point
            text_p->SetTextAlign(23); // Horizontally centered, vertically top-aligned
            text_p->DrawText(s_path[i], pressures[i] - offset_p, point_labels[i]);
        } else {
            // Place label above the point
            text_p->SetTextAlign(21); // Horizontally centered, vertically bottom-aligned
            text_p->DrawText(s_path[i], pressures[i] + offset_p, point_labels[i]);
        }
    }

    // Draw Temperature Plot
    c1->cd(2);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);

    // --- Manually set Y-axis range for temperature plot to add padding ---
    auto t_min_it = std::min_element(temperatures.begin(), temperatures.end());
    auto t_max_it = std::max_element(temperatures.begin(), temperatures.end());
    double t_min_val = *t_min_it;
    double t_max_val = *t_max_it;
    double t_range = t_max_val - t_min_val;
    // Add 10% padding, handle case where range is zero
    double t_padding = (t_range > 1e-9) ? 0.1 * t_range : 0.1 * t_max_val;
    gr_temperature->GetYaxis()->SetRangeUser(t_min_val - t_padding, t_max_val + t_padding);

    gr_temperature->Draw("ALP");
    gr_temperature->GetYaxis()->SetTitleOffset(1.5); // Give more space for axis title

    // Add labels to temperature plot
    double y_min_t = gr_temperature->GetYaxis()->GetXmin();
    double y_max_t = gr_temperature->GetYaxis()->GetXmax();
    double offset_t = (y_max_t - y_min_t) * 0.04; // 4% of the y-axis range

    TText *text_t = new TText();
    text_t->SetTextSize(0.025);
    for(size_t i = 0; i < s_path.size(); ++i) {
        // Alternate labels above and below the point to avoid overlap
        if (i % 2 == 0) {
            text_t->SetTextAlign(23); // Horizontally centered, vertically top-aligned
            text_t->DrawText(s_path[i], temperatures[i] - offset_t, point_labels[i]);
        } else {
            text_t->SetTextAlign(21); // Horizontally centered, vertically bottom-aligned
            text_t->DrawText(s_path[i], temperatures[i] + offset_t, point_labels[i]);
        }
    }

    c1->Update();

    std::cout << "\n\nFinished. Gradient plots are now displayed." << std::endl;
}
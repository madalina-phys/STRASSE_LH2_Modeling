// Standard and ROOT libraries
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include <iostream>

// Including the .C file makes the CryoSolver function available.
// It also transitively includes Definitions.h and Functions.h.
#include "CryoSolver.C"

void Plot_Dependencies()
{
    // --- Define parameter ranges for the study ---
    // Condenser Pressure Range
    double p_min = 2.5e4; // 0.3 bar
    double p_max = 10.5e4; // 1 bar
    int p_steps = 40;

    // Target Heat Load Range
    double q_min = 0.0;  // 10 W/m^2
    double q_max = 50.0; // 200 W/m^2
    int q_steps = 10;

    // --- Create TGraph2D objects to store the results ---
    TGraph2D *gr_filling_level = new TGraph2D();
    gr_filling_level->SetName("gr_filling_level");
    gr_filling_level->SetTitle("Filling Level (%) vs. P_cond (Pa) and q_target (W/m^2);P_cond (Pa);q_target (W/m^2);Filling Level (m)");

    TGraph2D *gr_mass_flow = new TGraph2D();
    gr_mass_flow->SetName("gr_mass_flow");
    gr_mass_flow->SetTitle("Mass Flow (kg/s) vs. P_cond (Pa) and q_target (W/m^2);P_cond (Pa);q_target (W/m^2);Mass Flow (kg/s)");

    TGraph2D *gr_P3 = new TGraph2D();
    gr_P3->SetName("gr_P3");
    gr_P3->SetTitle("P3 (Pa) vs. P_cond (Pa) and q_target (W/m^2);P_cond (Pa);q_target (W/m^2);P3 (Pa)");

    // --- Loop over parameters and run the solver ---
    int point_index = 0;
    for (int i = 0; i <= p_steps; ++i)
    {
        double p_current = p_min + (p_max - p_min) * i / p_steps;
        for (int j = 0; j <= q_steps; ++j)
        {
            double q_current = q_min + (q_max - q_min) * j / q_steps;

            std::cout << "\n======================================================\n";
            std::cout << "Running point " << point_index + 1 << "/" << (p_steps + 1) * (q_steps + 1) << ": ";
            std::cout << "P_cond = " << p_current << " Pa, q_target = " << q_current << " W/m^2" << std::endl;
            std::cout << "======================================================\n";

            // Call the solver. Assume the thermal holder heat load (q_th) is 0 for this study.
            // The debug flag is set to false to keep the output clean.
            CryoSolver(p_current, q_current, 5, false);

            // Store results. y_level_target and mass_flow_rate are global variables updated by CryoSolver.
            gr_filling_level->SetPoint(point_index, p_current, q_current, y_level_target/D_target*100); // Store as percentage of target diameter
            gr_mass_flow->SetPoint(point_index, p_current, q_current, mass_flow_rate);
            gr_P3->SetPoint(point_index, p_current, q_current, mass_flow_rate*(Cv_l*(T_x-T_cond)+L_v));
            point_index++;
        }
    }

    // --- Create canvases and draw the graphs ---
    gStyle->SetPalette(kRainBow);
    TCanvas *c1 = new TCanvas("c_all", "CryoSolver Dependencies", 1800, 700);
    c1->Divide(3, 1);

    c1->cd(1)->SetRightMargin(0.15);
    gr_filling_level->Draw("colz");
    c1->Update();

    c1->cd(2)->SetRightMargin(0.15);
    gr_mass_flow->Draw("colz");

    c1->Update();
    c1->cd(3)->SetRightMargin(0.15);
    gr_P3->Draw("colz");
    c1->Update();
    // --- Save results to a ROOT file ---
    TFile *outFile = new TFile("CryoSolver_Dependencies.root", "RECREATE"); 
    gr_filling_level->Write();
    gr_mass_flow->Write();
    outFile->Close();

    std::cout << "\n\nFinished. Plots are displayed and results are saved to CryoSolver_Dependencies.root" << std::endl;
}
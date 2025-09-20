#include "Definitions.h"
#include "Functions.h"

void CryoSolver(double P_cond, double q_target, double q_th, bool debug = false)
{
    T_cond = SaturationTemperature(P_cond); //update condensation temperature based on condensation pressure
    A_target = L_target * D_target * TMath::Pi(); //total area surrounding the target [m^2]

    //Initial guess
    T1 = T_cond; //initial guess for the temperature at the bottom of the supply [K]
    T2 = T_cond; //initial guess for the temperature at the end of the horizontal supply [K]
    T2_th = T_cond; //initial guess for the temperature at the target inlet [K]
    T_x = T_cond; //initial guess for the temperature at the liq-vap interface in the target [K]
    T3_th = T_cond; //initial guess for the temperature at the target outlet [K]
    T3 = T_cond; //initial guess for the temperature at the top of the target [K]
    T4 = T_cond; //initial guess for the temperature at the bottom of the return line [K]
    T5 = T_cond; //initial guess for the temperature at the top of the return line [K]
    T6 = T_cond; //initial guess for the temperature at the condenser outlet [K]
    P_x = P_cond; //initial guess for the pressure at the liq-vap interface in the target [Pa]
    m_dot_old = m_dot; //initial guess for the vapor mass flow rate [kg/s]
    x_old = x; //initial guess for the filling level in the target
    double error = 1, error_m_dot, error_x; //initial error
    int iter = 0; //iteration counter
    cout<<"Starting CryoSolver with P_cond: "<<P_cond<<" Pa, T_cond: "<<T_cond<<" K, q_target: "<<q_target<<" W/m^2, q_th: "<<q_th<<" W/m^2"<<endl;


    while ((error > 1E-6 && iter < 100000) || iter < 4)
    {
        if (debug) cout<<endl<<"--- Iteration "<<iter<<" ---"<<endl<<endl;

        // fluid properties as function of T and P
        
        rho_l = DensityLiquidHydrogen(T_cond); //update liquid density based on condensation temperature
        rho_v = DensityVaporHydrogen(T3); //update vapor density based on condensation temperature
        rho_v2 = DensityVaporHydrogen(T_cond); //update vapor density based on condensation temperature
        mu_l = ViscosityLiquidHydrogen(T_cond); //update liquid viscosity based on condensation temperature
        mu_v = ViscosityVaporHydrogen(T3); //update vapor viscosity based on condensation temperature
        mu_v2 = ViscosityVaporHydrogen(T_cond);
        k_l = ThermalConductivityLiquidHydrogen(T_cond); //update liquid thermal conductivity based on condensation temperature
        k_v = ThermalConductivityVaporHydrogen(T3); //update vapor thermal conductivity based on condensation temperature
        Cv_l = SpecificHeatLiquidHydrogen(T_cond); //update liquid specific heat capacity based on condensation temperature
        Cv_v = SpecificHeatVaporHydrogen(T3); //update vapor specific heat capacity based on condensation temperature
        L_v = LatentHeatHydrogen(P_x); //update latent heat based on pressure at the liq-vap interface in the target
        
        //Calculate friction factors
        friction_factor_l = FrictionFactor(rho_l, mu_l, D_pipe, m_dot, vel_l);
        friction_factor_v = FrictionFactor(rho_v, mu_v, D_pipe, m_dot, vel_v);
        friction_factor_v2 = FrictionFactor(rho_v2, mu_v2, D_pipe, m_dot, vel_v2);

        CalculatePressuresAndTemperatures(P_cond, T_cond, q_target, q_target_front, q_th, m_dot, x, y_level_cond,
                                       P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th,
                                       P3, T3, P4, T4, P5, T5, P6, T6,
                                       vel_l, vel_v,
                                       debug);
        if (debug) PrintIterationDebugInfo(iter, P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th, P3, T3, P4, T4, P5, T5, P6, T6, x, y, m_dot, vel_l, vel_v, error,  error_m_dot, error_x);
        //Calculate new filling level in the target
        x = FillingLevel(q_target, q_target_front, q_th, m_dot, 
                        rho_l, rho_v, 
                        T2_th, P2_th, T_x, P_x,
                        debug);
        if (x>D_target)x = D_target;
        //Calculate new condenser level based on filling level in the target
        y = LevelCondenser(x);
        if (debug) cout<<"\t[LevelCondenser] Calculated condenser level: "<<y<<" m"<<endl;
        //Update y_level_cond global variable
        y_level_cond = y;                                



        //Calculate new vapor mass flow rate
        m_dot = VaporMassFlowRate_Bisection(P_cond, T_cond, q_target, q_target_front, q_th, m_dot_old, x, y_level_cond,
                                       P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th,
                                       P3, T3, P4, T4, P5, T5, P6, T6,
                                       vel_l, vel_v,
                                       debug);
        
        CalculatePressuresAndTemperatures(P_cond, T_cond, q_target, q_target_front, q_th, m_dot, x, y_level_cond,
                                       P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th,
                                       P3, T3, P4, T4, P5, T5, P6, T6,
                                       vel_l, vel_v,
                                       debug);                  
        //Calculate error
        error_m_dot = fabs((m_dot - m_dot_old) / m_dot_old);
        error_x = fabs((x - x_old) / (x_old+0.01)); 
        error = error_m_dot + error_x;
        
        if (debug) PrintIterationDebugInfo(iter, P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th, P3, T3, P4, T4, P5, T5, P6, T6, x, y, m_dot, vel_l, vel_v, error,  error_m_dot, error_x);

        //Update old values
        m_dot_old = m_dot;
        x_old = x;
        iter++;
    }
    y_level_target = x; //update the global variable for the filling level in the target
    y_level_cond = y;
    mass_flow_rate = m_dot; //update the global variable for the vapor mass flow rate
    PrintIterationDebugInfo(iter, P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th, P3, T3, P4, T4, P5, T5, P6, T6, x, y, m_dot, vel_l, vel_v, error, error_m_dot, error_x);

}

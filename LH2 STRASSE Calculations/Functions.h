

double SaturationTemperature(double P){
    //parameters for LH2
    double A = 0.34317799;
    double B = -4.59895759;
    double C = 27.69524633;
    double T_sat = A*pow(log(P),2) + B*log(P) + C; //in K
    return T_sat;
}

double SaturationPressure(double T){
    //parameters for LH2
    // Coefficients for polynomial: p5*T^5 + p4*T^4 + ... + p0
    double p5 = 7.24480785e-03;
    double p4 = -9.50206160e-02;
    double p3 = 1.08001326e+02;
    double p2 = -3.73158679e+03;
    double p1 = 4.52369927e+04;
    double p0 = -1.90683817e+05;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double P_sat = (((((p5 * T + p4) * T + p3) * T + p2) * T + p1) * T + p0); //in Pa
    return P_sat;
}

double DensityLiquidHydrogen(double T){ 
    //parameters for LH2
    // Coefficients for 3rd degree polynomial: p3*T^3 + p2*T^2 + p1*T + p0
    // -2.44398187e-03  1.10790064e-01 -2.58840598e+00  9.81502143e+01 poly3
    double p3 = -2.44398187e-03;
    double p2 =  1.10790064e-01;
    double p1 = -2.58840598e+00;
    double p0 =  9.81502143e+01;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double rho_l = ((p3 * T + p2) * T + p1) * T + p0; //in kg/m^3
    return rho_l;
}

double DensityVaporHydrogen(double T){ 
    //parameters for LH2 
    // 1.66393764e-04 -1.18030170e-02  3.38343049e-01 -4.35410388e+00  2.08082954e+01   poly4
    double p4 =  1.66393764e-04;
    double p3 = -1.18030170e-02;
    double p2 =  3.38343049e-01;
    double p1 = -4.35410388e+00;
    double p0 =  2.08082954e+01;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double rho_v = ((((p4 * T + p3) * T + p2) * T + p1) * T + p0); //in kg/m^3
    return rho_v;
}

double ViscosityLiquidHydrogen(double T){       
    //parameters for LH2
    // -4.29608892e-09  3.43112703e-07 -9.75692114e-06  1.06001866e-04 poly3
    double p3 = -4.29608892e-09;
    double p2 =  3.43112703e-07;
    double p1 = -9.75692114e-06;
    double p0 =  1.06001866e-04;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double mu_l = ((p3 * T + p2) * T + p1) * T + p0; //in Pa*s
    return mu_l;
}

double ViscosityVaporHydrogen(double T){       
    //parameters for LH2
    // 1.05984092e-11 -7.82973846e-10  2.17374955e-08 -2.13433478e-07  1.11788367e-06 poly4
    double p4 =  1.05984092e-11;
    double p3 = -7.82973846e-10;
    double p2 =  2.17374955e-08;
    double p1 = -2.13433478e-07;
    double p0 =  1.11788367e-06;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double mu_v = ((((p4 * T + p3) * T + p2) * T + p1) * T + p0); //in Pa*s
    return mu_v;
}

double ThermalConductivityLiquidHydrogen(double T){ 
    // 1.45240215e-06 -2.89630453e-04  1.06463472e-02 -8.27115517e-03 poly
    double p3 = 1.45240215e-06;
    double p2 = -2.89630453e-04;
    double p1 =  1.06463472e-02;
    double p0 = -8.27115517e-03;   
    double k_l = (((p3 * T + p2) * T + p1) * T + p0); //in W/m*K
    return k_l;
}

double  ThermalConductivityVaporHydrogen(double T){ 
    //parameters for LH2
    // 4.65642861e-07 -3.60230977e-05  1.06085737e-03 -1.29629405e-02 poly
    double p3 =  4.65642861e-07;
    double p2 = -3.60230977e-05;
    double p1 =  1.06085737e-03;
    double p0 = -1.29629405e-02;
    double k_v = (((p3 * T + p2) * T + p1) * T + p0); //in W/m*K
    return k_v;

}

double SpecificHeatLiquidHydrogen(double T){       
    //parameters for LH2
    // 9.98115256e-02 -1.00173729e+01  3.99420289e+02 -7.87179864e+03   7.68975462e+04 -2.91664049e+05 poly5
    double p5 =  9.98115256e-02;
    double p4 = -1.00173729e+01;
    double p3 =  3.99420289e+02;
    double p2 = -7.87179864e+03;
    double p1 =  7.68975462e+04;
    double p0 = -2.91664049e+05;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double Cv_l = (((((p5 * T + p4) * T + p3) * T + p2) * T + p1) * T + p0); //in J/K*kg
    return Cv_l;
}

double SpecificHeatVaporHydrogen(double T){       
    //parameters for LH2
    // 1.48484086e-01 -1.49676287e+01  6.00385831e+02 -1.19278080e+04  1.17274475e+05 -4.45892740e+05 poly5
    double p5 =  1.48484086e-01;
    double p4 = -1.49676287e+01;
    double p3 =  6.00385831e+02;
    double p2 = -1.19278080e+04;
    double p1 =  1.17274475e+05;
    double p0 = -4.45892740e+05;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double Cv_v = (((((p5 * T + p4) * T + p3) * T + p2) * T + p1) * T + p0); //in J/K*kg
    return Cv_v;
}

double LatentHeatHydrogen(double P){       
    //parameters for LH2
    // 2.23767861e-21 -9.78934270e-16  1.67304839e-10 -1.44046249e-05  5.09861609e-01  4.47589140e+05 poly5
    double p5 =  2.23767861e-21;
    double p4 = -9.78934270e-16;
    double p3 =  1.67304839e-10;
    double p2 = -1.44046249e-05;
    double p1 =  5.09861609e-01;
    double p0 =  4.47589140e+05;
    // Using Horner's method for efficient and accurate polynomial evaluation
    double L_v = (((((p5 * P + p4) * P + p3) * P + p2) * P + p1) * P + p0); //in J/Kg
    return L_v;
}

double FrictionFactor(double rho, double mu, double D, double m_dot, double& v){
    //Calculate velocity
    double A = TMath::Pi()/4 * D*D; //cross-sectional area [m^2]
    v = m_dot/(rho*A); //velocity [m/s]

    //Calculate Reynolds number
    double Re = rho*v*D/mu;

    //Calculate friction factor using the Blasius correlation for turbulent flow
    double f;
    if (Re < 2300) {
        f = 64/Re; //laminar flow
        //cout<<"\t\tLaminar flow detected (Re = "<<Re<<"). Friction factor calculated using Hagen-Poiseuille equation."<<endl;
    } else {
        f = 0.3164*pow(Re,-0.25); //turbulent flow
        //cout<<"\t\tWarning: Turbulent flow detected (Re = "<<Re<<"). Friction factor calculated using Blasius correlation."<<endl;
    }

    f = f / D * rho * pow(v,2) / 2; //convert to pressure drop per unit length [Pa/m]
    return f;
}

double CalculateBoilingPressure(double m_dot, double rho_l_temp, double rho_v_temp, double P2_th_temp, double T2_th_temp, double T3_th_temp){ 
    if (m_dot <= 0) {
        return 0.0;
    }

    // 2. Calculate Mass Flux (G)
    double D_tmp = D_target/2; // Use the pipe diameter for mass flux calculation
    double G = m_dot / (TMath::Pi() * D_tmp * D_tmp/4);

    double dP_a = pow(G, 2) * (1/rho_v_temp - 1/rho_l_temp) + P2_th_temp*factor_boiling_pressure;
    //cout<<T3_th_temp<<".  "<<T2_th_temp<<endl;
    return dP_a;
}

void CalculatePressuresAndTemperatures(double P_cond, double T_cond, double q_tmp, double q_front_tmp, double q_th_tmp, double m, double x, double y_level_cond,
                                       double& P1, double& T1,
                                       double& P2, double& T2,
                                       double& P2_th, double& T2_th,
                                       double& P_x, double& T_x,
                                       double& P3_th, double& T3_th,
                                       double& P3, double& T3,
                                       double& P4, double& T4,
                                       double& P5, double& T5,
                                       double& P6, double& T6,
                                       double& vel_l, double& vel_v,
                                       bool debug){ 
    
    if (m<=0){
        friction_factor_l = 0;
        friction_factor_v = 0;
        friction_factor_v2 = 0;
        vel_l = 0;
        vel_v = 0;
        vel_v2 = 0;
    }
    else{
    friction_factor_l = FrictionFactor(rho_l, mu_l, D_pipe, m, vel_l);
    friction_factor_v = FrictionFactor(rho_v, mu_v, D_pipe, m, vel_v);
    friction_factor_v2 = FrictionFactor(rho_v2, mu_v2, D_pipe, m, vel_v2);
    }
    // The wetted surface area of a horizontal cylinder is L * D * acos(1 - 2x/D)
    // where x is the fill height from the bottom.
    A_target_l = L_target * D_target * TMath::ACos(1.0 - 2.0 * x / D_target);
    A_target_v = (TMath::Pi() * D_target * L_target) - A_target_l;

    A_target_front_l = D_target * D_target * TMath::ACos(1.0 - 2.0 * x / D_target) /4 - (D_target/2 -x)*D_target/2*TMath::Sin(1.0 - 2.0 * x / D_target);
    A_target_front_v = (TMath::Pi() * D_target * D_target / 4) - A_target_front_l;


    // State 1: Supply Bottom
    P1 = P_cond + rho_l * g * y_level_cond - friction_factor_l * y_supply; //pressure at the condenser outlet [Pa]
    T1 = T_cond;
    //T1 = (T_cond / y_supply + T2 / x_supply) / (1/y_supply + 1/x_supply); 

    // State 2: Target Inlet
    P2 = P1 - friction_factor_l * x_supply; //pressure at the target inlet [Pa]
    T2 = T1;
    //T2 = (T1 / x_supply + T2_th / L_th) / (1/x_supply + 1/L_th); 

    // State 2_th: Target Inlet (after target holder)
    P2_th = P2 + rho_l/2 * pow(vel_l,2) * pow((1 - pow(D_pipe/D_target*2,2)),2)/2;;
    T2_th = T2 + q_th_tmp * A_th / (m * Cv_l); //temperature at the end of the target holder pipe [K]          
    
    T_x = T2_th + (q_tmp * A_target_l + 2*q_front_tmp * A_target_front_l) / (m * Cv_l + m * vap_percentage * L_v);
    T3_th = T_x + (q_tmp * A_target_v + 2*q_front_tmp * A_target_front_v) / (m * Cv_v + m * (1-vap_percentage) * L_v); //temperature at the target outlet [K]

    double dP_boiling = CalculateBoilingPressure(m, rho_l, rho_v, P2_th, T2_th, T3_th);
    P_x = P2_th - rho_l * g * x + dP_boiling; //pressure at the liq-vap interface in the target [Pa]
    P3_th = P_x + rho_v * g * (D_target - x); //pressure at the target outlet [Pa]

    P3 = P3_th - rho_v/2 * pow(vel_v,2) * pow((1 - pow(D_pipe/D_target*2,2)),2)/2; //pressure at the target outlet [Pa]
    T3 = T3_th + q_th_tmp * A_th / (m * Cv_v); //temperature at the end of the target holder pipe [K]

    P4 = P3 - friction_factor_v * x_return; //pressure at the bottom of the return line [Pa]
    T4 = T3;
    //T4 = (T3 / x_return + T5 / y_return) / (1/x_return + 1/y_return); 

    P5 = P4 - friction_factor_v2 * y_return + rho_v * g * y_return; //pressure at the top of the return line [Pa]
    //T5 = T4;
    T5 = (T4*2 / y_return + T6 / dx_supply_return/2) / (2/y_return + 1/dx_supply_return/2);

    P6 = P5 - friction_factor_v2 * 2*dx_supply_return;// + rho_v/2 * pow(vel_v,2) * pow((1 - pow(D_pipe/D_condenser,2)),2)/2; //pressure at the condenser [Pa]
    T6 = T_cond;
}


double FillingLevel(double q_tmp, double q_front_tmp, double q_th_tmp, double m_dot, 
                    double rho_l, double rho_v, 
                    double T2_th_temp, double P2_th_temp, double& T_x, double& P_x,
                    bool debug){

    //Calculate filling level in the target
    double x = 0;
    while (x <= D_target) {

        // The wetted surface area of a horizontal cylinder is L * D * acos(1 - 2x/D)
        // where x is the fill height from the bottom.
        double A_l = L_target * D_target * TMath::ACos(1.0 - 2.0 * x / D_target);
        double A_v = (TMath::Pi() * D_target * L_target) - A_target_l;

        double A_front_l = D_target * D_target * TMath::ACos(1.0 - 2.0 * x / D_target) /4  - (D_target/2 -x)*D_target/2*TMath::Sin(1.0 - 2.0 * x / D_target);
        double A_front_v = (TMath::Pi() * D_target * D_target / 4) - A_front_l;

        T_x = T2_th_temp + (q_tmp * A_l + 2*q_front_tmp * A_front_l) / (m_dot * Cv_l + m_dot * vap_percentage * L_v); //temperature at the liq-vap interface in the target [K]
        double T3_th_temp = T_x + (q_tmp * A_v + 2*q_front_tmp * A_front_v) / (m_dot * Cv_v + m_dot * (1-vap_percentage) * L_v);

        dP_boiling = CalculateBoilingPressure(m_dot, rho_l, rho_v, P2_th_temp, T2_th_temp, T3_th_temp);
        P_x = P2 - rho_l * g * x + dP_boiling; //pressure at the liq-vap interface in the target [Pa]
        
        T_sat = SaturationTemperature(P_x);
        if (debug) cout<<"\t\tTrying filling level: "<<x<<" m, A_l: "<<A_l<<" m^2, P_x: "<<P_x<<" Pa, T_sat: "<<T_sat<<" K, T_x: "<<T_x<<" K"<<endl;
        if (T_x >= T_sat) {
            if (debug) cout<<"\t[FillingLevel] breaking at filling level: "<<x<<" m, A_l: "<<A_l<<" m^2, P_x: "<<P_x<<" Pa, T_sat: "<<T_sat<<" K, T_x: "<<T_x<<" K"<<endl;
            break;
        }
        else x += 0.00001; //increment filling level [m]
    }
    if (debug) cout<<"\t[FillingLevel] Calculated filling level: "<<x<<" m, P_x: "<<P_x<<" Pa, T_sat: "<<T_sat<<" K"<<endl;
    T_x = T_sat; //set T_x to the saturation temperature at the pressure P_x
    return x;
}

double LevelCondenser(double x){
    //Calculate the y level of the condenser based on the filling level in the target
    double y = y_supply + x * H_condenser / D_target/2+dMass; //y level of the top of the condenser [m]
    return y;
}

void PrintIterationDebugInfo(int iter,
                             double P1, double T1,
                             double P2, double T2,
                             double P2_th, double T2_th,
                             double P_x, double T_x,
                             double P3_th, double T3_th,
                             double P3, double T3,
                             double P4, double T4,
                             double P5, double T5,
                             double P6, double T6,
                             double x, double y, double m_dot,
                             double vel_l, double vel_v,
                             double error, double error_m_dot, double error_x)
{
    cout << "\tState Points for iteration: " << iter  << endl;
    cout << "\t  P1: " << P1 << " Pa, T1: " << T1 << " K" << endl;
    cout << "\t  P2: " << P2 << " Pa, T2: " << T2 << " K" << endl;
    cout << "\t  P2_th: " << P2_th << " Pa, T2_th: " << T2_th << " K (Target Inlet)" << endl;
    cout << "\t  P_x: " << P_x << " Pa, T_x: " << T_x << " K (L/V Interface)" << endl;
    cout << "\t  P3_th: " << P3_th << " Pa, T3_th: " << T3_th << " K (Target Outlet)" << endl;
    cout << "\t  P3: " << P3 << " Pa, T3: " << T3 << " K" << endl;
    cout << "\t  P4: " << P4 << " Pa, T4: " << T4 << " K" << endl;
    cout << "\t  P5: " << P5 << " Pa, T5: " << T5 << " K" << endl;
    cout << "\t  P6: " << P6 << " Pa, T6: " << T6 << " K" << endl;
    cout << "\tIteration Results:" << endl;
    cout << "\t  Filling Level (x):  ------------" << x << " m / " << x/D_target*100 << " % ------------" << endl;
    cout << "\t  Y Level Cond: "<<y<<" m"<<endl;
    cout << "\t  Vapor Mass Flow (m_dot): " << m_dot << " kg/s" << endl;
    cout << "\t  Liquid Velocity (vel_l): " << vel_l << " m/s" << endl;
    cout << "\t  Vapor Velocity (vel_v): " << vel_v << " m/s" << endl;
    cout << "\tConvergence:" << endl;
    cout << "\t  Total Error: " << error << endl;
    cout << "\t    - Mass Flow Error:   " << error_m_dot << endl;
    cout << "\t    - Filling Level Error: " << error_x << endl;
}


double VaporMassFlowRate_Bisection(double P_cond, double T_cond, double q_tmp, double q_front_tmp, double q_th_tmp, double m, double x, double y_level_cond,
                                       double& P1, double& T1,
                                       double& P2, double& T2,
                                       double& P2_th, double& T2_th,
                                       double& P_x, double& T_x,
                                       double& P3_th, double& T3_th,
                                       double& P3, double& T3,
                                       double& P4, double& T4,
                                       double& P5, double& T5,
                                       double& P6, double& T6,
                                       double& vel_l, double& vel_v,
                                       bool debug){
    
    double theta = TMath::Pi() - atan((D_target/2 - x)/(D_target/2)); //angle for the liquid phase
    double d = abs(sin(theta)*D_target); 
    double A_interface = 2*d*L_target; //cross-sectional area of the liquid phase [m^2], assuming 1 m length for pressure drop calculation
    // Bisection method variables
    double m_low = 0.0;
    double m_high = 0.01; // A reasonable upper bound, e.g., 100 g/s. Adjust if necessary.
    double m_mid;
    const int MAX_ITER = 1000;
    const double TOLERANCE = 1e-6;

    if (debug) {
        std::cout << "\t[VaporMassFlowRate] Starting bisection search for m_dot..." << std::endl;
    }

    CalculatePressuresAndTemperatures(P_cond, T_cond, q_tmp, q_front_tmp, q_th_tmp, m_low, x, y_level_cond,
                                       P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th,
                                       P3, T3, P4, T4, P5, T5, P6, T6,
                                       vel_l, vel_v,
                                       debug);
    
    double err_low = P6-P_cond;
    //cout<< "err_low: "<<err_low<<endl;
    // Check that our bounds actually bracket the solution
    
    CalculatePressuresAndTemperatures(P_cond, T_cond, q_tmp, q_front_tmp, q_th_tmp, m_high, x, y_level_cond,
                                       P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th,
                                       P3, T3, P4, T4, P5, T5, P6, T6,
                                       vel_l, vel_v,
                                       debug);
    double err_high = P6-P_cond;
    //cout<< "err_high: "<<err_high<<endl;

        if (err_low * err_high > 0) {
        std::cerr << "\tERROR in VaporMassFlowRate: Root is not bracketed by initial guesses." << std::endl;
        std::cerr << "\tError at m_low=0: " << err_low << ", Error at m_high=" << m_high << ": " << err_high << std::endl;
        return m_low; // Return a safe value
        }

        for (int i = 0; i < MAX_ITER; ++i) {
        m_mid = (m_low + m_high) / 2.0;
        CalculatePressuresAndTemperatures(P_cond, T_cond, q_tmp, q_front_tmp, q_th_tmp, m_mid, x, y_level_cond,
                                       P1, T1, P2, T2, P2_th, T2_th, P_x, T_x, P3_th, T3_th,
                                       P3, T3, P4, T4, P5, T5, P6, T6,
                                       vel_l, vel_v,
                                       debug);
        double err_mid = P6-P_cond;
        
        if (std::abs(err_mid) < TOLERANCE) { // <-- Typo: was error_mid
            break; // Convergence
        }

        if (std::signbit(err_mid) == std::signbit(err_low)) { // <-- Typo: was error_mid
             // If error at mid has same sign as error at low, move low to mid
            m_low = m_mid;
        } else {
            m_high = m_mid;
        }
        }

        if (debug) {
        std::cout << "\t[VaporMassFlowRate] Converged to: m_dot = " << m_mid << " kg/s" << std::endl;
        }

        return m_mid;
}

    
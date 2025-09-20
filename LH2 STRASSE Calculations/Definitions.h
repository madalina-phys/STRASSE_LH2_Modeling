double P_cond = 101325; //condensation pressure [Pa]
double T_cond = 20.238; //condensation temperature [K]
double T_amb = 293; //ambient temperature [K]
double q_target = 0.01*138; //heat load on the target [W/m^2]
double q_target_front = 0.5*138; //heat load on the target [W/m^2]
double q_th = 0.01*138; //heat load on the target [W/m^2] 
double D_target = 0.031; //diameter of the target pipe [m]
double L_target = 0.15; //length of the target pipe [m]
double A_target = L_target * D_target * TMath::Pi(); //total area surrounding the target [m^2]
double A_target_l, A_target_v; //area of the liquid and vapor phase in the target [m^2]
double A_target_front_l, A_target_front_v; //area of the liquid and vapor phase in the target [m^2]
double V_target = TMath::Pi() * pow((D_target/2),2) * L_target; //volume of the target [m^3]
double V_target_l, V_target_v; //volume of the liquid and vapor phase in the target [m^3]
double D_pipe = 0.004; //diameter of the supply/return pipe [m]
double D_th = 0.003; //diameter of the pipe connecting the target to the target holder [m]
double L_th = 0.05; //length of the pipe connecting the target to the target holder [m]
double A_th = TMath::Pi() * D_th * L_th;; //area of the target holder pipe [m^2]

double H_condenser = 0.05; //height of the condenser [m]
double D_condenser = 0.076; //diameter of the condenser [m]
double y_supply = 0.32; //y length of the supply line [m] / cond bottom
double y_cond_top = y_supply + H_condenser; //y level of the top of the condenser [m] / diameter condenser = 76mm / height condenser 35mm / height with the conic part 56mm
double y_return = y_cond_top - D_target; //y length of the return line [m]
double x_supply = 0.38; //x length of the supply line [m]
double dx_supply_return = 0.08; //horizontal distance between the supply and the return line [m]
double x_return = x_supply - dx_supply_return; //x length of the return line [m]

double y_level_cond = y_supply + H_condenser/2; //y level of liquid in the condenser [m]
double y_level_target = D_target/2; //y level of the target [m] - to be calculated by the CryoSolver 
double x = y_level_target;
double x_old = x;
double y = y_level_cond;
double y_old = y; 

double rho_l = 70.79; //liquid density [kg/m^3]
double rho_v = 1.34;//gas density [kg/m^3]
double rho_v2 = 1.34;//gas density [kg/m^3]
double mu_l = 13.92E-6; //liquid viscosity [Pa*s] = [kg/m*s]
double mu_v = 10.93E-6; //gas viscosity [Pa*s] = [kg/m*s]
double mu_v2 = 10.93E-6; //gas viscosity [Pa*s] = [kg/m*s]
double k_l = 1e-4; //liquid thermal conductivity [W/m*K]
double k_v = 1e-4; //gas thermal conductivity [W/m*K]
double L_v = 445400; //latent heat [J/Kg]
double Cv_l = 10300; //specific heat capacity [J/K*kg]
double Cv_v = 10300; //specific heat capacity [J/K*kg]
double g = 9.80665; //gravity [m/s^2]
double bool_heat_conductance = 1;

double m_dot = 0.001; //vapor mass flow rate [kg/s] - to be calculated by the CryoSolver
double m_dot_old = m_dot;
double mass_flow_rate = m_dot; //total mass flow rate [kg/s]
double T_sat; //saturation temperature in the target[K] - T_x
double T_x, P_x; //temperature and pressure at the liq-vap interface in the target [K],[Pa]
double friction_factor_l = 0; //friction factor for liquid phase
double friction_factor_v = 0; //friction factor for vapor phase
double friction_factor_v2 = 0; //friction factor for vapor phase
double vel_l, vel_v, vel_v2; //velocity of liquid and vapor phase in the pipe [m/s] l= liq, v=vap, v2 = vap cold
double vap_percentage = 0.1; //percentage of vapor in the mass flow rate [0-1]
double factor_boiling_pressure = 0.00;

double T1,P1; //temperature and pressure at the bottom of the supply line [K],[Pa]
double T2,P2; //temperature and pressure at the end of the supply line [K],[Pa]
double T2_th,P2_th; //temperature and pressure at the end of the target holder pipe [K],[Pa]
double T3,P3; //temperature and pressure at the top of the target [K],[Pa]
double T3_th,P3_th; //temperature and pressure at the top of the target holder pipe [K],[Pa]
double T4,P4; //temperature and pressure at the bottom of the return line [K],[Pa]
double T5,P5; //temperature and pressure at the top of the return line [K],[Pa]
double T6,P6; //temperature and pressure at the condenser [K],[Pa]

double A_interface; //area of the liquid-vapor interface in the target [m^2]
double dP_total, dP_expected; //total and expected pressure drop in the return line [Pa]
double dP_boiling = 0;
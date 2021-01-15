#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters
{

public:
    int lx, ly, ns, IterMax, MCNorm, RandomSeed;
    int n_orbs;
    string ModelType;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus, mus_Cluster, Fill, pi;
    double Total_Particles;
    double K1x, K1y;
    Mat_1_doub J_Hund, OnSiteE;
    double lambda_lattice;
    double k_const;
    double Disorder_Strength, RandomDisorderSeed;
    double Boltzman_constant;
    double BoundaryConnection;
    Mat_1_doub Temp_values;
    bool Read_Seed_from_file_;
    string Seed_file_name_;

    Matrix<double> hopping_NN_X;
    Matrix<double> hopping_NN_Y;
    Matrix<double> hopping_NNN_PXPY;
    Matrix<double> hopping_NNN_PXMY;
    double hopping_intracell;

    Mat_1_string MC_DOF;

    bool Cooling_;
    bool ED_;

    bool Metropolis_Algorithm;
    bool Heat_Bath_Algorithm;

    bool MC_on_theta, MC_on_phi, MC_on_theta_and_phi, MC_on_theta_and_phi_and_u, MC_on_moment_size, MC_on_local_density;

    bool fix_mu;

    double fixed_mu_value;

    /*
SavingMicroscopicStates=1
NoOfMicroscopicStates=50
      */
    bool Saving_Microscopic_States;
    int No_Of_Microscopic_States;

    double temp_max, beta_min;
    double temp_min, beta_max;
    double d_Temp;

    int Last_n_sweeps_for_measurement;
    int Measurement_after_each_m_sweeps;

    double temp, beta, Eav, maxmoment;
    double WindowSize, AccCount[2];
    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);

    double eta;
    double omega_min, omega_max, domega;
};

void Parameters::Initialize(string inputfile_)
{


    maxmoment = 10.0;
    string monte_carlo_variables_;
    double cooling_double;
    double Read_Seed_from_file_double;
    double metropolis_double;
    double Perform_HF_SC_calculation_double;
    double Simple_Mixing_double;
    double ED_double;
    int SavingMicroscopicStates_int;
    int no_of_temp_points;
    string temp_values_;
    string temp_string;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;

    monte_carlo_variables_ = matchstring2(inputfile_, "Monte_Carlo_variables");

    stringstream monte_carlo_variables_stream(monte_carlo_variables_);
    int no_of_MC_vars;
    monte_carlo_variables_stream >> no_of_MC_vars;
    assert(no_of_MC_vars < 6);
    MC_on_theta = false;
    MC_on_phi = false;
    MC_on_theta_and_phi =false;
    MC_on_moment_size = false;
    MC_on_local_density = false;
    MC_DOF.resize(no_of_MC_vars);

    for (int i = 0; i < no_of_MC_vars; i++)
    {
        monte_carlo_variables_stream >> temp_string;
        MC_DOF[i] = temp_string;

        if (temp_string == "phi")
        {
            MC_on_phi = true;
        }
        else if (temp_string == "theta")
        {
            MC_on_theta = true;
        }
        else if (temp_string == "theta_and_phi")
        {
            MC_on_theta_and_phi = true;
        }
        else if (temp_string == "theta_and_phi_and_u")
        {
            MC_on_theta_and_phi_and_u = true;
        }
        else if (temp_string == "moment_size")
        {
            MC_on_moment_size = true;
        }
        else if (temp_string == "local_density")
        {
            MC_on_local_density = true;
        }
        else
        {
            cout << "MC cannot be performed on " << temp_string << ", please choose only from phi theta theta_and_phi theta_and_phi_and_u moment_size local_density" << endl;
            assert(false);
        }
    }

    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));

    TBC_mx = int(matchstring(inputfile_, "TwistedBoundaryCond_mx"));
    n_orbs = int(matchstring(inputfile_, "N_Orbs"));
    J_Hund.resize(n_orbs);
    OnSiteE.resize(n_orbs);

    TBC_my = int(matchstring(inputfile_, "TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_, "TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_, "TBC_cellsY"));
    lx_cluster = int(matchstring(inputfile_, "Cluster_lx"));
    ly_cluster = int(matchstring(inputfile_, "Cluster_ly"));
    SavingMicroscopicStates_int = int(matchstring(inputfile_, "SavingMicroscopicStates"));
    fix_mu = matchstring(inputfile_, "Fix_mu");
    fixed_mu_value = double(matchstring(inputfile_, "fixed_mu_value")) * 1.0;
    BoundaryConnection = double(matchstring(inputfile_, "PBC"));

    assert(SavingMicroscopicStates_int == 1 ||
           SavingMicroscopicStates_int == 0);
    if (SavingMicroscopicStates_int == 1)
    {
        Saving_Microscopic_States = true;
    }
    else
    {
        Saving_Microscopic_States = false;
    }

    No_Of_Microscopic_States = int(matchstring(inputfile_, "NoOfMicroscopicStates"));

    ns = lx * ly;
    cout << "TotalNumberOfSites = " << ns << endl;

    Fill = matchstring(inputfile_, "Filling_per_site");
    Total_Particles = ns * Fill;
    cout << "TotalNumberOfParticles = " << Total_Particles << endl;


    IterMax = int(matchstring(inputfile_, "MaxMCsweeps"));
    MCNorm = 0.0; //matchstring(inputfile,"MCNorm")
    RandomSeed = matchstring(inputfile_, "RandomSeed");
    RandomDisorderSeed = matchstring(inputfile_, "RandomDisorderSeed");
    Disorder_Strength = matchstring(inputfile_, "Disorder_Strength");
    Boltzman_constant = matchstring(inputfile_, "Boltzman_constant");

    eta = matchstring(inputfile_, "eta");
    // omega_min = matchstring(inputfile_, "omega_min");
    // omega_max = matchstring(inputfile_, "omega_max");
    domega = matchstring(inputfile_, "domega");


    string J_Hund_str = "J_Hund";
    string OnSiteE_str = "OnSiteE";
    string J_Hund_str_temp, OnSiteE_str_temp;

    for(int n=0;n<n_orbs;n++){
        J_Hund_str_temp = J_Hund_str + to_string(n);
        OnSiteE_str_temp = OnSiteE_str + to_string(n);

        J_Hund[n] = matchstring(inputfile_, J_Hund_str_temp);
        OnSiteE[n] = matchstring(inputfile_, OnSiteE_str_temp);
    }

    //Hopping matrices -------------------


    hopping_NN_X.resize(n_orbs,n_orbs);
    hopping_NN_Y.resize(n_orbs,n_orbs);
    string Nearest_Neigh_Hopping_basis_X;
    string Nearest_Neigh_Hopping_basis_Y;

    string NN_X_str, NN_Y_str;
    for (int m=0;m<n_orbs;m++){

        NN_X_str = "Nearest_Neigh_Hopping_X_basis_row" + to_string(m);
        NN_Y_str = "Nearest_Neigh_Hopping_Y_basis_row" + to_string(m);
        Nearest_Neigh_Hopping_basis_X=matchstring2(inputfile_, NN_X_str);
        Nearest_Neigh_Hopping_basis_Y=matchstring2(inputfile_, NN_Y_str);

        stringstream stream_X(Nearest_Neigh_Hopping_basis_X);
        stringstream stream_Y(Nearest_Neigh_Hopping_basis_Y);

        for(int n=0;n<n_orbs;n++){
            stream_X >> hopping_NN_X(m,n);
            stream_Y >> hopping_NN_Y(m,n);

        }
    }



    //Next Nearest hopping------------
    hopping_NNN_PXPY.resize(n_orbs,n_orbs);
    hopping_NNN_PXMY.resize(n_orbs,n_orbs);
    //If needed read from input file

    //Hopping matrices done---------------




    lambda_lattice = matchstring (inputfile_, "lambda_lattice");
    K1x = matchstring(inputfile_, "K");
    K1y = K1x;
    cout << "K1x= " << K1x << endl;

    Dflag = 'N';

    metropolis_double = double(matchstring(inputfile_, "Metropolis_Algo"));
    if (metropolis_double == 1.0)
    {
        Metropolis_Algorithm = true;
        Heat_Bath_Algorithm = false;
    }
    else if (metropolis_double == 0.0)
    {
        Metropolis_Algorithm = false;
        Heat_Bath_Algorithm = true;
    }
    else
    {
        cout << "ERROR: Metropolis_Algo can be only 1 (true) or 0 (false)" << endl;
        assert(metropolis_double == 0.0);
    }


    cooling_double = double(matchstring(inputfile_, "Cooling"));
    if (cooling_double == 1.0)
    {
        Cooling_ = true;

        temp_min = double(matchstring(inputfile_, "Temperature_min"));
        temp_max = double(matchstring(inputfile_, "Temperature_max"));
        d_Temp = double(matchstring(inputfile_, "dTemperature"));
        beta_max = double(Boltzman_constant / temp_min);
        beta_min = double(Boltzman_constant / temp_max);

        no_of_temp_points = int( (( temp_max - temp_min )/(d_Temp)) + 1);
        Temp_values.resize(no_of_temp_points);
        for(int point_no=0;point_no<no_of_temp_points;point_no++){
            Temp_values[point_no] = temp_min + (point_no*(d_Temp));
        }

    }
    else if (cooling_double == 2.0)
    {
        Cooling_ = true;
        temp_values_ = matchstring2(inputfile_, "Temperature_Values");

        stringstream temp_values_stream(temp_values_);
        temp_values_stream>>no_of_temp_points;

        Temp_values.resize(no_of_temp_points);

        for(int point_no=0;point_no<no_of_temp_points;point_no++){
            temp_values_stream >> Temp_values[point_no];
        }

    }
    else if (cooling_double == 0.0)
    {
        Cooling_ = false;

        temp = double(matchstring(inputfile_, "Temperature")); // temperature in kelvin
        beta = double(1.0/(Boltzman_constant*temp));                         //Beta which is (T*k_b)^-1

        temp_min = temp;
        temp_max = temp;
        d_Temp = 10.0; //arbitrary positive number

        Temp_values.resize(1);
        Temp_values[0]=temp_min;
    }
    else
    {
        cout << "ERROR: Cooling can be only 1, 0 or 2" << endl;
        assert(cooling_double == 0.0);
    }

    ED_double = double(matchstring(inputfile_, "Perform_ED"));
    if (ED_double == 1.0)
    {
        ED_ = true;
        lx_cluster = lx;
        ly_cluster = ly;
    }
    else if (ED_double == 0.0)
    {
        ED_ = false;
    }
    else
    {
        cout << "ERROR: Perform_ED can be only 1 (true) or 0 (false)" << endl;
        assert(ED_double == 0.0);
    }

    Read_Seed_from_file_double = double(matchstring(inputfile_, "Read_Seed_from_file"));
    if (Read_Seed_from_file_double == 1.0)
    {
        Read_Seed_from_file_ = true;
    }
    else if (Read_Seed_from_file_double == 0.0)
    {
        Read_Seed_from_file_ = false;
    }
    else
    {
        cout << "ERROR:Read_Seed_from_file can be only 1 (true) or 0 (false)" << endl;
        assert(Read_Seed_from_file_double == 0.0);
    }

    Seed_file_name_ = matchstring2(inputfile_, "Seed_file_name");

    Last_n_sweeps_for_measurement = int(matchstring(inputfile_, "Last_n_sweeps_for_measurement"));
    Measurement_after_each_m_sweeps = int(matchstring(inputfile_, "Measurement_after_each_m_sweeps"));

    pi = 4.00 * atan(double(1.0));
    Eav = 0.0;
    AccCount[0] = 0;
    AccCount[1] = 0;

    WindowSize = double(0.01);
    mus = 0.25;
    hopping_intracell=-1.0;
    assert(n_orbs==3);
    cout << "____________________________________" << endl;
}

double Parameters::matchstring(string file, string match)
{
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass = false;
    while (std::getline(readFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass == false)
        {
            // ---------------------------------
            if (iss >> amount && test == match)
            {
                // cout << amount << endl;
                pass = true;
            }
            else
            {
                pass = false;
            }
            // ---------------------------------
            if (pass)
                break;
        }
    }
    if (pass == false)
    {
        string errorout = match;
        errorout += "= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file, string match)
{

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if (readFile.is_open())
    {
        while (!readFile.eof())
        {
            getline(readFile, line);

            if ((offset = line.find(match, 0)) != string::npos)
            {
                amount = line.substr(offset + match.length() + 1);
            }
        }
        readFile.close();
    }
    else
    {
        cout << "Unable to open input file while in the Parameters class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}

#endif

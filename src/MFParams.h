#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams
{
public:
    // Define Fields
    Mat_3_doub etheta, ephi;

    Mat_3_doub Sz, Sx, Sy;
    Mat_3_doub etheta_avg, ephi_avg;
    Mat_3_doub Moment_Size;
    Mat_3_doub Disorder;

    // Constructor
    MFParams(Parameters &Parameters__, Coordinates &Coordinates__, mt19937_64 &Generator1__, mt19937_64 &Generator2__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }

    double random1();
    double random2();
    void FieldThrow(int site, int orb, string mc_dof_type);
    void initialize();
    void Adjust_MCWindow();
    void Calculate_Fields_Avg();
    void Read_classical_DOFs(string filename);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_, ly_, ns_;

    uniform_real_distribution<double> dis1_; //for random fields
    uniform_real_distribution<double> dis2_; //for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);
};

void MFParams::Adjust_MCWindow()
{
    double ratio;
    ratio = Parameters_.AccCount[0] / (Parameters_.AccCount[0] + Parameters_.AccCount[1]);
    //cout<<"ratio= "<< ratio << "temp= "<<Parameters_.temp << endl;
    Parameters_.AccCount[0] = 0;
    Parameters_.AccCount[1] = 0;
    Parameters_.WindowSize *= abs(1.0 + 1.0 * (ratio - 0.5));
    if(Parameters_.WindowSize>10){
        Parameters_.WindowSize=10.0;
    }
    //Parameters_.WindowSize =0.2;
    cout << "Ratio: " << ratio << "  window size:  " << Parameters_.WindowSize << endl;
    return;
} // ----------

void MFParams::FieldThrow(int site, int orb, string mc_dof_type)
{
    int a, b;

    int Pi_multiple;

    double Pi = Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    a = Coordinates_.indx_cellwise(site);
    b = Coordinates_.indy_cellwise(site);

    //ANGLES
    if (mc_dof_type == "phi")
    {
        ephi[a][b][orb] += 2 * Pi * (random1() - 0.5) * MC_Window;

        Pi_multiple = ephi[a][b][orb]/Pi;


        if (ephi[a][b][orb] < 0.0)
        {
            ephi[a][b][orb] = -ephi[a][b][orb];
        }

        ephi[a][b][orb] = fmod(ephi[a][b][orb], 2.0 * Pi);


    }

    if (mc_dof_type == "theta")
    {
        etheta[a][b][orb] += Pi * (random1() - 0.5) * MC_Window;
        if (etheta[a][b][orb] < 0.0)
        {
            etheta[a][b][orb] = -etheta[a][b][orb];
        }

        etheta[a][b][orb] = fmod(etheta[a][b][orb],  Pi);

    }


    if (mc_dof_type == "theta_and_phi")
    {
        //phi
        ephi[a][b][orb] += 2 * Pi * (random1() - 0.5) * MC_Window;

        Pi_multiple = ephi[a][b][orb]/Pi;

        if (ephi[a][b][orb] < 0.0)
        {
            ephi[a][b][orb] = -ephi[a][b][orb];
        }

        ephi[a][b][orb] = fmod(ephi[a][b][orb], 2.0 * Pi);


        //theta
        etheta[a][b][orb] += Pi * (random1() - 0.5) * MC_Window;
        if (etheta[a][b][orb] < 0.0)
        {
            etheta[a][b][orb] = -etheta[a][b][orb];
        }

        etheta[a][b][orb] = fmod(etheta[a][b][orb],  Pi);
    }


} // ----------

double MFParams::random1()
{

    return dis1_(Generator1_);
}

double MFParams::random2()
{

    return dis2_(Generator2_);
}

void MFParams::initialize()
{

    bool Diagonal_ZigZag_Ising_alongZ=false;
    bool Diagonal_ZigZag_Ising_alongZ_rotatedby90deg=true;
    bool two_by_two_Plaquettes_Ising_alongZ=false;
    bool FM_state_Ising=false;
    bool AFM_state_Ising=false;
    lx_ = Coordinates_.lx_;
    ly_ = Coordinates_.ly_;

    // srand(Parameters_.RandomSeed);

    //    etheta_avg.resize(lx_, ly_);
    //    ephi_avg.resize(lx_, ly_);
    Disorder.resize(lx_);
    etheta.resize(lx_);
    ephi.resize(lx_);
    Moment_Size.resize(lx_);
    Sz.resize(lx_);
    Sx.resize(lx_);
    Sy.resize(lx_);
    for(int ix=0;ix<lx_;ix++){
        Disorder[ix].resize(ly_);
        etheta[ix].resize(ly_);
        ephi[ix].resize(ly_);
        Moment_Size[ix].resize(ly_);
        Sz[ix].resize(ly_);
        Sx[ix].resize(ly_);
        Sy[ix].resize(ly_);
        for(int iy=0;iy<ly_;iy++){
            Disorder[ix][iy].resize(3);
            etheta[ix][iy].resize(3);
            ephi[ix][iy].resize(3);
            Moment_Size[ix][iy].resize(3);
            Sz[ix][iy].resize(3);
            Sx[ix][iy].resize(3);
            Sy[ix][iy].resize(3);
        }
    }


    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file << "#seed=" << Parameters_.RandomDisorderSeed << " for mt19937_64 Generator is used" << endl;
    Disorder_conf_file << "#ix   iy   orb  Dis[ix,iy,orb]" << endl;

    ofstream Initial_MC_DOF_file("Initial_MC_DOF_values");

    Initial_MC_DOF_file << "#seed=" << Parameters_.RandomSeed << " for mt19937_64 Generator is used" << endl;
    Initial_MC_DOF_file << "#ix   iy  orb  Theta(x,y,orb)    Phi(x,y,orb)      Moment_Size(x,y,orb)" << endl;


    string temp_string;

    int spin_offset;
    int ix_, iy_, orb_;

    //Initialization
    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {
            for(int orb =0;orb<3;orb++){

                ephi[i][j][orb] = 0.0;
                etheta[i][j][orb]= 0.0;
            }
        }
    }


    if (Parameters_.Read_Seed_from_file_ == true)
    {
        cout<<"Configuration read from : '"<<Parameters_.Seed_file_name_<<"'"<<endl;
        ifstream Initial_Seed(Parameters_.Seed_file_name_);
        getline(Initial_Seed, temp_string);
        //cout << temp_string << endl;
        for (int ix = 0; ix < lx_; ix++)
        {
            for (int iy = 0; iy < ly_; iy++)
            {
                for(int orb =0;orb<3;orb++){
                    Initial_Seed >> ix_ >> iy_ >> orb_ >> etheta[ix][iy][orb] >> ephi[ix][iy][orb] >> Moment_Size[ix][iy][orb];
                    assert(ix_ == ix);
                    assert(iy_ == iy);
                    assert(orb_ ==orb);
                }
            }
        }
    }

    else
    {
        for (int j = 0; j < ly_; j++)
        {
            for (int i = 0; i < lx_; i++)
            {
                for(int orb =0;orb<3;orb++){

                    //RANDOM fields
                    if (Parameters_.MC_on_theta_and_phi == true)
                    {
                        ephi[i][j][orb] = 2.0 * random1() * PI;
                        etheta[i][j][orb] = random1() * PI;
                    }
                    else
                    {
                        if (Parameters_.MC_on_phi == true)
                        {
                            ephi[i][j][orb] = 2.0 * random1() * PI;
                        }

                        if (Parameters_.MC_on_theta == true)
                        {
                            etheta[i][j][orb] = random1() * PI;
                        }

                    }


                    Moment_Size[i][j][orb] = 1.0;
                }

            }
        }

        for (int j = 0; j < ly_; j++)
        {
            for (int i = 0; i < lx_; i++)
            {
                for(int orb =0;orb<3;orb++){
                    etheta[i][j][orb] += random1()*0.05;
                    ephi[i][j][orb] += random1()*0.05;
                    Initial_MC_DOF_file << i << setw(15) << j << setw(15) << orb<< setw(15) <<etheta[i][j][orb] << setw(15) << ephi[i][j][orb]
                                           << setw(15) << Moment_Size[i][j][orb] << endl;
                }
            }
        }
    }

    //RANDOM Disorder
    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {
            for(int orb =0;orb<3;orb++){
                Disorder[i][j][orb] = Parameters_.Disorder_Strength * ((2.0 * random2()) - 1.0);
                Disorder_conf_file << i << "  " << j << "  "<<orb<<"   "<< Disorder[i][j][orb] << endl;
            }
        }
        Disorder_conf_file << endl;
    }

} // ----------

void MFParams::Calculate_Fields_Avg()
{

    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {

            //            ephi_avg(i, j) = ephi_avg(i, j) + ephi(i, j);
            //            etheta_avg(i, j) = etheta_avg(i, j) + etheta(i, j);
        }
    }

} // ----------

void MFParams::Read_classical_DOFs(string filename)
{

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    getline (fl_in,tmp_str);

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            for(int orb =0;orb<3;orb++){
                fl_in >> tmp_double >> tmp_double >> tmp_double >> etheta[i][j][orb] >> ephi[i][j][orb]>> tmp_double;
            }
        }
    }

} // ----------

#endif

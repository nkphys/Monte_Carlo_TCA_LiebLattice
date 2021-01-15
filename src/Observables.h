#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

class Observables
{
public:
    Observables(Parameters &Parameters__, Coordinates &Coordinates__,
                MFParams &MFParams__, Hamiltonian &Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ncells_(Coordinates_.ncells_), nbasis_(Coordinates_.nbasis_), n_orbs_(Coordinates_.n_orbs_)
    {
        Initialize();
    }

    void Initialize();

    void OccDensity();
    void Calculate_Akw();
    void Calculate_Nw();
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    void TotalOccDensity();
    void DensityOfStates();
    void Calculate_OrbResolved_Nw();
    void SiSjFULL();
    double fermi_function(int n);

    void calculate_quantum_SiSj();
    void quantum_SiSjQ_Average();
    void quantum_SiSj_Average();

    void SiSjQ_Average();
    void SiSj_Average();
    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void OccDensity(int tlabel);
    void DOSprint(int tlabel);
    complex<double> SiSjQ(int i, int j);
    double SiSj(int i, int j);
    double Omega(int i);

    complex<double> SiSjQ_Mean(int i, int j);
    complex<double> SiSjQ_square_Mean(int i, int j);

    double SiSj_Mean(int i, int j);
    double SiSj_square_Mean(int i, int j);

    double BandWidth;
    Matrix<complex<double>> SiSjQ_, SiSjQ_Mean_, SiSjQ_square_Mean_;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    int lx_, ly_, ncells_, nbasis_;
    int n_orbs_;
    double dosincr_, tpi_;
    Matrix<double> SiSj_, dos;
    vector<double> sx_, sy_, sz_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;

    Mat_2_doub local_density;
    Mat_2_doub local_density_Mean;
    Mat_2_doub local_density_square_Mean;

    Matrix<complex<double>> quantum_SiSjQ_, quantum_SiSjQ_Mean_, quantum_SiSjQ_square_Mean_;
    Matrix<complex<double>> quantum_SiSj_, quantum_SiSj_Mean_, quantum_SiSj_square_Mean_;

    void calculate_local_density();
    void local_density_average();
};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/

void Observables::Calculate_Akw()
{

    //NEED TO CHANGE

    //---------Read from input file-----------------------//
    string fileout = "Akw.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.001;
    omega_min = -8.0;
    omega_max = 4.0;
    d_omega = 0.03;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Akw_out(fileout.c_str());

    int c1, c2;
    int j_posx, j_posy;
    int l_posx, l_posy;

    Mat_3_Complex_doub A_up_00, A_up_11, A_up_22;
    Mat_3_Complex_doub A_dn_00, A_dn_11, A_dn_22;
    A_up_00.resize(Coordinates_.ncells_);
    A_dn_00.resize(Coordinates_.ncells_);
    A_up_11.resize(Coordinates_.ncells_);
    A_dn_11.resize(Coordinates_.ncells_);
    A_up_22.resize(Coordinates_.ncells_);
    A_dn_22.resize(Coordinates_.ncells_);

    for (int i = 0; i < Coordinates_.ncells_; i++)
    {
        A_up_00[i].resize(Coordinates_.ncells_);
        A_dn_00[i].resize(Coordinates_.ncells_);
        A_up_11[i].resize(Coordinates_.ncells_);
        A_dn_11[i].resize(Coordinates_.ncells_);
        A_up_22[i].resize(Coordinates_.ncells_);
        A_dn_22[i].resize(Coordinates_.ncells_);

        for (int j = 0; j < Coordinates_.ncells_; j++)
        {
            A_up_00[i][j].resize(omega_index_max);
            A_dn_00[i][j].resize(omega_index_max);
            A_up_11[i][j].resize(omega_index_max);
            A_dn_11[i][j].resize(omega_index_max);
            A_up_22[i][j].resize(omega_index_max);
            A_dn_22[i][j].resize(omega_index_max);
        }
    }

    complex<double> Nup_check(0, 0);
    complex<double> Ndn_check(0, 0);

    for (int j = 0; j < Coordinates_.ncells_; j++)
    {
        j_posx = Coordinates_.indx_cellwise(j);
        j_posy = Coordinates_.indy_cellwise(j);

        for (int l = 0; l < Coordinates_.ncells_; l++)
        {
            l_posx = Coordinates_.indx_cellwise(l);
            l_posy = Coordinates_.indy_cellwise(l);

            for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
            {
                A_up_00[j][l][omega_ind] = zero_complex;
                A_dn_00[j][l][omega_ind] = zero_complex;
                A_up_11[j][l][omega_ind] = zero_complex;
                A_dn_11[j][l][omega_ind] = zero_complex;
                A_up_22[j][l][omega_ind] = zero_complex;
                A_dn_22[j][l][omega_ind] = zero_complex;

                for (int n = 0; n < Hamiltonian_.Ham_.n_row(); n++)
                {

                    //c= l + or1*ns_ + ns_*orbs_*spin;

                    //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                    c1 = Coordinates_.Nbasis(l_posx, l_posy, 0) + Coordinates_.nbasis_;
                    c2 = Coordinates_.Nbasis(j_posx, j_posy, 0) + Coordinates_.nbasis_;
                    A_dn_00[j][l][omega_ind] += conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c2, n) *
                                                Lorentzian(omega_min + (omega_ind * d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = Coordinates_.Nbasis(l_posx, l_posy, 0);
                    c2 = Coordinates_.Nbasis(j_posx, j_posy, 0);
                    A_up_00[j][l][omega_ind] += conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c2, n) *
                                                Lorentzian(omega_min + (omega_ind * d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = Coordinates_.Nbasis(l_posx, l_posy, 1) + Coordinates_.nbasis_;
                    c2 = Coordinates_.Nbasis(j_posx, j_posy, 1) + Coordinates_.nbasis_;
                    A_dn_11[j][l][omega_ind] += conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c2, n) *
                                                Lorentzian(omega_min + (omega_ind * d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = Coordinates_.Nbasis(l_posx, l_posy, 1);
                    c2 = Coordinates_.Nbasis(j_posx, j_posy, 1);
                    A_up_11[j][l][omega_ind] += conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c2, n) *
                                                Lorentzian(omega_min + (omega_ind * d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = Coordinates_.Nbasis(l_posx, l_posy, 2) + Coordinates_.nbasis_;
                    c2 = Coordinates_.Nbasis(j_posx, j_posy, 2) + Coordinates_.nbasis_;
                    A_dn_22[j][l][omega_ind] += conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c2, n) *
                                                Lorentzian(omega_min + (omega_ind * d_omega) - Hamiltonian_.eigs_[n], eta);

                    c1 = Coordinates_.Nbasis(l_posx, l_posy, 2);
                    c2 = Coordinates_.Nbasis(j_posx, j_posy, 2);
                    A_up_22[j][l][omega_ind] += conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c2, n) *
                                                Lorentzian(omega_min + (omega_ind * d_omega) - Hamiltonian_.eigs_[n], eta);
                }

                if (j == l)
                {
                    Nup_check += (A_up_00[j][l][omega_ind] + A_up_11[j][l][omega_ind] + A_up_22[j][l][omega_ind]) * d_omega;
                    Ndn_check += (A_dn_00[j][l][omega_ind] + A_dn_11[j][l][omega_ind] + A_dn_22[j][l][omega_ind]) * d_omega;
                }
            }
        }
    }

    cout << "Nup_check = " << Nup_check << endl;
    cout << "Ndn_check = " << Ndn_check << endl;

    complex<double> temp_up_00, temp_up_11, temp_up_22;
    complex<double> temp_dn_00, temp_dn_11, temp_dn_22;
    double kx, ky;
    int kx_i, ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------\Gamma to X-----------------
    ky_i = 0;
    for (kx_i = 0; kx_i <= (Parameters_.lx / 2); kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i = (Parameters_.lx / 2);
    for (ky_i = 1; ky_i <= (Parameters_.lx / 2); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i = (Parameters_.lx / 2) - 1;
    ky_i = (Parameters_.lx / 2) - 1;
    for (kx_i = (Parameters_.lx / 2) - 1; kx_i >= -1; kx_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    double k22_offset = 0;
    for (int k_point = 0; k_point < k_path.size(); k_point++)
    {

        kx_i = k_path[k_point].first;
        ky_i = k_path[k_point].second;
        kx = (2.0 * PI * kx_i) / (1.0 * Parameters_.lx);
        ky = (2.0 * PI * ky_i) / (1.0 * Parameters_.ly);

        for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
        {
            temp_up_00 = zero_complex;
            temp_dn_00 = zero_complex;
            temp_up_11 = zero_complex;
            temp_dn_11 = zero_complex;
            temp_up_22 = zero_complex;
            temp_dn_22 = zero_complex;

            for (int j = 0; j < Coordinates_.ncells_; j++)
            {
                j_posx = Coordinates_.indx_cellwise(j);
                j_posy = Coordinates_.indy_cellwise(j);

                for (int l = 0; l < Coordinates_.ncells_; l++)
                {
                    l_posx = Coordinates_.indx_cellwise(l);
                    l_posy = Coordinates_.indy_cellwise(l);

                    temp_up_00 += one_complex *
                                  exp(iota_complex * (kx * (j_posx - l_posx) +
                                                      ky * (j_posy - l_posy))) *
                                  A_up_00[j][l][omega_ind];

                    temp_dn_00 += one_complex *
                                  exp(iota_complex * (kx * (j_posx - l_posx) +
                                                      ky * (j_posy - l_posy))) *
                                  A_dn_00[j][l][omega_ind];

                    temp_up_11 += one_complex *
                                  exp(iota_complex * (kx * (j_posx - l_posx) +
                                                      ky * (j_posy - l_posy))) *
                                  A_up_11[j][l][omega_ind];

                    temp_dn_11 += one_complex *
                                  exp(iota_complex * (kx * (j_posx - l_posx) +
                                                      ky * (j_posy - l_posy))) *
                                  A_dn_11[j][l][omega_ind];

                    temp_up_22 += one_complex *
                                  exp(iota_complex * (kx * (j_posx - l_posx) +
                                                      ky * (j_posy - l_posy))) *
                                  A_up_22[j][l][omega_ind];

                    temp_dn_22 += one_complex *
                                  exp(iota_complex * (kx * (j_posx - l_posx) +
                                                      ky * (j_posy - l_posy))) *
                                  A_dn_22[j][l][omega_ind];
                }
            }
            //Use 1:6:7----for gnuplot
            file_Akw_out << k_point << "   " << kx_i << "   " << ky_i << "   " << (ky_i * Parameters_.lx) + kx_i << "    " << omega_min + (d_omega * omega_ind) << "   " << omega_ind << "    "
                         << temp_up_00.real() << "    " << temp_dn_00.real() << "    "
                         << temp_up_11.real() << "    " << temp_dn_11.real() << "    "
                         << temp_up_22.real() << "    " << temp_dn_22.real() << "    "
                         << temp_up_00.imag() << "     " << temp_dn_00.imag() << "    "
                         << temp_up_11.imag() << "     " << temp_dn_11.imag() << "    "
                         << temp_up_22.imag() << "     " << temp_dn_22.imag() << "    "
                         << endl;
        }
        file_Akw_out << endl;
    }
}

void Observables::Calculate_Nw()
{

    //---------Read from input file-----------------------//
    string fileout = "Nw_total_" + Parameters_.Seed_file_name_;
    double omega_min, omega_max, d_omega;
    double eta = Parameters_.eta;
    omega_min = Hamiltonian_.eigs_[0] - 5.0;
    omega_max = Hamiltonian_.eigs_[Hamiltonian_.eigs_.size() - 1] + 5.0;
    d_omega = Parameters_.domega;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Nw_out(fileout.c_str());

    double temp_val;

    for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
    {

        temp_val = 0.0;

        for (int n = 0; n < Hamiltonian_.Ham_.n_row(); n++)
        {

            temp_val +=
                Lorentzian(omega_min + (omega_ind * d_omega) -Hamiltonian_.eigs_[n], eta);
        }

        file_Nw_out << omega_min + (omega_ind * d_omega) << "     "
                    << (1.0 / Hamiltonian_.eigs_.size()) * temp_val
                    << "     " << endl;
    }

    file_Nw_out << "#mu = " << Parameters_.mus << endl;

    ofstream spectrum;
    spectrum.open("spectrum" + Parameters_.Seed_file_name_);

    for (int i=0; i < Hamiltonian_.eigs_.size(); i++){
        spectrum << i << "\t" << Hamiltonian_.eigs_[i] << endl;
    }
    spectrum.close();
}

void Observables::Calculate_OrbResolved_Nw()
{

    //NEED TO CHANGE

    //---------Read from input file-----------------------//
    string fileout = "NwOrbResolved_total.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.1;
    omega_min = -100;
    omega_max = 100.0;
    d_omega = 0.001;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Nw_out(fileout.c_str());

    double temp_val;

    for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
    {

        temp_val = 0.0;

        for (int n = 0; n < Hamiltonian_.Ham_.n_row(); n++)
        {

            temp_val += Lorentzian(omega_min + (omega_ind * d_omega) - Hamiltonian_.eigs_[n], eta);
        }

        file_Nw_out << omega_min + (omega_ind * d_omega) << "     " << temp_val << "     " << endl;
    }

    file_Nw_out << "#mu = " << Parameters_.mus << endl;
}

void Observables::Get_Non_Interacting_dispersion()
{
}

double Observables::Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}

void Observables::DensityOfStates()
{
    //-----------Calculate Bandwidth------//
    // double omega_min, omega_max;
    // omega_min = Hamiltonian_.eigs_[0] - 5.0;
    // omega_max - Hamiltonian

    //-----------------------------------//
    // string fileout = "";

} // ----------

void Observables::OccDensity()
{

} // ----------

void Observables::TotalOccDensity()
{

} // ----------

complex<double> Observables::SiSjQ(int i, int j) { return SiSjQ_(i, j); }

double Observables::SiSj(int i, int j) { return SiSj_(i, j); }

complex<double> Observables::SiSjQ_Mean(int i, int j) { return SiSjQ_Mean_(i, j); }

complex<double> Observables::SiSjQ_square_Mean(int i, int j) { return SiSjQ_square_Mean_(i, j); }

double Observables::SiSj_Mean(int i, int j) { return SiSj_Mean_(i, j); }

double Observables::SiSj_square_Mean(int i, int j) { return SiSj_square_Mean_(i, j); }

double Observables::fermi_function(int n)
{
    double value;
    value = 1.0 / (exp(Parameters_.beta * (Hamiltonian_.eigs_[n] - Parameters_.mus)) + 1.0);
    return value;
}

void Observables::calculate_quantum_SiSj()
{
    /*
    Matrix<complex<double>> F_u_u;
    Matrix<complex<double>> F_d_d;
    Matrix<complex<double>> F_u_d;
    Matrix<complex<double>> F_d_u;
    Matrix<complex<double>> omF_u_u;
    Matrix<complex<double>> omF_d_d;
    Matrix<complex<double>> omF_u_d;
    Matrix<complex<double>> omF_d_u;
    int nx, ny;
    int jx, jy;
    F_u_u.resize(ns_, ns_);
    F_d_d.resize(ns_, ns_);
    F_u_d.resize(ns_, ns_);
    F_d_u.resize(ns_, ns_);
    omF_u_u.resize(ns_, ns_);
    omF_d_d.resize(ns_, ns_);
    omF_u_d.resize(ns_, ns_);
    omF_d_u.resize(ns_, ns_);
    for (int i = 0; i < ns_; i++)
    {
        for (int j = 0; j < ns_; j++)
        {
            for (int n = 0; n < Hamiltonian_.eigs_.size(); n++)
            {
                F_u_u(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j, n) * fermi_function(n));
                F_d_d(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j + ns_, n) * fermi_function(n));
                F_u_d(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j + ns_, n) * fermi_function(n));
                F_d_u(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j, n) * fermi_function(n));
                omF_u_u(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j, n) * (1.0 - fermi_function(n)));
                omF_d_d(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j + ns_, n) * (1.0 - fermi_function(n)));
                omF_u_d(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j + ns_, n) * (1.0 - fermi_function(n)));
                omF_d_u(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j, n) * (1.0 - fermi_function(n)));
            }
        }
    }

    int i_;
    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {
            // i = Coordinates_.Nc(ix, iy);

            quantum_SiSj_(ix, iy) = 0.0;
            for (int j = 0; j < ns_; j++)
            {
                jx = Coordinates_.indx(j);
                jy = Coordinates_.indy(j);
                nx = (jx + ix) % lx_;
                ny = (jy + iy) % ly_;
                i_ = Coordinates_.Nc(nx, ny);
                quantum_SiSj_(ix, iy) += (
                0.25*(F_u_u(i_, i_) * F_u_u(j, j) + F_u_u(i_, j) * omF_u_u(j, i_)
                    - ( F_u_u(i_, i_) * F_d_d(j, j) + F_u_d(i_, j) * omF_d_u(j, i_) )
                    - ( F_d_d(i_, i_) * F_u_u(j, j) + F_d_u(i_, j) * omF_u_d(j, i_) )
                    + F_d_d(i_, i_) * F_d_d(j, j) + F_d_d(i_, j) * omF_d_d(j, i_))
                    + 0.5 * (F_u_d(i_, i_) * F_d_u(j, j) + F_u_u(i_, j) * omF_d_d(j, i_))
                    + 0.5 * (F_d_u(i_, i_) * F_u_d(j, j) + F_d_d(i_, j) * omF_u_u(j, i_))
                    ).real();
            }
            quantum_SiSj_(ix, iy) /= (ns_ * 1.0);
        }
    }

    //Fourier transform
    double phase, Cos_ij, Sin_ij;
    for (int qx = 0; qx < lx_; qx++)
    {
        for (int qy = 0; qy < ly_; qy++)
        {
            quantum_SiSjQ_(qx, qy) = complex<double>(0.0, 0.0);
            for (int xr = 0; xr < lx_; xr++)
            {
                for (int yr = 0; yr < ly_; yr++)
                {
                    phase = 2.0 * Parameters_.pi * (double(qx * xr) / double(lx_) + double(qy * yr) / double(ly_));
                    Cos_ij = cos(phase);
                    Sin_ij = sin(phase);
                    quantum_SiSjQ_(qx, qy) += quantum_SiSj_(xr, yr) * complex<double>(Cos_ij, Sin_ij);
                }
            }
            quantum_SiSjQ_(qx, qy) *= double(1.0 / (lx_ * ly_));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

    */
}

void Observables::calculate_local_density()
{
    int c1;
    complex<double> value = zero_complex;
    // cout <<"Parameter mus="<< Parameters_.mus<<endl;
    for (int i = 0; i < nbasis_; i++)
    {
        for (int sigma = 0; sigma < 2; sigma++)
        {
            local_density[i][sigma] = 0.0;
            c1 = i + (sigma * nbasis_);
            for (int n = 0; n < Hamiltonian_.eigs_.size(); n++)
            {
                local_density[i][sigma] += (conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c1, n) * fermi_function(n)).real();
            }

            // value += (conj(Hamiltonian_.Ham_(c1, 1)) * Hamiltonian_.Ham_(c1, 1));
        }
    }
}

void Observables::SiSjFULL()
{

    double Cos_ij, Sin_ij, ei, ai, ej, aj, phase;
    double sxi, syi, szi, sxj, syj, szj;
    int site_, site_p, ax, ay, ix, iy, iorb, jx, jy, jorb;

    for (int i = 0; i < 3 * lx_ * ly_; i++)
    {
        ix = Coordinates_.indx_basiswise_[i];
        iy = Coordinates_.indy_basiswise_[i];
        iorb = Coordinates_.indorb_basiswise_[i];
        ei = MFParams_.etheta[ix][iy][iorb];
        ai = MFParams_.ephi[ix][iy][iorb];
        sxi = MFParams_.Moment_Size[ix][iy][iorb] * cos(ai) * sin(ei);
        syi = MFParams_.Moment_Size[ix][iy][iorb] * sin(ai) * sin(ei);
        szi = MFParams_.Moment_Size[ix][iy][iorb] * cos(ei);

        for (int j = 0; j < 3 * lx_ * ly_; j++)
        {
            jx = Coordinates_.indx_basiswise_[j];
            jy = Coordinates_.indy_basiswise_[j];
            jorb = Coordinates_.indorb_basiswise_[j];
            ej = MFParams_.etheta[jx][jy][jorb];
            aj = MFParams_.ephi[jx][jy][jorb];
            sxj = MFParams_.Moment_Size[jx][jy][jorb] * cos(aj) * sin(ej);
            syj = MFParams_.Moment_Size[jx][jy][jorb] * sin(aj) * sin(ej);
            szj = MFParams_.Moment_Size[jx][jy][jorb] * cos(ej);

            SiSj_(i, j) = (sxi * sxj) + (syi * syj) + (szi * szj);

            //cout << xr << " "<< yr<< " "<<  SiSj_(xr,yr) << endl;
        }
    }

    //    for (int qx = 0; qx < lx_; qx++)
    //    {
    //        for (int qy = 0; qy < ly_; qy++)
    //        {
    //            SiSjQ_(qx, qy) = complex<double>(0.0, 0.0);
    //            for (int xr = 0; xr < lx_; xr++)
    //            {
    //                for (int yr = 0; yr < ly_; yr++)
    //                {
    //                    phase = 2.0 * Parameters_.pi * (double(qx * xr) / double(lx_) + double(qy * yr) / double(ly_));
    //                    Cos_ij = cos(phase);
    //                    Sin_ij = sin(phase);
    //                    SiSjQ_(qx, qy) += SiSj_(xr, yr) * complex<double>(Cos_ij, Sin_ij);
    //                }
    //            }
    //            SiSjQ_(qx, qy) *= double(1.0 / (lx_ * ly_));
    //            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
    //        }
    //    }

} // ----------

void Observables::SiSjQ_Average()
{

    for (int qx = 0; qx < lx_; qx++)
    {
        for (int qy = 0; qy < ly_; qy++)
        {
            SiSjQ_Mean_(qx, qy) += SiSjQ_(qx, qy);
            SiSjQ_square_Mean_(qx, qy) += (SiSjQ_(qx, qy) * SiSjQ_(qx, qy));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

} // ----------

void Observables::quantum_SiSjQ_Average()
{

    for (int qx = 0; qx < lx_; qx++)
    {
        for (int qy = 0; qy < ly_; qy++)
        {
            quantum_SiSjQ_Mean_(qx, qy) += quantum_SiSjQ_(qx, qy);
            quantum_SiSjQ_square_Mean_(qx, qy) += (quantum_SiSjQ_(qx, qy) * quantum_SiSjQ_(qx, qy));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }
}

void Observables::SiSj_Average()
{

    for (int x = 0; x < lx_ * ly_ * 3; x++)
    {
        for (int y = 0; y < ly_ * ly_ * 3; y++)
        {
            SiSj_Mean_(x, y) += SiSj_(x, y);
            SiSj_square_Mean_(x, y) += (SiSj_(x, y) * SiSj_(x, y));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

} // ----------

void Observables::quantum_SiSj_Average()
{

    for (int x = 0; x < lx_; x++)
    {
        for (int y = 0; y < ly_; y++)
        {
            quantum_SiSj_Mean_(x, y) += quantum_SiSj_(x, y);
            quantum_SiSj_square_Mean_(x, y) += (quantum_SiSj_(x, y) * quantum_SiSj_(x, y));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

} // ----------

void Observables::local_density_average()
{
    for (int i = 0; i < nbasis_; i++)
    {
        for (int sigma = 0; sigma < 2; sigma++)
        {
            local_density_Mean[i][sigma] += local_density[i][sigma];
            local_density_square_Mean[i][sigma] += pow(local_density[i][sigma], 2);
        }
    }
}

void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE)
{

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE) * (Curr_QuantE + CurrE);
}

void Observables::OccDensity(int tlabel)
{

} // ----------

void Observables::DOSprint(int tlabel)
{

} // ----------

void Observables::Initialize()
{

    complex<double> zero(0.0, 0.0);
    int space = 2 * ncells_ * n_orbs_;
    sx_.resize(ncells_ * n_orbs_);
    sy_.resize(ncells_ * n_orbs_);
    sz_.resize(ncells_ * n_orbs_);

    local_density.resize(nbasis_);
    local_density_Mean.resize(nbasis_);
    local_density_square_Mean.resize(nbasis_);
    for (int i = 0; i < local_density.size(); i++)
    {
        local_density[i].resize(2);
        local_density_Mean[i].resize(2);
        local_density_square_Mean[i].resize(2);
    }

    Nematic_order_mean_ = 0.0;
    Nematic_order_square_mean_ = 0.0;

    SiSj_.resize(lx_ * ly_ * 3, lx_ * ly_ * 3);
    SiSj_Mean_.resize(lx_ * ly_ * 3, lx_ * ly_ * 3);
    SiSj_square_Mean_.resize(lx_ * ly_ * 3, lx_ * ly_ * 3);

    SiSjQ_Mean_.resize(lx_, ly_);
    SiSjQ_square_Mean_.resize(lx_, ly_);
    SiSjQ_.resize(lx_, ly_);

    quantum_SiSj_.resize(lx_, ly_);
    quantum_SiSj_Mean_.resize(lx_, ly_);
    quantum_SiSj_square_Mean_.resize(lx_, ly_);

    quantum_SiSjQ_Mean_.resize(lx_, ly_);
    quantum_SiSjQ_square_Mean_.resize(lx_, ly_);
    quantum_SiSjQ_.resize(lx_, ly_);

    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {
            SiSjQ_Mean_(ix, iy) = zero;
            SiSjQ_square_Mean_(ix, iy) = zero;
        }
    }

} // ----------

double Observables::Omega(int i)
{
    return -20.0 + double(i) * dosincr_;
} // ----------

#endif // OBSERVABLES_H

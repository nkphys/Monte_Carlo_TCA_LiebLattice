#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "tensor_type.h"

#ifndef MCENGINE_H
#define MCENGINE_H

class MCEngine
{
public:
    MCEngine(Parameters &Parameters__, Coordinates &Coordinates__,
             MFParams &MFParams__, Hamiltonian &Hamiltonian__,
             Observables &Observables__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__),
          MFParams_(MFParams__), Hamiltonian_(Hamiltonian__),
          Observables_(Observables__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ncells_(Coordinates_.ncells_),
          nbasis_(Coordinates_.nbasis_),n_orbs_(Parameters_.n_orbs),
          ED_(Parameters_.ED_)
    {
    }

    void RUN_MC();
    double Prob(double muu, double mu_new);
    double ProbCluster(double muu, double mu_new);
    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    Hamiltonian &Hamiltonian_;
    Observables &Observables_;
    const int lx_, ly_, n_orbs_, ncells_, nbasis_;
    bool ED_;
};

/*
 * ***********
 *  Functions in Class MCEngine ------
 *  ***********
*/

void MCEngine::RUN_MC()
{

    complex<double> zero(0.0, 0.0);

    int MC_sweeps_used_for_Avg = Parameters_.Last_n_sweeps_for_measurement;
    int Gap_bw_sweeps = Parameters_.Measurement_after_each_m_sweeps;

    double PrevE, CurrE, P_new, P12, muu_prev;
    double Prob_check;
    double muu_prevCluster;
    double Prev_QuantE;
    double Curr_QuantECluster;
    double Prev_QuantECluster;
    int x, y, act;
    double saved_Params[6];

    string File_Out_progress;
    string File_Out_theta_phi;
    string File_Out_real_space_corr;
    string File_Out_q_space_corr;

    double temp_ = Parameters_.temp_max;

    double initial_mu_guess;
    int n_states_occupied_zeroT;

    //starting with a random guess

    //while (temp_ >= Parameters_.temp_min)
    for(int temp_point=0;temp_point<Parameters_.Temp_values.size();temp_point++){
        temp_ = Parameters_.Temp_values[temp_point];
        cout << "Temperature = " << temp_ << " is being done" << endl;
        Parameters_.temp = temp_;
        Parameters_.beta = double( 1.0/(Parameters_.Boltzman_constant*temp_));

        for (int ix = 0; ix < lx_*ly_*3; ix++)
        {
            for (int iy = 0; iy < lx_*ly_*3; iy++)
            {
//                Observables_.SiSjQ_Mean_(ix, iy) = zero;
//                Observables_.SiSjQ_square_Mean_(ix, iy) = zero;
                Observables_.SiSj_square_Mean_(ix, iy) = 0.0;
                Observables_.SiSj_Mean_(ix, iy) = 0.0;
//                for(int orb=0;orb<n_orbs_;orb++){
//                    Observables_.local_density_Mean[Coordinates_.Nbasis(ix, iy, orb)][0] = 0.0;
//                    Observables_.local_density_Mean[Coordinates_.Nbasis(ix, iy, orb)][1] = 0.0;
//                    Observables_.local_density_square_Mean[Coordinates_.Nbasis(ix, iy, orb)][0] = 0.0;
//                    Observables_.local_density_square_Mean[Coordinates_.Nbasis(ix, iy, orb)][1] = 0.0;
//                }
                //                Observables_.quantum_SiSjQ_(ix, iy) = zero;
                //                Observables_.quantum_SiSjQ_Mean_(ix, iy) = zero;
                //                Observables_.quantum_SiSjQ_square_Mean_(ix, iy) = zero;
                //                Observables_.quantum_SiSj_(ix, iy) = zero;
                //                Observables_.quantum_SiSj_Mean_(ix, iy) = zero;
                //                Observables_.quantum_SiSj_square_Mean_(ix, iy) = zero;
            }
        }
        Observables_.AVG_Total_Energy = 0.0;
        Observables_.AVG_Total_Energy_sqr = 0.0;


//        MFParams_.etheta_avg.fill(0.0);
//        MFParams_.ephi_avg.fill(0.0);

        char temp_char[50];
        sprintf(temp_char, "%.10f", temp_);

        File_Out_progress = "output_Temp" + string(temp_char) + ".txt";
        ofstream file_out_progress(File_Out_progress.c_str());

        File_Out_theta_phi = "ThetaPhi_Temp" + string(temp_char) + ".txt";
        ofstream File_Out_Theta_Phi(File_Out_theta_phi.c_str());

        File_Out_real_space_corr = "classical_real_space_corr" + string(temp_char) + ".txt";
        ofstream File_Out_Real_Space_Corr(File_Out_real_space_corr.c_str());

//        File_Out_q_space_corr = "classical_momentum_space_corr" + string(temp_char) + ".txt";
//        ofstream File_Out_Q_Space_Corr(File_Out_q_space_corr.c_str());

        file_out_progress << "Total " << Parameters_.IterMax << " sweeps are performed." << endl;
        file_out_progress << "First " << Parameters_.IterMax - (Gap_bw_sweeps * (MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg) << " sweeps are used for thermalization and every " << Gap_bw_sweeps + 1 << " in last " << Gap_bw_sweeps * (MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg << " sweeps are used for measurement." << endl;
        act = 1;

        Parameters_.WindowSize = 0.1; //2f + 0.003f*beta0 ;
        Parameters_.Eav = 0.0;
        Parameters_.MCNorm = 0;
        Parameters_.Dflag = 'N'; //N // flag to calculate only Eigenvalue
        //std::string name="Output/Conf_" + to_string(ltemp) + ".dat";
        //Parameters_.beta = double(11604.0/ (Parameters_.temp +20.0) );
        //cout << "TEMP  " << Parameters_.temp << endl;

        file_out_progress << "I_MC" << setw(15) << "S(0,1)" << setw(15) << "S(0,2)"
                          << setw(17) << "< N_total >"
                          << setw(15) << "E_CL" << setw(15) << "E_QM" << setw(15) << "E_Total" << setw(15) << "mu" << endl;

        PrevE = Hamiltonian_.GetCLEnergy();
        Hamiltonian_.InteractionsCreate();
        //        Hamiltonian_.Ham_.print();
        //cout << "Here 1"<<endl;
        // Hamiltonian_.Check_Hermiticity();
        //cout << "Here 2" << endl;
        Hamiltonian_.Diagonalize(Parameters_.Dflag);

        //        for(int n=0;n<2*Parameters_.ns;n++){
        //            cout<<n<<"  "<<Hamiltonian_.eigs_[n]<<endl;
        //        }
        //        assert(false);

        n_states_occupied_zeroT = Coordinates_.nbasis_*(Parameters_.Fill/(n_orbs_));
        if(!Parameters_.fix_mu){
            initial_mu_guess = 0.5 * (Hamiltonian_.eigs_[n_states_occupied_zeroT - 1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
        }
        else{
            initial_mu_guess=Parameters_.fixed_mu_value;
        }
        //initial_mu_guess=0.25;
        Parameters_.mus = Hamiltonian_.chemicalpotential(initial_mu_guess, (Parameters_.Fill/(n_orbs_*2.0)));
        Prev_QuantE = Hamiltonian_.E_QM();
        muu_prev = Parameters_.mus;
        Hamiltonian_.copy_eigs(1);
        cout << "Initial Classical Energy[Full System] = " << PrevE << endl;
        cout << "Initial Quantum Energy[Full System] = " << Prev_QuantE << endl;
        cout << "Initial Total Energy[Full System] = " << PrevE + Prev_QuantE << endl;
        cout << "Initial mu=" << muu_prev << endl;

        int Confs_used = 0;
        int measure_start = 0;
        muu_prevCluster = muu_prev;

        if (ED_)
        {
            Prev_QuantECluster = Prev_QuantE;
            Hamiltonian_.eigsCluster_saved_ = Hamiltonian_.eigs_saved_;
        }

        for (int count = 0; count < Parameters_.IterMax; count++)
        {

            for (int i = 0; i < ncells_; i++)
            { // For each site

                for(int orb=0;orb<n_orbs_;orb++)
                {

                for (int mc_dof = 0; mc_dof < Parameters_.MC_DOF.size(); mc_dof++)
                {

                    //***Before change*************//

                    if (ED_ == false)
                    {
                        //TCA is used
                        Parameters_.Dflag = 'N'; //N
                        PrevE = Hamiltonian_.GetCLEnergy();

                        Hamiltonian_.InteractionsClusterCreate(i);
                        Hamiltonian_.DiagonalizeCluster(Parameters_.Dflag);

                        //n_states_occupied_zeroT=Parameters_.Fill*Hamiltonian_.eigsCluster_.size();
                        //initial_mu_guess=0.5*(Hamiltonian_.eigsCluster_[n_states_occupied_zeroT-1] + HamiltonianCluster_.eigs_[n_states_occupied_zeroT])
                        muu_prevCluster = Hamiltonian_.chemicalpotentialCluster(muu_prevCluster, (Parameters_.Fill/(n_orbs_*2.0)));
                        Prev_QuantECluster = Hamiltonian_.E_QMCluster();

                        Hamiltonian_.copy_eigs_Cluster(1);
                    }
                    else
                    {
                        assert(ED_);
                    }

                    //*******************************//

                    x = Coordinates_.indx_cellwise(i);
                    y = Coordinates_.indy_cellwise(i);

                    saved_Params[0] = MFParams_.etheta[x][y][orb];
                    saved_Params[1] = MFParams_.ephi[x][y][orb];
                    saved_Params[2] = MFParams_.Moment_Size[x][y][orb];


                    MFParams_.FieldThrow(i,orb, Parameters_.MC_DOF[mc_dof]);
                    CurrE = Hamiltonian_.GetCLEnergy();

                    if (count < (Parameters_.IterMax - (Gap_bw_sweeps * (MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)))
                    {
                        Parameters_.Dflag = 'N'; //N
                    }
                    else
                    {
                        if (ED_)
                        {
                            Parameters_.Dflag = 'N'; //quantum observables are calculated separately.
                        }
                        else
                        {
                            Parameters_.Dflag = 'N'; //N
                        }
                    }

                    Hamiltonian_.InteractionsClusterCreate(i);

                    Hamiltonian_.DiagonalizeCluster(Parameters_.Dflag);

                    Parameters_.mus_Cluster = Hamiltonian_.chemicalpotentialCluster(muu_prevCluster, (Parameters_.Fill/(n_orbs_*2.0)));

                    Curr_QuantECluster = Hamiltonian_.E_QMCluster();

                    //Ratio of Quantum partition functions
                    /*P = [ Tr(exp(-beta(Hquant_new)))/Tr(exp(-beta(Hquant_old)))]*
                      [exp(-beta*E_classical(New)) / exp(-beta*E_classical(old))]
                     * [sin(Theta_i(New)) / sin(Theta_i(Old)) ]*/
                    /*exp(P12) = P
                  P12 = log (P)
                  */

                    //same mu-refrence is used, otherwise engine does not work properly
                    P_new = ProbCluster(muu_prevCluster*1.0, muu_prevCluster*1.0);
                    P12 = P_new - Parameters_.beta * ((CurrE) - (PrevE));
                    //P12 = - Parameters_.beta*((CurrE)-(PrevE));
                    //cout<<P12<<endl;
                    P12 += log((sin(MFParams_.etheta[x][y][orb]) / sin(saved_Params[0])));



                    //Heat bath algorithm [See page-129 of Prof. Elbio's Book]
                    //Heat bath algorithm works for small changes i.e. when P~1.0
                    //  if (Heat_Bath_Algo){
                    //     P =P/(1.0+P);
                    //  }
                    //Prob_check = exp(P12)/(1.0+exp(P12));

                    //Metropolis Algotithm
                    // if (Metropolis_Algo){
                    //    P=min(1.0,P);
                    // }
                    Prob_check = min(1.0,exp(P12));

                    /*
       * VON NEUMANN's REJECTING METHOD:
       * Random number < P -----> ACCEPT
       * Random number > P -----> REJECT
       */

                    //ACCEPTED
                    if (Prob_check > ( MFParams_.random1()) )
                    {
                        Parameters_.AccCount[0]++;
                        act = 1;
                        if (ED_)
                        {
                            PrevE = CurrE;
                            Prev_QuantECluster = Curr_QuantECluster;
                            Hamiltonian_.copy_eigs_Cluster(1);
                            muu_prevCluster = Parameters_.mus_Cluster;
                        }
                    }

                    //REJECTED
                    else
                    {
                        Parameters_.AccCount[1]++;
                        act = 0;
                        MFParams_.etheta[x][y][orb] = saved_Params[0];
                        MFParams_.ephi[x][y][orb] = saved_Params[1];
                        MFParams_.Moment_Size[x][y][orb] = saved_Params[2];
                    }

                    // if ((act == 1) && (count<1000)) {

                    //muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
                    //Parameters_.mus = Parameters_.mus*0.999 + muu*0.001;
                    //Parameters_.mus = muu;
                    //}

                } //MC_DOF_TYPE

            }//orb loop
            } // site loop

            //      if (act == 1) {
            //       muu = Hamiltonian_.chemicalpotential(Parameters_.mus,Parameters_.Fill);
            //       Parameters_.mus = Parameters_.mus*0.99f + muu*0.01f;
            //      }

            if ((count % 10 == 0))
            {
                MFParams_.Adjust_MCWindow();
            }

            if (count < (Parameters_.IterMax - (Gap_bw_sweeps * (MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg)))
            {
                if ((count % 10 == 0))
                {
                    Observables_.SiSjFULL();
                    file_out_progress << int(1.0 * count) << setw(20) << Observables_.SiSj(0, 1) << setw(16) << Observables_.SiSj(0, 2)
                                      << setw(16) << Hamiltonian_.ClusterDensity() << setw(16) << CurrE
                                      << setw(16) << Curr_QuantECluster << setw(16) << Curr_QuantECluster + CurrE << setw(15) << Parameters_.mus_Cluster << endl;
                }
            }
            //Average and Std. deviation is calculated is done
            else
            {

                /*********************/
                if (ED_ == false)
                {

                 /* No diagonalization for Full system is done
                    //TCA is used
                    Parameters_.Dflag = 'V';
                    Hamiltonian_.InteractionsCreate();
                    Hamiltonian_.Diagonalize(Parameters_.Dflag);

                    //n_states_occupied_zeroT=Parameters_.Fill*Hamiltonian_.eigsCluster_.size();
                    //initial_mu_guess=0.5*(Hamiltonian_.eigsCluster_[n_states_occupied_zeroT-1] + HamiltonianCluster_.eigs_[n_states_occupied_zeroT])
                    Parameters_.mus = Hamiltonian_.chemicalpotential(muu_prevCluster, Parameters_.Fill);
                */

                }
                else
                {
                    assert(ED_);
                    Parameters_.mus = Parameters_.mus_Cluster;
                    Hamiltonian_.eigs_ = Hamiltonian_.eigsCluster_;
                    Hamiltonian_.Ham_ = Hamiltonian_.HamCluster_;
                }
                //----------------------------------//

                if (measure_start == 0)
                {
                    measure_start++;
                    cout << "----------Measurement is started----------" << endl;
                    file_out_progress << "----------Measurement is started----------" << endl;
                    file_out_progress << "I_MC      Avg{S(pi,0)}    Avg{S(0,pi)}    std.dev{S(pi,0)}   std.dev{S(0,pi)}    Avg{q_S(pi,0)}    Avg{q_S(0,pi)}    std.dev{q_S(pi,0)}   std.dev{q_S(0,pi)}  Avg{S(1,0)}  Avg{S(0,1)}  std.dev{S(1,0)}  std.dev{S(0,1)}   Avg{E_classical}  std.dev{E_classical}   Avg.Filling" << endl;
                }
                int temp_count = count -
                        (Parameters_.IterMax - (Gap_bw_sweeps * (MC_sweeps_used_for_Avg - 1) + MC_sweeps_used_for_Avg));
                int zero_or_not = temp_count % (Gap_bw_sweeps + 1);
                if (zero_or_not == 0)
                {

                    Observables_.SiSjFULL();
//                    Observables_.SiSjQ_Average();
                    Observables_.SiSj_Average();
                    //Observables_.calculate_quantum_SiSj();
                    //Observables_.quantum_SiSjQ_Average();
                    //Observables_.quantum_SiSj_Average();
                    //Observables_.calculate_local_density();
                    //Observables_.local_density_average();

                    //Just Classical Energy
//                    Observables_.Total_Energy_Average(0.0, CurrE);
//                    MFParams_.Calculate_Fields_Avg();

                    if ((Parameters_.Saving_Microscopic_States == true) &&
                            (Confs_used < Parameters_.No_Of_Microscopic_States))
                    {

                        char Confs_char[50];
                        sprintf(Confs_char, "%d", Confs_used);
                        string File_Out_theta_phi_microState = "ThetaPhi_Temp" + string(temp_char) +
                                "MicroState" + string(Confs_char) + ".txt";
                        ofstream File_Out_Theta_Phi_MicroState(File_Out_theta_phi_microState.c_str());

                        File_Out_Theta_Phi_MicroState << "#x" << setw(15) << "y" << setw(15) << "orb"<<setw(15)<<"Theta(x,y,orb)" << setw(15) << "Phi(x,y,orb)"
                                                      << setw(15) << "Moment_Size(x,y,orb)"<< endl;
                        for (int ix = 0; ix < lx_; ix++)
                        {
                            for (int iy = 0; iy < ly_; iy++)
                            {
                                for(int iorb=0;iorb<n_orbs_;iorb++){
                                File_Out_Theta_Phi_MicroState << ix << setw(15) << iy << setw(15) <<  iorb << setw(15) << MFParams_.etheta[ix][iy][iorb] << setw(15) << MFParams_.ephi[ix][iy][iorb]
                                                              << setw(15) << MFParams_.Moment_Size[ix][iy][iorb]<< endl;
                                }
                            }
                        }

                    }

                    Confs_used = Confs_used + 1;

                    //double MC_steps_Avg_insitu = (1.0 + 1.0*(count - (Parameters_.IterMax - MC_steps_used_for_Avg)));

                    file_out_progress << int(1.0 * count) << setw(20) <<

                                         //-------------------------

                                         Observables_.SiSj_Mean(1, 0) / (Confs_used * 1.0)
                                      << setw(16) << Observables_.SiSj_Mean(0, 2) / (Confs_used * 1.0)
                                      << setw(16) << sqrt(((Observables_.SiSj_square_Mean(1, 0) / (Confs_used * 1.0)) - ((Observables_.SiSj_Mean(1, 0) * Observables_.SiSj_Mean(1, 0)) / (Confs_used * Confs_used * 1.0))))

                                      << setw(16) << sqrt(((Observables_.SiSj_square_Mean(0, 2) / (Confs_used * 1.0)) - ((Observables_.SiSj_Mean(0, 2) * Observables_.SiSj_Mean(0, 2)) / (Confs_used * Confs_used * 1.0))))
                                      << endl;
                }
            }

        } // Iter Loop

        file_out_progress << "Total " << Confs_used << " configurations were used were measurement" << endl;

        temp_ = temp_ - Parameters_.d_Temp;


//        File_Out_Q_Space_Corr << "qx" << setw(15) << "qy" << setw(15) << "<SSQ(qx,qy)>" << setw(15) << "sd(SSQ(qx,qy))" << endl;
//        // int temp_site_;
//        double qx, qy;
//        for (int ix = 0; ix < lx_; ix++)
//        {
//            for (int iy = 0; iy < ly_; iy++)
//            {
//                qx = 2 * 3.141593 * ix / (lx_ * 1.0);
//                qy = 2 * 3.141593 * iy / (ly_ * 1.0);

//                temp_site_ = Coordinates_.Ncell(ix, iy);
//                File_Out_Q_Space_Corr << qx << setw(15) << qy << setw(15) << Observables_.SiSjQ_Mean(ix, iy).real() / (Confs_used * 1.0)
//                                      << setw(15) << sqrt(((Observables_.SiSjQ_square_Mean(ix, iy).real() / (Confs_used * 1.0)) - ((Observables_.SiSjQ_Mean(ix, iy) * Observables_.SiSjQ_Mean(ix, iy)).real() / (Confs_used * Confs_used * 1.0)))) << endl;
//            }
//            File_Out_Q_Space_Corr << endl;
//        }

        string File_Out_theta_phi_microState0_toread = "ThetaPhi_Temp" + string(temp_char) +
                "MicroState0.txt";
        MFParams_.Read_classical_DOFs(File_Out_theta_phi_microState0_toread);
    } //Temperature loop

} // ---------

double MCEngine::Prob(double muu, double mu_new)
{

    double P = 0.0;
    double X, Y, X2;

    for (int i = 0; i < Hamiltonian_.eigs_.size(); i++)
    {
        X = Parameters_.beta * ((mu_new)-Hamiltonian_.eigs_[i]);
        Y = Parameters_.beta * ((muu)-Hamiltonian_.eigs_saved_[i]);
        //P += log(1 + exp(X)) - log(1 + exp(Y));

        if (X > 5)
        {
            P += X;
        }
        else if (fabs(X) < 0.001)
        {
            P += log(2.0 + X);
        }
        else if (X < -5)
        {
            P += exp(X);
        }
        else
        {
            P += log(1.0 + exp(X));
        }

        if (Y > 5)
        {
            P -= Y;
        }
        else if (fabs(Y) < 0.001)
        {
            P -= log(2.0 + Y);
        }
        else if (Y < -5)
        {
            P -= exp(Y);
        }
        else
        {
            P -= log(1.0 + exp(Y));
        }
    }

    return P;

} // ---------

double MCEngine::ProbCluster(double muu, double mu_new)
{

    double P = 0.0;
    double X, Y;
    int ns = (Parameters_.lx_cluster) * (Parameters_.ly_cluster);

    for (int i = 0; i < Hamiltonian_.eigsCluster_.size(); i++)
    {
        X = Parameters_.beta * ((mu_new)-Hamiltonian_.eigsCluster_[i]);
        Y = Parameters_.beta * ((muu)-Hamiltonian_.eigsCluster_saved_[i]);
        //P += log(1 + exp(X)) - log(1 + exp(Y));

        if (X > 15)
        {
            P += X;
        }
        else if (fabs(X) < 0.001)
        {
            P += log(2.0 + X);
        }
        else if (X < -15)
        {
            P += exp(X);
        }
        else
        {
            P += log(1.0 + exp(X));
        }

        if (Y > 15)
        {
            P -= Y;
        }
        else if (fabs(Y) < 0.001)
        {
            P -= log(2.0 + Y);
        }
        else if (Y < -15)
        {
            P -= exp(Y);
        }
        else
        {
            P -= log(1.0 + exp(Y));
        }
    }

    return P;

} // ---------

#endif // MCENGINE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;

#include "Matrix.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "MCEngine.h"
#include "random"


int main(int argc, char *argv[]) {
    if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }

    string ex_string_original =argv[0];
    cout<<ex_string_original<<endl;
    string ex_string;
    ex_string=ex_string_original.substr (ex_string_original.length() - 2);
    cout<<ex_string<<endl;

    string inputfile = argv[1];

    bool check_Non_Int=false;


    Parameters Parameters_;
    Parameters_.Initialize(inputfile);

    Coordinates Coordinates_(Parameters_.lx, Parameters_.ly, Parameters_.n_orbs);
    Coordinates CoordinatesCluster_(Parameters_.lx_cluster, Parameters_.ly_cluster, Parameters_.n_orbs);

    mt19937_64 Generator_(Parameters_.RandomSeed); //for random fields
    mt19937_64 Generator2_(Parameters_.RandomDisorderSeed); //for random disorder

    MFParams MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);

    Hamiltonian Hamiltonian_(Parameters_,Coordinates_,CoordinatesCluster_,MFParams_);


    Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);




    if(ex_string=="NI"){

        //Parameters_.J_HUND=0.0;
        Hamiltonian_.InteractionsCreate();
        //  Hamiltonian_.Ham_.print();
        // Hamiltonian_.Check_up_down_symmetry();
        //Hamiltonian_.Check_Hermiticity();
        Hamiltonian_.Diagonalize('V');
        int temp=Coordinates_.nbasis_*Parameters_.Fill*2.0;
        cout<<"mu for n=4 = "<<0.5*(Hamiltonian_.eigs_[temp-1] + Hamiltonian_.eigs_[temp])<<"   "<<
             Hamiltonian_.eigs_[temp-1]<<"   "<<Hamiltonian_.eigs_[temp]<<endl;
        Parameters_.mus=0.5*(Hamiltonian_.eigs_[temp-1] + Hamiltonian_.eigs_[temp]);
        double Quantum_E=Hamiltonian_.E_QM();
        double Classical_E=Hamiltonian_.GetCLEnergy();
        cout<<setprecision(9);
        cout<<"Total_Energy = "<<Quantum_E+Classical_E<<endl;
        //double mu = chemicalpotential(0.5, temp);
        // Observables_.Get_Non_Interacting_dispersion();
        //Hamiltonian_.Ham_.print();
        Observables_.Calculate_Akw();
        //Observables_.Calculate_Akw_at_w(mu);
        //Observables_.Calculate_Nw();

    }

    else if (ex_string=="sf"){

     cout<<setprecision(9);
     MCEngine MCEngine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);

     Observables_.Initialize();     // Set All Observables to zero

     MCEngine_.RUN_MC();      // Monte-Carlo Engine

     Observables_.Calculate_Nw();

    }

    else if (ex_string=="ve"){ //observe
        double Classical_E, QuantE;
        Classical_E= Hamiltonian_.GetCLEnergy();
        cout <<"Classical_E = "<<Classical_E<<endl;
        Hamiltonian_.InteractionsCreate();
        Hamiltonian_.Diagonalize(Parameters_.Dflag);



        int n_states_occupied_zeroT;
        double initial_mu_guess;
        n_states_occupied_zeroT = Coordinates_.nbasis_*2.0*(Parameters_.Fill/(Parameters_.n_orbs*2.0));

        if(!Parameters_.fixed_mu_value){
           initial_mu_guess = 0.5 * (Hamiltonian_.eigs_[n_states_occupied_zeroT - 1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
        }
        else{
            initial_mu_guess=Parameters_.fixed_mu_value;
        }

        Parameters_.mus = Hamiltonian_.chemicalpotential(initial_mu_guess, Parameters_.Fill);
        QuantE = Hamiltonian_.E_QM();
        cout <<"Quantum_E = "<<QuantE<<endl;

        Observables_.Calculate_Nw();

    }

    else {
        cout <<"Executable not present"<<endl;
    }





    cout << "--------THE END--------" << endl;
} // main

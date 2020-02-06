//
// Created by iskakoff on 19/07/16.
//
#include <iostream>
#include <fstream>

#include <edlib/EDParams.h>
#include "edlib/Hamiltonian.h"
#include "edlib/SzSymmetry.h"
#include "edlib/SOCRSStorage.h"
#include "edlib/CRSStorage.h"
#include "edlib/HubbardModel.h"
#include "edlib/GreensFunction.h"
#include "edlib/ChiLoc.h"
#include "edlib/HDF5Utils.h"
#include "edlib/SpinResolvedStorage.h"
#include "edlib/StaticObservables.h"
#include "edlib/MeshFactory.h"

int main(int argc, const char ** argv) {

  std::ofstream out("ED_output.txt");
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!


  

// Init MPI if enabled
#ifdef USE_MPI
  MPI_Init(&argc, (char ***) &argv);
  alps::mpi::communicator comm;
#endif
// Define and read model parameters
  alps::params params(argc, argv);
  EDLib::define_parameters(params);
  if(params.help_requested(std::cout)) {
    exit(0);
  }
  // open output file
  alps::hdf5::archive ar;
#ifdef USE_MPI
  if(!comm.rank())
#endif
    ar.open(params["OUTPUT_FILE"].as<std::string>(), "w");
// Start calculations
  try {
    // Construct Hamiltonian object
    typedef EDLib::SRSSIAMHamiltonian HType;
#ifdef USE_MPI
    HType ham(params, comm);
#else
    HType ham(params);
#endif
    // Diagonalize Hamiltonian
    ham.diag();

    // Save eigenvalues to HDF5 file
    EDLib::hdf5::save_eigen_pairs(ham, ar, "results");

    // Calculate Static Observables
    EDLib::StaticObservables<HType> so(params);
    std::map<std::string, std::vector<double>> observables = so.calculate_static_observables(ham);
    // so.print_static_observables(ham);
    #ifdef USE_MPI
    if(comm.rank() == 0) {
    #endif      
    EDLib::hdf5::save_static_observables(observables, ar, "results");
    #ifdef USE_MPI
    }
    #endif  
    // Construct Green's function object
    //EDLib::gf::GreensFunction < HType,  alps::gf::matsubara_positive_mesh, alps::gf::statistics::statistics_type> greensFunction(params, ham,alps::gf::statistics::statistics_type::FERMIONIC);
    // Compute and save Green's function
    //greensFunction.compute();
    //greensFunction.save(ar, "results");

    
    // Construct Green's function object
    //EDLib::gf::GreensFunction < HType,  alps::gf::real_frequency_mesh> greensFunctionR(params, ham);
    // Compute and save Green's function
    //greensFunctionR.compute();
    //greensFunctionR.save(ar, "results");


    // Init two particle Green's function object
    EDLib::gf::ChiLoc<HType, alps::gf::real_frequency_mesh> susc(params, ham);
    // Compute and save spin susceptibility
    susc.compute();
    susc.save(ar, "results");
    // Compute and save charge susceptibility
    //susc.compute<EDLib::gf::NOperator<double> >();
    //susc.save(ar, "results");

    #ifdef USE_MPI
    if(comm.rank() == 0) {
    #endif
      std::ofstream end_file;
      end_file.open("./end_file.dat");
      end_file << "ED is done!" << std::endl;
      end_file.close();
      #ifdef USE_MPI
          }
      #endif

  } catch (std::exception & e) {
#ifdef USE_MPI
    if(comm.rank() == 0) {
      std::cerr<<e.what();
      ar.close();
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
#else
    std::cerr<<e.what();
#endif



  }
#ifdef USE_MPI
  if(!comm.rank())
#endif
  ar.close();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}

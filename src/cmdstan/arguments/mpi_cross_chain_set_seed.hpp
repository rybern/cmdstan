#ifndef CMDSTAN_ARGUMENTS_MPI_CROSS_CHAIN_SET_SEED_HPP
#define CMDSTAN_ARGUMENTS_MPI_CROSS_CHAIN_SET_SEED_HPP

#include <cmdstan/arguments/categorical_argument.hpp>
#include <cmdstan/arguments/arg_diagnostic_file.hpp>
#include <cmdstan/arguments/arg_refresh.hpp>

#ifdef STAN_LANG_MPI
#include <stan/math/mpi/envionment.hpp>
#endif

namespace cmdstan {

  void mpi_cross_chain_set_seed(unsigned int& seed, int num_chains) {
#ifdef MPI_ADAPTED_WARMUP
    using stan::math::mpi::Session;
    using stan::math::mpi::Communicator;

    const Communicator& inter_comm = Session::inter_chain_comm(num_chains);
    const Communicator& intra_comm = Session::intra_chain_comm(num_chains);
    MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, MPI_COMM_STAN);
    seed += inter_comm.rank();
    MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, intra_comm.comm());
#endif
  }
}
#endif
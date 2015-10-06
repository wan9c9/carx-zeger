#include<string>
#include<iomanip>
#include<cstdlib>
#include<iostream>
#include<mpi.h>

//usage mpi-wrapper nRep prepare-cmd singleton-cmd final-cmd
//
main(int argc, char* argv[]) {
    int         my_rank;       /* rank of process      */
    int         p;             /* number of processes  */
    int         source;        /* rank of sender       */
    int         dest;          /* rank of receiver     */
    int         tag = 0;       /* tag for messages     */
    char        message[100];  /* storage for message  */
    MPI_Status  status;        /* return status for    */

    long unsigned nRep;
    if(argc <4)
      return 0;

    nRep = std::atoi(argv[1]);

    
    /* receive              */
    /* Start up MPI */
    MPI_Init(&argc, &argv);
    /* Find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(my_rank==0)
    {
      std::cout << "Total number of processes:" << p << std::endl;

      std::cout << "arguments provided:" << std::endl;
      for(unsigned i=0;i<argc;++i)
      {
        std::cout << "Arg"  << i << ": " << argv[i] << std::endl;
      }
			std::cout << "running preparation cmd: " << argv[2] << std::endl;
			system(argv[2]);
			std::cout << "running individual simulations." << std::endl;
    }
    
    for(unsigned i=0;i<nRep;++i)
    {
      if(i%p==my_rank)
      {
        try
        {
          std::string cmd;
          cmd = argv[3] + (" " + std::to_string(i));
          system(cmd.c_str());
					std::cout << "Rep: " << i << std::setw(8) << ",done, by " << my_rank << std::endl;
        }
        catch(...){
					std::cout << "Rep: " << i << std::setw(8) << ",ERROR, by " << my_rank << std::endl;
        }
      }
    }
    /* Shut down MPI */
    MPI_Finalize();
    if(argc>3)
      system(argv[4]);
    return 0;
} /* main */


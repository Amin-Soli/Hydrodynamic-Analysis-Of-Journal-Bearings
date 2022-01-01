#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#  include "mpi.h"
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif 

#include <Sacado.hpp>		
#include <cstdio>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"
#include "../classes/SimpleProblemInterface.H"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int 
main (int argc, char **argv)
{
  clock_t tStart = clock();
  
  srand (time(NULL));
  
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm CommWorld (MPI_COMM_WORLD);
#else
  Epetra_SerialComm CommWorld;
#endif

  if (CommWorld.MyPID() == 0)
    {
#ifdef HAVE_MPI
      Epetra_MpiComm Comm (MPI_COMM_SELF);
#else
      Epetra_SerialComm Comm;
#endif

	  cout.precision(12); 
      cout.setf(ios::fixed);
      
      const char *path="../inputData/bearingProperties.txt";
      
      tiltingPadBearings <F> myFObject(path);
      
      int n_pad = myFObject.padNumber();
      
      Epetra_Map Map ((n_pad + 2), 0, Comm);

      Epetra_Vector InitialGuess (Map);
      Epetra_Vector testInitialGuess (Map);
      
      for(int i = 0; i <n_pad ; i++)
         InitialGuess[i] = (0.002)*((double)rand() / (double)RAND_MAX ) - 0.001;
      
      InitialGuess[n_pad] = 0.0 ;  
      InitialGuess[n_pad + 1] = 0.4*((double)rand() / (double)RAND_MAX ) + 0.2;     
      
      F  *y , x[n_pad + 2];
      
	  for(int i = 0; i <(n_pad + 2); i++)
	  { 
		  x[i]=InitialGuess[i];
		  x[i].diff(i, (n_pad + 2));
	  }
		
      y = myFObject.solve(x);   

      RCP<SimpleProblemInterface> interface = 
	  rcp (new SimpleProblemInterface );
                
      /*******************  Trust Region Based Solver  **************************/
      
     RCP<ParameterList> params = parameterList ("NOX") ;
     params->set ("Nonlinear Solver", "Trust Region Based");
     ParameterList& printParams = params->sublist ("Printing");
     printParams.set ("MyPID", Comm.MyPID ()); 
     printParams.set ("Output Precision", 3);
     printParams.set ("Output Processor", 0);
     const bool verbose = false;
     if (verbose) {
	    printParams.set ("Output Information", 
			 NOX::Utils::OuterIteration + 
			 NOX::Utils::OuterIterationStatusTest + 
			 NOX::Utils::InnerIteration +
			 NOX::Utils::Parameters + 
			 NOX::Utils::Details + 
			 NOX::Utils::Warning);
      } else {
	    printParams.set ("Output Information", NOX::Utils::Warning);
      }   

      ParameterList& dirParams1 = params->sublist ("Direction");
      dirParams1.set ("Method", "Newton");
      
      ParameterList& dirParams2 = params->sublist ("Cauchy Direction");
      dirParams2.set ("Method", "Steepest Descent");
      
      ParameterList& newtonParams = dirParams1.sublist ("Newton");
      newtonParams.set ("Forcing Term Method", "Constant");
     
      ParameterList& CauchyParams = dirParams2.sublist ("Steepest Descent");
      CauchyParams.set ("Scaling Type", "Quadratic Model Min");
      
      ParameterList& variable = params->sublist ("Trust Region");
      variable.set ("Minimum Trust Region Radius", 1.0e-6);
      variable.set ("Maximum Trust Region Radius", 1.0e+10);
      variable.set ("Minimum Improvement Ratio", 1.0e-4);
      variable.set ("Contraction Trigger Ratio", 0.1);
      variable.set ("Expansion Trigger Ratio", 0.75);
      variable.set ("Contraction Factor", 0.25);
      variable.set ("Expansion Factor", 4.0);
      variable.set ("Recovery Step", 1.0);
      variable.set ("Use Ared/Pred Ratio Calculation", false);
      
      ParameterList& lsParams = newtonParams.sublist ("Linear Solver") ;
   
      lsParams.set ("Aztec Solver", "GMRES");  
      lsParams.set ("Max Iterations", 800);  
      lsParams.set ("Tolerance", 1e-4);
      lsParams.set ("Output Frequency", 50);    
      lsParams.set ("Aztec Preconditioner", "ilu"); 
      
      /*********************  Line Search Based Solver  **************************/
      
      RCP<ParameterList> params1 = parameterList ("NOX");

      params1->set ("Nonlinear Solver", "Line Search Based");
   
      ParameterList& printParams1 = params1->sublist ("Printing");
      printParams1.set ("MyPID", Comm.MyPID ()); 
      printParams1.set ("Output Precision", 3);
      printParams1.set ("Output Processor", 0);
      
      if (verbose) {
	    printParams1.set ("Output Information", 
			 NOX::Utils::OuterIteration + 
			 NOX::Utils::OuterIterationStatusTest + 
			 NOX::Utils::InnerIteration +
			 NOX::Utils::Parameters + 
			 NOX::Utils::Details + 
			 NOX::Utils::Warning);
      }
      
      else {
	    printParams1.set ("Output Information", NOX::Utils::Warning);
      }   

      ParameterList& searchParams = params1->sublist ("Line Search");
      searchParams.set ("Method", "Full Step");

      ParameterList& dirParams = params1->sublist ("Direction");
      dirParams.set ("Method", "Newton");

      ParameterList& newtonParams1 = dirParams.sublist ("Newton");
      newtonParams1.set ("Forcing Term Method", "Constant");
      
      ParameterList& lsParams1 = newtonParams1.sublist ("Linear Solver") ;
      
      lsParams1.set ("Aztec Solver", "GMRES");  
      lsParams1.set ("Max Iterations", 800);  
      lsParams1.set ("Tolerance", 1e-4);
      lsParams1.set ("Output Frequency", 50);    
      lsParams1.set ("Aztec Preconditioner", "ilu"); 
      
      /***************************************************************************/
      
      RCP<Epetra_CrsMatrix> A = rcp (new Epetra_CrsMatrix (Copy, Map, (n_pad + 2)));
      {
		std::vector<int> indices(n_pad + 2);
		std::vector<double> values(n_pad + 2);

		for(int i=0;i<(n_pad + 2);i++)
		   indices[i] = i; 
		  
		for(int i=0;i<(n_pad + 2);i++)
		{
		   for(int j=0;j<(n_pad + 2);j++)
			   values[j]=y[i].dx(j);
		   
		   A.get()->InsertGlobalValues (i, (n_pad + 2), &values[0], &indices[0]);
		}

		A.get()->FillComplete();
      }  
      
      delete []y;

      RCP<NOX::Epetra::Interface::Required> iReq = interface;
      RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;

      RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
	  rcp (new NOX::Epetra::LinearSystemAztecOO (printParams, lsParams,
						    iReq, iJac, A, InitialGuess));
						   
	  RCP<NOX::Epetra::LinearSystemAztecOO> linSys1 = 
	  rcp (new NOX::Epetra::LinearSystemAztecOO (printParams1, lsParams1,
						    iReq, iJac, A, InitialGuess));

     /******************* start algorithm by Trust Region Based Solver  **************************/
      
      NOX::Epetra::Vector noxInitGuess (InitialGuess, NOX::DeepCopy);
      RCP<NOX::Epetra::Group> group = 
	  rcp (new NOX::Epetra::Group (printParams, iReq, noxInitGuess, linSys));

      RCP<NOX::StatusTest::NormF> testNormF = 
	  rcp (new NOX::StatusTest::NormF (1.0e-6));

      RCP<NOX::StatusTest::MaxIters> testMaxIters = 
	  rcp (new NOX::StatusTest::MaxIters (1000));

      RCP<NOX::StatusTest::Combo> combo = 
	  rcp (new NOX::StatusTest::Combo (NOX::StatusTest::Combo::OR, 
					 testNormF, testMaxIters));

      RCP<NOX::Solver::Generic> solver = 
	  NOX::Solver::buildSolver (group, combo, params);

      NOX::StatusTest::StatusType status;
      int count_algorithm=0, num_correction=0; 
      cout<<endl;
      bool flag_switch_solver=true;
      bool flag_check_solution=false;
      double sol[n_pad + 2] ;
      double residual_solution=0 ;
      
      /***************************** start solution algorithm **********************************/
     
      cout<<endl<<"W= "<<myFObject.mass()<<endl<<endl;
      cout<<endl<<"initial guess : "<<endl<<endl;
      cout<<InitialGuess<<endl;
      cout<<"/************************************************************************/"<<endl;
     
      do
      {		 
		 cout<<"                              iterate( "<<count_algorithm<<" )                                  "<<endl;
        
         const NOX::Epetra::Group& intermediateGroup = 
               dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
         const Epetra_Vector& intermediateSolution =
               dynamic_cast<const NOX::Epetra::Vector&> (intermediateGroup.getX ()).getEpetraVector (); 
		 
         InitialGuess=intermediateSolution;
         
         checkSolution(InitialGuess,flag_check_solution,n_pad);
               
		 cout<<"solution before correction : "<<endl;
		 cout<<intermediateSolution<<endl;
		 cout<<endl;
		 
		 cout<<"solution after correction : "<<endl;
		 cout<<InitialGuess<<endl;
		 
		 if(flag_check_solution)
		 {	
			num_correction++ ;
			 			   
		    NOX::Epetra::Vector noxSol (InitialGuess, NOX::DeepCopy);
		    
            group = rcp (new NOX::Epetra::Group(printParams, iReq, noxSol, linSys));
	         
            solver = NOX::Solver::buildSolver (group, combo, params);
            
            flag_check_solution=false;           
            
	     }
	     
	     cout<<"number of correction solution at this moment: "<<num_correction;
	         
	     cout<<endl<<"..............................................................................."<<endl;
		 
		 status =solver->step();
		 
	     residual_solution = solver->getPreviousSolutionGroupPtr()->getNormF();
		 cout<<endl<<"2-Norm of Residual = "<< residual_solution <<endl<<endl;
		 
		 const NOX::Epetra::Group& testIntermediateGroup = 
               dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
         const Epetra_Vector& testIntermediateSolution =
               dynamic_cast<const NOX::Epetra::Vector&> (testIntermediateGroup.getX ()).getEpetraVector (); 
		 
         testInitialGuess=testIntermediateSolution;
         
         if( status == NOX::StatusTest::Unconverged && testSolution(InitialGuess,testInitialGuess,count_algorithm,residual_solution,n_pad) )
         {
			 num_correction++ ;
			 
		     NOX::Epetra::Vector noxSol (InitialGuess, NOX::DeepCopy);
		    
             group = rcp (new NOX::Epetra::Group(printParams, iReq, noxSol, linSys));
	         
             solver = NOX::Solver::buildSolver (group, combo, params);           
             
		 }
		 
		 if( status == NOX::StatusTest::Unconverged  && 
		     residual_solution <0 && flag_switch_solver 
		   )
                  
		 {	
			 		 			   
		    NOX::Epetra::Vector noxSol (InitialGuess, NOX::DeepCopy);
		    
            group = rcp (new NOX::Epetra::Group(printParams1, iReq, noxSol, linSys1));
	         
            solver = NOX::Solver::buildSolver (group, combo, params1);
            
            flag_switch_solver=false;
    
	     } 
		 		
		 count_algorithm++;
		 
	  }while(status == NOX::StatusTest::Unconverged);
	  
	   /***************************** end solution algorithm **********************************/
	 
	  cout<<endl<<"..........................  information after converge .........................."<<endl;
	 
	  cout<<endl<<"Nonlinear Iterations = "<<count_algorithm<<endl;
	  residual_solution = solver->getSolutionGroupPtr()->getNormF();
	  cout<<"2-Norm of Residual = "<< residual_solution <<endl<<endl;

      if (Comm.MyPID() == 0) {
	    cout << endl << "-- Parameter List From Solver --" << endl;
	    solver->getList ().print (cout);
      }

      const NOX::Epetra::Group& finalGroup = 
	  dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());

      const Epetra_Vector& finalSolution = 
	  dynamic_cast<const NOX::Epetra::Vector&> (finalGroup.getX ()).getEpetraVector ();
      
      if (Comm.MyPID() == 0) {
	    cout << "Computed solution: " << endl;
      }

      cout << finalSolution;  
      
      for(int i=0;i<(n_pad + 2);i++)
		  sol[i]=finalSolution[i];
		  
	  tiltingPadBearings <double> myDoubleObject(path);
		  
	  myDoubleObject.solve_withPrintPressure(sol);
	  
	  myDoubleObject.dynamicCoeffTotal(sol);
	  
	  cout<<endl<<"..........................................................................."<<endl<<endl;
      
      myDoubleObject.writeAllInformationOfBearing();      
  }

#ifdef HAVE_MPI
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  double time_program =(double)(clock() - tStart)/CLOCKS_PER_SEC;
  int min_program = time_program/60.0 ;
  double sec_program = time_program - min_program*60.0 ;
  cout<<endl<<endl<<"time of run program = "<<min_program<<" min , "<<
       sec_program<<" sec "<<endl;
  
  return EXIT_SUCCESS;
}

#include "tiltingPadBearings.H"
typedef Sacado::Fad::DFad<double> F; 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Classes  ///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SimpleProblemInterface : 
  public NOX::Epetra::Interface::Required,
  public NOX::Epetra::Interface::Jacobian
{
public:

  // The constructor accepts an initial guess and the exact solution
  // vector (which we know because we created the example).  We make
  // deep copies of each.
  SimpleProblemInterface ()
  {}

  // Destructor.
  ~SimpleProblemInterface() {}

  // Compute f := F(x), where x is the input vector and f the output
  // vector.
  bool 
  computeF (const Epetra_Vector & x, 
	    Epetra_Vector & f,
	    NOX::Epetra::Interface::Required::FillType F)
  {
	
	const char *path="../inputData/bearingProperties.txt";
	
	tiltingPadBearings <double> myObject(path);
	
	int n_pad = myObject.padNumber();
	
	double r[n_pad + 2], *forces;
	
	for(int i=0;i<(n_pad + 2);i++)
		r[i]=x[i];
		 
   	forces = myObject.solve(r);
   	
   	for(int i=0;i<(n_pad + 2);i++)
   	  f[i]=forces[i];
      
    delete []forces;  

    return true;
  }

  bool 
  computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac)
  {
	  
	  const char *path="../inputData/bearingProperties.txt";
	
	  tiltingPadBearings <F> myObject(path);
	
	  int n_pad = myObject.padNumber();
	  
	  F  *y , r[n_pad + 2];
	  
	  for(int i = 0; i <(n_pad + 2); i++) 
	  {
		 r[i]=x[i];
		 r[i].diff(i, (n_pad + 2));
	  }
		
    y = myObject.solve(r);
    
    Epetra_CrsMatrix* J = dynamic_cast<Epetra_CrsMatrix*>(&Jac);

    if (J == NULL) {
      std::ostringstream os;
      os << "*** Problem_Interface::computeJacobian() - The supplied "
	 << "Epetra_Operator object is NOT an Epetra_CrsMatrix! ***";
      throw std::runtime_error (os.str());
    }

    std::vector<int> indices(n_pad + 2);
    std::vector<double> values(n_pad + 2);
    
    for(int i=0;i<(n_pad + 2);i++)
      indices[i] = i; 
      
    for(int i=0;i<(n_pad + 2);i++)
    {
       for(int j=0;j<(n_pad + 2);j++)
		   values[j]=y[i].dx(j);
	   
	   J->ReplaceGlobalValues (i,(n_pad + 2), &values[0], &indices[0]);
    }
    
    delete []y;  
		  
    return true;
  }

  bool 
  computePrecMatrix (const Epetra_Vector & x, 
		     Epetra_RowMatrix & M) 
  {
    throw std::runtime_error ("*** SimpleProblemInterface does not implement "
			      "computing an explicit preconditioner from an "
			      "Epetra_RowMatrix ***");
  }  

  bool 
  computePreconditioner (const Epetra_Vector & x, 
			 Epetra_Operator & O)
  {
    throw std::runtime_error ("*** SimpleProblemInterface does not implement "
			      "computing an explicit preconditioner from an "
			      "Epetra_Operator ***");
  }  

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// functions  ///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void  checkSolution (Epetra_Vector &x, bool &flag, int n_pad)
{
    double min=-0.002 , max=0.002 ; 

    for(int i = 0; i <n_pad ; i++)
       if(x[i]<=min || x[i]>=max)
       {
          x[i] = (max-min)*((double)rand() / (double)RAND_MAX ) + min ;
          flag=true;
	   }
	     
    if(x[n_pad]<-0.1 || x[n_pad]>0.1)
    {
       x[n_pad] = 0.0 ;
       flag=true;  
    }    
    
    if(x[n_pad + 1]<0.0 || x[n_pad + 1]>1.0)
    {
       x[n_pad + 1] = 0.4*((double)rand() / (double)RAND_MAX ) + 0.2;
       x[n_pad] = 0.0;
       flag=true;  
    }
   
   return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool testSolution(Epetra_Vector &x,Epetra_Vector &y,int count_algorithm,double residual_solution, int n_pad)
{
	bool flag=false;
	bool flag_count=true;
	static int count=0;

	for(int i = 0; i <(n_pad + 2) ; i++)
	   if(x[i]!=y[i])
	      flag_count=false;
	      
	if(flag_count)
	   count++;
	   
	if(count>1 || (count_algorithm % 50 ==0 && count_algorithm!=0 && residual_solution>100.0) || residual_solution>1.0e5  )
	  {
         for(int i = 0; i <n_pad ; i++)
            x[i] = (0.004)*((double)rand() / (double)RAND_MAX ) - 0.002;
      
            x[n_pad] = 0.0; 
            x[n_pad + 1] = 0.6*((double)rand() / (double)RAND_MAX ) + 0.2;
		 
		 count=0;
		 flag=true;
	  }
	   
	return flag ;
	
}

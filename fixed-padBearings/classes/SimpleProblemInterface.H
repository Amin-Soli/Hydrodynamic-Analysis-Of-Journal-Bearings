#include "fixedPadBearings.H"
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
	
	fixedPadBearings <double> myObject(path);
	
	double r[2], *forces;
	
	for(int i=0;i<2;i++)
		r[i]=x[i];
		 
   	forces = myObject.solve(r);
   	
   	for(int i=0;i<2;i++)
   	  f[i]=forces[i];
      
    delete []forces;

    return true;
  }

  bool 
  computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac)
  {
	
	const char *path="../inputData/bearingProperties.txt";
	
	fixedPadBearings <F> myObject(path);
	
	F  *y , r[2];
	
	for(int i = 0; i < 2; i++) 
	{
		r[i]=x[i];
		r[i].diff(i, 2);
	}
		
    y = myObject.solve(r);
    
    Epetra_CrsMatrix* J = dynamic_cast<Epetra_CrsMatrix*>(&Jac);

    if (J == NULL) {
      std::ostringstream os;
      os << "*** Problem_Interface::computeJacobian() - The supplied "
	 << "Epetra_Operator object is NOT an Epetra_CrsMatrix! ***";
      throw std::runtime_error (os.str());
    }

    std::vector<int> indices(2);
    std::vector<double> values(2);
    
    for(int i=0;i<2;i++)
      indices[i] = i; 
      
    for(int i=0;i<2;i++)
    {
       for(int j=0;j<2;j++)
		   values[j]=y[i].dx(j);
	   
	   J->ReplaceGlobalValues (i, 2, &values[0], &indices[0]);
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

void  checkSolution (Epetra_Vector &x, bool &flag )
{
   double pi=acos(-1);
   double min=0 , max=pi ; 
   
   if(x[0]<=0 || x[0]>=1)
   {
	 x[0]= 0.4*((double)rand() / (double)RAND_MAX ) + 0.2; 
	 flag=true;
   }
   
   if(x[1]<=min || x[1]>=max)
   {
	 x[1]= (max-min)*((double)rand() / (double)RAND_MAX ) + min; 
	 flag=true;
   }
   
   return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool testSolution(Epetra_Vector &x,Epetra_Vector &y)
{
	double pi=acos(-1);
	bool flag=false;
	static int count=0;
	
	if(x[0]==y[0] && x[1]==y[1])
	   count++;
	   
	if(count>1)
	  {
		 x[0]= 0.4*((double)rand() / (double)RAND_MAX ) + 0.2; ;
		 x[1]= pi*((double)rand() / (double)RAND_MAX ) ;
		 count=0;
		 flag=true;
	  }
	   
	return flag ;
	
}

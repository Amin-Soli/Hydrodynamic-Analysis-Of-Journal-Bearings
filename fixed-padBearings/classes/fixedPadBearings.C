#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>		
#include <cstdio>
#include <string>
using namespace std;

/********************************** constructor ******************************************/

template <class T>
fixedPadBearings<T>::fixedPadBearings (const char *path)
{	 
   readOfBearingProperties(path);  
   
   informationOfSolution = new T[4];
   
   informationOfDynamicCoeff_dimensionLess = new double[8];
   
   pi=acos(-1);
   
   for(int i=0;i<numberOfGrooves;i++)
   {
	   centerGrooves[i]=centerGrooves[i]*pi/180.0;
	   magnitudGrooves[i]=magnitudGrooves[i]*pi/180.0;
   }	
    
   u = N*2*pi*r/60.0 ;
    
   W = -1.0*(1.0*r/c)*(1.0*r/c)*v*2*r*L*N/(sommerfeld*60) ;
    
   dx=2*pi/(n-1);
   
   dy=L/(m-1); 
	
   p_in =0 ; 

   
}

/*****************************************************************************************/
/******************************* member functions ****************************************/
/*****************************************************************************************/

template <class T>
void fixedPadBearings<T>::solve_withPrintPressure (T *x)
{
     int i,j,count_groove;
	 T  e_p, cop_p, error_p , force_y=0 , force_x=0 ;
     T *functions;
     T **p;
     bool flag_grooves;
     p = new T *[m];
     for(i = 0; i <m; i++)
     p[i] = new T[n];
     
	 functions=new T[2];
	 
	 for (j=0;j<m;j++)
       for (i=0;i<n;i++)
          p[j][i]=p_in;

     if(typeOfBearing=="Elliptical")
     {
		 phi_L=atan( x[0]*sin(x[1])/(x[0]*cos(x[1])+ preload));
		 phi_U=atan( x[0]*sin(x[1])/(x[0]*cos(x[1])- preload));
		 epsilon_L=sqrt(x[0]*x[0] + preload*preload - 2*x[0]*preload*cos(pi-x[1]));
		 epsilon_U=sqrt(x[0]*x[0] + preload*preload - 2*x[0]*preload*cos(-1*x[1]));
	 }

   do
    {
        e_p=0;

      for (j=0;j<m;j++)
		for (i=0;i<n;i++)
        {
            if(j==0 || j==m-1)
                p[j][i]=p_in;


            else
            {

               flag_grooves=true;

               for(count_groove=0;(count_groove<numberOfGrooves && flag_grooves);count_groove++)
                  if( i*dx>=(centerGrooves[count_groove] - magnitudGrooves[count_groove]/2.0) &&
                      i*dx<=(centerGrooves[count_groove] + magnitudGrooves[count_groove]/2.0)    ) 
                  {
                      p[j][i]=p_in;
                      flag_grooves=false;
                  }
                  
                  

               if(i==0 && flag_grooves)
               {

                   p[j][i]=(
                              ( pow(h(i+1,x),3)/(12.0*v )- pow(h(n-2,x),3)/(12.0*v ) )*(p[j][i+1]-p[j][n-2])/(4.0*dx*dx)
                              + pow(h(i,x),3)*(p[j][i+1]+ p[j][n-2])/(12.0*v*dx*dx )
                              + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                              +(pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                              - u*r*(h(i+1,x)-h(n-2,x))/(4.0*dx)
                           )
                           /
                           (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                   if(p[j][i]<0)
                       p[j][i]=p_in;

               }


               else if(i==n-1 && flag_grooves)
                   p[j][i]=p[j][0];


               else if(i==n-2 && flag_grooves)
               {
                   p[j][i]=(
                             ( pow(h(0,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p[j][0]-p[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p[j][0]+ p[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                             +(pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                             - u*r*(h(0,x)-h(i-1,x))/(4.0*dx)
                           )
                           /
                           (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                    if(p[j][i]<0)
                        p[j][i]=p_in;

               }

               else if(flag_grooves)
               {

                  cop_p=p[j][i];

                  p[j][i]=(
                            ( pow(h(i+1,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                            + pow(h(i,x),3)*(p[j][i+1]+ p[j][i-1])/(12.0*v*dx*dx )
                            + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                            +(pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            - u*r*(h(i+1,x)-h(i-1,x))/(4.0*dx)
                          )
                          /
                          (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );


                  if(p[j][i]<0)
                      p[j][i]=p_in;

                  error_p=fabs(cop_p-p[j][i]);

                  if(error_p>e_p)
                     e_p=error_p;                  

	           }


            }

        }

	}while(e_p>errorLimit);
	
	string filename= "pressureDistributionOnMidplaneOfBearing.txt";
	ofstream file(filename.c_str());
	file<<"theta (degree)"<<'\t'<<'\t'<<"presuure (Pa)"<<endl<<endl;
	for(j=0;j<n;j++)
	         file<<j*dx*180.0/pi<<'\t'<<'\t'<<'\t'<<p[m/2][j]<<endl;
	file.close();

	for (j=0;j<m;j++)
        for(i=0;i<n;i++)
        {
			force_x=force_x+(r*dx*dy*p[j][i]*cos(i*dx) );
            force_y=force_y+(r*dx*dy*p[j][i]*sin(i*dx) ); 
		} 
		
	functions[0]= force_x - W ; 
    
    functions[1]= force_y ; 
    
    informationOfSolution[0]= x[0];
    informationOfSolution[1]= x[1];
    
    for(i=0;i<2;i++)
	     informationOfSolution[i+2]=functions[i];	   
	
	 for(i = 0; i < m; ++i) 
         delete[] p[i];   
     
     delete[] p;
	 
	 delete[] functions; 
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
T* fixedPadBearings<T>::solve (T *x)
{
     int i,j,count_groove;
	 T  e_p, cop_p, error_p , force_y=0 , force_x=0 ;
     T *functions;
     T **p;
     bool flag_grooves;
     p = new T *[m];
     for(i = 0; i <m; i++)
     p[i] = new T[n];
     
	 functions=new T[2];
	 
	 for (j=0;j<m;j++)
       for (i=0;i<n;i++)
          p[j][i]=p_in;
          
     if(typeOfBearing=="Elliptical")
     {
		 phi_L=atan( x[0]*sin(x[1])/(x[0]*cos(x[1])+ preload));
		 phi_U=atan( x[0]*sin(x[1])/(x[0]*cos(x[1])- preload));
		 epsilon_L=sqrt(x[0]*x[0] + preload*preload - 2*x[0]*preload*cos(pi-x[1]));
		 epsilon_U=sqrt(x[0]*x[0] + preload*preload - 2*x[0]*preload*cos(-1*x[1]));
	 }

   do
   {
      e_p=0;

      for (j=0;j<m;j++)
		for (i=0;i<n;i++)
        {
            if(j==0 || j==m-1)
                p[j][i]=p_in;


            else
            {

               flag_grooves=true;

               for(count_groove=0;(count_groove<numberOfGrooves && flag_grooves);count_groove++)
                  if( i*dx>=(centerGrooves[count_groove] - magnitudGrooves[count_groove]/2.0) &&
                      i*dx<=(centerGrooves[count_groove] + magnitudGrooves[count_groove]/2.0)    ) 
                  {
                      p[j][i]=p_in;
                      flag_grooves=false;
                  }
                  
                  

               if(i==0 && flag_grooves)
               {

                   p[j][i]=(
                              ( pow(h(i+1,x),3)/(12.0*v )- pow(h(n-2,x),3)/(12.0*v ) )*(p[j][i+1]-p[j][n-2])/(4.0*dx*dx)
                              + pow(h(i,x),3)*(p[j][i+1]+ p[j][n-2])/(12.0*v*dx*dx )
                              + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                              +(pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                              - u*r*(h(i+1,x)-h(n-2,x))/(4.0*dx)
                           )
                           /
                           (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                   if(p[j][i]<0)
                       p[j][i]=p_in;

               }


               else if(i==n-1 && flag_grooves)
                   p[j][i]=p[j][0];


               else if(i==n-2 && flag_grooves)
               {
                   p[j][i]=(
                             ( pow(h(0,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p[j][0]-p[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p[j][0]+ p[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                             +(pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                             - u*r*(h(0,x)-h(i-1,x))/(4.0*dx)
                           )
                           /
                           (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                    if(p[j][i]<0)
                        p[j][i]=p_in;

               }

               else if(flag_grooves)
               {

                  cop_p=p[j][i];

                  p[j][i]=(
                            ( pow(h(i+1,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                            + pow(h(i,x),3)*(p[j][i+1]+ p[j][i-1])/(12.0*v*dx*dx )
                            + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                            +(pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            - u*r*(h(i+1,x)-h(i-1,x))/(4.0*dx)
                          )
                          /
                          (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );


                  if(p[j][i]<0)
                      p[j][i]=p_in;

                  error_p=fabs(cop_p-p[j][i]);

                  if(error_p>e_p)
                     e_p=error_p;                  

	           }


            }

        }

	}while(e_p>errorLimit);

	for (j=0;j<m;j++)
        for(i=0;i<n;i++)
        {
			force_x=force_x+(r*dx*dy*p[j][i]*cos(i*dx) );
            force_y=force_y+(r*dx*dy*p[j][i]*sin(i*dx) ); 
		} 
		
	functions[0]= force_x - W ; 
    
    functions[1]= force_y ;    
	
	 for(i = 0; i < m; ++i) 
         delete[] p[i];   
     
     delete[] p;
	 
	 return functions; 
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
void fixedPadBearings<T>::writeAllInformationOfBearing()
{
   cout<<"\n"<<"bearing of properties : "<<"\n"<<"\n";
   cout << "type of bearing = " << typeOfBearing << "\n";
   cout << "value of clearance = " << c << "\n";
   cout << "value of length = " << L << "\n";
   cout << "value of radius = " << r << "\n";
   cout << "value of viscosity = " << v << "\n";
   cout << "value of rotationOfSpeed (rpm) = " << N << "\n";	
   cout << "value of sommerfeld = " << sommerfeld << "\n";	
   cout << "value of preload = " << preload << "\n";
   cout << "value of maximum error for calculating pressure = " << errorLimit << "\n";
   cout << "is grooves = " << isGrooves << "\n"; 
   
   if(isGrooves == "Yes")
   {
	   cout << "number of grooves = " << numberOfGrooves << "\n";
	   for(int i=0;i<numberOfGrooves;i++)
	   {
		  cout << "center of groove("<<(i+1)<<") = " << (centerGrooves[i]*180.0/pi) << "\n";
		  cout << "magnitud of groove("<<(i+1)<<") = " << (magnitudGrooves[i]*180.0/pi) << "\n";
	   }
	   
   }
     	
   cout << "value of numberOfNodsInTetaDirection = " << n << "\n";
   cout << "value of numberOfNodsInLengthDirection = " << m << "\n";   
   cout<<"\n"<<"...................................................................."<<"\n";
   
   cout<<"\n"<<"all information of solution : "<<"\n"<<"\n";
   cout << "epsilon = " << informationOfSolution[0] << "\n";
   cout << "loading angle = " << informationOfSolution[1]*180.0/pi << "\n";
   cout << "shaftWeight = " << W << "\n";
   cout << "Fy - shaftWeight = " << informationOfSolution[2] << "\n";
   cout << "Fx = " << informationOfSolution[3] << "\n";       
   cout<<"\n"<<"...................................................................."<<"\n";
  
   cout<<"\n"<<"DimensionLess Of Dynamic Coeffs : "<<"\n"<<"\n";
   cout << "Kxx = " << informationOfDynamicCoeff_dimensionLess[0] << "\n";
   cout << "Kxy = " << informationOfDynamicCoeff_dimensionLess[1] << "\n";
   cout << "Kyx = " << informationOfDynamicCoeff_dimensionLess[2] << "\n";
   cout << "Kyy = " << informationOfDynamicCoeff_dimensionLess[3] << "\n";
   cout << "Cxx = " << informationOfDynamicCoeff_dimensionLess[4] << "\n";
   cout << "Cxy = " << informationOfDynamicCoeff_dimensionLess[5] << "\n";
   cout << "Cyx = " << informationOfDynamicCoeff_dimensionLess[6] << "\n";
   cout << "Cyy = " << informationOfDynamicCoeff_dimensionLess[7] << "\n" << "\n";

}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
void fixedPadBearings<T>::readOfBearingProperties(const char *path)
{
    ifstream file(path);

	string searchForClearance = "clearance";
	string searchForLength = "length";
	string searchForRadius = "radius";
	string searchForViscosity = "viscosity";
	string searchForTypeOfBearing = "typeOfBearing";
	string searchForRotationOfSpeed = "rotationOfSpeed";
	string searchForSommerfeld = "sommerfeld";
	string searchForPreload = "preload";
	string searchForIsGrooves = "isGrooves";
	string searchForNumberOfNodsInTetaDirection = "numberOfNodsInTetaDirection";
	string searchForNumberOfNodsInLengthDirection = "numberOfNodsInLengthDirection";
	string searchForNumberOfGrooves = "numberOfGrooves";
	string searchForErrorLimit = "maximumErrorForCalculatingPressure";
		

	string lineOfText;

	for(;;)
	{	
		getline(file, lineOfText);
			  
		if (file.eof()) break;
		
		if (lineOfText.find(searchForClearance, 0) != string::npos)
		  file >> c;
		  	  
		if (lineOfText.find(searchForLength, 0) != string::npos)
		  file >> L;      	  	  
		
		if (lineOfText.find(searchForRadius, 0) != string::npos)
		  file >> r;		  	  	  
		
		if (lineOfText.find(searchForViscosity, 0) != string::npos)
		  file >> v;
		 
		if (lineOfText.find(searchForTypeOfBearing, 0) != string::npos)
		  file >> typeOfBearing;

		if (lineOfText.find(searchForRotationOfSpeed, 0) != string::npos)
		  file >> N;
		
		if (lineOfText.find(searchForSommerfeld, 0) != string::npos)
		  file >> sommerfeld;
		
		if (lineOfText.find(searchForPreload, 0) != string::npos)
		  file >> preload;
		
		if (lineOfText.find(searchForNumberOfNodsInTetaDirection, 0) != string::npos)
		  file >> n;		  	  
		
		if (lineOfText.find(searchForNumberOfNodsInLengthDirection, 0) != string::npos)
		  file >> m;
		  
		if (lineOfText.find(searchForIsGrooves, 0) != string::npos)
		  file >> isGrooves;
		  
		if (lineOfText.find(searchForNumberOfGrooves, 0) != string::npos)
		  file >> numberOfGrooves;
		  
		if (lineOfText.find(searchForErrorLimit, 0) != string::npos)
		  file >> errorLimit;
    }
    
    file.close();
		  
		if(isGrooves=="Yes")
		{	
			magnitudGrooves = new double[numberOfGrooves];
			centerGrooves = new double[numberOfGrooves];
		    
			string searchForCenterGrooves;
			string searchForMagnitudGrooves;
			
			for(int i=0;i<numberOfGrooves;i++)
			{
				stringstream ss;
				ss<<(i+1);
				string str= ss.str();
			
				searchForCenterGrooves = "centerGroove(" + str + ")";
				searchForMagnitudGrooves = "magnitudGroove(" + str + ")"; 
			
				ifstream file(path);
				string lineOfText;
			
			    for(;;)
			    {
					getline(file, lineOfText);
			  
					if (file.eof()) break;
				
					if (lineOfText.find(searchForCenterGrooves, 0) != string::npos)
						file >> centerGrooves[i];	
	    
					if (lineOfText.find(searchForMagnitudGrooves, 0) != string::npos)
						file >> magnitudGrooves[i];
			    }
			    
			    file.close();	

			}
		
		}			  	  							
	
	
	
	if(isGrooves=="No")
	   numberOfGrooves=0;
	   
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
double fixedPadBearings<T>::mass()
{
	return W;
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
T fixedPadBearings<T>::h(int i, T *x)
{
	if(typeOfBearing=="Elliptical")
    {
        if( i*dx<=pi/2.0 || i*dx>=3*pi/2.0 )  
           return c*(1+ epsilon_U*cos(i*dx-pi-phi_U) ) ;

        else
           return c*(1 + epsilon_L*cos(i*dx-phi_L) ) ;
    }

    if(typeOfBearing=="Plain")
        return c*(1 + x[0]*cos(i*dx- x[1]));
	
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
void fixedPadBearings<T>::dynamicCoeffTotal(double *x)
{
	 int i,j,count_groove;
	 double  Kxx=0, Kxy=0, Kyx=0, Kyy=0, Cxx=0, Cxy=0, Cyx=0, Cyy=0 ,Fx=0 , 
	         e_p, cop_p, error_p ,e_p_x, cop_p_x, error_p_x ,
	         e_p_y, cop_p_y, error_p_y ,e_p_dx ,cop_p_dx,
	         error_p_dx ,e_p_dy , cop_p_dy, error_p_dy ;
	          
     double *dyCoeff;
     dyCoeff=new double[8];
     double **p, **p_x, **p_y, **p_dx, **p_dy ;
     bool flag_grooves;
     p = new double *[m];
     p_x = new double *[m];
     p_y = new double *[m];
     p_dx = new double *[m];
     p_dy = new double *[m];
     
     for(i = 0; i <m; i++)
     {
        p[i] = new double[n];
        p_x[i] = new double[n];
        p_y[i] = new double[n];
        p_dx[i] = new double[n];
        p_dy[i] = new double[n];
	 }
	
     if(typeOfBearing=="Elliptical")
     {
		 phi_L=atan( x[0]*sin(x[1])/(x[0]*cos(x[1])+ preload));
		 phi_U=atan( x[0]*sin(x[1])/(x[0]*cos(x[1])- preload));
		 epsilon_L=sqrt(x[0]*x[0] + preload*preload - 2*x[0]*preload*cos(pi-x[1]));
		 epsilon_U=sqrt(x[0]*x[0] + preload*preload - 2*x[0]*preload*cos(-1*x[1]));
	 }
     
     for (j=0;j<m;j++)
       for (i=0;i<n;i++)
       {
          p[j][i]=p_in;
          p_x[j][i]=p_in;
          p_y[j][i]=p_in;
          p_dx[j][i]=p_in;
          p_dy[j][i]=p_in;
	   }
	   
    do
    {
        e_p=0;        
        e_p_x=0;      
        e_p_y=0; 
        e_p_dx=0; 
        e_p_dy=0;
        
      for (j=0;j<m;j++)
		for (i=0;i<n;i++)
        {
            if(j==0 || j==m-1)
            {
                p[j][i]=p_in;
                p_x[j][i]=p_in;
                p_y[j][i]=p_in;
                p_dx[j][i]=p_in;
                p_dy[j][i]=p_in;
			}


            else
            {

               flag_grooves=true;

               for(count_groove=0;(count_groove<numberOfGrooves && flag_grooves);count_groove++)
                  if( i*dx>=(centerGrooves[count_groove] - magnitudGrooves[count_groove]/2.0) &&
                      i*dx<=(centerGrooves[count_groove] + magnitudGrooves[count_groove]/2.0)    ) 
                  {
                      p[j][i]=p_in;
                      p_x[j][i]=p_in;
                      p_y[j][i]=p_in;
                      p_dx[j][i]=p_in;
                      p_dy[j][i]=p_in;
                      
                      flag_grooves=false;
                  }
                  
                  

               if(i==0 && flag_grooves)
               {

                   p[j][i]=(
                              ( pow(h(i+1,x),3)/(12.0*v )- pow(h(n-2,x),3)/(12.0*v ) )*(p[j][i+1]-p[j][n-2])/(4.0*dx*dx)
                              + pow(h(i,x),3)*(p[j][i+1]+ p[j][n-2])/(12.0*v*dx*dx )
                              + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                              + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                              - u*r*(h(i+1,x)-h(n-2,x))/(4.0*dx)
                           )
                           /
                           (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                  
                  p_x[j][i]=(
                             ( pow(h(i+1,x),3)/(12.0*v )- pow(h(n-2,x),3)/(12.0*v ) )*(p_x[j][i+1]-p_x[j][n-2])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p_x[j][i+1]+ p_x[j][n-2])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p_x[j+1][i]+ p_x[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_x[j+1][i]-p_x[j-1][i])/(4.0*dy*dy)
                             + u*r*sin(i*dx)/2.0 + ( 3*cos((i+1)*dx)*pow(h(i+1,x),2)/(12.0*v ) - 3*cos((n-2)*dx)*pow(h(n-2,x),2)/(12.0*v ) ) * (p[j][i+1]-p[j][n-2])/(4.0*dx*dx)
                             + 3*cos(i*dx)*pow(h(i,x),2)*(p[j][n-2]-2*p[j][i]+p[j][i+1])/(12.0*v*dx*dx )
                             + 3*r*r*cos(i*dx)*pow(h(i,x),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             + (3*cos(i*dx)*pow(h(i,x),2)/(12.0*v ) - 3*cos(i*dx)*pow(h(i,x),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                  
                  p_y[j][i]=(
                             ( pow(h(i+1,x),3)/(12.0*v )- pow(h(n-2,x),3)/(12.0*v ) )*(p_y[j][i+1]-p_y[j][n-2])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p_y[j][i+1]+ p_y[j][n-2])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p_y[j+1][i]+ p_y[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_y[j+1][i]-p_y[j-1][i])/(4.0*dy*dy)
                             - u*r*cos(i*dx)/2.0 + ( 3*sin((i+1)*dx)*pow(h(i+1,x),2)/(12.0*v ) - 3*sin((n-2)*dx)*pow(h(n-2,x),2)/(12.0*v ) ) * (p[j][i+1]-p[j][n-2])/(4.0*dx*dx)
                             + 3*sin(i*dx)*pow(h(i,x),2)*(p[j][n-2]-2*p[j][i]+p[j][i+1])/(12.0*v*dx*dx )
                             + 3*r*r*sin(i*dx)*pow(h(i,x),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             + (3*sin(i*dx)*pow(h(i,x),2)/(12.0*v ) - 3*sin(i*dx)*pow(h(i,x),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                  
                  
                  p_dx[j][i]=(
                               ( pow(h(i+1,x),3)/(12.0*v )- pow(h(n-2,x),3)/(12.0*v ) )*(p_dx[j][i+1]-p_dx[j][n-2])/(4.0*dx*dx)
                               + pow(h(i,x),3)*(p_dx[j][i+1]+ p_dx[j][n-2])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x),3)*(p_dx[j+1][i]+ p_dx[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_dx[j+1][i]-p_dx[j-1][i])/(4.0*dy*dy)
                               - r*r*cos(i*dx)
                              )
                              /
                              (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                              
                  
                  p_dy[j][i]=(
                               ( pow(h(i+1,x),3)/(12.0*v )- pow(h(n-2,x),3)/(12.0*v ) )*(p_dy[j][i+1]-p_dy[j][n-2])/(4.0*dx*dx)
                               + pow(h(i,x),3)*(p_dy[j][i+1]+ p_dy[j][n-2])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x),3)*(p_dy[j+1][i]+ p_dy[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_dy[j+1][i]-p_dy[j-1][i])/(4.0*dy*dy)
                               - r*r*sin(i*dx)
                              )
                              /
                              (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                              
                              
                  if(p[j][i]<0)
                  {
                     p[j][i]=p_in;
                     p_x[j][i]=p_in;
                     p_y[j][i]=p_in;
                     p_dx[j][i]=p_in;
                     p_dy[j][i]=p_in;
                  }

               }


               else if(i==n-1 && flag_grooves)
               {
                   p[j][i]=p[j][0];
                   p_x[j][i]=p_x[j][0];
                   p_y[j][i]=p_y[j][0];
                   p_dx[j][i]=p_dx[j][0];
                   p_dy[j][i]=p_dy[j][0];
               }


               else if(i==n-2 && flag_grooves)
               {
                   p[j][i]=(
                             ( pow(h(0,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p[j][0]-p[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p[j][0]+ p[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                             +(pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                             - u*r*(h(0,x)-h(i-1,x))/(4.0*dx)
                           )
                           /
                           (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                 
                 
                  p_x[j][i]=(
                             ( pow(h(0,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_x[j][0]-p_x[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p_x[j][0]+ p_x[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p_x[j+1][i]+ p_x[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_x[j+1][i]-p_x[j-1][i])/(4.0*dy*dy)
                             + u*r*sin(i*dx)/2.0 + ( 3*cos((0)*dx)*pow(h(0,x),2)/(12.0*v ) - 3*cos((i-1)*dx)*pow(h(i-1,x),2)/(12.0*v ) ) * (p[j][0]-p[j][i-1])/(4.0*dx*dx)
                             + 3*cos(i*dx)*pow(h(i,x),2)*(p[j][i-1]-2*p[j][i]+p[j][0])/(12.0*v*dx*dx )
                             + 3*r*r*cos(i*dx)*pow(h(i,x),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             + (3*cos(i*dx)*pow(h(i,x),2)/(12.0*v ) - 3*cos(i*dx)*pow(h(i,x),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                  
                  p_y[j][i]=(
                             ( pow(h(0,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_y[j][0]-p_y[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p_y[j][0]+ p_y[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p_y[j+1][i]+ p_y[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_y[j+1][i]-p_y[j-1][i])/(4.0*dy*dy)
                             - u*r*cos(i*dx)/2.0 + ( 3*sin((0)*dx)*pow(h(0,x),2)/(12.0*v ) - 3*sin((i-1)*dx)*pow(h(i-1,x),2)/(12.0*v ) ) * (p[j][0]-p[j][i-1])/(4.0*dx*dx)
                             + 3*sin(i*dx)*pow(h(i,x),2)*(p[j][i-1]-2*p[j][i]+p[j][0])/(12.0*v*dx*dx )
                             + 3*r*r*sin(i*dx)*pow(h(i,x),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             + (3*sin(i*dx)*pow(h(i,x),2)/(12.0*v ) - 3*sin(i*dx)*pow(h(i,x),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                  
                  
                  p_dx[j][i]=(
                               ( pow(h(0,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_dx[j][0]-p_dx[j][i-1])/(4.0*dx*dx)
                               + pow(h(i,x),3)*(p_dx[j][0]+ p_dx[j][i-1])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x),3)*(p_dx[j+1][i]+ p_dx[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_dx[j+1][i]-p_dx[j-1][i])/(4.0*dy*dy)
                               - r*r*cos(i*dx)
                              )
                              /
                              (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                              
                  
                  p_dy[j][i]=(
                               ( pow(h(0,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_dy[j][0]-p_dy[j][i-1])/(4.0*dx*dx)
                               + pow(h(i,x),3)*(p_dy[j][0]+ p_dy[j][i-1])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x),3)*(p_dy[j+1][i]+ p_dy[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_dy[j+1][i]-p_dy[j-1][i])/(4.0*dy*dy)
                               - r*r*sin(i*dx)
                              )
                              /
                              (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                              
                 
                  if(p[j][i]<0)
                  {
                     p[j][i]=p_in;
                     p_x[j][i]=p_in;
                     p_y[j][i]=p_in;
                     p_dx[j][i]=p_in;
                     p_dy[j][i]=p_in;
                  }

               }

               else if(flag_grooves)
               {

                  cop_p=p[j][i];
                  cop_p_x=p_x[j][i];
                  cop_p_y=p_y[j][i];
                  cop_p_dx=p_dx[j][i];
                  cop_p_dy=p_dy[j][i];

                  p[j][i]=(
                            ( pow(h(i+1,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                            + pow(h(i,x),3)*(p[j][i+1]+ p[j][i-1])/(12.0*v*dx*dx )
                            + r*r*pow(h(i,x),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                            + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            - u*r*(h(i+1,x)-h(i-1,x))/(4.0*dx)
                          )
                          /
                          (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                 
                  
                  p_x[j][i]=(
                             ( pow(h(i+1,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_x[j][i+1]-p_x[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p_x[j][i+1]+ p_x[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p_x[j+1][i]+ p_x[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_x[j+1][i]-p_x[j-1][i])/(4.0*dy*dy)
                             + u*r*sin(i*dx)/2.0 + ( 3*cos((i+1)*dx)*pow(h(i+1,x),2)/(12.0*v ) - 3*cos((i-1)*dx)*pow(h(i-1,x),2)/(12.0*v ) ) * (p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                             + 3*cos(i*dx)*pow(h(i,x),2)*(p[j][i-1]-2*p[j][i]+p[j][i+1])/(12.0*v*dx*dx )
                             + 3*r*r*cos(i*dx)*pow(h(i,x),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             +(3*cos(i*dx)*pow(h(i,x),2)/(12.0*v ) - 3*cos(i*dx)*pow(h(i,x),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );

                  
                  p_y[j][i]=(
                             ( pow(h(i+1,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_y[j][i+1]-p_y[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x),3)*(p_y[j][i+1]+ p_y[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x),3)*(p_y[j+1][i]+ p_y[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_y[j+1][i]-p_y[j-1][i])/(4.0*dy*dy)
                             - u*r*cos(i*dx)/2.0 + ( 3*sin((i+1)*dx)*pow(h(i+1,x),2)/(12.0*v ) - 3*sin((i-1)*dx)*pow(h(i-1,x),2)/(12.0*v ) ) * (p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                             + 3*sin(i*dx)*pow(h(i,x),2)*(p[j][i-1]-2*p[j][i]+p[j][i+1])/(12.0*v*dx*dx )
                             + 3*r*r*sin(i*dx)*pow(h(i,x),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             + (3*sin(i*dx)*pow(h(i,x),2)/(12.0*v ) - 3*sin(i*dx)*pow(h(i,x),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                  
                  
                  p_dx[j][i]=(
                               ( pow(h(i+1,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_dx[j][i+1]-p_dx[j][i-1])/(4.0*dx*dx)
                               + pow(h(i,x),3)*(p_dx[j][i+1]+ p_dx[j][i-1])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x),3)*(p_dx[j+1][i]+ p_dx[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_dx[j+1][i]-p_dx[j-1][i])/(4.0*dy*dy)
                               - r*r*cos(i*dx)
                              )
                              /
                              (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                              
                  
                  p_dy[j][i]=(
                               ( pow(h(i+1,x),3)/(12.0*v )- pow(h(i-1,x),3)/(12.0*v ) )*(p_dy[j][i+1]-p_dy[j][i-1])/(4.0*dx*dx)
                               + pow(h(i,x),3)*(p_dy[j][i+1]+ p_dy[j][i-1])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x),3)*(p_dy[j+1][i]+ p_dy[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x),3)/(12.0*v )- pow(h(i,x),3)/(12.0*v ) )*r*r*(p_dy[j+1][i]-p_dy[j-1][i])/(4.0*dy*dy)
                               - r*r*sin(i*dx)
                              )
                              /
                              (2.0*pow(h(i,x),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x),3)/(12.0*v*dy*dy) );
                              
                                          
                  if(p[j][i]<0)
                  {
                     p[j][i]=p_in;
                     p_x[j][i]=p_in;
                     p_y[j][i]=p_in;
                     p_dx[j][i]=p_in;
                     p_dy[j][i]=p_in;
                  }

                  error_p=fabs(cop_p-p[j][i]);                  
                  error_p_x=fabs(cop_p_x-p_x[j][i]);                
                  error_p_y=fabs(cop_p_y-p_y[j][i]);                 
                  error_p_dx=fabs(cop_p_dx-p_dx[j][i]);                 
                  error_p_dy=fabs(cop_p_dy-p_dy[j][i]);


                  if(error_p>e_p)
                     e_p=error_p;
                     
                  if(error_p_x>e_p_x)
                     e_p_x=error_p_x;
                     
                  if(error_p_y>e_p_y)
                     e_p_y=error_p_y;
                     
                  if(error_p_dx>e_p_dx)
                     e_p_dx=error_p_dx;
                     
                  if(error_p_dy>e_p_dy)
                     e_p_dy=error_p_dy; 

	           }


            }

        }

	}while(      
	        e_p>errorLimit || e_p_x>errorLimit  || e_p_y>errorLimit || 
	        e_p_dx>errorLimit  || e_p_dy>errorLimit   
	      );

    for (j=0;j<m;j++)
      for(i=0;i<n;i++)
      {
		  Fx=Fx + (r*dx*dy*p[j][i]*cos(i*dx) );
          Kxx=Kxx + (r*dx*dy*p_x[j][i]*cos(i*dx) );
          Kyx=Kyx + (r*dx*dy*p_x[j][i]*sin(i*dx) );
          Kxy=Kxy + (r*dx*dy*p_y[j][i]*cos(i*dx) );
          Kyy=Kyy + (r*dx*dy*p_y[j][i]*sin(i*dx) );
          Cxx=Cxx + (r*dx*dy*p_dx[j][i]*cos(i*dx) );
          Cxy=Cxy + (r*dx*dy*p_dx[j][i]*sin(i*dx) );
          Cyx=Cyx + (r*dx*dy*p_dy[j][i]*cos(i*dx) );
          Cyy=Cyy + (r*dx*dy*p_dy[j][i]*sin(i*dx) );
      } 

      dyCoeff[0]= Kxx*c/Fx ;
      dyCoeff[1]= Kxy*c/Fx ;
      dyCoeff[2]= Kyx*c/Fx ;
      dyCoeff[3]= Kyy*c/Fx ;
      dyCoeff[4]= Cxx*c*u/(r*Fx) ;
      dyCoeff[5]= Cxy*c*u/(r*Fx) ;
      dyCoeff[6]= Cyx*c*u/(r*Fx) ;
      dyCoeff[7]= Cyy*c*u/(r*Fx) ;
      
      for(i=0;i<8;i++)
          informationOfDynamicCoeff_dimensionLess[i]=dyCoeff[i];
      
      delete[] dyCoeff;
}

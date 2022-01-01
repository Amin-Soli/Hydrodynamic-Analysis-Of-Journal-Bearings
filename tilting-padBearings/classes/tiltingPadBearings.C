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
tiltingPadBearings<T>::tiltingPadBearings (const char *path)
{
   readOfBearingProperties(path);  
   
   solution = new T[numberOfPads+2];
   
   informationOfSolution = new T[numberOfPads+2];
   
   informationOfDynamicCoeff_dimensionLess = new double[9];
	
   teta_eachPad_from_Xaxis = new double[numberOfPads];
	
   tetaPivot_eachPad_from_Xaxis = new double[numberOfPads];

   teta_eachPad_from_Xaxis[0]= x_1 ;
    
   tetaPivot_eachPad_from_Xaxis[0]= xp_1 ;
    
   for(int i=1;i<numberOfPads;i++)
   {
	  teta_eachPad_from_Xaxis[i]=teta_eachPad_from_Xaxis[i-1] + teta_pad + (360- numberOfPads*teta_pad)/numberOfPads ; 
	  tetaPivot_eachPad_from_Xaxis[i]=tetaPivot_eachPad_from_Xaxis[i-1] + teta_pad + (360- numberOfPads*teta_pad)/numberOfPads ;
   }
           
   pi=acos(-1);
   
   for(int i=0;i<numberOfPads;i++)
   {
	  teta_eachPad_from_Xaxis[i]=teta_eachPad_from_Xaxis[i]*pi/180.0 ;
	  tetaPivot_eachPad_from_Xaxis[i]=tetaPivot_eachPad_from_Xaxis[i]*pi/180.0 ;
   }
     
   teta_pad = teta_pad*pi/180.0 ;
    
   u = N*2*pi*r/60.0 ;
    
   Omega = u/r ;

   dx=teta_pad/(n-1);
    
   dy=L/(m-1); 
	
   p_in = 0 ; 
   
   frequency = frequency_dimensionless*Omega;
   
   real_eigenValue = real_eigenValue_dimensionless*Omega ;
}

/*****************************************************************************************/
/******************************* member functions ****************************************/
/*****************************************************************************************/

template <class T>
void tiltingPadBearings<T>::solve_withPrintPressure (T *x)
{
     int i ,j ;
	 T force_x=0, force_y=0 ;
     T *functions;
     T **p;
     p = new T *[m];
     for(i = 0; i <m; i++)
     p[i] = new T[n];
     
	 functions=new T[numberOfPads+2];
	 
	 for(i=0;i<numberOfPads;i++)
	 {
		 calc_p(x,p,teta_eachPad_from_Xaxis[i],i);
		 
		 stringstream ss;
		 ss<<(i+1);
		 string str= ss.str();
	     string filename= "pad" + str + "_pressureDistributionOnMidplaneOfBearing.txt" ;
	     ofstream file(filename.c_str());
	     	file<<"theta (degree)"<<'\t'<<'\t'<<"presuure (Pa)"<<endl<<endl;
		 for(j=0;j<n;j++)
	         file<<j*dx*180.0/pi<<'\t'<<'\t'<<'\t'<<p[m/2][j]<<endl;
	     file.close(); 
	     
		 functions[i] = momentom_pad(p);
		 force_x = force_x + Fx_pad(p,teta_eachPad_from_Xaxis[i]);
		 force_y = force_y + Fy_pad(p,teta_eachPad_from_Xaxis[i]);		 
	 }
	 
	 functions[numberOfPads]=force_y ;
	 
	 functions[numberOfPads+1]= force_x - W ;
	 
	 for(i=0;i<(numberOfPads+2);i++)
	     solution[i]=x[i];
	 	 
	 for(i=0;i<(numberOfPads+2);i++)
	     informationOfSolution[i]=functions[i];	 
	
	 for(i = 0; i < m; ++i) 
         delete[] p[i];   
     
     delete[] p;
	 
	 delete[] functions; 
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
T* tiltingPadBearings<T>::solve (T *x)
{
     int i ,j ;
	 T force_x=0, force_y=0 ;
     T *functions;
     T **p;
     p = new T *[m];
     for(i = 0; i <m; i++)
     p[i] = new T[n];
     
	 functions=new T[numberOfPads+2];
	 
	 for(i=0;i<numberOfPads;i++)
	 {
		 calc_p(x,p,teta_eachPad_from_Xaxis[i],i);
		 functions[i] = momentom_pad(p);
		 force_x = force_x + Fx_pad(p,teta_eachPad_from_Xaxis[i]);
		 force_y = force_y + Fy_pad(p,teta_eachPad_from_Xaxis[i]);		 
	 }
	 
	 functions[numberOfPads]=force_y ;
	 
	 functions[numberOfPads+1]= force_x - W ;
	 	 
	 for(i=0;i<(numberOfPads+2);i++)
	     informationOfSolution[i]=functions[i];	 
	
	 for(i = 0; i < m; ++i) 
         delete[] p[i];   
     
     delete[] p;
	 
	 return functions; 
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
void  tiltingPadBearings<T>::dynamicCoeffTotal(double *x)
{
     double **dyCoeff_pad ;
     double *total_dyCoeff ;
     double sum;
     
     dyCoeff_pad = new double *[numberOfPads];
     
     for(int i = 0; i < numberOfPads; i++)
        dyCoeff_pad[i] = new double[9];

     total_dyCoeff=new double[9];
	 	 
	 for(int i=0;i< numberOfPads; i++)
	     dyCoeff_pad[i]=dynamicCoeffEachPad(x,teta_eachPad_from_Xaxis[i],tetaPivot_eachPad_from_Xaxis[i],i);
 
	 /****************** total dyCoeff **************************/
	 
	 for(int i=0;i<9;i++)
	 {
		sum=0.0; 
	    for(int j=0;j<numberOfPads;j++)
	       sum = sum + dyCoeff_pad[j][i];
	    
	    total_dyCoeff[i]=sum;   	     
	 }
	 /****************** non dimension *************************/
	 
	 total_dyCoeff[1] = (total_dyCoeff[1]*c)/total_dyCoeff[0] ;
	 total_dyCoeff[2] = (total_dyCoeff[2]*c)/total_dyCoeff[0] ;
	 total_dyCoeff[3] = (total_dyCoeff[3]*c)/total_dyCoeff[0] ;
	 total_dyCoeff[4] = (total_dyCoeff[4]*c)/total_dyCoeff[0] ;
	 total_dyCoeff[5] = (total_dyCoeff[5]*c*u)/(r*total_dyCoeff[0]) ;
	 total_dyCoeff[6] = (total_dyCoeff[6]*c*u)/(r*total_dyCoeff[0]) ;
	 total_dyCoeff[7] = (total_dyCoeff[7]*c*u)/(r*total_dyCoeff[0]) ;
	 total_dyCoeff[8] = (total_dyCoeff[8]*c*u)/(r*total_dyCoeff[0]) ;	
	 
	 for(int i=0;i<9;i++)
	     informationOfDynamicCoeff_dimensionLess[i]=total_dyCoeff[i];
	  
	 for(int i = 0; i < numberOfPads; ++i) 
         delete[] dyCoeff_pad[i];   
     
     delete[] dyCoeff_pad;
	 
	 delete[] total_dyCoeff;
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
void tiltingPadBearings<T>::calc_p(T *x, T **p, double &pad_angel_from_teta , int &k)
{
   T e_p, cop_p, error_p ;
   int i,j ;
	
   for (j=0;j<m;j++)
     for (i=0;i<n;i++)
         p[j][i]=p_in;

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

             if(i==0 || i==n-1)
                p[j][i]=p_in;


             else
              {

                cop_p=p[j][i];

                p[j][i]=(
                          ( pow(h(i+1,x,pad_angel_from_teta,k),3) - pow(h(i-1,x,pad_angel_from_teta,k),3) )*(p[j][i+1]-p[j][i-1])/(48.0*v*dx*dx)
                          + pow(h(i,x,pad_angel_from_teta,k),3)*(p[j][i+1]+ p[j][i-1])/(12.0*v*dx*dx )
                          + r*r*pow(h(i,x,pad_angel_from_teta,k),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                          - u*r*(h(i+1,x,pad_angel_from_teta,k)-h(i-1,x,pad_angel_from_teta,k))/(4.0*dx)
                        )
                        /
                        (2.0*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dx*dx) + 2.0*r*r*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dy*dy) );

                if(p[j][i]< p_in)
                   p[j][i]= p_in;

                 error_p=fabs(cop_p-p[j][i]);

                 if(error_p>e_p)
                    e_p=error_p;
                   

	            }


            }

        }

	}while( e_p>errorLimit );
	  
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
T tiltingPadBearings<T>::h(int i, T *x, double &pad_angel_from_teta , int &k)
{	
    return c*( 1 - preload*cos(teta_pad/2.0-i*dx) + x[numberOfPads+1]*cos( (i*dx + pad_angel_from_teta) - x[numberOfPads] )

           + x[k]*r*sin(teta_pad/2.0 - i*dx)/c  );	
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
T tiltingPadBearings<T>::momentom_pad(T **p)
{	
     T momentom_=0 ;
     
	 for (int j=0;j<m;j++)
        for(int i=0;i<n;i++)
		   momentom_= momentom_ + r*dx*dy*r*sin(teta_pad/2.0-i*dx)*p[j][i];
		   
     return momentom_;	
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
T tiltingPadBearings<T>::Fx_pad(T **p,double &pad_angel_from_teta)
{	
    T force=0;
     
	for (int j=0;j<m;j++)
        for(int i=0;i<n;i++)
            force=force+(r*dx*dy*p[j][i]*cos(i*dx+pad_angel_from_teta) );

    return force;
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
T tiltingPadBearings<T>::Fy_pad(T **p,double &pad_angel_from_teta)
{
     T force=0;

	 for (int j=0;j<m;j++)
        for(int i=0;i<n;i++)
            force=force+(r*dx*dy*p[j][i]*sin(i*dx+pad_angel_from_teta) );

     return force;
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
void tiltingPadBearings<T>::writeAllInformationOfBearing()
{
   cout<<"\n"<<"bearing of properties : "<<"\n"<<"\n";
   cout << "value of clearance = " << c << "\n";
   cout << "value of length = " << L << "\n";
   cout << "value of radius = " << r << "\n";
   cout << "value of viscosity = " << v << "\n";
   cout << "value of rotationOfSpeed (rpm) = " << N << "\n";	
   cout << "value of preload = " << preload << "\n";	
   cout << "value of magnitudOfAnglePad (degree) = " << teta_pad*180.0/pi << "\n";
   cout << "value of numberOfPads = " << numberOfPads << "\n";
   cout << "value of numberOfNodsInTetaDirection = " << n << "\n";
   cout << "value of numberOfNodsInLengthDirection = " << m << "\n";
   cout << "value of padNumberOneTetaFromXaxis (degree) = " << x_1 << "\n";
   cout << "value of pivotPadNumberOneTetaFromXaxis (degree) = " << xp_1 << "\n";
   cout << "value of maximum error for calculating pressure = " << errorLimit << "\n";
   cout << "value of instabilityDimensionlessFrequency = " << frequency_dimensionless << "\n";
   cout << "value of realTermDimensionlessEigenValueForInstability = " << real_eigenValue_dimensionless << "\n";   
   cout<<"\n"<<"...................................................................."<<"\n";
   cout<<"\n"<<"all information of solution : "<<"\n"<<"\n";
   
   for(int i=0;i<numberOfPads;i++)
       cout << "tilting angle of pad ("<<i+1<<") = " << solution[i] << "\n"; 
       
   cout << "loading angle = " << solution[numberOfPads]*180.0/pi << "\n";    
   cout << "epsilon = " << solution[numberOfPads + 1] << "\n";
   cout << "shaftWeight = " << W << "\n";
   cout << "Fy - shaftWeight = " << informationOfSolution[numberOfPads+1] << "\n";
   cout << "Fx = " << informationOfSolution[numberOfPads] << "\n";
  
   for(int i=0;i<numberOfPads;i++)
       cout << "momentom of pad ("<<i+1<<") = " << informationOfSolution[i] << "\n";  
       
   cout<<"\n"<<"...................................................................."<<"\n";
   cout<<"\n"<<"DimensionLess Of Dynamic Coeffs : "<<"\n"<<"\n";
   cout << "Kxx = " << informationOfDynamicCoeff_dimensionLess[1] << "\n";
   cout << "Kxy = " << informationOfDynamicCoeff_dimensionLess[2] << "\n";
   cout << "Kyx = " << informationOfDynamicCoeff_dimensionLess[3] << "\n";
   cout << "Kyy = " << informationOfDynamicCoeff_dimensionLess[4] << "\n";
   cout << "Cxx = " << informationOfDynamicCoeff_dimensionLess[5] << "\n";
   cout << "Cxy = " << informationOfDynamicCoeff_dimensionLess[6] << "\n";
   cout << "Cyx = " << informationOfDynamicCoeff_dimensionLess[7] << "\n";
   cout << "Cyy = " << informationOfDynamicCoeff_dimensionLess[8] << "\n";
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
void tiltingPadBearings<T>::readOfBearingProperties(const char *path)
{
    ifstream file(path);

	string searchForClearance = "clearance";
	string searchForLength = "length";
	string searchForRadius = "radius";
	string searchForViscosity = "viscosity";
	string searchForRotationOfSpeed = "rotationOfSpeed";
	string searchForShaftWeight = "ShaftWeight";
	string searchForPreload = "preload";
	string searchForMagnitudOfAngelPad = "magnitudOfAnglePad";
	string searchForPadNumberOneTetaFromXaxis = "padNumberOneTetaFromXaxis";
    string searchForPivotPadNumberOneTetaFromXaxis = "pivotPadNumberOneTetaFromXaxis";
	string searchForNumberOfPads = "numberOfPads";
	string searchForNumberOfNodsInTetaDirection = "numberOfNodsInTetaDirection";
	string searchForNumberOfNodsInLengthDirection = "numberOfNodsInLengthDirection";
	string searchForInstabilityDimensionlessFrequency = "instabilityDimensionlessFrequency";
	string searchForRealTermDimensionlessEigenValueForInstability = "realTermDimensionlessEigenValueForInstability";	
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
		
		if (lineOfText.find(searchForRotationOfSpeed, 0) != string::npos)
		  file >> N;
		
		if (lineOfText.find(searchForShaftWeight, 0) != string::npos)
		  file >> W;
		
		if (lineOfText.find(searchForPreload, 0) != string::npos)
		  file >> preload;
		
		if (lineOfText.find(searchForMagnitudOfAngelPad, 0) != string::npos)
		  file >> teta_pad;
		
		if (lineOfText.find(searchForNumberOfPads, 0) != string::npos)
		  file >> numberOfPads;
		
		if (lineOfText.find(searchForNumberOfNodsInTetaDirection, 0) != string::npos)
		  file >> n;		  	  
		
		if (lineOfText.find(searchForNumberOfNodsInLengthDirection, 0) != string::npos)
		  file >> m;
	    
	    if (lineOfText.find(searchForPadNumberOneTetaFromXaxis, 0) != string::npos)
		  file >> x_1;
		
	    if (lineOfText.find(searchForPivotPadNumberOneTetaFromXaxis, 0) != string::npos)
		  file >> xp_1;
		  
	    if (lineOfText.find(searchForInstabilityDimensionlessFrequency , 0) != string::npos)
		  file >> frequency_dimensionless;
		
	    if (lineOfText.find(searchForRealTermDimensionlessEigenValueForInstability, 0) != string::npos)
		  file >> real_eigenValue_dimensionless;
		  
		if (lineOfText.find(searchForErrorLimit, 0) != string::npos)
		  file >> errorLimit;
				  	  							
	}
	
	file.close();
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
int tiltingPadBearings<T>::padNumber()
{
	return numberOfPads;
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
double tiltingPadBearings<T>::mass()
{
	return W;
}

/************************************************************************************************************************************/
/************************************************************************************************************************************/

template <class T>
double*  tiltingPadBearings<T>::dynamicCoeffEachPad(double *x, double &pad_angel_from_teta, double &pivot_angel , int &k )
{
     bool flag = false ;
	 int i,j;
	 double  Kxx=0, Kxy=0, Kyx=0, Kyy=0, Cxx=0, Cxy=0, Cyx=0, Cyy=0 ,Fx=0 , coeff_p =0 , coeff_q=0 ,
	         Kzz=0, Kze=0, Kez=0, Kee=0, Czz=0, Cze=0, Cez=0, Cee=0 , Kzz_new=0 , Czz_new=0 ,
	         e_p,cop_p, error_p ,e_p_x,cop_p_x, error_p_x ,
	         e_p_y,cop_p_y, error_p_y ,e_p_dx,cop_p_dx,
	         error_p_dx ,e_p_dy,cop_p_dy, error_p_dy  ;
	          
     double *dyCoeff;
     dyCoeff=new double[9];
     
     double **p, **p_x, **p_y, **p_dx, **p_dy ;
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

               if(i==0 || i==n-1)
               {
                     p[j][i]=p_in;
                     p_x[j][i]=p_in;
                     p_y[j][i]=p_in;
                     p_dx[j][i]=p_in;
                     p_dy[j][i]=p_in;
               }

               else 
               {

                  cop_p=p[j][i];
                  cop_p_x=p_x[j][i];
                  cop_p_y=p_y[j][i];
                  cop_p_dx=p_dx[j][i];
                  cop_p_dy=p_dy[j][i];

                  p[j][i]=(
                            ( pow(h(i+1,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i-1,x,pad_angel_from_teta,k),3)/(12.0*v ) )*(p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                            + pow(h(i,x,pad_angel_from_teta,k),3)*(p[j][i+1]+ p[j][i-1])/(12.0*v*dx*dx )
                            + r*r*pow(h(i,x,pad_angel_from_teta,k),3)*(p[j+1][i]+ p[j-1][i])/(12.0*v*dy*dy )
                            + (pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v ) )*r*r*(p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            - u*r*(h(i+1,x,pad_angel_from_teta,k)-h(i-1,x,pad_angel_from_teta,k))/(4.0*dx)
                          )
                          /
                          (2.0*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dy*dy) );
                 
                  
                  p_x[j][i]=(
                             ( pow(h(i+1,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i-1,x,pad_angel_from_teta,k),3)/(12.0*v ) )*(p_x[j][i+1]-p_x[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x,pad_angel_from_teta,k),3)*(p_x[j][i+1]+ p_x[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x,pad_angel_from_teta,k),3)*(p_x[j+1][i]+ p_x[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v ) )*r*r*(p_x[j+1][i]-p_x[j-1][i])/(4.0*dy*dy)
                             + u*r*sin(i*dx + pad_angel_from_teta)/2.0 + ( 3*cos((i+1)*dx + pad_angel_from_teta)*pow(h(i+1,x,pad_angel_from_teta,k),2)/(12.0*v ) - 3*cos((i-1)*dx + pad_angel_from_teta)*pow(h(i-1,x,pad_angel_from_teta,k),2)/(12.0*v ) ) * (p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                             + 3*cos(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)*(p[j][i-1]-2*p[j][i]+p[j][i+1])/(12.0*v*dx*dx )
                             + 3*r*r*cos(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             +(3*cos(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)/(12.0*v ) - 3*cos(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dy*dy) );

                  
                  p_y[j][i]=(
                             ( pow(h(i+1,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i-1,x,pad_angel_from_teta,k),3)/(12.0*v ) )*(p_y[j][i+1]-p_y[j][i-1])/(4.0*dx*dx)
                             + pow(h(i,x,pad_angel_from_teta,k),3)*(p_y[j][i+1]+ p_y[j][i-1])/(12.0*v*dx*dx )
                             + r*r*pow(h(i,x,pad_angel_from_teta,k),3)*(p_y[j+1][i]+ p_y[j-1][i])/(12.0*v*dy*dy )
                             + (pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v ) )*r*r*(p_y[j+1][i]-p_y[j-1][i])/(4.0*dy*dy)
                             - u*r*cos(i*dx + pad_angel_from_teta)/2.0 + ( 3*sin((i+1)*dx + pad_angel_from_teta)*pow(h(i+1,x,pad_angel_from_teta,k),2)/(12.0*v ) - 3*sin((i-1)*dx + pad_angel_from_teta)*pow(h(i-1,x,pad_angel_from_teta,k),2)/(12.0*v ) ) * (p[j][i+1]-p[j][i-1])/(4.0*dx*dx)
                             + 3*sin(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)*(p[j][i-1]-2*p[j][i]+p[j][i+1])/(12.0*v*dx*dx )
                             + 3*r*r*sin(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)*(p[j-1][i]-2*p[j][i]+p[j+1][i])/(12.0*v*dy*dy )
                             + (3*sin(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)/(12.0*v ) - 3*sin(i*dx + pad_angel_from_teta)*pow(h(i,x,pad_angel_from_teta,k),2)/(12.0*v ) ) *r*r* (p[j+1][i]-p[j-1][i])/(4.0*dy*dy)
                            )
                            /
                            (2.0*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dy*dy) );
                  
                  
                  p_dx[j][i]=(
                               ( pow(h(i+1,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i-1,x,pad_angel_from_teta,k),3)/(12.0*v ) )*(p_dx[j][i+1]-p_dx[j][i-1])/(4.0*dx*dx)
                               + pow(h(i,x,pad_angel_from_teta,k),3)*(p_dx[j][i+1]+ p_dx[j][i-1])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x,pad_angel_from_teta,k),3)*(p_dx[j+1][i]+ p_dx[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v ) )*r*r*(p_dx[j+1][i]-p_dx[j-1][i])/(4.0*dy*dy)
                               - r*r*cos(i*dx + pad_angel_from_teta)
                              )
                              /
                              (2.0*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dy*dy) );
                              
                  
                  p_dy[j][i]=(
                               ( pow(h(i+1,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i-1,x,pad_angel_from_teta,k),3)/(12.0*v ) )*(p_dy[j][i+1]-p_dy[j][i-1])/(4.0*dx*dx)
                               + pow(h(i,x,pad_angel_from_teta,k),3)*(p_dy[j][i+1]+ p_dy[j][i-1])/(12.0*v*dx*dx )
                               + r*r*pow(h(i,x,pad_angel_from_teta,k),3)*(p_dy[j+1][i]+ p_dy[j-1][i])/(12.0*v*dy*dy )
                               + (pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v )- pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v ) )*r*r*(p_dy[j+1][i]-p_dy[j-1][i])/(4.0*dy*dy)
                               - r*r*sin(i*dx + pad_angel_from_teta)
                              )
                              /
                              (2.0*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dx*dx)+ 2.0*r*r*pow(h(i,x,pad_angel_from_teta,k),3)/(12.0*v*dy*dy) );
                              
                                          
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
	        e_p>errorLimit || e_p_x>errorLimit  || e_p_y>errorLimit 
	                      || e_p_dx>errorLimit  || e_p_dy>errorLimit 	          
	      );
	      
	for (j=0;(j<m && flag == false);j++)
       for(i=0;i<n;i++)
          if(p[j][i] != 0)
          {
             flag = true;
             break;
		  } 
            
    if(flag == false)
    {
		for(i=0;i<9;i++)
		    dyCoeff[i]=0;
		    
		return dyCoeff;
	}
		          

    for (j=0;j<m;j++)
      for(i=0;i<n;i++)
      {
		  Fx=Fx + (r*dx*dy*p[j][i]*cos(i*dx + pad_angel_from_teta) );
          Kxx=Kxx + (r*dx*dy*p_x[j][i]*cos(i*dx + pad_angel_from_teta) );
          Kyx=Kyx + (r*dx*dy*p_x[j][i]*sin(i*dx + pad_angel_from_teta) );
          Kxy=Kxy + (r*dx*dy*p_y[j][i]*cos(i*dx + pad_angel_from_teta) );
          Kyy=Kyy + (r*dx*dy*p_y[j][i]*sin(i*dx + pad_angel_from_teta) );
          Cxx=Cxx + (r*dx*dy*p_dx[j][i]*cos(i*dx + pad_angel_from_teta) );
          Cxy=Cxy + (r*dx*dy*p_dx[j][i]*sin(i*dx + pad_angel_from_teta) );
          Cyx=Cyx + (r*dx*dy*p_dy[j][i]*cos(i*dx + pad_angel_from_teta) );
          Cyy=Cyy + (r*dx*dy*p_dy[j][i]*sin(i*dx + pad_angel_from_teta) );
      } 
      
    Kzz = Kxx*cos(pivot_angel)*cos(pivot_angel) + Kyy*sin(pivot_angel)*sin(pivot_angel) + (Kxy + Kyx)*cos(pivot_angel)*sin(pivot_angel);
    Kze = Kxy*cos(pivot_angel)*cos(pivot_angel) - Kyx*sin(pivot_angel)*sin(pivot_angel) + (Kyy - Kxx)*cos(pivot_angel)*sin(pivot_angel);
    Kez = Kyx*cos(pivot_angel)*cos(pivot_angel) - Kxy*sin(pivot_angel)*sin(pivot_angel) + (Kyy - Kxx)*cos(pivot_angel)*sin(pivot_angel);	
    Kee = Kyy*cos(pivot_angel)*cos(pivot_angel) + Kxx*sin(pivot_angel)*sin(pivot_angel) - (Kxy + Kyx)*cos(pivot_angel)*sin(pivot_angel);
    Czz = Cxx*cos(pivot_angel)*cos(pivot_angel) + Cyy*sin(pivot_angel)*sin(pivot_angel) + (Cxy + Cyx)*cos(pivot_angel)*sin(pivot_angel);    
    Cze = Cxy*cos(pivot_angel)*cos(pivot_angel) - Cyx*sin(pivot_angel)*sin(pivot_angel) + (Cyy - Cxx)*cos(pivot_angel)*sin(pivot_angel);
    Cez = Cyx*cos(pivot_angel)*cos(pivot_angel) - Cxy*sin(pivot_angel)*sin(pivot_angel) + (Cyy - Cxx)*cos(pivot_angel)*sin(pivot_angel);        	
    Cee = Cyy*cos(pivot_angel)*cos(pivot_angel) + Cxx*sin(pivot_angel)*sin(pivot_angel) - (Cxy + Cyx)*cos(pivot_angel)*sin(pivot_angel);
    
    coeff_p = (Kee + real_eigenValue*Cee)/( (Kee + real_eigenValue*Cee)*(Kee + real_eigenValue*Cee) + frequency*frequency*Cee*Cee ) ;
    coeff_q = (frequency*Cee)/( (Kee + real_eigenValue*Cee)*(Kee + real_eigenValue*Cee) + frequency*frequency*Cee*Cee ) ;
   
    Kzz_new = Kzz + Cez*(  (Kze + real_eigenValue*Cze)*(real_eigenValue*coeff_p - 
                           frequency*coeff_q) + frequency*Cze*(real_eigenValue*coeff_q + frequency*coeff_p)  )
                  - (Kez + real_eigenValue*Cez)*(  (Kze + real_eigenValue*Cze)*(coeff_p + real_eigenValue*coeff_q/frequency)
                                                  + frequency*Cze*(coeff_q - real_eigenValue*coeff_p/frequency)  ) ;
       
    Czz_new = Czz - Cez*(  (Kze + real_eigenValue*Cze)*coeff_p + frequency*Cze*coeff_q  )
                  + (Kez + real_eigenValue*Cez)/frequency *(  (Kze + real_eigenValue*Cze)*coeff_q
                                                              - frequency*Cze*coeff_p  )  ;
     
    dyCoeff[0] = Fx ;
    dyCoeff[1] = Kzz_new*cos(pivot_angel)*cos(pivot_angel) ;
    dyCoeff[2] = Kzz_new*cos(pivot_angel)*sin(pivot_angel) ;
    dyCoeff[3] = Kzz_new*cos(pivot_angel)*sin(pivot_angel) ;
    dyCoeff[4] = Kzz_new*sin(pivot_angel)*sin(pivot_angel) ;
    dyCoeff[5] = Czz_new*cos(pivot_angel)*cos(pivot_angel) ;
    dyCoeff[6] = Czz_new*cos(pivot_angel)*sin(pivot_angel) ;
    dyCoeff[7] = Czz_new*cos(pivot_angel)*sin(pivot_angel) ;
    dyCoeff[8] = Czz_new*sin(pivot_angel)*sin(pivot_angel) ;
    
    for(i = 0; i < m; ++i) 
        delete[] p[i];   
     
    delete[] p;
     
    return dyCoeff;
}

template <class T>
class fixedPadBearings {
	
  private:
  
    double  c , L , r , v , N , sommerfeld , preload ,          
             
            u  , W , dx , dy , pi , p_in , errorLimit;
            
                  
    double *magnitudGrooves ; 
    
    double *centerGrooves ;
    
    T *informationOfSolution ;
    
    T phi_U , phi_L, epsilon_L, epsilon_U ;
    
    double *informationOfDynamicCoeff_dimensionLess ;
    
    int n , m , numberOfGrooves  ;
    
    string typeOfBearing , isGrooves ;
    
  public:
    
    fixedPadBearings (const char *path );

    T* solve (T *x);
    
    void solve_withPrintPressure (T *x);
    
    void dynamicCoeffTotal(double *x);
    
    double mass();
    
    T h(int i, T *x);
    
    void readOfBearingProperties (const char *path);
    
    void writeAllInformationOfBearing();
};

#include "fixedPadBearings.C"

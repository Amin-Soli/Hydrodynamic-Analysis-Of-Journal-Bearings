template <class T>
class tiltingPadBearings {
	
  private:
  
    double  c , L , r , v , N , preload , teta_pad  , x_1 , xp_1 , 
    
            frequency_dimensionless , real_eigenValue_dimensionless , frequency ,
             
            real_eigenValue , u  , W , Omega , dx , dy , pi , p_in , errorLimit ;
                  
    double *teta_eachPad_from_Xaxis ; 
    
    double *tetaPivot_eachPad_from_Xaxis ;
    
    T *solution, *informationOfSolution ;
    
    double *informationOfDynamicCoeff_dimensionLess ;
    
    int n , m , numberOfPads ;
    
  public:
    
    tiltingPadBearings (const char *path);

    T* solve (T *x);
    
    void solve_withPrintPressure (T *x);
    
    void dynamicCoeffTotal(double *x);
    
    double* dynamicCoeffEachPad(double *x, double &pad_angel_from_teta, double &pivot_angel , int &k );
    
    int padNumber();
    
    double mass();
    
    void calc_p(T *x, T **p, double &pad_angel_from_teta , int &k);
    
    T h(int i, T *x, double &pad_angel_from_teta , int &k);
    
    T momentom_pad(T **p);
    
    T Fx_pad(T **p,double &pad_angel_from_teta);
    
    T Fy_pad(T **p, double &pad_angel_from_teta);
    
    void readOfBearingProperties (const char *path);
    
    void writeAllInformationOfBearing();
};

#include "tiltingPadBearings.C"

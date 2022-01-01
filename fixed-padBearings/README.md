
# This program has been written in the form of Object-Oriented Programming in C++ (Class).

# Note: This program is executable only in linux.

# This program calculates the static solution in two dimensions for two types of fixed-pad journal bearings called Elliptical and Plain journal bearings using the NOX optimization package in the Trilinos library. The considered algorithm uses Newton optimization with the Trust region method.

# The algorithm is:
  1. Make initial guess for unknown variables (X vector)  
     X[0] = eccentricity ratio  
     X[1] = loading angle
  2. Solve the Reynolds equations using finite difference method and calculate pressure distributions on the bearing
  3. Obtain forces exerted on the shaft
  4. Calculate the Residuals vector which is defined as:  
     Residuals[0] = (summation of forces applied on the shaft in the y direction) - shaft weight  
     Residuals[1] = summation of forces applied on the shaft in the x direction
  5. Calculate stepSize (alpha) using NOX optimization algorithm with Trust Region method
  6. Obtain new X vector as:  
     X_new = X_old + alpha
  7. Return to step 2
  * This algorithm is repeated until the norm of Residuals is smaller than 1.0e-6

# After calculating eccentricity ratio and loading angle, this program solves dynamic equations in order to calculate dynamic coefficients (damping and stiffness coefficients).  

# This program consists of three folders: 
  * One folder which is named "classes", includes all required classes to solve the problem
  * One folder which is named "solver", includes your main code file ("fixedPadBearings.cpp") and executable file "trilinosExample" 
  * Another folder which is named "inputData", includes "bearingProperties.txt" file to enter required data for solving the problem

# Note: "fixedPadBearings.cpp" file in "solver" folder is already compiled, and the executable file "trilinosExample" exists there.

# Note: If you want to compile "fixedPadBearings.cpp" file again, you should first install Trilinos library (enable NOX package and its dependencies) and then compile it using cmake.

# If you want to run this code, all you need to do, is:
  1. Go to "inputData" folder 
  2. Open "bearingProperties.txt" file and enter inlet parameters
  3. Open terminal
  4. go to directory of "solver" folder
  5. type "./trilinosExample"

class PIDController {
  
  // setpoint and current values;
  float setpoint;

  
  // Proportional, Integral and Differential controller constants
  // Set to zero to turn off that controller.
  float kP, kI, kD;
  
  //  Integral and error terms;
  private float integral;
  private float prev_error;
  
  // Time constants for Integration and Differentiation (seperate time values)
  // set iT = 0 to turn off.  Never set dT=0 (divide by zero for Differentiation). 
  float iT;
  float dT;



  // PID Constructor
  PIDController(){
    // setpoint and current values;
    setpoint=0;
    // Proportional, Integral and Differential controller constants
    kP=kI=kD=0;
    //  Integral and error terms;
    integral=0;
    prev_error=0;;
    // Time constants for Integration and Differentiation (seperate time values)
    iT=0.01;
    dT=0.01;
  }



  float PIDOutputValue(float measured_value){
    float error;
    float derivative;
    float output;



    error=(setpoint-measured_value)/(setpoint+Float.MIN_VALUE);
    
    integral=integral+error*(iT+Float.MIN_VALUE);
    if (integral > 1.0)
      integral=1.0;
    else if (integral < 0.0)
      integral=0.0;
    
    derivative=(error - prev_error)/(dT+Float.MIN_VALUE);
    
    output= kP*error + kI*integral + kD*derivative;
    
    prev_error=error;
    
    if (output > 1.0)
      output=1.0;
    else if (output < 0.0)
      output=0.0;
      
    return(output);
    
  }
  
  
  


  
}

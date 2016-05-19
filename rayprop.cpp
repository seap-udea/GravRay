#include <gravray.cpp>
using namespace std;

int main(int argc,char* argv[])
{
  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();

  ////////////////////////////////////////////////////
  //GET ARGUMENTS
  ////////////////////////////////////////////////////
  /*
    Funcion: 

      Propagate the orbit of an object in time starting at a given
      position and time and including in the integration the combined
      gravitational effect.

    Arguments are: 

      object, latitude (degrees), longitude (degrees), elevation
      (meters), date

    Date format:

       MM/DD/CCYY HH:MM:SS.dcm UTC-L

    Example:
    
       ./rayprop.exe  "07/19/2015 00:00:00.000 UTC-5"
  */
}

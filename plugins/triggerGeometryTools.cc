#include "L1Trigger/phase2L1EcalTimingTP/interface/triggerGeometryTools.hh"

triggerGeometryTools::triggerGeometryTools(){};



float triggerGeometryTools::getRecoPhi(int iphi){
    float towerPhiMap[72]=                        
      {0.044, 0.131, 0.218, 0.305, 0.393, 0.480, 0.567, 0.654, 0.742, 0.829, 0.916, 1.004, 1.091, 1.178, 1.265, 1.353, 1.440, 1.527, 1.614, 1.702, 1.789, 1.876, 1.963, 2.051, 2.138, 2.225, 2.313, 2.400, 2.487, 2.574, 2.662, 2.749, 2.836, 2.923, 3.011, 3.098,
       -3.098, -3.011, -2.923, -2.836, -2.749, -2.662, -2.574, -2.487, -2.400, -2.313, -2.225, -2.138, -2.051, -1.963, -1.876, -1.789, -1.702, -1.614, -1.527, -1.440, -1.353, -1.265, -1.178, -1.091, -1.004, -0.916, -0.829, -0.742, -0.654, -0.567, -0.480, -0.393, -0.305, -0.218, -0.131, -0.044};
    /*{-0.1305, -0.0435,  0.0435, 0.1308, 0.2178, 0.3048, 
      0.3918,  0.4788, 0.5658, 0.6528, 0.7398, 0.8268, 
      0.9138, 1.0008, 1.0878, 1.1748, 1.2618, 1.3488, 
      1.4358, 1.5228, 1.6098, 1.6968, 1.7838, 1.8708, 
      1.9578, 2.0448, 2.1318, 2.2188, 2.3058, 2.3928, 
      2.4798, 2.5668, 2.6538, 2.7408, 2.8278, 2.9148, 
      3.0018, 3.0888, -3.0885, -3.0015, -2.9145, -2.8275, 
      -2.7405, -2.6535, -2.5665, -2.4795, -2.3925, -2.3055, 
      -2.2185, -2.1315, -2.0445, -1.9575, -1.8705, -1.7835, 
      -1.6965, -1.6095, -1.5225, -1.4355, -1.3485, -1.2615, 
      -1.1745, -1.0875, -1.0005, -0.9135, -0.8265, -0.7395, 
      -0.6525, -0.5655, -0.4785, -0.3915, -0.3045, -0.2175 };*/
    float temp = towerPhiMap[iphi-1];
    
    return temp;
};

float triggerGeometryTools::getRecoEta(int ieta, short zside){
    //from tower eta -28 to 28
  float towerEtaMap[57]=   { 
      -2.913, //-2.739, 
      -2.565, -2.391, 2.217, 	//switch to smaller trigger towers here    
      -2.0445, -1.9575, -1.8705, -1.7835, -1.6965, 
      -1.6095, -1.5225, -1.4355, -1.3485, -1.2615, 
    -1.1745, -1.0875, -1.0005, -0.9135, -0.8265, 
      -0.7395, -0.6525, -0.5655, -0.4785, -0.3915, 
      -0.3045, -0.2175, -0.1305, -0.0435, 0.0435, 
      0.1305, 0.2175, 0.3045, 0.3915, 0.4785, 
      0.5655, 0.6525, 0.7395, 0.8265, 0.9135, 
      1.0005, 1.0875, 1.1745, 1.2615, 1.3485, 
      1.4355, 1.5225, 1.6095, 1.6965, 1.7835, 
      1.8705, 1.9575, 2.0445, 2.217, 2.391, 
      2.565, //2.739,
      2.913 
    };
    
    
    
    float eta = -999;
    if(ieta<0 || ieta>(28*2)){
      std::cout<<"Error!!! towereta out of bounds in triggerGeometryTools.h "<<std::endl;
      std::cout<<"ieta "<<ieta<<std::endl;
      exit(0);
    }
    
    if(zside == 1)
      eta = towerEtaMap[ieta];
    else if(zside == -1)
      eta = towerEtaMap[ieta];
    else{
      std::cout<<"Error!!! zside out of bounds in triggerGeometryTools.h "<<std::endl;
      std::cout<<"zside "<<zside<<std::endl;
      exit(0);
    }
    return eta;

};



int triggerGeometryTools::TPGEtaRange(int ieta)
{
  int iEta = 0;
  // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
  // And we have mapped [-28, -1], [1, 28] onto [0, 55]   
  if(ieta < 0)
    iEta = ieta + 28;
  else if(ieta > 0)
    iEta = ieta + 27;
  if(ieta==0){
    std::cout<<"Error! ieta is 0, ieta: "<<ieta<<" iEta "<<iEta<<std::endl;
    exit(0);
  }
  return iEta;
}

  int convertGenEta(double inputEta) 
  {
    const double tpgEtaValues[27] = {
      0.087,      
      0.174, // HB and inner HE bins are 0.348 wide
      0.261,
      0.348,
      0.522,
      0.609,
      0.696,
      0.783,
      0.870,
      0.957,
      1.044,
      1.131,
      1.218,
      1.305,
      1.392,
      1.479,
      1.566,
      1.653,
      1.74,
      1.848,
      1.956, // Last two HE bins are 0.432 and 0.828 wide
      2.064,
      2.172,
      2.379,
      2.586,
      2.793,
      3
      //IGNORING HF
      //3.250, // HF bins are 0.5 wide
      //3.750,
      //4.250,
      //4.750
    };

    for (int n=1; n<29; n++){
      //std::cout<<"inputEta "<<inputEta<< " n "<< n <<" tpgEtaValues[n-1] "<< tpgEtaValues[n-1] << " abs(inputEta)<tpgEtaValues[n-1]"<<std::endl;
      if (std::fabs(inputEta)<tpgEtaValues[n-1]) {
	//Positive eta is >28 negative eta is 0 to 27
	if(inputEta>0){ return n + 28;}
	else{ return n;}
	break;
      }
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputeta: "<<inputEta<<std::endl;
    return -9;
  }

  //-pi < phi <= +pi,
  int convertGenPhi(double inputPhi)
  {
    double posPhi[36];
    for(int n = 0; n < 36; n++)
      posPhi[n] = (0.087) * n + 0.0435;
    double negPhi[36];
    for(int n = 0; n < 36; n++)
      negPhi[n] = -3.14159 + 0.087 * n - 0.0435;

    //1 to 36 is 0 to pi
    if( 3.1416 > inputPhi && inputPhi >= 0){

      for(int n = 1; n < 36; n++){
	//std::cout<<"inputPhi "<<inputPhi<< " posPhi[n-1] "<< posPhi[n-1] << " n "<<n<<std::endl;
	if(inputPhi <= posPhi[n-1]){
	  int tpgPhi = n;
	  return tpgPhi;
	}
      }
    }

    //37 to 72 is -pi to 0
    else if(-3.1416 < inputPhi && inputPhi < 0){
      for(int n = 1; n < 36; n++)
	if(inputPhi < negPhi[n-1]){
	  int tpgPhi = n + 36;
	  return tpgPhi;
	}
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputphi: "<<inputPhi<<std::endl;
    return -9;
  }

//DEFINE_FWK_MODULE(triggerGeometryTools);

#include <iostream>
#include <math.h>
#include <random>
#include <cstdlib>
#include <fstream>
#include <chrono> 
#include <ctime>
using namespace std;

int N_EVENTS = 10000;
int N_CONVERGED;
int N_CRYSTALS = 674;
double step_size = 0.0035;//alpha
double error = 3;
double maxiterations = 100;
double MaxFunction = 100;
double dcmin = 0.1;

std::vector<double> CalibrationConstants;


struct CrystalList{
	std::vector<double> crystal_energy;
	std::vector<int> crystal_number;
   	CrystalList();
     	CrystalList(std::vector<double> _crystal_energy, std::vector<int> _crystal_number) : crystal_energy(_crystal_energy), crystal_number(_crystal_number) {};
};

struct Event{
	double track_energy;
	int cluster_size;
	CrystalList crystal_list;
	Event();
	Event(double _energy, double _size) : track_energy(_energy), cluster_size(_size){};
	Event(double _energy, double _size, CrystalList _list ) : track_energy(_energy), cluster_size(_size), crystal_list(_list){};
	
};


double F(Event event, int event_number, int iteration, std::vector<double> constants){
  
	int sum = 0;
        for(unsigned int c=0;c < event.cluster_size;c++){
		int cry_i = event.crystal_list.crystal_number[c];
                sum+=constants[cry_i]*event.crystal_list.crystal_energy[c];
	}
    	F+= pow((sum - event.track_energy)*(1/error) ,2);
    
	return F;
}


std::vector<double> IncrementalGradientDescent(Event event, int j, std::vector<double> constants, std::vector<double> seed_constants, double FVAL){
    double Loss = FVAL;
    double old_c ;
    double new_c ;
    double dc;
    double dFdVm;
    double dVmdc;
    bool converged =false;
    int k=0;
	
    double InitLoss =F(event,j,0, constants);
    std::vector<double> previous_constants = constants;
    
    while(converged == false and k < 100){
       
        for(unsigned int m=0; m<event.cluster_size;m++){

                    int Cm = event.crystal_list.crystal_number[m];//Get crystal number
                    double prediction =0;
                    for(unsigned int i=0;i<event.cluster_size;i++){
                        int Ci = event.crystal_list.crystal_number[i];//Get crystal number
                        prediction +=constants[Ci]*event.crystal_list.crystal_energy[i];  //Find the coveriance matrix A terms  
                     }
                    dFdVm = 2*(prediction -event.track_energy);
                    dVmdc = event.crystal_list.crystal_energy[m];
              
                    new_c = (constants[Cm] - step_size*(1/event.cluster_size)*(1/pow(error,2))*constants[Cm]*dFdVm*dVmdc);
                    dc = abs(new_c - constants[Cm]); 
                    if(dc < 0.1 and (new_c - seed_constants[m])<0.1){
                        constants[Cm] = new_c; 
                    }
             }
            Loss = F(event, j, k, constants);
            if((Loss <= InitLoss) and (Loss < MaxFunction)){
                converged =true;
                N_CONVERGED +=1;
            }

         k++;
       
    	}
    
    	if(converged == true){
          for(unsigned int m=0; m<event.cluster_size;m++){
                int Ci = event.crystal_list.crystal_number[m];
                if (m==event.cluster_size-1){Loss = F(event,j, k,constants);}
          }
        
        return constants;
    	}else {
        return previous_constants;
    	}
}
int main(){
    
     	ofstream outputfile;
     	outputfile.open("Info.csv");
     	std::vector<Event> epoch;
     	double offset; //used to set constant offest
	std::vector<double> RawCalibrationResults; //from other calib sources
	std::vector<double> offset_vector;
     
     	for(int c=0;c<N_CRYSTALS;c++){
		CalibrationConstants.push_back(0);
     	}
    
	 std::random_device rd;
     	 std::mt19937 mt(rd());
	 std::normal_distribution<double> te(46.,3.); 
	 std::normal_distribution<double> cs(4.,1);
	 std::uniform_real_distribution<double> cn(0, N_CRYSTALS);
	 std::uniform_real_distribution<double> randoff(0.1, 1.0);
    
	 for(int c=0;c<N_CRYSTALS;c++){
		
        	auto const off = randoff(mt);
        	offset_vector.push_back(off);
		RawCalibrationResults.push_back(off*0.9);//imagine that the raw constants are only 90% accurate (assume same for all crystals)
	 }
    
         auto start = chrono::high_resolution_clock::now();
	 for(size_t n=0;n<N_EVENTS;n++){
         
		if(n==0){
            		CalibrationConstants = RawCalibrationResults; //seed with raw inputs
       		}
         
		auto const track_energy = te(mt);
		auto const size = cs(mt);
        	int cluster_size = round(size);
		
       		std::vector<double> crystal_energy;
		std::vector<int> crystal_number;
		
		for(int m=0;m<cluster_size; m++){
            		int C_number = round(cn(mt));
			crystal_number.push_back(C_number);
            		offset = offset_vector[C_number];
            		crystal_energy.push_back((1/offset)*track_energy/(cluster_size));
			
		}
       
		CrystalList crystal_list(crystal_energy,crystal_number);
		Event event(track_energy,cluster_size,crystal_list);
        
		CalibrationConstants= IncrementalGradientDescent(event, n, CalibrationConstants, RawCalibrationResults, FVAL);
		auto currenttime = chrono::high_resolution_clock::now();
		auto eventduration = chrono::duration_cast<chrono::microseconds>(currenttime- start); 
		outputfile<<n<<","<<eventduration.count()<<endl;
		epoch.push_back(event);
     	}
    
	for(int i =0 ;i<N_CRYSTALS;i++){
			std::cout<<"constant for crystal "<<i<<" is "<<CalibrationConstants[i]<<"True Offset is "<<offset_vector[i]<<" Residuals "<<CalibrationConstants[i]-offset_vector[i]<<std::endl;
          
	}
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start); 
    cout<<"NEvents Processed "<<N_EVENTS<<" NEVents converged "<<N_CONVERGED<<"Time "<<duration.count()<<endl;
    return 0;
}

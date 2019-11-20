#include <iostream>
#include <math.h>
#include <random>
#include <cstdlib>
#include <fstream>
#include <chrono> 
#include <ctime>
using namespace std;

double FVAL;
int N_EVENTS = 10000;
int N_CONVERGED;
int N_CRYSTALS = 674;
int N_TOTALHITS = 0;
double Loss=0;
double step_size = 0.0035;
double error =3;
double maxiterations = 100;
double MaxFunction = 100;
    
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


double FullF(std::vector<Event> event_list, std::vector<double> constants, double F){
	int sum = 0;
    for(unsigned int e = 0; e< N_EVENTS;e++){
        Event event = event_list[e];
        for(unsigned int c=0;c < event.cluster_size;c++){
            int cry_i = event.crystal_list.crystal_number[c];
            sum+=constants[cry_i]*event.crystal_list.crystal_energy[c];
        }
    F+= pow((sum - event.track_energy)*(1/error) ,2);
    }  
    
	return F;
}


double F(Event event, int event_number, int iteration, std::vector<double> constants, double F, bool is_end){
  
	int sum = 0;
    for(unsigned int c=0;c < event.cluster_size;c++){
		int cry_i = event.crystal_list.crystal_number[c];
        sum+=constants[cry_i]*event.crystal_list.crystal_energy[c];
	}
    F+= pow((sum - event.track_energy)*(1/error) ,2);
    if(is_end ){
        
        cout<<event_number<<","<<iteration<<","<<F<<endl;
        
    }
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
    
    double InitLoss =F(event,j,0, constants,FVAL,false);
    std::vector<double> previous_constants = constants;
    int k=0;
    double Etrk = event.track_energy; //Get the trackers output   
    while(converged == false and k < 100){
       
        for(unsigned int m=0; m<event.cluster_size;m++){

                    int Cm = event.crystal_list.crystal_number[m];//Get crystal number
                    old_c = constants[Cm]; //Find previous constant
                    double Vm=event.crystal_list.crystal_energy[m]; //Get crystal m's energy in this event and this hit
                    double prediction =0;
                    for(unsigned int i=0;i<event.cluster_size;i++){
                        int Ci = event.crystal_list.crystal_number[i];//Get crystal number
                        prediction +=constants[Ci]*event.crystal_list.crystal_energy[i];  //Find the coveriance matrix A terms  
                     }
                    dFdVm = 2*(prediction -Etrk);
                    dVmdc = Vm;
              
                    new_c = (old_c - step_size*(1/pow(error,2))*constants[Cm]*dFdVm*dVmdc);
                    dc = abs(new_c - old_c); 
                    if(dc < 0.1 and (new_c - seed_constants[m])<0.1){
                        constants[Cm] = new_c; 
                    }
             }
            Loss = F(event, j, k, constants, Loss, false);
            if((Loss <= InitLoss) and (Loss < MaxFunction)){
                converged =true;
                N_CONVERGED +=1;
            }

         k++;
       
    }
    
    if(converged ==true){
          for(unsigned int m=0; m<event.cluster_size;m++){
                int Ci = event.crystal_list.crystal_number[m];//Get crystal number
                //cout<<"updated "<<previous_constants[Ci]<<" "<<constants[Ci]<<" with k iterations "<<k<<endl;
                if (m==event.cluster_size-1){Loss = F(event,j, k,constants, Loss,true);}
          }
        
        return constants;
    }else {
        return previous_constants;
    }
    
}

int main(){
    
     ofstream outputfile;
     outputfile.open("Info.csv");
     
     double offset = 0.5;
     std::vector<double> RawCalibrationResults;
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
     std::vector<Event> event_list;
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
        event_list.push_back(event);
     }
    
	for(int i =0 ;i<N_CRYSTALS;i++){
			std::cout<<"constant for crystal "<<i<<" is "<<CalibrationConstants[i]<<"True Offset is "<<offset_vector[i]<<" Residuals "<<CalibrationConstants[i]-offset_vector[i]<<std::endl;
            //outputfile<<i<<","<<CalibrationConstants[i]<<","<<offset_vector[i]<<","<<CalibrationConstants[i]-offset_vector[i]<<std::endl;
            
		}
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start); 
    double endLoss = FullF(event_list, CalibrationConstants, FVAL);
    cout<<"NEvents Processed "<<N_EVENTS<<" NEVents converged "<<N_CONVERGED<<"Time "<<duration.count()<<" Final Loss function "<<endLoss<<endl;
    
	return 0;
}

  

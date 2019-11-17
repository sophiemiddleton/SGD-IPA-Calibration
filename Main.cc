#include <iostream>
#include <math.h>
#include <random>
#include <cstdlib>

using namespace std;

double FVAL =0;
int N_EVENTS = 10000;
int N_CONVERGED;
int N_CRYSTALS = 674;
int N_TOTALHITS = 0;
double Loss=0;
double LearningRate = 0.001;
double error =3;
double d = 0.1;//decay
double maxiterations = 1000;

    
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


double F(Event event, std::vector<double> constants, double F){
	int sum = 0;
    for(unsigned int c=0;c < event.cluster_size;c++){
		int cry_i = event.crystal_list.crystal_number[c];
        sum+=constants[cry_i]*event.crystal_list.crystal_energy[c];
	}
    F+= pow((sum - event.track_energy)*(1/pow(error,2)) ,2);
    cout<<"Test f "<<F<<endl;
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
    cout<<"init F "<<FVAL<<endl;
    double InitLoss =F(event, constants,FVAL);
    std::vector<double> previous_constants = constants;
    //constraints ??? TODO!!!
    int k=0;
    double Etrk = event.track_energy; //Get the trackers output   
    while(converged == false and k < 100){
        cout<<"-------"<<k<<"---------"<<endl;
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
              
                    new_c = (old_c - LearningRate*(1/pow(error,2))*dFdVm*dVmdc);
                    cout<<"old C = "<<old_c<<" new C =  "<<new_c<<endl;
                    dc = abs(new_c - old_c); 
                    if(dc < 0.1 and (new_c - seed_constants[m])<0.1){
                        constants[Cm] = new_c; 
                    }
             }
            Loss = F(event, constants, Loss);

            if((Loss <= InitLoss and Loss < 100) ){
                converged =true;
                N_CONVERGED +=1;
            }

        k++;
       
    }
    cout<<"end F "<<endl;
    FVAL = F(event, constants, FVAL);
    if(converged ==true){
          for(unsigned int m=0; m<event.cluster_size;m++){
                int Ci = event.crystal_list.crystal_number[m];//Get crystal number
                cout<<"updated "<<previous_constants[Ci]<<" "<<constants[Ci]<<" with k iterations "<<k<<endl;
          }
        return constants;
    }else {
        return previous_constants;
    }
    
}

int main(){
    
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
		RawCalibrationResults.push_back(off*0.9);
	 }
	 for(size_t n=0;n<N_EVENTS;n++){
		if(n==0){
            CalibrationConstants = RawCalibrationResults;
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
        cout<<"FVAL before "<<FVAL<<endl;
        CalibrationConstants= IncrementalGradientDescent(event, n, CalibrationConstants, RawCalibrationResults, FVAL);
        cout<<"FVAL after "<<FVAL<<endl;
      
     }
    
	for(int i =0 ;i<N_CRYSTALS;i++){
			std::cout<<"constant for crystal "<<i<<" is "<<CalibrationConstants[i]<<"True Offset is "<<offset_vector[i]<<" Residuals "<<CalibrationConstants[i]-offset_vector[i]<<std::endl;
		}
    cout<<"NEvents Processed "<<N_EVENTS<<" NEVents converged "<<N_CONVERGED<<"Final Loss function FVAL "<<FVAL<<endl;
	return 0;
}

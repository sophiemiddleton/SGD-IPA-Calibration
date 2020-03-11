// To run :  g++ -std=c++11 SGDiterativev0.cpp -o a.out

#include <iostream>
#include <math.h>
#include <random>
#include <cstdlib>
#include <fstream>
#include <chrono> 
#include <ctime>
using namespace std;

double FVAL;
int N_EVENTS;
int N_CONVERGED = 0;
int N_CRYSTALS = 674;
int N_TOTALHITS = 0;
double Loss = 0;
double step_size = 0.00035;
double error =1;
double MaxIterations = 1000;
double MaxFunction = 100;
double dcmin = 0.00001;
double tolerance = 0.1;

std::vector<double> CalibrationConstants;


struct CrystalList{
	std::vector<double> crystal_energy;
	std::vector<int> crystal_number;
   	CrystalList();
     	CrystalList(std::vector<double> _crystal_energy, std::vector<int> _crystal_number) : crystal_energy(_crystal_energy), crystal_number(_crystal_number) {};
};

struct Event{
    unsigned int EventNumber;
	double track_energy;
	int cluster_size;
	CrystalList crystal_list;
	Event();
	Event(unsigned int _n, double _energy, double _size) : EventNumber(_n), track_energy(_energy), cluster_size(_size){};
	Event(unsigned int _n, double _energy, double _size, CrystalList _list ) : EventNumber(_n), track_energy(_energy), cluster_size(_size), crystal_list(_list){};
	
};

std::vector<Event> FakeDateMaker(std::vector<double> RawCalibrationResults, std::vector<double> offset_vector){
    
	std::vector<Event> event_list;

	double offset;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::normal_distribution<double> te(46.,3.); 
	std::normal_distribution<double> cs(4.,1);
	std::uniform_real_distribution<double> cn(0, N_CRYSTALS);
	std::uniform_real_distribution<double> randoff(0.1, 1.0);

	for(size_t n=0;n<N_EVENTS;n++){
        
        	std::cout<<"[In FakeDataMaker()] Event Loop ..."<<n<<std::endl;
		
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
         
		Event event(n, track_energy,cluster_size,crystal_list);
        
        	event_list.push_back(event);
     }
    
    return event_list;
  
}

void SetOffsetVector(std::vector<double> &RawCalibrationResults, std::vector<double> &offset_vector){
    
	std::vector<Event> event_list;

	double offset;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> randoff(0.1, 1.0);

	for(int c=0;c<N_CRYSTALS;c++){

		std::cout<<"[In OffsetMaker()] Finding Raw Values ..."<<std::endl;
		auto const off = randoff(mt);
		offset_vector.push_back(off);
		RawCalibrationResults.push_back(off*0.8);
	 }
}

CrystalList FillCrystals(const char *CrystalsFile, unsigned int &n, int eventN){
	std::cout<<"[In ReadInputData()] : opening "<<CrystalsFile<<std::endl;
	
    	FILE *fC = fopen(CrystalsFile, "r");
	if ( fC == NULL) {
		std::cout<<"[In ReadDataInput()] : Error: Cannot open Cry files "<<std::endl;
		exit(1);
	}
	
    	float EoP, E, P, CryE;
    	int eventT, runT, eventC, runC, CryId;

    	std::vector<double> crystal_energy;
	std::vector<int> crystal_number;
    	n=0;
	while(fscanf(fC, "%i,%i,%i,%f\n", &eventC, &runC, &CryId, &CryE)!=EOF){
            if (eventC == eventN){
                crystal_energy.push_back(CryE);
                crystal_number.push_back(CryId);
                cout<<eventN<<" "<<CryId<<endl;
                n++;
            }
	}
    	CrystalList crystal_list(crystal_energy, crystal_number);
	fclose(fC);
    	return crystal_list;
}


std::vector<Event> BuildEventsFromData(const char *CrystalsFile, const char *TrackFile){
	std::cout<<"[In ReadInputData()] : opening "<<TrackFile<<std::endl;
	std::vector<Event> event;
	FILE *fT = fopen(TrackFile, "r");
   
	if (fT == NULL ) {
		std::cout<<"[In ReadDataInput()] : Error: Cannot open Track files "<<std::endl;
		exit(1);
	}
	
    	float EoP, E, P;
    	int eventT, runT;
    
	while(fscanf(fT, "%i,%i,%f,%f,%f\n", &eventT, &runT, &EoP, &E, &P)!=EOF){
            unsigned int clSize;
            CrystalList crystal_list = FillCrystals(CrystalsFile, clSize, eventT);
	    Event e(eventT, E, clSize ,crystal_list);
            event.push_back(e);
	}
	std::cout<<"[In ReadInputData()] : closing "<<TrackFile<<" and "<<CrystalsFile<<std::endl;
    	fclose(fT);
	return event;
}

unsigned int GetLines (const char *filename){
    ifstream file(filename);
    unsigned int lines=0;
    std::string line;
    while (getline(file, line)){
	   lines++;
    }
    std::cout<<"[In GetLines()] File has "<<lines<<" Lines"<<std::endl;
    return lines;
}

double FullF(std::vector<Event> event_list, std::vector<double> constants){
    int sum = 0;
    double F = 0;
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

double F(Event event, std::vector<double> constants){
    double F = 0;
    int sum = 0;
        for(unsigned int c=0;c < event.cluster_size;c++){
		int cry_i = event.crystal_list.crystal_number[c];
        sum += constants[cry_i]*event.crystal_list.crystal_energy[c];
	}
    F+= pow((sum - event.track_energy)*(1/error) ,2);
	return F;
}

std::vector<double> SGD(Event event, int j, std::vector<double> constants, double& FVAL){
	std::cout<<"[In SGD ()] Beginning ..."<<std::endl;
	N_EVENTS =  GetLines("Tracks.csv");
	double Loss = FVAL;
	double old_c ;
	double new_c ;
	double dc;
	double dFdCi;
	bool converged =false;

	double InitLoss = F(event, constants);
	std::cout<<"[In SGD()] Initial Loss is "<<InitLoss<<std::endl;
	std::vector<double> previous_constants = constants;
	unsigned int k = 0;
	double Etrk = event.track_energy;  

	while(converged == false and k < MaxIterations){
       
        for(unsigned int m=0; m<event.cluster_size;m++){
            int Cm = event.crystal_list.crystal_number[m];
            old_c = constants[Cm]; 
            double Vm=event.crystal_list.crystal_energy[m]; 
            double prediction =0;
            
            for(unsigned int i=0;i<event.cluster_size;i++){
                int Ci = event.crystal_list.crystal_number[i];
                prediction +=constants[Ci]*event.crystal_list.crystal_energy[i];
             }
            
            dFdCi = 2*Vm*(1/pow(error,2))*(prediction -Etrk);

            new_c = (old_c - step_size*constants[Cm]*dFdCi);

            dc = abs(new_c - old_c); 
           
            if(dc > dcmin ){
                constants[Cm] = new_c; 
            }
        }
        Loss = F(event, constants);
        
        if((Loss <  InitLoss) and (Loss < MaxFunction)){
            converged =true;
            N_CONVERGED +=1;
        }

         k++;
       
    }
    
    if(converged ==true){
        
          for(unsigned int m=0; m<event.cluster_size;m++){
              
                int Ci = event.crystal_list.crystal_number[m];
                cout<<"[In SGD() ] Updated "<<Ci<<" from "<<previous_constants[Ci]<<" to "<<constants[Ci]<<" with k iterations "<<k<<endl;
                if (m==event.cluster_size-1) Loss = F(event,constants); 
          }
        
        return constants;
    } else {
        return previous_constants;
    } 
}

void Randomize(std::vector<Event> &EventList){
    std::random_shuffle (EventList.begin(), EventList.end());
}

int main(){
	bool fake = false;
	std::cout<<"[In Main()] Beginning ..."<<std::endl;
	ofstream outputfile, TrackFile, CrystalsFile;
	outputfile.open("SGDv1.csv");

	outputfile<<"cryId,reco,true,Residual"<<std::endl;

	std::vector<double> offset_vector;
	std::vector<double> RawCalibrationResults;

	for(int c=0;c<N_CRYSTALS;c++){

		CalibrationConstants.push_back(0);

	}
    
 	SetOffsetVector(RawCalibrationResults, offset_vector);
     	std::vector<Event> event_list;
     	if(fake) event_list = FakeDateMaker(RawCalibrationResults, offset_vector);
    	if(!fake){
        
		event_list = BuildEventsFromData("Crystals.csv","Tracks.csv");
		std::cout<<"Found "<<event_list.size()<<" Events "<<std::endl;
    	 }
	auto start = chrono::high_resolution_clock::now();
	CalibrationConstants = RawCalibrationResults;

	for(unsigned int l = 0; l < 1 ; l++){
		int N_CONVERGED = 0;
		for(auto const& event : event_list){

		Randomize(event_list);
		CalibrationConstants = SGD(event, event.EventNumber, CalibrationConstants,  FVAL);

         	}
       
     	}
	for(int i =0 ;i<N_CRYSTALS;i++){
        
        std::cout<<"constant for crystal "<<i<<" is "<<CalibrationConstants[i]<<"True Offset is "<<offset_vector[i]<<" Residuals "<<CalibrationConstants[i]-offset_vector[i]<<std::endl;
        
        outputfile<<i<<","<<CalibrationConstants[i]<<","<<offset_vector[i]<<","<<CalibrationConstants[i]-offset_vector[i]<<std::endl;
            
    	}
    
	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(end - start); 

	double endLoss = FullF(event_list, CalibrationConstants);
	cout<<"NEvents Processed "<<N_EVENTS<<" NEVents converged "<<N_CONVERGED<<"Time "<<duration.count()<<" Final Loss function "<<endLoss<<endl;

	return 0;
	}

#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <fstream>
#include <random>
#include <ctime>
//build:
//under HW4 folder
//g++ -I../Eigen HW4.cpp -o HW4

//execute:
//HW4 problem_number
//for example, to solve Problem_1_1: HW4 1
#include "GA.h"

double Problem_1_1(const Eigen::VectorXd &x){
    return x(0)*sin(x(0))+x(1)*sin(5*x(1));
}

double Problem_1_2_a(const Eigen::VectorXd &x){    
    return -exp(0.2*sqrt(pow(x(0)-1,2)+pow(x(1)-1,2))+cos(2*x(0))+sin(2*x(0)));
}

double Problem_1_2_b(const Eigen::VectorXd &x){
    double dotx=x.dot(x);
    return 0.5+(sin(dotx)-0.5)/(1+0.1*dotx);
}

Eigen::VectorXd Decode_1_2_b(const Eigen::VectorXd &x){
    Eigen::VectorXd r(x.size());
    for(int i=0;i<x.size();i++){
        r(i)=tan(x(i));
        //atan project infinite range to +-pi/2, so use tan to decode
    }
    return r;
}
int main(int argc, char * argv[])
{
    int p=3;
    if(argc>=2){
        p=atoi(argv[1]);
    }
    if(p==1){
        std::cout<<"============="<<std::endl;
        std::cout<<"=Problem_1_1="<<std::endl;
        std::cout<<"============="<<std::endl;
        int Dim=2;
        int ChromoNum=100;
        Eigen::VectorXd MinX(Dim),MaxX(Dim);
        MinX<<0,4;
        MaxX<<10,6;
        
        double CrossOverRate=0.5;
        double CrossOverAlpha=0.386;
        double MutationRate=0.2;
        double MutationParam=20;
        
        
        GA ga(&Problem_1_1, Dim, ChromoNum,
                MinX,MaxX,
                CrossOverRate, CrossOverAlpha,
                MutationRate, MutationParam);
        ga.Solve(20,true);
        ga.ShowResult();
    }
    else if(p==2){
        std::cout<<"==============="<<std::endl;
        std::cout<<"=Problem_1_2_a="<<std::endl;
        std::cout<<"==============="<<std::endl;
        int Dim=2;
        int ChromoNum=100;
        Eigen::VectorXd MinX(Dim),MaxX(Dim);
        MinX<<-5,-5;
        MaxX<<5,5;
        
        double CrossOverRate=0.5;
        double CrossOverAlpha=0.386;
        double MutationRate=0.2;
        double MutationParam=20;
        
        
        GA ga(&Problem_1_2_a, Dim, ChromoNum,
                MinX,MaxX,
                CrossOverRate, CrossOverAlpha,
                MutationRate, MutationParam);
        ga.Solve(20,true);
        ga.ShowResult();
    }
    else{
        std::cout<<"==============="<<std::endl;
        std::cout<<"=Problem_1_2_b="<<std::endl;
        std::cout<<"==============="<<std::endl;
        int Dim=2;
        int ChromoNum=100;
        Eigen::VectorXd MinX(Dim),MaxX(Dim);
        MinX<<-M_PI/2,-M_PI/2;
        MaxX<<M_PI/2,M_PI/2;
        
        double CrossOverRate=0.5;
        double CrossOverAlpha=0.386;
        double MutationRate=0.2;
        double MutationParam=20;
        
        
        GA ga(&Problem_1_2_b, Dim, ChromoNum,
                MinX,MaxX,
                CrossOverRate, CrossOverAlpha,
                MutationRate, MutationParam);
        ga.SetDecoder(Decode_1_2_b);
        ga.Solve(20,true);
        ga.ShowResult();
    }
}
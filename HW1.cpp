#include <iostream>
#include <cmath>
#include <iomanip>

//g++ -o HW1 HW1.cpp
float F_7_2(float x){
    return x*x+4*cos(x);
}
float dF_7_2(float x){
    return 2*x-4*sin(x);
}
float ddF_7_2(float x){
    return 2-4*cos(x);
}


float GoldenSectionSearch(float (*F)(float),float a,float b,float threshold,int max_iteration,bool vb){
    float p=0.618;
    float q=1-p;
    float last_a=a;
    float last_b=b;
    float a_=p*last_a+q*last_b;
    float b_=p*last_b+q*last_a;
    
    float fa_=(*F)(a_);
    float fb_=(*F)(b_);
    if(vb){
        std::cout<<"GSS"<<std::endl;
        std::cout<<"================================================"<<std::endl;
        std::cout<<"Iter\tak\tbk\tF(ak)\tF(bk)\tNewInterval"<<std::endl;
        std::cout<<"================================================"<<std::endl;
        std::cout<<"1\t"<<last_a<<"\t"<<last_b<<"\t"<<fa_<<"\t"<<fb_<<"\t"<<"["<<a_<<","<<b_<<"]"<<std::endl;
    }
    if(fabs(b_-a_)<threshold){
        return (a_+b_)/2;
    }
    
    for(int i=2;i<max_iteration;i++){
        if(fa_>fb_){
            last_a=a_;
            a_=b_;
            fa_=fb_;
            b_=p*last_b+q*last_a;
            fb_=(*F)(b_);
        }
        else{
            last_b=b_;
            b_=a_;
            fb_=fa_;
            a_=p*last_a+q*last_b;
            fa_=(*F)(a_);
        }
        if(vb)std::cout<<i<<"\t"<<last_a<<"\t"<<last_b<<"\t"<<fa_<<"\t"<<fb_<<"\t"<<"["<<a_<<","<<b_<<"]"<<std::endl;
        if(fabs(b_-a_)<threshold){
            if(vb)std::cout<<"x*="<<(a_+b_)/2<<std::endl<<std::endl;
            return (a_+b_)/2;
        }
    }
    if(vb)std::cout<<"GG"<<std::endl;
    return (a+b)/2;

}

float Newton(float (*F)(float),float (*dF)(float),float (*ddF)(float),float initial,float threshold,float learning_rate,int max_iteration,bool vb){
    float x=initial;
    float dx=(*dF)(x)/(*ddF)(x);
    if(vb){
        std::cout<<"Newton\tinitial="<<x<<std::endl;
        std::cout<<"================================================"<<std::endl;
        std::cout<<"Iter\tx\tdx\tF(x)"<<std::endl;
        std::cout<<"================================================"<<std::endl;
        std::cout<<"0\t"<<x<<"\t"<<dx<<"\t"<<(*F)(x)<<std::endl;
    }
    for(int i=0;i<max_iteration;i++){
        x=x-learning_rate*dx;
        dx=(*dF)(x)/(*ddF)(x);
        std::cout<<i+1<<"\t"<<x<<"\t"<<dx<<"\t"<<(*F)(x)<<std::endl;
        if(fabs(dx)<threshold){
            if(vb)std::cout<<"x*="<<x<<std::endl<<std::endl;
            return x;
        }
    }
    if(vb)std::cout<<"GG"<<std::endl;
    return x;
}

float GoldenSectionSearchLearningRate(float (*F)(float),float a,float b,float threshold,int max_iteration,bool vb,float x,float dx){
    float p=0.618;
    float q=1-p;
    float last_a=a;
    float last_b=b;
    float a_=p*last_a+q*last_b;
    float b_=p*last_b+q*last_a;
    
    float fa_=(*F)(x-a_*dx);
    float fb_=(*F)(x-b_*dx);
    if(vb){
        std::cout<<"GSS"<<std::endl;
        std::cout<<"================================================"<<std::endl;
        std::cout<<"Iter\tak\tbk\tF(ak)\tF(bk)\tNewInterval"<<std::endl;
        std::cout<<"================================================"<<std::endl;
        std::cout<<"1\t"<<last_a<<"\t"<<last_b<<"\t"<<fa_<<"\t"<<fb_<<"\t"<<"["<<a_<<","<<b_<<"]"<<std::endl;
    }
    if(fabs(b_-a_)<threshold){
        return (a_+b_)/2;
    }
    
    for(int i=2;i<max_iteration;i++){
        if(fa_>fb_){
            last_a=a_;
            a_=b_;
            fa_=fb_;
            b_=p*last_b+q*last_a;
            fb_=(*F)(x-b_*dx);
        }
        else{
            last_b=b_;
            b_=a_;
            fb_=fa_;
            a_=p*last_a+q*last_b;
            fa_=(*F)(x-a_*dx);
        }
        if(vb)std::cout<<i<<"\t"<<last_a<<"\t"<<last_b<<"\t"<<fa_<<"\t"<<fb_<<"\t"<<"["<<a_<<","<<b_<<"]"<<std::endl;
        if(fabs(b_-a_)<threshold){
            if(vb)std::cout<<"x*="<<(a_+b_)/2<<std::endl<<std::endl;
            return (a_+b_)/2;
        }
    }
    if(vb)std::cout<<"GG"<<std::endl;
    return (a+b)/2;

}


float NewtonGSS(float (*F)(float),float (*dF)(float),float (*ddF)(float),
    float Newton_initial,float Newton_threshold,
    float LR_lower,float LR_upper,float LR_threshold,
    int max_iteration){

    float x=Newton_initial;
    float dx=(*dF)(x)/(*ddF)(x);
    float learning_rate=GoldenSectionSearchLearningRate(F,LR_lower,LR_upper,LR_threshold,max_iteration,false,x,dx);
    
    std::cout<<"NewtonGSS\tinitial="<<x<<std::endl;
    std::cout<<"================================================"<<std::endl;
    std::cout<<"Iter\tx\tdx\tF(x)\tLR"<<std::endl;
    std::cout<<"================================================"<<std::endl;
    std::cout<<"0\t"<<x<<"\t"<<dx<<"\t"<<(*F)(x)<<"\t"<<learning_rate<<std::endl;
    for(int i=0;i<max_iteration;i++){
        x=x-learning_rate*dx;
        dx=(*dF)(x)/(*ddF)(x);
        learning_rate=GoldenSectionSearchLearningRate(F,LR_lower,LR_upper,LR_threshold,max_iteration,false,x,dx);
        std::cout<<i+1<<"\t"<<x<<"\t"<<dx<<"\t"<<(*F)(x)<<"\t"<<learning_rate<<std::endl;
        if(fabs(dx)<Newton_threshold){
            std::cout<<"x*="<<x<<std::endl<<std::endl;
            return x;
        }
    }
    std::cout<<"GG"<<std::endl;
    return x;
}

int main(){
    std::cout << std::setprecision(4);
    std::cout << std::fixed;
    std::cout<<"7.2.b"<<std::endl;
    float x=GoldenSectionSearch(&F_7_2,1,2,0.05,10,true);
    
    std::cout<<"7.2.d"<<std::endl;
    x=Newton(&F_7_2,&dF_7_2,&ddF_7_2,1.0,0.05,1,10,true);
    x=Newton(&F_7_2,&dF_7_2,&ddF_7_2,1.1,0.05,1,10,true);
    
    std::cout<<"7.2.extra"<<std::endl;
    NewtonGSS(&F_7_2,&dF_7_2,&ddF_7_2,
        1,0.05,
        0,1,0.2,
        10);
    NewtonGSS(&F_7_2,&dF_7_2,&ddF_7_2,
        1.1,0.05,
        0,1,0.2,
        10);
    return 0;
}

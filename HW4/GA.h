#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <random>
class GA{
    public:
        Eigen::MatrixXd chromo;
        Eigen::VectorXd minx;
        Eigen::VectorXd maxx;
        int dim;
        int chromonum;
        double crossoverrate;
        double alpha;//0~1 : higher -> offspring vary more
        double mutationrate;
        double mp;//about 20~100 : lower -> expand mutation range
        
        bool usedecoder;
        Eigen::VectorXd (*decode)(const Eigen::VectorXd &x);
        
        double (*f)(const Eigen::VectorXd &x);
        std::default_random_engine re;
        
        GA(double (*F)(const Eigen::VectorXd &x),int Dim,int ChromoNum,
            const Eigen::VectorXd MinX,const Eigen::VectorXd MaxX,
            double CrossOverRate,double CrossOverAlpha,
            double MutationRate,double MutationParam);
        void Init();
        void Solve(int Iteration,bool verbose);
        void Mutation();
        void Crossover();
        void Selection();
        void ShowResult();
        void SetDecoder(Eigen::VectorXd (*Decode)(const Eigen::VectorXd &x));
};


GA::GA(double (*F)(const Eigen::VectorXd &x),int Dim,int ChromoNum,
            const Eigen::VectorXd MinX,const Eigen::VectorXd MaxX,
            double CrossOverRate,double CrossOverAlpha,
            double MutationRate,double MutationParam){
    dim=Dim;
    chromonum=ChromoNum;
    crossoverrate=CrossOverRate;
    mutationrate=MutationRate;
    minx=MinX;
    maxx=MaxX;
    f=F;
    alpha=CrossOverAlpha;
    mp=MutationParam;
    Init();
    usedecoder=false;
}

void GA::SetDecoder(Eigen::VectorXd (*Decode)(const Eigen::VectorXd &x)){
    decode=Decode;
    usedecoder=true;
}

void GA::Init(){
    chromo=Eigen::MatrixXd(dim,chromonum);
    // to obtain one chromo: chromo.col(j)
    for(int i=0;i<dim;i++){
        std::uniform_real_distribution<double> rd(minx(i),maxx(i));
        for(int j=0;j<chromonum;j++){
            chromo(i,j) = rd(re);
        }
    }
    
}
void GA::Solve(int Iteration,bool verbose){
    for(int i=0;i<Iteration;i++){
        Selection();

        Crossover();

        Mutation();
        
        if(verbose){
            std::cout<<"iteration="<<i<<std::endl;
            ShowResult();
        }

    }
}

void GA::Selection(){
    std::uniform_int_distribution<> indexgen(0, chromonum-1);
    for(int i=0;i<chromonum;i++){
        int idx=indexgen(re);
        //minimize
        double vidx,vi;
        if(usedecoder){
            vidx=(*f)((*decode)(chromo.col(idx)));
            vi=(*f)((*decode)(chromo.col(i)));
        }
        else{
            vidx=(*f)(chromo.col(idx));
            vi=(*f)(chromo.col(i));
        }
        if(vidx<vi){
            chromo.col(i)=chromo.col(idx);
        }
    }
}

void GA::Crossover(){
    //Blend crossover alpha 
    std::uniform_int_distribution<> indexgen(0, chromonum-1);
    std::uniform_real_distribution<double> rd(0,1);
    for(int i=0;i<chromonum;i++){
        if(rd(re)>crossoverrate)continue;
        int idx1,idx2;
        idx1=indexgen(re);
        idx2=indexgen(re);
        
        Eigen::VectorXd x1=chromo.col(idx1);
        Eigen::VectorXd x2=chromo.col(idx2);
        Eigen::VectorXd d(dim);
        Eigen::VectorXd min12(dim); //min max of x1 x2
        Eigen::VectorXd max12(dim);
        
        for(int j=0;j<dim;j++){
            if(x1(j)<x2(j)){
                min12(j)=x1(j);
                max12(j)=x2(j);
            }
            else{
                min12(j)=x2(j);
                max12(j)=x1(j);
            }
        }
        d=max12-min12;
        min12=min12-alpha*d;
        max12=max12+alpha*d;
        
        for(int j=0;j<dim;j++){
            std::uniform_real_distribution<double> rdof(min12(j),max12(j));
            chromo(j,idx1)=rdof(re);
            if(chromo(j,idx1)<minx(j)){
                chromo(j,idx1)=minx(j);
            }
            else if(chromo(j,idx1)>maxx(j)){
                chromo(j,idx1)=maxx(j);
            }
            chromo(j,idx2)=rdof(re);
            if(chromo(j,idx2)<minx(j)){
                chromo(j,idx2)=minx(j);
            }
            else if(chromo(j,idx2)>maxx(j)){
                chromo(j,idx2)=maxx(j);
            }
        }
    }
    
}
void GA::Mutation(){
    //Polynomial Mutation
    std::uniform_int_distribution<> indexgen(0, chromonum-1);
    std::uniform_real_distribution<double> rd(0,1);
    for(int i=0;i<chromonum;i++){
        if(rd(re)>mutationrate)continue;
        int idx=indexgen(re);
        for(int j=0;j<dim;j++){
            double u=rd(re);
            if(u>0.5){//upper
                double deltaU=1-pow(2.0*(1.0-u),1.0/(1.0+mp));
                chromo(j,idx)=chromo(j,idx)+deltaU*(maxx(j)-chromo(j,idx));
            }
            else{//lower
                double deltaL=pow(2*u,1.0/(1.0+mp))-1;
                chromo(j,idx)=chromo(j,idx)+deltaL*(chromo(j,idx)-minx(j));
            }
        }
    }
}

void GA::ShowResult(){
    double m=usedecoder?(*f)((*decode)(chromo.col(0))) : (*f)(chromo.col(0));
    int idx=0;
    for(int i=1;i<chromonum;i++){
        double v=usedecoder?(*f)((*decode)(chromo.col(i))) : (*f)(chromo.col(i));

        if(v<m){
            m=v;
            idx=i;
        }
    }
    std::cout<<"min value:"<<m<<std::endl;
    if(usedecoder){
        std::cout<<"\tx:"<<(*decode)(chromo.col(idx)).transpose()<<std::endl;
    }
    else{
        std::cout<<"\tx:"<<chromo.col(idx).transpose()<<std::endl;
    }
    std::cout<<"------------------------------------------"<<std::endl;
}

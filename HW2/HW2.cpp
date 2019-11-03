#include <iostream>
#include <cmath>
//under HW2 folder
//g++ -I../Eigen HW2.cpp -o HW2
#include <Eigen/Core>

Eigen::MatrixXd GramSchmit(const Eigen::MatrixXd &Q){
    int n = Q.rows();
    Eigen::MatrixXd CD = Eigen::MatrixXd::Identity(n,n);//CD=conjugate directions
    if(n != Q.cols()){return CD;}
    for(int i=0;i<n;i++){
        Eigen::VectorXd projection = Eigen::VectorXd::Zero(n);
        for(int j=0;j<i;j++){
            Eigen::VectorXd Qj=Q*CD.col(j);
            double scalar = CD.col(i).dot(Qj)/CD.col(j).dot(Qj);
            projection = projection + scalar * CD.col(j);
        }
        CD.col(i) = CD.col(i) - projection;
    }
    return CD;
}

double Quadratic(const Eigen::VectorXd &x,const Eigen::MatrixXd &Q, const Eigen::VectorXd &b){
    return x.dot(Q*x)+x.dot(b);
}
Eigen::VectorXd dQuadratic(const Eigen::VectorXd &x,const Eigen::MatrixXd &Q, const Eigen::VectorXd &b){
    return Q*x+b;
}

Eigen::VectorXd ConjugateGradientQuadratic(
    double (*F)(const Eigen::VectorXd &,const Eigen::MatrixXd &, const Eigen::VectorXd &),
    Eigen::VectorXd (*dF)(const Eigen::VectorXd &,const Eigen::MatrixXd &, const Eigen::VectorXd &),
    const Eigen::MatrixXd &Q,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &initial)
{
    
    int n = initial.size();
    Eigen::VectorXd x = initial;
    Eigen::MatrixXd CD = GramSchmit(Q);
    for(int i=0;i<n;i++){
        Eigen::VectorXd g = (*dF)(x,Q,b);
        if(g.norm()<0.01){return x;}
        Eigen::VectorXd d = CD.col(i);
        double a = - g.dot(d)/d.dot(Q*d);
        x = x + a*d; 
    }
    return x;
}


int main(){
    Eigen::Matrix2d Q;
    Q << 5, -3, -3, 2;
    Eigen::Vector2d b;
    b << 0, -1;
    Eigen::Vector2d initial;
    // xQx/2+bx
    initial<<0,0;
    Eigen::Vector2d x = ConjugateGradientQuadratic(&Quadratic,&dQuadratic,Q,b,initial);
    std::cout<<"x* = "<<std::endl<<x;
}
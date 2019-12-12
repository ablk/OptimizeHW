#include <iostream>
#include <cmath>
//under HW3 folder
//g++ -I../Eigen HW3.cpp -o HW3
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <fstream>
int main(int argc, char * argv[])
{
    std::fstream fin("denoise.txt", std::ios_base::in);
    std::vector<double> corrupted_data;
    
    double lambda=0.1;
    if(argc>=2){
        lambda=atof(argv[1]);
    }
    
    float d;
    while (fin >> d)
    {
        corrupted_data.push_back(d);
    }
    Eigen::VectorXd Xcor=Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(corrupted_data.data(), corrupted_data.size());
    //std::cout<<Xcor<<std::endl;
    //std::cout<<Xcor.size();
    int Dim=Xcor.size();
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(Dim*2-2);
    for(int i=0;i<Dim-1;i++)
    {
        tripletList.push_back(Eigen::Triplet<double>(i,i,-1));
        tripletList.push_back(Eigen::Triplet<double>(i,i+1,1));
    }
    
    Eigen::SparseMatrix<double> D(Dim-1,Dim);
    D.setFromTriplets(tripletList.begin(), tripletList.end());
    
    Eigen::SparseMatrix<double> I(Dim,Dim);
    I.setIdentity();
    Eigen::SparseMatrix<double> A=I+lambda*D.transpose()*D;
    
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(A);
    Eigen::VectorXd X=solver.solve(Xcor);
    
    std::ostringstream ss;
    ss << "filtered" <<lambda<<".txt";
    std::string foutname(ss.str());
    
    std::fstream fout(foutname, std::ios_base::out);
    for(int i=0;i<X.size();i++){
        fout<<X(i)<<std::endl;
    }
    return 0;
}
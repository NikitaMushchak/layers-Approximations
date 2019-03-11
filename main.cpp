#include <iostream>
#include <vector>
#include "ai.hh"



bool desc (double i,double j) { return (i<j); }
void ApproximateLayers(std::vector<std::vector<double> >&layers){

        std::size_t num;
        double H ;
        double dx;

        for(std::size_t i = 1 ; i < layers.size(); ++i){
                if(layers[i][0]*layers[i][1] < 0.){
                    H = layers[i][1] - layers[i][0];
                    num = i;
                    break;
                }

        }
        dx = 0.09 * H;
        double sigma = layers[num][2];
        //std::cout<<"sigma = "<<sigma<<std::endl;
        std::cout<<"Heigth of centeral layer = "<<H<<"    number = "<<num<<std::endl;
        std::vector<std::size_t> sloi;

        double h;
        double hnew;
        if(H < 11.2){
            std::cout<<"H less than 11.2!!!!"<<std::endl;
            //H = 11.2;

            std::size_t it;
            it = num ;

            sloi.push_back(it);

            h = layers[it][1] - layers[it][0];
            it--;

            while( 11.2/2. > std::abs(layers[it][1])){
                //ai::printMarker();
                //std::cout<<"std::abs(layers[it][0]))= "<<std::abs(layers[it][0])<<std::endl;

                hnew = layers[it][1] - layers[it][0];

                sigma = (sigma * h + hnew * layers[it][2])/(h+hnew);
                h = h+hnew;
                //std::cout<<"h = "<<h<<"    h new = "<<hnew<<"  sigma = "<<sigma<<std::endl;

                sloi.push_back(it);
                it--;
            }


            it = num+1;

            while(11.2/2. > layers[it][0]){
                //ai::printMarker();
                //std::cout<<"std::abs(layers[it][0]))= "<<std::abs(layers[it][0])<<std::endl;

                hnew = layers[it][1] - layers[it][0];

                sigma = (sigma * h + hnew * layers[it][2])/(h+hnew);
                h += hnew;
                //std::cout<<"h = "<<h<<"    h new = "<<hnew<<"  sigma = "<<sigma<<std::endl;

                sloi.push_back(it);
                it++;

            }
            H = 11.2;
            dx = 0.09 * H;
        }

        if(H > 11.2){

        std::size_t it;
        it = num ;

        sloi.push_back(it);

        h = layers[it][1] - layers[it][0];
        it--;

        //while(  11.2/2. > std::abs(layers[it][1])){
            //ai::printMarker();
            //std::cout<<"std::abs(layers[it][0]))= "<<std::abs(layers[it][0])<<std::endl;

            hnew = layers[it][1] - layers[it][0];

            sigma = (sigma * h + hnew * layers[it][2])/(h+hnew);
            h = h+hnew;
            //std::cout<<"h = "<<h<<"    h new = "<<hnew<<"  sigma = "<<sigma<<std::endl;

            sloi.push_back(it);
            it--;
        //}


        it = num+1;

        //while(11.2/2. > layers[it][0]){
            //ai::printMarker();
            //std::cout<<"std::abs(layers[it][0]))= "<<std::abs(layers[it][0])<<std::endl;

            hnew = layers[it][1] - layers[it][0];

            sigma = layers[num][2];
            h += hnew;
            std::cout<<"h = "<<h<<"    h new = "<<hnew<<"  sigma = "<<sigma<<std::endl;

            sloi.push_back(it);
            it++;

        //}
    }

        ai::printVector(sloi);

        std::size_t min = ai::min(sloi);
        std::size_t max = ai::max(sloi);
        ai::printMarker();
        std::vector<std::vector<double> > layers1;

        for(size_t i = 0; i < layers.size(); i++ ){
            layers1.push_back(std::vector<double>{layers[i][0] ,layers[i][1], layers[i][2], layers[i][3], layers[i][4], layers[i][5]} );
            if(i == min){
                layers1.push_back(std::vector<double>{-H/2. ,H/2., sigma, layers[i][3], layers[i][4], layers[i][5]});
                i = max-1;
            }

        }
        // Обрабатываем пересечение слоев


        for(size_t i = 0 ; i < layers1.size();i++){
            if (layers1[i][0]*layers1[i][1] < 0.) num = i;
        }

        for(size_t i = 0; i < layers1.size()-1; ++i){
            if(layers1[i][1] > layers1[i+1][0]){

                if(i!=num)
                    layers1[i][1] = layers1[i+1][0];

                if(i==num)
                    layers1[i+1][0]= layers1[i][1];

            }
            // if(std::abs(layers[i][1]) > std::abs(layers[i+1][0]) ){
            //     layers1[i+1][0]= layers1[i][1];
            // }
        }
        std::vector<std::vector<double> > mesh;
        //Обсчтываем все сверху продуктивного

        //это прокуктивный слой
        mesh.push_back(std::vector<double>{layers1[num][0], layers1[num][1],

        layers1[num][2], layers1[num][3], layers1[num][4], layers1[num][5]} );
        double up = layers1[num][0];

        while(std::abs(layers1[0][0]) > std::abs(up)){
            mesh.push_back(std::vector<double>{up-dx, up, 0, 0 ,0 , 0});
            up-=dx;
        }

    std::cout<<"STEP = "<<dx<<std::endl;
        // std::cout<<"mesh:"<<std::endl;
        // ai::printMatrix(mesh);

        std::vector<std::vector<double> > mesh1;
        for(size_t i = mesh.size()-1 ; i>0; --i){
            mesh1.push_back(std::vector<double>{mesh[i][0], mesh[i][1], mesh[i][2], mesh[i][3], mesh[i][4], mesh[i][5] });
        }
        mesh1.push_back(std::vector<double>{layers1[num][0], layers1[num][1],

        layers1[num][2], layers1[num][3], layers1[num][4], layers1[num][5]});




         up = layers1[num][1];
         //std::cout<<"up = "<<up<<std::endl;
         size_t f = layers1.size()-1;
         //ai::printMarker();
        while(std::abs(layers1[f][1]) > std::abs(up)){
            //ai::printMarker();
            mesh1.push_back(std::vector<double>{up, up+dx, 0, 0 ,0 , 0});
            up+=dx;
        }
        mesh.clear();



        double eps = pow(10, -8);

        //идем по строкам mesh1
        // for(size_t i = 0 ;i<mesh1.size(); ++i){
        //     if(mesh1[i][2] < eps){
        //         for(size_t j = 0 ; j < layers1.size(); ++j){
        //             if((std::abs(mesh1[i][0]) <= std::abs(layers1[j][0])) && (std::abs(mesh1[i][1]) >= std::abs(layers1[j][1])) )
        //                 mesh1[i][2] = layers1[j][2];
        //             }
        //         }
        // }
        size_t middle;
        for(size_t i = 0 ; i < mesh1.size(); ++i){
            if(mesh1[i][0]*mesh1[i][1] < 0.) {
                middle = i;
                break;
            }
        }
        std::cout<<"middle = "<<middle<<std::endl;
        double rastlayer;
        double rastmesh;
        double dl;
        double hadd;


        size_t j = middle -1;


        for(size_t i = num-1; i >= 0; --i){  //Идем по слоям layers1 i - итератор
            std::cout<<" "<<std::endl;
            std::cout<<"iteration = "<<i<<std::endl;
            h = layers1[i][1]-layers1[i][0];     //j - итератор по mesh1
            std::cout<<"h curent layer = "<<h<<std::endl;
            dl = std::ceil(h / dx );
            std::cout<<"Number of dx in layer  = "<<dl<<std::endl;
            if((int)dl == 1 ){
                std::cout<<"number of dx in layer = "<<dl<<std::endl;
                std::cout<<"j = "<<j<<std::endl;
                std::cout<<"i = "<<i<<std::endl;
                ai::printMarker();
                mesh1[j][2]+=h*layers1[i][2];
                mesh1[j][3]+=h*layers1[i][3];
                mesh1[j][4]+=h*layers1[i][4];
                mesh1[j][5]+=h*layers1[i][5];
                std::cout<<"mesh1 = "<<mesh1[j][2]<<std::endl;
                if(i>=0){
                    //ai::printMarker();
                    i--;
                    mesh1[j][2]+=(dx-h)*layers1[i][2];
                    mesh1[j][3]+=(dx-h)*layers1[i][3];
                    mesh1[j][4]+=(dx-h)*layers1[i][4];
                    mesh1[j][5]+=(dx-h)*layers1[i][5];

                    mesh1[j][2]/=dx;
                    mesh1[j][3]/=dx;
                    mesh1[j][4]/=dx;
                    mesh1[j][5]/=dx;
                }
                std::cout<<"mesh1 = "<<mesh1[j][2]<<std::endl;
                i++;
                j--;
            }

            if((int)dl > 1){
                int iter = 0;
                std::cout<<"j = "<<j<<std::endl;
                std::cout<<"i = "<<i<<std::endl;
                double dh;

                while(iter < (int)dl){
                    if(j>=0){
                        ai::printMarker();
                        // dh =
                        mesh1[j][2] = layers1[i][2];
                        mesh1[j][3] = layers1[i][3];
                        mesh1[j][4] = layers1[i][4];
                        mesh1[j][5] = layers1[i][5];
                        if(j==0){break;}

                        j--;

                    }
                    iter++;
                    std::cout<<"j = "<<j<<"   iter = "<<iter<<std::endl;
                    //if((int)dl == 1) j--;
                }
                // j++;
                //i++;
            }
            if (i==0){break;}
            //if (j==0){break;}
        }
    // std::cout<<"Up is fineshed"<<std::endl;
    // std::cout<<" "<<std::endl;
    // std::cout<<" "<<std::endl;
    // std::cout<<" "<<std::endl;
         j = middle + 1;


        for(size_t i = num+1; i < layers1.size(); ++i){  //Идем по слоям layers1 i - итератор
            // std::cout<<" "<<std::endl;
            // std::cout<<"iteration = "<<i<<std::endl;
            h = layers1[i][1]-layers1[i][0];     //j - итератор по mesh1
            // std::cout<<"h curent layer = "<<h<<std::endl;
            dl = std::floor(h / dx + 0.51);
            // std::cout<<"number of dx in layer = "<<dl<<std::endl;
            if( (int)dl == 1){
                // std::cout<<"i = " <<i<<"    j = "<<j<<std::endl;

                mesh1[j][2]+=h*layers1[i][2];
                mesh1[j][3]+=h*layers1[i][3];
                mesh1[j][4]+=h*layers1[i][4];
                mesh1[j][5]+=h*layers1[i][5];
                // std::cout<<"h= "<< h<<"*"<<"layers1[i][2]"<<layers1[i][2]<<" = "<<mesh1[j][2]<<std::endl;
                //std::cout<<"layers1[i][2] = "<<layers1[i][2]<<std::endl;

                if( i < layers1.size() ){
                    i++;
                    mesh1[j][2]+=(dx-h)*layers1[i][2];
                    mesh1[j][3]+=(dx-h)*layers1[i][3];
                    mesh1[j][4]+=(dx-h)*layers1[i][4];
                    mesh1[j][5]+=(dx-h)*layers1[i][5];

                    // std::cout<<"dx - h= "<< dx-h<<"*"<<"layers1[i][2] = "<<layers1[i][2]<<" = "<<(dx-h)*layers1[i][2]<<std::endl;
                }
                // std::cout<<"mesh1 = "<<mesh1[j][2]<<std::endl;
                i--;
                j++;
            }

            if((int)dl > 1){
                int iter = 0;
                // std::cout<<"j = "<<j<<std::endl;
                // std::cout<<"i = "<<i<<std::endl;
                // std::cout<<"number of layers = "<<(int)dl<<std::endl;

                h = 0;
                while(iter <= (int)dl){
                    if( j < mesh1.size() ) {
                        mesh1[j][2] = layers1[i][2];
                        mesh1[j][3] = layers1[i][3];
                        mesh1[j][4] = layers1[i][4];
                        mesh1[j][5] = layers1[i][5];

                            if (j==mesh1.size() - 1){break;}

                        j++;
                    }
                    iter++;
                    //std::cout<<"j = "<<j<<"   iter = "<<iter<<std::endl;
                }
                if (i==layers1.size()-1){break;}
                if (j==mesh1.size()-1){break;}
                //j++;
                //i++;
            }
        }

        std::vector<std::vector<double> > mesh2;
        double start , end;
        // ai::printMarker();
         for(size_t i = 0 ; i < mesh1.size()-1; ++i){
            // ai::printMarker();
             sigma = mesh1[i][2];
             end = mesh1[i][1];
             start = mesh1[i][0];
            while(mesh1[i][2] == mesh1[i+1][2]){
                //std::cout<<" i = "<<i<<std::endl;
                //std::cout<<"Sigma ="<<sigma<<std::endl;
                sigma = mesh1[i+1][2];
                end = mesh1[i+1][1];
                    i++;

                    if( i == mesh1.size()-1 ){ break;}

            }
            mesh2.push_back(std::vector<double>{start,  end, sigma, mesh1[i][3], mesh1[i][4], mesh1[i][5]} );


         }

         for(size_t i = 0 ; i < mesh2.size();++i){
             if(mesh2[i][1]*mesh2[i][0]<0.){
                 middle = i;
                 break;
             }
         }

         for(size_t i = middle ; i>1 ;--i){
            if(mesh2[i-1][2]-mesh2[i][2]>6.){
                mesh2[i-1][2]= 6. + mesh[i][2];
                std::cout<<">6"<<std::endl;
            }
         }

         for(size_t i = middle ; i<mesh2.size()-1 ;++i){
            if(mesh2[i+1][2]-mesh2[i][2]>6.){
                mesh2[i+1][2]= 6. + mesh[i][2];
                std::cout<<">6"<<std::endl;
            }
         }
         std::cout<<"mesh2:"<<std::endl;
         std::cout<<"middle = "<<middle<<std::endl;
         // layers.resize(mesh2.size());
         // for(size_)

         ai::printMatrix(mesh2);
         ai::saveMatrix("mesh2",mesh2);

         layers.clear();
         std::vector<std::vector<double> > layerss;
         // for(size_t i =0 ; i < mesh2.size(); ++i){
         //     layers[i].resize(mesh2.size());
         // }
         ai::printMarker();
         for(size_t i = mesh2.size() ; i > 0 ; ){
             --i;
             layerss.push_back(std::vector<double>{mesh2[i][0], mesh2[i][1], mesh2[i][2],mesh2[i][3],mesh2[i][4],mesh2[i][5]} );
             std::cout<<"iter = "<<i<<std::endl;
         }
         ai::saveMatrix("layerss", layerss);
        std::cout<<"mesh1:"<<std::endl;
        ai::printMatrix(mesh1);

        //ai::saveMatrix("mesh",mesh1);

        std::cout<<"New layers"<<std::endl;
        ai::printMatrix(layers1);

    }
    void ParseFromLayers(std::string filename, std::vector<std::vector<double> >& layers){
            std::vector<std::vector<double> > la;
            //ai::printMarker();
            //std::cout<<"la:"<<std::endl;
            ai::parseFileInMatrix(filename,' ',la);
            //ai::printMatrix(la);
            //ai::printMarker();
            bool srt = 1;
            if(la[1][0] < 0. ){
                std::cout<<"WRONG LAYER SORT!!!!"<<std::endl;
                srt = 0;
            }
        //ai::printMarker();
            if(srt == 1){
                for(size_t i = la.size() ; i > 1 ; ){
                    --i;
                    layers.push_back(std::vector<double>{la[i][0], la[i][1], la[i][2],la[i][3],la[i][4],la[i][5]} );
                    //std::cout<<"iter = "<<i<<std::endl;
                }
            }
            if(srt == 0){
                for(size_t i = 1; i < la.size() ;++i ){
                    layers.push_back(std::vector<double>{la[i][0], la[i][1], la[i][2],la[i][3],la[i][4],la[i][5]} );
                    std::cout<<"iter = "<<i<<std::endl;
                }
            }
            la.clear();
            ai::printMarker();
            std::cout<<"layers:"<<std::endl;
            ai::printMatrix(layers);
    }
void ParseGeomechanics(std::string filename, std::vector<std::vector<double> >&layers){
    std::vector<std::vector<double> > la;
    //ai::printMarker();
    std::cout<<"la:"<<std::endl;
    ai::parseFileInMatrix(filename,' ',la);
    int productiv = la.size();
    productiv/=2;
    double h = la[productiv][0] + la[productiv][1]/2.;

    std::cout<<"productive layer num = "<<productiv<<std::endl;
    std::cout<<"h = "<<h<<std::endl;
    for(size_t i = 0; i < la.size(); ++i){
        la[i][0] -= h;
        la[i][1] += la[i][0];


    }
    //ai::printMatrix(la);
    for(size_t i = 0; i<la.size(); ++i){//   start      finish    stress   Young      poisson   Carter
        layers.push_back(std::vector<double>{la[i][0], la[i][1], la[i][5], la[i][3], la[i][4], la[i][6] });
    }
 la.clear();
 ai::printMatrix(layers);
}

int main(){

std::vector<std::vector<double> > layers;
std::string filename = "simpleLayers2.txt";

 // ParseFromLayers(filename, layers);
 ParseGeomechanics(filename, layers);

 ApproximateLayers(layers);

    return 1 ;
}

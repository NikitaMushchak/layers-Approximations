#include <iostream>
#include <vector>
#include "ai.hh"



bool desc (double i,double j) { return (i<j); }



int main(){


    std::string filename = "layers2.txt";

    std::vector<std::vector<double> > layers;

    ai::parseFileInMatrix(filename,' ',layers);

    //ai::printMatrix(layers);

    std::size_t num;
    double H ;
    double dx;

    for(std::size_t i = 1 ; i < layers.size(); ++i){
            if(layers[i][0]*layers[i][1] < 0.){
                H = layers[i][1] - layers[i][0];
                num = i;
            }

    }
    double sigma = layers[num][2];
    //std::cout<<"sigma = "<<sigma<<std::endl;
    std::vector<std::size_t> sloi;

    double h;
    double hnew;
    if(H < 11.2){

        //H = 11.2;

        std::size_t it;
        it = num ;

        sloi.push_back(it);

        h = layers[it][1] - layers[it][0];
        it--;

        while(  11.2/2. > std::abs(layers[it][1])){
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
        dx = 1;
    }


    std::size_t min = ai::min(sloi);
    std::size_t max = ai::max(sloi);

    std::vector<std::vector<double> > layers1;

    for(size_t i = 1; i < layers.size(); i++ ){
        layers1.push_back(std::vector<double>{layers[i][0] ,layers[i][1], layers[i][2], layers[i][3], layers[i][4], layers[i][5]} );
        if(i == min){
            layers1.push_back(std::vector<double>{-H/2. ,H/2., sigma, layers[i][3], layers[i][4], layers[i][5]});
            i = max-1;
        }

    }
    // Обрабатываем пересечение слоев

    // for(std::size_t i = 0 ; i < layers1.size()-1;++i){
    //         if( abs(layers1[i][1])>abs(layers1[i+1][0]) ){
    //
    //             layers1[i+1][0] = layers1[i][1];
    //         }
    //         else{
    //             layers1[i][1] = layers1[i+1][0];
    //         }
    // }
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


    for(size_t i = num-1; i > 0; --i){  //Идем по слоям layers1 i - итератор
        //std::cout<<" "<<std::endl;
        //std::cout<<"iteration = "<<i<<std::endl;
        h = layers1[i][1]-layers1[i][0];     //j - итератор по mesh1
        //std::cout<<"h curent layer = "<<h<<std::endl;
        dl = std::ceil(h / dx);
        //std::cout<<"number of dx in layer = "<<dl<<std::endl;
        if(dl < eps){
            mesh1[j][2]+=h*layers1[i][2];
            //std::cout<<"mesh1 = "<<mesh1[j][2]<<std::endl;
            i--;
            mesh1[j][2]+=(dx-h)*layers1[i][2];
            //std::cout<<"mesh1 = "<<mesh1[j][2]<<std::endl;
            i++;
            j--;
        }

        if(dl > eps){
            int iter = 0;
            //std::cout<<"j = "<<j<<std::endl;
            //std::cout<<"i = "<<i<<std::endl;

            h = 0;
            while(iter < (int)dl){
                mesh1[j][2] = layers1[i][2];
                j--;
                iter++;
                //std::cout<<"j = "<<j<<"   iter = "<<iter<<std::endl;
            }
            //j++;
            //i++;
        }
    }

     j = middle +1;


    for(size_t i = num+1; i < layers1.size(); ++i){  //Идем по слоям layers1 i - итератор
        std::cout<<" "<<std::endl;
        std::cout<<"iteration = "<<i<<std::endl;
        h = layers1[i][1]-layers1[i][0];     //j - итератор по mesh1
        std::cout<<"h curent layer = "<<h<<std::endl;
        dl = std::ceil(h / dx);
        std::cout<<"number of dx in layer = "<<dl<<std::endl;
        if(dl < eps){
            mesh1[j][2]+=h*layers1[i][2];
            //std::cout<<"mesh1 = "<<mesh1[j][2]<<std::endl;
            i++;
            mesh1[j][2]+=(dx-h)*layers1[i][2];
            //std::cout<<"mesh1 = "<<mesh1[j][2]<<std::endl;
            i--;
            j++;
        }

        if(dl > eps){
            int iter = 0;
            std::cout<<"j = "<<j<<std::endl;
            std::cout<<"i = "<<i<<std::endl;

            h = 0;
            while(iter < (int)dl){
                if(j<mesh1.size())
                mesh1[j][2] = layers1[i][2];
                j++;
                iter++;
                std::cout<<"j = "<<j<<"   iter = "<<iter<<std::endl;
            }
            //j++;
            //i++;
        }
    }
    // double dist;
    // size_t j = num-1;
    // size_t it;
    //
    // for(size_t i = middle-1 ; i > 0 ; i--){ //идеем по mesh i
    //     rastlayer = 0;
    //     rastmesh = mesh1[i][1] - mesh1[i][0];
    //     std::cout<<"                                   "<<std::endl;
    //     std::cout<<" "<<std::endl;
    //     dist = 0;
    //     // for(size_t j = num-1 ; j > 0; j-- ){ //идем по layers1 j
    //     while(j > 0){      //идем по layers1 j- индекс
    //             //j--;
    //
    //             it=0;
    //             bool q = 1;
    //
    //             rastlayer += layers1[j][1] - layers1[j][0];
    //             dist = (dx - rastlayer) > 0 ? dx-rastlayer :dist;
    //             std::cout<<"j = "<<j<<std::endl;
    //             std::cout<<"dist = "<<dist<<std::endl;
    //             if(rastlayer < rastmesh){
    //                 std::cout<<"rastlayer < rastmesh"<<std::endl;
    //                 std::cout<<"rastlayer = "<<rastlayer<<"   rastmesh = "<<rastmesh<<std::endl;
    //                 mesh1[i][2] += layers1[j][2]*(layers1[j][1] - layers1[j][0]);
    //                 std::cout<<"mesh1["<<i<<"][2] = "<<mesh1[i][2]<<std::endl;
    //
    //             }
    //             if(rastlayer >= rastmesh){
    //                 std::cout<<"rastlayer >= rastmesh"<<std::endl;
    //                 std::cout<<"rastlayer = "<<rastlayer<<std::endl;
    //                 std::cout<<"layers1[j][1]-layers1[j][0] = "<<layers1[j][1]-layers1[j][0]<<std::endl;
    //                 dist = dist < eps ? dx : dist;
    //                 mesh1[i][2] += layers1[j][2]*(dist);
    //                 it=1;
    //                 mesh1[i][2]/=dx;
    //                 q = 0;
    //
    //                 std::cout<<"j = "<<j<<std::endl;
    //                 j--;
    //
    //                 std::cout<<"  # j  = "<<j<<std::endl;
    //                 std::cout<<"mesh1["<<i<<"][2] = "<<mesh1[i][2]<<std::endl;
    //                 break;
    //             }
    //             // if(!q){
    //             //     j--;
    //             //     break;
    //             //
    //             // }
    //             j--;
    //     }
    // }
    // for(size_t i = num-1 ; i > num -4 ; i--){ //идеем по layers1 i
    //     double rastlayer = layers1[i][1] - layers1[i][0];
    //     double rastmesh ;
    //     for(size_t j = middle-1 ; j > 0; j-- ){ //идем по mesh1 j
    //
    //             rastmesh += mesh1[j][1] - mesh1[j][0];
    //             if( rastlayer <= rastmesh ){
    //
    //                 mesh1[j][2] = layers1[i][2];
    //                 j--;
    //                 break;
    //             }
    //
    //
    //     }
    // }

    std::cout<<"mesh1:"<<std::endl;
    ai::printMatrix(mesh1);
    //num- номер стоки с продуктивным слоем
    //Запоминаем максимальную и минимальную высоту относительно нуля

   //
   //  double up = layers1[0][0];
   //  double down = layers1[layers1.size()-1][1];
   //
   //
   //  std::cout<<"UP = "<<up<<"    DOWN = "<<down<<std::endl;
   //
   //
   //  //делаем сетку
   //  std::vector<double> mesh;
   //  double coord = dx/2. ;
   //
   //  mesh.push_back(coord);
   //
   //  while(coord <  down){
   //
   //      coord += dx;
   //
   //      mesh.push_back(coord);
   //
   //  }
   //
   //  coord = -dx/2.;
   //
   //  mesh.push_back(coord);
   //  while(std::abs(coord)< std::abs(up)){
   //      coord -= dx;
   //
   //      mesh.push_back(coord);
   //
   //  }
   //  //Сортируем сетку
   //  std::sort(mesh.begin() ,mesh.end(), desc);
   //  //Ячейки
   //  std::vector<std::vector<double> >cells;
   //
   //  std::cout<<"Mesh:"<<std::endl;
   //  //ai::printVector(mesh);
   //  for(size_t i = 0 ; i < mesh.size() - 1 ;i++){
   //      //ai::printMarker();
   //      cells.push_back(std::vector<double>{mesh[i], mesh[i+1]} );
   //  }
   //  //mesh.clear();
   //  //ai::printMarker();
   //  std::cout<<"Cells"<<std::endl;
   //  ai::printMatrix(cells);
   // // матрица хранящая ячейки
   //  std::vector<std::vector<size_t> > sl;//для хранения слоев
   //
   //  for(size_t i = 0; i < layers1.size(); i++){
   //      for(size_t j = 0; j < cells.size(); ++j){
   //          if(layers1[i][0]<=cells[j][0] && layers1[i][1]>=cells[j][1])
   //          sl.push_back(std::vector<size_t>{i,j});
   //      }
   //  }
    // ai::printVector(sloi);
    std::cout<<"New layers"<<std::endl;
    ai::printMatrix(layers1);

    // std::cout<<"sl"<<std::endl;
    // ai::printMatrix(sl);
    return 1 ;
}

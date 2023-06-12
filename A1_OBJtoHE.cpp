#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include "HEMesh.h"

using namespace std;

typedef Vec<double, 3> Point;
typedef typename std::list<Edge<double>>::iterator EdgeIterator;
typedef typename std::list<Vert<double>>::iterator VertIterator;
typedef typename std::list<Face<double>>::iterator FaceIterator;

int main(){
    FILE *fp = fopen("test.txt","r");
    double x , y , z;
    int v1 , v2 , v3;
    char typ;
    list<Vert<double>> my_vert_list;
    my_vert_list.clear();
    list<Edge<double>> my_edge_list;
    my_edge_list.clear();
    list<Face<double>> my_face_list;
    my_face_list.clear();
    vector<T_IDX> vect_index;
    vect_index.clear();

    while(fscanf(fp,"%c%lf%lf%lf",&typ,&x,&y,&z)>0 ){
        if(typ=='\n'){continue;}
        //cout <<"typ:"<< typ<<"  x:" << x <<"    y:" << y <<"    z:"<< z<<endl;
        if(typ=='v'){
            my_vert_list.push_back(Vert<double>(Point(x, y, z)));
        }
        else {
            v1 = int(x-1);
            v2 = int(y-1);
            v3 = int(z-1);
            vect_index.push_back(v1);vect_index.push_back(v2);vect_index.push_back(v3);

            my_edge_list.push_back(Edge<double>());
            my_edge_list.push_back(Edge<double>());
            my_edge_list.push_back(Edge<double>());
            my_face_list.push_back(Face<double>());
        }
    }
    T_IDX num_verts = my_vert_list.size();
    T_IDX num_faces = my_face_list.size();
    T_IDX num_edges = num_faces * 3;
    T_IDX i_tmp;

    EdgeIterator edge_iterators[num_edges];
    i_tmp = 0;
    for (EdgeIterator edge_it_tmp = my_edge_list.begin(); edge_it_tmp != my_edge_list.end(); edge_it_tmp++, i_tmp++) {
        edge_iterators[i_tmp] = edge_it_tmp;
    }

    i_tmp = 0;
    VertIterator vert_iterators[num_verts];
    for (VertIterator vert_it = my_vert_list.begin(); vert_it != my_vert_list.end(); vert_it++, i_tmp++) {
        vert_iterators[i_tmp] = vert_it;
    }

    FaceIterator face_iterators[num_faces];
    i_tmp = 0;
    for (FaceIterator face_it = my_face_list.begin(); face_it != my_face_list.end(); face_it++, i_tmp++) {
        face_iterators[i_tmp] = face_it;
    }

    /// @note finding the other half of the edge
    map<T_IDX,vector<pair<T_IDX,T_IDX>>> HE_pair;
    //the key is index of the vex being pointed to, 
    //the value is pairs of vexes that point to it and the index of the edge
    HE_pair.clear();
    for(T_IDX i = 0 ; i<num_edges ; i++){
        T_IDX v_index = vect_index[i];
        T_IDX next_e_index;
        T_IDX next_v_index;
        T_IDX pair_e_index = 0;
        next_e_index= ( (i+1)%3 ? i+1 : i-2 );
        next_v_index = ((i+1)%3 ? vect_index[i+1] : vect_index[i-2]);
        int flag = -1;
        auto map_iter = HE_pair.find(v_index);
        if(map_iter != HE_pair.end()){
            flag = 0;
            for(int j = 0 ; j<map_iter->second.size();j++){
                if(next_v_index == map_iter->second[j].first){
                    //the other half exist
                    pair_e_index = map_iter->second[j].second;
                    flag = 1 ; break;
                }
            }
        }
        auto map_iter_next = HE_pair.find(next_v_index);
        if(map_iter_next != HE_pair.end())
            {HE_pair[next_v_index].push_back(make_pair(v_index,i));}//key exists but no right value
        else{
            //set up new pair with next_v_index as key
            flag = -1;
            vector<pair<T_IDX,T_IDX>> tmp_pair;
            tmp_pair.push_back(make_pair(v_index,i));
            HE_pair.insert(make_pair(next_v_index,tmp_pair));
            tmp_pair.clear();
        }
        edge_iterators[i] ->setIncidencyAndAdjacency(
            vert_iterators[v_index], face_iterators[i/3], edge_iterators[pair_e_index], edge_iterators[next_e_index]);
        edge_iterators[pair_e_index] ->setPairEdge(edge_iterators[i]);
        cout << i <<"   " <<flag<< ": v_index:" <<v_index<<"  face_index:"<<i/3<<"   pair_e_index:"<<pair_e_index<<"  next_e_index:"<<next_e_index<<endl;
    }
    cout<<"face:"<<endl;
    for(T_IDX i = 0 ; i<num_faces ; i++){
        face_iterators[i]->m_edge_it = edge_iterators[i*3];
        cout<<i<<": "<<i*3<<endl;
    }
    cout<<"vert:"<<endl;
    for(T_IDX i = 0 ; i<num_verts ; i++){
        T_IDX first_e_index;
        first_e_index = HE_pair[i][0].second;
        vert_iterators[i]->m_edge_it = edge_iterators[first_e_index];
        cout<<i<<": "<<first_e_index<<endl;
    }
    for (FaceIterator face_it2 = my_face_list.begin(); face_it2 != my_face_list.end(); face_it2++) {
        face_it2->m_edge_it->printEdge();
    }
    Mesh<double> MyMesh = Mesh<double>(my_vert_list,my_edge_list,my_face_list);
    MyMesh.exportOBJ("A1_out.txt");
    system("pause");
    return 0;
}
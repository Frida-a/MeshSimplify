#ifndef __HEMESH_H__
#define __HEMESH_H__

#include "matrix.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <cassert>
#include <fstream>
#include <list>


template<typename DType> struct Vert;
template<typename DType> struct Edge;
template<typename DType> struct Face;
template<typename DType> class  Mesh;


template<typename DType>
struct Vert {
    typedef Vec<DType, 3> Point;
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;
    
    T_IDX m_id = -1;
    Point m_coords;
    EdgeIterator m_edge_it; // incident edge

    Vert(Point coords): m_coords(coords) {}
};

template<typename DType>
struct Edge {
    typedef Vec<DType, 3> Point;
    typedef typename std::list<Vert<DType>>::iterator VertIterator;
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;
    typedef typename std::list<Face<DType>>::iterator FaceIterator;

    EdgeIterator m_pair_it;     // pair half-edge
    EdgeIterator m_next_it;     // next half-edge
    VertIterator m_vert_it;     // incident vert (which is the origin of this)
    FaceIterator m_face_it;     // incident face on the left side
    
    Edge() {}
    inline void setIncidencyAndAdjacency(VertIterator vert_it, FaceIterator face_it, EdgeIterator pair_it, EdgeIterator next_it) {
        m_vert_it = vert_it;    m_face_it = face_it;    m_pair_it = pair_it;    m_next_it = next_it;
    }
    inline void setPairEdge(EdgeIterator pair_it){
        m_pair_it = pair_it;
    }
    inline void setFace(FaceIterator face_it){
        m_face_it = face_it;
    }
    inline void setNext(EdgeIterator next_it){
        m_next_it = next_it;
    }
    inline void printEdge(){std::cout<<"id of Edge:" << (*m_vert_it).m_id<<" "<<std::endl;}
};


template<typename DType>
struct Face {
    typedef Vec<DType, 2> Point;
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;
    
    EdgeIterator m_edge_it;

    Face() {}
};

template<typename DType>
class Mesh {
    typedef Vec<DType, 3> Point;
    typedef typename std::list<Vert<DType>>::iterator VertIterator;
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;
    typedef typename std::list<Face<DType>>::iterator FaceIterator;

    std::list<Vert<DType>> vert_list;
    std::list<Edge<DType>> edge_list;
    std::list<Face<DType>> face_list;
public:
    Mesh(std::list<Vert<DType>> t_vert_list,std::list<Edge<DType>> t_edge_list,std::list<Face<DType>> t_face_list) {
        vert_list = t_vert_list; edge_list = t_edge_list; face_list = t_face_list;
    }
    Mesh(){}
    void initTetrahedron();
    void initFromFile(const char* inputf);
    void singleEdgeCollaspe(EdgeIterator target_it);
    void linearlyIndexVerts(T_IDX index_start);
    void exportOBJ(std::string obj_fn);
};


template<typename DType>
void Mesh<DType>::initTetrahedron() {
    /// @brief initializes this mesh to a tetrahedron, corresponding to the following OBJ file:
    // ---------OBJ example---------
    // v 0 0 0
    // v 1 0 0
    // v 0 1 0
    // v 0 0 1
    // f 1 3 2
    // f 1 2 4
    // f 2 3 4
    // f 1 4 3
    // -----------------------------
    vert_list.clear();
    edge_list.clear();
    face_list.clear();


    T_IDX num_verts = 4;
    T_IDX num_faces = 4;
    T_IDX num_edges = 12;   // half edges
    T_IDX i_tmp = 0;

    /// @note create 4 vertices with coordinates
    vert_list.push_back(Vert<DType>(Point(0, 0, 0)));
    vert_list.push_back(Vert<DType>(Point(1, 0, 0)));
    vert_list.push_back(Vert<DType>(Point(0, 1, 0)));
    vert_list.push_back(Vert<DType>(Point(0, 0, 1)));
    // temporary records for random access
    i_tmp = 0;
    VertIterator vert_iterators[num_verts];
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++, i_tmp++) {
        vert_iterators[i_tmp] = vert_it;
    }

    /// @note create 4 faces, to be connected later
    for (T_IDX i = 0; i < num_faces; i++) {
        face_list.push_back(Face<DType>());
    }
    FaceIterator face_iterators[num_faces];
    i_tmp = 0;
    for (FaceIterator face_it = face_list.begin(); face_it != face_list.end(); face_it++, i_tmp++) {
        face_iterators[i_tmp] = face_it;
    }

    /// @note create 12 half edges, to be connected later
    for (T_IDX i = 0; i < num_edges; i++) {
        edge_list.push_back(Edge<DType>());
    }
    EdgeIterator edge_iterators[num_edges];
    i_tmp = 0;
    for (EdgeIterator edge_it_tmp = edge_list.begin(); edge_it_tmp != edge_list.end(); edge_it_tmp++, i_tmp++) {
        edge_iterators[i_tmp] = edge_it_tmp;
    }

    /// @note connect half-edges to its neighbors
    // format: setIncAndAdj(vert, face, pair, next)
    // face 0
    edge_iterators[0 ]->setIncidencyAndAdjacency(vert_iterators[0], face_iterators[0], edge_iterators[11], edge_iterators[1 ]); // h0: v0 -> v2
    edge_iterators[1 ]->setIncidencyAndAdjacency(vert_iterators[2], face_iterators[0], edge_iterators[6 ], edge_iterators[2 ]); // h1: v2 -> v1
    edge_iterators[2 ]->setIncidencyAndAdjacency(vert_iterators[1], face_iterators[0], edge_iterators[3 ], edge_iterators[0 ]); // h2: v1 -> v0
    // face 1
    edge_iterators[3 ]->setIncidencyAndAdjacency(vert_iterators[0], face_iterators[1], edge_iterators[2 ], edge_iterators[4 ]); // h3: v0 -> v1
    edge_iterators[4 ]->setIncidencyAndAdjacency(vert_iterators[1], face_iterators[1], edge_iterators[8 ], edge_iterators[5 ]); // h4: v1 -> v3
    edge_iterators[5 ]->setIncidencyAndAdjacency(vert_iterators[3], face_iterators[1], edge_iterators[9 ], edge_iterators[3 ]); // h5: v3 -> v0
    // face 2
    edge_iterators[6 ]->setIncidencyAndAdjacency(vert_iterators[1], face_iterators[2], edge_iterators[1 ], edge_iterators[7 ]); // h6: v1 -> v2
    edge_iterators[7 ]->setIncidencyAndAdjacency(vert_iterators[2], face_iterators[2], edge_iterators[10], edge_iterators[8 ]); // h7: v2 -> v3
    edge_iterators[8 ]->setIncidencyAndAdjacency(vert_iterators[3], face_iterators[2], edge_iterators[4 ], edge_iterators[6 ]); // h8: v3 -> v1
    // face 3
    edge_iterators[9 ]->setIncidencyAndAdjacency(vert_iterators[0], face_iterators[3], edge_iterators[5 ], edge_iterators[10]); // h9: v0 -> v3
    edge_iterators[10]->setIncidencyAndAdjacency(vert_iterators[3], face_iterators[3], edge_iterators[7 ], edge_iterators[11]); // h10: v3 -> v2
    edge_iterators[11]->setIncidencyAndAdjacency(vert_iterators[2], face_iterators[3], edge_iterators[0], edge_iterators[9 ]); // h11: v2 -> v0

    /// @note connect faces to incident half-edges
    face_iterators[0]->m_edge_it = edge_iterators[0];   // f0: v0 -> v2 -> v1
    face_iterators[1]->m_edge_it = edge_iterators[3];   // f1: v0 -> v1 -> v3
    face_iterators[2]->m_edge_it = edge_iterators[6];   // f2: v1 -> v2 -> v3
    face_iterators[3]->m_edge_it = edge_iterators[9];   // f3: v0 -> v3 -> v2

    /// @note connect verts to incident half-edges
    vert_iterators[0]->m_edge_it = edge_iterators[2];
    vert_iterators[1]->m_edge_it = edge_iterators[3];
    vert_iterators[2]->m_edge_it = edge_iterators[6];
    vert_iterators[3]->m_edge_it = edge_iterators[9];
}

template<typename DType>
void Mesh<DType>::linearlyIndexVerts(T_IDX index_start) {
    T_IDX v_counter = index_start;
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++, v_counter++) {
        vert_it->m_id = v_counter;
    }
}

template<typename DType>
void Mesh<DType>::initFromFile(const char* inputf){
    FILE *fp = fopen(inputf,"r");
    double x , y , z;
    int v1 , v2 , v3;
    char typ;

    vert_list.clear();
    edge_list.clear();
    face_list.clear();
    //vector index from faces in order
    std::vector<T_IDX> vect_index;
    vect_index.clear();
    while(fscanf(fp,"%c%lf%lf%lf",&typ,&x,&y,&z)>0 ){
        if(typ=='\n'){continue;}
        //cout <<"typ:"<< typ<<"  x:" << x <<"    y:" << y <<"    z:"<< z<<endl;
        if(typ=='v'){
            vert_list.push_back(Vert<double>(Point(x, y, z)));
        }
        else {
            v1 = int(x-1);
            v2 = int(y-1);
            v3 = int(z-1);
            vect_index.push_back(v1);vect_index.push_back(v2);vect_index.push_back(v3);

            edge_list.push_back(Edge<double>());
            edge_list.push_back(Edge<double>());
            edge_list.push_back(Edge<double>());
            face_list.push_back(Face<double>());
        }
    }
    T_IDX num_verts = vert_list.size();
    T_IDX num_faces = face_list.size();
    T_IDX num_edges = num_faces * 3;
    T_IDX i_tmp;

    EdgeIterator edge_iterators[num_edges];
    i_tmp = 0;
    for (EdgeIterator edge_it_tmp = edge_list.begin(); edge_it_tmp != edge_list.end(); edge_it_tmp++, i_tmp++) {
        edge_iterators[i_tmp] = edge_it_tmp;
    }

    i_tmp = 0;
    VertIterator vert_iterators[num_verts];
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++, i_tmp++) {
        vert_iterators[i_tmp] = vert_it;
    }

    FaceIterator face_iterators[num_faces];
    i_tmp = 0;
    for (FaceIterator face_it = face_list.begin(); face_it != face_list.end(); face_it++, i_tmp++) {
        face_iterators[i_tmp] = face_it;
    }

    /// @note finding the other half of the edge
    std::map<T_IDX,std::vector<std::pair<T_IDX,T_IDX>>> HE_pair;
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
            {HE_pair[next_v_index].push_back(std::make_pair(v_index,i));}//key exists but no right value
        else{
            //set up new pair with next_v_index as key
            flag = -1;
            std::vector<std::pair<T_IDX,T_IDX>> tmp_pair;
            tmp_pair.push_back(std::make_pair(v_index,i));
            HE_pair.insert(std::make_pair(next_v_index,tmp_pair));
            tmp_pair.clear();
        }
        edge_iterators[i] ->setIncidencyAndAdjacency(
            vert_iterators[v_index], face_iterators[i/3], edge_iterators[pair_e_index], edge_iterators[next_e_index]);
        edge_iterators[pair_e_index] ->setPairEdge(edge_iterators[i]);
        //std::cout << i <<"   " <<flag<< ": v_index:" <<v_index<<"  face_index:"<<i/3<<"   pair_e_index:"<<pair_e_index<<"  next_e_index:"<<next_e_index<<std::endl;
    }
    //std::cout<<"face:"<<std::endl;
    for(T_IDX i = 0 ; i<num_faces ; i++){
        face_iterators[i]->m_edge_it = edge_iterators[i*3];
        //std::cout<<i<<": "<<i*3<<std::endl;
    }
    //std::cout<<"vert:"<<std::endl;
    for(T_IDX i = 0 ; i<num_verts ; i++){
        T_IDX first_e_index;
        first_e_index = HE_pair[i][0].second;
        vert_iterators[i]->m_edge_it = edge_iterators[first_e_index];
        //std::cout<<i<<": "<<first_e_index<<std::endl;
    }
    //for (FaceIterator face_it2 = face_list.begin(); face_it2 != face_list.end(); face_it2++) {
        //face_it2->m_edge_it->printEdge();}

}

template<typename DType>
void Mesh<DType>::singleEdgeCollaspe(EdgeIterator target_it){
    VertIterator vert_it1 = target_it->m_vert_it;
    VertIterator vert_it2 = target_it->m_pair_it->m_vert_it;
    Point NewPt;
    for(int i = 0 ; i < 3 ; i++){
        NewPt.coords[i]= (vert_it1->m_coords[i] + vert_it2->m_coords[i]);
    }
    //add new vert,edges,faces
    Vert<DType> Newvert = Vert(NewPt);
    VertIterator vert_nearby;
    EdgeIterator edge_in, edge_out,edge_in_before,edge_out_before;//in view of new vert
    FaceIterator face_left;
    EdgeIterator edge_it = vert_it1->m_edge_it;
    Point publicPt1,publicPt2;
    //first point
    do{
        if(edge_it == target_it ){
            publicPt1 = vert_nearby->m_coords;
            edge_it = edge_it->m_pair_it->m_next_it;
            publicPt2 = edge_it->m_pair_it->m_vert_it->m_coords;
            continue;}

        vert_nearby = edge_it->m_pair_it->m_vert_it;
        face_left->m_edge_it = edge_in_before;
        if(edge_it != vert_it1->m_edge_it){
            edge_out_before = edge_out;
            edge_out_before->setFace(face_left);
            edge_out_before->m_next_it->setFace(face_left);//bind surrounding edges with new faces
            edge_list.push_back((*edge_out_before));//stored after the connection is built
        }
        edge_in->setIncidencyAndAdjacency(vert_nearby,face_left,edge_out,edge_out);//next_it to be connected
        edge_out->setIncidencyAndAdjacency(Newvert,face_left,edge_in,edge_it->m_next_it);//face to be connected
        edge_in_before->setNext(edge_out);

        edge_list.push_back((*edge_in_before));
        face_list.push_back((*face_left));
        edge_in_before = edge_in;
        
        edge_it = edge_it->m_pair_it->m_next_it;
    }while (edge_it != vert_it1->m_edge_it);
    face_left->m_edge_it = edge_in;
    edge_out->setFace(face_left);
    edge_list.push_back(edge_out);//stored after the connection is built
    //second point
    edge_it = vert_it2->m_edge_it;
    do{
        vert_nearby = edge_it->m_pair_it->m_vert_it;
        if(edge_it == target_it->m_pair_it || vert_nearby->m_coords == publicPt1.coords
             ||vert_nearby->m_coords == publicPt2.coords){
            edge_it = edge_it->m_pair_it->m_next_it;
            continue;}

        face_left->m_edge_it = edge_in;
        if(edge_it != vert_it2->m_edge_it){
            edge_out->setFace(face_left);
            edge_out->m_next_it->setFace(face_left);//bind surrounding edges with new faces
            edge_list.push_back(edge_out);//stored after the connection is built
        }
        EdgeIterator edge_next = edge_out;
        edge_in->setIncidencyAndAdjacency(vert_nearby,face_left,edge_out,edge_out);//next_it to be connected
        edge_out->setIncidencyAndAdjacency(Newvert,face_left,edge_in,edge_it->m_next_it);//face to be connected
        edge_in_before = edge_in;
        face_list.push_back(face_left);
        edge_list.push_back(edge_in);

        edge_it = edge_it->m_pair_it->m_next_it;
    }while (edge_it != vert_it2->m_edge_it);
    face_left->m_edge_it = edge_in;
    edge_out->setFace(face_left);
    edge_list.push_back(edge_out);//stored after the connection is built

    vert_list.push_back(Newvert);

    //delete points,edges,faces
        //first point
    edge_it = vert_it1->m_edge_it;
    do{
        face_list.erase(edge_it->m_face_it);
        EdgeIterator tmp_edge = edge_it;
        edge_it = edge_it->m_pair_it->m_next_it;
        edge_list.erase(tmp_edge->m_pair_it);
        edge_list.erase(tmp_edge);
    }while (edge_it != vert_it1->m_edge_it);
    //second point
    edge_it = vert_it2->m_edge_it;
    do{
        if(edge_it!=target_it->m_next_it){
            face_list.erase(edge_it->m_face_it);//otherwise the face has been erased
        }
        EdgeIterator tmp_edge = edge_it;
        edge_it = edge_it->m_pair_it->m_next_it;
        edge_list.erase(tmp_edge->m_pair_it);
        edge_list.erase(tmp_edge);
    }while (edge_it != vert_it2->m_edge_it);

    vert_list.erase(vert_it1);
    vert_list.erase(vert_it2);
    
}

template<typename DType>
void Mesh<DType>::exportOBJ(std::string obj_fn) {

    // file stream
    std::ofstream fout(obj_fn);
    
    // assign ids to verts and output their coordinates
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++) {
        Point coords = vert_it->m_coords;
        fout << "v " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }

    // OBJ subscripts starts from 1
    linearlyIndexVerts(1);

    //for (FaceIterator face_it2 = face_list.begin(); face_it2 != face_list.end(); face_it2++) {
        //face_it2->m_edge_it->printEdge();}

    // output faces
    for (FaceIterator face_it = face_list.begin(); face_it != face_list.end(); face_it++) {
        T_IDX vid_0, vid_1, vid_2;

        EdgeIterator edge_it = face_it->m_edge_it;
        vid_0 = edge_it->m_vert_it->m_id;  edge_it = edge_it->m_next_it;
        vid_1 = edge_it->m_vert_it->m_id;  edge_it = edge_it->m_next_it;
        vid_2 = edge_it->m_vert_it->m_id;

        fout << "f " << vid_0 << " " << vid_1 << " " << vid_2 << "\n";
    }
}




#endif
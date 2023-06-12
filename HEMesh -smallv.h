#ifndef __HEMESH_H__
#define __HEMESH_H__

#include "matrix.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <unordered_map>
#include <cstddef>
#include <functional>
#include <algorithm>
#include <vector>

#include <cassert>
#include <fstream>
#include <list>
#include <iterator>


template<typename DType> struct Vert;
template<typename DType> struct Edge;
template<typename DType> struct Face;
template<typename DType> class  Mesh;

#define INF 0x3f3f3f3f;

template<typename DType>
struct Vert {
    typedef Vec<DType, 3> Point;
    typedef Mat<DType,1,4> PlaneParam;
    typedef Mat<DType,4,4> VertQ;
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;

    T_IDX m_id = -1;
    Point m_coords;
    EdgeIterator m_edge_it; // incident edge
    bool val = 1; //0:invalid
    VertQ m_vq;

    Vert(Point coords): m_coords(coords) {}
    void inval(){val = 0;}
    bool isval(){   return val;}
    Point get_coords(){ return m_coords;}
    void calcuVq();
};

template<typename DType>
struct Edge {
    typedef Vec<DType, 3> Point;
    typedef Mat<DType,4,4> VertQ;
    typedef typename std::list<Vert<DType>>::iterator VertIterator;
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;
    typedef typename std::list<Face<DType>>::iterator FaceIterator;

    EdgeIterator m_pair_it;     // pair half-edge
    EdgeIterator m_next_it;     // next half-edge
    VertIterator m_vert_it;     // incident vert (which is the origin of this)
    FaceIterator m_face_it;     // incident face on the left side
    bool val = 1; //0:invalid
    Point best_pt;
    DType cost;

    Edge() {}
    inline void setIncidencyAndAdjacency(VertIterator vert_it, FaceIterator face_it, EdgeIterator pair_it, EdgeIterator next_it) {
        m_vert_it = vert_it;    m_face_it = face_it;    m_pair_it = pair_it;    m_next_it = next_it;
    }
    inline void setIncidency(VertIterator vert_it, FaceIterator face_it, EdgeIterator pair_it) {
        m_vert_it = vert_it;    m_face_it = face_it;    m_pair_it = pair_it;
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
    void inval(){val = 0;}
    bool isval(){   return val;}
    void calcu_bestpt();
};


template<typename DType>
struct Face {
    typedef Vec<DType, 3> Point;
    typedef Mat<DType,1,4> PlaneParam;
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;
    
    EdgeIterator m_edge_it;
    bool val = 1; //0:invalid
    PlaneParam m_pp;
    Point p_norm;//p_norm : the vertical vector

    void inval(){val = 0;}
    bool isval(){   return val;}
    void calcu_Param();
    Face() {}
};
template<typename DType>
struct Edge_cost{
    typedef typename std::list<Edge<DType>>::iterator EdgeIterator;
    EdgeIterator m_edge_half;
    DType m_cost;

    Edge_cost() {}
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
    std::vector<Edge_cost<DType>> m_edge_cost;
public:
    Mesh(std::list<Vert<DType>> t_vert_list,std::list<Edge<DType>> t_edge_list,std::list<Face<DType>> t_face_list) {
        vert_list = t_vert_list; edge_list = t_edge_list; face_list = t_face_list;
    }
    Mesh(){}
    void initTetrahedron();
    void initFromFile(const char* inputf);
    void singleEdgeCollaspe(EdgeIterator target_it);
    void insertion(FaceIterator target_face);
    void EdgeFlip(EdgeIterator target_edge);
    inline EdgeIterator searchFirstEdge(){  return edge_list.begin();}
    inline FaceIterator searchFirstFace(){  return face_list.begin();}
    void preCalcu();
    int EdgeCollaspe(EdgeIterator target_it , Point NewPt);
    void meshSimplfy(int target_face_num);
    void intialHeapContainer();
    void printInfo();
    void linearlyIndexVerts(T_IDX index_start);
    void exportOBJ(std::string obj_fn);
};


struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

template<typename DType>
bool Edge_cost_cmp(const Edge_cost<DType> a, const Edge_cost<DType> b){
    return a.m_cost > b.m_cost;
}



template<typename DType>
void Face<DType>:: calcu_Param(){
    EdgeIterator edge_it;

    Point p1 = this->m_edge_it->m_vert_it->m_coords;
    edge_it = m_edge_it->m_next_it;
    Point p2 = edge_it->m_vert_it->m_coords;
    edge_it = edge_it->m_next_it;
    Point p3 = edge_it->m_vert_it->m_coords;

    Point ptmp = p2-p1;
    p_norm = ptmp.cross(p3-p1);
    p_norm = p_norm.normalized();
    DType d = - p_norm.dot(p1);

    m_pp(0,0) = p_norm[0];
    m_pp(0,1) = p_norm[1];
    m_pp(0,2) = p_norm[2];
    m_pp(0,3) = d;

    //std::cout<<m_pp<<std::endl;
}

template<typename DType>
void Vert<DType>:: calcuVq(){
    for(int i = 0 ; i<4 ; i++){
        for(int j = 0; j<4 ; j++){
            m_vq (i , j) = 0;
        }
    }
    EdgeIterator edge_it;
    VertQ tmpVQ;
    PlaneParam tmp_pp;
    edge_it = m_edge_it;
    do{
        tmp_pp = edge_it->m_face_it->m_pp;
        m_vq += tmp_pp.transpose().matmul(tmp_pp);
        //std::cout<<m_vq<<std::endl;
        edge_it = edge_it->m_next_it->m_pair_it;
    }while(edge_it != m_edge_it);
    //std::cout<<m_vq<<std::endl;
}

template<typename DType>
void Edge<DType>:: calcu_bestpt(){
    VertQ q_sum;
    VertIterator vert_it = m_vert_it;
    q_sum = vert_it->m_vq + m_pair_it->m_vert_it->m_vq;
    Mat<DType,3,3> mat_a,a_invert, a_cofact;
    Mat<DType,3,1> mat_b;
    
    mat_a(0,0) = q_sum(0,0);    mat_a(0,1) = q_sum(0,1);    mat_a(0,2) = q_sum(0,2);
    mat_a(1,0) = q_sum(0,1);    mat_a(1,1) = q_sum(1,1);    mat_a(1,2) = q_sum(1,2);
    mat_a(2,0) = q_sum(0,2);    mat_a(2,1) = q_sum(1,2);    mat_a(2,2) = q_sum(2,2);

    mat_b(0,0) = q_sum(0,3);    mat_b(1,0) = q_sum(1,3);    mat_b(2,0) = q_sum(2,3);

    //std::cout<<mat_a<<std::endl;
    DType del_a = mat_a.det();
    if(del_a != 0){
        for(int i = 0 ; i< 3 ; i++){
            for(int j = 0 ; j< 3 ;j++){
                Mat<DType,2,2> minor;
                int a =0;   int b = 0;
                for(int m = 0 ;m<3 ; m++){
                    if(m != i){
                        for(int n = 0 ; n<3 ; n++){
                            if(n != j){
                            //int a =0;   int b = 0;
                                minor(a,b) = mat_a(m,n);
                                b++;
                            }
                        }
                        b = 0;
                        a++;
                    }

                }
                a_cofact(i,j) = ((i+j)%2? -1*minor.det() : minor.det());
            }
        }

        a_invert = a_cofact.transpose()/del_a;
        //std::cout<<"Invert test:    "<<a_invert<<std::endl;
        auto v = a_invert.matmul(mat_b);
        //std::cout << v <<std::endl;
        best_pt[0] = -v(0,0);  best_pt[1] = -v(1,0);  best_pt[2] = -v(2,0);
        Mat<DType,4,1> v_hat;
        v_hat(0,0) = best_pt[0];    v_hat(1,0) = best_pt[1];    v_hat(2,0) = best_pt[2];
        v_hat(3,0) = 1;
        auto res = v_hat.transpose().matmul(q_sum).matmul(v_hat);
        cost = res(0,0);
        //std::cout<<cost<<std::endl;
    }else{//discard if uninvertible
        cost = INF;
        std::cout<<"uninvertible  "<<cost<<std::endl;
    }
}

template<typename DType>
void Mesh<DType>::preCalcu(){
    for (FaceIterator face_it = face_list.begin(); face_it != face_list.end(); face_it++) {
        if(face_it->isval()){
            //Point coords = face_it->m_edge_it->m_vert_it->m_coords;
            //std::cout << "f " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
            face_it->calcu_Param();
        }
    }
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++) {
        if(vert_it->isval()){
            vert_it->calcuVq();
        }
    }
    for (EdgeIterator edge_it = edge_list.begin(); edge_it != edge_list.end(); edge_it++) {

        edge_it->calcu_bestpt();
    }
}

template<typename DType>
void Mesh<DType>::intialHeapContainer(){
    m_edge_cost.clear();
    EdgeIterator edge_it;
    for (FaceIterator face_it = face_list.begin(); face_it != face_list.end(); face_it++) {
        if(face_it->isval()){
            Edge_cost<DType> tmpEcs;
            edge_it = face_it->m_edge_it;
            if(edge_it->isval()){
                tmpEcs.m_cost = edge_it->cost;
                tmpEcs.m_edge_half = edge_it;
                m_edge_cost.push_back(tmpEcs);
            }

            edge_it = edge_it->m_next_it;
            if(edge_it->isval()){
                tmpEcs.m_cost = edge_it->cost;
                tmpEcs.m_edge_half = edge_it;
                m_edge_cost.push_back(tmpEcs);
            }

            edge_it = edge_it->m_next_it;
            if(edge_it->isval()){
                tmpEcs.m_cost = edge_it->cost;
                tmpEcs.m_edge_half = edge_it;
                m_edge_cost.push_back(tmpEcs);
            }
        }
    }
}

template<typename DType>
int Mesh<DType>::EdgeCollaspe(EdgeIterator target_it , Point NewPt){
    VertIterator vert_it1 = target_it->m_vert_it;
    VertIterator vert_it2 = target_it->m_pair_it->m_vert_it;
    //test
    VertIterator vert_temp;
    vert_temp = vert_it1;
    vert_it1 = vert_it2;
    vert_it2 = vert_temp;

    //add new vert,edges,faces
    Vert<DType> Newvert = Vert<DType>(NewPt);
    vert_list.push_back(Newvert);
    VertIterator Newvert_it ;
    Newvert_it = prev(vert_list.end());//add new vert, edge to be connected

    std::vector<EdgeIterator> tmpEdgeIts;
    std::vector<FaceIterator> tmpFaceIts, oldFaceIts;
    std::vector<VertIterator> tmpVertIts;
    EdgeIterator edge_it;
    VertIterator vert_nearby;
    FaceIterator face_it;

    tmpVertIts.clear();
    tmpFaceIts.clear();
    tmpEdgeIts.clear();
    oldFaceIts.clear();

    T_IDX i_tmp = 0;
    //first point
    edge_it = target_it->m_pair_it->m_next_it;
    do{
        //std::cout << edge_it->m_vert_it->m_id<<std::endl;
        //locate surrounding verts
        vert_nearby = edge_it->m_pair_it->m_vert_it;
        tmpVertIts.push_back(vert_nearby);
        //add new edges
        edge_list.push_back(Edge<DType>());//in
        edge_list.push_back(Edge<DType>());//out
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts[i_tmp*2] = prev(prev(edge_list.end()));
        tmpEdgeIts[i_tmp*2 + 1] = prev(edge_list.end());
        tmpEdgeIts[i_tmp*2 + 1] ->setNext(edge_it->m_next_it);//set next edge for out edges
        //add new faces
        face_list.push_back(Face<DType>());
        tmpFaceIts.push_back(face_it);
        tmpFaceIts[i_tmp] = prev(face_list.end());

        oldFaceIts.push_back(edge_it->m_face_it);

        i_tmp++;
        //std::cout<<"EC1";
        edge_it = edge_it->m_pair_it->m_next_it;

        if(i_tmp>50){ 
            //target_it->inval();  
            return 1;}

    }while (edge_it != target_it);   
    tmpEdgeIts[1] ->m_next_it = tmpEdgeIts[1] ->m_next_it ->m_pair_it ->m_next_it;
    oldFaceIts.push_back(target_it->m_face_it);
    //oldFaceIts.push_back(target_it->m_next_it->m_face_it);
    //second point
    edge_it = target_it->m_next_it;
    edge_it = edge_it->m_pair_it->m_next_it;
    do{
        //std::cout << edge_it->m_vert_it->m_id<<std::endl;
        //locate surrounding verts
        vert_nearby = edge_it->m_pair_it->m_vert_it;
        tmpVertIts.push_back(vert_nearby);
        //add new edges
        edge_list.push_back(Edge<DType>());//in
        edge_list.push_back(Edge<DType>());//out
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts[i_tmp*2] = prev(prev(edge_list.end()));
        tmpEdgeIts[i_tmp*2 + 1] = prev(edge_list.end());
        tmpEdgeIts[i_tmp*2 + 1] ->setNext(edge_it->m_next_it);
        //std::cout<<"set next:  "<<tmpEdgeIts[i_tmp*2+1]->m_next_it->m_vert_it->m_id <<std::endl;
        //add new faces
        face_list.push_back(Face<DType>());
        tmpFaceIts.push_back(face_it);
        tmpFaceIts[i_tmp] = prev(face_list.end());

        oldFaceIts.push_back(edge_it->m_face_it);

        i_tmp++;
        edge_it = edge_it->m_pair_it->m_next_it;
        //std::cout<<"ec2";
        if(i_tmp>100){  
            //target_it->inval(); 
            return 1;}
        }while (edge_it->m_pair_it->m_next_it != target_it->m_pair_it);//the last vert has already been discovered
    oldFaceIts.push_back(edge_it->m_face_it);

    //whether the new face will filp

    for(int i = 0 ; i < tmpVertIts.size() ; i++){
        Point p_norm_tmp,p_norm_tmp2;
        Point p1 = tmpVertIts[i]->m_coords;
        Point p2 = tmpVertIts[(i+1)%(tmpVertIts.size())]->m_coords;
        Point pnew = Newvert_it->m_coords;
        Point pold = vert_it1->m_coords;

        Point ptmp = p2-p1;
        p_norm_tmp = ptmp.cross(pnew-p1);
        p_norm_tmp = p_norm_tmp.normalized();

        p_norm_tmp2 = ptmp.cross(pold-p1);
        p_norm_tmp2 = p_norm_tmp2.normalized();

        if(p_norm_tmp.dot(p_norm_tmp2) < 0 ){
            //target_it->inval();
            std::cout<<"bad face";
            target_it->inval();
            for(int j = 0 ;j<tmpFaceIts.size() ; j++){
                face_list.pop_back();
                edge_list.pop_back();
                edge_list.pop_back();
            }
            vert_list.pop_back();
            return 1;
        }
    }
    for(int i = 0 ; i < tmpVertIts.size() ; i++){
        Point p_norm_tmp,p_norm_tmp2;
        Point p1 = tmpVertIts[i]->m_coords;
        Point p2 = tmpVertIts[(i+1)%(tmpVertIts.size())]->m_coords;
        Point pnew = Newvert_it->m_coords;
        Point pold = vert_it2->m_coords;//第二个顶点

        Point ptmp = p2-p1;
        p_norm_tmp = ptmp.cross(pnew-p1);
        p_norm_tmp = p_norm_tmp.normalized();

        p_norm_tmp2 = ptmp.cross(pold-p1);
        p_norm_tmp2 = p_norm_tmp2.normalized();

        if(p_norm_tmp.dot(p_norm_tmp2) < 0 ){
            //target_it->inval();
            std::cout<<"bad face";
            target_it->inval();
            for(int j = 0 ;j<tmpFaceIts.size() ; j++){
                face_list.pop_back();
                edge_list.pop_back();
                edge_list.pop_back();
            }
            vert_list.pop_back();
            return 1;
        }
    }



    //initalize edges and faces
    T_IDX tmp_vert_num = tmpVertIts.size();
    //std::cout << tmpEdgeIts.size()<<" "<<tmp_vert_num<<"    "<<tmpFaceIts.size()<<"    ";
    for(T_IDX tmp_i = 0 ; tmp_i<tmp_vert_num ; tmp_i++){
        tmpEdgeIts[tmp_i*2]->setIncidencyAndAdjacency(tmpVertIts[tmp_i],tmpFaceIts[tmp_i],tmpEdgeIts[tmp_i*2+1],tmpEdgeIts[(tmp_i*2+3+tmp_vert_num*2) % (tmp_vert_num*2)]);
        tmpEdgeIts[tmp_i*2+1]->setIncidency(Newvert_it,tmpFaceIts[(tmp_i+tmp_vert_num-1)%tmp_vert_num],tmpEdgeIts[tmp_i*2]);//next_it already set
        tmpFaceIts[tmp_i]->m_edge_it = tmpEdgeIts[tmp_i*2];
        //std::cout<<"tmp_i:  "<<tmp_i<<" ";
        tmpVertIts[tmp_i]->m_edge_it = tmpEdgeIts[tmp_i*2+1];
        tmpEdgeIts[tmp_i*2+1]->m_next_it->m_face_it = tmpFaceIts[(tmp_i+tmp_vert_num-1)%tmp_vert_num];//the faces & next_it for surrounding edges have changed
        tmpEdgeIts[tmp_i*2+1]->m_next_it->m_next_it = tmpEdgeIts[(tmp_i*2-2+tmp_vert_num*2)% (tmp_vert_num*2)];
    }
 
    Newvert_it->m_edge_it = tmpEdgeIts[0];

    //whether the new face fliped
    /*
    for(int i = 0 ;i<tmpFaceIts.size();i++){
        tmpFaceIts[i]->calcu_Param();
        if(tmpFaceIts[i]->p_norm.dot(oldFaceIts[0]->p_norm)<0 ){
            //std::cout<<tmpFaceIts[i]->p_norm[0]<<" "<<tmpFaceIts[i]->p_norm[1]<<" "<<tmpFaceIts[i]->p_norm[2]<<std::endl;
            //std::cout<<oldFaceIts[0]->p_norm[0]<<" "<<oldFaceIts[0]->p_norm[1]<<" "<<oldFaceIts[0]->p_norm[2]<<std::endl;
            target_it->inval();
            std::cout<<"here"<<std::endl;
            //delete newly added faces ,edges, verts
            Newvert_it->inval();

            for(int j = 0 ; j < tmpFaceIts.size();j++){
                tmpFaceIts[j]->inval();
            }
            return 1;
        }
    }
    */
    


    //delete points,edges,faces
    //std::cout<<"old face num   "<<oldFaceIts.size()<<std::endl;
    for(int i = 0 ; i<oldFaceIts.size(); i++){
        oldFaceIts[i]->inval();
    }

    target_it->m_pair_it->inval();
    target_it->inval();
    vert_it1->inval();
    vert_it2->inval();

    //maintain the heap
    
    for(int i = 0 ;i<tmpFaceIts.size();i++){
        tmpFaceIts[i]->calcu_Param();
    }
    
    Newvert_it->calcuVq();
    for(int i = 0 ; i<tmpVertIts.size();i++){
        tmpVertIts[i]->calcuVq();
    }
    edge_it = Newvert_it->m_edge_it->m_pair_it;
    do{
        edge_it->calcu_bestpt();
        edge_it->m_pair_it->calcu_bestpt();//and their pair !
        edge_it = edge_it->m_pair_it ->m_next_it;
    }while(edge_it != Newvert_it->m_edge_it->m_pair_it);
    for(int i = 0 ; i<tmpVertIts.size() ; i++){
        edge_it = tmpVertIts[i]->m_edge_it->m_pair_it;
        do{
            edge_it->calcu_bestpt();
            edge_it->m_pair_it->calcu_bestpt();//and their pair !
            edge_it = edge_it->m_pair_it ->m_next_it;
        }while(edge_it != tmpVertIts[i]->m_edge_it->m_pair_it);
    }
    //std::cout<<"NewPT " << NewPt[0]<<" "<<NewPt[1]<<" "<<NewPt[2]<<std::endl;
    return 0;
}

template<typename DType>
void Mesh<DType>::meshSimplfy(int target_face_num){
    preCalcu();

    //do edge collapse
    T_IDX face_num = face_list.size();
    Edge_cost<DType> tmpEc,tmpEc_before;
    int flag = 0;
    while(face_num > target_face_num){
        std::cout<<"face_num    "<<face_num<<std::endl;

        intialHeapContainer();
        make_heap(m_edge_cost.begin(),m_edge_cost.end(),Edge_cost_cmp<DType>);
        pop_heap(m_edge_cost.begin(),m_edge_cost.end(),Edge_cost_cmp<DType>);
        tmpEc_before = tmpEc;
        tmpEc = m_edge_cost.back();
        m_edge_cost.pop_back();
        if(tmpEc_before.m_edge_half == tmpEc.m_edge_half->m_pair_it && flag==0){  continue;}//skip pair edge
        //std::cout<<"min cost is :  "<<tmpEc.m_cost<<std::endl;
        flag = 0;

        EdgeIterator target_it = tmpEc.m_edge_half;
        VertIterator vert_it1 = target_it->m_vert_it;
        VertIterator vert_it2 = target_it->m_pair_it->m_vert_it;

        std::vector<VertIterator> tmpVertIts;
        std::vector<EdgeIterator> tmpEdgeIts;

        EdgeIterator edge_it;
        VertIterator vert_nearby;
        FaceIterator face_it;

        //check commen verts
        tmpVertIts.clear();
        T_IDX i_tmp = 0;
        T_IDX i_tmp2 = 0;

        //for first point
        edge_it = target_it->m_pair_it->m_next_it;
        do{
            //locate commen verts
            vert_nearby = edge_it->m_pair_it->m_vert_it;
            //add new edges
            //std::cout<<"first";
            tmpVertIts.push_back(vert_nearby);
            edge_it = edge_it->m_pair_it->m_next_it;
            i_tmp2++;
            if(i_tmp2>50){  target_it->inval();   break;}
        }while (edge_it != target_it);   
        //second point
        tmpEdgeIts.clear();
        edge_it = target_it->m_next_it;
        edge_it = edge_it->m_pair_it->m_next_it;
        do{
            //locate commen verts
            vert_nearby = edge_it->m_pair_it->m_vert_it;
            for(int i = 0  ; i<tmpVertIts.size() ; i++){
                if(vert_nearby == tmpVertIts[i]){
                    tmpEdgeIts.push_back(edge_it);
                }
            }
            //std::cout<<"second";
            i_tmp2++;
            if(i_tmp2>50){ target_it->inval(); break;}
            edge_it = edge_it->m_pair_it->m_next_it;
        }while (edge_it->m_pair_it->m_next_it != target_it->m_pair_it); 
        if(i_tmp2>50){ continue;}

        if(tmpEdgeIts.size()>0){ target_it->inval(); flag = 1;}
        //std::cout<<"flag:   "<<flag <<std::endl;
        if(flag == 0){
            flag = EdgeCollaspe(target_it,target_it->best_pt);
            if(flag == 0){
                face_num = face_num-2;//判断本次循环确实缩边了
            }
            i_tmp++;
            //std::cout<<i_tmp<<std::endl;
        }
    }
}

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
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++) {
        if(vert_it->isval()){
            vert_it->m_id = v_counter;
             v_counter++;
        }
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
            //std::cout<<"typ:"<<int(typ)<<"  x"<<x<<"    y"<<y<<"    z"<<z<<std::endl;
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
    //std::map<T_IDX,std::vector<std::pair<T_IDX,T_IDX>>> HE_pair;
    
    std::unordered_map<std::pair<T_IDX,T_IDX>,T_IDX,pair_hash> HE_pair;
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

        edge_iterators[i] ->setIncidencyAndAdjacency(
            vert_iterators[v_index], face_iterators[i/3], edge_iterators[pair_e_index], edge_iterators[next_e_index]);

        auto map_iter = HE_pair.find(std::make_pair(v_index,next_v_index));
        if(HE_pair.empty()){
            HE_pair.insert(std::make_pair(std::make_pair(next_v_index,v_index),i));
        }
        else if(map_iter != HE_pair.end()){
            flag = 0;
            //the other half exist
            pair_e_index = map_iter->second;
            edge_iterators[pair_e_index] ->setPairEdge(edge_iterators[i]);
            edge_iterators[i] ->setPairEdge(edge_iterators[pair_e_index]);
            edge_iterators[i] ->m_vert_it->m_edge_it = edge_iterators[pair_e_index];
            vert_iterators[v_index]->m_edge_it = edge_iterators[pair_e_index];
        }
        else{
            flag = 1;
            HE_pair.insert(std::make_pair(std::make_pair(next_v_index,v_index),i));
        }
        //std::cout << i <<"   " <<flag<< ": v_index:" <<v_index<<"  face_index:"<<i/3<<"   pair_e_index:"<<pair_e_index<<"  next_e_index:"<<next_e_index<<std::endl;
    }
    //std::cout<<"face:"<<std::endl;
    for(T_IDX i = 0 ; i<num_faces ; i++){
        face_iterators[i]->m_edge_it = edge_iterators[i*3];
        //std::cout<<i<<": "<<i*3<<std::endl;
    }
    //std::cout<<"vert:"<<std::endl;
    /*
    for(T_IDX i = 0 ; i<num_verts ; i++){
        T_IDX first_e_index;
        first_e_index = HE_pair[i].second;
        vert_iterators[i]->m_edge_it = edge_iterators[first_e_index];
        //std::cout<<i<<": "<<first_e_index<<std::endl;
    }
    */
    //for (FaceIterator face_it2 = face_list.begin(); face_it2 != face_list.end(); face_it2++) {
        //face_it2->m_edge_it->printEdge();}

}

template<typename DType>
void Mesh<DType>::singleEdgeCollaspe(EdgeIterator target_it){
    VertIterator vert_it1 = target_it->m_vert_it;
    VertIterator vert_it2 = target_it->m_pair_it->m_vert_it;
    //test
    VertIterator vert_temp;
    vert_temp = vert_it1;
    vert_it1 = vert_it2;
    vert_it2 = vert_temp;

    Point NewPt;
    NewPt= (vert_it1->get_coords()+ vert_it2->get_coords())/2;
    //add new vert,edges,faces
    Vert<DType> Newvert = Vert<DType>(NewPt);
    vert_list.push_back(Newvert);
    VertIterator Newvert_it ;
    Newvert_it = prev(vert_list.end());//add new vert, edge to be connected

    std::vector<EdgeIterator> tmpEdgeIts;
    std::vector<FaceIterator> tmpFaceIts;
    std::vector<VertIterator> tmpVertIts;
    EdgeIterator edge_it;
    VertIterator vert_nearby;
    FaceIterator face_it;

    tmpVertIts.clear();
    tmpFaceIts.clear();
    tmpEdgeIts.clear();

    T_IDX i_tmp = 0;
    //first point
    edge_it = target_it->m_pair_it->m_next_it;
    do{
        //std::cout << edge_it->m_vert_it->m_id<<std::endl;
        //locate surrounding verts
        vert_nearby = edge_it->m_pair_it->m_vert_it;
        tmpVertIts.push_back(vert_nearby);
        //add new edges
        edge_list.push_back(Edge<DType>());//in
        edge_list.push_back(Edge<DType>());//out
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts[i_tmp*2] = prev(prev(edge_list.end()));
        tmpEdgeIts[i_tmp*2 + 1] = prev(edge_list.end());
        tmpEdgeIts[i_tmp*2 + 1] ->setNext(edge_it->m_next_it);//set next edge for out edges
        //add new faces
        face_list.push_back(Face<DType>());
        tmpFaceIts.push_back(face_it);
        tmpFaceIts[i_tmp] = prev(face_list.end());

        i_tmp++;
        edge_it = edge_it->m_pair_it->m_next_it;

    }while (edge_it != target_it);   
    tmpEdgeIts[1] ->m_next_it = tmpEdgeIts[1] ->m_next_it ->m_pair_it ->m_next_it;
    //second point
    edge_it = target_it->m_next_it;
    edge_it = edge_it->m_pair_it->m_next_it;
    do{
        std::cout << edge_it->m_vert_it->m_id<<std::endl;
        //locate surrounding verts
        vert_nearby = edge_it->m_pair_it->m_vert_it;
        tmpVertIts.push_back(vert_nearby);
        //add new edges
        edge_list.push_back(Edge<DType>());//in
        edge_list.push_back(Edge<DType>());//out
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts.push_back(edge_it);
        tmpEdgeIts[i_tmp*2] = prev(prev(edge_list.end()));
        tmpEdgeIts[i_tmp*2 + 1] = prev(edge_list.end());
        tmpEdgeIts[i_tmp*2 + 1] ->setNext(edge_it->m_next_it);
        //std::cout<<"set next:  "<<tmpEdgeIts[i_tmp*2+1]->m_next_it->m_vert_it->m_id <<std::endl;
        //add new faces
        face_list.push_back(Face<DType>());
        tmpFaceIts.push_back(face_it);
        tmpFaceIts[i_tmp] = prev(face_list.end());

        i_tmp++;
        edge_it = edge_it->m_pair_it->m_next_it;

        }while (edge_it->m_pair_it->m_next_it != target_it->m_pair_it);//the last vert has already been discovered

    //initalize edges and faces
    T_IDX tmp_vert_num = tmpVertIts.size();
    //std::cout << tmpEdgeIts.size()<<" "<<tmp_vert_num<<"    "<<tmpFaceIts.size()<<"    ";
    for(T_IDX tmp_i = 0 ; tmp_i<tmp_vert_num ; tmp_i++){
        tmpEdgeIts[tmp_i*2]->setIncidencyAndAdjacency(tmpVertIts[tmp_i],tmpFaceIts[tmp_i],tmpEdgeIts[tmp_i*2+1],tmpEdgeIts[(tmp_i*2+3+tmp_vert_num*2) % (tmp_vert_num*2)]);
        tmpEdgeIts[tmp_i*2+1]->setIncidency(Newvert_it,tmpFaceIts[(tmp_i+tmp_vert_num-1)%tmp_vert_num],tmpEdgeIts[tmp_i*2]);//next_it already set
        tmpFaceIts[tmp_i]->m_edge_it = tmpEdgeIts[tmp_i*2];
        //std::cout<<"tmp_i:  "<<tmp_i<<" ";
        tmpVertIts[tmp_i]->m_edge_it = tmpEdgeIts[tmp_i*2+1];
        tmpEdgeIts[tmp_i*2+1]->m_next_it->m_pair_it->m_face_it = tmpFaceIts[tmp_i];//the faces & next_it for surrounding edges have changed
        tmpEdgeIts[tmp_i*2+1]->m_next_it->m_next_it = tmpEdgeIts[(tmp_i*2-2+tmp_vert_num*2)% (tmp_vert_num*2)];
    }
    //tmpEdgeIts[5]->m_next_it = tmpEdgeIts[(2+tmp_vert_num*2)% (tmp_vert_num*2)];

    Newvert_it->m_edge_it = tmpEdgeIts[0];
    //vert_list.back().m_edge_it = tmpEdgeIts[1];
    //Newvert.m_edge_it = tmpEdgeIts[1];

    //delete points,edges,faces
    //first point
    
    edge_it = target_it->m_pair_it;
    edge_it = edge_it->m_pair_it->m_next_it;
    do{
        //face_list.erase(edge_it->m_face_it);
        edge_it->m_face_it->inval();
        //std::cout<<tmpEdgeIts[1]->m_next_it->m_vert_it->m_id;
        EdgeIterator tmp_edge = edge_it;
        edge_it = edge_it->m_pair_it->m_next_it;
        //edge_list.erase(tmp_edge->m_pair_it);
        //tmp_edge->m_pair_it->inval();
        //edge_list.erase(tmp_edge);
        //tmp_edge->inval();
    }while (edge_it != target_it->m_pair_it);

    //second point
    //std::cout<<"pl2"<<std::endl;
    edge_it = target_it->m_pair_it->m_next_it;
    do{
        //face_list.erase(edge_it->m_face_it);
        edge_it->m_face_it->inval();
        //std::cout<<tmpEdgeIts[1]->m_next_it->m_vert_it->m_id;
        EdgeIterator tmp_edge = edge_it;
        edge_it = edge_it->m_pair_it->m_next_it;
        //edge_list.erase(tmp_edge->m_pair_it);
        //tmp_edge->m_pair_it->inval();
        //edge_list.erase(tmp_edge);
        //tmp_edge->inval();
        
    }while (edge_it != target_it);
    //edge_list.erase(target_it->m_pair_it);
    target_it->m_pair_it->inval();
    //edge_list.erase(target_it);
    target_it->inval();
    //std::cout<<"pl3"<<std::endl;
    //vert_list.erase(vert_it1);
    vert_it1->inval();
    //vert_list.erase(vert_it2);
    vert_it2->inval();
    
    //for debug
    /*
    linearlyIndexVerts(1);
    for(int i = 0 ; i<tmp_vert_num ;i++){
        std::cout<<"Edge:"<<tmpEdgeIts[i*2]->m_vert_it->m_id <<" "<<tmpEdgeIts[i*2]->m_pair_it->m_vert_it->m_id<<"  "<<tmpEdgeIts[i*2]->m_next_it->m_vert_it->m_id<<std::endl;
        std::cout<<"Edge:"<<tmpEdgeIts[i*2+1]->m_vert_it->m_id <<" "<<tmpEdgeIts[i*2+1]->m_pair_it->m_vert_it->m_id<<"  "<<tmpEdgeIts[i*2+1]->m_next_it->m_vert_it->m_id<<std::endl;
        std::cout<<"Face:"<<tmpFaceIts[i]->m_edge_it->m_vert_it->m_id<<std::endl;
    }
    */
    
}

template<typename DType>
void Mesh<DType>::insertion(FaceIterator target_face){
    EdgeIterator tmpEdgeIts[9];
    FaceIterator tmpFaceIts[3];
    EdgeIterator edge_it;
    edge_it = target_face->m_edge_it;
    for(int i = 0; i<3 ;i++){
        tmpEdgeIts[i] = edge_it;
        edge_it = edge_it -> m_next_it;
    }
    VertIterator v0 = tmpEdgeIts[0]->m_vert_it;
    VertIterator v1 = tmpEdgeIts[1]->m_vert_it;
    VertIterator v2 = tmpEdgeIts[2]->m_vert_it;
    Point NewPt = (v0->get_coords()+v1->get_coords()+v2->get_coords())/3;

    //add new vert,edges,faces
    Vert<DType> Newvert = Vert<DType>(NewPt);
    vert_list.push_back(Newvert);
    VertIterator Newvert_it ;
    Newvert_it = prev(vert_list.end());//add new vert, edge to be connected

    edge_list.push_back(Edge<DType>());
    edge_list.push_back(Edge<DType>());
    edge_list.push_back(Edge<DType>());
    edge_list.push_back(Edge<DType>());
    edge_list.push_back(Edge<DType>());
    edge_list.push_back(Edge<DType>());
    tmpEdgeIts[6] = prev(prev(prev(edge_list.end())));
    tmpEdgeIts[7] = prev(prev(edge_list.end()));
    tmpEdgeIts[8] = prev(edge_list.end());
    tmpEdgeIts[3] = prev(prev(prev(tmpEdgeIts[6])));
    tmpEdgeIts[4] = prev(prev(tmpEdgeIts[6]));
    tmpEdgeIts[5] = prev(tmpEdgeIts[6]);

    face_list.push_back(Face<DType>());
    face_list.push_back(Face<DType>());
    face_list.push_back(Face<DType>());
    tmpFaceIts[0] = prev(prev(prev(face_list.end())));
    tmpFaceIts[1] = prev(prev(face_list.end()));
    tmpFaceIts[2] = prev(face_list.end());

    //setIncidencyAndAdjacency
    tmpEdgeIts[3]->setIncidencyAndAdjacency(tmpEdgeIts[0]->m_vert_it , tmpFaceIts[2] , tmpEdgeIts[6] , tmpEdgeIts[8]);
    tmpEdgeIts[4]->setIncidencyAndAdjacency(tmpEdgeIts[1]->m_vert_it , tmpFaceIts[0] , tmpEdgeIts[7] , tmpEdgeIts[6]);
    tmpEdgeIts[5]->setIncidencyAndAdjacency(tmpEdgeIts[2]->m_vert_it , tmpFaceIts[1] , tmpEdgeIts[8] , tmpEdgeIts[7]);
    tmpEdgeIts[6]->setIncidencyAndAdjacency(Newvert_it, tmpFaceIts[0] , tmpEdgeIts[3],tmpEdgeIts[0]);
    tmpEdgeIts[7]->setIncidencyAndAdjacency(Newvert_it, tmpFaceIts[1] , tmpEdgeIts[4],tmpEdgeIts[1]);
    tmpEdgeIts[8]->setIncidencyAndAdjacency(Newvert_it, tmpFaceIts[2] , tmpEdgeIts[5],tmpEdgeIts[2]);

    tmpFaceIts[0]->m_edge_it = tmpEdgeIts[0];
    tmpFaceIts[1]->m_edge_it = tmpEdgeIts[1];
    tmpFaceIts[2]->m_edge_it = tmpEdgeIts[2];
    Newvert_it->m_edge_it = tmpEdgeIts[3];

    tmpEdgeIts[0]->m_next_it = tmpEdgeIts[4];
    tmpEdgeIts[1]->m_next_it = tmpEdgeIts[5];
    tmpEdgeIts[2]->m_next_it = tmpEdgeIts[3];
    tmpEdgeIts[0]->m_face_it = tmpFaceIts[0];
    tmpEdgeIts[1]->m_face_it = tmpFaceIts[1];
    tmpEdgeIts[2]->m_face_it = tmpFaceIts[2];

    //delete face
    target_face->inval();

}

template<typename DType>
void Mesh<DType>::EdgeFlip(EdgeIterator target_edge){
    EdgeIterator tmpEdgeIts[6];
    FaceIterator tmpFaceIts[2];
    EdgeIterator edge_it;
    tmpEdgeIts[0] = target_edge->m_next_it;
    tmpEdgeIts[1] = target_edge->m_next_it->m_next_it;
    tmpEdgeIts[2] = target_edge->m_pair_it->m_next_it;
    tmpEdgeIts[3] = target_edge->m_pair_it->m_next_it->m_next_it;

    VertIterator v0 , v1 , v2 , v3;
    v0 = tmpEdgeIts[0]->m_vert_it;
    v1 = tmpEdgeIts[1]->m_vert_it;
    v2 = tmpEdgeIts[2]->m_vert_it;
    v3 = tmpEdgeIts[3]->m_vert_it;

    //add new faces & edges
    edge_list.push_back(Edge<DType>());
    edge_list.push_back(Edge<DType>());
    tmpEdgeIts[4] = prev(prev(edge_list.end()));
    tmpEdgeIts[5] = prev(edge_list.end());

    face_list.push_back(Face<DType>());
    face_list.push_back(Face<DType>());
    tmpFaceIts[0] = prev(prev(face_list.end()));
    tmpFaceIts[1] = prev(face_list.end());

    //set Incidency And Adjacency
    tmpEdgeIts[4]->setIncidencyAndAdjacency(v1,tmpFaceIts[0],tmpEdgeIts[5],tmpEdgeIts[3]);
    tmpEdgeIts[5]->setIncidencyAndAdjacency(v3,tmpFaceIts[1],tmpEdgeIts[4],tmpEdgeIts[1]);

    tmpFaceIts[0]->m_edge_it = tmpEdgeIts[0];
    tmpFaceIts[1]->m_edge_it = tmpEdgeIts[1];

    v0->m_edge_it = tmpEdgeIts[3];
    v1->m_edge_it = tmpEdgeIts[0];
    v2->m_edge_it = tmpEdgeIts[1];
    v3->m_edge_it = tmpEdgeIts[2];

    tmpEdgeIts[0]->m_face_it = tmpFaceIts[0];
    tmpEdgeIts[1]->m_face_it = tmpFaceIts[1];
    tmpEdgeIts[2]->m_face_it = tmpFaceIts[1];
    tmpEdgeIts[3]->m_face_it = tmpFaceIts[0];

    tmpEdgeIts[0]->m_next_it = tmpEdgeIts[4];
    tmpEdgeIts[1]->m_next_it = tmpEdgeIts[2];
    tmpEdgeIts[2]->m_next_it = tmpEdgeIts[5];
    tmpEdgeIts[3]->m_next_it = tmpEdgeIts[0];

    //delete old faces
    target_edge->m_face_it->inval();
    target_edge->m_pair_it->m_face_it->inval();

}





template<typename DType>
void Mesh<DType>::exportOBJ(std::string obj_fn) {

    // file stream
    std::ofstream fout(obj_fn);
    
    // assign ids to verts and output their coordinates
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++) {
        if(vert_it->isval()){
            Point coords = vert_it->m_coords;
            fout << "v " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
        }
    }

    // OBJ subscripts starts from 1
    linearlyIndexVerts(1);

    //for (FaceIterator face_it2 = face_list.begin(); face_it2 != face_list.end(); face_it2++) {
        //face_it2->m_edge_it->printEdge();}

    // output faces
    for (FaceIterator face_it = face_list.begin(); face_it != face_list.end(); face_it++) {
        if(face_it->isval()){
            T_IDX vid_0, vid_1, vid_2;

            EdgeIterator edge_it = face_it->m_edge_it;
            vid_0 = edge_it->m_vert_it->m_id;  edge_it = edge_it->m_next_it;
            vid_1 = edge_it->m_vert_it->m_id;  edge_it = edge_it->m_next_it;
            vid_2 = edge_it->m_vert_it->m_id;

            fout << "f " << vid_0 << " " << vid_1 << " " << vid_2 << "\n";
        }
    }
}


template<typename DType>
void Mesh<DType>::printInfo(){
    for (VertIterator vert_it = vert_list.begin(); vert_it != vert_list.end(); vert_it++) {
        Point coords = vert_it->m_coords;
        std::cout << "v " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }
    linearlyIndexVerts(1);
    for (FaceIterator face_it = face_list.begin(); face_it != face_list.end(); face_it++) {
        T_IDX vid_0, vid_1, vid_2,vid_00, vid_11, vid_22;

        EdgeIterator edge_it = face_it->m_edge_it;
        vid_0 = edge_it->m_vert_it->m_id;  vid_00 = edge_it->m_pair_it->m_vert_it->m_id;
        edge_it = edge_it->m_next_it;
        vid_1 = edge_it->m_vert_it->m_id;  vid_11 = edge_it->m_pair_it->m_vert_it->m_id;
        edge_it = edge_it->m_next_it;
        vid_2 = edge_it->m_vert_it->m_id;  vid_22 = edge_it->m_pair_it->m_vert_it->m_id;

        std::cout << "f " << vid_0 << " " << vid_1 << " " << vid_2 << "\n";
        std::cout << "pf " << vid_00 << " " << vid_11 << " " << vid_22 << "\n";//pair edge
    }

}

#endif
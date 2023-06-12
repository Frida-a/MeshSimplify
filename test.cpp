#include <iostream>
#include <unordered_map>
#include <utility>

typedef std::pair<std::string,std::string> pair;

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

int main()
{
    std::unordered_map<pair,int,pair_hash> HE_pair ;
    

    auto map_iter = HE_pair.find(1,2);
        if(map_iter != HE_pair.end()){
            int flag = 0;
            //the other half exist

            break;
        }
        else{
            flag = 1;
            std::vector<std::pair<T_IDX,T_IDX>> tmp_pair;
            tmp_pair.push_back(std::make_pair(next_v_index,v_index));
            HE_pair.insert(std::make_pair(tmp_pair,i));
            tmp_pair.clear();
        }
    return 0;
}
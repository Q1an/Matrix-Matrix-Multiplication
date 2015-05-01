#include <cstdlib>
#include <iostream>
#include <list>

const int N = 3000;

void quicksort(std::list<int> &data);

int main(){

    std::list<int> data;

    std::list<int>::iterator it;
    for(int i=0; i<N; i++){
        data.push_back(rand());
    }

    /*
    for(it=data.begin(); it!=data.end(); it++){
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    */

    quicksort(data);

    /*
    for(it=data.begin(); it!=data.end(); it++){
        std::cout << *it << " ";
    }
    std::cout << std::endl;
    */

    std::cout << "Quick sort completed." << std::endl;

    return 0;
}

void quicksort(std::list<int> &data){
    
    if(data.size()<=1){
        return;
    }

    int pivot;
    std::list<int> greater;
    std::list<int> less;

    pivot=data.front();
    data.pop_front();

    std::list<int>::iterator it;
    for(it=data.begin(); it!=data.end(); it++){
        if((*it)>pivot){
            greater.push_back(*it);
        }
        else{
            less.push_back(*it);
        }
    }

    quicksort(greater);
    quicksort(less);

    data.clear();
    for(it=less.begin(); it!=less.end(); it++){
        data.push_back(*it);
    }

    data.push_back(pivot);

    for(it=greater.begin(); it!=greater.end(); it++){
        data.push_back(*it);
    }

    return;
}

